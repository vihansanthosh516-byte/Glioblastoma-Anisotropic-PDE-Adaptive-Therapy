#!/usr/bin/env python3
"""Script 133: H-3 locked UCSF test, run ONCE after the freeze (configs/h3_frozen.yaml + git tag h3-frozen).

Amendments 9 and 11. Reads the frozen config, checks the SHA-256 of the final logits (script 132) and of the input pack
(script 131), then opens the ALPTDG labels for the first time. Per patient (PREDICT-GBM array order; ALPTDG arrays are
flipped in x and y as script 124): standard plan = EDT <= 15 mm from the pre-op core inside the aligned pre-op brain;
model plan = predictgbm_io.model_plan of the upsampled sigmoid logits at the standard-plan volume.
Primary (A9.4, A11.2): enhancing coverage (time1 label 3), final minus standard, patients with label 3 > 0 at time1;
patient bootstrap 95% CI, two-sided Wilcoxon p; success = lower CI bound > 0.
Secondary: all recurrence (time1 labels 1-3); 90% volume vs standard (margin -0.01); time2 targets (exploratory);
context: the released U-Net (all 15 checkpoints, none trained on UCSF; script 125) at equal volume.
Output: output/predictgbm/h3_final/ucsf_test_results.json, ucsf_test_pairs.csv
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import subprocess
import sys
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
import yaml
from scipy.stats import wilcoxon

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "h3_final"
FROZEN = PROJECT_ROOT / "configs" / "h3_frozen.yaml"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
ZIP = _Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402
from src import patient_stats as ps  # noqa: E402


def sha(p):
    return hashlib.sha256(_Path(p).read_bytes()).hexdigest()


def is_frozen():
    if not FROZEN.exists():
        return False
    r = subprocess.run(["git", "tag", "-l", "h3-frozen"], cwd=PROJECT_ROOT, capture_output=True, text=True)
    return r.stdout.strip() == "h3-frozen"


def up1(a):
    u = np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)
    return np.concatenate([u, u[:, :, -1:]], 2)


def label(z, sid, tp):
    a = np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_{tp}_seg.nii.gz"))).dataobj)
    return np.ascontiguousarray(np.rint(a).astype(np.uint8)[::-1, ::-1])


def main():
    if not is_frozen():
        raise SystemExit("not frozen: configs/h3_frozen.yaml and tag h3-frozen are required")
    cfg = yaml.safe_load(FROZEN.read_text())
    lf = OUT / "ucsf_logits_2mm.npz"
    if sha(lf) != cfg["final_model"]["ucsf_logits_sha256"]:
        raise SystemExit("logits hash does not match the frozen config")
    if sha(PROJECT_ROOT / "kaggle_run" / "pack_ucsf" / "ucsf_h3_inputs_2mm.npz") != cfg["inputs"]["pack_sha256"]:
        raise SystemExit("input pack hash does not match the frozen config")
    L = np.load(lf)
    U = np.load(H3 / "unet_ucsf_logits.npz")
    z = zipfile.ZipFile(ZIP)
    rows = []
    for key in sorted(L.files):
        sid = key.replace("ucsf_", "")
        x = np.load(H3 / "inputs" / f"{sid}.npz")
        seg, brain = np.rint(x["seg"]).astype(np.uint8), x["brain"] > 0
        std = pg.standard_plan(seg, brain)
        k = int(std.sum())
        score = 1.0 / (1.0 + np.exp(-up1(L[key].astype(np.float32))))
        plan = pg.model_plan(score, seg, brain, k=k)
        plan90 = pg.topk_plan(score, int(round(0.9 * k)), brain)
        uplan = pg.topk_plan(U[sid].astype(np.float32), k, brain)
        r = {"pid": sid, "k": k}
        for tp in ("time1", "time2"):
            lab = label(z, sid, tp)
            re_, ra = (lab == 3) & brain, np.isin(lab, (1, 2, 3)) & brain
            r[f"{tp}_n_enh"] = int(re_.sum())
            r[f"{tp}_standard_enh"], r[f"{tp}_final_enh"] = pg.coverage(re_, std), pg.coverage(re_, plan)
            r[f"{tp}_final90_enh"], r[f"{tp}_unet_enh"] = pg.coverage(re_, plan90), pg.coverage(re_, uplan)
            r[f"{tp}_standard_all"], r[f"{tp}_final_all"] = pg.coverage(ra, std), pg.coverage(ra, plan)
        rows.append(r)
        print(sid, round(r["time1_standard_enh"], 3) if r["time1_n_enh"] else "-", round(r["time1_final_enh"], 3) if r["time1_n_enh"] else "-", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "ucsf_test_pairs.csv", index=False)

    def contrast(d, a, b):
        v = (d[a] - d[b]).dropna().to_numpy()
        s = ps.summarize_delta(v)
        s["wilcoxon_p_two_sided"] = float(wilcoxon(v).pvalue) if np.any(v != 0) else 1.0
        return s

    p1 = df[df["time1_n_enh"] > 0]
    p2 = df[df["time2_n_enh"] > 0]
    res = {"script": "133_h3_ucsf_test_score", "frozen_config_sha256": sha(FROZEN), "n_test": int(len(df)),
           "n_time1_enhancing": int(len(p1)), "n_time2_enhancing": int(len(p2)),
           "means_time1": {c: float(p1[c].mean()) for c in p1.columns if c.startswith("time1_") and c.endswith("_enh") and "n_" not in c},
           "PRIMARY_final_minus_standard_enh_time1": contrast(p1, "time1_final_enh", "time1_standard_enh"),
           "secondary": {"all_recurrence_time1": contrast(df[df["time1_standard_all"].notna()], "time1_final_all", "time1_standard_all"),
                         "volume90_vs_standard100_time1": contrast(p1, "time1_final90_enh", "time1_standard_enh"),
                         "time2_enh_exploratory": contrast(p2, "time2_final_enh", "time2_standard_enh"),
                         "context_unet_minus_standard_time1": contrast(p1, "time1_unet_enh", "time1_standard_enh")}}
    lo = res["secondary"]["volume90_vs_standard100_time1"]["mean_ci95"][0]
    res["secondary"]["volume90_noninferior_margin_minus_0.01"] = bool(lo > -0.01)
    res["PRIMARY_success_lower_ci_above_0"] = bool(res["PRIMARY_final_minus_standard_enh_time1"]["mean_ci95"][0] > 0)
    res["weakest_points"] = ["Time1 enhancement may be residual or treatment-related, not only recurrence (no RANO label).",
                             "Rigid pre-op to time1 alignment; resection shift not modelled.",
                             "UCSF tissue = dev-mean atlas (A11.4); single centre; n with enhancing tumour may be small."]
    (OUT / "ucsf_test_results.json").write_text(json.dumps(res, indent=1, default=float))
    p = res["PRIMARY_final_minus_standard_enh_time1"]
    print(json.dumps(res["means_time1"], indent=1))
    print("PRIMARY", round(p["mean"], 4), p["mean_ci95"], "p", p["wilcoxon_p_two_sided"], "n", p["n"])


if __name__ == "__main__":
    main()
