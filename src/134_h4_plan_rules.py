#!/usr/bin/env python3
"""Script 134: H-4 step 1, plan-rule candidates on the pooled development set (Amendment 12). POST HOC: LUMIERE, RHUH
and UCSF were locked tests (H-2, H-3) and are re-used here as development data.

Maps (no retraining): TUM = C1 out-of-fold logits (script 115, colab_h3/out); LUMIERE + RHUH = H-2 frozen logits
(output/predictgbm/final/test_logits_2mm.npz); UCSF = H-3 frozen logits (output/predictgbm/h3_final/ucsf_logits_2mm.npz).
Rules (A12.3): R0 = model top-k; G5 / G10 = core + 5 / 10 mm always treated, model fills the rest; f = fraction of the
15 mm standard-plan volume. Endpoints vs the 15 mm plan at f = 1: enhancing and all-recurrence coverage.
Output: output/predictgbm/h4/plan_rules.json, plan_rules_pairs.csv
"""
from __future__ import annotations

import gzip
import json
import os
import sys
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "h4"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
OOF = _Path(os.environ.get("H3_OOF_DIR", r"C:\Users\vihan\data\colab_h3\out"))
ZIP = _Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402
from src import patient_stats as ps  # noqa: E402

FRACS = (0.7, 0.8, 0.9, 1.0)
RULES = {"R0": 0, "G5": 5, "G10": 10}


def up1(a):
    u = np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)
    return np.concatenate([u, u[:, :, -1:]], 2)


def plan(score, d, brain, guard, k):
    if guard == 0:
        return pg.topk_plan(score, k, brain)
    g = (d <= guard) & brain
    if g.sum() >= k:
        return pg.topk_plan(1.0 / (1.0 + d), k, brain)   # positive score: topk_plan treats <= 0 as missing
    return g | pg.topk_plan(np.where(g, -np.inf, score), k - int(g.sum()), brain & ~g)


def score_case(pid, ds, logit, seg, brain, rec_e, rec_a):
    d = distance_transform_edt(~pg.core(seg))
    std = (d <= pg.CTV_MARGIN) & brain
    k = int(std.sum())
    s = 1.0 / (1.0 + np.exp(-up1(logit.astype(np.float32))))
    s = s[:seg.shape[0], :seg.shape[1], :seg.shape[2]]
    r = {"pid": pid, "dataset": ds, "n_enh": int(rec_e.sum()), "standard_enh": pg.coverage(rec_e, std),
         "standard_all": pg.coverage(rec_a, std)}
    for name, gd in RULES.items():
        for f in FRACS:
            p = plan(s, d, brain, gd, int(round(f * k)))
            r[f"{name}_{f}_enh"], r[f"{name}_{f}_all"] = pg.coverage(rec_e, p), pg.coverage(rec_a, p)
    return r


def main():
    rows = []
    tum = np.load(OOF / "oof_logits_hybrid_mri_rank_mu_s0.npz")
    for pid in sorted(set(pg.dev_ids()) & set(tum.files)):
        c = pg.load_case(pid, with_recurrence=True)
        if not pg.core(c["seg"]).any():
            continue
        rows.append(score_case(pid, "tum", tum[pid], c["seg"], c["brain"], pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])))
        print(pid, flush=True)
    h2 = np.load(PROJECT_ROOT / "output" / "predictgbm" / "final" / "test_logits_2mm.npz")
    for pid in sorted(h2.files):
        c = pg.load_case(pid, with_recurrence=True)
        ds = "rhuh" if pid.upper().startswith("RHUH") else "lumiere"
        rows.append(score_case(pid, ds, h2[pid], c["seg"], c["brain"], pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])))
        print(pid, flush=True)
    h3 = np.load(PROJECT_ROOT / "output" / "predictgbm" / "h3_final" / "ucsf_logits_2mm.npz")
    z = zipfile.ZipFile(ZIP)
    for key in sorted(h3.files):
        sid = key.replace("ucsf_", "")
        x = np.load(H3 / "inputs" / f"{sid}.npz")
        seg, brain = np.rint(x["seg"]).astype(np.uint8), x["brain"] > 0
        lab = np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_time1_seg.nii.gz"))).dataobj)
        lab = np.rint(lab).astype(np.uint8)[::-1, ::-1]
        rows.append(score_case(sid, "ucsf", h3[key], seg, brain, (lab == 3) & brain, np.isin(lab, (1, 2, 3)) & brain))
        print(sid, flush=True)
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "plan_rules_pairs.csv", index=False)
    pe = df[df["n_enh"] > 0]
    res = {"script": "134_h4_plan_rules", "declared": "Amendment 12", "label": "POST HOC: H-2/H-3 test sets re-used as dev",
           "n": {"all": int(len(df)), "with_enh": int(len(pe))}, "by_dataset_n": df["dataset"].value_counts().to_dict(), "rules": {}}
    for name in RULES:
        for f in FRACS:
            key = f"{name}_{f}"
            e = ps.summarize_delta((pe[f"{key}_enh"] - pe["standard_enh"]).to_numpy())
            a = ps.summarize_delta((df[f"{key}_all"] - df["standard_all"]).to_numpy())
            per = {g: {"enh": float((d[f"{key}_enh"] - d["standard_enh"])[d["n_enh"] > 0].mean()),
                       "all": float((d[f"{key}_all"] - d["standard_all"]).mean())} for g, d in df.groupby("dataset")}
            res["rules"][key] = {"enh": e, "all": a, "per_dataset_mean": per}
    eq = [k for k in res["rules"] if k.endswith("_1.0") and res["rules"][k]["all"]["mean_ci95"][0] > 0]
    res["A12.4a_equal_volume_pick"] = max(eq, key=lambda k: res["rules"][k]["enh"]["mean_ci95"][0]) if eq else None
    ok = [k for k, v in res["rules"].items() if v["enh"]["mean_ci95"][0] > -0.01 and v["all"]["mean_ci95"][0] > -0.01
          and all(p["enh"] > -0.01 and p["all"] > -0.01 for p in v["per_dataset_mean"].values())]
    res["A12.4b_reduced_volume_ok"] = sorted(ok, key=lambda k: (float(k.split("_")[1]), list(RULES).index(k.split("_")[0])))
    (OUT / "plan_rules.json").write_text(json.dumps(res, indent=1, default=float))
    for k, v in res["rules"].items():
        print(f"{k:8s} enh {v['enh']['mean']:+.4f} {[round(x, 4) for x in v['enh']['mean_ci95']]}  all {v['all']['mean']:+.4f} {[round(x, 4) for x in v['all']['mean_ci95']]}")
    print("pick", res["A12.4a_equal_volume_pick"], "reduced ok", res["A12.4b_reduced_volume_ok"])


if __name__ == "__main__":
    main()
