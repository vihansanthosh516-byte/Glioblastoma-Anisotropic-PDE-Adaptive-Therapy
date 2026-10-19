#!/usr/bin/env python3
"""Script 126: H-3 development scoring of candidates (Amendment 9, A9.3 and A9.5). Header originally said C2 / C4 / C5; the saved
output contains arms C1, C2, C2_C3, C1C3, C3 and unet, and Amendment 11 dropped C4/C5 (the declared comparison is C1/C2/C3).

Members are out-of-fold maps per TUM dev patient:
  - learned arms of script 115: <UCSF_H3_DIR>/../colab_h3/out/oof_logits_<arm>.npz (2 mm logits, upsampled to 1 mm);
  - "unet": the released PREDICT-GBM U-Net, leak-free (only the 3 checkpoints whose fold lists the patient as TEST;
    script 125 --stage dev).
Each map is rank-normalised inside the brain (script 120 rule); an ensemble is the mean of its members' ranks.
C5 (uncertainty margin): mean rank + 1.0 x standard deviation of the members' ranks (lambda fixed, not tuned).
Plan = top-k inside the brain at the 15 mm standard-plan volume (script 117 rank_plan).
Endpoint as A9.4: enhancing-recurrence coverage minus the standard plan; patients with enhancing recurrence.
Output: output/predictgbm/h3/dev_ensembles.json, dev_ensembles_pairs.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm" / "h3"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
OOF = _Path(os.environ.get("H3_OOF_DIR", r"C:\Users\vihan\data\colab_h3\out"))
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402


def _mod(name, file):
    s = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "src" / file)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


s112 = _mod("s112", "112_predictgbm_models_dev.py")
s117 = _mod("s117", "117_predictgbm_volume_curves.py")
s120 = _mod("s120", "120_predictgbm_ensembles.py")


def main():
    from src import patient_stats as ps
    ap = argparse.ArgumentParser()
    ap.add_argument("--members", required=True, help="comma list of script-115 arms; 'unet' = leak-free U-Net")
    ap.add_argument("--arms", required=True, help="semicolon list name=member+member[+...][:c5]")
    a = ap.parse_args()
    mem = a.members.split(",")
    Z = {m: np.load(H3 / "unet_dev_logits.npz") if m == "unet" else np.load(OOF / f"oof_logits_{m}.npz") for m in mem}
    arms = {}
    for spec in a.arms.split(";"):
        name, rest = spec.split("=")
        c5 = rest.endswith(":c5")
        arms[name] = (rest.replace(":c5", "").split("+"), c5)
    rows = []
    for pid in sorted(pg.dev_ids()):
        if not all(pid in Z[m].files for m in mem):
            continue
        c = pg.load_case(pid, with_recurrence=True)
        re_ = pg.rec_enhancing(c["rec"])
        if re_.sum() == 0 or not pg.core(c["seg"]).any():
            continue
        std = pg.standard_plan(c["seg"], c["brain"])
        k = int(std.sum())
        R = {}
        for m in mem:
            v = Z[m][pid].astype(np.float32)
            R[m] = s120.ranks(v if m == "unet" else s112.up1(v), c["brain"])
        r = {"pid": pid, "standard_enh": pg.coverage(re_, std)}
        for name, (ms, c5) in arms.items():
            st = np.stack([R[m] for m in ms])
            s = st.mean(0) + (st.std(0) if c5 else 0)
            r[name] = pg.coverage(re_, s117.rank_plan(s, k, c["brain"]))
        rows.append(r)
        print(pid, {n: round(r[n], 3) for n in arms}, flush=True)
    df = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / "dev_ensembles_pairs.csv", index=False)
    res = {"script": "126_h3_dev_ensembles", "declared": "Amendment 9 (A9.3, A9.5)", "split": "TUM development, out-of-fold",
           "n_patients": int(len(df)), "standard_enh_mean": float(df["standard_enh"].mean()), "arms": {}}
    for n in arms:
        res["arms"][n] = {"members": arms[n][0], "c5": arms[n][1], "enh_mean": float(df[n].mean()),
                          "delta_vs_standard": ps.summarize_delta((df[n] - df["standard_enh"]).to_numpy())}
    res["weakest_points"] = ["Development patients only; the dev number of the selected arm is optimistic.",
                             "U-Net member: trained on TUM patients of other folds and on LUMIERE/RHUH/others (not UCSF)."]
    (OUT / "dev_ensembles.json").write_text(json.dumps(res, indent=1, default=float))
    for n, v in res["arms"].items():
        d = v["delta_vs_standard"]
        print(f"{n:28s} {v['enh_mean']:.4f}  d {d['mean']:+.4f} {[round(x, 4) for x in d['mean_ci95']]}")


if __name__ == "__main__":
    main()
