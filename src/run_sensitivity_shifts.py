#!/usr/bin/env python3
"""Exploratory sensitivity on H-2 locked test outputs (does NOT touch frozen model).

Reads frozen test logits + recurrence via the locked loader (freeze already in
place), rebuilds standard/final plans exactly as script 123, then reports:
  - coverage under +2 / +5 voxel (mm) shifts of the recurrence mask,
  - distant-component split (>30 mm from pre-op core surface, A7.4),
  - local-only coverage contrast.
Writes output/predictgbm/final/sensitivity_shifts.csv + .json (new files).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt, label as ndlabel
from scipy.ndimage import shift as ndshift

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from src import patient_stats as ps  # noqa: E402

OUT = PROJECT_ROOT / "output" / "predictgbm" / "final"


def up1(a):
    u = np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)
    return np.concatenate([u, u[:, :, -1:]], 2)


def cov(m, p):
    return 1.0 if m.sum() == 0 else float((m & p).sum() / m.sum())


def main():
    z = np.load(OUT / "test_logits_2mm.npz")
    df = pd.read_csv(OUT / "test_pairs.csv")
    use = df[df["n_rec_enh"] > 0].reset_index(drop=True)
    rows = []
    for i, r in use.iterrows():
        pid = r["pid"]
        c = pg.load_case(pid, with_recurrence=True)
        std = pg.standard_plan(c["seg"], c["brain"])
        k = int(std.sum())
        score = 1.0 / (1.0 + np.exp(-up1(z[pid].astype(np.float32))))
        plan = pg.model_plan(score, c["seg"], c["brain"], k=k)
        re = pg.rec_enhancing(c["rec"])
        lab, nlab = ndlabel(re)
        dcore = distance_transform_edt(~pg.core(c["seg"]))
        distant = np.zeros_like(re, bool)
        for li in range(1, nlab + 1):
            comp = lab == li
            cent = np.argwhere(comp).mean(axis=0).round().astype(int)
            cent = np.clip(cent, 0, np.array(re.shape) - 1)
            if dcore[tuple(cent)] > 30:
                distant |= comp
        local = re & ~distant
        row = {
            "pid": pid, "nlab": int(nlab),
            "ndist_vox": int(distant.sum()), "nrec": int(re.sum()),
            "base_S": cov(re, std), "base_F": cov(re, plan),
            "local_S": cov(local, std), "local_F": cov(local, plan),
        }
        for s in (2, 5):
            rs = ndshift(re.astype(np.float32), (s, 0, 0), order=0).astype(bool)
            row[f"shift{s}_S"] = cov(rs, std)
            row[f"shift{s}_F"] = cov(rs, plan)
            row[f"shift{s}_d"] = row[f"shift{s}_F"] - row[f"shift{s}_S"]
        row["base_d"] = row["base_F"] - row["base_S"]
        row["local_d"] = row["local_F"] - row["local_S"]
        rows.append(row)
        if (i + 1) % 10 == 0:
            print(f"{i + 1}/{len(use)}", flush=True)
    dd = pd.DataFrame(rows)
    dd.to_csv(OUT / "sensitivity_shifts.csv", index=False)

    def summ(x):
        s = ps.summarize_delta(np.asarray(x, float))
        return {kk: s[kk] for kk in
                ("n", "mean", "mean_ci95", "median", "pct_improved",
                 "pct_worse", "wilcoxon_p_two_sided") if kk in s}

    out = {
        "script": "sensitivity_shifts (exploratory, frozen outputs only)",
        "n": len(dd),
        "primary_base": summ(dd["base_d"]),
        "shift2": summ(dd["shift2_d"]),
        "shift5": summ(dd["shift5_d"]),
        "local_only": summ(dd["local_d"]),
        "distant_patients": int((dd["ndist_vox"] > 0).sum()),
        "distant_only_patients": int(((dd["ndist_vox"] == dd["nrec"]) & (dd["nrec"] > 0)).sum()),
        "distant_vox_frac": float(dd["ndist_vox"].sum() / dd["nrec"].sum()),
        "mean_base_S": float(dd["base_S"].mean()),
        "mean_base_F": float(dd["base_F"].mean()),
        "mean_local_S": float(dd["local_S"].mean()),
        "mean_local_F": float(dd["local_F"].mean()),
    }
    (OUT / "sensitivity_shifts.json").write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
