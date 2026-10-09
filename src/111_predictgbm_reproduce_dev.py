#!/usr/bin/env python3
"""Script 111: reproduce the PREDICT-GBM evaluation on DEVELOPMENT patients only (TUM; Amendment 6 A6.6, A7.6).

Checks that src/predictgbm_io.py rebuilds the released standard plan and model plans voxel for voxel, then reports
development coverage for the standard plan and each published model. Test patients are never read here (the loader
lock would refuse their recurrence anyway). This gives the bar our models must clear on development data.
Output: output/predictgbm/dev_reproduce.json, dev_coverage_pairs.csv (+ manifest)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ids = pg.dev_ids()
    rows, plan_mismatch = [], {}
    for k, pid in enumerate(ids):
        c = pg.load_case(pid, with_recurrence=True)
        std = pg.standard_plan(c["seg"], c["brain"])
        rel_std = pg._load(pid, "standard_plan.nii.gz") > 0
        plan_mismatch.setdefault("standard", []).append(int((std ^ rel_std).sum()))
        re, ra = pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])
        r = {"pid": pid, "n_rec_enh": int(re.sum()), "n_rec_all": int(ra.sum()), "plan_voxels": int(std.sum()),
             "standard_enh": pg.coverage(re, std), "standard_all": pg.coverage(ra, std)}
        for m in pg.PUBLISHED_MODELS:
            f = pg.DATA / pid / f"{m}_pred.nii.gz"
            if not f.exists():
                continue
            plan = pg.model_plan(pg._load(pid, f"{m}_pred.nii.gz"), c["seg"], c["brain"])
            rel = pg.DATA / pid / f"{m}_plan.nii.gz"
            if rel.exists():
                relp = pg._load(pid, f"{m}_plan.nii.gz") > 0
                plan_mismatch.setdefault(m, []).append(int((plan ^ relp).sum()))
                r[f"{m}_released_enh"] = pg.coverage(re, relp)   # plan file as released (what the paper scored)
            r[f"{m}_enh"] = pg.coverage(re, plan)
            r[f"{m}_all"] = pg.coverage(ra, plan)
        rows.append(r)
        if (k + 1) % 20 == 0:
            print(f"{k + 1}/{len(ids)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "dev_coverage_pairs.csv", index=False)
    use = df[df["n_rec_enh"] > 0]
    res = {"script": "111_predictgbm_reproduce_dev", "split": "development (TUM) only", "n_patients": int(len(df)),
           "n_with_enhancing_recurrence": int(len(use)), "n_boot": ps.DEFAULT_N_BOOT,
           "plan_reproduction_voxel_mismatch": {m: {"n": len(v), "n_exact": int(sum(x == 0 for x in v)), "max": int(max(v))}
                                                for m, v in plan_mismatch.items()},
           "coverage_enh_mean": {}, "delta_vs_standard_enh": {}}
    for col in ["standard_enh"] + [f"{m}{x}_enh" for m in pg.PUBLISHED_MODELS for x in ("", "_released")]:
        if col in use:
            res["coverage_enh_mean"][col] = float(use[col].mean())
            if col != "standard_enh":
                d = (use[col] - use["standard_enh"]).dropna().to_numpy()
                res["delta_vs_standard_enh"][col] = ps.summarize_delta(d)
    res["weakest_points"] = ["Development patients only; published models may have been tuned on some of them.",
                             "Coverage is geometric, not dose; retrospective.",
                             "Rebuilt U-Net / LMI / GlioMap plans differ from the released plan files (U-Net scores are raw values "
                             "outside 0-1 and are clipped as in evaluate.py); *_released_enh scores the released files."]
    (OUT / "dev_reproduce.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("predictgbm_dev_111", OUT / "dev_reproduce.manifest.json", script="src/111_predictgbm_reproduce_dev.py",
                       seed=ps.DEFAULT_SEED, config={"ctv_margin": pg.CTV_MARGIN, "models": pg.PUBLISHED_MODELS},
                       inputs=[], dataset="PREDICT-GBM (TUM development)", patient_split="TUM dev; LUMIERE+RHUH test (locked)",
                       primary_endpoint="enhancing-recurrence coverage")
    print(json.dumps({k: res[k] for k in ("n_patients", "n_with_enhancing_recurrence", "plan_reproduction_voxel_mismatch",
                                          "coverage_enh_mean")}, indent=1, default=float))


if __name__ == "__main__":
    main()
