#!/usr/bin/env python3
"""Script 120: GRAND_PLAN 15 Z1 (= 13 X3) ensembles, PREDICT-GBM development (TUM).

Declared in GRAND_PLAN section 15 before this script ran. Ten arms are scored on the same development patients and
best_at_equal_volume is chosen by the same-set mean, so it is optimistic (published maps may also have seen TUM). Per patient, each map is rank-normalised inside the brain
(rank / n, ties averaged); an ensemble is the mean of its members' ranks. Plans = top-k inside the brain at
k = fraction x standard-plan volume (fractions as script 117). Members: GliODIL, LOTI (sbtc) (released maps, clipped as
the benchmark), U-Net raw (released prediction ranked by raw value), M1 (our isotropic tissue model at lambda 2 mm,
script 115 m1_density, upsampled), Euclidean (minus distance from the core).
Arms (declared): a GliODIL+LOTI; b GliODIL+U-Net; c GliODIL+LOTI+U-Net; d c+M1; e c+Euclidean; plus each single member.
Output: output/predictgbm/ensembles.json, ensembles_pairs.csv
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt
from scipy.stats import rankdata

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402


def _mod(name, file):
    s = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "src" / file)
    m = importlib.util.module_from_spec(s)
    s.loader.exec_module(m)
    return m


s112 = _mod("s112", "112_predictgbm_models_dev.py")
s115 = _mod("s115", "115_predictgbm_unet_dev.py")
s117 = _mod("s117", "117_predictgbm_volume_curves.py")

FRACS = [0.8, 0.9, 1.0]
ARMS = {"gliodil": ["gliodil"], "sbtc": ["sbtc"], "unet_raw": ["unet_raw"], "m1": ["m1"], "euclid": ["euclid"],
        "a_gliodil_sbtc": ["gliodil", "sbtc"], "b_gliodil_unet": ["gliodil", "unet_raw"],
        "c_gliodil_sbtc_unet": ["gliodil", "sbtc", "unet_raw"], "d_c_plus_m1": ["gliodil", "sbtc", "unet_raw", "m1"],
        "e_c_plus_euclid": ["gliodil", "sbtc", "unet_raw", "euclid"]}


def ranks(score, brain):
    r = np.zeros(score.shape, np.float32)
    v = score[brain].astype(np.float64)
    r[brain] = rankdata(v, method="average") / v.size
    return r


def main():
    from src import patient_stats as ps
    rows = []
    ids = pg.dev_ids()
    for i, pid in enumerate(ids):
        c = pg.load_case(pid, with_recurrence=True)
        core = pg.core(c["seg"])
        re = pg.rec_enhancing(c["rec"])
        if not core.any() or re.sum() == 0:
            continue
        std = pg.standard_plan(c["seg"], c["brain"])
        need = {"gliodil", "sbtc", "unet_raw"}
        if not all((pg.DATA / pid / f"{m.replace('_raw', '')}_pred.nii.gz").exists() for m in need):
            continue
        maps = {"euclid": -distance_transform_edt(~core)}
        for m in ("gliodil", "sbtc"):
            maps[m] = np.clip(pg.resample_to(pg._load(pid, f"{m}_pred.nii.gz").astype(np.float32), c["seg"].shape), 0, 1)
        maps["unet_raw"] = pg.resample_to(pg._load(pid, "unet_pred.nii.gz").astype(np.float32), c["seg"].shape)
        core2 = s112.to2(core) >= 0.5
        maps["m1"] = s112.up1(s115.m1_density(c, core2, std)) if core2.any() else -distance_transform_edt(~core)
        R = {k: ranks(v, c["brain"]) for k, v in maps.items()}
        r = {"pid": pid, "standard_enh": pg.coverage(re, std)}
        for arm, mem in ARMS.items():
            s = np.mean([R[m] for m in mem], 0)
            for f in FRACS:
                r[f"{arm}@{f:g}"] = pg.coverage(re, s117.rank_plan(s, int(round(f * std.sum())), c["brain"]))
        rows.append(r)
        print(f"{i + 1}/{len(ids)} {pid}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "ensembles_pairs.csv", index=False)
    res = {"script": "120_predictgbm_ensembles", "declared": "GRAND_PLAN 15 Z1", "split": "development (TUM)",
           "n_patients": int(len(df)), "standard_enh_mean": float(df["standard_enh"].mean()), "arms": {}}
    for arm in ARMS:
        res["arms"][arm] = {f"{f:g}": {"mean": float(df[f"{arm}@{f:g}"].mean()),
                                       "delta_vs_standard": ps.summarize_delta((df[f"{arm}@{f:g}"] - df["standard_enh"]).to_numpy())}
                            for f in FRACS}
    best = max(ARMS, key=lambda a: res["arms"][a]["1"]["mean"])
    res["best_at_equal_volume"] = best
    res["weakest_points"] = ["Development patients; published members may have seen TUM.",
                             "Ensemble choice among 10 arms on the same patients: the best arm's dev number is optimistic."]
    (OUT / "ensembles.json").write_text(json.dumps(res, indent=1, default=float))
    for arm in ARMS:
        d = res["arms"][arm]
        print(f"{arm:22s}", {f: round(d[f]["mean"], 4) for f in d}, "d@1", round(d["1"]["delta_vs_standard"]["mean"], 4),
              [round(x, 4) for x in d["1"]["delta_vs_standard"]["mean_ci95"]])


if __name__ == "__main__":
    main()
