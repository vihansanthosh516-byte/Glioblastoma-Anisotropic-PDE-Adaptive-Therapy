#!/usr/bin/env python3
"""Script 102: subgroup, catastrophic-error and residual analysis of the baseline ladder (Phase 3.2, 3.6; plan A1.13-A1.15).

Runs on the CSV outputs of script 98 only (no PDE, light on CPU). The PDE arms join once script 100 has finished.
Subgroups use ONLY information available before the forecast (plan A1.13, A1.14):
  prior-growth stratum  c = (V_t - V_prev) / V_prev from the two scans before the forecast (native-resolution core voxels):
                        shrink c < -0.10, stable -0.10 <= c <= +0.10, growth c > +0.10, no_history (first forecast or V_prev = 0)
  radiotherapy timing   days from the end of radiotherapy to the INPUT scan: pre_or_during (<= 0), early (1-84),
                        late (> 84), unknown (no radiation dates). Analysis-defined; not a RANO category.
  horizon bin           short < 60 d, medium 60-120 d, long > 120 d (pre-set)
  diagnosis             GBM vs other (plan A2.1)
Delta = Dice(method) - Dice(persistence), ONE value per patient per stratum (patient mean over the stratum's pairs),
patient bootstrap (2,000 draws; the registered primary run uses 10,000).
Catastrophic error (analysis-defined research threshold, not clinical): Dice < 0.1 on a pair whose persistence Dice >= 0.1;
sensitivity 0.05 and 0.2. Cluster bootstrap by patient.
Residuals: Spearman between the geometric-baseline gain and covariates known at forecast time; the observed growth
ln(V_out/V_in) is shown for explanation only and is outcome-defined.

Input : output/baseline_ladder_pairs.csv, data/manifests/*.csv
Output: output/subgroups.json, output/subgroups_pairs.csv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
MAN = PROJECT_ROOT / "data" / "manifests"
sys.path.insert(0, str(PROJECT_ROOT))

from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

METHODS = ["geometric_train_rate", "last_rate", "linear", "gompertz", "oracle_volume_matched"]
N_BOOT = 2000


def growth_stratum(c):
    if c is None or not np.isfinite(c):
        return "no_history"
    return "shrink" if c < -0.10 else ("stable" if c <= 0.10 else "growth")


def rt_stratum(rel):
    if rel is None or not np.isfinite(rel):
        return "unknown_rt"
    return "pre_or_during" if rel <= 0 else ("early_1_84" if rel <= 84 else "late_gt84")


def horizon(dt):
    return "short_lt60" if dt < 60 else ("medium_60_120" if dt <= 120 else "long_gt120")


def build_table():
    lad = pd.read_csv(OUT / "baseline_ladder_pairs.csv")
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")[["patient_id", "tp_in", "tp_out", "day_in"]]
    pats = pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)
    pats = pats[pats["dataset"] == "MU"].set_index("patient_id")
    scans = pd.read_csv(MAN / "scan_manifest.csv", low_memory=False)
    scans = scans[scans["dataset"] == "MU"].sort_values(["patient_id", "day"])
    prev_vol = {}
    for pid, g in scans.groupby("patient_id"):
        g = g.reset_index(drop=True)
        for i in range(len(g)):
            prev_vol[(pid, int(g.loc[i, "scan_number"]))] = (
                (float(g.loc[i, "core_voxels_native"]), float(g.loc[i - 1, "core_voxels_native"])) if i > 0 else (float(g.loc[i, "core_voxels_native"]), np.nan))
    per = lad[lad["method"] == "persistence"][["patient_id", "tp_in", "tp_out", "dt_days", "in_primary", "is_gbm", "v_in", "v_out", "dice", "has_prior_scan"]]
    per = per.rename(columns={"dice": "dice_persistence"}).merge(pairs, on=["patient_id", "tp_in", "tp_out"], how="left")
    c = []
    for r in per.itertuples(index=False):
        vt, vp = prev_vol.get((r.patient_id, int(r.tp_in)), (np.nan, np.nan))
        c.append((vt - vp) / vp if np.isfinite(vp) and vp > 0 else np.nan)
    per["prior_growth_c"] = c
    per["growth_stratum"] = per["prior_growth_c"].map(growth_stratum)
    per["rt_end_day"] = per["patient_id"].map(pats["radiation_end_day"])
    per["rt_rel_days"] = per["day_in"] - per["rt_end_day"]
    per["rt_stratum"] = per["rt_rel_days"].map(rt_stratum)
    per["horizon"] = per["dt_days"].map(horizon)
    per["log_obs_growth"] = np.log((per["v_out"] + 1) / (per["v_in"] + 1))
    per["log_v_in"] = np.log(per["v_in"] + 1)
    wide = lad.pivot_table(index=["patient_id", "tp_in", "tp_out"], columns="method", values="dice").reset_index()
    return per.merge(wide, on=["patient_id", "tp_in", "tp_out"], how="left")


def stratum_block(df):
    out = {"n_pairs": int(len(df)), "n_patients": int(df["patient_id"].nunique()),
           "persistence_mean_dice": float(ps.patient_means(df, "dice_persistence").mean()) if len(df) else None, "methods": {}}
    for m in METHODS:
        d = df.copy()
        d["delta"] = d[m] - d["dice_persistence"]
        pm = ps.patient_means(d, "delta")
        if len(pm) < 8:
            continue
        s = ps.summarize_delta(pm.to_numpy(), n_boot=N_BOOT)
        out["methods"][m] = {k: s[k] for k in ("n", "mean", "mean_ci95", "median", "pct_improved", "pct_worse", "broadly_positive")}
    return out


def catastrophic(df, thr):
    sub = df[df["dice_persistence"] >= 0.1]
    res = {"threshold": thr, "n_pairs_eligible": int(len(sub))}
    for m in ["persistence"] + METHODS:
        col = "dice_persistence" if m == "persistence" else m
        flag = (sub[col] < thr).astype(float)
        t = sub.assign(flag=flag)
        lo, hi = ps.cluster_bootstrap(t, lambda d: d["flag"].mean(), n_boot=N_BOOT)
        res[m] = {"rate": float(flag.mean()), "ci95": [lo, hi]}
    return res


def main():
    df = build_table()
    df.to_csv(OUT / "subgroups_pairs.csv", index=False)
    prim = df[df["in_primary"]]
    res = {"script": "102_subgroups", "n_primary_pairs": int(len(prim)), "n_patients": int(prim["patient_id"].nunique()), "strata": {}}
    for kind in ("growth_stratum", "rt_stratum", "horizon"):
        res["strata"][kind] = {k: stratum_block(g) for k, g in prim.groupby(kind)}
        res["strata"][kind + "_gbm_only"] = {k: stratum_block(g) for k, g in prim[prim["is_gbm"] == True].groupby(kind)}  # noqa: E712
    res["strata"]["diagnosis"] = {"gbm": stratum_block(prim[prim["is_gbm"] == True]), "non_gbm": stratum_block(prim[prim["is_gbm"] == False])}  # noqa: E712
    res["catastrophic_error"] = {str(t): catastrophic(prim, t) for t in (0.1, 0.05, 0.2)}
    prim = prim.assign(delta_geo=prim["geometric_train_rate"] - prim["dice_persistence"])
    cov = {}
    for name in ("log_v_in", "dt_days", "prior_growth_c", "rt_rel_days", "log_obs_growth"):
        sub = prim.dropna(subset=[name, "delta_geo"])
        if len(sub) < 20:
            continue
        rho = float(stats.spearmanr(sub[name], sub["delta_geo"]).statistic)
        lo, hi = ps.cluster_bootstrap(sub, lambda d, n=name: stats.spearmanr(d[n], d["delta_geo"]).statistic, n_boot=N_BOOT)
        cov[name] = {"spearman_with_geometric_gain": rho, "ci95_cluster_bootstrap": [lo, hi], "n_pairs": int(len(sub)),
                     "known_before_forecast": name != "log_obs_growth"}
    res["residual_covariates"] = cov
    g = ps.patient_means(prim.assign(delta=prim["delta_geo"]), "delta")
    res["geometric_gain_patient_quantiles"] = {str(q): float(g.quantile(q)) for q in (0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)}
    (OUT / "subgroups.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("subgroups_102", OUT / "subgroups.manifest.json", script="src/102_subgroups.py", seed=ps.DEFAULT_SEED,
                       config={"n_boot": N_BOOT}, inputs=[OUT / "baseline_ladder_pairs.csv", MAN / "forecast_pairs_mu.csv",
                                                           MAN / "patient_manifest.csv", MAN / "scan_manifest.csv"],
                       dataset="MU-Glioma-Post core", patient_split="data/manifests/split_mu.csv", primary_endpoint="n/a (subgroup)")
    for kind in ("growth_stratum", "rt_stratum", "horizon"):
        print("\n", kind)
        for k, b in res["strata"][kind].items():
            gm = b["methods"].get("geometric_train_rate")
            print(f"  {k:16s} pat {b['n_patients']:3d} pers {b['persistence_mean_dice']:.3f} " + (
                f"geo {gm['mean']:+.3f} [{gm['mean_ci95'][0]:+.3f},{gm['mean_ci95'][1]:+.3f}] med {gm['median']:+.3f}" if gm else "(n<8)"))
    print("\ncatastrophic (thr 0.1):", {m: round(v["rate"], 3) for m, v in res["catastrophic_error"]["0.1"].items() if isinstance(v, dict)})
    print("covariates:", {k: round(v["spearman_with_geometric_gain"], 3) for k, v in cov.items()})


if __name__ == "__main__":
    main()
