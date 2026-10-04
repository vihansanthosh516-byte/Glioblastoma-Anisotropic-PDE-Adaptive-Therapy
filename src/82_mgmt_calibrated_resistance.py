#!/usr/bin/env python3
"""
Script 82: MGMT-calibrated resistance analysis (Track C negative 3, real-data follow-up).
=========================================================================================
Script 76 found that adaptive therapy gains nothing from resistance, using an assumed
resistant fraction f_r0 in {0.001, 0.01, 0.1}. Here the resistant fraction is tied to data:
MGMT promoter methylation, the clinical marker of temozolomide sensitivity, from the
MU-Glioma-Post clinical sheet.

PRE-SPECIFIED ANALYSIS (fixed before the first run; not changed after seeing data)
---------------------------------------------------------------------------------
Hypothesis  MGMT-unmethylated (MGMT = 0) tumours progress faster than MGMT-methylated
            (MGMT = 1) tumours.
Cohort      MU Glioma Post clinical sheet: Primary Diagnosis == "GBM", MGMT in {0, 1}
            (codes 2 indeterminate, 3/4 not assessed are excluded). Subgroups: MGMT = 0 and
            MGMT = 1, nothing else.
Event       Progression == 1, time = "Time to First Progression (Days)" (from diagnosis /
            resection). Progression == 0 is censored at the last known day: the later of the
            last MRI day in output/mu_glioma_cohort.json and the day of death if recorded.
            Patients with no usable time are listed, not imputed.
PRIMARY     median time to progression by MGMT status (Kaplan-Meier, censoring-aware, 95% CI),
            and the two-sided log-rank p-value. Also reported: the median among progressors
            only (the biased quantity that gave 167 vs 182 days in the unrestricted sheet),
            Cox hazard ratio (MGMT 0 vs 1) unadjusted and adjusted for age, and restricted
            mean TTP to day 365 (RMST_365) with bootstrap CI (2000 resamples, seed 82).
SECONDARY   Does a resistance model calibrated to the MGMT difference make adaptive therapy
            beat non-adaptive therapy?
  Calibration  Anchor f_r0(MGMT = 1) = 0.01 (script 76's primary cell, Strobl 2021). Find
            f_r0(MGMT = 0) on the grid {0.01, 0.03, 0.1, 0.3, 0.6} (log-linear interpolation)
            at which the model's RMST_365 ratio, Stupp schedule, pooled script-76 evaluation
            sets (111 patients), cost 0.25, eps 0, script-59 seed, equals the observed
            RMST_365 ratio MGMT 0 / MGMT 1. Three calibrations: observed ratio, and the
            lower and upper bootstrap bounds of the observed ratio. If the observed ratio is
            >= 1 the calibrated f_r0 equals the anchor (no excess resistance in the data).
  Test      Script-76 primary-cell settings (cost 0.25, eps 0, script-59 seed, 365-day
            window) at each calibrated f_r0 (cells cached from script 76 where identical).
  FLIPS POSITIVE iff, on real_test (the only real-patient set), at the point-estimate
            calibrated f_r0 for MGMT = 0, either AT50 or AT80 has (a) CI lower bound > 0 for
            the difference vs the same-window paced arm AND (b) resistance-attributable days
            (difference-in-differences vs the f_r0 = 0 control) > 0. All sets, cells and arms
            are reported either way. If it does not flip, the result stays a limitation.
Output: output/mgmt_resistance/{results.json, km_mgmt_ttp.png, cells/}
"""
from __future__ import annotations

import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUT_DIR = PROJECT_ROOT / "output" / "mgmt_resistance"
CLIN = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post_ClinicalData-July2025.xlsx"
COHORT_JSON = PROJECT_ROOT / "output" / "mu_glioma_cohort.json"
SEED = 82
ANCHOR = 0.01
F_GRID = (0.01, 0.03, 0.1, 0.3, 0.6)
RMST_T = 365.0
PRIMARY_CELL = {"cost": 0.25, "eps": 0.0, "occupancy": "seed59", "window": 365}


def clinical_table() -> tuple[pd.DataFrame, list]:
    x = pd.read_excel(CLIN, sheet_name="MU Glioma Post")
    cohort = {p["patient_id"]: p for p in json.loads(COHORT_JSON.read_text())}
    num = lambda c: pd.to_numeric(x[c], errors="coerce")  # noqa: E731
    df = pd.DataFrame({"id": x["Patient_ID"], "dx": x["Primary Diagnosis"], "mgmt": num("MGMT methylation"),
                       "prog": num("Progression"), "t_prog": num("Time to First Progression (Days)"),
                       "death_day": num("Number of days from Diagnosis to death (Days)"),
                       "age": num("Age at diagnosis")})
    last_scan = []
    for pid in df["id"]:
        days = [t.get("day_from_diagnosis") for t in cohort.get(pid, {}).get("timepoints", [])]
        days = [float(d) for d in days if d is not None and np.isfinite(float(d))]
        last_scan.append(max(days) if days else np.nan)
    df["last_scan_day"] = last_scan
    df = df[(df["dx"] == "GBM") & df["mgmt"].isin([0, 1])].copy()
    excluded = []
    t, e = [], []
    for _, r in df.iterrows():
        if r["prog"] == 1 and np.isfinite(r["t_prog"]) and r["t_prog"] > 0:
            t.append(r["t_prog"]); e.append(1)
        elif r["prog"] == 0:
            c = np.nanmax([r["last_scan_day"], r["death_day"]]) if np.isfinite([r["last_scan_day"], r["death_day"]]).any() else np.nan
            t.append(c); e.append(0)
        else:
            t.append(np.nan); e.append(np.nan)
    df["time"], df["event"] = t, e
    bad = df[~(df["time"] > 0)]
    excluded = [(r["id"], "no usable time") for _, r in bad.iterrows()]
    return df[df["time"] > 0].copy(), excluded


def km_block(df: pd.DataFrame) -> dict:
    from lifelines import CoxPHFitter, KaplanMeierFitter
    from lifelines.statistics import logrank_test
    from lifelines.utils import median_survival_times, restricted_mean_survival_time
    g0, g1 = df[df["mgmt"] == 0], df[df["mgmt"] == 1]
    out = {"n": {"mgmt0": len(g0), "mgmt1": len(g1)}, "events": {"mgmt0": int(g0["event"].sum()),
                                                                 "mgmt1": int(g1["event"].sum())}}
    fits = {}
    for name, g in (("mgmt0", g0), ("mgmt1", g1)):
        k = KaplanMeierFitter().fit(g["time"], g["event"], label=name)
        ci = median_survival_times(k.confidence_interval_)
        fits[name] = k
        out[f"km_median_{name}"] = float(k.median_survival_time_)
        out[f"km_median_{name}_ci95"] = [float(ci.iloc[0, 0]), float(ci.iloc[0, 1])]
        out[f"rmst365_{name}"] = float(restricted_mean_survival_time(k, t=RMST_T))
        prog = g[g["event"] == 1]["time"]
        out[f"median_progressors_only_{name}"] = float(prog.median())
    lr = logrank_test(g0["time"], g1["time"], g0["event"], g1["event"])
    out["logrank_p"] = float(lr.p_value)
    out["logrank_chi2"] = float(lr.test_statistic)
    cox = CoxPHFitter().fit(df[["time", "event", "mgmt"]].assign(mgmt0=(df["mgmt"] == 0).astype(float))
                            [["time", "event", "mgmt0"]], "time", "event")
    s = cox.summary.loc["mgmt0"]
    out["cox_hr_mgmt0_vs_mgmt1"] = {"hr": float(s["exp(coef)"]), "ci95": [float(s["exp(coef) lower 95%"]),
                                                                           float(s["exp(coef) upper 95%"])],
                                    "p": float(s["p"])}
    d2 = df.dropna(subset=["age"])[["time", "event", "age"]].assign(mgmt0=(df["mgmt"] == 0).astype(float))
    d2 = d2.dropna()
    cox2 = CoxPHFitter().fit(d2, "time", "event")
    s2 = cox2.summary.loc["mgmt0"]
    out["cox_adjusted_age"] = {"n": len(d2), "hr": float(s2["exp(coef)"]),
                               "ci95": [float(s2["exp(coef) lower 95%"]), float(s2["exp(coef) upper 95%"])],
                               "p": float(s2["p"])}
    rng = np.random.default_rng(SEED)
    ratios = []
    for _ in range(2000):
        r0, r1 = g0.sample(len(g0), replace=True, random_state=int(rng.integers(1 << 31))), \
            g1.sample(len(g1), replace=True, random_state=int(rng.integers(1 << 31)))
        a = restricted_mean_survival_time(KaplanMeierFitter().fit(r0["time"], r0["event"]), t=RMST_T)
        b = restricted_mean_survival_time(KaplanMeierFitter().fit(r1["time"], r1["event"]), t=RMST_T)
        ratios.append(a / b)
    out["rmst365_ratio_mgmt0_over_mgmt1"] = float(out["rmst365_mgmt0"] / out["rmst365_mgmt1"])
    out["rmst365_ratio_ci95"] = [float(np.percentile(ratios, 2.5)), float(np.percentile(ratios, 97.5))]
    return out, fits


def plot_km(fits: dict, out: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    for name, c in (("mgmt0", "firebrick"), ("mgmt1", "navy")):
        fits[name].plot_survival_function(ax=ax, color=c, ci_alpha=0.15)
    ax.set_xlabel("days from diagnosis")
    ax.set_ylabel("progression-free probability")
    ax.set_title(f"MGMT 0 (n={out['n']['mgmt0']}) vs 1 (n={out['n']['mgmt1']}), log-rank p = {out['logrank_p']:.3f}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "km_mgmt_ttp.png", dpi=150)
    plt.close(fig)


def calibrate(obs_ratio: float, s76, pooled: dict) -> dict:
    """Model RMST_365 ratio vs f_r0 (anchor 0.01), then invert at obs_ratio (clamped to the grid)."""
    import torch
    rho, kill = pooled["rho"], pooled["kill"]
    sched = s76.schedule(torch.tensor([s76.STUPP_SCHEDULE] * len(rho)))
    rm = {}
    for f in F_GRID:
        sim = s76.TwoPopSim(rho.tolist(), kill, f, PRIMARY_CELL["cost"], PRIMARY_CELL["eps"], PRIMARY_CELL["occupancy"])
        rm[f] = float(s76.rollout(sim, sched, 90)["ttp"].mean())
    ratio = {f: rm[f] / rm[ANCHOR] for f in F_GRID}
    fs, rs = np.log(np.array(F_GRID)), np.array([ratio[f] for f in F_GRID])
    order = np.argsort(rs)  # ratio decreases with f_r0; interpolate on a monotone grid
    return {"model_rmst365": rm, "model_ratio_vs_anchor": ratio, "_fs": fs, "_rs": rs, "_order": order}


def invert(cal: dict, target: float) -> float:
    fs, rs = cal["_fs"], cal["_rs"]
    if target >= rs[0]:
        return float(np.exp(fs[0]))            # no excess resistance needed
    if target <= rs.min():
        return float(np.exp(fs[int(np.argmin(rs))]))   # beyond the grid: clamp to the most resistant value
    return float(np.exp(np.interp(target, rs[::-1], fs[::-1])))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "cells").mkdir(exist_ok=True)
    df, excluded = clinical_table()
    km, fits = km_block(df)
    plot_km(fits, km)
    print(json.dumps(km, indent=2))

    sp = spec_from_file_location("s76", PROJECT_ROOT / "src" / "76_resistance_adaptive_equal_budget.py")
    s76 = module_from_spec(sp)
    sp.loader.exec_module(s76)
    import torch
    sets = s76.s68.reaction_sets()
    pooled = {"rho": torch.cat([s["rho"] for s in sets.values()]),
              "kill": np.concatenate([np.asarray(s["kill"], float) for s in sets.values()])}
    cal = calibrate(km["rmst365_ratio_mgmt0_over_mgmt1"], s76, pooled)
    obs = km["rmst365_ratio_mgmt0_over_mgmt1"]
    lo, hi = km["rmst365_ratio_ci95"]
    f_point, f_low_ratio, f_high_ratio = invert(cal, obs), invert(cal, hi), invert(cal, lo)
    f_cal = {"point": f_point, "at_ratio_upper_ci": f_low_ratio, "at_ratio_lower_ci": f_high_ratio}
    print("calibrated f_r0 (MGMT0):", f_cal, "model ratios:", cal["model_ratio_vs_anchor"])

    # adaptive comparison at calibrated f_r0 (cells reuse script-76 cells when the config matches)
    results = {}
    fixed_cache = {}
    needed = {"anchor_mgmt1": ANCHOR, **{f"mgmt0_{k}": v for k, v in f_cal.items()}}
    ctrl_cfg = {"f_r0": 0.0, "cost": 0.0, "eps": 0.0, "occupancy": "seed59", "window": 365}
    ctrl = json.loads(s76.cell_path(ctrl_cfg).read_text())
    for tag, f in needed.items():
        cfg = {"f_r0": float(f), **PRIMARY_CELL}
        p76 = s76.cell_path(cfg)
        cell_file = OUT_DIR / "cells" / f"{tag}.json"
        if cell_file.exists():
            cell = json.loads(cell_file.read_text())
        elif p76.exists():
            cell = json.loads(p76.read_text())
        else:
            cell = s76.run_cell(7000 + len(results), cfg, sets, fixed_cache)
        cell_file.write_text(json.dumps(cell))
        rows = {}
        for set_name, st in cell["sets"].items():
            for at in ("AT50", "AT80"):
                a, a0 = st["arms"][at], ctrl["sets"][set_name]["arms"][at]
                rows[f"{set_name}:{at}"] = {
                    "vs_paced_days": a["vs_paced"]["mean_diff_days"], "vs_paced_ci95": a["vs_paced"]["ci95"],
                    "vs_paced_control_days": a0["vs_paced"]["mean_diff_days"],
                    "resistance_attributable_days": a["vs_paced"]["mean_diff_days"] - a0["vs_paced"]["mean_diff_days"],
                    "vs_best_nonadaptive_days": a["vs_best_nonadaptive"]["mean_diff_days"],
                    "vs_best_nonadaptive_ci95": a["vs_best_nonadaptive"]["ci95"],
                    "resistant_fraction_day365_median": a["resistant_fraction_day365_median"]}
        results[tag] = {"f_r0": float(f), "rows": rows}
    pt = results["mgmt0_point"]["rows"]
    flips = any(pt[f"real_test:{at}"]["vs_paced_ci95"][0] > 0 and pt[f"real_test:{at}"]["resistance_attributable_days"] > 0
                for at in ("AT50", "AT80"))
    s76_primary = {s: json.loads((PROJECT_ROOT / "output" / "resistance_adaptive" / "evaluate.json").read_text())["summary"][s]["rows"]
                   for s in ("real_test", "cohort64")}
    primary76 = {s: [r for r in rows if r["f_r0"] == 0.01 and r["cost"] == 0.25 and r["eps"] == 0.0
                     and r["occupancy"] == "seed59" and r["window"] == 365] for s, rows in s76_primary.items()}
    res = {"script": "82_mgmt_calibrated_resistance", "cohort_excluded": excluded, "clinical": km,
           "calibration": {"anchor_f_r0_mgmt1": ANCHOR, "grid": F_GRID, "observed_rmst365_ratio": obs,
                           "observed_ratio_ci95": [lo, hi],
                           "model_rmst365": cal["model_rmst365"], "model_ratio_vs_anchor": cal["model_ratio_vs_anchor"],
                           "f_r0_mgmt0": f_cal},
           "adaptive_comparison": results, "flips_positive_per_prespecified_rule": bool(flips),
           "script76_primary_cell_rows": primary76}
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k not in ("clinical", "script76_primary_cell_rows")}, indent=2))
    print(f"[saved] {OUT_DIR / 'results.json'}  flips_positive={flips}")


if __name__ == "__main__":
    main()
