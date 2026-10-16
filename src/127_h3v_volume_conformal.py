#!/usr/bin/env python3
"""Script 127: H-3v volume forecast with split-conformal intervals, MU development (Amendment 10 A10.1-A10.3).

Target: y = ln(V_out / V_in) of the core (labels 1+3, 2 mm) on MU primary pairs.
Arms (out-of-fold on manifest folds, refit per outer fold on outer-train only):
  persistence (0), geometric_train_rate (script-98 v_pred, no refit), MixedLM (primary).
Split-conformal: outer-train patients split 75/25 proper-train/calibration
(seed 20261004, per-fold RandomState(20261004+f)); q80/q90 from calibration
patient-mean |residuals|; coverage on the outer-fold test patients.
Populations: GBM-only (primary, 105 as analysed) + all-eligible (sensitivity).
Exclusions: v_in<=0 or v_out<=0 (ln undefined), listed in output.
Output: output/h3v_dev.json, output/h3v_dev_pairs.csv
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from src import patient_stats as ps  # noqa: E402

SEED = 20261004
CAL_FRAC = 0.25
# Measured floor, ledger 2ab (input-volume tertile edges in mm3; MU voxels are
# 8 mm3 at 2 mm, so convert voxel counts to mm3 before binning).
FLOOR_EDGES = (2059.0, 7614.0)
FLOORS = (1.03, 0.88, 0.35)


def floor_for_vin(v_mm3: float) -> float:
    if v_mm3 < FLOOR_EDGES[0]:
        return FLOORS[0]
    if v_mm3 < FLOOR_EDGES[1]:
        return FLOORS[1]
    return FLOORS[2]


def load_pairs() -> pd.DataFrame:
    b = pd.read_csv(PROJECT_ROOT / "output" / "baseline_ladder_pairs.csv")
    p = b[(b["in_primary"] == True)].copy()  # noqa: E712
    pers = p[p["method"] == "persistence"][
        ["patient_id", "tp_in", "tp_out", "fold", "dt_days", "is_gbm",
         "has_prior_scan", "v_in", "v_out"]].copy()
    geo = p[p["method"] == "geometric_train_rate"][
        ["patient_id", "tp_in", "tp_out", "v_pred"]].copy()
    df = pers.merge(geo, on=["patient_id", "tp_in", "tp_out"], how="left")
    f = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "forecast_pairs_mu.csv",
                    usecols=["patient_id", "tp_in", "tp_out", "day_in"])
    df = df.merge(f, on=["patient_id", "tp_in", "tp_out"], how="left")
    m = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "patient_manifest.csv")
    mu = m[m["dataset"] == "MU"][["patient_id", "radiation_end_day", "age", "sex",
                                  "mgmt_code"]].copy()
    df = df.merge(mu, on="patient_id", how="left")
    return df


def design(df: pd.DataFrame, rt_median: float) -> pd.DataFrame:
    d = df.copy()
    d["y"] = np.log(d["v_out"] / d["v_in"])
    d["yhat_persist"] = 0.0
    d["yhat_geo"] = np.log(d["v_pred"] / d["v_in"])
    d["logVin"] = np.log(d["v_in"])
    d["rt_known"] = d["radiation_end_day"].notna().astype(float)
    d["daysRT"] = (d["day_in"] - d["radiation_end_day"]).where(
        d["radiation_end_day"].notna(), rt_median)
    d["sexM"] = (d["sex"] == "Male").astype(float)
    d["mgmt_meth"] = d["mgmt_code"].map({1: 1.0, 0: 0.0})
    d["mgmt_miss"] = d["mgmt_meth"].isna().astype(float)
    d["mgmt_meth"] = d["mgmt_meth"].fillna(0.0)
    d["const"] = 1.0
    return d


FIXED = ["const", "logVin", "dt_days", "daysRT", "rt_known", "age", "sexM",
         "mgmt_meth", "mgmt_miss"]


def fit_mixed(train: pd.DataFrame):
    import statsmodels.api as sm
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            md = sm.MixedLM(train["y"], train[FIXED], groups=train["patient_id"],
                            exog_re=train[["const", "dt_days"]])
            res = md.fit(method="lbfgs", maxiter=200, disp=0)
            if not np.all(np.isfinite(res.params)):
                raise RuntimeError("non-finite params")
            return ("mixedlm", res)
        except Exception as e:  # noqa: BLE001 - declared fallback, logged
            ols = sm.OLS(train["y"], train[FIXED]).fit()
            return (f"ols_fallback:{type(e).__name__}", ols)


def predict(fit, df: pd.DataFrame) -> np.ndarray:
    _, res = fit
    return np.asarray(res.predict(df[FIXED]))


def patient_means(df: pd.DataFrame, col: str) -> pd.Series:
    return df.groupby("patient_id")[col].mean()


def main() -> None:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--write-freeze", action="store_true",
                    help="refit MixedLM on ALL MU primary pairs and write "
                         "configs/h3v_frozen.yaml (no LUMIERE read)")
    args = ap.parse_args()
    if args.write_freeze:
        write_freeze()
        return
    run_dev()


def run_dev() -> None:
    df = load_pairs()
    bad = df[(df["v_in"] <= 0) | (df["v_out"] <= 0)].copy()
    df = df[(df["v_in"] > 0) & (df["v_out"] > 0)].copy()
    rt_median = float(
        (df.drop_duplicates("patient_id")
           .assign(drt=lambda t: t["day_in"] - t["radiation_end_day"])
           ["drt"]).median())
    # NOTE: day_in here is per-pair; use first-pair day per patient for the median
    firsts = df.sort_values(["patient_id", "day_in"]).drop_duplicates("patient_id")
    rt_median = float(((firsts["day_in"] - firsts["radiation_end_day"])
                       .dropna()).median())
    df = design(df, rt_median)

    out: dict = {"script": "127_h3v_volume_conformal", "seed": SEED,
                 "rt_median_days": rt_median,
                 "excluded_nonpositive_volume": bad[
                     ["patient_id", "tp_in", "tp_out", "v_in", "v_out"]
                 ].to_dict("records"),
                 "pops": {}}
    all_rows = []
    for pop, flag in (("gbm_only", True), ("all_eligible", False)):
        sub = df[df["is_gbm"] == True].copy() if flag else df.copy()  # noqa: E712
        folds = sorted(sub["fold"].dropna().unique())
        prow = []
        for f in folds:
            te = sub[sub["fold"] == f]
            tr = sub[sub["fold"] != f]
            pats = np.array(sorted(tr["patient_id"].unique()))
            rng = np.random.RandomState(SEED + int(f))
            rng.shuffle(pats)
            ncal = max(1, int(round(len(pats) * CAL_FRAC)))
            cal_pats, pro_pats = set(pats[:ncal]), set(pats[ncal:])
            pro, cal = tr[tr["patient_id"].isin(pro_pats)], tr[
                tr["patient_id"].isin(cal_pats)]
            fit = fit_mixed(pro)
            for name, frame in (("proper", pro), ("cal", cal), ("test", te)):
                frame = frame.copy()
                frame["yhat_mix"] = predict(fit, frame)
                frame["pop"], frame["outer_fold"], frame["split_role"] = pop, int(f), name
                frame["fit_kind"] = fit[0]
                prow.append(frame)
            # conformal quantiles from calibration patient-mean |residuals|
            cal_r = cal.copy()
            cal_r["yhat_mix"] = predict(fit, cal_r)
            cpat = (cal_r.assign(res=(cal_r["y"] - cal_r["yhat_mix"]).abs())
                    .groupby("patient_id")["res"].mean())
            q80, q90 = float(cpat.quantile(0.80)), float(cpat.quantile(0.90))
            out.setdefault("quantiles", {})[f"{pop}|fold{int(f)}"] = {
                "q80": q80, "q90": q90, "n_cal_patients": int(len(cpat)),
                "fit_kind": fit[0]}
        P = pd.concat(prow)
        all_rows.append(P)
        # out-of-fold test predictions only
        T = P[P["split_role"] == "test"].copy()
        T["ae_mix"] = (T["y"] - T["yhat_mix"]).abs()
        T["ae_persist"] = T["y"].abs()
        T["ae_geo"] = (T["y"] - T["yhat_geo"]).abs()
        pm = pd.DataFrame({
            "mix": patient_means(T, "ae_mix"),
            "persist": patient_means(T, "ae_persist"),
            "geo": patient_means(T, "ae_geo"),
            "ymean": T.groupby("patient_id")["y"].mean(),
            "yhatmean": T.groupby("patient_id")["yhat_mix"].mean(),
        })
        pm["d_mix_persist"] = pm["mix"] - pm["persist"]
        pm["d_mix_geo"] = pm["mix"] - pm["geo"]
        res = {"n_patients": int(len(pm)),
               "mean_ae": {k: float(pm[k].mean()) for k in ("mix", "persist", "geo")},
               "mix_minus_persist": ps.summarize_delta(pm["d_mix_persist"].to_numpy()),
               "mix_minus_geo": ps.summarize_delta(pm["d_mix_geo"].to_numpy())}
        # coverage: patient interval from fold quantiles
        cov = {}
        for level, qk in (("c80", "q80"), ("c90", "q90")):
            hit = []
            for pid, r in pm.iterrows():
                f = int(T[T["patient_id"] == pid]["outer_fold"].iloc[0])
                q = out["quantiles"][f"{pop}|fold{f}"][qk]
                hit.append(abs(float(r["ymean"]) - float(r["yhatmean"])) <= q)
            cov[level] = {"coverage": float(np.mean(hit)), "n": int(len(hit))}
        res["conformal_oof"] = cov
        # abstain subgroup (descriptive): |yhat| below measured floor
        T["floor"] = T["v_in"].apply(lambda v: floor_for_vin(float(v) * 8.0))
        T["abstain"] = T["yhat_mix"].abs() < T["floor"]
        pa = T.groupby("patient_id")["abstain"].mean()
        for lab, sel in (("forecast", pa[pa < 0.5].index),
                         ("abstain_majority", pa[pa >= 0.5].index)):
            dd = pm.loc[pm.index.intersection(sel), "d_mix_persist"].to_numpy()
            res[f"abstain_{lab}"] = (
                ps.summarize_delta(dd) if len(dd) else {"n": 0})
        out["pops"][pop] = res
    P = pd.concat(all_rows)
    P.to_csv(PROJECT_ROOT / "output" / "h3v_dev_pairs.csv", index=False)
    (PROJECT_ROOT / "output" / "h3v_dev.json").write_text(json.dumps(out, indent=1))
    print(json.dumps({k: out["pops"][k]["mean_ae"] for k in out["pops"]}, indent=1))
    print(json.dumps({k: {kk: out["pops"][k][kk]["mean_ci95"]
                          for kk in ("mix_minus_persist",)} for k in out["pops"]},
                     indent=1))


def write_freeze() -> None:
    """Refit MixedLM on ALL MU primary pairs (GBM-only, A5.1) and write the
    frozen config. Reads MU only; no LUMIERE access."""
    import yaml
    df = load_pairs()
    df = df[(df["v_in"] > 0) & (df["v_out"] > 0)].copy()
    firsts = df.sort_values(["patient_id", "day_in"]).drop_duplicates("patient_id")
    rt_median = float(((firsts["day_in"] - firsts["radiation_end_day"])
                       .dropna()).median())
    df = design(df, rt_median)
    gbm = df[df["is_gbm"] == True].copy()  # noqa: E712
    kind, res = fit_mixed(gbm)
    # pooled calibration quantiles from the dev run (already on disk)
    dev = json.loads((PROJECT_ROOT / "output" / "h3v_dev.json").read_text())
    qs = [v for k, v in dev["quantiles"].items() if k.startswith("gbm_only|")]
    q80 = float(np.mean([v["q80"] for v in qs]))
    q90 = float(np.mean([v["q90"] for v in qs]))
    cfg = {"model": "mixedlm_gbm_all_pairs",
           "script": "127_h3v_volume_conformal --write-freeze",
           "fixed_effects": FIXED,
           "params": {k: float(v) for k, v in zip(FIXED, np.asarray(res.params))},
           "fit_kind": kind,
           "rt_median_days": rt_median,
           "q80": q80, "q90": q90,
           "floor_edges_mm3": list(FLOOR_EDGES), "floors": list(FLOORS),
           "lumiere_rules": {
               "daysRT_missing": "RTmiss=1, daysRT=<rt_median_days>",
               "mgmt_map": {"methylated": 1.0, "not methylated": 0.0,
                            "other_or_missing": "0.0 with mgmt_miss=1"},
               "sex_map": {"male": 1.0, "female": 0.0},
               "volume": "contrast-enhancing CT1 VoxelVolume, HD-GLIO-AUTO "
                         "primary / DeepBraTumIA sensitivity"}}
    (PROJECT_ROOT / "configs" / "h3v_frozen.yaml").write_text(yaml.dump(cfg))
    print("wrote configs/h3v_frozen.yaml", kind)


if __name__ == "__main__":
    main()
