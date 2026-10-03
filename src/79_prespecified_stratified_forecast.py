#!/usr/bin/env python3
"""
Script 79: Pre-specified stratification forecast (Track B negative 1 follow-up).
================================================================================
Problem: in run_improved_aniso.py (scan 1 -> scan 2, 152 patients) 53% of tumours shrank or
stayed the same. A growth-only model cannot forecast shrinkage and loses to no-change.
Restricting the headline to tumours that GREW is outcome selection (it conditions on the
target) and is not defensible. A defensible version stratifies on information available
BEFORE the forecast. This script does that with a pre-specified rule.

PRE-SPECIFIED DESIGN (fixed before the first run; not changed after seeing data)
-------------------------------------------------------------------------------
Cohort    MU-Glioma-Post patients with >= 3 scans, ordered scan days, scan-volume > 0 on the
          first three scans (output/mu_glioma_cohort.json; the gate counted 108).
Forecast  scan index 1 (second scan) -> scan index 2 (third scan). The history scan 0 -> 1 is
          the only extra information. Forecast pipeline = run_improved_aniso.forecast_patient,
          unchanged (same atlas, same parameter grid, same treatment kill term, same
          5-fold patient-level CV selection of (rho, d, r), all Dice out-of-fold).
Classifier  Predict, before the forecast, whether the tumour grows between scans 1 and 2
          (target: 2 mm voxel count at scan 2 > scan 1, the existing `grew` definition).
          Features (all known at scan 1): log(v1/v0)/(day1-day0) from cohort scan volumes;
          log v1; interval day1-day0; forecast interval day2-day1; days from radiation end to
          scan 1; days from surgery to scan 1; indicators: radiation recorded, TMZ recorded,
          scan 1 during TMZ, scan 1 during radiation. Missing numeric values -> training-fold
          median. Model: standardised L2 logistic regression, C = 1 (no tuning), 5-fold
          patient-level CV (the folds of the forecast CV), decision threshold 0.5.
          Reported: out-of-fold accuracy and AUC against two baselines: always-predict-grow,
          and persistence (grew iff scan 0 -> 1 grew).
Strata    predicted-to-grow (out-of-fold probability >= 0.5) vs predicted-not-to-grow.
Arms      anisotropic, iso_same, iso_homog, no_change (all as in run_improved_aniso).
          Secondary, shape-only (uses the observed scan-2 volume, so not a forecast):
          volume-matched arms and uniform dilation, on the predicted-to-grow stratum.
Primary endpoints (Holm-corrected over the two):
          P1  out-of-fold Dice, anisotropic vs no_change, predicted-to-grow stratum
          P2  out-of-fold Dice, anisotropic vs iso_same, predicted-to-grow stratum
Reported in full (intention-to-forecast): all patients, every arm, mean Dice with bootstrap
          CI, plus strata by prediction and (descriptive only, outcome-selected) by observed
          growth. The predicted-to-grow stratum is NOT the headline for "forecast beats
          no-change"; the all-patient table is.
Stages    forecast (resumable per-patient cache), analyze.

Output: output/forecast_stratified/{forecast_grid.json, results.json}
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
import run_improved_aniso as ria  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "forecast_stratified"
CACHE_DIR = OUT_DIR / "cache"
GRID_JSON = OUT_DIR / "forecast_grid.json"
RESULTS_JSON = OUT_DIR / "results.json"
N_FOLDS = ria.N_FOLDS
PROB_THRESHOLD = 0.5


def f(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def stratified_pairs():
    cohort = json.loads(ria.COHORT_JSON.read_text())
    pairs, excluded = [], []
    for p in cohort:
        tp = p["timepoints"]
        if len(tp) < 3:
            continue
        days = [f(t["day_from_diagnosis"]) for t in tp[:3]]
        vols = [f(t["volume_mm3"]) for t in tp[:3]]
        if any(d is None for d in days) or not (days[0] < days[1] < days[2]) \
                or any(v is None or v <= 0 for v in vols):
            excluded.append((p["patient_id"], "missing/unordered day or non-positive volume in first 3 scans"))
            continue
        pairs.append({"patient_id": p["patient_id"], "tp1": tp[1]["number"], "tp2": tp[2]["number"],
                      "dt_days": days[2] - days[1], "t1_day": days[1], "treatment": p.get("treatment_schedule"),
                      "hist": {"day0": days[0], "day1": days[1], "v0": vols[0], "v1": vols[1]}})
    return pairs, excluded


def worker(pair: dict) -> dict:
    cache = CACHE_DIR / f"forecast_{pair['patient_id']}.json"
    if cache.exists():
        rec = json.loads(cache.read_text())
        if rec.get("atlas_sha") == ria._ATLAS_SHA and rec.get("dt_days") == pair["dt_days"] \
                and rec.get("tp1") == pair["tp1"]:
            return rec
    t0 = time.time()
    rec = ria.forecast_patient(pair, ria._ATLAS)
    rec["hist"] = pair["hist"]
    rec["seconds"] = round(time.time() - t0, 1)
    rec["atlas_sha"] = ria._ATLAS_SHA
    tmp = cache.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.replace(cache)
    return rec


def stage_forecast(limit: int | None) -> None:
    import multiprocessing as mp
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    pairs, excluded = stratified_pairs()
    if limit:
        pairs = pairs[:limit]
    recs = {}
    with mp.get_context("spawn").Pool(ria.FORECAST_WORKERS, initializer=ria._forecast_init) as pool:
        for k, rec in enumerate(pool.imap_unordered(worker, pairs)):
            recs[rec["patient_id"]] = rec
            best = max(rec["dice"]["grid"].values()) if "grid" in rec.get("dice", {}) else None
            print(f"[{k + 1}/{len(pairs)}] {rec['patient_id']} dt={rec['dt_days']:.0f}d "
                  f"no_change={rec.get('dice', {}).get('no_change')} best_grid={best} {rec['seconds']}s", flush=True)
    GRID_JSON.write_text(json.dumps({"rhos": ria.RHOS, "ds_mm2_per_day": ria.DS, "sharpen": ria.SHARPEN,
                                     "u_visible": ria.U_VISIBLE, "atlas_sha": ria._atlas_sha(),
                                     "excluded": excluded,
                                     "patients": [recs[p["patient_id"]] for p in pairs]}, indent=1))


def features(rec: dict) -> list:
    h, t = rec["hist"], rec.get("treatment") or {}
    d1 = rec["t1_day"]
    rs, re_ = f(t.get("radiation_start_day")), f(t.get("radiation_end_day"))
    ts, te = f(t.get("tmz_start_day")), f(t.get("tmz_end_day"))
    sg = f(t.get("surgery_day"))
    nan = float("nan")
    return [math.log(h["v1"] / h["v0"]) / (h["day1"] - h["day0"]), math.log(h["v1"]),
            h["day1"] - h["day0"], rec["dt_days"],
            d1 - re_ if re_ is not None else nan, d1 - sg if sg is not None else nan,
            float(rs is not None and re_ is not None), float(ts is not None and te is not None),
            float(ts is not None and te is not None and ts <= d1 <= te),
            float(rs is not None and re_ is not None and rs <= d1 <= re_)]


def oof_classifier(recs, y):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.preprocessing import StandardScaler
    X = np.array([features(r) for r in recs], float)
    fold = np.random.default_rng(51).permutation(len(recs)) % N_FOLDS   # same folds as cv_select
    prob = np.zeros(len(recs))
    for k in range(N_FOLDS):
        tr, te = fold != k, fold == k
        med = np.nanmedian(X[tr], axis=0)
        med = np.where(np.isnan(med), 0.0, med)
        Xtr, Xte = np.where(np.isnan(X[tr]), med, X[tr]), np.where(np.isnan(X[te]), med, X[te])
        sc = StandardScaler().fit(Xtr)
        m = LogisticRegression(C=1.0, max_iter=1000).fit(sc.transform(Xtr), y[tr])
        prob[te] = m.predict_proba(sc.transform(Xte))[:, 1]
    pred = prob >= PROB_THRESHOLD
    persistence = X[:, 0] > 0
    return prob, pred, {
        "oof_accuracy": float((pred == y).mean()), "oof_auc": float(roc_auc_score(y, prob)),
        "always_grow_accuracy": float(y.mean()), "persistence_accuracy": float((persistence == y).mean()),
        "n_predicted_grow": int(pred.sum()), "n_predicted_not_grow": int((~pred).sum()),
        "ppv_predicted_grow": float(y[pred].mean()) if pred.any() else None,
        "feature_names": ["log_vol_rate_scan0to1", "log_v1", "interval_0to1", "interval_1to2",
                          "days_since_rt_end", "days_since_surgery", "rt_recorded", "tmz_recorded",
                          "scan1_during_tmz", "scan1_during_rt"]}


def analyze() -> None:
    g = json.loads(GRID_JSON.read_text())
    recs = [r for r in g["patients"] if "skip" not in r]
    skipped = [(r["patient_id"], r["skip"]) for r in g["patients"] if "skip" in r]
    y = np.array([r["grew"] for r in recs])
    prob, pred, clf = oof_classifier(recs, y)

    oof, chosen = ria.cv_select(recs, "dice")
    oof["no_change"] = np.array([r["dice"]["no_change"] for r in recs])
    arms = ["anisotropic", "iso_same", "iso_homog", "no_change"]

    def summ(x):
        return {"mean": float(np.mean(x)), "mean_ci95": ria.boot_ci(x), "median": float(np.median(x)),
                "n": int(len(x))}

    def block(mask):
        return {a: summ(oof[a][mask]) for a in arms} if mask.sum() > 0 else None

    def tests(mask):
        if mask.sum() <= 5:
            return None
        return {"anisotropic_vs_no_change": ria.paired(oof["anisotropic"][mask], oof["no_change"][mask]),
                "anisotropic_vs_iso_same": ria.paired(oof["anisotropic"][mask], oof["iso_same"][mask])}

    all_mask = np.ones(len(recs), bool)
    pg = tests(pred)
    adj = ria.holm({k: v["wilcoxon_p"] for k, v in pg.items()}) if pg else {}
    for k in adj:
        pg[k]["wilcoxon_p_holm"] = adj[k]

    sec = None
    idx = np.where(pred)[0]
    if len(idx) > 5:
        sub = [recs[i] for i in idx]
        vm, vm_chosen = ria.cv_select(sub, "dice_volume_matched")
        vm["uniform_dilation"] = np.array([r["dice_volume_matched"]["uniform_dilation"] for r in sub])
        sec = {"n": len(sub), "note": "uses observed scan-2 volume: shape test, not a forecast",
               "arms": {a: summ(v) for a, v in vm.items()},
               "anisotropic_vs_iso_same": ria.paired(vm["anisotropic"], vm["iso_same"]),
               "anisotropic_vs_uniform_dilation": ria.paired(vm["anisotropic"], vm["uniform_dilation"])}

    res = {
        "design": "scan index 1 -> 2 forecast; pre-specified growth classifier on scan0->1 history, "
                  "treatment timing and interval; run_improved_aniso.forecast_patient unchanged",
        "n_scored": len(recs), "n_grew": int(y.sum()), "n_shrank_or_same": int((~y).sum()),
        "excluded": g["excluded"], "skipped": skipped,
        "interval_days_median": float(np.median([r["dt_days"] for r in recs])),
        "classifier": clf,
        "intention_to_forecast_all_patients": block(all_mask),
        "all_patients_tests": tests(all_mask),
        "stratum_predicted_to_grow": {"dice": block(pred), "tests_primary": pg},
        "stratum_predicted_not_to_grow": {"dice": block(~pred), "tests": tests(~pred)},
        "descriptive_outcome_selected_grew": {"dice": block(y), "tests": tests(y),
                                              "note": "selected on the scan-2 outcome; not a defensible headline"},
        "descriptive_outcome_selected_shrank_or_same": {"dice": block(~y), "tests": tests(~y)},
        "secondary_shape_only_predicted_to_grow": sec,
        "chosen_parameters_per_fold": chosen,
        "per_patient": [{"patient_id": r["patient_id"], "grew": bool(y[i]), "prob_grow_oof": float(prob[i]),
                         "predicted_grow": bool(pred[i]), **{f"dice_{a}": float(oof[a][i]) for a in arms}}
                        for i, r in enumerate(recs)]}
    RESULTS_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_patient", "chosen_parameters_per_fold")},
                     indent=2))
    print(f"[saved] {RESULTS_JSON}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["forecast", "analyze"])
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    stage_forecast(a.limit) if a.stage == "forecast" else analyze()


if __name__ == "__main__":
    main()
