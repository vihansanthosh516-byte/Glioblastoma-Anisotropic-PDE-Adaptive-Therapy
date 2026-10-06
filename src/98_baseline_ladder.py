#!/usr/bin/env python3
"""Script 98: baseline ladder on MU-Glioma-Post core masks (Phase 2.2; analysis_plan_v1.md A1.4).

Forecast the next scan's core mask (labels 1 and 3, 2 mm grid, as script 81) with no PDE:
  persistence            S(t)
  geometric_train_rate   S(t) grown/shrunk (signed distance inside brain) to V_t exp(g dt), g = median log-growth rate
                         of the TRAINING-fold patients (manifest folds), so no information from the test patient
  last_rate              V_t exp(g_last dt), g_last from the previous pair of scans (persistence if none)
  linear                 V_t + slope dt, slope from the previous pair, floored at 0 (persistence if none)
  gompertz               K fixed at the training-fold 95th percentile volume, V0 = first scan, rate fitted on scans up to t
                         (needs >= 3 scans; persistence otherwise)
  oracle_volume_matched  S(t) reshaped to the TRUE next-scan volume. NOT a baseline: an upper bound on volume-only methods.
Exponential extrapolation is last_rate (reported once). Every arm shares the crop, 2 mm grid, brain domain and scorer.

Input : data/manifests/{forecast_pairs_mu,patient_manifest,split_mu}.csv, MU masks, MU T1n (brain mask)
Output: output/baseline_ladder.json, output/baseline_ladder_pairs.csv, output/baseline_ladder.manifest.json
Primary population per analysis_plan_v1.md A1.2 and A2.1 (all eligible; GBM-only reported next to it).
Delta = metric(method) - metric(persistence), aggregated to ONE value per patient, patient-level bootstrap.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage, optimize

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
MAN = PROJECT_ROOT / "data" / "manifests"
MU_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
sys.path.insert(0, str(PROJECT_ROOT))

import run_improved_aniso as ria  # noqa: E402
from src import forecast_metrics as fm  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

CORE = (1, 3)
SPACING = 2.0
METHODS = ["persistence", "geometric_train_rate", "last_rate", "linear", "gompertz", "oracle_volume_matched"]
BINS = [("short_lt60", 0, 60), ("medium_60_120", 60, 120.0001), ("long_gt120", 120.0001, 1e9)]
GRID2 = (120, 120, 77)


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def core_mask(pid, tp):
    p = MU_DIR / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_tumorMask.nii.gz"
    return ria.to_2mm(np.isin(np.asarray(nib.load(str(p)).dataobj), CORE)) >= 0.5


def brain_mask(pid, tp):
    p = MU_DIR / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_brain_t1n.nii.gz"
    return ria.to_2mm(np.asarray(nib.load(str(p)).dataobj) > 0) >= 0.5


def signed_domain_distance(S, brain):
    """Signed distance to the boundary of S (negative inside), +inf outside brain | S. Computed once per pair."""
    sd = ndimage.distance_transform_edt(~S) - ndimage.distance_transform_edt(S)
    # tie-break: voxels at equal distance are ordered by distance to the centroid of S (1e-4 mm per mm), so the
    # volume is matched exactly instead of overshooting by a whole distance shell
    c = np.argwhere(S).mean(0) if S.any() else np.zeros(3)
    grids = np.meshgrid(*[np.arange(n) for n in S.shape], indexing="ij", sparse=True)
    dc = np.sqrt(sum((g - ci) ** 2 for g, ci in zip(grids, c)))
    return np.where(brain | S, sd + 1e-4 * dc, np.inf)


def reshape_to_volume(S, sd, n_vox):
    """Mask with n_vox voxels: S grown (or shrunk) along the signed distance inside the brain domain."""
    n = int(round(max(n_vox, 0)))
    if n <= 0 or not S.any():
        return np.zeros_like(S)
    n = min(n, int(np.isfinite(sd).sum()))
    thr = np.partition(sd.ravel(), n - 1)[n - 1]
    return sd <= thr


def gompertz_forecast(days, vols, K, day_out):
    """Predicted volume at day_out from scans (days, vols), K fixed, V0 = first scan; None if not fittable."""
    if len(days) < 3 or vols[0] <= 0:
        return None
    t = np.asarray(days, float) - days[0]
    lv = np.log(np.asarray(vols, float) + 1.0)
    lK, l0 = np.log(K + 1.0), lv[0]

    def sse(a):
        return float(((lK + (l0 - lK) * np.exp(-a * t) - lv) ** 2).sum())

    r = optimize.minimize_scalar(sse, bounds=(1e-4, 0.2), method="bounded")
    return float(np.exp(lK + (l0 - lK) * np.exp(-r.x * (day_out - days[0]))) - 1.0)


def main(limit=None):
    t0 = time.time()
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")
    pats = pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)
    pats = pats[pats["dataset"] == "MU"].set_index("patient_id")
    split = pd.read_csv(MAN / "split_mu.csv").set_index("patient_id")["fold"]
    pairs["fold"] = pairs["patient_id"].map(split)
    pairs["is_gbm"] = pairs["patient_id"].map(pats["is_gbm"])
    pairs["patient_eligible"] = pairs["patient_id"].map(pats["eligible"]).astype(bool)
    if limit:
        keep = pairs["patient_id"].drop_duplicates().head(limit)
        pairs = pairs[pairs["patient_id"].isin(keep)]

    # pass 1: 2 mm core masks (packed) and volumes for every scan used
    packed, vol, day = {}, {}, {}
    for pid, g in pairs.groupby("patient_id"):
        for tp, d in list(zip(g["tp_in"], g["day_in"])) + list(zip(g["tp_out"], g["day_out"])):
            if (pid, tp) not in packed:
                m = core_mask(pid, tp)
                packed[(pid, tp)] = np.packbits(m.ravel())
                vol[(pid, tp)] = int(m.sum())
                day[(pid, tp)] = d
    log(f"pass 1: {len(packed)} masks ({time.time() - t0:.0f}s)")

    def unpack(k):
        return np.unpackbits(packed[k], count=int(np.prod(GRID2))).astype(bool).reshape(GRID2)

    pairs["v_in"] = [vol[(p, a)] for p, a in zip(pairs["patient_id"], pairs["tp_in"])]
    pairs["v_out"] = [vol[(p, b)] for p, b in zip(pairs["patient_id"], pairs["tp_out"])]
    pairs["log_rate"] = np.log((pairs["v_out"] + 1.0) / (pairs["v_in"] + 1.0)) / pairs["dt_days"]

    # training-fold constants (only primary pairs of patients in other folds)
    prim = pairs[pairs["in_primary"] & (pairs["v_in"] > 0)]
    g_train = {int(f): float(prim.loc[prim["fold"] != f, "log_rate"].median()) for f in sorted(split.unique())}
    K_train = {}
    for f in (int(x) for x in sorted(split.unique())):
        vs = [vol[(p, tp)] for (p, tp) in vol if split.get(p, -1) != f and vol[(p, tp)] > 0]
        K_train[f] = float(np.percentile(vs, 95))
    log(f"g_train {g_train}  K_train {K_train}")

    rows = []
    for i, r in enumerate(pairs.itertuples(index=False)):
        if not r.in_primary and not (r.input_core_nonempty and r.dt_days > 0):
            continue
        pid, f = r.patient_id, int(r.fold)
        S = unpack((pid, r.tp_in))
        T = unpack((pid, r.tp_out))
        if not S.any():
            continue
        sd = signed_domain_distance(S, brain_mask(pid, r.tp_in))
        # prior scans of this patient strictly before the input scan
        prior = sorted([(day[(pid, tp)], tp) for (p, tp) in packed if p == pid and day[(pid, tp)] < r.day_in])
        hist_days = [d for d, _ in prior] + [r.day_in]
        hist_vols = [vol[(pid, tp)] for _, tp in prior] + [r.v_in]
        V = float(r.v_in)
        forecasts, flags = {"persistence": S}, {}
        forecasts["geometric_train_rate"] = reshape_to_volume(S, sd, V * np.exp(g_train[f] * r.dt_days))
        if len(prior) >= 1:
            d0, v0 = hist_days[-2], hist_vols[-2]
            gl = np.log((V + 1) / (v0 + 1)) / (r.day_in - d0)
            slope = (V - v0) / (r.day_in - d0)
            forecasts["last_rate"] = reshape_to_volume(S, sd, (V + 1) * np.exp(gl * r.dt_days) - 1)
            forecasts["linear"] = reshape_to_volume(S, sd, max(V + slope * r.dt_days, 0.0))
        else:
            forecasts["last_rate"] = S
            forecasts["linear"] = S
            flags["last_rate"] = flags["linear"] = True
        gv = gompertz_forecast(hist_days, hist_vols, K_train[f], r.day_out)
        if gv is None:
            forecasts["gompertz"] = S
            flags["gompertz"] = True
        else:
            forecasts["gompertz"] = reshape_to_volume(S, sd, gv)
        forecasts["oracle_volume_matched"] = reshape_to_volume(S, sd, float(T.sum()))
        for m, F in forecasts.items():
            met = fm.all_metrics(F, T, SPACING)
            rows.append({"patient_id": pid, "tp_in": r.tp_in, "tp_out": r.tp_out, "method": m, "fold": f,
                         "dt_days": r.dt_days, "in_primary": bool(r.in_primary), "has_prior_scan": bool(r.has_prior_scan),
                         "n_prior_scans": len(prior), "is_gbm": bool(r.is_gbm),
                         "fallback_to_persistence": bool(flags.get(m, False)),
                         "v_in": r.v_in, "v_out": r.v_out, "v_pred": int(F.sum()), **met})
        if i % 40 == 0:
            log(f"pair {i}/{len(pairs)} ({time.time() - t0:.0f}s)")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "baseline_ladder_pairs.csv", index=False)
    log(f"forecasts done ({time.time() - t0:.0f}s)")
    summarise(res, pairs, g_train, K_train)


def horizon_bin(dt):
    for name, lo, hi in BINS:
        if lo <= dt < hi:
            return name


def summarise(res: pd.DataFrame, pairs, g_train, K_train):
    res = res.copy()
    res["bin"] = res["dt_days"].map(horizon_bin)
    out = {"script": "98_baseline_ladder", "target": "core labels 1,3 at 2 mm", "g_train_per_day": g_train,
           "K_train_voxels": K_train, "n_pairs_primary": int(res[res["in_primary"] & (res["method"] == "persistence")].shape[0])}

    def deltas(df, metric, method):
        a = ps.patient_means(df[df["method"] == method], metric)
        b = ps.patient_means(df[df["method"] == "persistence"], metric)
        j = pd.concat([a, b], axis=1, keys=["m", "p"]).dropna()
        return (j["m"] - j["p"]).to_numpy()

    strata = {
        "all_eligible_primary_pairs": res[res["in_primary"]],
        "gbm_only_primary_pairs": res[res["in_primary"] & res["is_gbm"]],
        "first_pair_only_primary": res[res["in_primary"] & ~res["has_prior_scan"]],
        "with_prior_scan_primary": res[res["in_primary"] & res["has_prior_scan"]],
    }
    for b, _, _ in BINS:
        strata[f"horizon_{b}"] = res[res["in_primary"] & (res["bin"] == b)]
    out["strata"] = {}
    for sname, df in strata.items():
        block = {"n_patients": int(df["patient_id"].nunique()), "methods": {}}
        pers = ps.patient_means(df[df["method"] == "persistence"], "dice")
        block["persistence_mean_dice"] = float(pers.mean()) if len(pers) else None
        for m in METHODS[1:]:
            sub = df[df["method"] == m]
            if sub.empty:
                continue
            d = deltas(df, "dice", m)
            if len(d) < 5:
                continue
            block["methods"][m] = {
                "dice_delta_vs_persistence": ps.summarize_delta(d, n_boot=2000),
                "mean_dice": float(ps.patient_means(sub, "dice").mean()),
                "median_abs_log_vol_error": float(np.nanmedian(np.abs(sub["log_vol_ratio"]))),
                "hd95_mean_mm": float(np.nanmean(ps.patient_means(sub, "hd95_mm"))),
                "centroid_mean_mm": float(np.nanmean(ps.patient_means(sub, "centroid_mm"))),
                "fallback_share": float(sub["fallback_to_persistence"].mean())}
        out["strata"][sname] = block
    # consistency with script 81 (same masks, first pair, any interval): persistence Dice per patient
    cons = []
    for f in (OUT / "forecast_labels_core" / "cache").glob("forecast_*.json"):
        c = json.loads(f.read_text())
        if "dice" in c and "no_change" in c["dice"]:
            cons.append((c["patient_id"], c["dice"]["no_change"]))
    mine = res[(res["method"] == "persistence") & (~res["has_prior_scan"])].groupby("patient_id")["dice"].first()
    s81 = pd.Series(dict(cons))
    j = pd.concat([s81, mine], axis=1, keys=["s81", "s98"]).dropna()
    out["consistency_with_script81_persistence"] = {
        "n_patients_compared": int(len(j)), "max_abs_diff": float((j["s81"] - j["s98"]).abs().max()) if len(j) else None,
        "mean_s81": float(j["s81"].mean()) if len(j) else None, "mean_s98": float(j["s98"].mean()) if len(j) else None}
    (OUT / "baseline_ladder.json").write_text(json.dumps(out, indent=1, default=float))
    write_run_manifest("baseline_ladder_98", OUT / "baseline_ladder.manifest.json", script="src/98_baseline_ladder.py",
                       seed=ps.DEFAULT_SEED, config={"methods": METHODS, "bins": BINS, "spacing_mm": SPACING},
                       inputs=[MAN / "forecast_pairs_mu.csv", MAN / "patient_manifest.csv", MAN / "split_mu.csv"],
                       dataset="MU-Glioma-Post", patient_split="data/manifests/split_mu.csv",
                       primary_endpoint="Dice delta vs persistence (per patient)")
    log("summary written")


if __name__ == "__main__":
    lim = int(sys.argv[1]) if len(sys.argv) > 1 else None
    main(lim)
