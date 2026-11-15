"""Track B: visit-informed growth hybrid (VIG), one locked run on LUMIERE.

Protocol: paper/trackB_vig_lumiere_protocol_2026-11-14.md (locked before scoring).
Fit on all development pairs (MU + UCSD + UCSF). Apply once to LUMIERE.
`--inventory` prints input-side counts only (no outcome is summarised).
"""
from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.ensemble import HistGradientBoostingRegressor

import track_b_mu_tcks as t
import track_b_v8_ucsd_external as ext
import track_b_v8_ucsf_alptdg as uf

ROOT = Path(__file__).resolve().parent.parent
LUM = ROOT / "data" / "external" / "lumiere"
SEG_CSV = LUM / "LUMIERE-pyradiomics-deepbratumia-features.csv"
RATE_CSV = LUM / "LUMIERE-ExpertRating-v202211.csv"
DEMO_CSV = LUM / "LUMIERE-Demographics_Pathology.csv"
OUT = ROOT / "output" / "track_b_vig_lumiere_20261114"
SEED = 20261114
MIN_PATIENTS = 25
NBOOT = 50000
COLS = ["patient_id", "day_in", "dt_days", "y_in", "y_out", "core_vox_in", "core_vox_out",
        "r_hist", "dt_prior", "phase", "age", "mgmt", "gbm", "r_obs"]
STARTS = [[np.log(0.01), 0.0, np.log(60)], [np.log(0.02), -0.1, np.log(30)], [np.log(0.005), 0.0, np.log(120)]]


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def dev_pairs():
    mu, _ = t.load()
    u = ext.ucsd_pairs()
    f, _ = uf.build(zipfile.ZipFile(uf.ZIP))
    out = []
    for name, d in (("MU", mu), ("UCSD", u), ("UCSF", f)):
        d = d[COLS].copy()
        d["patient_id"] = name + "_" + d.patient_id.astype(str)
        d["cohort"] = name
        out.append(d)
    return pd.concat(out, ignore_index=True)


def lumiere_scans():
    use = ["Patient", "Time point", "Label name", "Sequence", "diagnostics_Mask-original_Spacing", "diagnostics_Mask-original_VoxelNum"]
    d = pd.read_csv(SEG_CSV, usecols=use)
    d = d[d.Sequence == "CT1"]
    sp = d["diagnostics_Mask-original_Spacing"].dropna().unique()
    assert set(sp) <= {"(1.0, 1.0, 1.0)"}, f"unexpected voxel spacing {sp}"
    v = d.set_index(["Patient", "Time point", "Label name"])["diagnostics_Mask-original_VoxelNum"].unstack("Label name").fillna(0.0)
    s = v.reset_index()
    s["core"] = s["Necrosis"] + s["Contrast-enhancing"]
    parts = s["Time point"].str.split("-")
    s["week"] = parts.str[1].astype(int)
    s["sub"] = parts.str[2].fillna("0").astype(int)
    r = pd.read_csv(RATE_CSV)
    r.columns = ["Patient", "Time point", "lt3", "nm", "rating", "why"]
    pre = r[r.rating.astype(str).str.strip() == "Pre-Op"]
    keys = set(zip(pre["Patient"], pre["Time point"]))          # a scan is pre-operative if any rating row says so
    s["pre_op"] = [(a, b) in keys for a, b in zip(s["Patient"], s["Time point"])]
    return s


def lumiere_pairs():
    s = lumiere_scans()
    n_scans, n_pre = len(s), int(s.pre_op.sum())
    s = s[~s.pre_op].sort_values(["Patient", "week", "sub"])
    age = pd.read_csv(DEMO_CSV).set_index("Patient")["Age at surgery (years)"]
    rows = []
    for pid, g in s.groupby("Patient"):
        g = g.reset_index(drop=True)
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            r_hist = dt_prior = np.nan
            if i > 0:
                pv = g.iloc[i - 1]
                if pv.core > 0 and a.week > pv.week:
                    dt_prior = 7.0 * (a.week - pv.week)
                    r_hist = (np.log(a.core + 1) - np.log(pv.core + 1)) / dt_prior
            rows.append({"patient_id": "LUM_" + pid, "day_in": 7.0 * a.week, "dt_days": 7.0 * (b.week - a.week),
                         "core_vox_in": a.core, "core_vox_out": b.core, "r_hist": r_hist, "dt_prior": dt_prior,
                         "phase": "other", "age": pd.to_numeric(age.get(pid), errors="coerce"), "mgmt": np.nan, "gbm": 1.0})
    d = pd.DataFrame(rows)
    n_all = len(d)
    d = d[(d.dt_days >= 14) & (d.dt_days <= 365) & (d.core_vox_in > 0)].reset_index(drop=True)
    d["y_in"], d["y_out"] = np.log(d.core_vox_in + 1.0), np.log(d.core_vox_out + 1.0)
    d["r_obs"] = (d.y_out - d.y_in) / d.dt_days
    d["cohort"] = "LUMIERE"
    return d, {"n_scans_with_mask": n_scans, "n_pre_op_excluded": n_pre, "n_consecutive_pairs": n_all, "n_eligible_pairs": int(len(d))}


# ---- VIG: growth law x visit weight + history smoothing -------------------------------------
def vig_curve(th, dt):
    rho, dq, tau = np.exp(th[0]), th[1], np.exp(th[2])
    w = np.exp(-dt / tau)                       # weight of the "growth suspected" regime given the scan interval
    return w * rho * dt + (1 - w) * dq          # exponential growth in that regime, constant drift in the other


def last_change(d):
    return (d.r_hist * d.dt_prior).fillna(0.0).values


def fit_vig(tr):
    ch, dt = (tr.y_out - tr.y_in).values, tr.dt_days.values
    loss = lambda th: float(np.mean(np.abs(ch - vig_curve(th, dt))))
    fits = [minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-6, "fatol": 1e-9, "maxiter": 10000}) for x0 in STARTS]
    th = min(fits, key=lambda r: r.fun).x
    has = tr.r_hist.notna().values
    res, lc = (ch - vig_curve(th, dt))[has], last_change(tr)[has]
    X = np.c_[np.maximum(lc, 0.0), np.minimum(lc, 0.0)]      # pull-back after a rise, pull-back after a fall
    grid, lam = np.linspace(0.0, 0.6, 31), np.zeros(2)
    for _ in range(4):                                        # coordinate grid search, absolute loss
        for j in range(2):
            cand = []
            for g in grid:
                b = lam.copy(); b[j] = g
                cand.append(np.mean(np.abs(res + X @ b)))
            lam[j] = grid[int(np.argmin(cand))]
    return th, lam


def gompertz_l1(tr):
    f = lambda th, y, dt: th[1] - (th[1] - y) * np.exp(-th[0] * dt)
    loss = lambda th: float(np.mean(np.abs(f(th, tr.y_in.values, tr.dt_days.values) - tr.y_out.values)))
    x = min((minimize(loss, x0, method="Nelder-Mead", options={"xatol": 1e-7, "fatol": 1e-9, "maxiter": 10000})
             for x0 in ([1e-3, tr.y_in.max() + 0.5], [1e-4, 9.0], [5e-3, 12.0])), key=lambda r: r.fun).x
    return lambda te: f(x, te.y_in.values, te.dt_days.values)


def fit_predict(tr, te):
    out = {}
    r_pop = float(tr.r_obs.median())
    out["persistence"] = te.y_in.values
    out["two_point_exp"] = (te.y_in + te.r_hist.fillna(0.0) * te.dt_days).values
    out["pop_rate"] = (te.y_in + r_pop * te.dt_days).values
    t2, s2 = t.eb_params(tr, r_pop)
    out["eb_random_slope"] = np.asarray(t.shrink_pred(te, r_pop, t2, s2)[0])
    out["gompertz"] = t.gompertz_fit(tr)(te)
    out["gbr"] = HistGradientBoostingRegressor(random_state=0).fit(t.features(tr), tr.y_out).predict(t.features(te))
    out["gompertz_l1"] = gompertz_l1(tr)(te)
    out["gbr_l1"] = HistGradientBoostingRegressor(loss="absolute_error", random_state=0).fit(t.features(tr), tr.y_out).predict(t.features(te))
    th, lam = fit_vig(tr)
    out["vig_nohist"] = te.y_in.values + vig_curve(th, te.dt_days.values)
    lc = last_change(te)
    out["vig"] = out["vig_nohist"] - (lam[0] * np.maximum(lc, 0.0) + lam[1] * np.minimum(lc, 0.0))
    par = {"rho_per_day": float(np.exp(th[0])), "drift": float(th[1]), "tau_days": float(np.exp(th[2])), "lambda_after_rise": float(lam[0]), "lambda_after_fall": float(lam[1]), "pop_rate_per_day": r_pop,
           "curve_at_days": {str(k): float(vig_curve(th, np.array([float(k)]))[0]) for k in (20, 30, 45, 60, 90, 120, 180, 300)}}
    return out, par


def paired(diff, rng):
    idx = rng.integers(0, len(diff), (NBOOT, len(diff)))
    b = diff[idx].mean(1)
    flips = rng.choice([-1.0, 1.0], (NBOOT, len(diff)))
    perm = np.abs((flips * diff).mean(1))
    return {"mean": float(diff.mean()), "ci95": [float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5))],
            "sign_flip_p": float((np.sum(perm >= abs(diff.mean())) + 1) / (NBOOT + 1)), "sd": float(diff.std(ddof=1))}


def run():
    OUT.mkdir(parents=True, exist_ok=True)
    if (OUT / "result.json").exists():
        raise SystemExit("locked test already run; refusing to overwrite")
    dev = dev_pairs()
    lum, inv = lumiere_pairs()
    n_pat = int(lum.patient_id.nunique())
    res = {"script": "src/track_b_vig_lumiere.py", "protocol": "paper/trackB_vig_lumiere_protocol_2026-11-14.md",
           "label": "one locked run; LUMIERE was exposed before for other targets, so it is not an untouched cohort",
           "dev": {c: {"pairs": int(len(g)), "patients": int(g.patient_id.nunique())} for c, g in dev.groupby("cohort")},
           "lumiere_inventory": inv, "n_pairs": int(len(lum)), "n_patients": n_pat,
           "input_sha256": {p.name: sha(p) for p in (SEG_CSV, RATE_CSV, DEMO_CSV)}}
    if n_pat < MIN_PATIENTS:
        res["status"] = "no-go: fewer than 25 patients"
    else:
        preds, par = fit_predict(dev, lum)
        y = lum.y_out.values
        pat = pd.DataFrame({k: np.abs(v - y) for k, v in preds.items()}).assign(pid=lum.patient_id.values).groupby("pid").mean()
        rng = np.random.default_rng(SEED)
        comps = ["persistence", "two_point_exp", "pop_rate", "eb_random_slope", "gompertz", "gbr", "gompertz_l1", "gbr_l1"]
        vs = {b: paired((pat.vig - pat[b]).values, rng) for b in comps}
        vs_nohist = {b: paired((pat.vig_nohist - pat[b]).values, rng) for b in comps}
        e = pd.DataFrame({k: np.abs(v - y) for k, v in preds.items()})
        dtb = pd.cut(lum.dt_days, [0, 45, 75, 105, 400]).astype(str)
        res.update({"status": "scored", "vig_parameters": par, "patient_mae_log_volume": pat.mean().to_dict(),
                    "vig_minus": vs, "success": bool(all(v["ci95"][1] < 0 for v in vs.values())),
                    "secondary": {"vig_nohist_minus": vs_nohist,
                                  "pair_mae_by_interval": {k: {"n": int(len(g)), "persistence": float(g.persistence.mean()), "pop_rate": float(g.pop_rate.mean()), "vig": float(g.vig.mean())} for k, g in e.groupby(dtb.values)},
                                  "n_pairs_with_history": int(lum.r_hist.notna().sum()), "n_zero_outcome_pairs": int((lum.core_vox_out == 0).sum())}})
        pd.DataFrame(preds).assign(patient_id=lum.patient_id, dt_days=lum.dt_days, y_out=y).to_csv(OUT / "predictions.csv", index=False)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(json.dumps(res, indent=1))


def inventory():
    lum, inv = lumiere_pairs()
    print(json.dumps({**inv, "n_patients": int(lum.patient_id.nunique()), "dt_days_quantiles": lum.dt_days.quantile([0, .25, .5, .75, 1]).to_dict(),
                      "share_pairs_with_history": float(lum.r_hist.notna().mean()), "y_in_quantiles": lum.y_in.quantile([0, .25, .5, .75, 1]).round(2).to_dict(),
                      "pairs_per_patient_median": float(lum.groupby("patient_id").size().median())}, indent=1))


if __name__ == "__main__":
    inventory() if "--inventory" in sys.argv else run()
