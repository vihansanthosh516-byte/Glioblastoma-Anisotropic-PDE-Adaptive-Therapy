"""Track B TCKS vs baselines on MU-Glioma-Post core volumes. Exploratory.

Protocol: paper/trackB_mu_tcks_protocol_2026-11-14.md (locked before scoring).
Reads manifest volumes and clinical dates only; no image is read.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
PAIRS = ROOT / "data" / "manifests" / "forecast_pairs_mu.csv"
SPLIT = ROOT / "data" / "manifests" / "split_mu.csv"
CLIN = ROOT / "data" / "tcia" / "MU-Glioma-Post_ClinicalData-July2025.xlsx"
OUT = ROOT / "output" / "track_b_mu_tcks_20261114"
PHASES = ["post_progression", "chemoradiation", "adjuvant", "other"]
SEED = 20261114


def num(s):
    return pd.to_numeric(s, errors="coerce")


def load():
    p = pd.read_csv(PAIRS)
    p = p.merge(pd.read_csv(SPLIT), on="patient_id", how="left")
    p["y_in"] = np.log(p.core_vox_in + 1.0)
    p["y_out"] = np.log(p.core_vox_out + 1.0)
    p = p.sort_values(["patient_id", "day_in"]).reset_index(drop=True)
    # prior-scan rate from the previous pair of the same patient (any window)
    p["r_hist"], p["dt_prior"] = np.nan, np.nan
    for _, g in p.groupby("patient_id"):
        idx = list(g.index)
        for a, b in zip(idx[:-1], idx[1:]):
            if p.at[a, "tp_out"] == p.at[b, "tp_in"] and p.at[a, "core_vox_in"] > 0:
                dtp = p.at[b, "day_in"] - p.at[a, "day_in"]
                if dtp > 0:
                    p.at[b, "r_hist"] = (p.at[b, "y_in"] - p.at[a, "y_in"]) / dtp
                    p.at[b, "dt_prior"] = dtp
    c = pd.read_excel(CLIN, sheet_name="MU Glioma Post")
    c.columns = [str(x).strip() for x in c.columns]
    col = lambda key: next(x for x in c.columns if key.lower() in x.lower())
    k = pd.DataFrame({
        "patient_id": c["Patient_ID"].astype(str).str.strip(),
        "rt_s": num(c[col("Radiation Therapy Start")]),
        "rt_e": num(c[col("Radiation Therapy end")]),
        "ad_s": num(c[col("Diagnosis to Starting Additional Therapy")]),
        "ad_e": num(c[col("Diagnosis to Complete Additional Therapy")]),
        "prog": num(c[col("Diagnosis to date of First Progression")]),
        "age": num(c["Age at diagnosis"]),
        "mgmt": num(c["MGMT methylation"]),
        "gbm": (c["Primary Diagnosis"].astype(str).str.strip() == "GBM").astype(float),
    })
    p = p.merge(k, on="patient_id", how="left")
    m = (p.day_in + p.day_out) / 2.0
    ph = np.full(len(p), "other", dtype=object)
    ph[(p.ad_s <= m) & (m <= p.ad_e)] = "adjuvant"
    ph[(p.rt_s <= m) & (m <= p.rt_e)] = "chemoradiation"
    ph[p.prog <= p.day_in] = "post_progression"
    p["phase"] = ph
    p["r_obs"] = (p.y_out - p.y_in) / p.dt_days
    prim = p[p.in_primary & p.input_core_nonempty].reset_index(drop=True)
    return prim, p


def eb_params(tr, center):
    """Method-of-moments tau2 (true rate spread) and s2 (log-volume noise)."""
    dev = tr.r_obs - center
    short = tr[tr.dt_days < 60]
    s2 = float(np.var(short.y_out - short.y_in) / 2.0) if len(short) > 5 else float(np.var(tr.y_out - tr.y_in) / 2)
    noise_rate = np.mean(2 * s2 / tr.dt_days ** 2)
    tau2 = max(float(np.var(dev)) - noise_rate, 1e-12)
    return tau2, s2


def shrink_pred(te, center_te, tau2, s2):
    w = np.where(te.r_hist.notna(), tau2 / (tau2 + s2 / np.square(te.dt_prior.fillna(1.0))), 0.0)
    r = w * te.r_hist.fillna(0.0) + (1 - w) * center_te
    return te.y_in + r * te.dt_days, w


def gompertz_fit(tr):
    def f(th, y, dt):
        a, ymax = th
        return ymax - (ymax - y) * np.exp(-a * dt)
    res = least_squares(lambda th: f(th, tr.y_in.values, tr.dt_days.values) - tr.y_out.values,
                        x0=[1e-3, tr.y_in.max() + 0.5], bounds=([-0.1, 0], [0.1, 20]))
    return lambda te: f(res.x, te.y_in.values, te.dt_days.values)


def features(d):
    X = pd.DataFrame({"y_in": d.y_in, "dt": d.dt_days, "r_hist": d.r_hist,
                      "has_prior": d.r_hist.notna().astype(float), "age": d.age,
                      "mgmt": d.mgmt, "gbm": d.gbm})
    for ph in PHASES:
        X["ph_" + ph] = (d.phase == ph).astype(float)
    return X


def run():
    d, _ = load()
    preds = {k: np.full(len(d), np.nan) for k in
             ["persistence", "two_point_exp", "pop_rate", "eb_random_slope", "gompertz", "gbr", "tcks"]}
    wts = np.full(len(d), np.nan)
    for f in sorted(d.fold.unique()):
        tr, te_m = d[d.fold != f], d.fold == f
        te = d[te_m]
        r_pop = float(tr.r_obs.median())
        preds["persistence"][te_m] = te.y_in
        preds["two_point_exp"][te_m] = te.y_in + te.r_hist.fillna(0.0) * te.dt_days
        preds["pop_rate"][te_m] = te.y_in + r_pop * te.dt_days
        t2, s2 = eb_params(tr, r_pop)
        preds["eb_random_slope"][te_m] = shrink_pred(te, r_pop, t2, s2)[0]
        preds["gompertz"][te_m] = gompertz_fit(tr)(te)
        g = HistGradientBoostingRegressor(random_state=0).fit(features(tr), tr.y_out)
        preds["gbr"][te_m] = g.predict(features(te))
        r_ph = tr.groupby("phase").r_obs.median().reindex(PHASES).fillna(r_pop)
        c_tr, c_te = tr.phase.map(r_ph).values, te.phase.map(r_ph).values
        t2p, s2p = eb_params(tr, c_tr)
        yh, w = shrink_pred(te, c_te, t2p, s2p)
        preds["tcks"][te_m], wts[te_m] = yh, w
    err = {k: np.abs(v - d.y_out.values) for k, v in preds.items()}
    pat = pd.DataFrame({k: v for k, v in err.items()}).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    base = [k for k in mae if k != "tcks"]
    best = min(base, key=mae.get)
    diff = (pat.tcks - pat[best]).values
    rng = np.random.default_rng(SEED)
    boots = np.array([diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(5000)])
    ci = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
    prior = d.r_hist.notna().values
    pat_p = pd.DataFrame({k: v[prior] for k, v in err.items()}).assign(pid=d.patient_id.values[prior]).groupby("pid").mean()
    grow = (d.core_vox_out > 1.25 * d.core_vox_in).astype(int).values
    auc = {k: float(roc_auc_score(grow, v - d.y_in.values)) for k, v in preds.items() if k != "persistence"}
    vs_all = {}
    for b in base:
        dd = (pat.tcks - pat[b]).values
        bb = np.array([dd[rng.integers(0, len(dd), len(dd))].mean() for _ in range(5000)])
        vs_all[b] = {"mean_diff": float(dd.mean()), "ci95": [float(np.percentile(bb, 2.5)), float(np.percentile(bb, 97.5))]}
    res = {
        "script": "src/track_b_mu_tcks.py", "protocol": "paper/trackB_mu_tcks_protocol_2026-11-14.md",
        "label": "exploratory (reused MU-Glioma-Post cohort)",
        "n_pairs": int(len(d)), "n_patients": int(d.patient_id.nunique()),
        "n_pairs_with_prior": int(prior.sum()),
        "phase_counts": d.phase.value_counts().to_dict(),
        "patient_mae_log_volume": mae, "best_baseline": best,
        "tcks_minus_best": {"mean": float(diff.mean()), "ci95": ci},
        "success": bool(ci[1] < 0), "tcks_vs_each_baseline": vs_all,
        "secondary_mae_pairs_with_prior": pat_p.mean().to_dict(),
        "secondary_auc_growth25": auc,
        "tcks_mean_weight_with_prior": float(np.nanmean(wts[prior])),
        "fisher_kpp_pde": "not rerun; prior MU mask result did not beat persistence (output/mu_fitted_response/results.json)",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame(preds).assign(patient_id=d.patient_id, y_out=d.y_out, phase=d.phase, fold=d.fold).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
