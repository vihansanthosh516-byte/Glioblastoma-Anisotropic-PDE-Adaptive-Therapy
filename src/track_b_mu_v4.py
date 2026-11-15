"""Track B v4: median-loss boosted latent model. Exploratory, try 4.

Amendment: paper/trackB_mu_v4_amendment_2026-11-14.md (locked before scoring).
Baseline predictions are read unchanged from the v2 run output.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score

import track_b_mu_tcks as t

ROOT = Path(__file__).resolve().parent.parent
V2 = ROOT / "output" / "track_b_mu_tcks_v2_20261114" / "predictions.csv"
OUT = ROOT / "output" / "track_b_mu_v4_20261114"
TRY = 4
BASE = ["persistence", "two_point_exp", "pop_rate", "eb_random_slope", "gompertz", "gbr"]


def latent(df, tr):
    r_pop = float(tr.r_obs.median())
    r_ph = tr.groupby("phase").r_obs.median().reindex(t.PHASES).fillna(r_pop)
    c_tr = tr.phase.map(r_ph).values
    t2, s2 = t.eb_params(tr, c_tr)
    c = df.phase.map(r_ph).values
    has = df.r_hist.notna().values
    dtp = df.dt_prior.fillna(1.0).values
    x_prior = df.y_in.values - df.r_hist.fillna(0.0).values * dtp + c * dtp
    P = s2 + t2 * dtp ** 2
    return np.where(has, x_prior + P / (P + s2) * (df.y_in.values - x_prior), df.y_in.values)


def design(df, xh):
    X = pd.DataFrame({"x_hat": xh, "y_in": df.y_in.values, "dt": df.dt_days.values,
                      "has_prior": df.r_hist.notna().astype(float).values})
    for ph in t.PHASES[:-1]:
        X["dt_" + ph] = df.dt_days.values * (df.phase.values == ph)
        X["is_" + ph] = (df.phase.values == ph).astype(float)
    return X


def boot_ci(diff, rng, level):
    b = np.array([diff[rng.integers(0, len(diff), len(diff))].mean() for _ in range(5000)])
    a = (1 - level) / 2 * 100
    return [float(np.percentile(b, a)), float(np.percentile(b, 100 - a))]


def run():
    d, _ = t.load()
    v2 = pd.read_csv(V2)
    assert len(v2) == len(d) and (v2.patient_id.values == d.patient_id.values).all()
    pred = np.full(len(d), np.nan)
    for f in sorted(d.fold.unique()):
        trm, tem = (d.fold != f).values, (d.fold == f).values
        tr, te = d[trm], d[tem]
        q = HistGradientBoostingRegressor(loss="absolute_error", random_state=0)
        q.fit(design(tr, latent(tr, tr)), tr.y_out.values)
        pred[tem] = q.predict(design(te, latent(te, tr)))
    y = d.y_out.values
    errs = {k: np.abs(v2[k].values - y) for k in BASE}
    errs["v4"] = np.abs(pred - y)
    pat = pd.DataFrame(errs).assign(pid=d.patient_id.values).groupby("pid").mean()
    mae = pat.mean().to_dict()
    best = min(BASE, key=mae.get)
    level = 1 - 0.05 / TRY
    rng = np.random.default_rng(t.SEED)
    diff = (pat.v4 - pat[best]).values
    ci = boot_ci(diff, rng, level)
    grow = (d.core_vox_out > 1.25 * d.core_vox_in).astype(int).values
    res = {
        "script": "src/track_b_mu_v4.py", "amendment": "paper/trackB_mu_v4_amendment_2026-11-14.md",
        "label": "exploratory, post hoc, try 4 on reused MU-Glioma-Post",
        "n_pairs": int(len(d)), "n_patients": int(len(pat)),
        "patient_mae_log_volume": mae, "best_baseline": best,
        "v4_minus_best": {"mean": float(diff.mean()), "ci_level": level, "ci": ci},
        "success": bool(ci[1] < 0),
        "secondary_auc_growth25_v4": float(roc_auc_score(grow, pred - d.y_in.values)),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "result.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    pd.DataFrame({"patient_id": d.patient_id, "fold": d.fold, "y_out": y, "v4": pred}).to_csv(OUT / "predictions.csv", index=False)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    run()
