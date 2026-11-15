"""Track B hybrid candidate 'stack': physics-informed L1 stacking with a monotone size prior.

Protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md
Families joined: mechanistic growth (Gompertz, phase-rate latent state) + statistical learners
(median regression, L1 gradient boosting) + logistic collapse gate (v8 style).
Fitting uses tr only. Meta-weight settings are chosen by inner grouped CV inside tr.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, QuantileRegressor
from sklearn.model_selection import GroupKFold

import track_b_mu_tcks as t
import track_b_mu_v3 as v3

INNER_K = 5
ALPHAS = [0.0, 0.001, 0.01]
K_BASE = 5  # gompertz, phase_rate_latent, pop_rate, median_regression, l1_gbm


def _X(df, comp, xh):
    D = v3.design(df, xh).values.astype(float)
    C = comp.apply(pd.to_numeric, errors="coerce").values.astype(float)
    return np.column_stack([D, C])


def _prep(Xtr, Xte):
    med = np.nanmedian(Xtr, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    Xtr = np.where(np.isnan(Xtr), med, Xtr)
    Xte = np.where(np.isnan(Xte), med, Xte)
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0)
    sd = np.where(sd == 0, 1.0, sd)
    return (Xtr - mu) / sd, (Xte - mu) / sd


def _base_preds(fit, ap, cf, ca):
    """Fit all base learners on fit (rows of tr), predict rows of ap. Returns (len(ap), K_BASE)."""
    r_pop = float(fit.r_obs.median())
    r_ph = fit.groupby("phase").r_obs.median().reindex(t.PHASES).fillna(r_pop)
    c_ap = ap.phase.map(r_ph).values.astype(float)
    dt_ap = ap.dt_days.values.astype(float)
    cols = []
    # 1. Gompertz growth forecast
    cols.append(np.asarray(t.gompertz_fit(fit)(ap), float))
    # 2. Phase-rate latent forecast: smoothed input size + phase median rate * dt
    xh_ap = v3.latent(ap, fit)
    cols.append(xh_ap + c_ap * dt_ap)
    # 3. Population-rate forecast
    cols.append(ap.y_in.values.astype(float) + r_pop * dt_ap)
    # 4 and 5 need design + composition features
    xh_fit = v3.latent(fit, fit)
    Xf = _X(fit, cf, xh_fit)
    Xa = _X(ap, ca, xh_ap)
    Zf, Za = _prep(Xf, Xa)
    yf = fit.y_out.values.astype(float)
    collapse = (fit.core_vox_out < 0.05 * fit.core_vox_in).values
    # 4. v8-style median regression on non-collapse rows
    q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Zf[~collapse], yf[~collapse])
    cols.append(np.asarray(q.predict(Za), float))
    # 5. L1 gradient-boosted model on design + composition
    g = HistGradientBoostingRegressor(loss="absolute_error", random_state=0).fit(Xf, yf)
    cols.append(np.asarray(g.predict(Xa), float))
    return np.column_stack(cols)


def _fit_meta(P, y, alpha):
    """L1 loss, weights >= 0 summing to 1, free intercept.
    alpha penalises L1 distance of weights from uniform (plain L1 of w is constant under sum=1)."""
    K = P.shape[1]
    u = np.full(K, 1.0 / K)

    def obj(x):
        w, b = x[:K], x[K]
        return np.mean(np.abs(y - (b + P @ w))) + alpha * np.sum(np.abs(w - u))

    x0 = np.concatenate([u, [float(np.median(y - P @ u))]])
    cons = [{"type": "eq", "fun": lambda x: np.sum(x[:K]) - 1.0}]
    bounds = [(0.0, 1.0)] * K + [(None, None)]
    res = minimize(obj, x0, method="SLSQP", bounds=bounds, constraints=cons,
                   options={"maxiter": 1000, "ftol": 1e-10})
    w = np.clip(res.x[:K], 0.0, None)
    w = w / w.sum()
    return w, float(res.x[K])


def _choose_alpha(P, y, groups):
    best_a, best_m = ALPHAS[0], np.inf
    for a in ALPHAS:
        errs = []
        for tri, vai in GroupKFold(n_splits=INNER_K).split(P, y, groups):
            w, b = _fit_meta(P[tri], y[tri], a)
            errs.append(np.abs(y[vai] - (b + P[vai] @ w)))
        m = float(np.concatenate(errs).mean())
        if m < best_m:
            best_a, best_m = a, m
    return best_a


def fit_predict(tr, te, comp_tr, comp_te):
    tr = tr.reset_index(drop=True)
    te = te.reset_index(drop=True)
    comp_tr = comp_tr.reset_index(drop=True)
    comp_te = comp_te.reset_index(drop=True)
    y = tr.y_out.values.astype(float)
    groups = tr.patient_id.values

    # out-of-fold base predictions inside tr (inner grouped 5-fold by patient)
    P_oof = np.full((len(tr), K_BASE), np.nan)
    for a, b in GroupKFold(n_splits=INNER_K).split(tr, y, groups):
        P_oof[b] = _base_preds(tr.iloc[a].reset_index(drop=True), tr.iloc[b].reset_index(drop=True),
                               comp_tr.iloc[a].reset_index(drop=True), comp_tr.iloc[b].reset_index(drop=True))
    assert np.isfinite(P_oof).all()

    # base predictions for te from models fit on all of tr
    P_te = _base_preds(tr, te, comp_tr, comp_te)

    alpha = _choose_alpha(P_oof, y, groups)
    w, b0 = _fit_meta(P_oof, y, alpha)
    stacked = b0 + P_te @ w

    # v8 collapse gate
    collapse = (tr.core_vox_out < 0.05 * tr.core_vox_in).values
    if collapse.any() and (~collapse).any():
        Xg_tr = _X(tr, comp_tr, v3.latent(tr, tr))
        Xg_te = _X(te, comp_te, v3.latent(te, tr))
        Zt, Ze = _prep(Xg_tr, Xg_te)
        clf = LogisticRegression(max_iter=1000).fit(Zt, collapse)
        p = clf.predict_proba(Ze)[:, 1]
        y_col = float(np.median(y[collapse]))
        return np.where(p > 0.5, y_col, stacked)
    return stacked
