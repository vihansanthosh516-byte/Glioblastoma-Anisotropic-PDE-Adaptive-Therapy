"""Track B hybrid candidate 'compartment' (one fixed method).

Protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md
Families joined in one method:
  1. Mechanistic three-compartment kinetic model, linear ODE per phase:
       E = enhancing core voxels, N = necrotic core voxels, F = edema reservoir.
       dE/dt = (rho_ph - k) E + s F ;  dN/dt = k E - c N ; F held constant over the interval.
       Initial E, N from the input scan (core_vox_in * enh_frac, core_vox_in * nec_frac);
       F = expm1(log_edema). Forecast y = log(E + N + 1) at day_out (dt_days ahead).
       Fit: L1 of log-volume on tr only (Powell, bounded).
  2. L1 residual learner on input-scan composition features + ODE forecast
       (QuantileRegressor or HistGradientBoosting absolute_error), fit on non-collapse tr rows.
  3. v8-style collapse gate (logistic; collapse = core_vox_out < 0.05 * core_vox_in).
       p > 0.5 -> median tr y_out of collapse rows.
Settings: residual learner chosen by inner GroupKFold on patient_id inside tr only.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, QuantileRegressor
from sklearn.model_selection import GroupKFold

import track_b_mu_tcks as t

PH = list(t.PHASES)
SEED = t.SEED
N_SUB = 40
LEARNERS = ["qr", "hgb"]
BOUNDS = [(-0.05, 0.1)] * len(PH) + [(0.0, 0.1), (0.0, 0.1), (0.0, 0.05)]
X0 = np.array([0.01] * len(PH) + [0.01, 0.01, 1e-4])


def _states(df, comp):
    core = df.core_vox_in.values.astype(float)
    E = core * comp.enh_frac.fillna(0.0).values
    N = core * comp.nec_frac.fillna(0.0).values
    F = np.clip(np.expm1(comp.log_edema.fillna(0.0).values), 0.0, None)
    return E, N, F


def _rho_vec(df, rho):
    idx = np.array([PH.index(p) if p in PH else PH.index("other") for p in df.phase.values])
    return rho[idx]


def _simulate(theta, df, comp):
    rho = theta[: len(PH)]
    k, c, s = theta[len(PH):]
    E, N, F = _states(df, comp)
    r = _rho_vec(df, rho)
    h = df.dt_days.values.astype(float) / N_SUB

    def f(e, n):
        return (r - k) * e + s * F, k * e - c * n

    for _ in range(N_SUB):
        k1e, k1n = f(E, N)
        k2e, k2n = f(E + 0.5 * h * k1e, N + 0.5 * h * k1n)
        k3e, k3n = f(E + 0.5 * h * k2e, N + 0.5 * h * k2n)
        k4e, k4n = f(E + h * k3e, N + h * k3n)
        E = E + h / 6.0 * (k1e + 2 * k2e + 2 * k3e + k4e)
        N = N + h / 6.0 * (k1n + 2 * k2n + 2 * k3n + k4n)
    return np.log(np.clip(E, 0.0, None) + np.clip(N, 0.0, None) + 1.0)


def _fit_ode(tr, comp):
    y = tr.y_out.values.astype(float)

    def loss(th):
        return float(np.mean(np.abs(_simulate(th, tr, comp) - y)))

    res = minimize(loss, X0, method="Powell", bounds=BOUNDS,
                   options={"maxfev": 3000, "xtol": 1e-4, "ftol": 1e-6})
    return res.x


def _feat(df, comp, ode):
    X = comp.select_dtypes(include=[np.number]).copy().reset_index(drop=True)
    X["y_in"] = df.y_in.values
    X["dt"] = df.dt_days.values
    X["ode"] = ode
    for ph in PH:
        X["ph_" + ph] = (df.phase.values == ph).astype(float)
    return X


def _core(tr, te, comp_tr, comp_te, theta, learner):
    ode_tr = _simulate(theta, tr, comp_tr)
    ode_te = _simulate(theta, te, comp_te)
    Xtr, Xte = _feat(tr, comp_tr, ode_tr), _feat(te, comp_te, ode_te)
    keep = Xtr.columns[Xtr.notna().any()]
    Xtr, Xte = Xtr[keep], Xte[keep]
    med = Xtr.median()
    Xtr, Xte = Xtr.fillna(med).fillna(0.0), Xte.fillna(med).fillna(0.0)
    mu, sd = Xtr.mean(), Xtr.std().replace(0, 1.0).fillna(1.0)
    Ztr, Zte = ((Xtr - mu) / sd).values, ((Xte - mu) / sd).values

    y = tr.y_out.values.astype(float)
    c_tr = tr.core_vox_out.values < 0.05 * tr.core_vox_in.values
    if c_tr.any() and (~c_tr).any():
        clf = LogisticRegression(max_iter=1000).fit(Ztr, c_tr)
        p = clf.predict_proba(Zte)[:, 1]
    else:
        p = np.zeros(len(te))
    nc = ~c_tr
    rtr = y[nc] - ode_tr[nc]
    if learner == "qr":
        m = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Ztr[nc], rtr)
    else:
        m = HistGradientBoostingRegressor(loss="absolute_error", max_iter=150, learning_rate=0.05,
                                          max_depth=3, min_samples_leaf=10,
                                          random_state=SEED).fit(Ztr[nc], rtr)
    base = ode_te + m.predict(Zte)
    y_col = float(np.median(y[c_tr])) if c_tr.any() else 0.0
    return np.where(p > 0.5, y_col, base)


def fit_predict(tr, te, comp_tr, comp_te):
    tr = tr.reset_index(drop=True)
    te = te.reset_index(drop=True)
    comp_tr = comp_tr.reset_index(drop=True)
    comp_te = comp_te.reset_index(drop=True)

    gkf = GroupKFold(n_splits=min(4, tr.patient_id.nunique()))
    errs = {L: [] for L in LEARNERS}
    for itr, iva in gkf.split(tr, groups=tr.patient_id.values):
        a, b = tr.iloc[itr].reset_index(drop=True), tr.iloc[iva].reset_index(drop=True)
        ca, cb = comp_tr.iloc[itr].reset_index(drop=True), comp_tr.iloc[iva].reset_index(drop=True)
        th = _fit_ode(a, ca)
        for L in LEARNERS:
            p = _core(a, b, ca, cb, th, L)
            errs[L].append(np.abs(p - b.y_out.values))
    best = min(LEARNERS, key=lambda L: float(np.mean(np.concatenate(errs[L]))))

    th = _fit_ode(tr, comp_tr)
    return _core(tr, te, comp_tr, comp_te, th, best)
