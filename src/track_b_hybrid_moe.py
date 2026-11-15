"""Track B hybrid candidate 'moe': mechanistic-expert gated weighted median + v8 collapse override.

Protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md
Experts (each predicts y_out): (0) persistence, (1) population median rate, (2) Gompertz fit,
(3) v3 latent size + phase median rate, (4) v8-style median regression on non-collapse rows.
Gate: multinomial logistic on standardized comp + y_in + dt + phase dummies, trained to predict
which expert had the smallest absolute error on inner out-of-fold predictions (inner 5-fold, grouped by patient).
Forecast: weighted median of expert predictions with gate probabilities as weights.
Override: v8 collapse logistic (collapse = core_vox_out < 0.05*core_vox_in); p > 0.5 -> median tr y_out of collapse rows.
Fits use training folds only.
"""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression, QuantileRegressor
from sklearn.model_selection import GroupKFold

import track_b_mu_tcks as t
import track_b_mu_v3 as v3

K = 5  # number of experts
INNER_K = 5
COLLAPSE_FRAC = 0.05


def _num(comp):
    return pd.DataFrame(comp).apply(pd.to_numeric, errors="coerce").to_numpy(float)


def _standardize(Xtr, Xte):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        mu = np.nanmean(Xtr, axis=0)
        sd = np.nanstd(Xtr, axis=0)
    mu = np.where(np.isfinite(mu), mu, 0.0)
    sd = np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)
    Ztr = (np.where(np.isnan(Xtr), mu, Xtr) - mu) / sd
    Zte = (np.where(np.isnan(Xte), mu, Xte) - mu) / sd
    return Ztr, Zte


def _design_comp(df, ref, comp):
    X = v3.design(df, v3.latent(df, ref)).to_numpy(float)
    return np.column_stack([X, comp])


def _experts(A, B, cA, cB):
    """Fit experts on A, predict y_out for B. Returns (len(B), K)."""
    r_pop = float(A.r_obs.median())
    r_ph = A.groupby("phase").r_obs.median().reindex(t.PHASES).fillna(r_pop)
    y_in = B.y_in.values.astype(float)
    dt = B.dt_days.values.astype(float)
    c_b = B.phase.map(r_ph).values.astype(float)
    E = np.zeros((len(B), K))
    E[:, 0] = y_in
    E[:, 1] = y_in + r_pop * dt
    E[:, 2] = t.gompertz_fit(A)(B)
    E[:, 3] = v3.latent(B, A) + c_b * dt
    nc = A.core_vox_out.values >= COLLAPSE_FRAC * A.core_vox_in.values
    Za, Zb = _standardize(_design_comp(A, A, cA), _design_comp(B, A, cB))
    reg = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(Za[nc], A.y_out.values[nc])
    E[:, 4] = reg.predict(Zb)
    return np.where(np.isfinite(E), E, y_in[:, None])


def _gate_feats(df, comp):
    return np.column_stack([comp, df.y_in.values, df.dt_days.values]), \
        np.column_stack([(df.phase.values == ph).astype(float) for ph in t.PHASES])


def _gate_probs(A, cA, oof, B, cB):
    lab = np.argmin(np.abs(oof - A.y_out.values[:, None]), axis=1)
    Na, Pa = _gate_feats(A, cA)
    Nb, Pb = _gate_feats(B, cB)
    Ga, Gb = _standardize(Na, Nb)
    Ga, Gb = np.column_stack([Ga, Pa]), np.column_stack([Gb, Pb])
    P = np.zeros((len(B), K))
    classes = np.unique(lab)
    if len(classes) == 1:
        P[:, classes[0]] = 1.0
        return P
    clf = LogisticRegression(C=1.0, max_iter=5000).fit(Ga, lab)
    P[:, clf.classes_] = clf.predict_proba(Gb)
    return P


def _wmedian(V, W):
    order = np.argsort(V, axis=1)
    Vs = np.take_along_axis(V, order, 1)
    Ws = np.take_along_axis(W, order, 1)
    cw = np.cumsum(Ws, axis=1)
    idx = np.argmax(cw >= 0.5 * cw[:, -1:] - 1e-12, axis=1)
    return Vs[np.arange(len(V)), idx]


def fit_predict(tr, te, comp_tr, comp_te):
    tr = tr.reset_index(drop=True)
    te = te.reset_index(drop=True)
    ctr, cte = _num(comp_tr), _num(comp_te)

    # out-of-fold expert predictions inside tr (inner 5-fold grouped by patient)
    oof = np.zeros((len(tr), K))
    for a, b in GroupKFold(n_splits=INNER_K).split(np.zeros(len(tr)), groups=tr.patient_id.values):
        A = tr.iloc[a].reset_index(drop=True)
        B = tr.iloc[b].reset_index(drop=True)
        oof[b] = _experts(A, B, ctr[a], ctr[b])

    P_te = _gate_probs(tr, ctr, oof, te, cte)
    E_te = _experts(tr, te, ctr, cte)
    med = _wmedian(E_te, P_te)

    # v8 collapse override
    collapse_tr = tr.core_vox_out.values < COLLAPSE_FRAC * tr.core_vox_in.values
    Zc_tr, Zc_te = _standardize(_design_comp(tr, tr, ctr), _design_comp(te, tr, cte))
    if 0 < collapse_tr.sum() < len(collapse_tr):
        p = LogisticRegression(max_iter=1000).fit(Zc_tr, collapse_tr).predict_proba(Zc_te)[:, 1]
    else:
        p = np.full(len(te), float(collapse_tr.mean()))
    y_col = float(np.median(tr.y_out.values[collapse_tr])) if collapse_tr.any() else float(np.median(tr.y_out.values))
    return np.where(p > 0.5, y_col, med)
