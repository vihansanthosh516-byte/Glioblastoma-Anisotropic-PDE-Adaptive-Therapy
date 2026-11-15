"""Track B hybrid candidate 'pkpd' (protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md).

Families joined into one method:
  1. Mechanistic population model in log-volume: Gompertz growth minus a phase-dependent
     kill rate (phases from track_b_mu_tcks.PHASES), fit by robust soft-L1 least squares on tr.
  2. Median (L1) gradient-boosted residual learner on comp features + latent size x_hat
     (v3.latent) + mechanistic forecast; corrects the mechanistic forecast.
  3. Collapse gate (v8-style): logistic on standardized features; if p > 0.5 predict the
     median training y_out of collapse rows.
Setting choice (GBR depth/size) by inner GroupKFold on tr only, grouped by patient_id.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
import track_b_mu_tcks as t  # noqa: E402
import track_b_mu_v3 as v3  # noqa: E402

PH = t.PHASES
SEED = t.SEED
LR = 0.05
GRID = [{"max_depth": 2, "n_estimators": 150}, {"max_depth": 3, "n_estimators": 150}]


def _phase_idx(df):
    m = {p: i for i, p in enumerate(PH)}
    return np.array([m.get(p, len(PH) - 1) for p in df.phase.values])


def _mech(th, y_in, dt, pi):
    a, ymax = th[0], th[1]
    k = th[2:2 + len(PH)]
    growth = ymax - (ymax - y_in) * np.exp(-a * dt)
    return growth - k[pi] * dt


def _mech_fit(tr):
    y0 = tr.y_in.values.astype(float)
    y1 = tr.y_out.values.astype(float)
    dt = tr.dt_days.values.astype(float)
    pi = _phase_idx(tr)
    lo = [-0.1, 0.0] + [-0.5] * len(PH)
    hi = [0.1, 25.0] + [0.5] * len(PH)
    x0 = [0.0, max(float(y0.max()) + 0.5, 1.0)] + [0.0] * len(PH)
    res = least_squares(lambda th: _mech(th, y0, dt, pi) - y1, x0=x0, bounds=(lo, hi),
                        loss="soft_l1", f_scale=0.1)
    return res.x


def _mech_pred(th, df):
    return _mech(th, df.y_in.values.astype(float), df.dt_days.values.astype(float), _phase_idx(df))


def _design(df, comp, xh, ym):
    X = v3.design(df, xh).reset_index(drop=True)
    X["y_mech"] = ym
    return pd.concat([X, comp.reset_index(drop=True)], axis=1)


def _predict(tr, comp_tr, te, comp_te, setting):
    tr = tr.reset_index(drop=True)
    te = te.reset_index(drop=True)
    th = _mech_fit(tr)
    ym_tr, ym_te = _mech_pred(th, tr), _mech_pred(th, te)
    Xtr = _design(tr, comp_tr, v3.latent(tr, tr), ym_tr)
    Xte = _design(te, comp_te, v3.latent(te, tr), ym_te)
    med = Xtr.median()
    Xtr, Xte = Xtr.fillna(med).fillna(0.0), Xte.fillna(med).fillna(0.0)

    collapse = (tr.core_vox_out < 0.05 * tr.core_vox_in).values
    nc = ~collapse
    resid = tr.y_out.values - ym_tr
    gbr = GradientBoostingRegressor(loss="absolute_error", learning_rate=LR, random_state=SEED,
                                    **setting).fit(Xtr.values[nc], resid[nc])
    hyb = ym_te + gbr.predict(Xte.values)

    mu = Xtr.mean()
    sd = Xtr.std().replace(0, 1.0).fillna(1.0)
    Ztr, Zte = ((Xtr - mu) / sd).values, ((Xte - mu) / sd).values
    if 0 < collapse.sum() < len(collapse):
        clf = LogisticRegression(max_iter=1000).fit(Ztr, collapse)
        p = clf.predict_proba(Zte)[:, 1]
        y_col = float(np.median(tr.y_out.values[collapse]))
        return np.where(p > 0.5, y_col, hyb)
    return hyb


def _select(tr, comp_tr):
    tr = tr.reset_index(drop=True)
    comp_tr = comp_tr.reset_index(drop=True)
    gkf = GroupKFold(n_splits=3)
    best, best_s = GRID[0], np.inf
    for s in GRID:
        errs = []
        for a, b in gkf.split(tr, groups=tr.patient_id.values):
            trI, vaI = tr.iloc[a].reset_index(drop=True), tr.iloc[b].reset_index(drop=True)
            p = _predict(trI, comp_tr.iloc[a], vaI, comp_tr.iloc[b], s)
            errs.append(np.abs(p - vaI.y_out.values))
        sc = float(np.concatenate(errs).mean())
        if sc < best_s:
            best, best_s = s, sc
    return best


def fit_predict(tr, te, comp_tr, comp_te):
    setting = _select(tr, comp_tr)
    return np.asarray(_predict(tr, comp_tr, te, comp_te, setting), dtype=float)
