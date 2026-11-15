"""Track B hybrid candidate 'switch': regime-switching state-space forecast with a
multinomial regime classifier. Median-type output (MAE score).

Protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md
Families joined: (1) state-space / Kalman latent-size filter per regime,
(2) multinomial logistic regime classifier on input-scan features.
Fit uses training folds only. The one setting (classifier C) is chosen by inner
GroupKFold CV on patient_id inside the training folds.
Output: median of the regime mixture of Laplace predictive components.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GroupKFold

import track_b_mu_tcks as t
import track_b_mu_v3 as v3

STABLE, GROWTH, COLLAPSE = 0, 1, 2
N_REG = 3
LOG_GROWTH = 0.22          # regime threshold on log-ratio y_out - y_in
COLLAPSE_FRAC = 0.05       # collapse if core_vox_out < 0.05 * core_vox_in
C_GRID = [0.01, 0.1, 1.0]  # inner-CV grid for classifier regularisation
MIN_N = 5                  # min rows to estimate a regime rate, else pooled rate
B_FLOOR = 1e-3             # floor on Laplace scale
BISECT_IT = 100


def regime_labels(df):
    lr = df.y_out.values - df.y_in.values
    lab = np.full(len(df), STABLE, dtype=int)
    lab[lr > LOG_GROWTH] = GROWTH
    lab[df.core_vox_out.values < COLLAPSE_FRAC * df.core_vox_in.values] = COLLAPSE
    return lab


def _latent_k(df, ck, t2, s2):
    """Kalman-style filtered size at day_in under regime drift ck (same gain as v3.latent)."""
    has = df.r_hist.notna().values
    dtp = df.dt_prior.fillna(1.0).values
    y_in = df.y_in.values
    x_prior = y_in - df.r_hist.fillna(0.0).values * dtp + ck * dtp
    P = s2 + t2 * dtp ** 2
    x = x_prior + P / (P + s2) * (y_in - x_prior)
    return np.where(has, x, y_in)


def _std_stats(tr, comp, xh):
    cc = comp.apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    keep = cc.columns[cc.notna().any()]
    cc = cc[keep]
    med = cc.median()
    filled = cc.fillna(med)
    mu_c = filled.mean()
    sd_c = filled.std().replace(0, 1.0).fillna(1.0)
    num = np.column_stack([tr.y_in.values, tr.dt_days.values, xh])
    return {"keep": keep, "med": med, "mu_c": mu_c, "sd_c": sd_c,
            "mu_n": num.mean(0), "sd_n": np.where(num.std(0) > 0, num.std(0), 1.0)}


def _phase_dummies(df):
    return np.column_stack([(df.phase.values == p).astype(float) for p in t.PHASES[:-1]])


def _cls_X(df, comp, xh, st):
    cc = comp[st["keep"]].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan)
    cz = (cc.fillna(st["med"]) - st["mu_c"]) / st["sd_c"]
    num = np.column_stack([df.y_in.values, df.dt_days.values, xh])
    nz = (num - st["mu_n"]) / st["sd_n"]
    return np.column_stack([cz.values, nz, _phase_dummies(df)])


def _mix_median(P, mus, bs):
    """Median of sum_k P_k Laplace(mus_k, bs_k), per row, by bisection on the CDF."""
    bmax = bs.max()
    lo = mus.min(axis=1) - 60 * bmax - 1.0
    hi = mus.max(axis=1) + 60 * bmax + 1.0

    def cdf(y):
        z = (y[:, None] - mus) / bs[None, :]
        F = np.where(z < 0, 0.5 * np.exp(np.minimum(z, 0)), 1 - 0.5 * np.exp(-np.maximum(z, 0)))
        return (P * F).sum(axis=1)

    for _ in range(BISECT_IT):
        mid = 0.5 * (lo + hi)
        below = cdf(mid) < 0.5
        lo = np.where(below, mid, lo)
        hi = np.where(below, hi, mid)
    return 0.5 * (lo + hi)


def _predict_mix(tr, te, comp_tr, comp_te, C):
    lab = regime_labels(tr)
    r = tr.r_obs.values
    r_pop = float(np.median(r))
    c_reg = np.array([float(np.median(r[lab == k])) if (lab == k).sum() >= MIN_N else r_pop
                      for k in range(N_REG)])
    t2, s2 = t.eb_params(tr, c_reg[lab])
    off, bs = np.zeros(N_REG), np.zeros(N_REG)
    for k in range(N_REG):
        mu_tr = _latent_k(tr, c_reg[k], t2, s2) + c_reg[k] * tr.dt_days.values
        e = tr.y_out.values - mu_tr
        off[k] = float(np.median(e))
        bs[k] = max(float(np.mean(np.abs(e - off[k]))), B_FLOOR)
    xh_tr, xh_te = v3.latent(tr, tr), v3.latent(te, tr)
    st = _std_stats(tr, comp_tr, xh_tr)
    clf = LogisticRegression(C=C, max_iter=2000).fit(_cls_X(tr, comp_tr, xh_tr, st), lab)
    P = np.zeros((len(te), N_REG))
    P[:, clf.classes_] = clf.predict_proba(_cls_X(te, comp_te, xh_te, st))
    mus = np.column_stack([_latent_k(te, c_reg[k], t2, s2) + c_reg[k] * te.dt_days.values + off[k]
                           for k in range(N_REG)])
    return _mix_median(P, mus, bs)


def _choose_C(tr, comp_tr):
    gkf = GroupKFold(n_splits=4)
    groups = tr.patient_id.values
    best_C, best_mae = C_GRID[0], np.inf
    for C in C_GRID:
        rows = []
        for a, b in gkf.split(tr, groups=groups):
            tra, teb = tr.iloc[a].reset_index(drop=True), tr.iloc[b].reset_index(drop=True)
            m = _predict_mix(tra, teb, comp_tr.iloc[a].reset_index(drop=True),
                             comp_tr.iloc[b].reset_index(drop=True), C)
            rows.append(pd.DataFrame({"pid": teb.patient_id.values,
                                      "err": np.abs(m - teb.y_out.values)}))
        mae = float(pd.concat(rows).groupby("pid").err.mean().mean())
        if mae < best_mae:
            best_C, best_mae = C, mae
    return best_C


def fit_predict(tr, te, comp_tr, comp_te):
    C = _choose_C(tr, comp_tr)
    return np.asarray(_predict_mix(tr, te, comp_tr, comp_te, C), dtype=float)
