"""Track B hybrid candidate: segmentation-noise-aware forecasting (deconv).

Protocol: paper/trackB_hybrid_search_protocol_2026-11-14.md
Families joined: (1) measurement-noise regression + shrinkage state estimate (errors-in-variables),
(2) previous-scan latent state filter (v3.latent), (3) gradient-boosted L1 residual learner,
(4) logistic collapse gate (v8 style). Fit on tr only; choices by inner GroupKFold on tr.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, QuantileRegressor
from sklearn.model_selection import GroupKFold

import track_b_mu_tcks as t
import track_b_mu_v3 as v3

NOISE_COLS = ["y_in", "compact", "enh_frac", "surf_ratio"]
SIZE_COLS = ["compact", "enh_frac", "surf_ratio", "log_whole"]
LAM_GRID = [0.0, 0.25, 0.5, 0.75, 1.0]
GAMMA_GRID = [0.0, 0.5, 1.0]
SHORT_DT = 60.0
MIN_SHORT = 30


def _frame(df, comp):
    F = pd.DataFrame({"y_in": df.y_in.values.astype(float)})
    return pd.concat([F, comp.reset_index(drop=True)], axis=1)


def _qr(X, y):
    return QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(X, y)


def _fit_stage(tr, F):
    mu, sd = F.mean(), F.std().replace(0, 1.0)
    Z = ((F - mu) / sd).values
    zi = {c: i for i, c in enumerate(F.columns)}
    ncol = [zi[c] for c in NOISE_COLS]
    scol = [zi[c] for c in SIZE_COLS]
    ad = np.abs(tr.y_out.values - tr.y_in.values)
    short = tr.dt_days.values < SHORT_DT
    noise = _qr(Z[short][:, ncol], ad[short]) if short.sum() >= MIN_SHORT else None
    const = float(np.median(ad))
    size = _qr(Z[:, scol], tr.y_in.values)
    m = size.predict(Z[:, scol])
    sig = _sig_from(noise, const, Z, ncol)
    tau2 = max(float(np.var(tr.y_in.values - m) - np.mean(sig ** 2)), 0.02)
    rp = tr.groupby("phase").r_obs.median()
    return dict(mu=mu, sd=sd, zi=zi, ncol=ncol, scol=scol, noise=noise, const=const,
                size=size, tau2=tau2, rp=rp, rall=float(tr.r_obs.median()))


def _sig_from(noise, const, Z, ncol):
    if noise is None:
        return np.full(Z.shape[0], const)
    return np.maximum(noise.predict(Z[:, ncol]), 0.02)


def _apply(st, df, F, lat, lam):
    Z = ((F - st["mu"]) / st["sd"]).values
    m = st["size"].predict(Z[:, st["scol"]])
    sig = _sig_from(st["noise"], st["const"], Z, st["ncol"])
    y = df.y_in.values.astype(float)
    w = st["tau2"] / (st["tau2"] + sig ** 2)
    shrunk = m + w * (y - m)
    dn = (1 - lam) * shrunk + lam * lat
    rate = df.phase.map(st["rp"]).fillna(st["rall"]).values
    base = dn + rate * df.dt_days.values
    return dn, sig, base


def _resid_X(df, F, dn, sig, base):
    X = F.reset_index(drop=True).copy()
    X["dt"] = df.dt_days.values
    X["has_prior"] = df.r_hist.notna().astype(float).values
    X["denoised"] = dn
    X["sigma"] = sig
    X["base"] = base
    for ph in t.PHASES:
        X["ph_" + ph] = (df.phase.values == ph).astype(float)
    return X


def _gbm(X, y):
    return HistGradientBoostingRegressor(loss="absolute_error", random_state=0, max_iter=150,
                                         learning_rate=0.05, max_depth=3, min_samples_leaf=15).fit(X, y)


def _select(tr, F):
    cv = GroupKFold(n_splits=min(5, tr.patient_id.nunique()))
    err_l = np.zeros(len(LAM_GRID))
    splits = list(cv.split(tr, groups=tr.patient_id.values))
    for a, b in splits:
        trI, teI = tr.iloc[a].reset_index(drop=True), tr.iloc[b].reset_index(drop=True)
        Fa, Fb = F.iloc[a].reset_index(drop=True), F.iloc[b].reset_index(drop=True)
        st = _fit_stage(trI, Fa)
        lat = v3.latent(teI, trI)
        for k, lam in enumerate(LAM_GRID):
            _, _, base = _apply(st, teI, Fb, lat, lam)
            err_l[k] += np.abs(base - teI.y_out.values).mean()
    lam = LAM_GRID[int(np.argmin(err_l))]
    err_g = np.zeros(len(GAMMA_GRID))
    for a, b in splits:
        trI, teI = tr.iloc[a].reset_index(drop=True), tr.iloc[b].reset_index(drop=True)
        Fa, Fb = F.iloc[a].reset_index(drop=True), F.iloc[b].reset_index(drop=True)
        st = _fit_stage(trI, Fa)
        dn_a, sg_a, bs_a = _apply(st, trI, Fa, v3.latent(trI, trI), lam)
        g = _gbm(_resid_X(trI, Fa, dn_a, sg_a, bs_a), trI.y_out.values - bs_a)
        dn_b, sg_b, bs_b = _apply(st, teI, Fb, v3.latent(teI, trI), lam)
        r = g.predict(_resid_X(teI, Fb, dn_b, sg_b, bs_b))
        for k, gm in enumerate(GAMMA_GRID):
            err_g[k] += np.abs(bs_b + gm * r - teI.y_out.values).mean()
    gamma = GAMMA_GRID[int(np.argmin(err_g))]
    return lam, gamma


def _gate_X(df, dn, F, mu, sd):
    D = pd.DataFrame({"denoised": dn, "y_in": df.y_in.values, "dt": df.dt_days.values,
                      "has_prior": df.r_hist.notna().astype(float).values})
    for ph in t.PHASES[:-1]:
        D["dt_" + ph] = df.dt_days.values * (df.phase.values == ph)
    X = pd.concat([D, F.reset_index(drop=True).drop(columns=["y_in"])], axis=1)
    return X


def fit_predict(tr, te, comp_tr, comp_te):
    tr = tr.reset_index(drop=True)
    te = te.reset_index(drop=True)
    Ftr_raw = _frame(tr, comp_tr)
    Fte_raw = _frame(te, comp_te)
    med = Ftr_raw.median()
    Ftr = Ftr_raw.fillna(med).fillna(0.0)
    Fte = Fte_raw.fillna(med).fillna(0.0)

    lam, gamma = _select(tr, Ftr)

    st = _fit_stage(tr, Ftr)
    dn_tr, sg_tr, bs_tr = _apply(st, tr, Ftr, v3.latent(tr, tr), lam)
    dn_te, sg_te, bs_te = _apply(st, te, Fte, v3.latent(te, tr), lam)
    g = _gbm(_resid_X(tr, Ftr, dn_tr, sg_tr, bs_tr), tr.y_out.values - bs_tr)
    pred = bs_te + gamma * g.predict(_resid_X(te, Fte, dn_te, sg_te, bs_te))

    collapse_tr = (tr.core_vox_out.values < 0.05 * tr.core_vox_in.values)
    if 0 < collapse_tr.sum() < len(tr):
        Xg_tr = _gate_X(tr, dn_tr, Ftr, None, None)
        Xg_te = _gate_X(te, dn_te, Fte, None, None)
        gmu, gsd = Xg_tr.mean(), Xg_tr.std().replace(0, 1.0)
        clf = LogisticRegression(max_iter=1000).fit((Xg_tr - gmu) / gsd, collapse_tr)
        p = clf.predict_proba((Xg_te - gmu) / gsd)[:, 1]
        y_col = float(np.median(tr.y_out.values[collapse_tr]))
        pred = np.where(p > 0.5, y_col, pred)

    return np.asarray(pred, dtype=float)
