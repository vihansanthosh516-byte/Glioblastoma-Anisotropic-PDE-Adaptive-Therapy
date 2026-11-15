"""Exploration 3 on DEVELOPMENT cohorts only (MU, UCSD, UCSF). Post hoc.

Candidate hybrids built on the scan-interval signal:
  ldt   : median regression of change on log(dt)                      (reference, not a hybrid)
  gv    : growth law x visit weight, fitted with absolute loss
  mix   : joint mixture (growth regime / quiet regime) of size change and scan interval,
          EM fit, forecast = median of the mixture with a vanish point mass
  +h    : same with a mean-reversion term on the last observed change
Leave-one-cohort-out (train on two cohorts, score the third). LUMIERE is not read here.
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm
from sklearn.linear_model import QuantileRegressor, LogisticRegression

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from explore1 import pool, pmae

GRID = np.linspace(-12, 8, 4001)


def lap_pdf(x, loc, b):
    return np.exp(-np.abs(x - loc) / b) / (2 * b)


def lap_cdf(x, loc, b):
    z = (x - loc) / b
    return np.where(z < 0, 0.5 * np.exp(z), 1 - 0.5 * np.exp(-z))


def wmedian(v, w):
    o = np.argsort(v); cw = np.cumsum(w[o])
    return float(v[o][np.searchsorted(cw, cw[-1] / 2.0)])


def fit_mix(tr, n_iter=200):
    """EM for the joint model of (change, log interval) with two regimes. Pairs with y_out > 0 only."""
    g = tr[tr.core_vox_out > 0]
    d, dt, ldt, s = g.ch.values, g.dt_days.values, np.log(g.dt_days.values), g.y_in.values - 8.0
    c = float(np.polyfit(s, np.log(np.abs(d) + 0.05), 1)[0])          # size slope of the noise scale
    sc = np.exp(c * s)
    r = np.where(d > 0.2, 0.8, 0.2)
    th = None
    for _ in range(n_iter):
        pi = float(r.mean())
        mG, mQ = np.average(ldt, weights=r), np.average(ldt, weights=1 - r)
        vG = max(np.average((ldt - mG) ** 2, weights=r), 0.02); vQ = max(np.average((ldt - mQ) ** 2, weights=1 - r), 0.02)
        rho = max(wmedian(d / dt, r * dt / sc), 0.0)
        dq = wmedian(d, (1 - r) / sc)
        bG = max(float(np.sum(r * np.abs(d - rho * dt) / sc) / r.sum()), 0.05)
        bQ = max(float(np.sum((1 - r) * np.abs(d - dq) / sc) / (1 - r).sum()), 0.05)
        a = pi * lap_pdf(d, rho * dt, bG * sc) * norm.pdf(ldt, mG, np.sqrt(vG))
        b = (1 - pi) * lap_pdf(d, dq, bQ * sc) * norm.pdf(ldt, mQ, np.sqrt(vQ))
        r_new = a / (a + b + 1e-300)
        new = (pi, mG, vG, mQ, vQ, rho, dq, bG, bQ, c)
        if th is not None and np.max(np.abs(np.array(new) - np.array(th))) < 1e-7:
            th = new; break
        th, r = new, r_new
    van = LogisticRegression(C=1.0).fit(np.c_[tr.y_in.values, np.log(tr.dt_days.values)], (tr.core_vox_out == 0).astype(int))
    return th, van


def pred_mix(te, th, van, lam=0.0):
    pi, mG, vG, mQ, vQ, rho, dq, bG, bQ, c = th
    dt, ldt, sc = te.dt_days.values, np.log(te.dt_days.values), np.exp(c * (te.y_in.values - 8.0))
    a = pi * norm.pdf(ldt, mG, np.sqrt(vG)); b = (1 - pi) * norm.pdf(ldt, mQ, np.sqrt(vQ))
    w = a / (a + b)
    p0 = van.predict_proba(np.c_[te.y_in.values, ldt])[:, 1]
    out = np.zeros(len(te))
    for i in range(len(te)):
        F = p0[i] * (GRID >= -te.y_in.values[i]) + (1 - p0[i]) * (w[i] * lap_cdf(GRID, rho * dt[i], bG * sc[i]) + (1 - w[i]) * lap_cdf(GRID, dq, bQ * sc[i]))
        out[i] = GRID[np.searchsorted(F, 0.5)]
    return out - lam * te.lastch.values, w


def fit_gv(tr):
    d, dt = tr.ch.values, tr.dt_days.values
    def f(th):
        rho, dq, tau, kap = np.exp(th[0]), th[1], np.exp(th[2]), np.exp(th[3])
        w = 1 / (1 + (dt / tau) ** kap)
        return np.mean(np.abs(d - (w * rho * dt + (1 - w) * dq)))
    best = min((minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 4000}) for x0 in
                ([np.log(0.01), 0.0, np.log(60), 0.7], [np.log(0.005), -0.1, np.log(90), 1.0], [np.log(0.02), 0.0, np.log(40), 0.3])), key=lambda r: r.fun)
    return best.x


def pred_gv(te, th):
    rho, dq, tau, kap = np.exp(th[0]), th[1], np.exp(th[2]), np.exp(th[3])
    w = 1 / (1 + (te.dt_days.values / tau) ** kap)
    return w * rho * te.dt_days.values + (1 - w) * dq


def ci(diff, level=0.95, n=20000, seed=20261114):
    rng = np.random.default_rng(seed)
    b = diff[rng.integers(0, len(diff), (n, len(diff)))].mean(1)
    a = (1 - level) / 2 * 100
    return round(float(np.percentile(b, a)), 3), round(float(np.percentile(b, 100 - a)), 3)


def main():
    p = pool()
    p["ch"] = p.y_out - p.y_in
    p["lastch"] = (p.r_hist * p.dt_prior).fillna(0.0)
    p["ldt"] = np.log(p.dt_days / 60.0)
    cohorts = ["MU", "UCSD", "UCSF"]
    rows, allpat = [], {}
    for te_c in cohorts:
        tr, te = p[p.cohort != te_c].reset_index(drop=True), p[p.cohort == te_c].reset_index(drop=True)
        y, pid = te.y_out.values, te.patient_id.values
        preds = {"persistence": te.y_in.values,
                 "pop_rate": te.y_in.values + float(((tr.ch) / tr.dt_days).median()) * te.dt_days.values}
        q = QuantileRegressor(quantile=0.5, alpha=0.0, solver="highs").fit(tr[["ldt"]].values, tr.ch.values)
        preds["ldt"] = te.y_in.values + q.predict(te[["ldt"]].values)
        tg = fit_gv(tr); preds["gv"] = te.y_in.values + pred_gv(te, tg)
        th, van = fit_mix(tr)
        m0, w = pred_mix(te, th, van)
        preds["mix"] = te.y_in.values + m0
        hist = tr[tr.r_hist.notna()]
        lam = float(np.clip(-wmedian(((hist.ch - pred_mix(hist, th, van)[0]) / hist.lastch.replace(0, np.nan)).fillna(0).values, np.abs(hist.lastch.values) + 1e-9), 0, 0.5)) if len(hist) > 20 else 0.0
        preds["mix+h"] = te.y_in.values + pred_mix(te, th, van, lam)[0]
        pat = pd.DataFrame({k: np.abs(v - y) for k, v in preds.items()}).assign(pid=pid).groupby("pid").mean()
        allpat[te_c] = pat
        print(f"\n== test {te_c}: n patients {len(pat)} | mix params pi={th[0]:.2f} mG={np.exp(th[1]):.0f}d mQ={np.exp(th[3]):.0f}d rho={th[5]:.4f}/d dq={th[6]:.3f} bG={th[7]:.2f} bQ={th[8]:.2f} c={th[9]:.2f} lam={lam:.2f}"
              f" | gv rho={np.exp(tg[0]):.4f} dq={tg[1]:.3f} tau={np.exp(tg[2]):.0f} kap={np.exp(tg[3]):.2f} | mean growth weight w={w.mean():.2f}")
        print("   patient MAE", pat.mean().round(3).to_dict())
        for k in ["ldt", "gv", "mix", "mix+h"]:
            print(f"   {k:6s} minus persistence {float((pat[k] - pat.persistence).mean()):+.3f} CI {ci((pat[k] - pat.persistence).values)} | minus pop_rate {float((pat[k] - pat.pop_rate).mean()):+.3f} CI {ci((pat[k] - pat.pop_rate).values)} | sd diff {float((pat[k] - pat.persistence).std()):.3f}")
    allp = pd.concat(allpat.values())
    print("\n== pooled over the three held-out cohorts (n patients", len(allp), ")")
    for k in ["ldt", "gv", "mix", "mix+h"]:
        print(f"   {k:6s} minus persistence {float((allp[k] - allp.persistence).mean()):+.3f} CI {ci((allp[k] - allp.persistence).values)} | minus pop_rate {float((allp[k] - allp.pop_rate).mean()):+.3f} CI {ci((allp[k] - allp.pop_rate).values)}")


if __name__ == "__main__":
    main()
