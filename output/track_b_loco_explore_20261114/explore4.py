"""Exploration 4 on DEVELOPMENT cohorts only (MU, UCSD, UCSF). Post hoc.

Refines the growth-law x visit-weight hybrid ("gv"):
  change = w(dt) * rho * dt + (1 - w(dt)) * dq,   w(dt) = 1 / (1 + (dt / tau)^kappa)
Variants: free rho; rho fixed to 0.014/day (untreated GBM growth, Stensjoen 2015, value from memory);
exponential visit weight; and an added mean-reversion term on the last change (history).
Leave-one-cohort-out and correction curves. LUMIERE is not read here.
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from explore1 import pool
from explore3 import ci

RHO_LIT = 0.014
DTS = np.array([20, 30, 45, 60, 90, 120, 180, 300.0])


def curve(th, dt, kind):
    if kind == "free":
        rho, dq, tau, kap = np.exp(th[0]), th[1], np.exp(th[2]), np.exp(th[3])
        w = 1 / (1 + (dt / tau) ** kap)
    elif kind == "fix":
        rho, dq, tau, kap = RHO_LIT, th[0], np.exp(th[1]), np.exp(th[2])
        w = 1 / (1 + (dt / tau) ** kap)
    else:  # exponential weight
        rho, dq, tau = np.exp(th[0]), th[1], np.exp(th[2])
        w = np.exp(-dt / tau)
    return w * rho * dt + (1 - w) * dq


STARTS = {"free": [[np.log(0.01), 0.0, np.log(60), 0.7], [np.log(0.005), -0.1, np.log(90), 1.0], [np.log(0.02), 0.0, np.log(40), 0.3]],
          "fix": [[0.0, np.log(60), 0.7], [-0.1, np.log(90), 1.0], [0.0, np.log(40), 0.3], [0.0, np.log(120), 0.0]],
          "exp": [[np.log(0.01), 0.0, np.log(60)], [np.log(0.02), -0.1, np.log(30)], [np.log(0.005), 0.0, np.log(120)]]}


def fit(tr, kind):
    d, dt = tr.ch.values, tr.dt_days.values
    f = lambda th: np.mean(np.abs(d - curve(th, dt, kind)))
    return min((minimize(f, x0, method="Nelder-Mead", options={"xatol": 1e-5, "fatol": 1e-8, "maxiter": 6000}) for x0 in STARTS[kind]), key=lambda r: r.fun).x


def fit_lam(tr, base):
    h = tr[tr.r_hist.notna()]
    if len(h) < 30:
        return 0.0
    res, lc = (h.ch - base[tr.r_hist.notna().values]).values, h.lastch.values
    grid = np.linspace(0, 0.5, 51)
    return float(grid[np.argmin([np.mean(np.abs(res + l * lc)) for l in grid])])


def main():
    p = pool()
    p["ch"] = p.y_out - p.y_in
    p["lastch"] = (p.r_hist * p.dt_prior).fillna(0.0)
    cohorts = ["MU", "UCSD", "UCSF"]
    kinds = ["free", "fix", "exp"]
    pats = {c: {} for c in cohorts}
    for te_c in cohorts:
        tr, te = p[p.cohort != te_c].reset_index(drop=True), p[p.cohort == te_c].reset_index(drop=True)
        y, pid = te.y_out.values, te.patient_id.values
        preds = {"persistence": te.y_in.values, "pop_rate": te.y_in.values + float((tr.ch / tr.dt_days).median()) * te.dt_days.values}
        for k in kinds:
            th = fit(tr, k)
            preds[k] = te.y_in.values + curve(th, te.dt_days.values, k)
            lam = fit_lam(tr, curve(th, tr.dt_days.values, k))
            preds[k + "+h"] = preds[k] - lam * te.lastch.values
            print(f"test {te_c} {k:4s} curve at dt {DTS.astype(int).tolist()}: {np.round(curve(th, DTS, k), 3).tolist()} lam={lam:.2f}")
        pats[te_c] = pd.DataFrame({k: np.abs(v - y) for k, v in preds.items()}).assign(pid=pid).groupby("pid").mean()
    allp = pd.concat(pats.values())
    print("\nLOCO patient MAE minus persistence [95% CI] | minus pop_rate [95% CI]")
    for k in [a + b for a in kinds for b in ("", "+h")]:
        line = f"{k:7s}"
        for c in cohorts + ["POOLED"]:
            pt = allp if c == "POOLED" else pats[c]
            d1, d2 = (pt[k] - pt.persistence).values, (pt[k] - pt.pop_rate).values
            line += f" | {c} {d1.mean():+.3f} {ci(d1)} ; {d2.mean():+.3f} {ci(d2)}"
        print(line)
    th = fit(p, "free"); print("\nall-development fit, free:", np.round([np.exp(th[0]), th[1], np.exp(th[2]), np.exp(th[3])], 4).tolist(), "curve", np.round(curve(th, DTS, "free"), 3).tolist())
    th = fit(p, "fix"); print("all-development fit, fix :", np.round([th[0], np.exp(th[1]), np.exp(th[2])], 4).tolist(), "curve", np.round(curve(th, DTS, "fix"), 3).tolist())
    # history term inside MU only (5-fold patient CV), where most history pairs are
    g = p[p.cohort == "MU"].reset_index(drop=True)
    rng = np.random.default_rng(0); pid = g.patient_id.unique(); fold = g.patient_id.map(dict(zip(pid, rng.permutation(len(pid)) % 5))).values
    e0, e1 = np.zeros(len(g)), np.zeros(len(g))
    for k in range(5):
        tr, te = g[fold != k], g[fold == k]
        th = fit(tr, "fix"); lam = fit_lam(tr.reset_index(drop=True), curve(th, tr.dt_days.values, "fix"))
        b = te.y_in.values + curve(th, te.dt_days.values, "fix")
        e0[fold == k], e1[fold == k] = np.abs(b - te.y_out.values), np.abs(b - lam * te.lastch.values - te.y_out.values)
    pm = pd.DataFrame({"fix": e0, "fix+h": e1, "pers": (g.ch).abs().values}).assign(pid=g.patient_id.values).groupby("pid").mean()
    print("MU 5-fold CV: fix minus persistence", round(float((pm.fix - pm.pers).mean()), 3), ci((pm.fix - pm.pers).values), "| fix+h minus fix", round(float((pm["fix+h"] - pm.fix).mean()), 3), ci((pm["fix+h"] - pm.fix).values))


if __name__ == "__main__":
    main()
