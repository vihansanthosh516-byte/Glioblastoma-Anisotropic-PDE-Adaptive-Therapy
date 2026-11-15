"""Exploration 5 on DEVELOPMENT cohorts only. Post hoc. History part of the hybrid.

Uses MU + UCSD (the cohorts with earlier scans). 5-fold patient CV, repeated 5 times.
Base = interval curve (exp form) fitted on the training folds (all three cohorts' training pairs).
Variants of the history term on top of the base. LUMIERE is not read here.
"""
import sys, warnings
from pathlib import Path
import numpy as np, pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
from explore1 import pool
from explore3 import ci
from explore4 import fit, curve


def l1_fit(res, X, grid):
    """coordinate grid search for coefficients b >= 0 minimising mean |res + X b|"""
    b = np.zeros(X.shape[1])
    for _ in range(4):
        for j in range(X.shape[1]):
            cand = []
            for g in grid:
                bb = b.copy(); bb[j] = g
                cand.append(np.mean(np.abs(res + X @ bb)))
            b[j] = grid[int(np.argmin(cand))]
    return b


def design(d, kind):
    lc = d.lastch.values
    if kind == "h1":
        return np.c_[lc]
    if kind == "h2":      # different pull after a rise and after a fall
        return np.c_[np.maximum(lc, 0), np.minimum(lc, 0)]
    if kind == "h3":      # different pull for small and large lesions
        s = (d.y_in.values < 8.0)
        return np.c_[lc * s, lc * ~s]
    if kind == "h4":      # stronger pull when the earlier interval was short
        s = (d.dt_prior.fillna(999).values < 60)
        return np.c_[lc * s, lc * ~s]
    raise ValueError(kind)


def main():
    p = pool()
    p["ch"] = p.y_out - p.y_in
    p["lastch"] = (p.r_hist * p.dt_prior).fillna(0.0)
    kinds = ["h1", "h2", "h3", "h4"]
    grid = np.linspace(0, 0.6, 31)
    out = []
    for rep in range(5):
        rng = np.random.default_rng(rep)
        pid = p.patient_id.unique(); fold = p.patient_id.map(dict(zip(pid, rng.permutation(len(pid)) % 5))).values
        pred = {k: np.zeros(len(p)) for k in ["base"] + kinds}
        coefs = {k: [] for k in kinds}
        for f in range(5):
            tr, te = p[fold != f], p[fold == f]
            th = fit(tr, "exp")
            btr, bte = curve(th, tr.dt_days.values, "exp"), curve(th, te.dt_days.values, "exp")
            pred["base"][fold == f] = te.y_in.values + bte
            has = tr.r_hist.notna().values
            for k in kinds:
                b = l1_fit((tr.ch.values - btr)[has], design(tr[has], k), grid)
                coefs[k].append(b)
                pred[k][fold == f] = te.y_in.values + bte - design(te, k) @ b
        e = pd.DataFrame({k: np.abs(v - p.y_out.values) for k, v in pred.items()})
        e["pers"] = p.ch.abs().values; e["pid"] = p.patient_id.values; e["cohort"] = p.cohort.values; e["has"] = p.r_hist.notna().values
        for c in ["MU", "UCSD"]:
            g = e[(e.cohort == c) & e.has].groupby("pid").mean(numeric_only=True)
            out.append({"rep": rep, "cohort": c, "n_pat": len(g), **{k + "-base": float((g[k] - g.base).mean()) for k in kinds}, "base-pers": float((g.base - g.pers).mean())})
        if rep == 0:
            print("coefficients, rep 0, per fold:", {k: np.round(np.mean(v, axis=0), 2).tolist() for k, v in coefs.items()})
            g = e[e.has & e.cohort.isin(["MU", "UCSD"])].groupby("pid").mean(numeric_only=True)
            for k in kinds:
                print(f"   rep 0, patients with history (n={len(g)}): {k} minus base {float((g[k] - g.base).mean()):+.4f} CI {ci((g[k] - g.base).values)}")
    o = pd.DataFrame(out)
    print(o.groupby("cohort").mean(numeric_only=True).drop(columns="rep").round(4).to_string())


if __name__ == "__main__":
    main()
