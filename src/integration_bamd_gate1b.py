"""BAMD Gate 1b: does biopsy-site MRI predict biopsy RNA pathway scores in held-out patients?

Protocol: paper/integration_protocol_bamd_gate1b_2026-11-12.md. Data: Hu et al. 2023 public
table (data/external/hu2023/MultiReg_imaging_pathclass.csv). No survival data.

  python -I src/integration_bamd_gate1b.py [--perms 2000]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA = PROJECT_ROOT / "data" / "external" / "hu2023" / "MultiReg_imaging_pathclass.csv"
OUT = PROJECT_ROOT / "output" / "integration_bamd"
MRI = ["T2", "RCBV", "STDRCBV", "MD", "FA", "SPGRC", "EPI"]
TARGETS = ["GPM", "MTC", "NEU", "PPR"]
ALPHAS = (0.1, 1.0, 10.0, 100.0)
RHO_MIN, LEVEL, SEED = 0.25, 0.05, 20261116


def load(path=DATA):
    d = pd.read_csv(path)
    d["CE"] = (d["MRI"] == "CE").astype(float)
    keep = d[MRI + TARGETS + ["CE", "PatientID.x"]].dropna()
    return keep.reset_index(drop=True), len(d) - len(keep)


def ridge(Xtr, Ytr, Xte, a):
    mu, sd = Xtr.mean(0), Xtr.std(0)
    sd = np.where(sd > 0, sd, 1)
    A, B = (Xtr - mu) / sd, (Xte - mu) / sd
    ym = Ytr.mean(0)
    W = np.linalg.solve(A.T @ A + a * np.eye(A.shape[1]), A.T @ (Ytr - ym))
    return B @ W + ym


def lopo_predict(X, Y, groups):
    pred = np.zeros_like(Y)
    pats = np.unique(groups)
    for p in pats:
        te, tr = groups == p, groups != p
        # inner LOPO to choose alpha
        inner = np.unique(groups[tr])
        best, best_e = ALPHAS[0], np.inf
        for a in ALPHAS:
            e = 0.0
            for q in inner:
                ite, itr = tr & (groups == q), tr & (groups != q)
                e += ((ridge(X[itr], Y[itr], X[ite], a) - Y[ite]) ** 2).sum()
            if e < best_e:
                best, best_e = a, e
        pred[te] = ridge(X[tr], Y[tr], X[te], best)
    return pred


def rhos(pred, Y):
    return np.array([spearmanr(pred[:, j], Y[:, j]).statistic for j in range(Y.shape[1])])


def holm(p):
    p = np.asarray(p)
    adj, run = np.empty_like(p), 0.0
    for r, i in enumerate(np.argsort(p)):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        adj[i] = run
    return adj


def gate(X, Y, groups, strata, perms, rng):
    obs = rhos(lopo_predict(X, Y, groups), Y)
    null = np.empty((perms, Y.shape[1]))
    for k in range(perms):
        Yp = Y.copy()
        for s in np.unique(strata):
            idx = np.nonzero(strata == s)[0]
            Yp[idx] = Y[rng.permutation(idx)]
        null[k] = rhos(lopo_predict(X, Yp, groups), Yp)
    p = (1 + (null >= obs).sum(0)) / (1 + perms)
    return obs, p, holm(p)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=2000)
    args = ap.parse_args(argv)
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    d, dropped = load()
    groups, strata = d["PatientID.x"].astype(str).to_numpy(), d["CE"].to_numpy()
    Y = d[TARGETS].to_numpy(float)
    res = {"n_biopsies": len(d), "n_patients": int(len(np.unique(groups))), "dropped_missing": int(dropped),
           "data_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
           "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "perms": args.perms}
    for name, cols in (("primary_with_CE", MRI + ["CE"]), ("secondary_without_CE", MRI)):
        obs, p, ph = gate(d[cols].to_numpy(float), Y, groups, strata, args.perms, rng)
        res[name] = {"rho": dict(zip(TARGETS, obs.tolist())), "p": dict(zip(TARGETS, p.tolist())),
                     "holm_p": dict(zip(TARGETS, ph.tolist())),
                     "passed": [t for t, r, q in zip(TARGETS, obs, ph) if r >= RHO_MIN and q < LEVEL]}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "gate1b.json").write_text(json.dumps(res, indent=2), encoding="utf-8")
    print(json.dumps(res, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
