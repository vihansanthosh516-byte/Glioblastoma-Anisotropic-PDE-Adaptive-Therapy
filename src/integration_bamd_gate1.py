"""BAMD step 3: Gate 1 (outcome-free). Can MRI features predict RNA cell-state scores?

Protocol: paper/integration_protocol_bamd_2026-11-12.md. Reads MRI features and TCGA
U133A expression only; never reads survival. Writes output/integration_bamd/gate1.json
and, for passed states, the frozen map and residual PCA (frozen_map.json).

  python -I src/integration_bamd_gate1.py [--perms 1000]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from integration_bamd_features import feature_names  # noqa: E402

OUT = PROJECT_ROOT / "output" / "integration_bamd"
EXPR = PROJECT_ROOT / "data" / "external" / "tcga_gbm_u133a" / "HT_HG-U133A.gz"
MODULES = PROJECT_ROOT / "data" / "external" / "neftel2019" / "IDHwt.GBM.MetaModules.tsv"
STATES = {"MES": ("MESlike1", "MESlike2"), "AC": ("AClike",), "OPC": ("OPClike",), "NPC": ("NPClike1", "NPClike2")}
ALPHAS = (0.1, 1.0, 10.0, 100.0, 1000.0)
RHO_MIN, ALPHA_LEVEL, N_PCS = 0.25, 0.05, 3
SEED = 20261112


def load_modules(path=MODULES) -> dict:
    t = pd.read_csv(path, sep="\t")
    return {s: sorted({g for c in cols for g in t[c].dropna().astype(str)}) for s, cols in STATES.items()}


def pick_samples(columns) -> dict:
    """One sample per patient: primary tumor -01 first, then the lowest sample code."""
    best = {}
    for c in columns:
        pid, code = c[:12], c[13:15]
        rank = (0 if code == "01" else 1, code, c)
        if pid not in best or rank < best[pid][0]:
            best[pid] = (rank, c)
    return {pid: v[1] for pid, v in best.items()}


def state_scores(expr: pd.DataFrame, modules: dict) -> pd.DataFrame:
    """expr: genes x patients (log scale). Mean module z-score minus mean all-gene z-score."""
    z = expr.sub(expr.mean(axis=1), axis=0).div(expr.std(axis=1).replace(0, np.nan), axis=0).dropna()
    background = z.mean(axis=0)
    out = {}
    for s, genes in modules.items():
        present = [g for g in genes if g in z.index]
        if len(present) < 10:
            raise ValueError(f"{s}: only {len(present)} module genes on the array")
        out[s] = z.loc[present].mean(axis=0) - background
    return pd.DataFrame(out)


def ridge_fit(X, Y, alpha):
    mu, sd = X.mean(0), X.std(0)
    sd = np.where(sd > 0, sd, 1.0)
    Xs = (X - mu) / sd
    ym = Y.mean(0)
    W = np.linalg.solve(Xs.T @ Xs + alpha * np.eye(X.shape[1]), Xs.T @ (Y - ym))
    return {"mu": mu, "sd": sd, "W": W, "ym": ym}


def ridge_predict(m, X):
    return ((X - m["mu"]) / m["sd"]) @ m["W"] + m["ym"]


def impute(train, test):
    med = np.nanmedian(train, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    return np.where(np.isnan(train), med, train), np.where(np.isnan(test), med, test)


def folds(n, k, rng):
    idx = rng.permutation(n)
    return [idx[i::k] for i in range(k)]


def choose_alpha(X, Y, rng, k=5):
    errs = []
    fs = folds(len(X), k, rng)
    for a in ALPHAS:
        e = 0.0
        for te in fs:
            tr = np.setdiff1d(np.arange(len(X)), te)
            Xtr, Xte = impute(X[tr], X[te])
            e += ((ridge_predict(ridge_fit(Xtr, Y[tr], a), Xte) - Y[te]) ** 2).sum()
        errs.append(e)
    return ALPHAS[int(np.argmin(errs))]


def cv_rho(X, Y, rng, repeats=10, k=5) -> np.ndarray:
    """Mean over repeats of per-state Spearman rho of out-of-fold predictions."""
    rhos = []
    for _ in range(repeats):
        pred = np.zeros_like(Y)
        for te in folds(len(X), k, rng):
            tr = np.setdiff1d(np.arange(len(X)), te)
            a = choose_alpha(X[tr], Y[tr], rng)
            Xtr, Xte = impute(X[tr], X[te])
            pred[te] = ridge_predict(ridge_fit(Xtr, Y[tr], a), Xte)
        rhos.append([spearmanr(pred[:, j], Y[:, j]).statistic for j in range(Y.shape[1])])
    return np.mean(rhos, axis=0)


def holm(p):
    p = np.asarray(p)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[i]))
        adj[i] = running
    return adj


def freeze(X, Y, passed_idx, rng):
    """Map on all anchor patients + residual PCA, as plain lists for JSON."""
    med = np.nanmedian(X, axis=0)
    med = np.where(np.isnan(med), 0.0, med)
    Xi = np.where(np.isnan(X), med, X)
    a = choose_alpha(Xi, Y[:, passed_idx], rng)
    m = ridge_fit(Xi, Y[:, passed_idx], a)
    Xs = (Xi - m["mu"]) / m["sd"]
    P = Xs @ m["W"]                                    # biology-anchored scores (centered)
    B = np.linalg.lstsq(P, Xs, rcond=None)[0]          # project X on predicted states
    R = Xs - P @ B
    _, _, Vt = np.linalg.svd(R - R.mean(0), full_matrices=False)
    return {"alpha": a, "impute_median": med.tolist(), "mu": m["mu"].tolist(), "sd": m["sd"].tolist(),
            "W": m["W"].tolist(), "ym": m["ym"].tolist(), "B": B.tolist(), "R_mean": R.mean(0).tolist(),
            "pcs": Vt[:N_PCS].tolist()}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=1000)
    args = ap.parse_args(argv)
    rng = np.random.Generator(np.random.PCG64DXSM(SEED))
    feats = pd.read_csv(OUT / "features.csv")
    tcga = feats[feats["collection"] == "TCGA-GBM"].copy()
    with gzip.open(EXPR, "rt") as fh:
        expr = pd.read_csv(fh, sep="\t", index_col=0)
    samples = pick_samples(expr.columns)
    tcga = tcga[tcga["patient_id"].isin(samples)].drop_duplicates("patient_id").sort_values("patient_id")
    cols = [samples[p] for p in tcga["patient_id"]]
    S = state_scores(expr[cols], load_modules())
    X = tcga[feature_names()].to_numpy(float)
    Y = S.to_numpy(float)
    observed = cv_rho(X, Y, rng)
    null = np.array([cv_rho(X, Y[rng.permutation(len(Y))], rng, repeats=1) for _ in range(args.perms)])
    p = (1 + (null >= observed).sum(0)) / (1 + args.perms)
    p_holm = holm(p)
    passed = [s for j, s in enumerate(STATES) if observed[j] >= RHO_MIN and p_holm[j] < ALPHA_LEVEL]
    result = {"n_anchor": len(tcga), "states": list(STATES), "cv_rho": dict(zip(STATES, observed.tolist())),
              "perm_p": dict(zip(STATES, p.tolist())), "holm_p": dict(zip(STATES, p_holm.tolist())),
              "passed": passed, "rule": f"rho >= {RHO_MIN} and Holm p < {ALPHA_LEVEL}", "perms": args.perms,
              "null_note": "permutation null uses 1 CV repeat per shuffle (observed uses 10); conservative for p",
              "survival_read": False,
              "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    if passed:
        frozen = freeze(X, Y, [list(STATES).index(s) for s in passed], rng)
        frozen.update({"passed_states": passed, "feature_names": feature_names()})
        body = json.dumps(frozen).encode()
        (OUT / "frozen_map.json").write_bytes(body)
        result["frozen_map_sha256"] = hashlib.sha256(body).hexdigest()
    (OUT / "gate1.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("n_anchor", "cv_rho", "holm_p", "passed")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
