"""
Track A composition-aware core/periphery audit (Claude, plan trackA_feasibility_findings_2026-10-25.md, option B).
================================================================================================================
Pure functions on arrays (no file I/O). Leave-one-donor-out logistic models:

  celltype_only   one-hot cell type -> core vs periphery   (pure composition)
  genes           top-variance genes (chosen on training donors only)
  genes_celltype  genes + one-hot cell type
  within_type     genes-only model trained and scored inside one cell type, so the cell mix
                  cannot help; the question is whether same-type cells differ by region

Each model is scored per held-out donor with ROC AUC. A label-shuffle control permutes labels
within each training donor.
"""
from __future__ import annotations

import itertools

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


def top_var_genes(X, n):
    m = X.mean(0)
    v = (X * X).mean(0) - m * m
    return np.argsort(-v)[:n]


def one_hot(types, levels):
    return (types[:, None] == np.asarray(levels)[None, :]).astype(np.float32)


def _fit_score(Xtr, ytr, Xte, C=0.1, seed=0):
    mu, sd = Xtr.mean(0), Xtr.std(0) + 1e-6
    clf = LogisticRegression(C=C, max_iter=2000, random_state=seed)
    clf.fit((Xtr - mu) / sd, ytr)
    return clf.predict_proba((Xte - mu) / sd)[:, 1]


def _auc(y, s):
    return float(roc_auc_score(y, s)) if 0 < y.sum() < y.size else float("nan")


def lodo_audit(X, y, donor, ctype, n_genes=2000, min_cells=50, within_types=(), shuffle_seed=None, C=0.1):
    """Leave-one-donor-out audit. X dense cells x genes; y 1 = periphery; returns per-donor AUCs."""
    levels = sorted(set(ctype))
    oh = one_hot(ctype, levels)
    out = {}
    rng = np.random.default_rng(shuffle_seed) if shuffle_seed is not None else None
    for d in sorted(set(donor)):
        tr, te = donor != d, donor == d
        ytr = y[tr].copy()
        if rng is not None:                      # shuffle labels within each training donor
            for dd in set(donor[tr]):
                m = donor[tr] == dd
                ytr[m] = rng.permutation(ytr[m])
        g = top_var_genes(X[tr], n_genes)
        row = {"n_test": int(te.sum())}
        row["celltype_only"] = _auc(y[te], _fit_score(oh[tr], ytr, oh[te], C))
        row["genes"] = _auc(y[te], _fit_score(X[tr][:, g], ytr, X[te][:, g], C))
        Xg = np.hstack([X[:, g], oh])
        row["genes_celltype"] = _auc(y[te], _fit_score(Xg[tr], ytr, Xg[te], C))
        for t in within_types:
            mt_tr, mt_te = tr & (ctype == t), te & (ctype == t)
            ok = (min(np.sum(ytr[ctype[tr] == t] == 0), np.sum(ytr[ctype[tr] == t] == 1)) >= min_cells
                  and min(np.sum(y[mt_te] == 0), np.sum(y[mt_te] == 1)) >= min_cells)
            if not ok:
                row[f"within_{t}"] = float("nan")
                continue
            gt = top_var_genes(X[mt_tr], n_genes)
            row[f"within_{t}"] = _auc(y[mt_te], _fit_score(X[mt_tr][:, gt], ytr[ctype[tr] == t],
                                                           X[mt_te][:, gt], C))
        out[d] = row
    return out


def random_split_auc(X, y, n_genes=2000, seed=0, C=0.1):
    """Optimistic reference: 5-fold random cell split (donors leak across folds)."""
    aucs = []
    for tr, te in StratifiedKFold(5, shuffle=True, random_state=seed).split(X, y):
        g = top_var_genes(X[tr], n_genes)
        aucs.append(_auc(y[te], _fit_score(X[tr][:, g], y[tr], X[te][:, g], C)))
    return float(np.mean(aucs))


def sign_flip_p(diffs):
    """Exact one-sided sign-flip p-value for mean(diffs) > 0 (all 2^n sign patterns)."""
    d = np.asarray([x for x in diffs if np.isfinite(x)], float)
    if d.size == 0:
        return float("nan")
    obs = d.mean()
    flips = np.array(list(itertools.product((1, -1), repeat=d.size)))
    return float(np.mean((flips * d).mean(1) >= obs - 1e-12))
