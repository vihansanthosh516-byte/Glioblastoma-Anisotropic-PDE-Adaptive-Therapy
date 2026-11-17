"""Per-nucleus program scores for the two fixed procedures.

Procedure R (reference): Tirosh et al. 2016 / Neftel et al. 2019 control-gene score.
    E = log2(CPM/10 + 1); relative expression Er = E - c (c = per-gene mean over the reference pool);
    SC(cell) = mean_{g in module} Er - mean_{g in control} Er, controls = for each module gene, 100 genes
    drawn from the same of 30 equal-size bins of aggregate expression log2(mean CPM + 1) over the pool;
    the control set is the unique union minus module genes.
Procedure A (alternative): UCell (Andreatta & Carmona 2021) Mann-Whitney U score on within-nucleus ranks,
    ranks above max_rank set to max_rank + 1, ties averaged: score = 1 - U / (n * max_rank),
    U = sum(rank) - n(n+1)/2.
MES-like per nucleus (Neftel 2019) = max(score MESlike1, score MESlike2), separately for R and A.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.stats import rankdata

N_BINS = 30
N_CTRL = 100
UCELL_MAX_RANK = 1500


def log_cpm(counts: sp.spmatrix) -> sp.csc_matrix:
    """genes x cells counts -> log2(CPM/10 + 1), sparse (zeros stay zero)."""
    m = sp.csc_matrix(counts, dtype=np.float64, copy=True)
    lib = np.asarray(m.sum(axis=0)).ravel()
    if np.any(lib <= 0):
        raise ValueError("nucleus with zero total counts")
    m.data = np.log2(m.data / np.repeat(lib, np.diff(m.indptr)) * 1e6 / 10.0 + 1.0)
    return m


def aggregate_expression(counts: sp.spmatrix) -> np.ndarray:
    m = sp.csc_matrix(counts, dtype=np.float64)
    lib = np.asarray(m.sum(axis=0)).ravel()
    cpm = m @ sp.diags(1e6 / lib)
    return np.log2(np.asarray(cpm.mean(axis=1)).ravel() + 1.0)


def expression_bins(agg: np.ndarray, n_bins: int = N_BINS) -> np.ndarray:
    """Equal-size bins on the rank of aggregate expression (ties broken by gene order)."""
    order = np.argsort(agg, kind="stable")
    bins = np.empty(len(agg), dtype=int)
    bins[order] = np.arange(len(agg)) * n_bins // len(agg)
    return bins


def control_genes(module_idx: np.ndarray, bins: np.ndarray, rng: np.random.Generator, n_ctrl: int = N_CTRL) -> np.ndarray:
    module_idx = np.asarray(module_idx)
    chosen = set()
    for g in module_idx:
        pool = np.flatnonzero(bins == bins[g])
        chosen.update(rng.choice(pool, size=min(n_ctrl, len(pool)), replace=False).tolist())
    chosen -= set(module_idx.tolist())
    return np.array(sorted(chosen), dtype=int)


def tirosh_score(logexpr: sp.spmatrix, center: np.ndarray, module_idx: np.ndarray, ctrl_idx: np.ndarray) -> np.ndarray:
    m = sp.csr_matrix(logexpr)
    mod = np.asarray(m[module_idx].mean(axis=0)).ravel() - center[module_idx].mean()
    ctl = np.asarray(m[ctrl_idx].mean(axis=0)).ravel() - center[ctrl_idx].mean()
    return mod - ctl


def ucell_score(counts: sp.spmatrix, module_idx: np.ndarray, max_rank: int = UCELL_MAX_RANK) -> np.ndarray:
    m = sp.csc_matrix(counts)
    n_genes = m.shape[0]
    module_idx = np.asarray(module_idx)
    n = len(module_idx)
    in_mod = np.zeros(n_genes, dtype=bool)
    in_mod[module_idx] = True
    pos = np.full(n_genes, -1)
    pos[module_idx] = np.arange(n)
    out = np.empty(m.shape[1])
    u_max = n * max_rank
    for j in range(m.shape[1]):
        lo, hi = m.indptr[j], m.indptr[j + 1]
        idx, val = m.indices[lo:hi], m.data[lo:hi]
        nnz = hi - lo
        ranks_nz = rankdata(-val, method="average")
        zero_rank = (nnz + 1 + n_genes) / 2.0   # average rank shared by all zero genes
        r = np.full(n, zero_rank, dtype=float)
        hit = in_mod[idx]
        r[pos[idx[hit]]] = ranks_nz[hit]
        r = np.minimum(r, max_rank + 1)
        u = r.sum() - n * (n + 1) / 2.0
        out[j] = 1.0 - u / u_max
    return out


def mes_like(scores_by_module: dict[str, np.ndarray]) -> np.ndarray:
    return np.maximum(scores_by_module["MESlike1"], scores_by_module["MESlike2"])
