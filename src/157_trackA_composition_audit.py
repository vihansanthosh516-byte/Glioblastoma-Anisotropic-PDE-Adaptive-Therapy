#!/usr/bin/env python3
"""
Script 157: Track A composition-aware core/periphery audit, leave-one-donor-out.
===============================================================================
Question: how much of core-vs-periphery prediction comes only from the cell-type mix, and do
cells of the SAME type differ between core and periphery in held-out donors?

Data: output/01_filtered_three_class.h5ad (UCSC multiomic-gbm atlas; X holds NORMALIZED values,
not raw counts). The 7 GBM donors with both Core and Periphery cells (script 151); SNU33 and
SNU33citeseq merged. Up to 3,000 cells per donor and region (seeded subsample). Cell types are
the atlas authors' `Major.celltype` labels, not re-derived here.

Exploratory: 7 donors, so per-donor results and exact sign-flip p-values (minimum 1/128).
Output stays LOCAL (UCSC-derived; public repo): output/trackA_composition_audit/.

Run:  python src/157_trackA_composition_audit.py
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackA_composition_audit as ca  # noqa: E402

H5AD = PROJECT_ROOT / "output" / "01_filtered_three_class.h5ad"
OUT_DIR = PROJECT_ROOT / "output" / "trackA_composition_audit"
DONORS = {"SNU27": "SNU27", "SNU32": "SNU32", "SNU33": "SNU33", "SNU33citeseq": "SNU33",
          "SNU34citeseq": "SNU34citeseq", "SNU35": "SNU35", "SNU46": "SNU46", "SNU51": "SNU51"}
WITHIN = ("microglia", "oligodendrocyte", "macrophage", "granulocyte", "tumor cell")


def _cat(obs, key):
    """Decode an AnnData categorical column; code -1 becomes explicit ``None``."""
    g = obs[key]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
    codes = np.asarray(g["codes"][:])
    if codes.ndim != 1 or not np.issubdtype(codes.dtype, np.integer):
        raise ValueError(f"categorical codes for {key!r} must be a one-dimensional integer array")
    invalid = (codes < -1) | (codes >= len(cats))
    if invalid.any():
        bad = codes[invalid]
        raise ValueError(f"categorical codes for {key!r} contain out-of-range values: {bad[:5].tolist()}")
    values = np.full(codes.shape, None, dtype=object)
    valid = codes >= 0  # AnnData uses -1 as the standard missing-category code.
    values[valid] = cats[codes[valid]]
    return values


def load(max_per_group, seed):
    with h5py.File(H5AD, "r") as f:
        obs = f["obs"]
        sample, cls, ctype = _cat(obs, "ID1"), _cat(obs, "class"), _cat(obs, "Major.celltype")
        rng = np.random.default_rng(seed)
        keep = []
        donor_all = np.array([DONORS.get(s, "") for s in sample], dtype=object)
        missing_sample_rows = sum(value is None for value in sample)
        missing_class_rows = sum(value is None for value in cls)
        missing_celltype_in_analysis_rows = 0
        for d in sorted(set(DONORS.values())):
            for c in ("Core", "Periphery"):
                idx = np.flatnonzero((donor_all == d) & (cls == c))
                valid_celltype = np.array([value is not None for value in ctype[idx]], dtype=bool)
                missing_celltype_in_analysis_rows += int((~valid_celltype).sum())
                idx = idx[valid_celltype]
                if idx.size > max_per_group:
                    idx = np.sort(rng.choice(idx, max_per_group, replace=False))
                keep.append(idx)
        rows = np.sort(np.concatenate(keep))
        if rows.size == 0:
            raise ValueError("no Core/Periphery rows with known donor and cell type are eligible")
        X = f["X"]
        shape = tuple(X.attrs["shape"])
        M = sp.csr_matrix((X["data"][:], X["indices"][:], X["indptr"][:]), shape=shape)[rows]
        gi = f["var"]["_index"]
        raw = gi["values"][:] if isinstance(gi, h5py.Group) else gi[:]     # nullable-string encoding
        genes = np.array([g.decode() if isinstance(g, bytes) else g for g in raw])
    diagnostics = {"missing_sample_id_rows": missing_sample_rows,
                   "missing_class_label_rows": missing_class_rows,
                   "eligible_core_periphery_rows_excluded_missing_celltype": missing_celltype_in_analysis_rows}
    return (M.toarray().astype(np.float32), (cls[rows] == "Periphery").astype(int),
            donor_all[rows], ctype[rows], genes, diagnostics)


# Region and library are confounded (core and periphery are separate libraries per donor), so a
# within-type signal could be handling/dissociation stress rather than biology. Flag such genes.
STRESS = {"FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "IER2", "IER3", "ATF3", "HSPA1A", "HSPA1B",
          "HSPA8", "HSPB1", "DNAJB1", "HSP90AA1", "DUSP1", "ZFP36", "NR4A1", "KLF2", "KLF4", "KLF6"}
MG_HOMEOSTATIC = {"P2RY12", "TMEM119", "CX3CR1", "SALL1", "CSF1R", "SELPLG", "GPR34", "OLFML3"}


def top_genes(X, y, genes, n_genes, n_top=25):
    """Descriptive all-donor logistic coefficient ranking; not differential expression."""
    from sklearn.linear_model import LogisticRegression
    g = ca.top_var_genes(X, n_genes)
    Z = (X[:, g] - X[:, g].mean(0)) / (X[:, g].std(0) + 1e-6)
    coef = LogisticRegression(C=0.1, max_iter=2000).fit(Z, y).coef_[0]
    order = np.argsort(coef)
    up = [str(genes[g[i]]) for i in order[::-1][:n_top]]
    down = [str(genes[g[i]]) for i in order[:n_top]]

    def flag(lst):
        return {"stress": [x for x in lst if x in STRESS],
                "ribo_mito": [x for x in lst if x.startswith(("RPL", "RPS", "MT-"))],
                "homeostatic": [x for x in lst if x in MG_HOMEOSTATIC]}

    return {"positive_logistic_coefficient_genes": up,
            "negative_logistic_coefficient_genes": down,
            "positive_coefficient_gene_flags": flag(up),
            "negative_coefficient_gene_flags": flag(down)}


def summarize(res, keys):
    out = {}
    for k in keys:
        v = np.array([r[k] for r in res.values()], float)
        out[k] = {"mean": float(np.nanmean(v)) if np.isfinite(v).any() else None,
                  "n_donors": int(np.isfinite(v).sum()),
                  "per_donor": {d: (None if not np.isfinite(r[k]) else r[k]) for d, r in res.items()}}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-per-group", type=int, default=3000)
    ap.add_argument("--n-genes", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=20261025)
    ap.add_argument("--exclude-technical", action="store_true")
    args = ap.parse_args()
    out_dir = OUT_DIR / ("exclude_technical" if args.exclude_technical else "")
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    X, y, donor, ctype, genes, input_missingness = load(args.max_per_group, args.seed)
    t_load = time.time() - t0
    n_excluded = 0
    if args.exclude_technical:
        # drop genes that often reflect handling or contamination, not cell state: immediate-early
        # stress, ribosomal, mitochondrial, hemoglobin, nuclear lncRNAs (sensitivity analysis)
        tech = np.array([g in STRESS or g.startswith(("RPL", "RPS", "MT-", "MTRNR", "HBA", "HBB", "HBD"))
                         or g in {"MALAT1", "NEAT1"} for g in genes])
        X = X[:, ~tech]
        genes = genes[~tech]
        n_excluded = int(tech.sum())
    mg = ctype == "microglia"
    mg_genes = top_genes(X[mg], y[mg], genes, args.n_genes)

    keys = ["celltype_only", "genes", "genes_celltype"] + [f"within_{t}" for t in WITHIN]
    res = ca.lodo_audit(X, y, donor, ctype, n_genes=args.n_genes, within_types=WITHIN)
    shuf = ca.lodo_audit(X, y, donor, ctype, n_genes=args.n_genes, within_types=("microglia",), shuffle_seed=7)
    rand = ca.random_split_auc(X, y, n_genes=args.n_genes, seed=7)

    def diffs(a, b):
        return [res[d][a] - (res[d][b] if isinstance(b, str) else b) for d in res]

    s = summarize(res, keys)
    g, c = s["genes"]["mean"], s["celltype_only"]["mean"]
    tests = {
        "genes_minus_celltype_only": {"mean": float(np.mean(diffs("genes", "celltype_only"))),
                                      "p_one_sided": ca.sign_flip_p(diffs("genes", "celltype_only"))},
        "genes_celltype_minus_celltype_only": {
            "mean": float(np.mean(diffs("genes_celltype", "celltype_only"))),
            "p_one_sided": ca.sign_flip_p(diffs("genes_celltype", "celltype_only"))},
    }
    for t in WITHIN:
        dd = diffs(f"within_{t}", 0.5)
        tests[f"within_{t}_minus_0.5"] = {"mean": float(np.nanmean(dd)) if np.isfinite(dd).any() else None,
                                          "p_one_sided": ca.sign_flip_p(dd)}

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    counts = {d: {"n_core": int(((donor == d) & (y == 0)).sum()), "n_periphery": int(((donor == d) & (y == 1)).sum()),
                  "by_type_periphery": {t: int(((donor == d) & (y == 1) & (ctype == t)).sum()) for t in WITHIN},
                  "by_type_core": {t: int(((donor == d) & (y == 0) & (ctype == t)).sum()) for t in WITHIN}}
              for d in sorted(set(donor))}
    result = {
        "script": "157_trackA_composition_audit",
        "label": ("EXPLORATORY, 7 donors, leave-one-donor-out. X is normalized (not raw counts). Cell types are "
                  "atlas labels. Not a diagnostic or clinical result."),
        "git_commit_at_run": commit, "input": str(H5AD.relative_to(PROJECT_ROOT)),
        "settings": {"max_per_group": args.max_per_group, "n_genes": args.n_genes, "seed": args.seed,
                     "model": "L2 logistic regression, C = 0.1, standardized features",
                     "within_min_cells_per_side": 50,
                     "exclude_technical": args.exclude_technical, "n_genes_excluded": n_excluded},
        "cells_used": counts,
        "input_missingness": input_missingness,
        "auc": s,
        "celltype_only_to_gene_excess_auc_ratio": (c - 0.5) / (g - 0.5) if g and g > 0.5 else None,
        "tests": tests,
        "shuffle_control": summarize(shuf, ["celltype_only", "genes", "within_microglia"]),
        "random_cell_split_auc_genes": rand,
        "donor_vs_random_optimism_gap": rand - g,
        "microglia_top_genes_all_donors_descriptive": mg_genes,
        "runtime_s": {"load": t_load, "total": time.time() - t0},
    }
    (out_dir / "audit_results.json").write_text(json.dumps(result, indent=2, default=float, allow_nan=False))
    print(json.dumps({k: result[k] for k in ("celltype_only_to_gene_excess_auc_ratio", "tests",
                                             "random_cell_split_auc_genes")}, indent=1, default=float, allow_nan=False))
    print({k: v["mean"] for k, v in s.items()})


if __name__ == "__main__":
    main()
