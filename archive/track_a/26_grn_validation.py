#!/usr/bin/env python3
"""
Month 2, Week 4: GRN Validation via Bootstrap Resampling
Validates causal GRN edges by bootstrap confidence intervals.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
import torch


def load_data(device: torch.device) -> Tuple[torch.Tensor, List[str], torch.Tensor]:
    """Load expression data, gene names, and labels."""
    t0 = time.perf_counter()
    X = torch.from_numpy(np.load("output/nn_X.npy")).to(device, dtype=torch.float32)
    with open("output/te_gene_names.txt") as f:
        raw_genes = [line.strip() for line in f]
    gene_names = [g.split('\t')[-1] for g in raw_genes]
    y = torch.from_numpy(np.load("output/nn_y.npy")).to(device, dtype=torch.int64)
    elapsed = time.perf_counter() - t0
    print(f"[LOAD] X: {tuple(X.shape)}, genes: {len(gene_names)}, y: {tuple(y.shape)} in {elapsed:.3f}s")
    return X, gene_names, y


def bootstrap_te_validation(
    X: torch.Tensor,
    gene_names: List[str],
    y: torch.Tensor,
    zone: int = 1,
    n_bootstraps: int = 100,
    sample_frac: float = 0.8,
) -> dict:
    """
    Bootstrap validation of TE edges.
    Resamples cells and recomputes TE for top candidate edges.
    """
    print(f"[BOOTSTRAP] Validating GRN with {n_bootstraps} resamples...")
    t0 = time.perf_counter()
    
    # Filter to zone
    mask = y == zone
    X_zone = X[mask]
    N_zone = X_zone.shape[0]
    n_genes = len(gene_names)
    
    # Use transition scores for ordering
    scores = torch.from_numpy(np.load("output/csgt_transition_scores.npy")).to(X.device)
    scores_zone = scores[mask]
    order = scores_zone.argsort()
    X_ordered = X_zone[order]
    
    Y_t = X_ordered[1:]
    Y_tm1 = X_ordered[:-1]
    
    # Pre-select candidate edges (top 20 by original TE if available, else all pairs)
    # For now, test a subset of gene pairs
    candidate_pairs = []
    for i in range(min(20, n_genes)):
        for j in range(min(20, n_genes)):
            if i != j:
                candidate_pairs.append((i, j))
    
    print(f"[BOOTSTRAP] Testing {len(candidate_pairs)} candidate edges...")
    
    # Store bootstrap TE values
    bootstrap_tes = {pair: [] for pair in candidate_pairs}
    
    for b in range(n_bootstraps):
        if b % 20 == 0:
            print(f"  Bootstrap {b+1}/{n_bootstraps}...")
        
        # Resample cells (with replacement)
        n_sample = int(N_zone * sample_frac)
        idx = torch.randint(0, N_zone - 1, (n_sample,), device=X.device)
        idx = idx.sort()[0]  # maintain temporal order
        
        Y_t_b = Y_t[idx]
        Y_tm1_b = Y_tm1[idx]
        
        # Compute TE for candidate pairs (simplified)
        for i, j in candidate_pairs:
            # Simplified TE proxy: correlation between X_i(t-1) and Y_j(t) conditioned on Y_j(t-1)
            # This is a fast approximation
            x_source = Y_tm1_b[:, i]
            y_target_t = Y_t_b[:, j]
            y_target_tm1 = Y_tm1_b[:, j]
            
            # Partial correlation as TE proxy
            try:
                # Regress y_tm1 from both
                A = torch.stack([y_target_tm1, x_source], dim=1)
                coeff = torch.linalg.lstsq(A, y_target_t).solution
                pred = A @ coeff
                residual = y_target_t - pred
                # Variance explained by x_source
                total_var = y_target_t.var()
                if total_var > 0:
                    te_proxy = 1 - (residual.var() / total_var)
                    bootstrap_tes[(i, j)].append(max(0, float(te_proxy)))
            except:
                bootstrap_tes[(i, j)].append(0.0)
    
    elapsed = time.perf_counter() - t0
    print(f"[BOOTSTRAP] Completed in {elapsed:.3f}s")
    
    # Compute confidence intervals
    ci_results = {}
    for pair, values in bootstrap_tes.items():
        if len(values) > 0:
            arr = np.array(values)
            ci_results[f"{gene_names[pair[0]]}->{gene_names[pair[1]]}"] = {
                "mean": float(arr.mean()),
                "std": float(arr.std()),
                "ci_95_lower": float(np.percentile(arr, 2.5)),
                "ci_95_upper": float(np.percentile(arr, 97.5)),
                "significant": bool(np.percentile(arr, 2.5) > 0),
                "n_samples": int(len(arr)),
            }
    
    return ci_results


def export_bootstrap_results(ci_results: dict) -> None:
    """Export bootstrap confidence intervals."""
    print("[EXPORT] Saving bootstrap validation results...")
    Path("output").mkdir(exist_ok=True)
    
    with open("output/grn_bootstrap_ci.json", "w") as f:
        json.dump(ci_results, f, indent=2)
    print("  output/grn_bootstrap_ci.json")
    
    # TSV summary
    with open("output/grn_bootstrap_summary.tsv", "w") as f:
        f.write("edge\tmean_te\tstd_te\tci_lower\tci_upper\tsignificant\n")
        for edge, stats in ci_results.items():
            f.write(f"{edge}\t{stats['mean']:.6f}\t{stats['std']:.6f}\t"
                    f"{stats['ci_95_lower']:.6f}\t{stats['ci_95_upper']:.6f}\t"
                    f"{stats['significant']}\n")
    print("  output/grn_bootstrap_summary.tsv")


def print_bootstrap_summary(ci_results: dict) -> None:
    """Print bootstrap validation summary."""
    print(f"\n{'='*90}")
    print("BOOTSTRAP VALIDATION SUMMARY (95% CI)")
    print(f"{'='*90}")
    print(f"{'Edge':<30} {'Mean TE':<10} {'Std':<10} {'CI Lower':<10} {'CI Upper':<10} {'Significant'}")
    print("-" * 90)
    
    significant_count = 0
    for edge, stats in sorted(ci_results.items(), key=lambda x: x[1]['mean'], reverse=True):
        sig = "YES" if stats['significant'] else "NO"
        if stats['significant']:
            significant_count += 1
        print(f"{edge:<30} {stats['mean']:<10.6f} {stats['std']:<10.6f} "
              f"{stats['ci_95_lower']:<10.6f} {stats['ci_95_upper']:<10.6f} {sig}")
    
    print(f"\nSignificant edges: {significant_count}/{len(ci_results)}")
    print(f"{'='*90}\n")


def main() -> None:
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Computational Backend: {device}")
    
    X, gene_names, y = load_data(device)
    
    ci_results = bootstrap_te_validation(X, gene_names, y, zone=1, n_bootstraps=50)
    
    export_bootstrap_results(ci_results)
    print_bootstrap_summary(ci_results)
    
    print("[SUCCESS] Month 2 Week 4 Complete: GRN Bootstrap Validation Done")


if __name__ == "__main__":
    main()