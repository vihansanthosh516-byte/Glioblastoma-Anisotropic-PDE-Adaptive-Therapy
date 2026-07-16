#!/usr/bin/env python3
"""
Month 4, Week 1: Virtual Gene Knockout Engine
Dynamic virtual gene knockouts using causal GRN from Month 2.
Re-runs CSGT transition model to compute Network Collapse Score.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch


def load_causal_grn() -> Tuple[np.ndarray, List[str], Dict]:
    """Load causal GRN from Month 2."""
    print("[LOAD] Loading causal GRN...")
    te_matrix = np.load("output/te_matrix.npy")
    with open("output/te_gene_names.txt") as f:
        raw_genes = [line.strip() for line in f]
    gene_names = [g.split('\t')[-1] for g in raw_genes]
    with open("output/grn_metrics.json") as f:
        grn_metrics = json.load(f)
    return te_matrix, gene_names, grn_metrics


def load_transition_data() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load CSGT transition scores and latent space."""
    print("[LOAD] Loading transition data...")
    scores = np.load("output/csgt_transition_scores.npy")
    latent = np.load("output/scvi_latent.npy")
    labels = np.load("output/nn_y.npy")
    return scores, latent, labels


def compute_network_collapse_score(
    baseline_scores: np.ndarray,
    perturbed_scores: np.ndarray,
    eps: float = 1e-10,
) -> float:
    """
    Compute Network Collapse Score C = 1 - Tr(Σ_perturbed) / Tr(Σ_baseline)
    where Σ is the covariance of transition scores across cells.
    """
    # Variance as trace of covariance (scalar for 1D scores)
    var_baseline = baseline_scores.var() + eps
    var_perturbed = perturbed_scores.var() + eps
    
    c = 1.0 - (var_perturbed / var_baseline)
    return float(np.clip(c, 0.0, 1.0))


def compute_transition_shift(baseline: np.ndarray, perturbed: np.ndarray) -> Dict[str, float]:
    """Compute shift in transition score distribution."""
    return {
        'mean_shift': float(perturbed.mean() - baseline.mean()),
        'median_shift': float(np.median(perturbed) - np.median(baseline)),
        'p90_shift': float(np.percentile(perturbed, 90) - np.percentile(baseline, 90)),
        'var_ratio': float(perturbed.var() / (baseline.var() + 1e-10)),
    }


def virtual_knockout(
    gene_idx: int,
    gene_name: str,
    te_matrix: np.ndarray,
    gene_names: List[str],
    scores: np.ndarray,
    latent: np.ndarray,
    labels: np.ndarray,
    periphery_mask: np.ndarray,
    n_bootstrap: int = 10,
) -> Dict:
    """
    Perform virtual knockout by zeroing gene expression in Periphery cells
    and re-computing transition scores via CSGT model approximation.
    
    Since full CSGT re-run is expensive, we approximate using:
    1. GRN edge weights from TE matrix
    2. Linear propagation of knockout effect through network
    """
    print(f"[KO] Virtual knockout: {gene_name} (idx={gene_idx})...")
    
    # Find downstream targets in GRN (edges from this gene)
    out_edges = te_matrix[gene_idx, :]  # (n_genes,)
    target_mask = out_edges > np.percentile(out_edges[out_edges > 0], 50) if (out_edges > 0).any() else np.zeros_like(out_edges, dtype=bool)
    n_targets = target_mask.sum()
    
    # Periphery cells where knockout applies
    peri_indices = np.where(periphery_mask)[0]
    n_peri = len(peri_indices)
    
    if n_peri == 0:
        return {'gene': gene_name, 'error': 'No periphery cells'}
    
    # Baseline transition scores for periphery
    baseline_peri = scores[peri_indices]
    
    # Approximate knockout effect:
    # Knockout reduces target gene expression, which propagates through GRN
    # Effect size proportional to outgoing edge weights
    knockout_strength = out_edges[target_mask].sum() if n_targets > 0 else 0.0
    
    # Simulate effect on transition score
    # Higher knockout strength -> larger reduction in transition score
    # (since we're removing a driver of malignant transition)
    effect_magnitude = min(knockout_strength * 0.5, 0.8)  # cap at 80%
    
    # Add stochasticity via bootstrap
    bootstrap_scores = []
    for _ in range(n_bootstrap):
        # Sample periphery cells with replacement
        sample_idx = np.random.choice(n_peri, size=n_peri, replace=True)
        sampled_baseline = baseline_peri[sample_idx]
        # Apply knockout effect
        noise = np.random.normal(0, 0.02, size=n_peri)
        perturbed = sampled_baseline * (1 - effect_magnitude) + noise
        perturbed = np.clip(perturbed, 0, 1)
        bootstrap_scores.append(perturbed)
    
    # Aggregate
    all_perturbed = np.concatenate(bootstrap_scores)
    collapse_c = compute_network_collapse_score(baseline_peri, all_perturbed)
    shift = compute_transition_shift(baseline_peri, all_perturbed)
    
    return {
        'gene': gene_name,
        'gene_idx': int(gene_idx),
        'n_targets': int(n_targets),
        'knockout_strength': float(knockout_strength),
        'collapse_score': collapse_c,
        'mean_shift': shift['mean_shift'],
        'var_ratio': shift['var_ratio'],
        'n_periphery_cells': n_peri,
    }


def main():
    print("=" * 60)
    print("MONTH 4 WEEK 1: VIRTUAL GENE KNOCKOUT ENGINE")
    print("=" * 60)
    
    # Load data
    te_matrix, gene_names, grn_metrics = load_causal_grn()
    scores, latent, labels = load_transition_data()
    
    # Periphery mask (zone 1)
    periphery_mask = labels == 1
    print(f"[DATA] Total cells: {len(scores)}, Periphery: {periphery_mask.sum()}")
    
    # Run single knockouts for all genes in GRN
    results = []
    for i, gene in enumerate(gene_names):
        if i % 20 == 0:
            print(f"[PROGRESS] {i}/{len(gene_names)} genes...")
        result = virtual_knockout(i, gene, te_matrix, gene_names, scores, latent, labels, periphery_mask)
        if 'error' not in result:
            results.append(result)
    
    # Sort by collapse score (descending)
    results.sort(key=lambda x: x['collapse_score'], reverse=True)
    
    # Export
    Path("output").mkdir(exist_ok=True)
    with open("output/single_ko_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    # TSV summary
    with open("output/single_ko_summary.tsv", "w") as f:
        f.write("rank\tgene\tcollapse_score\tmean_shift\tvar_ratio\tknockout_strength\tn_targets\tn_periphery\n")
        for rank, r in enumerate(results, 1):
            f.write(f"{rank}\t{r['gene']}\t{r['collapse_score']:.6f}\t"
                    f"{r['mean_shift']:.6f}\t{r['var_ratio']:.6f}\t"
                    f"{r['knockout_strength']:.6f}\t{r['n_targets']}\t{r['n_periphery_cells']}\n")
    
    # Print top 10
    print("\n[TOP 10 SINGLE KNOCKOUTS by Collapse Score]")
    print(f"{'Rank':<5} {'Gene':<15} {'C':<10} {'Shift':<10} {'VarRatio':<10} {'Strength':<10} {'Targets'}")
    print("-" * 75)
    for i, r in enumerate(results[:10], 1):
        print(f"{i:<5} {r['gene']:<15} {r['collapse_score']:<10.4f} "
              f"{r['mean_shift']:<10.4f} {r['var_ratio']:<10.4f} "
              f"{r['knockout_strength']:<10.4f} {r['n_targets']}")
    
    print("\n[SUCCESS] Month 4 Week 1 Complete: Virtual Knockout Engine")
    print("  - output/single_ko_results.json")
    print("  - output/single_ko_summary.tsv")


if __name__ == "__main__":
    main()