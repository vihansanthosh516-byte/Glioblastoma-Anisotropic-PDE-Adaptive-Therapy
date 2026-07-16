#!/usr/bin/env python3
"""
Month 4, Week 2: Combinatorial Drug Screen
Dual gene knockout screening for synergistic therapeutic combinations.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch


def load_data() -> Tuple[np.ndarray, List[str], np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load all required data."""
    print("[LOAD] Loading data...")
    te_matrix = np.load("output/te_matrix.npy")
    with open("output/te_gene_names.txt") as f:
        raw_genes = [line.strip() for line in f]
    gene_names = [g.split('\t')[-1] for g in raw_genes]
    scores = np.load("output/csgt_transition_scores.npy")
    labels = np.load("output/nn_y.npy")
    single_ko = json.load(open("output/single_ko_results.json"))
    return te_matrix, gene_names, scores, labels, single_ko


def compute_bliss_synergy(
    effect_a: float,
    effect_b: float,
    effect_ab: float,
) -> float:
    """
    Bliss independence synergy score.
    Expected combined effect: E_ab = E_a + E_b - E_a * E_b
    Synergy = Observed - Expected
    """
    expected = effect_a + effect_b - effect_a * effect_b
    synergy = effect_ab - expected
    return synergy


def compute_loewe_additivity(
    effect_a: float,
    effect_b: float,
    effect_ab: float,
) -> float:
    """
    Loewe additivity: combination index CI = D_a/D_a* + D_b/D_b*
    where D_a* is dose of A alone giving effect_ab.
    Simplified for binary knockouts.
    """
    if effect_a + effect_b == 0:
        return 0.0
    # CI < 1 = synergy, CI > 1 = antagonism
    ci = (effect_a / (effect_ab + 1e-10)) + (effect_b / (effect_ab + 1e-10))
    return 1.0 - ci  # Positive = synergy


def virtual_dual_knockout(
    gene_i: int,
    gene_j: int,
    gene_names: List[str],
    te_matrix: np.ndarray,
    scores: np.ndarray,
    labels: np.ndarray,
    single_ko_map: Dict[str, Dict],
    periphery_mask: np.ndarray,
) -> Dict:
    """
    Simulate dual knockout using single KO effects and GRN topology.
    
    Effects combine via:
    1. Direct edge between genes (from TE matrix)
    2. Shared downstream targets
    3. Bliss independence for independent pathways
    """
    gene_i_name = gene_names[gene_i]
    gene_j_name = gene_names[gene_j]
    
    # Get single KO effects
    ko_i = single_ko_map.get(gene_i_name, {})
    ko_j = single_ko_map.get(gene_j_name, {})
    
    effect_i = ko_i.get('collapse_score', 0.0)
    effect_j = ko_j.get('collapse_score', 0.0)
    
    # Direct interaction strength
    direct_interaction = max(te_matrix[gene_i, gene_j], te_matrix[gene_j, gene_i])
    
    # Shared targets in GRN
    targets_i = set(np.where(te_matrix[gene_i, :] > 0.01)[0])
    targets_j = set(np.where(te_matrix[gene_j, :] > 0.01)[0])
    shared_targets = len(targets_i & targets_j)
    
    # Combine effects
    # Base: Bliss independence
    bliss_expected = effect_i + effect_j - effect_i * effect_j
    
    # Synergy boost from direct interaction and shared targets
    synergy_boost = 0.0
    if direct_interaction > 0.02:
        synergy_boost += direct_interaction * 0.5
    if shared_targets > 2:
        synergy_boost += min(shared_targets * 0.02, 0.15)
    
    # Combined effect
    combined_effect = min(bliss_expected + synergy_boost, 0.95)
    
    # Bliss synergy score
    bliss_synergy = combined_effect - bliss_expected
    
    # Loewe synergy
    loewe_synergy = compute_loewe_additivity(effect_i, effect_j, combined_effect)
    
    return {
        'gene_a': gene_i_name,
        'gene_b': gene_j_name,
        'idx_a': int(gene_i),
        'idx_b': int(gene_j),
        'effect_a': float(effect_i),
        'effect_b': float(effect_j),
        'combined_effect': float(combined_effect),
        'bliss_synergy': float(bliss_synergy),
        'loewe_synergy': float(loewe_synergy),
        'direct_interaction': float(direct_interaction),
        'shared_targets': int(shared_targets),
    }


def main():
    print("=" * 60)
    print("MONTH 4 WEEK 2: COMBINATORIAL DRUG SCREEN")
    print("=" * 60)
    
    # Load data
    te_matrix, gene_names, scores, labels, single_ko = load_data()
    periphery_mask = labels == 1
    print(f"[DATA] Genes: {len(gene_names)}, Periphery cells: {periphery_mask.sum()}")
    print(f"[DATA] Single KO results: {len(single_ko)}")
    
    # Map single KO results
    single_ko_map = {r['gene']: r for r in single_ko}
    
    # Screen top genes by single KO collapse score
    top_genes = [r['gene'] for r in sorted(single_ko, key=lambda x: x['collapse_score'], reverse=True)[:50]]
    top_indices = [gene_names.index(g) for g in top_genes if g in gene_names]
    print(f"[SCREEN] Testing top {len(top_indices)} genes -> {len(top_indices)*(len(top_indices)-1)//2} pairs")
    
    # Run combinatorial screen
    results = []
    total_pairs = len(top_indices) * (len(top_indices) - 1) // 2
    pair_count = 0
    
    for i_idx, i in enumerate(top_indices):
        for j in top_indices[i_idx+1:]:
            pair_count += 1
            if pair_count % 500 == 0:
                print(f"[PROGRESS] {pair_count}/{total_pairs} pairs...")
            
            result = virtual_dual_knockout(i, j, gene_names, te_matrix, scores, labels, single_ko_map, periphery_mask)
            results.append(result)
    
    # Sort by combined effect
    results.sort(key=lambda x: x['combined_effect'], reverse=True)
    
    # Also sort by synergy
    results_by_bliss = sorted(results, key=lambda x: x['bliss_synergy'], reverse=True)
    results_by_loewe = sorted(results, key=lambda x: x['loewe_synergy'], reverse=True)
    
    # Export
    Path("output").mkdir(exist_ok=True)
    
    with open("output/dual_ko_results.json", "w") as f:
        json.dump({
            'by_combined_effect': results,
            'by_bliss_synergy': results_by_bliss,
            'by_loewe_synergy': results_by_loewe,
        }, f, indent=2)
    
    # TSV summary
    with open("output/dual_ko_summary.tsv", "w") as f:
        f.write("rank\tgene_a\tgene_b\tcombined_effect\teffect_a\teffect_b\tbliss_synergy\tloewe_synergy\tdirect_interaction\tshared_targets\n")
        for rank, r in enumerate(results[:100], 1):
            f.write(f"{rank}\t{r['gene_a']}\t{r['gene_b']}\t{r['combined_effect']:.6f}\t"
                    f"{r['effect_a']:.6f}\t{r['effect_b']:.6f}\t{r['bliss_synergy']:.6f}\t"
                    f"{r['loewe_synergy']:.6f}\t{r['direct_interaction']:.6f}\t{r['shared_targets']}\n")
    
    # Print top 10 by combined effect
    print("\n[TOP 10 DUAL KNOCKOUTS by Combined Effect]")
    print(f"{'Rank':<5} {'Gene A':<15} {'Gene B':<15} {'Combined':<10} {'Bliss':<10} {'Loewe':<10} {'Direct':<10} {'Shared'}")
    print("-" * 95)
    for i, r in enumerate(results[:10], 1):
        print(f"{i:<5} {r['gene_a']:<15} {r['gene_b']:<15} {r['combined_effect']:<10.4f} "
              f"{r['bliss_synergy']:<10.4f} {r['loewe_synergy']:<10.4f} "
              f"{r['direct_interaction']:<10.4f} {r['shared_targets']}")
    
    # Print top 10 by synergy
    print("\n[TOP 10 by Bliss Synergy]")
    for i, r in enumerate(results_by_bliss[:10], 1):
        print(f"  {i}. {r['gene_a']} + {r['gene_b']}: Bliss={r['bliss_synergy']:.4f}, Combined={r['combined_effect']:.4f}")
    
    print("\n[SUCCESS] Month 4 Week 2 Complete: Combinatorial Drug Screen")
    print("  - output/dual_ko_results.json")
    print("  - output/dual_ko_summary.tsv")


if __name__ == "__main__":
    main()