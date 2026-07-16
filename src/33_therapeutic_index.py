#!/usr/bin/env python3
"""
Month 4, Week 3: Therapeutic Index Calculation
Computes TI = C_tumor / C_healthy for all single and dual knockouts.
Identifies combinations that collapse tumor transition while sparing healthy cells.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import torch


def load_all_data() -> Dict:
    """Load all required data."""
    print("[LOAD] Loading all datasets...")
    
    # Month 2: GRN
    te_matrix = np.load("output/te_matrix.npy")
    with open("output/te_gene_names.txt") as f:
        raw_genes = [line.strip() for line in f]
    gene_names = [g.split('\t')[-1] for g in raw_genes]
    
    # Month 1: Transition scores, latent, labels
    scores = np.load("output/csgt_transition_scores.npy")
    latent = np.load("output/scvi_latent.npy")
    labels = np.load("output/nn_y.npy")
    
    # Month 4: Single and dual KO results
    single_ko = json.load(open("output/single_ko_results.json"))
    dual_ko_data = json.load(open("output/dual_ko_results.json"))
    dual_ko = dual_ko_data['by_combined_effect']
    
    # Zone masks
    healthy_mask = labels == 0
    periphery_mask = labels == 1
    core_mask = labels == 2
    
    print(f"[DATA] Cells: Healthy={healthy_mask.sum()}, Periphery={periphery_mask.sum()}, Core={core_mask.sum()}")
    print(f"[DATA] Single KO: {len(single_ko)}, Dual KO: {len(dual_ko)}")
    
    return {
        'te_matrix': te_matrix,
        'gene_names': gene_names,
        'scores': scores,
        'labels': labels,
        'single_ko': single_ko,
        'dual_ko': dual_ko,
        'healthy_mask': healthy_mask,
        'periphery_mask': periphery_mask,
        'core_mask': core_mask,
    }


def compute_zone_transition_stats(scores: np.ndarray, labels: np.ndarray) -> Dict[str, Dict]:
    """Compute transition score statistics per zone."""
    zones = {'healthy': 0, 'periphery': 1, 'core': 2}
    stats = {}
    for name, zone_id in zones.items():
        mask = labels == zone_id
        if mask.any():
            s = scores[mask]
            stats[name] = {
                'n': int(mask.sum()),
                'mean': float(s.mean()),
                'std': float(s.std()),
                'median': float(np.median(s)),
                'p90': float(np.percentile(s, 90)),
                'p95': float(np.percentile(s, 95)),
            }
        else:
            stats[name] = {'n': 0}
    return stats


def estimate_healthy_impact(
    gene_name: str,
    gene_idx: int,
    te_matrix: np.ndarray,
    gene_names: List[str],
    healthy_scores: np.ndarray,
) -> float:
    """
    Estimate impact of knockout on healthy cells.
    Healthy cells should have minimal transition activity.
    Impact = disruption of healthy GRN stability.
    """
    # Get outgoing edges from this gene
    out_edges = te_matrix[gene_idx, :]
    
    # Healthy impact proportional to edge weights to stable genes
    # Genes with high expression in healthy that would be disrupted
    # Approximate: sum of outgoing edges * baseline transition score
    baseline = healthy_scores.mean() if len(healthy_scores) > 0 else 0.0
    impact = float(out_edges.sum() * baseline * 0.1)  # scaling factor
    
    return min(impact, 1.0)


def compute_therapeutic_index(
    tumor_collapse: float,
    healthy_impact: float,
    eps: float = 1e-6,
) -> float:
    """
    Therapeutic Index: TI = C_tumor / C_healthy
    Higher TI = better (more tumor collapse, less healthy disruption)
    """
    if healthy_impact < eps:
        return float('inf') if tumor_collapse > 0 else 0.0
    return tumor_collapse / healthy_impact


def main():
    print("=" * 60)
    print("MONTH 4 WEEK 3: THERAPEUTIC INDEX CALCULATION")
    print("=" * 60)
    
    data = load_all_data()
    te_matrix = data['te_matrix']
    gene_names = data['gene_names']
    scores = data['scores']
    labels = data['labels']
    single_ko = data['single_ko']
    dual_ko = data['dual_ko']
    healthy_mask = data['healthy_mask']
    periphery_mask = data['periphery_mask']
    core_mask = data['core_mask']
    
    # Baseline stats per zone
    zone_stats = compute_zone_transition_stats(scores, labels)
    for zone, stats in zone_stats.items():
        print(f"  {zone}: n={stats['n']}, mean={stats['mean']:.4f}, std={stats['std']:.4f}")
    
    healthy_scores = scores[healthy_mask]
    tumor_scores = scores[periphery_mask | core_mask]
    
    # Map gene to index
    gene_to_idx = {name: i for i, name in enumerate(gene_names)}
    
    # --- SINGLE KO THERAPEUTIC INDEX ---
    print("\n[TI] Computing Single KO Therapeutic Indices...")
    single_ti = []
    
    for r in single_ko:
        gene = r['gene']
        idx = gene_to_idx.get(gene, -1)
        if idx == -1:
            continue
        
        tumor_c = r['collapse_score']
        healthy_impact = estimate_healthy_impact(gene, idx, te_matrix, gene_names, healthy_scores)
        ti = compute_therapeutic_index(tumor_c, healthy_impact)
        
        single_ti.append({
            'gene': gene,
            'tumor_collapse': tumor_c,
            'healthy_impact': healthy_impact,
            'therapeutic_index': ti if ti != float('inf') else 999.0,
            'is_infinite': ti == float('inf'),
        })
    
    # Sort by TI
    single_ti.sort(key=lambda x: x['therapeutic_index'], reverse=True)
    
    # --- DUAL KO THERAPEUTIC INDEX ---
    print("[TI] Computing Dual KO Therapeutic Indices...")
    dual_ti = []
    
    for r in dual_ko:
        gene_a = r['gene_a']
        gene_b = r['gene_b']
        idx_a = gene_to_idx.get(gene_a, -1)
        idx_b = gene_to_idx.get(gene_b, -1)
        
        if idx_a == -1 or idx_b == -1:
            continue
        
        tumor_c = r['combined_effect']
        
        # Healthy impact for combination (additive approximation)
        impact_a = estimate_healthy_impact(gene_a, idx_a, te_matrix, gene_names, healthy_scores)
        impact_b = estimate_healthy_impact(gene_b, idx_b, te_matrix, gene_names, healthy_scores)
        healthy_impact = min(impact_a + impact_b, 1.0)
        
        ti = compute_therapeutic_index(tumor_c, healthy_impact)
        
        dual_ti.append({
            'gene_a': gene_a,
            'gene_b': gene_b,
            'tumor_collapse': tumor_c,
            'healthy_impact': healthy_impact,
            'therapeutic_index': ti if ti != float('inf') else 999.0,
            'bliss_synergy': r['bliss_synergy'],
            'loewe_synergy': r['loewe_synergy'],
            'is_infinite': ti == float('inf'),
        })
    
    dual_ti.sort(key=lambda x: x['therapeutic_index'], reverse=True)
    
    # --- EXPORT ---
    Path("output").mkdir(exist_ok=True)
    
    # Single KO TI
    with open("output/single_ko_ti.json", "w") as f:
        json.dump(single_ti, f, indent=2)
    with open("output/single_ko_ti.tsv", "w") as f:
        f.write("rank\tgene\ttumor_collapse\thealthy_impact\ttherapeutic_index\tinfinite\n")
        for rank, r in enumerate(single_ti, 1):
            f.write(f"{rank}\t{r['gene']}\t{r['tumor_collapse']:.6f}\t{r['healthy_impact']:.6f}\t"
                    f"{r['therapeutic_index']:.2f}\t{r['is_infinite']}\n")
    
    # Dual KO TI
    with open("output/dual_ko_ti.json", "w") as f:
        json.dump(dual_ti, f, indent=2)
    with open("output/dual_ko_ti.tsv", "w") as f:
        f.write("rank\tgene_a\tgene_b\ttumor_collapse\thealthy_impact\ttherapeutic_index\tbliss_synergy\tloewe_synergy\tinfinite\n")
        for rank, r in enumerate(dual_ti, 1):
            f.write(f"{rank}\t{r['gene_a']}\t{r['gene_b']}\t{r['tumor_collapse']:.6f}\t"
                    f"{r['healthy_impact']:.6f}\t{r['therapeutic_index']:.2f}\t"
                    f"{r['bliss_synergy']:.6f}\t{r['loewe_synergy']:.6f}\t{r['is_infinite']}\n")
    
    # --- PRINT SUMMARIES ---
    print("\n[TOP 10 SINGLE KO by Therapeutic Index]")
    print(f"{'Rank':<5} {'Gene':<15} {'Tumor C':<10} {'Healthy I':<12} {'TI':<10}")
    print("-" * 55)
    for i, r in enumerate(single_ti[:10], 1):
        ti_str = "INF" if r['is_infinite'] else f"{r['therapeutic_index']:.2f}"
        print(f"{i:<5} {r['gene']:<15} {r['tumor_collapse']:<10.4f} {r['healthy_impact']:<12.4f} {ti_str:<10}")
    
    print("\n[TOP 10 DUAL KO by Therapeutic Index]")
    print(f"{'Rank':<5} {'Gene A':<15} {'Gene B':<15} {'Tumor C':<10} {'Healthy I':<12} {'TI':<10} {'Bliss':<8}")
    print("-" * 85)
    for i, r in enumerate(dual_ti[:10], 1):
        ti_str = "INF" if r['is_infinite'] else f"{r['therapeutic_index']:.2f}"
        print(f"{i:<5} {r['gene_a']:<15} {r['gene_b']:<15} {r['tumor_collapse']:<10.4f} "
              f"{r['healthy_impact']:<12.4f} {ti_str:<10} {r['bliss_synergy']:<8.4f}")
    
    # Zone preservation check
    print("\n[ZONE PRESERVATION] Healthy zone transition stats:")
    print(f"  Baseline mean: {zone_stats['healthy']['mean']:.4f}")
    print(f"  Target: < 10% shift in mean")
    
    # Compute predicted shift for top combos
    print("\n[PREDICTED HEALTHY SHIFT] Top 5 dual combos:")
    for r in dual_ti[:5]:
        shift = r['healthy_impact'] * zone_stats['healthy']['mean']
        print(f"  {r['gene_a']}+{r['gene_b']}: shift={shift:.4f} "
              f"({shift/zone_stats['healthy']['mean']*100:.1f}% of baseline)")
    
    print("\n[SUCCESS] Month 4 Week 3 Complete: Therapeutic Index Calculation")
    print("  - output/single_ko_ti.json, .tsv")
    print("  - output/dual_ko_ti.json, .tsv")


if __name__ == "__main__":
    main()