#!/usr/bin/env python3
"""
Month 2, Week 2: Causal GRN Builder
Constructs directed gene regulatory network from Transfer Entropy matrix.
Identifies master switches via out-degree centrality.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import List, Tuple

import networkx as nx
import numpy as np


def load_data() -> Tuple[np.ndarray, List[str], np.ndarray]:
    """Load TE matrix, gene names, and cell labels."""
    t0 = time.perf_counter()
    
    te_matrix = np.load("output/te_matrix.npy")
    
    with open("output/te_gene_names.txt") as f:
        raw_genes = [line.strip() for line in f]
    gene_names = [g.split('\t')[-1] for g in raw_genes]
    
    labels = np.load("output/nn_y.npy")
    
    elapsed = time.perf_counter() - t0
    print(f"[LOAD] TE matrix: {te_matrix.shape}, genes: {len(gene_names)}, labels: {labels.shape} in {elapsed:.3f}s")
    return te_matrix, gene_names, labels


def threshold_te_matrix(
    te_matrix: np.ndarray,
    method: str = "percentile",
    percentile: float = 95,
) -> Tuple[np.ndarray, float]:
    """
    Apply threshold to TE matrix to filter noise.
    
    Args:
        te_matrix: (N, N) transfer entropy values
        method: "percentile" or "sigma"
        percentile: percentile cutoff (e.g., 95 for top 5%)
    
    Returns:
        Binary adjacency matrix, threshold value
    """
    print(f"[THRESHOLD] Applying {method} threshold (percentile={percentile})...")
    
    # Get upper triangle (exclude diagonal)
    mask = np.triu(np.ones_like(te_matrix, dtype=bool), k=1)
    vals = te_matrix[mask]
    
    if method == "percentile":
        threshold = np.percentile(vals, percentile)
    elif method == "sigma":
        threshold = vals.mean() + 2 * vals.std()
    else:
        raise ValueError(f"Unknown method: {method}")
    
    # Binarize
    adj_matrix = (te_matrix >= threshold).astype(np.int32)
    np.fill_diagonal(adj_matrix, 0)
    
    n_edges = adj_matrix.sum()
    print(f"[THRESHOLD] Threshold: {threshold:.6f}, Edges retained: {n_edges}/{te_matrix.size}")
    return adj_matrix, threshold


def build_causal_grn(
    adj_matrix: np.ndarray,
    gene_names: List[str],
) -> nx.DiGraph:
    """Build directed graph from adjacency matrix."""
    print("[GRAPH] Building directed causal GRN...")
    
    G = nx.DiGraph()
    
    # Add nodes
    for i, name in enumerate(gene_names):
        G.add_node(i, name=name)
    
    # Add edges
    rows, cols = np.where(adj_matrix > 0)
    for i, j in zip(rows, cols):
        G.add_edge(i, j, weight=1.0)
    
    print(f"[GRAPH] Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")
    return G


def compute_master_switches(G: nx.DiGraph, gene_names: List[str], top_k: int = 10) -> List[dict]:
    """
    Compute out-degree centrality to identify master switches.
    
    Returns list of dicts with gene info and out-degree.
    """
    print(f"[CENTRALITY] Computing out-degree centrality (top {top_k})...")
    
    out_degrees = dict(G.out_degree())
    
    # Sort by out-degree
    sorted_nodes = sorted(out_degrees.items(), key=lambda x: x[1], reverse=True)
    
    master_switches = []
    for i, (node_id, out_deg) in enumerate(sorted_nodes[:top_k]):
        name = gene_names[node_id]
        in_deg = G.in_degree(node_id)
        total_deg = G.degree(node_id)
        
        # Get downstream targets
        targets = [gene_names[t] for t in G.successors(node_id)]
        
        master_switches.append({
            "rank": i + 1,
            "gene_id": int(node_id),
            "gene_name": name,
            "out_degree": int(out_deg),
            "in_degree": int(in_deg),
            "total_degree": int(total_deg),
            "downstream_targets": targets,
            "n_targets": len(targets),
        })
    
    return master_switches


def export_grn_assets(
    G: nx.DiGraph,
    gene_names: List[str],
    master_switches: List[dict],
    adj_matrix: np.ndarray,
    threshold: float,
) -> None:
    """Export all GRN artifacts."""
    print("[EXPORT] Saving GRN assets...")
    Path("output").mkdir(exist_ok=True)
    
    # GraphML for Cytoscape
    nx.write_graphml(G, "output/causal_grn.graphml")
    print("  output/causal_grn.graphml")
    
    # Master switches TSV
    with open("output/master_switches.tsv", "w") as f:
        f.write("rank\tgene_id\tgene_name\tout_degree\tin_degree\ttotal_degree\tn_targets\tdownstream_targets\n")
        for ms in master_switches:
            targets_str = ",".join(ms["downstream_targets"])
            f.write(f"{ms['rank']}\t{ms['gene_id']}\t{ms['gene_name']}\t"
                    f"{ms['out_degree']}\t{ms['in_degree']}\t{ms['total_degree']}\t"
                    f"{ms['n_targets']}\t{targets_str}\n")
    print("  output/master_switches.tsv")
    
    # Adjacency matrix
    np.save("output/grn_adjacency.npy", adj_matrix)
    print("  output/grn_adjacency.npy")
    
    # Network metrics
    metrics = {
        "n_nodes": int(G.number_of_nodes()),
        "n_edges": int(G.number_of_edges()),
        "density": float(nx.density(G)),
        "threshold": float(threshold),
        "avg_out_degree": float(np.mean(list(dict(G.out_degree()).values()))),
        "avg_in_degree": float(np.mean(list(dict(G.in_degree()).values()))),
        "max_out_degree": int(max(dict(G.out_degree()).values())),
        "max_in_degree": int(max(dict(G.in_degree()).values())),
        "master_switches": master_switches,
    }
    
    with open("output/grn_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print("  output/grn_metrics.json")


def print_master_switch_summary(master_switches: List[dict]) -> None:
    """Print summary table of master switches."""
    print(f"\n{'='*100}")
    print("TOP MASTER SWITCHES (by Out-Degree Centrality)")
    print(f"{'='*100}")
    print(f"{'Rank':<5} {'Gene':<15} {'Out-Deg':<8} {'In-Deg':<8} {'Total':<8} {'Targets':<15} {'Downstream Targets'}")
    print("-" * 100)
    
    for ms in master_switches:
        targets_str = ", ".join(ms["downstream_targets"][:5])
        if len(ms["downstream_targets"]) > 5:
            targets_str += f" +{len(ms['downstream_targets']) - 5} more"
        print(f"{ms['rank']:<5} {ms['gene_name']:<15} {ms['out_degree']:<8} "
              f"{ms['in_degree']:<8} {ms['total_degree']:<8} {ms['n_targets']:<15} {targets_str}")
    
    print(f"{'='*100}\n")


def main() -> None:
    print("=" * 60)
    print("MONTH 2 WEEK 2: CAUSAL GRN BUILDER")
    print("=" * 60)
    
    # Load data
    te_matrix, gene_names, labels = load_data()
    
    # Threshold TE matrix
    adj_matrix, threshold = threshold_te_matrix(te_matrix, method="percentile", percentile=95)
    
    # Build causal GRN
    G = build_causal_grn(adj_matrix, gene_names)
    
    # Identify master switches
    master_switches = compute_master_switches(G, gene_names, top_k=10)
    
    # Export
    export_grn_assets(G, gene_names, master_switches, adj_matrix, threshold)
    
    # Summary
    print_master_switch_summary(master_switches)
    
    print("[SUCCESS] Month 2 Week 2 Complete: Causal GRN Built")
    print("  - Network: output/causal_grn.graphml")
    print("  - Master switches: output/master_switches.tsv")
    print("  - Metrics: output/grn_metrics.json")


if __name__ == "__main__":
    main()