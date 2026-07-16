#!/usr/bin/env python3
"""
Month 2, Week 1: Transfer Entropy Engine — PATCHED VERSION
GPU-accelerated Kraskov-Stögbauer-Grassberger (KSG) Mutual Information Estimator
for directed Transfer Entropy T_{X→Y} between TFs and target programs.

FIXES APPLIED:
1. Numerical jitter (ε ~ N(0, 1e-6)) injected to break coordinate ties
2. Time-lag reconstruction uses phenotypic velocity vectors (Month 1) instead of simple index shift
3. KSG estimator computes true conditional information flow via PyTorch tensor operations
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List, Tuple

import numpy as np
import torch
import torch.nn.functional as F


def load_data(device: torch.device) -> Tuple[torch.Tensor, np.ndarray, torch.Tensor]:
    """Load gene expression, gene names, and spatial labels."""
    t0 = time.perf_counter()
    X = torch.from_numpy(np.load("output/nn_X.npy")).to(device, dtype=torch.float32)  # (15000, 2500)
    gene_names = np.loadtxt("output/nn_gene_names.tsv", dtype=str, delimiter='\t')
    y = torch.from_numpy(np.load("output/nn_y.npy")).to(device, dtype=torch.int64)  # (15000,)
    elapsed = time.perf_counter() - t0
    print(f"[LOAD] X {tuple(X.shape)}, gene_names {gene_names.shape}, y {tuple(y.shape)} in {elapsed:.3f}s")
    print(f"[LOAD] Zones: {torch.unique(y).cpu().numpy()}")
    return X, gene_names, y


def select_top_genes(
    X: torch.Tensor,
    gene_names: np.ndarray,
    y: torch.Tensor,
    n_genes: int = 100,
    zone: int = 1,
) -> Tuple[torch.Tensor, List[str], np.ndarray]:
    """Select top N highly variable genes in the Periphery zone."""
    print(f"[SELECT] Selecting top {n_genes} variable genes in zone {zone}...")
    t0 = time.perf_counter()
    
    mask = y == zone
    X_zone = X[mask]  # (N_zone, 2500)
    
    gene_var = X_zone.var(dim=0)  # (2500,)
    _, top_idx = torch.topk(gene_var, n_genes)
    top_idx_cpu = top_idx.cpu().numpy()
    
    X_selected = X[:, top_idx_cpu]  # (15000, 100)
    selected_names = gene_names[top_idx_cpu + 1].tolist()  # +1 for header 'gene'
    
    elapsed = time.perf_counter() - t0
    print(f"[SELECT] Selected {len(selected_names)} genes in {elapsed:.3f}s")
    print(f"[SELECT] Top 10: {selected_names[:10]}")
    print(f"[SELECT] Variance range: [{gene_var[top_idx_cpu].min():.4f}, {gene_var[top_idx_cpu].max():.4f}]")
    
    return X_selected, selected_names, top_idx_cpu


def inject_jitter(
    X: torch.Tensor,
    eps: float = 1e-6,
    seed: int = 42,
) -> torch.Tensor:
    """
    Inject tiny Gaussian noise to break coordinate ties for KSG estimator.
    This prevents distance degeneracy in high dimensions.
    """
    print(f"[JITTER] Injecting N(0, {eps}) noise...")
    gen = torch.Generator(device=X.device)
    gen.manual_seed(seed)
    noise = torch.normal(0, eps, size=X.shape, generator=gen, device=X.device, dtype=X.dtype)
    return X + noise


def load_velocity_field(device: torch.device) -> Tuple[torch.Tensor, torch.Tensor]:
    """Load phenotypic velocity vectors and magnitudes from Month 1."""
    print("[VELOCITY] Loading Month 1 phenotypic velocity field...")
    velocity = torch.from_numpy(np.load("output/phenotypic_velocity.npy")).to(device, dtype=torch.float32)
    magnitude = torch.from_numpy(np.load("output/phenotypic_velocity_magnitude.npy")).to(device, dtype=torch.float32)
    print(f"[VELOCITY] Velocity: {tuple(velocity.shape)}, Magnitude: {tuple(magnitude.shape)}")
    print(f"[VELOCITY] Magnitude range: [{magnitude.min():.4f}, {magnitude.max():.4f}]")
    return velocity, magnitude


def construct_pseudotime_neighbors(
    X_zone: torch.Tensor,          # (N_zone, n_genes) - expression in zone
    velocity: torch.Tensor,        # (15000, 32) - 32D latent velocity
    magnitude: torch.Tensor,       # (15000,) - velocity magnitude
    zone_mask: torch.Tensor,       # (15000,) - boolean mask for zone
    dt: float = 0.05,              # pseudotime step
    k_neighbors: int = 10,         # neighbors for local field estimation
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Construct pseudo-temporal pairs (Y_t, Y_{t-dt}, X_{t-dt}) using velocity field.
    
    Instead of shifting by index (which assumes perfect ordering), we use the 
    phenotypic velocity vectors to predict where each cell moves in expression space.
    
    For each cell i in the zone:
    - Y_{t-dt} = X_zone[i] (current expression)
    - Y_t ≈ X_zone[i] + v_i * dt (predicted future expression)
    - X_{t-dt} = X_zone[i] (source gene expression at current time)
    
    But since we only observe steady-state snapshots, we find nearest neighbors
    in the velocity direction to estimate the conditional distributions.
    """
    print(f"[PSEUDOTIME] Constructing pseudo-temporal pairs with dt={dt}, k={k_neighbors}...")
    t0 = time.perf_counter()
    
    N_zone = X_zone.shape[0]
    n_genes = X_zone.shape[1]
    device = X_zone.device
    
    # Get velocity for zone cells
    v_zone = velocity[zone_mask]  # (N_zone, 32)
    mag_zone = magnitude[zone_mask]  # (N_zone,)
    
    # Normalize velocity direction
    v_norm = v_zone / (mag_zone.unsqueeze(1) + 1e-8)  # (N_zone, 32)
    
    # Project velocity direction to gene expression space using PCA loadings
    # We'll use a simpler approach: find neighbors along velocity direction
    # in latent space, then map to expression space
    
    # Load latent coordinates for zone cells
    latent = torch.from_numpy(np.load("output/scvi_latent.npy")).to(device, dtype=torch.float32)
    latent_zone = latent[zone_mask]  # (N_zone, 32)
    
    # For each cell, find k nearest neighbors along velocity direction
    # Use cosine similarity with velocity vector to find forward neighbors
    knn_indices = torch.zeros(N_zone, k_neighbors, dtype=torch.int64, device=device)
    
    # Chunked k-NN search
    chunk_size = 500
    for i in range(0, N_zone, chunk_size):
        end = min(i + chunk_size, N_zone)
        chunk_latent = latent_zone[i:end]
        
        # Directional similarity: dot product with velocity direction
        # Higher dot = more aligned with velocity direction = "future" cells
        dir_sim = chunk_latent @ v_norm.T  # (chunk, N_zone)
        
        # We want neighbors with positive dot product (forward in velocity direction)
        # and close in Euclidean distance
        # Combine: penalize distance, reward alignment
        euclid = torch.cdist(chunk_latent, latent_zone, p=2)
        
        # Score: distance - λ * alignment (smaller = better future neighbor)
        # Only consider positive alignment
        align_mask = dir_sim > 0
        score = euclid.clone()
        score[~align_mask] = float('inf')
        
        # Top k by score
        _, idx = torch.topk(score, k_neighbors, dim=1, largest=False)
        knn_indices[i:end] = idx
        
        del dir_sim, euclid, score, align_mask
        if i % 5000 == 0:
            torch.cuda.empty_cache()
    
    # Build pseudo-temporal pairs
    # Y_{t-dt} = current cell expression (N_zone, n_genes)
    Y_tm1 = X_zone.clone()
    
    # Y_t = mean of k forward neighbors' expression
    Y_t = torch.zeros_like(X_zone)
    for i in range(N_zone):
        neighbors = X_zone[knn_indices[i]]  # (k, n_genes)
        Y_t[i] = neighbors.mean(dim=0)
    
    # X_{t-dt} = current cell expression (same as Y_{t-dt} for source)
    X_tm1 = X_zone.clone()
    
    # Filter out cells with no valid forward neighbors (velocity ~ 0)
    valid = mag_zone > 1e-4
    valid_count = valid.sum().item()
    print(f"[PSEUDOTIME] Valid cells with velocity > 1e-4: {valid_count}/{N_zone}")
    
    elapsed = time.perf_counter() - t0
    print(f"[PSEUDOTIME] Completed in {elapsed:.3f}s")
    
    return Y_t[valid], Y_tm1[valid], X_tm1[valid]


def ksg_mutual_information(
    x: torch.Tensor,
    y: torch.Tensor,
    k: int = 3,
) -> torch.Tensor:
    """
    KSG Mutual Information Estimator (Kraskov et al., 2004) on GPU.
    
    I(X;Y) = ψ(k) + ψ(N) - ⟨ψ(n_x + 1) + ψ(n_y + 1)⟩
    
    Args:
        x: (N, d_x) first variable
        y: (N, d_y) second variable
        k: number of neighbors
    
    Returns:
        MI estimate (scalar tensor)
    """
    N, dx = x.shape
    _, dy = y.shape
    device = x.device
    
    # Joint space
    joint = torch.cat([x, y], dim=1)  # (N, dx + dy)
    
    # Chunked distance computation to manage VRAM
    chunk_size = 500
    epsilons = torch.zeros(N, device=device)
    
    for i in range(0, N, chunk_size):
        end = min(i + chunk_size, N)
        chunk = joint[i:end]
        
        # Pairwise distances in joint space
        dists = torch.cdist(chunk, joint, p=2)  # (chunk, N)
        
        # Exclude self (diagonal)
        if i == 0:
            # Only for first chunk, diagonal elements are in first 'end' columns
            dists[:, :end].fill_diagonal_(float('inf'))
        else:
            # Diagonal is at column i:end for rows i:end
            local_diag = torch.arange(end - i, device=device)
            dists[local_diag, local_diag + i] = float('inf')
        
        # k-th nearest neighbor distance (epsilon)
        knn_dists, _ = torch.topk(dists, k, dim=1, largest=False)
        epsilons[i:end] = knn_dists[:, -1]
        
        del dists, knn_dists
        if i % 5000 == 0:
            torch.cuda.empty_cache()
    
    # Count neighbors within epsilon in marginal spaces
    # x space
    dists_x = torch.cdist(x, x, p=2)
    dists_x.fill_diagonal_(float('inf'))
    nx = (dists_x <= epsilons.unsqueeze(1)).sum(dim=1).float()
    
    # y space
    dists_y = torch.cdist(y, y, p=2)
    dists_y.fill_diagonal_(float('inf'))
    ny = (dists_y <= epsilons.unsqueeze(1)).sum(dim=1).float()
    
    del dists_x, dists_y
    torch.cuda.empty_cache()
    
    # Digamma function
    def digamma(z: torch.Tensor) -> torch.Tensor:
        return torch.log(z) - 1/(2*z) - 1/(12*z*z)
    
    psi_k = digamma(torch.tensor(float(k), device=device))
    psi_N = digamma(torch.tensor(float(N), device=device))
    psi_nx = digamma(nx + 1)
    psi_ny = digamma(ny + 1)
    
    mi = psi_k + psi_N - (psi_nx + psi_ny).mean()
    
    return mi.clamp(min=0)


def transfer_entropy(
    X: torch.Tensor,
    gene_names: List[str],
    y: torch.Tensor,
    velocity: torch.Tensor,
    magnitude: torch.Tensor,
    zone: int = 1,
    k: int = 3,
    dt: float = 0.05,
) -> Tuple[np.ndarray, List[Tuple[str, str, float]]]:
    """
    Compute Transfer Entropy T_{X→Y} for all TF→Target pairs in Periphery.
    
    T_{X→Y} = I(Y_t; X_{t-1} | Y_{t-1}) 
    = H(Y_t | Y_{t-1}) - H(Y_t | Y_{t-1}, X_{t-1})
    = I(Y_t; Y_{t-1}, X_{t-1}) - I(Y_t; Y_{t-1})
    
    Using pseudo-temporal pairs constructed from velocity field.
    """
    print(f"[TE] Computing Transfer Entropy for zone {zone}...")
    t0 = time.perf_counter()
    
    # Filter to Periphery
    mask = y == zone
    X_zone = X[mask]  # (N_zone, n_genes)
    N_zone = X_zone.shape[0]
    n_genes = X_zone.shape[1]
    
    print(f"[TE] Zone {zone} cells: {N_zone}, genes: {n_genes}")
    
    # Construct pseudo-temporal pairs using velocity field
    Y_t, Y_tm1, X_tm1 = construct_pseudotime_neighbors(
        X_zone, velocity, magnitude, mask, dt=dt, k_neighbors=5
    )
    
    N_eff = Y_t.shape[0]
    print(f"[TE] Effective pseudo-temporal samples: {N_eff}")
    
    if N_eff < 50:
        print("[TE] WARNING: Too few valid samples, falling back to index shift")
        # Fallback: simple index shift by transition score ordering
        scores = torch.from_numpy(np.load("output/csgt_transition_scores.npy")).to(X.device)
        scores_zone = scores[mask]
        order = scores_zone.argsort()
        X_ordered = X_zone[order]
        Y_t = X_ordered[1:]
        Y_tm1 = X_ordered[:-1]
        X_tm1 = X_ordered[:-1]
        N_eff = Y_t.shape[0]
    
    # Compute TE matrix: (n_genes, n_genes) where TE[i,j] = T_{gene_i → gene_j}
    te_matrix = torch.zeros(n_genes, n_genes, device=X.device)
    
    batch_size = 5
    
    for i in range(0, n_genes, batch_size):
        i_end = min(i + batch_size, n_genes)
        print(f"[TE] Processing source genes {i}-{i_end}/{n_genes}...")
        
        for j in range(n_genes):
            # Target gene j
            Y_t_j = Y_t[:, j:j+1]        # (N_eff, 1)
            Y_tm1_j = Y_tm1[:, j:j+1]    # (N_eff, 1)
            
            for ii in range(i, i_end):
                if ii == j:
                    te_matrix[ii, j] = 0
                    continue
                
                # Source gene ii at t-1
                X_tm1_ii = X_tm1[:, ii:ii+1]  # (N_eff, 1)
                
                # T_{X→Y} = I(Y_t; X_{t-1} | Y_{t-1})
                # = I(Y_t; Y_{t-1}, X_{t-1}) - I(Y_t; Y_{t-1})
                
                joint_cond = torch.cat([Y_tm1_j, X_tm1_ii], dim=1)  # (N_eff, 2)
                
                mi_full = ksg_mutual_information(Y_t_j, joint_cond, k=k)
                mi_cond = ksg_mutual_information(Y_t_j, Y_tm1_j, k=k)
                
                te = mi_full - mi_cond
                te_matrix[ii, j] = te.clamp(min=0)
    
    elapsed = time.perf_counter() - t0
    print(f"[TE] Transfer Entropy matrix computed in {elapsed:.3f}s")
    
    # Convert to numpy
    te_np = te_matrix.cpu().numpy()
    
    # Find top directional links
    top_links = []
    flat_idx = np.argsort(te_np.flatten())[::-1]
    for idx in flat_idx[:20]:
        i, j = np.unravel_index(idx, te_np.shape)
        if te_np[i, j] > 1e-10:  # Only meaningful links
            top_links.append((gene_names[i], gene_names[j], float(te_np[i, j])))
    
    return te_np, top_links


def export_results(
    te_matrix: np.ndarray,
    gene_names: List[str],
    top_links: List[Tuple[str, str, float]],
) -> None:
    """Save TE matrix and report top links."""
    print("[EXPORT] Saving Transfer Entropy matrix...")
    np.save("output/te_matrix.npy", te_matrix.astype(np.float32))
    
    # Save gene names mapping
    with open("output/te_gene_names.txt", "w") as f:
        for i, name in enumerate(gene_names):
            f.write(f"{i}\t{name}\n")
    
    # Report top links
    print("\n[TOP LINKS] Strongest directional causal links in Periphery:")
    for i, (src, tgt, val) in enumerate(top_links[:15]):
        print(f"  {i+1:2d}. {src:15s} → {tgt:15s}  TE = {val:.6f}")
    
    if not top_links:
        print("  (no significant links found)")
    
    # Statistics
    print(f"\n[STATS] TE Matrix: {te_matrix.shape}")
    print(f"[STATS] Range: [{te_matrix.min():.6f}, {te_matrix.max():.6f}]")
    print(f"[STATS] Mean: {te_matrix.mean():.6f}, Std: {te_matrix.std():.6f}")
    print(f"[STATS] Non-zero entries: {(te_matrix > 1e-10).sum()}/{te_matrix.size}")
    
    print(f"\n[EXPORT] Matrix saved: output/te_matrix.npy ({te_matrix.shape})")
    print(f"[EXPORT] Gene names: output/te_gene_names.txt")


def main() -> None:
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Computational Backend Initialized on: {device}")
    if device.type == 'cuda':
        print(f"   GPU: {torch.cuda.get_device_name(0)}")
        print(f"   VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
    
    Path("output").mkdir(exist_ok=True)
    
    # Load data
    X, gene_names, y = load_data(device)
    
    # Select top variable genes in Periphery
    X_sel, sel_names, sel_idx = select_top_genes(X, gene_names, y, n_genes=100, zone=1)
    
    # Inject jitter to break ties
    X_sel = inject_jitter(X_sel, eps=1e-6)
    
    # Load velocity field from Month 1
    velocity, magnitude = load_velocity_field(device)
    
    # Compute Transfer Entropy with pseudo-temporal reconstruction
    te_matrix, top_links = transfer_entropy(X_sel, sel_names, y, velocity, magnitude, zone=1, k=3, dt=0.05)
    
    # Export
    export_results(te_matrix, sel_names, top_links)
    
    # Memory report
    if device.type == 'cuda':
        print(f"\n[GPU] Peak memory allocated: {torch.cuda.max_memory_allocated() / 1e9:.2f} GB")
        print(f"[GPU] Peak memory reserved: {torch.cuda.max_memory_reserved() / 1e9:.2f} GB")
    
    print("\n[SUCCESS] Month 2 Week 1 Complete: Transfer Entropy Engine Executed")


if __name__ == "__main__":
    main()