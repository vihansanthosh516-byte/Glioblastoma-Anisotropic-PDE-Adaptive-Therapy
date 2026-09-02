#!/usr/bin/env python3
"""
Month 1: Biophysical Velocity Fields & Waddington Landscape Physics
Phenotypic Velocity Field Computation - GPU Accelerated

Computes local directional velocity vectors for 15,000 cells in 32D VAE latent space
using k-NN gradients of CSGT transition scores. Projects to 25x25 UMAP grid for
publication quiver visualization.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import torch
from matplotlib.colors import ListedColormap


def load_data(device: torch.device) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Load precomputed datasets directly to GPU memory."""
    t0 = time.perf_counter()
    latent = torch.from_numpy(np.load("output/scvi_latent.npy")).to(device, dtype=torch.float32)
    scores = torch.from_numpy(np.load("output/csgt_transition_scores.npy")).to(device, dtype=torch.float32)
    labels = torch.from_numpy(np.load("output/nn_y.npy")).to(device, dtype=torch.int64)
    elapsed = time.perf_counter() - t0
    print(f"[LOAD] Loaded latent {tuple(latent.shape)}, scores {tuple(scores.shape)}, labels {tuple(labels.shape)} in {elapsed:.3f}s")
    print(f"[LOAD] Latent range: [{latent.min():.4f}, {latent.max():.4f}], Scores range: [{scores.min():.4f}, {scores.max():.4f}]")
    print(f"[LOAD] Labels unique: {torch.unique(labels).cpu().numpy()}")
    return latent, scores, labels


def compute_velocity_vectors(
    latent: torch.Tensor,
    scores: torch.Tensor,
    k: int = 30,
    eps: float = 1e-8,
) -> Tuple[torch.Tensor, torch.Tensor]:
    """
    Compute phenotypic velocity vectors on GPU.
    
    For each cell, finds k-NN in latent space, computes gradient toward
    neighbors with higher transition scores, normalizes to unit vectors.
    
    Returns:
        velocity: (N, 32) velocity vectors in latent space
        magnitude: (N,) velocity magnitudes before normalization
    """
    N, D = latent.shape
    print(f"[VELOCITY] Computing {k}-NN on {N} cells x {D}D latent space...")
    t0 = time.perf_counter()

    # Full pairwise distance matrix on GPU: (N, N)
    # Use chunked computation to avoid OOM on 15k x 15k x 4 bytes = ~900MB
    chunk_size = 3000
    knn_indices = torch.empty(N, k, dtype=torch.int64, device=latent.device)
    
    for i in range(0, N, chunk_size):
        end = min(i + chunk_size, N)
        chunk = latent[i:end]
        # (chunk, N) distances
        dist_chunk = torch.cdist(chunk, latent, p=2)
        # Exclude self (diagonal) by setting to inf
        if i == 0:
            dist_chunk[:, :end].fill_diagonal_(float('inf'))
        # Top k smallest distances
        _, idx = torch.topk(dist_chunk, k, dim=1, largest=False)
        knn_indices[i:end] = idx
        del dist_chunk
        torch.cuda.empty_cache()

    elapsed = time.perf_counter() - t0
    print(f"[VELOCITY] k-NN search complete in {elapsed:.3f}s")

    print(f"[VELOCITY] Computing velocity vectors...")
    t0 = time.perf_counter()

    # Gather neighbor scores: (N, k)
    neighbor_scores = scores[knn_indices]  # (N, k)
    # Gather neighbor latent coords: (N, k, D)
    neighbor_latent = latent[knn_indices]  # (N, k, D)
    
    # Score differences: neighbor_score - self_score
    score_diff = neighbor_scores - scores.unsqueeze(1)  # (N, k)
    
    # Latent differences: neighbor - self
    latent_diff = neighbor_latent - latent.unsqueeze(1)  # (N, k, D)
    
    # Weight by positive score differences only (gradient toward higher transition)
    weights = torch.clamp(score_diff, min=0.0)  # (N, k)
    
    # Weighted sum of latent differences
    weighted_diff = latent_diff * weights.unsqueeze(-1)  # (N, k, D)
    velocity = weighted_diff.sum(dim=1)  # (N, D)
    
    # Magnitude before normalization
    magnitude = velocity.norm(dim=1)  # (N,)
    
    # Normalize to unit vectors (avoid division by zero)
    velocity = velocity / (magnitude.unsqueeze(1) + eps)
    
    elapsed = time.perf_counter() - t0
    print(f"[VELOCITY] Vector computation complete in {elapsed:.3f}s")
    print(f"[VELOCITY] Magnitude stats: mean={magnitude.mean():.6f}, max={magnitude.max():.6f}, min={magnitude.min():.6f}")
    print(f"[VELOCITY] Zero-magnitude cells: {(magnitude < eps).sum().item()}/{N}")
    
    return velocity, magnitude


def project_to_umap_grid(
    latent: torch.Tensor,
    velocity: torch.Tensor,
    labels: torch.Tensor,
    grid_size: int = 25,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Project latent coordinates to 2D UMAP-like grid and bin velocities.
    
    Uses PCA to 2D as a deterministic proxy for UMAP projection,
    then bins into uniform grid and averages velocities per cell.
    """
    print(f"[GRID] Projecting to 2D and binning to {grid_size}x{grid_size} grid...")
    t0 = time.perf_counter()

    # Move to CPU for projection and gridding
    latent_cpu = latent.cpu().numpy()
    velocity_cpu = velocity.cpu().numpy()
    labels_cpu = labels.cpu().numpy()

    # PCA to 2D for consistent projection
    from sklearn.decomposition import PCA
    pca = PCA(n_components=2, random_state=42)
    coords_2d = pca.fit_transform(latent_cpu)  # (N, 2)
    print(f"[GRID] PCA explained variance ratio: {pca.explained_variance_ratio_.sum():.4f}")

    # Normalize to [0, 1] range
    coords_min = coords_2d.min(axis=0)
    coords_max = coords_2d.max(axis=0)
    coords_norm = (coords_2d - coords_min) / (coords_max - coords_min + 1e-8)

    # Grid bin indices
    bin_x = (coords_norm[:, 0] * (grid_size - 1e-6)).astype(int)
    bin_y = (coords_norm[:, 1] * (grid_size - 1e-6)).astype(int)

    # Accumulate velocities and counts per grid cell
    grid_vx = np.zeros((grid_size, grid_size), dtype=np.float32)
    grid_vy = np.zeros((grid_size, grid_size), dtype=np.float32)
    grid_count = np.zeros((grid_size, grid_size), dtype=np.int32)
    grid_label = np.full((grid_size, grid_size), -1, dtype=np.int32)
    grid_mag = np.zeros((grid_size, grid_size), dtype=np.float32)

    np.add.at(grid_vx, (bin_y, bin_x), velocity_cpu[:, 0])
    np.add.at(grid_vy, (bin_y, bin_x), velocity_cpu[:, 1])
    np.add.at(grid_count, (bin_y, bin_x), 1)
    np.add.at(grid_mag, (bin_y, bin_x), np.linalg.norm(velocity_cpu, axis=1))
    # Majority label per bin
    for i in range(len(labels_cpu)):
        grid_label[bin_y[i], bin_x[i]] = labels_cpu[i]

    # Average velocities
    mask = grid_count > 0
    grid_vx[mask] /= grid_count[mask]
    grid_vy[mask] /= grid_count[mask]
    grid_mag[mask] /= grid_count[mask]

    # Grid centers for quiver
    xx, yy = np.meshgrid(np.arange(grid_size), np.arange(grid_size))
    grid_centers_x = (xx + 0.5) / grid_size
    grid_centers_y = (yy + 0.5) / grid_size

    # Map back to original coordinate scale for plotting
    grid_centers_x = grid_centers_x * (coords_max[0] - coords_min[0]) + coords_min[0]
    grid_centers_y = grid_centers_y * (coords_max[1] - coords_min[1]) + coords_min[1]

    elapsed = time.perf_counter() - t0
    print(f"[GRID] Gridding complete in {elapsed:.3f}s")
    print(f"[GRID] Occupied bins: {mask.sum()}/{grid_size*grid_size}")
    print(f"[GRID] Mean grid velocity magnitude: {grid_mag[mask].mean():.6f}")

    return (
        coords_2d, labels_cpu,
        grid_centers_x, grid_centers_y,
        grid_vx, grid_vy, grid_mag, grid_label, mask
    )


def plot_phenotypic_flux(
    coords_2d: np.ndarray,
    labels: np.ndarray,
    grid_x: np.ndarray,
    grid_y: np.ndarray,
    grid_vx: np.ndarray,
    grid_vy: np.ndarray,
    grid_mag: np.ndarray,
    grid_label: np.ndarray,
    mask: np.ndarray,
    output_path: Path,
) -> None:
    """Render publication-quality quiver plot of phenotypic velocity field."""
    print(f"[PLOT] Rendering quiver plot to {output_path}...")
    t0 = time.perf_counter()

    fig, ax = plt.subplots(figsize=(10, 10), dpi=300)
    
    # Zone colors: 0=Healthy(blue), 1=Periphery(orange), 2=Core(red)
    zone_colors = ['#1f77b4', '#ff7f0e', '#d62728']
    zone_labels = ['Healthy', 'Periphery', 'Core']
    cmap = ListedColormap(zone_colors)
    
    # Scatter cells colored by zone (subsample for visibility)
    n_cells = coords_2d.shape[0]
    idx = np.random.choice(n_cells, size=min(5000, n_cells), replace=False)
    scatter = ax.scatter(
        coords_2d[idx, 0], coords_2d[idx, 1],
        c=labels[idx], cmap=cmap, s=1, alpha=0.4, rasterized=True
    )
    
    # Quiver plot on grid
    q = ax.quiver(
        grid_x[mask], grid_y[mask],
        grid_vx[mask], grid_vy[mask],
        grid_mag[mask],
        cmap='hot', scale=25, width=0.003,
        headwidth=3, headlength=4, alpha=0.9,
        clim=(0, grid_mag[mask].max())
    )
    
    # Colorbar for velocity magnitude
    cbar = plt.colorbar(q, ax=ax, shrink=0.8, aspect=20)
    cbar.set_label('Velocity Magnitude', fontsize=12)
    
    # Legend for zones
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=c, label=l, alpha=0.6) 
        for c, l in zip(zone_colors, zone_labels)
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, framealpha=0.9)
    
    ax.set_title('Tumor Phenotypic Flux: Velocity Field of Cellular Transformation', fontsize=14, pad=15)
    ax.set_xlabel('PCA 1 (Latent Space)', fontsize=12)
    ax.set_ylabel('PCA 2 (Latent Space)', fontsize=12)
    ax.set_aspect('equal')
    ax.grid(True, alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    
    elapsed = time.perf_counter() - t0
    print(f"[PLOT] Saved in {elapsed:.3f}s")


def main() -> None:
    # ─── CUDA Gating ──────────────────────────────────────────────────────
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f'Computational Backend Initialized on: {device}')
    if device.type == 'cuda':
        print(f'   GPU: {torch.cuda.get_device_name(0)}')
        print(f'   VRAM: {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB')
        torch.cuda.empty_cache()
    
    # ─── Load Data ────────────────────────────────────────────────────────
    latent, scores, labels = load_data(device)
    
    # ─── Compute Velocity Vectors (GPU) ───────────────────────────────────
    velocity, magnitude = compute_velocity_vectors(latent, scores, k=30)
    
    # ─── Grid Projection & Quantification ─────────────────────────────────
    (coords_2d, labels_cpu,
     grid_x, grid_y, grid_vx, grid_vy, grid_mag, grid_label, mask) = project_to_umap_grid(
        latent, velocity, labels, grid_size=25
    )
    
    # ─── Publication Graphics ─────────────────────────────────────────────
    output_path = Path('output/tumor_phenotypic_flux.png')
    output_path.parent.mkdir(exist_ok=True)
    plot_phenotypic_flux(
        coords_2d, labels_cpu,
        grid_x, grid_y, grid_vx, grid_vy, grid_mag, grid_label, mask,
        output_path
    )
    
    # ─── Export velocity data for downstream months ───────────────────────
    np.save('output/phenotypic_velocity.npy', velocity.cpu().numpy())
    np.save('output/phenotypic_velocity_magnitude.npy', magnitude.cpu().numpy())
    np.save('output/phenotypic_velocity_pca2d.npy', coords_2d)
    print(f'[EXPORT] Saved velocity vectors to output/phenotypic_velocity.npy')
    print(f'[EXPORT] Saved velocity magnitudes to output/phenotypic_velocity_magnitude.npy')
    print(f'[EXPORT] Saved PCA 2D coords to output/phenotypic_velocity_pca2d.npy')
    
    print('\n[SUCCESS] Month 1 Foundation Complete: Phenotypic Velocity Field Computed')


if __name__ == '__main__':
    main()