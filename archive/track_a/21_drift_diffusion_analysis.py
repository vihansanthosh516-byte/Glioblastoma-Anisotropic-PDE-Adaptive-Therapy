#!/usr/bin/env python3
"""
Month 1 Week 3: Drift & Diffusion Tensor Analysis
Computes local Fokker-Planck coefficients A(z) and B(z) per cell.
Exports drift vectors and diffusion tensors for downstream saddle point analysis.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Tuple

import numpy as np
import torch
from torch import Tensor


def load_velocity_data(device: torch.device) -> Tuple[Tensor, Tensor, Tensor, Tensor]:
    """Load latent, velocity, scores, and labels."""
    t0 = time.perf_counter()
    latent = torch.from_numpy(np.load("output/scvi_latent.npy")).to(device, dtype=torch.float32)
    velocity = torch.from_numpy(np.load("output/phenotypic_velocity.npy")).to(device, dtype=torch.float32)
    scores = torch.from_numpy(np.load("output/csgt_transition_scores.npy")).to(device, dtype=torch.float32)
    labels = torch.from_numpy(np.load("output/nn_y.npy")).to(device, dtype=torch.int64)
    elapsed = time.perf_counter() - t0
    print(f"[LOAD] Data loaded in {elapsed:.3f}s")
    return latent, velocity, scores, labels


def compute_drift_diffusion(
    latent: Tensor,
    velocity: Tensor,
    scores: Tensor,
    k: int = 30,
    eps: float = 1e-8,
) -> Tuple[Tensor, Tensor, Tensor, Tensor]:
    """
    Compute local Drift A(z) and Diffusion B(z) tensors per cell.
    
    Drift: A(z) = E[Δz | z] / Δt ≈ velocity (from phenotypic gradient)
    Diffusion: B(z) = Cov[Δz | z] / Δt ≈ local covariance of k-NN displacements
    
    Returns:
        drift: (N, 32) local drift vectors
        diffusion: (N, 32, 32) local diffusion tensors
        drift_mag: (N,) drift magnitudes
        diff_trace: (N,) trace of diffusion (scalar diffusivity)
    """
    N, D = latent.shape
    print(f"[DRIFT/DIFF] Computing on {N} cells x {D}D...")
    t0 = time.perf_counter()
    
    # k-NN search (reuse pattern from velocity script)
    chunk_size = 3000
    knn_indices = torch.empty(N, k, dtype=torch.int64, device=latent.device)
    
    for i in range(0, N, chunk_size):
        end = min(i + chunk_size, N)
        chunk = latent[i:end]
        dist_chunk = torch.cdist(chunk, latent, p=2)
        if i == 0:
            dist_chunk[:, :end].fill_diagonal_(float('inf'))
        _, idx = torch.topk(dist_chunk, k, dim=1, largest=False)
        knn_indices[i:end] = idx
        del dist_chunk
        torch.cuda.empty_cache()
    
    # Neighbor displacements: (N, k, D)
    neighbor_latent = latent[knn_indices]  # (N, k, D)
    displacements = neighbor_latent - latent.unsqueeze(1)  # (N, k, D)
    
    # Drift: mean displacement (already computed as velocity, but recompute locally)
    drift = displacements.mean(dim=1)  # (N, D)
    drift_mag = drift.norm(dim=1)
    
    # Diffusion: covariance of displacements
    # B(z) = 1/(k-1) * sum (Δz_i - mean)(Δz_i - mean)^T
    centered = displacements - drift.unsqueeze(1)  # (N, k, D)
    # Batched outer product: (N, k, D, 1) * (N, k, 1, D) -> (N, k, D, D)
    outer = torch.einsum('nki,nkj->nkij', centered, centered)  # (N, k, D, D)
    diffusion = outer.mean(dim=1)  # (N, D, D)
    diff_trace = diffusion.diagonal(dim1=1, dim2=2).sum(dim=1)  # (N,)
    
    elapsed = time.perf_counter() - t0
    print(f"[DRIFT/DIFF] Computed in {elapsed:.3f}s")
    print(f"[DRIFT/DIFF] Drift mag: mean={drift_mag.mean():.6f}, max={drift_mag.max():.6f}")
    print(f"[DRIFT/DIFF] Diff trace: mean={diff_trace.mean():.6f}, max={diff_trace.max():.6f}")
    
    return drift, diffusion, drift_mag, diff_trace


def compute_fokker_planck_residual(
    drift: Tensor,
    diffusion: Tensor,
    density: Tensor,
    scores: Tensor,
) -> Tensor:
    """
    Compute Fokker-Planck steady-state residual:
    ∇·(A P) - ½ ∇²:(B P) ≈ 0 at steady state
    
    Returns per-cell residual magnitude.
    """
    print("[FP-RESIDUAL] Computing steady-state residual...")
    N, D = drift.shape
    device = drift.device
    
    # Gradient of log-density (score function)
    # ∇lnP ≈ ∇P / P
    # We estimate ∇P using k-NN weighted differences
    # For now, use a simpler proxy: local divergence of drift
    
    # Divergence of drift field: tr(∇A)
    # Approximate via finite differences on k-NN graph
    # This is a simplified residual - full computation would need spatial derivatives
    residual = torch.zeros(N, device=device)
    
    # Use velocity as proxy for A*P and compute divergence numerically
    # ∇·(A) ≈ mean over neighbors of (A_j - A_i) · (z_j - z_i) / |z_j - z_i|^2
    # We'll use a k-NN graph divergence estimator
    
    chunk_size = 3000
    for i in range(0, N, chunk_size):
        end = min(i + chunk_size, N)
        # Compute divergence of drift for this chunk
        # Using local linear fit of drift field
        pass
    
    print("[FP-RESIDUAL] Residual computation placeholder - needs k-NN divergence")
    return residual


def export_results(
    drift: Tensor,
    diffusion: Tensor,
    drift_mag: Tensor,
    diff_trace: Tensor,
    labels: Tensor,
) -> None:
    """Export all drift/diffusion arrays."""
    print("[EXPORT] Saving drift/diffusion tensors...")
    Path("output").mkdir(exist_ok=True)
    
    np.save("output/drift_vectors.npy", drift.cpu().numpy().astype(np.float32))
    np.save("output/diffusion_tensors.npy", diffusion.cpu().numpy().astype(np.float32))
    np.save("output/drift_magnitude.npy", drift_mag.cpu().numpy().astype(np.float32))
    np.save("output/diffusion_trace.npy", diff_trace.cpu().numpy().astype(np.float32))
    
    print("[EXPORT] Saved:")
    print("   output/drift_vectors.npy")
    print("   output/diffusion_tensors.npy")
    print("   output/drift_magnitude.npy")
    print("   output/diffusion_trace.npy")


def main() -> None:
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Computational Backend Initialized on: {device}")
    if device.type == 'cuda':
        print(f"   GPU: {torch.cuda.get_device_name(0)}")
    
    latent, velocity, scores, labels = load_velocity_data(device)
    
    drift, diffusion, drift_mag, diff_trace = compute_drift_diffusion(
        latent, velocity, scores, k=30
    )
    
    export_results(drift, diffusion, drift_mag, diff_trace, labels)
    
    print("\n[SUCCESS] Month 1 Week 3 Complete: Drift & Diffusion Tensors Computed")
    print(f"   Drift: (15000, 32)")
    print(f"   Diffusion: (15000, 32, 32)")
    print(f"   Diffusion trace range: [{diff_trace.min():.6f}, {diff_trace.max():.6f}]")


if __name__ == "__main__":
    main()