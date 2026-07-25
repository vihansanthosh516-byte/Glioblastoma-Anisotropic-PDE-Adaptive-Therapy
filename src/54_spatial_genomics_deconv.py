#!/usr/bin/env python3
"""
Phase 1: Spatial Genomics Deconvolution -> PDE Parameter Fields

Maps canonical Neftel/Suva GBM cellular states (MES, AC, NPC, OPC)
to biophysical PDE parameters: diffusion coefficient D(x,y) and
proliferation rate rho(x,y).

References:
- Neftel et al., Cell 2019 (GBM cellular states)
- Suva et al., Cell 2014 (single-cell RNA-seq in glioma)
- IVY GAP spatial transcriptomics for validation
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter


# --------------------------------------------------------------------------- #
# Configuration (biophysically grounded bounds)
# --------------------------------------------------------------------------- #
GRID_SIZE = 128
DOMAIN_MM = 100.0  # 100 mm x 100 mm spatial domain

# Neftel state-to-parameter weights (literature-informed)
# MES/AC: invasive -> high D; NPC/OPC: proliferative -> high rho
MES_WEIGHT_D = 1.0
AC_WEIGHT_D = 0.8
NPC_WEIGHT_D = 0.1
OPC_WEIGHT_D = 0.1

MES_WEIGHT_RHO = 0.1
AC_WEIGHT_RHO = 0.1
NPC_WEIGHT_RHO = 1.0
OPC_WEIGHT_RHO = 0.8

# Physical bounds (GBM literature)
D_MIN, D_MAX = 0.01, 0.5      # mm^2/day
RHO_MIN, RHO_MAX = 0.005, 0.12  # day^-1

# Spatial smoothing for realistic tissue gradients
GAUSSIAN_SIGMA = 1.5  # voxels


# --------------------------------------------------------------------------- #
# Core Functions
# --------------------------------------------------------------------------- #
def generate_cell_fractions(grid_size: int = GRID_SIZE) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Generate synthetic spatial fractions for the four Neftel GBM states.

    Uses overlapping Gaussian centers to create realistic spatial heterogeneity
    with biomechanically plausible interfaces (MES/AC at periphery, NPC/OPC core).

    Returns:
        Tuple of (MES, AC, NPC, OPC) arrays, each shape (grid_size, grid_size),
        normalized so fractions sum to 1 at every pixel.
    """
    x = np.linspace(-1, 1, grid_size)
    y = np.linspace(-1, 1, grid_size)
    X, Y = np.meshgrid(x, y)

    # MES: invasive leading edge (upper-right)
    mes = np.exp(-((X - 0.4) ** 2 + (Y - 0.4) ** 2) / 0.15)

    # AC: astrocyte-like, also periphery but shifted (upper-left)
    ac = np.exp(-((X + 0.4) ** 2 + (Y - 0.4) ** 2) / 0.15)

    # NPC: neural progenitor, core region (lower-right)
    npc = np.exp(-((X - 0.2) ** 2 + (Y + 0.3) ** 2) / 0.12)

    # OPC: oligodendrocyte progenitor, core region (lower-left)
    opc = np.exp(-((X + 0.2) ** 2 + (Y + 0.3) ** 2) / 0.12)

    # Add subtle gradients for realism
    mes += 0.1 * np.exp(-((X + 0.1) ** 2 + (Y + 0.1) ** 2) / 0.5)
    ac  += 0.1 * np.exp(-((X - 0.1) ** 2 + (Y + 0.1) ** 2) / 0.5)
    npc += 0.1 * np.exp(-((X + 0.1) ** 2 + (Y - 0.1) ** 2) / 0.5)
    opc += 0.1 * np.exp(-((X - 0.1) ** 2 + (Y - 0.1) ** 2) / 0.5)

    # Normalize so fractions sum to 1 at every pixel
    total = mes + ac + npc + opc + 1e-12
    mes, ac, npc, opc = mes / total, ac / total, npc / total, opc / total

    return mes, ac, npc, opc


def compute_pde_fields(
    mes: np.ndarray,
    ac: np.ndarray,
    npc: np.ndarray,
    opc: np.ndarray,
    sigma: float = GAUSSIAN_SIGMA,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Map cell-state fractions to PDE parameter fields D(x,y) and rho(x,y).

    Applies Gaussian smoothing to reflect tissue-level averaging and
    rescales to biophysical ranges.
    """
    # Diffusion: driven by invasive states (MES, AC)
    D_raw = (MES_WEIGHT_D * mes + AC_WEIGHT_D * ac +
             NPC_WEIGHT_D * npc + OPC_WEIGHT_D * opc)
    D = gaussian_filter(D_raw, sigma=sigma)
    D = D_MIN + (D_MAX - D_MIN) * (D - D.min()) / (D.max() - D.min() + 1e-12)

    # Proliferation: driven by progenitor states (NPC, OPC)
    rho_raw = (MES_WEIGHT_RHO * mes + AC_WEIGHT_RHO * ac +
               NPC_WEIGHT_RHO * npc + OPC_WEIGHT_RHO * opc)
    rho = gaussian_filter(rho_raw, sigma=sigma)
    rho = RHO_MIN + (RHO_MAX - RHO_MIN) * (rho - rho.min()) / (rho.max() - rho.min() + 1e-12)

    return D, rho


def compute_metrics(D: np.ndarray, rho: np.ndarray) -> Dict:
    """Calculate summary statistics for PDE parameter fields."""
    return {
        "D": {
            "mean": float(np.mean(D)),
            "std": float(np.std(D)),
            "min": float(np.min(D)),
            "max": float(np.max(D)),
            "median": float(np.median(D)),
        },
        "rho": {
            "mean": float(np.mean(rho)),
            "std": float(np.std(rho)),
            "min": float(np.min(rho)),
            "max": float(np.max(rho)),
            "median": float(np.median(rho)),
        },
        "domain": {
            "grid_size": int(GRID_SIZE),
            "domain_mm": DOMAIN_MM,
            "dx_mm": DOMAIN_MM / GRID_SIZE,
        },
    }


def save_outputs(
    mes: np.ndarray,
    ac: np.ndarray,
    npc: np.ndarray,
    opc: np.ndarray,
    D: np.ndarray,
    rho: np.ndarray,
    metrics: Dict,
    out_dir: Path,
) -> None:
    """Save figure and JSON metrics to output directory."""
    out_dir.mkdir(parents=True, exist_ok=True)

    # 2x3 panel figure
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), constrained_layout=True)
    cmap_state = "magma"
    cmap_pde = "viridis"

    state_data = [
        ("MES fraction", mes, cmap_state),
        ("AC fraction", ac, cmap_state),
        ("NPC fraction", npc, cmap_state),
        ("OPC fraction", opc, cmap_state),
        (r"$D(x,y)\ \text{[mm}^2/\text{day]}$", D, cmap_pde),
        (r"$\rho(x,y)\ \text{[day}^{-1}\text{]}$", rho, cmap_pde),
    ]

    for ax, (title, field, cmap) in zip(axes.flat, state_data):
        im = ax.imshow(field, cmap=cmap, origin="lower",
                       extent=[0, DOMAIN_MM, 0, DOMAIN_MM])
        ax.set_title(title, fontsize=11, fontweight="bold")
        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        cbar = fig.colorbar(im, ax=ax, shrink=0.8)
        cbar.ax.tick_params(labelsize=8)

    fig.suptitle("Phase 1: Spatial Genomics Deconvolution -> PDE Parameter Fields",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.savefig(out_dir / "phase1_spatial_genomics_deconv.png",
                dpi=300, bbox_inches="tight")
    plt.close(fig)

    # JSON metrics
    with open(out_dir / "phase1_deconv_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[SAVE] Figure -> {out_dir / 'phase1_spatial_genomics_deconv.png'}")
    print(f"[SAVE] Metrics -> {out_dir / 'phase1_deconv_metrics.json'}")


def deconvolve_spatial_genomics() -> Tuple[Dict, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Main pipeline: generate cell fractions, compute PDE fields, return all arrays + metrics.
    """
    print("[PHASE 1] Generating Neftel cell-state fractions...")
    mes, ac, npc, opc = generate_cell_fractions(GRID_SIZE)

    # Sanity check: fractions sum to 1
    total = mes + ac + npc + opc
    assert np.allclose(total, 1.0, atol=1e-6), "Fractions do not sum to 1"
    print(f"[CHECK] Fractions sum to 1.0 (max deviation: {np.max(np.abs(total - 1.0)):.2e})")

    print("[PHASE 1] Computing PDE parameter fields...")
    D, rho = compute_pde_fields(mes, ac, npc, opc)

    metrics = compute_metrics(D, rho)

    print(f"[METRICS] D: mean={metrics['D']['mean']:.4f}, std={metrics['D']['std']:.4f} mm^2/day")
    print(f"[METRICS] rho: mean={metrics['rho']['mean']:.5f}, std={metrics['rho']['std']:.5f} day^-1")

    return metrics, mes, ac, npc, opc, D, rho


# --------------------------------------------------------------------------- #
# Entry Point
# --------------------------------------------------------------------------- #
def main() -> None:
    os.makedirs("output", exist_ok=True)

    metrics, mes, ac, npc, opc, D, rho = deconvolve_spatial_genomics()
    save_outputs(mes, ac, npc, opc, D, rho, metrics, Path("output"))

    print("[PHASE 1] Complete.")


if __name__ == "__main__":
    main()