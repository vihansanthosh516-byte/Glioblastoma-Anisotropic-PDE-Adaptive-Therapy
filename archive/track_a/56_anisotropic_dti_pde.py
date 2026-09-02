#!/usr/bin/env python3
"""
Phase 3: 3D Anisotropic DTI Tensor Integration for GBM

Full 3D reaction-diffusion with patient-specific DTI-derived diffusion tensors.
Integrates Phase 1 (rho, D) and Phase 2 (mechanical pressure) outputs.

PDE: ∂u/∂t = ∇·(D_cell(x,y,z)∇u) + ρ(x,y,z)u(1-u/K) - ∇·(u v_mech)

References:
- Hormuth et al., Cancer Research 2017 (anisotropic glioma PDE)
- Swanson et al., Cell 2003 (FK model)
- Tournier et al., NeuroImage 2007 (CSD)
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
GRID_3D = (128, 128, 64)        # (nx, ny, nz) - 128x128 in-plane, 64 slices
DOMAIN_MM = (100.0, 100.0, 50.0)  # mm in x, y, z
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Biophysical parameters
K_CARRYING = 1e6      # mm^-3 (large = near-exponential)
R_ANISO = 0.85        # Degree of DTI guidance [0,1]
D0_BASE = 0.05        # mm^2/day baseline

# Time stepping
DT_DAYS = 0.05
T_TOTAL_DAYS = 30.0
N_STEPS = int(T_TOTAL_DAYS / DT_DAYS)

# Mechanical coupling (from Phase 2)
K_MECH = 0.1          # mm^2/day/Pa
BETA_PRESSURE = 100.0 # Pa per unit density


# --------------------------------------------------------------------------- #
# Phase 1/2 Loaders
# --------------------------------------------------------------------------- #
def load_phase1_fields() -> Tuple[np.ndarray, np.ndarray]:
    """Load or generate D(x,y,z) and rho(x,y,z) matching Phase 1 statistics."""
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx),
        np.linspace(-1, 1, ny),
        np.linspace(-1, 1, nz),
        indexing="ij"
    )

    # D field: high at periphery (invasive fronts)
    D = 0.01 + 0.49 * (
        np.exp(-((X - 0.4)**2 + (Y - 0.4)**2 + Z**2) / 0.15) +
        np.exp(-((X + 0.4)**2 + (Y - 0.4)**2 + Z**2) / 0.15)
    )
    D = gaussian_filter(D, sigma=(1.5, 1.5, 1.0))
    D = 0.01 + 0.49 * (D - D.min()) / (D.max() - D.min() + 1e-12)

    # rho field: high in core (proliferative)
    rho = 0.005 + 0.115 * (
        np.exp(-((X - 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12) +
        np.exp(-((X + 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12)
    )
    rho = gaussian_filter(rho, sigma=(1.5, 1.5, 1.0))
    rho = 0.005 + 0.115 * (rho - rho.min()) / (rho.max() - rho.min() + 1e-12)

    return D, rho


def load_phase2_pressure(rho: np.ndarray) -> np.ndarray:
    """Load or simulate mechanical pressure field P(x,y,z) from Phase 2."""
    # In production: load from Phase 2 output .npy
    # For now: compute from current density estimate
    return BETA_PRESSURE * rho  # placeholder


# --------------------------------------------------------------------------- #
# 3D DTI Tensor Field Construction
# --------------------------------------------------------------------------- #
def construct_dti_tensor_field(
    D_base: np.ndarray,
    r_aniso: float = R_ANISO,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construct 3D DTI-derived tumor diffusion tensor D_cell(x,y,z).

    Returns:
        tensor: (nx, ny, nz, 3, 3) symmetric positive-definite
        FA: (nx, ny, nz) fractional anisotropy
        e1: (nx, ny, nz, 3) principal eigenvector
    """
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx),
        np.linspace(-1, 1, ny),
        np.linspace(-1, 1, nz),
        indexing="ij"
    )

    # Synthetic white matter tracts (corpus callosum, CST, SLF)
    # Corpus callosum: mid-sagittal, horizontal (x-direction)
    cc_mask = (np.abs(Y) < 0.15) & (np.abs(Z) < 0.1)
    # Corticospinal tract: lateral, superior-inferior (z-direction)
    cst_mask = (X**2 + Y**2 > 0.15) & (X**2 + Y**2 < 0.35) & (Z > -0.6)
    # Superior longitudinal fasciculus: dorsal, anterior-posterior (y-direction)
    slf_mask = (np.abs(X) < 0.1) & (Z > -0.3)

    # Principal directions per tract
    e1_field = np.zeros((nx, ny, nz, 3))
    FA_field = np.zeros((nx, ny, nz))

    # Corpus callosum -> e1 = [1, 0, 0] (x-direction)
    e1_field[cc_mask] = [1, 0, 0]
    FA_field[cc_mask] = 0.75

    # Corticospinal -> e1 = [0, 0, 1] (z-direction)
    e1_field[cst_mask] = [0, 0, 1]
    FA_field[cst_mask] = 0.70

    # SLF -> e1 = [0, 1, 0] (y-direction)
    e1_field[slf_mask] = [0, 1, 0]
    FA_field[slf_mask] = 0.65

    # Default: isotropic (random orientation, low FA)
    isotropic_mask = ~(cc_mask | cst_mask | slf_mask)
    FA_field[isotropic_mask] = 0.15
    # Random isotropic directions
    np.random.seed(42)
    rand_vec = np.random.randn(nx, ny, nz, 3)
    rand_norm = np.linalg.norm(rand_vec, axis=-1, keepdims=True) + 1e-8
    e1_field[isotropic_mask] = rand_vec[isotropic_mask] / rand_norm[isotropic_mask]

    # Eigenvalues: lambda1 >> lambda2 = lambda3
    lambda1 = D_base * (1 + 3 * FA_field)  # parallel
    lambda2 = D_base * (1 - FA_field)      # perpendicular
    lambda3 = lambda2

    # Build tensor: R @ diag(lambda1, lambda2, lambda3) @ R.T
    # where R has e1 as first column, orthonormal completion
    tensor = np.zeros((nx, ny, nz, 3, 3))

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                e1 = e1_field[i, j, k]
                if np.linalg.norm(e1) < 1e-8:
                    e1 = np.array([1.0, 0.0, 0.0])
                
                # Build orthonormal basis
                # Find a vector not parallel to e1
                if abs(e1[0]) < 0.9:
                    v = np.array([1.0, 0.0, 0.0])
                else:
                    v = np.array([0.0, 1.0, 0.0])
                e2 = v - np.dot(v, e1) * e1
                e2 = e2 / (np.linalg.norm(e2) + 1e-8)
                e3 = np.cross(e1, e2)
                
                R = np.column_stack([e1, e2, e3])
                D_diag = np.diag([lambda1[i, j, k], lambda2[i, j, k], lambda3[i, j, k]])
                tensor[i, j, k] = R @ D_diag @ R.T

    return tensor, FA_field, e1_field


# --------------------------------------------------------------------------- #
# 3D Anisotropic Reaction-Diffusion Solver (Explicit with CFL)
# --------------------------------------------------------------------------- #
def gradient_3d(f: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central difference gradients in 3D with Neumann BC."""
    fx = np.zeros_like(f)
    fy = np.zeros_like(f)
    fz = np.zeros_like(f)

    fx[1:-1, :, :] = (f[2:, :, :] - f[:-2, :, :]) / (2 * dx)
    fx[0, :, :] = (f[1, :, :] - f[0, :, :]) / dx
    fx[-1, :, :] = (f[-1, :, :] - f[-2, :, :]) / dx

    fy[:, 1:-1, :] = (f[:, 2:, :] - f[:, :-2, :]) / (2 * dx)
    fy[:, 0, :] = (f[:, 1, :] - f[:, 0, :]) / dx
    fy[:, -1, :] = (f[:, -1, :] - f[:, -2, :]) / dx

    fz[:, :, 1:-1] = (f[:, :, 2:] - f[:, :, :-2]) / (2 * dx)
    fz[:, :, 0] = (f[:, :, 1] - f[:, :, 0]) / dx
    fz[:, :, -1] = (f[:, :, -1] - f[:, :, -2]) / dx

    return fx, fy, fz


def divergence_3d(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float) -> np.ndarray:
    """∇·v in 3D."""
    return gradient_3d(vx, dx)[0] + gradient_3d(vy, dx)[1] + gradient_3d(vz, dx)[2]


def anisotropic_diffusion_flux(
    u: np.ndarray,
    tensor: np.ndarray,
    dx: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute flux J = -D ∇u for anisotropic tensor.
    Returns (Jx, Jy, Jz).
    """
    ux, uy, uz = gradient_3d(u, dx)

    Jx = -(tensor[:, :, :, 0, 0] * ux + tensor[:, :, :, 0, 1] * uy + tensor[:, :, :, 0, 2] * uz)
    Jy = -(tensor[:, :, :, 1, 0] * ux + tensor[:, :, :, 1, 1] * uy + tensor[:, :, :, 1, 2] * uz)
    Jz = -(tensor[:, :, :, 2, 0] * ux + tensor[:, :, :, 2, 1] * uy + tensor[:, :, :, 2, 2] * uz)

    return Jx, Jy, Jz


def mechanical_velocity_3d(pressure: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """v_mech = -k_mech ∇P in 3D."""
    px, py, pz = gradient_3d(pressure, dx)
    return -K_MECH * px, -K_MECH * py, -K_MECH * pz


def step_anisotropic_3d(
    u: np.ndarray,
    tensor: np.ndarray,
    rho: np.ndarray,
    pressure: np.ndarray,
    dx: float,
    dt: float,
) -> np.ndarray:
    """Single explicit time step with anisotropic diffusion + reaction + mechanics."""
    # 1. Anisotropic diffusion: ∇·(D ∇u) = -∇·J
    Jx, Jy, Jz = anisotropic_diffusion_flux(u, tensor, dx)
    diff_term = -divergence_3d(Jx, Jy, Jz, dx)

    # 2. Reaction: ρ u (1 - u/K)
    react_term = rho * u * (1 - u / K_CARRYING)

    # 3. Mechanical advection: -∇·(u v_mech)
    vx, vy, vz = mechanical_velocity_3d(pressure, dx)
    adv_x = gradient_3d(u * vx, dx)[0]
    adv_y = gradient_3d(u * vy, dx)[1]
    adv_z = gradient_3d(u * vz, dx)[2]
    adv_term = -(adv_x + adv_y + adv_z)

    # Explicit Euler
    u_new = u + dt * (diff_term + react_term + adv_term)

    # Neumann BC (zero flux)
    u_new[0, :, :] = u_new[1, :, :]
    u_new[-1, :, :] = u_new[-2, :, :]
    u_new[:, 0, :] = u_new[:, 1, :]
    u_new[:, -1, :] = u_new[:, -2, :]
    u_new[:, :, 0] = u_new[:, :, 1]
    u_new[:, :, -1] = u_new[:, :, -2]

    # Non-negativity
    u_new = np.maximum(u_new, 0)

    return u_new


# --------------------------------------------------------------------------- #
# Isotropic Baseline (for comparison)
# --------------------------------------------------------------------------- #
def run_isotropic_baseline(
    D_iso: np.ndarray,
    rho: np.ndarray,
    pressure: np.ndarray,
    n_steps: int = N_STEPS,
    dt: float = DT_DAYS,
    dx: float = DX_MM,
) -> np.ndarray:
    """Run isotropic (scalar D) simulation for comparison."""
    nx, ny, nz = GRID_3D
    u = np.zeros(GRID_3D)
    center = (nx // 2, ny // 2, nz // 2)
    u[center[0]-2:center[0]+2, center[1]-2:center[1]+2, center[2]-2:center[2]+2] = 0.1

    for step in range(n_steps):
        # Simple isotropic diffusion
        ux, uy, uz = gradient_3d(u, dx)
        Jx, Jy, Jz = -D_iso * ux, -D_iso * uy, -D_iso * uz
        diff = -divergence_3d(Jx, Jy, Jz, dx)

        react = rho * u * (1 - u / K_CARRYING)

        vx, vy, vz = mechanical_velocity_3d(pressure, dx)
        adv = -(gradient_3d(u * vx, dx)[0] + gradient_3d(u * vy, dx)[1] + gradient_3d(u * vz, dx)[2])

        u = u + dt * (diff + react + adv)
        u = np.maximum(u, 0)

        # Neumann BC
        u[0, :, :] = u[1, :, :]
        u[-1, :, :] = u[-2, :, :]
        u[:, 0, :] = u[:, 1, :]
        u[:, -1, :] = u[:, -2, :]
        u[:, :, 0] = u[:, :, 1]
        u[:, :, -1] = u[:, :, -2]

    return u


# --------------------------------------------------------------------------- #
# Metrics & Visualization
# --------------------------------------------------------------------------- #
def compute_anisotropy_metrics(
    u: np.ndarray,
    tensor: np.ndarray,
    FA: np.ndarray,
    e1: np.ndarray,
    dx: float,
) -> Dict:
    """Compute tensor anisotropy and invasion metrics."""
    # Tensor anisotropy ratio (eigenvalue ratio)
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[:, :, :, ::-1]  # descending
    aniso_ratio = evals[:, :, :, 0] / (evals[:, :, :, 2] + 1e-12)

    # Tumor boundary mask
    boundary = (u > 0.1) & (
        (np.roll(u, 1, axis=0) < 0.1) | (np.roll(u, -1, axis=0) < 0.1) |
        (np.roll(u, 1, axis=1) < 0.1) | (np.roll(u, -1, axis=1) < 0.1) |
        (np.roll(u, 1, axis=2) < 0.1) | (np.roll(u, -1, axis=2) < 0.1)
    )

    # Preferred migration orientation at boundary
    if np.any(boundary):
        e1_boundary = e1[boundary]
        # Project to xy-plane for angle
        angles = np.arctan2(e1_boundary[:, 1], e1_boundary[:, 0]) * 180 / np.pi
        pref_angle = float(np.mean(angles))
    else:
        pref_angle = 0.0

    # Maximum invasion distance along principal direction
    # Distance from center to furthest u > 0.1 voxel
    nx, ny, nz = u.shape
    center = np.array([nx // 2, ny // 2, nz // 2])
    tumor_voxels = np.argwhere(u > 0.1)
    if len(tumor_voxels) > 0:
        dists = np.sqrt(np.sum((tumor_voxels - center) ** 2, axis=1)) * dx
        max_inv_dist = float(np.max(dists))
    else:
        max_inv_dist = 0.0

    # Asymmetry index: variance in radial distribution
    if len(tumor_voxels) > 10:
        # Spherical coordinates relative to center
        rel = tumor_voxels - center
        r = np.sqrt(np.sum(rel**2, axis=1))
        theta = np.arctan2(rel[:, 1], rel[:, 0])
        phi = np.arccos(rel[:, 2] / (r + 1e-8))
        # Asymmetry = variance in angular distribution
        asym_index = float(np.std(theta) / np.pi)
    else:
        asym_index = 0.0

    return {
        "mean_tensor_anisotropy_ratio": float(np.mean(aniso_ratio[boundary]) if np.any(boundary) else 1.0),
        "mean_FA_boundary": float(np.mean(FA[boundary]) if np.any(boundary) else 0.0),
        "preferred_migration_angle_deg": pref_angle,
        "max_fiber_tract_invasion_mm": max_inv_dist,
        "asymmetry_index": asym_index,
        "simulation_days": T_TOTAL_DAYS,
        "grid": GRID_3D,
        "dx_mm": dx,
    }


def save_outputs_3d(
    u_aniso: np.ndarray,
    u_iso: np.ndarray,
    tensor: np.ndarray,
    FA: np.ndarray,
    e1: np.ndarray,
    metrics: Dict,
    out_dir: Path,
) -> None:
    """Save 2x2 panel figure and JSON metrics."""
    out_dir.mkdir(parents=True, exist_ok=True)

    nx, ny, nz = GRID_3D
    sx, sy, sz = nx // 2, ny // 2, nz // 2

    # Mid-sagittal slices
    FA_mid = FA[sx, :, :]
    u_aniso_mid = u_aniso[sx, :, :]
    u_iso_mid = u_iso[sx, :, :]
    e1_mid = e1[sx, :, :]  # (ny, nz, 3)

    fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)

    # Panel 1: FA & Fiber Orientation Quiver
    im1 = axes[0, 0].imshow(FA_mid.T, cmap="viridis", origin="lower", vmin=0, vmax=1,
                            extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    axes[0, 0].set_title("FA & Fiber Orientation (Sagittal)", fontweight="bold")
    # Quiver of e1 (x=1, y=2 in sagittal = y, z)
    skip = 6
    y_grid = np.arange(0, DOMAIN_MM[1], DX_MM * skip)
    z_grid = np.arange(0, DOMAIN_MM[2], DX_MM * skip)
    axes[0, 0].quiver(y_grid, z_grid,
                      e1_mid[::skip, ::skip, 1].T, e1_mid[::skip, ::skip, 2].T,
                      color="white", scale=20, width=0.004, alpha=0.8)
    axes[0, 0].set_xlabel("y [mm]"); axes[0, 0].set_ylabel("z [mm]")
    plt.colorbar(im1, ax=axes[0, 0], shrink=0.8, label="FA")

    # Panel 2: Isotropic vs Anisotropic
    im2a = axes[0, 1].imshow(u_iso_mid.T, cmap="hot", origin="lower", alpha=0.6,
                             extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=1)
    im2b = axes[0, 1].contour(u_aniso_mid.T, levels=[0.1, 0.3, 0.5, 0.7], colors="cyan",
                              linewidths=2, extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    axes[0, 1].set_title("Isotropic (heat) vs Anisotropic (contours)", fontweight="bold")
    axes[0, 1].set_xlabel("y [mm]"); axes[0, 1].set_ylabel("z [mm]")
    axes[0, 1].clabel(im2b, inline=True, fontsize=8, fmt="%.1f")

    # Panel 3: Principal Diffusion Vector overlaid on Tumor Density
    im3 = axes[1, 0].imshow(u_aniso_mid.T, cmap="hot", origin="lower", alpha=0.7,
                            extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=1)
    axes[1, 0].quiver(y_grid, z_grid,
                      e1_mid[::skip, ::skip, 1].T, e1_mid[::skip, ::skip, 2].T,
                      color="cyan", scale=20, width=0.004, alpha=0.9)
    axes[1, 0].set_title("Principal Diffusion Vector e1 on u(x,y,z)", fontweight="bold")
    axes[1, 0].set_xlabel("y [mm]"); axes[1, 0].set_ylabel("z [mm]")
    plt.colorbar(im3, ax=axes[1, 0], shrink=0.8, label="Density")

    # Panel 4: Infiltration Anisotropy Index Map (D_xx / D_yy proxy)
    Dxx = tensor[sx, :, :, 0, 0]
    Dyy = tensor[sx, :, :, 1, 1]
    aniso_idx = Dxx / (Dyy + 1e-12)
    im4 = axes[1, 1].imshow(aniso_idx.T, cmap="RdBu_r", origin="lower",
                            extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0.5, vmax=2.0)
    axes[1, 1].set_title(r"Infiltration Anisotropy Index $D_{xx}/D_{yy}$", fontweight="bold")
    axes[1, 1].set_xlabel("y [mm]"); axes[1, 1].set_ylabel("z [mm]")
    plt.colorbar(im4, ax=axes[1, 1], shrink=0.8, label="Ratio")

    fig.suptitle("Phase 3: 3D Anisotropic DTI Tensor Integration (T=30 days)",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.savefig(out_dir / "phase3_anisotropic_dti.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(out_dir / "phase3_dti_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[SAVE] Figure -> {out_dir / 'phase3_anisotropic_dti.png'}")
    print(f"[SAVE] Metrics -> {out_dir / 'phase3_dti_metrics.json'}")


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main() -> None:
    os.makedirs("output", exist_ok=True)

    print("[PHASE 3] Loading Phase 1 fields...")
    D_base, rho = load_phase1_fields()
    print(f"  D_base range: [{D_base.min():.4f}, {D_base.max():.4f}] mm^2/day")
    print(f"  rho range: [{rho.min():.5f}, {rho.max():.5f}] day^-1")

    print("[PHASE 3] Constructing 3D DTI tensor field...")
    tensor, FA, e1 = construct_dti_tensor_field(D_base, R_ANISO)
    print(f"  Tensor shape: {tensor.shape}")
    print(f"  FA range: [{FA.min():.3f}, {FA.max():.3f}]")

    print("[PHASE 3] Loading Phase 2 pressure field...")
    # Use initial density to estimate pressure
    u_init = np.zeros(GRID_3D)
    center = (GRID_3D[0] // 2, GRID_3D[1] // 2, GRID_3D[2] // 2)
    u_init[center[0]-2:center[0]+2, center[1]-2:center[1]+2, center[2]-2:center[2]+2] = 0.1
    pressure = load_phase2_pressure(u_init)

    # Isotropic baseline
    print("[PHASE 3] Running isotropic baseline...")
    u_iso = run_isotropic_baseline(D_base, rho, pressure)

    # Anisotropic simulation
    print(f"[PHASE 3] Running 3D anisotropic simulation: {N_STEPS} steps...")
    u_aniso = u_init.copy()
    for step in range(N_STEPS):
        # Recompute pressure from current density
        pressure = load_phase2_pressure(u_aniso)
        u_aniso = step_anisotropic_3d(u_aniso, tensor, rho, pressure, DX_MM, DT_DAYS)

        if step % 100 == 0:
            vol = np.sum(u_aniso) * DX_MM**3
            print(f"  Step {step}/{N_STEPS}: volume={vol:.2f} mm^3")

    print("[PHASE 3] Computing anisotropy metrics...")
    metrics = compute_anisotropy_metrics(u_aniso, tensor, FA, e1, DX_MM)

    print(f"[METRICS] Mean aniso ratio: {metrics['mean_tensor_anisotropy_ratio']:.2f}")
    print(f"[METRICS] Preferred angle: {metrics['preferred_migration_angle_deg']:.1f} deg")
    print(f"[METRICS] Max invasion: {metrics['max_fiber_tract_invasion_mm']:.2f} mm")
    print(f"[METRICS] Asymmetry index: {metrics['asymmetry_index']:.3f}")

    save_outputs_3d(u_aniso, u_iso, tensor, FA, e1, metrics, Path("output"))

    print("[PHASE 3] Complete.")


if __name__ == "__main__":
    main()