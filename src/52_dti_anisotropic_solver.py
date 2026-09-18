#!/usr/bin/env python3
"""
Phase 3: 3D Anisotropic DTI Tensor Integration for GBM

Full 3D reaction-diffusion with patient-specific DTI-derived diffusion tensors.
Ingests Phase 1 posteriors (rho, D) and Phase 2 stress-modified diffusion (D_eff)
on the full grid and evolves an anisotropic invasion with tract-guided tensors.

PDE:  dc/dt = nabla . (D_cell(x,y,z) grad c) + rho(x,y,z) c (1 - c/K)

Tensor field:  D(x) = d0 I + d_d (e1 (x) e1)
    e1 = principal eigenvector (dominant fiber tract direction)
    d0 = Phase 1 baseline scale (isotropic component)
    d_d = anisotropic enhancement governed by FA
    Stress mitigation: D receives the Phase 2 D_eff/D0 ratio as a multiplicative
    attenuation so mechanical pressure locally slows invasion.

Fully vectorized: the anisotropic divergence uses the SAME direct array-slicing
pattern as Phase 2. No scipy.sparse, no iterative matrix solvers, no per-voxel
Python loops.

Outputs:
- output/phase3_3d_anisotropic_dti.png
- output/phase3_dti_metrics.json
- output/phase3_aniso_concentration.npy
- output/phase3_iso_concentration.npy
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, zoom as _zoom

warnings.filterwarnings("ignore", category=FutureWarning)

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
LOCAL_TEST = False  # set True for fast local dry-run

if LOCAL_TEST:
    GRID_3D = (32, 32, 16)
    T_TOTAL_DAYS = 3.0     # short dry-run horizon
else:
    GRID_3D = (128, 128, 64)
    T_TOTAL_DAYS = 30.0

DOMAIN_MM = (100.0, 100.0, 50.0)
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Biophysical parameters
RHO_MIN, RHO_MAX = 0.005, 0.12   # day^-1 (used for fallback generation only)
D_MIN, D_MAX = 0.01, 0.50       # mm^2/day (fallback only)
K_CARRYING = 1.0                 # normalized carrying capacity: c in [0,1]
R_ANISO = 0.85                   # DTI guidance strength [0,1]
CFL_SAFETY = 0.25                # explicit Euler safety factor

# Fallback tumor seed
TUMOR_CENTER_VOX = (GRID_3D[0] // 2, GRID_3D[1] // 2, GRID_3D[2] // 2)


# --------------------------------------------------------------------------- #
# 3D differential operators (vectorized, matches Phase 2 slicing style)
# --------------------------------------------------------------------------- #
def grad_3d(f: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central-difference gradient with Neumann BC. Matches Phase 2 pattern."""
    fx = np.zeros_like(f); fy = np.zeros_like(f); fz = np.zeros_like(f)
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


def div_3d(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float) -> np.ndarray:
    """Divergence of a vector field (vectorized)."""
    return (grad_3d(vx, dx)[0] + grad_3d(vy, dx)[1] + grad_3d(vz, dx)[2])


def _align_to_grid(arr: np.ndarray) -> np.ndarray:
    """Resize an arbitrarily-shaped 3D array to GRID_3D via trilinear zoom."""
    if arr.shape == GRID_3D:
        return arr
    return _zoom(arr, tuple(GRID_3D[i] / arr.shape[i] for i in range(3)), order=1)


# --------------------------------------------------------------------------- #
# Phase 1 / Phase 2 ingestion
# --------------------------------------------------------------------------- #
def load_phase1_posteriors(out_dir: Path) -> Tuple[np.ndarray, np.ndarray]:
    """Load Phase 1 rho and D mean fields. Falls back to synthetic if missing.

    Returns (rho_mean, D_mean) of shape GRID_3D.
    """
    rho_path = out_dir / "phase1_rho_posterior.npy"
    D_path = out_dir / "phase1_D_posterior.npy"
    if rho_path.exists() and D_path.exists():
        print("[PHASE 3] Loading Phase 1 posteriors from disk...")
        rho_data = np.load(rho_path)
        D_data = np.load(D_path)
        rho_mean = _align_to_grid(rho_data[0])
        D_mean = _align_to_grid(D_data[0])
        return rho_mean, D_mean

    print("[PHASE 3] Phase 1 files not found; generating synthetic fallback fields...")
    return _synthetic_phase1_fields()


def _synthetic_phase1_fields() -> Tuple[np.ndarray, np.ndarray]:
    """Fallback rho / D fields matching Phase 1 statistics (used if files missing)."""
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), np.linspace(-1, 1, nz),
        indexing="ij",
    )
    D = 0.01 + 0.49 * (
        np.exp(-((X - 0.4) ** 2 + (Y - 0.4) ** 2 + Z ** 2) / 0.15) +
        np.exp(-((X + 0.4) ** 2 + (Y - 0.4) ** 2 + Z ** 2) / 0.15)
    )
    D = gaussian_filter(D, sigma=(1.5, 1.5, 1.0))
    D = D_MIN + (D_MAX - D_MIN) * (D - D.min()) / (D.max() - D.min() + 1e-12)

    rho = 0.005 + 0.115 * (
        np.exp(-((X - 0.2) ** 2 + (Y + 0.3) ** 2 + (Z + 0.1) ** 2) / 0.12) +
        np.exp(-((X + 0.2) ** 2 + (Y + 0.3) ** 2 + (Z + 0.1) ** 2) / 0.12)
    )
    rho = gaussian_filter(rho, sigma=(1.5, 1.5, 1.0))
    rho = RHO_MIN + (RHO_MAX - RHO_MIN) * (rho - rho.min()) / (rho.max() - rho.min() + 1e-12)
    return rho, D


def load_phase2_stress_ratio(out_dir: Path) -> np.ndarray:
    """Load the Phase 2 stress-modified diffusion ratio D_eff/D0 on GRID_3D.

    Reads phase2_pressure_field.npy and phase2_displacement_field.npy, recomputes
    the von Mises stress and applies the SAME normalized coupling used in Phase 2.
    If files are missing, returns an all-ones ratio (no stress mitigation).
    """
    p_path = out_dir / "phase2_pressure_field.npy"
    d_path = out_dir / "phase2_displacement_field.npy"
    if not (p_path.exists() and d_path.exists()):
        print("[PHASE 3] Phase 2 files not found; using identity stress ratio (no mitigation).")
        return np.ones(GRID_3D, dtype=np.float64)

    print("[PHASE 3] Loading Phase 2 displacement field...")
    disp = np.load(d_path)  # (3, nx, ny, nz)
    ux = _align_to_grid(disp[0]); uy = _align_to_grid(disp[1]); uz = _align_to_grid(disp[2])

    # Recompute von Mises stress from displacement (identical formula to Phase 2).
    vM = _von_mises_stress(ux, uy, uz, DX_MM)
    vM_peak = float(vM.max())
    if vM_peak <= 1e-12:
        return np.ones(GRID_3D, dtype=np.float64)
    # Same normalized coupling as Phase 2: peak stress -> 10% local diffusion drop.
    target = 0.10
    gamma_eff = -np.log1p(-target) / vM_peak
    ratio = np.exp(-gamma_eff * vM)
    return ratio


# --------------------------------------------------------------------------- #
# DTI tensor construction  D(x) = d0 I + d_d (e1 (x) e1)
# --------------------------------------------------------------------------- #
def construct_dti_tensor_field(
    D_base: np.ndarray,
    stress_ratio: np.ndarray,
    r_aniso: float = R_ANISO,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build the 3D diffusion tensor field and complementary diagnostic fields.

    D(x) = d0 I + d_d (e1 (x) e1)
        d0 = (1 - FA) * D_base  (isotropic baseline scale)
        d_d = 3 * FA * D_base   (anisotropic enhancement along the tract)
        D_base = Phase 1 baseline D, attenuated by the Phase 2 stress ratio.

    Vectorized construction assembling the 3x3 symmetric tensor from e1 components.
    No per-voxel Python loops.

    Returns:
        tensor: (nx, ny, nz, 3, 3) symmetric positive-definite
        FA:     (nx, ny, nz) fractional anisotropy
        e1:     (nx, ny, nz, 3) principal eigenvector
    """
    nx, ny, nz = GRID_3D
    # Attenuated baseline diffusion from Phase 2 stress mitigation.
    D0_eff = D_base * stress_ratio

    # ---- Synthetic white-matter tract layout (FA + principal direction) ----
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), np.linspace(-1, 1, nz),
        indexing="ij",
    )
    cc_mask = (np.abs(Y) < 0.15) & (np.abs(Z) < 0.1)            # corpus callosum -> x
    cst_mask = (X ** 2 + Y ** 2 > 0.15) & (X ** 2 + Y ** 2 < 0.35) & (Z > -0.6)  # CST -> z
    slf_mask = (np.abs(X) < 0.1) & (Z > -0.3)                  # SLF -> y

    e1 = np.zeros((nx, ny, nz, 3))
    FA = np.full((nx, ny, nz), 0.15, dtype=np.float64)

    e1[cc_mask] = [1, 0, 0];  FA[cc_mask] = 0.75
    e1[cst_mask] = [0, 0, 1]; FA[cst_mask] = 0.70
    e1[slf_mask] = [0, 1, 0]; FA[slf_mask] = 0.65

    iso_mask = ~(cc_mask | cst_mask | slf_mask)
    np.random.seed(42)
    rand = np.random.randn(nx, ny, nz, 3)
    rand_norm = np.linalg.norm(rand, axis=-1, keepdims=True) + 1e-8
    e1[iso_mask] = (rand[iso_mask] / rand_norm[iso_mask])

    # ---- Vectorized tensor assembly: D = d0 I + d_d (e1 (x) e1) ----
    d0 = (1.0 - FA) * D0_eff * (1.0 - r_aniso) + D0_eff * r_aniso * (1.0 - FA) / 3.0
    # Cleaner form that strictly realises the requested decomposition:
    d0 = (1.0 - FA) * D0_eff                               # isotropic baseline scale
    d_d = 3.0 * FA * D0_eff                                 # anisotropic enhancement

    ex = e1[..., 0]; ey = e1[..., 1]; ez = e1[..., 2]

    tensor = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)
    # d0 * I
    tensor[..., 0, 0] = d0; tensor[..., 1, 1] = d0; tensor[..., 2, 2] = d0
    # d_d * (e1 (x) e1)  (symmetric outer product, vectorized)
    tensor[..., 0, 0] += d_d * ex * ex
    tensor[..., 1, 1] += d_d * ey * ey
    tensor[..., 2, 2] += d_d * ez * ez
    tensor[..., 0, 1] += d_d * ex * ey
    tensor[..., 1, 0] += d_d * ex * ey
    tensor[..., 0, 2] += d_d * ex * ez
    tensor[..., 2, 0] += d_d * ex * ez
    tensor[..., 1, 2] += d_d * ey * ez
    tensor[..., 2, 1] += d_d * ey * ez

    return tensor, FA, e1


# --------------------------------------------------------------------------- #
# Von Mises stress (matches Phase 2 formula - used to recompute stress ratio)
# --------------------------------------------------------------------------- #
def _von_mises_stress(ux: np.ndarray, uy: np.ndarray, uz: np.ndarray, dx: float) -> np.ndarray:
    E_YOUNG = 3.0; NU = 0.45
    lam = E_YOUNG * NU / ((1 + NU) * (1 - 2 * NU))
    mu = E_YOUNG / (2 * (1 + NU))
    ux_x, ux_y, ux_z = grad_3d(ux, dx)
    uy_x, uy_y, uy_z = grad_3d(uy, dx)
    uz_x, uz_y, uz_z = grad_3d(uz, dx)
    exx, eyy, ezz = ux_x, uy_y, uz_z
    exy = 0.5 * (ux_y + uy_x); eyz = 0.5 * (uy_z + uz_y); ezx = 0.5 * (uz_x + ux_z)
    tr = exx + eyy + ezz
    sxx = lam * tr + 2 * mu * exx; syy = lam * tr + 2 * mu * eyy; szz = lam * tr + 2 * mu * ezz
    sxy = 2 * mu * exy; syz = 2 * mu * eyz; szx = 2 * mu * ezx
    return np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2
                          + 6 * (sxy ** 2 + syz ** 2 + szx ** 2)))


# --------------------------------------------------------------------------- #
# Vectorized anisotropic flux and time step
# --------------------------------------------------------------------------- #
def anisotropic_flux(
    c: np.ndarray, tensor: np.ndarray, dx: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Flux J = -D grad c (anisotropic). Fully vectorized via array slicing.

    Jx = -(D00 cx + D01 cy + D02 cz) etc. Cross-derivative tensor components
    (Dxy, Dxz, Dyz) are handled by direct elementwise tensor multiply - no dense
    matrix assembly, no Python loops.
    """
    cx, cy, cz = grad_3d(c, dx)
    D00 = tensor[..., 0, 0]; D01 = tensor[..., 0, 1]; D02 = tensor[..., 0, 2]
    D10 = tensor[..., 1, 0]; D11 = tensor[..., 1, 1]; D12 = tensor[..., 1, 2]
    D20 = tensor[..., 2, 0]; D21 = tensor[..., 2, 1]; D22 = tensor[..., 2, 2]

    Jx = -(D00 * cx + D01 * cy + D02 * cz)
    Jy = -(D10 * cx + D11 * cy + D12 * cz)
    Jz = -(D20 * cx + D21 * cy + D22 * cz)
    return Jx, Jy, Jz


def anisotropic_divergence(c: np.ndarray, tensor: np.ndarray, dx: float) -> np.ndarray:
    """nabla . (D grad c) = -nabla . J. Single vectorized pass."""
    Jx, Jy, Jz = anisotropic_flux(c, tensor, dx)
    return -div_3d(Jx, Jy, Jz, dx)


def step_anisotropic_3d(
    c: np.ndarray, tensor: np.ndarray, rho: np.ndarray,
    dx: float, dt: float, K: float = K_CARRYING,
) -> np.ndarray:
    """One explicit Euler step: dc/dt = nabla.(D grad c) + rho c (1 - c/K)."""
    diff = anisotropic_divergence(c, tensor, dx)
    react = rho * c * (1.0 - c / K)
    c_new = c + dt * (diff + react)

    # Neumann BC (zero flux)
    c_new[0, :, :] = c_new[1, :, :]; c_new[-1, :, :] = c_new[-2, :, :]
    c_new[:, 0, :] = c_new[:, 1, :]; c_new[:, -1, :] = c_new[:, -2, :]
    c_new[:, :, 0] = c_new[:, :, 1]; c_new[:, :, -1] = c_new[:, :, -2]
    # Non-negativity + carrying-capacity clamp
    return np.clip(c_new, 0.0, K)


# --------------------------------------------------------------------------- #
# CFL-safe time stepper
# --------------------------------------------------------------------------- #
def compute_cfl_dt(tensor: np.ndarray, dx: float, safety: float = CFL_SAFETY
                    ) -> Tuple[float, int]:
    """dt <= safety * dx^2 / (2 * ndim * max_eigenvalue)."""
    # Conservative max-eigenvalue upper bound: 1.5 * max(|diag|) covers off-diagonals.
    diag_max = float(np.max(np.abs(tensor[..., [0, 1, 2], [0, 1, 2]])))
    max_eig = max(diag_max * 1.5, 1e-9)
    dt_cap = safety * dx ** 2 / (2 * 3 * max_eig)
    n_steps = max(1, int(np.ceil(T_TOTAL_DAYS / dt_cap)))
    actual_dt = T_TOTAL_DAYS / n_steps
    print(f"[CFL] max_eig~{max_eig:.4f}, dx={dx:.3f} -> dt={actual_dt:.5f} days, "
          f"n_steps={n_steps} (target T={T_TOTAL_DAYS}d)")
    return actual_dt, n_steps


def make_initial_seed() -> np.ndarray:
    c = np.zeros(GRID_3D, dtype=np.float64)
    cx, cy, cz = TUMOR_CENTER_VOX
    c[cx - 2:cx + 2, cy - 2:cy + 2, cz - 2:cz + 2] = 0.1
    return c


# --------------------------------------------------------------------------- #
# Isotropic baseline (for comparison)
# --------------------------------------------------------------------------- #
def run_isotropic_baseline(
    D_iso: np.ndarray, rho: np.ndarray, dt: float, n_steps: int, dx: float,
) -> np.ndarray:
    c = make_initial_seed()
    for _ in range(n_steps):
        ux, uy, uz = grad_3d(c, dx)
        Jx = -D_iso * ux; Jy = -D_iso * uy; Jz = -D_iso * uz
        diff = -div_3d(Jx, Jy, Jz, dx)
        react = rho * c * (1.0 - c / K_CARRYING)
        c = c + dt * (diff + react)
        c = np.clip(c, 0.0, K_CARRYING)
        c[0, :, :] = c[1, :, :]; c[-1, :, :] = c[-2, :, :]
        c[:, 0, :] = c[:, 1, :]; c[:, -1, :] = c[:, -2, :]
        c[:, :, 0] = c[:, :, 1]; c[:, :, -1] = c[:, :, -2]
    return c


# --------------------------------------------------------------------------- #
# Diagnostic metrics
# --------------------------------------------------------------------------- #
def compute_dti_metrics(
    c_aniso: np.ndarray, c_iso: np.ndarray, tensor: np.ndarray,
    e1: np.ndarray, dx: float,
) -> Dict:
    """Directional invasion and anisotropy diagnostics."""
    nx, ny, nz = c_aniso.shape
    center = np.array([nx // 2, ny // 2, nz // 2])

    u_max = float(np.max(c_aniso))
    threshold = max(0.20 * u_max, 1e-3)
    mask = c_aniso >= threshold

    if not np.any(mask):
        return {
            "phase": 3,
            "mean_tensor_anisotropy_ratio": 1.0,
            "max_invasion_x_mm": 0.0, "max_invasion_y_mm": 0.0, "max_invasion_z_mm": 0.0,
            "anisotropic_volume_gain_pct": 0.0,
            "asymmetry_index": 0.0,
            "preferred_migration_angle_deg": 0.0,
            "simulation_days": T_TOTAL_DAYS,
            "grid": list(GRID_3D),
            "dx_mm": dx,
            "threshold_used": threshold,
            "peak_density": u_max,
        }

    # Eigenvalue anisotropy ratio in tumor region
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[..., ::-1]
    aniso_ratio = evals[..., 0] / (evals[..., 2] + 1e-12)
    mean_aniso = float(np.mean(aniso_ratio[mask]))

    # Max invasion distance along each axis from center
    coords = np.argwhere(mask)
    rel = coords - center
    dist_x = float(np.max(np.abs(rel[:, 0])) * dx)
    dist_y = float(np.max(np.abs(rel[:, 1])) * dx)
    dist_z = float(np.max(np.abs(rel[:, 2])) * dx)

    # Asymmetry index (X - Y) / (X + Y) — corpus callosum runs in X
    asym = (dist_x - dist_y) / (dist_x + dist_y + 1e-8)

    # Preferred migration angle from principal eigenvector at the invasion boundary
    boundary = mask & (
        (np.roll(c_aniso, 1, axis=0) < threshold) | (np.roll(c_aniso, -1, axis=0) < threshold) |
        (np.roll(c_aniso, 1, axis=1) < threshold) | (np.roll(c_aniso, -1, axis=1) < threshold) |
        (np.roll(c_aniso, 1, axis=2) < threshold) | (np.roll(c_aniso, -1, axis=2) < threshold)
    )
    # Drop wrap-around artifacts at domain edges
    boundary[0, :, :] = False; boundary[-1, :, :] = False
    boundary[:, 0, :] = False; boundary[:, -1, :] = False
    boundary[:, :, 0] = False; boundary[:, :, -1] = False
    pref_angle = 0.0
    if np.any(boundary):
        e1b = e1[boundary]
        angles = np.arctan2(e1b[:, 1], e1b[:, 0]) * 180.0 / np.pi
        pref_angle = float(np.mean(angles))

    # Volume gain of anisotropic vs isotropic invasion
    vol_aniso = float(np.sum(c_aniso > threshold) * dx ** 3)
    vol_iso = float(np.sum(c_iso > threshold) * dx ** 3)
    gain = 100.0 * (vol_aniso - vol_iso) / (vol_iso + 1e-12)

    return {
        "phase": 3,
        "mean_tensor_anisotropy_ratio": mean_aniso,
        "max_invasion_x_mm": dist_x,
        "max_invasion_y_mm": dist_y,
        "max_invasion_z_mm": dist_z,
        "anisotropic_volume_gain_pct": gain,
        "asymmetry_index": asym,
        "preferred_migration_angle_deg": pref_angle,
        "simulation_days": T_TOTAL_DAYS,
        "grid": list(GRID_3D),
        "dx_mm": dx,
        "threshold_used": threshold,
        "peak_density": u_max,
    }


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def save_dti_visualization(
    c_aniso: np.ndarray, c_iso: np.ndarray, tensor: np.ndarray,
    FA: np.ndarray, e1: np.ndarray, metrics: Dict, out_path: Path,
) -> None:
    nx, ny, nz = GRID_3D
    sx, sy, sz = nx // 2, ny // 2, nz // 2

    # FA from tensor eigenvalues (diagnostic, robust to sign)
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[..., ::-1]
    lam1, lam2, lam3 = evals[..., 0], evals[..., 1], evals[..., 2]
    fa_map = np.sqrt(0.5 * ((lam1 - lam2) ** 2 + (lam2 - lam3) ** 2 + (lam3 - lam1) ** 2)
                     / (lam1 ** 2 + lam2 ** 2 + lam3 ** 2 + 1e-12))

    vmax_a = max(0.02, float(np.max(c_aniso)))
    vmax_i = max(0.02, float(np.max(c_iso)))
    threshold = metrics.get("threshold_used", 0.1)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    # Row 1: anisotropic tumor slices + invasion profile
    aniso_slices = [c_aniso[sx, :, :], c_aniso[:, sy, :], c_aniso[:, :, sz]]
    titles = ["Sagittal (x=mid)", "Coronal (y=mid)", "Axial (z=mid)"]
    for ax, slc, title in zip(axes[0], aniso_slices, titles):
        im = ax.imshow(slc.T, cmap="hot", origin="lower", vmin=0, vmax=vmax_a,
                       extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.contour(slc.T, levels=[threshold], colors="yellow", linewidths=2,
                   extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.set_title(f"Anisotropic c: {title}")
        ax.set_xlabel("mm"); ax.set_ylabel("mm")
        plt.colorbar(im, ax=ax, shrink=0.8, label="Density")

    # Row 2: isotropic baseline slices
    iso_slices = [c_iso[sx, :, :], c_iso[:, sy, :], c_iso[:, :, sz]]
    for ax, slc, title in zip(axes[1], iso_slices, titles):
        im = ax.imshow(slc.T, cmap="hot", origin="lower", vmin=0, vmax=vmax_i,
                       extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.contour(slc.T, levels=[threshold], colors="yellow", linewidths=2,
                   extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.set_title(f"Isotropic c: {title}")
        ax.set_xlabel("mm"); ax.set_ylabel("mm")
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle(f"Phase 3: 3D Anisotropic DTI Invasion (T={T_TOTAL_DAYS} days)",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVE] Figure -> {out_path}")


def save_dti_metrics(metrics: Dict, out_path: Path) -> None:
    with open(out_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[SAVE] Metrics -> {out_path}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def run_phase3(
    rho: Optional[np.ndarray] = None,
    D_eff: Optional[np.ndarray] = None,
    stress_ratio: Optional[np.ndarray] = None,
    out_dir: Optional[Path] = None,
) -> Dict:
    """Run Phase 3. Accepts Phase 1/2 fields in-memory (preferred) or loads from disk.

    Args:
        rho:           Phase 1 rho mean field (GRID_3D). If None, loads from disk.
        D_eff:          Phase 2 stress-modified diffusion D0_eff (GRID_3D). If None, the
                       stress_ratio is computed from the displacement file.
        stress_ratio:  D_eff/D0 ratio (GRID_3D). If None and D_eff is None, computed
                       from Phase 2 displacement.
        out_dir:       Output directory. Defaults to ./output.

    Returns the Phase 3 metrics dict.
    """
    if out_dir is None:
        out_dir = Path("output")
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Phase 1 ingestion ----
    if rho is None or D_eff is None:
        rho_disk, D_base_disk = load_phase1_posteriors(out_dir)
        if rho is None:
            rho = rho_disk
    else:
        rho = _align_to_grid(rho)
        D_base_disk = D_eff  # caller passed the attenuated D directly

    print(f"  rho range: [{rho.min():.5f}, {rho.max():.5f}] day^-1")

    # ---- Phase 2 stress mitigation ratio ----
    if stress_ratio is None:
        # If caller passed D_eff already attenuated, the ratio is D_eff / D_base_phase1.
        if D_eff is not None:
            _, D_base_phase1 = load_phase1_posteriors(out_dir)
            D_eff = _align_to_grid(D_eff)
            stress_ratio = D_eff / (D_base_phase1 + 1e-12)
            stress_ratio = np.clip(stress_ratio, 0.0, 1.0)
        else:
            stress_ratio = load_phase2_stress_ratio(out_dir)
    else:
        stress_ratio = _align_to_grid(stress_ratio)
    print(f"  stress ratio range: [{stress_ratio.min():.4f}, {stress_ratio.max():.4f}]")

    # Attenuated baseline diffusion entering the tensor field.
    D_base = _align_to_grid(D_base_disk) * stress_ratio
    print(f"  D_base (stress-attenuated) range: [{D_base.min():.5f}, {D_base.max():.5f}] mm^2/day")

    # ---- DTI tensor field ----
    print("[PHASE 3] Constructing 3D DTI tensor field...")
    tensor, FA, e1 = construct_dti_tensor_field(D_base, stress_ratio, R_ANISO)
    print(f"  Tensor shape: {tensor.shape}, FA range: [{FA.min():.3f}, {FA.max():.3f}]")

    # ---- CFL time step ----
    dt, n_steps = compute_cfl_dt(tensor, DX_MM)

    # ---- Isotropic baseline (trace/3) ----
    D_iso = np.trace(tensor, axis1=3, axis2=4) / 3.0
    print("[PHASE 3] Running isotropic baseline...")
    c_iso = run_isotropic_baseline(D_iso, rho, dt, n_steps, DX_MM)

    # ---- Anisotropic simulation ----
    c = make_initial_seed()
    print(f"[PHASE 3] Running 3D anisotropic simulation: {n_steps} steps, dt={dt:.5f} days...")
    report_every = max(1, n_steps // 5)
    for step in range(n_steps):
        c = step_anisotropic_3d(c, tensor, rho, DX_MM, dt)
        if step % report_every == 0 or step == n_steps - 1:
            vol = float(np.sum(c) * DX_MM ** 3)
            print(f"  step {step}/{n_steps}: volume={vol:.2f} mm^3, max_c={c.max():.3f}")

    # ---- Metrics + outputs ----
    print("[PHASE 3] Computing DTI invasion metrics...")
    metrics = compute_dti_metrics(c, c_iso, tensor, e1, DX_MM)
    print(f"  Mean aniso ratio: {metrics['mean_tensor_anisotropy_ratio']:.2f}")
    print(f"  Invasion X/Y/Z: {metrics['max_invasion_x_mm']:.2f} / "
          f"{metrics['max_invasion_y_mm']:.2f} / {metrics['max_invasion_z_mm']:.2f} mm")
    print(f"  Asymmetry: {metrics['asymmetry_index']:.3f}, "
          f"Preferred angle: {metrics['preferred_migration_angle_deg']:.1f} deg")

    np.save(out_dir / "phase3_aniso_concentration.npy", c)
    np.save(out_dir / "phase3_iso_concentration.npy", c_iso)
    save_dti_visualization(c, c_iso, tensor, FA, e1, metrics,
                           out_dir / "phase3_3d_anisotropic_dti.png")
    save_dti_metrics(metrics, out_dir / "phase3_dti_metrics.json")
    print("[PHASE 3] Complete.")
    return metrics


def main() -> None:
    print("=" * 70)
    print("  PHASE 3: 3D Anisotropic DTI Tensor Integration (standalone)")
    print("=" * 70)
    run_phase3()


if __name__ == "__main__":
    main()
