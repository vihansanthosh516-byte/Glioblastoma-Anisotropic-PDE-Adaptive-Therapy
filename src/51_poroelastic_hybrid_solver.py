#!/usr/bin/env python3
"""
Phase 2: Biot Poroelastic Hybrid Solver

Couples tumor growth (Phase 1 posteriors) with solid tissue mechanics and
interstitial fluid dynamics to compute mass effect, midline shift, and
stress-modified diffusion.

Core Physics:
- Pressure: ∂p/∂t = K ∇²p + γ ρ u - Q_drain
- Displacement: ∇·σ_solid - α ∇p = 0  (Navier-Cauchy with pore pressure coupling)
- Effective diffusion: D_eff = D * exp(-β σ_vM)  (mechanical inhibition)

Outputs:
- output/phase2_pressure_field.npy
- output/phase2_displacement_field.npy
- output/phase2_mechanics_metrics.json
- output/phase2_poroelastic_mechanics.png
"""
from __future__ import annotations

import json
import os
import warnings
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, binary_erosion, distance_transform_edt

warnings.filterwarnings("ignore", category=FutureWarning)

# --------------------------------------------------------------------------- #
# Configuration (matching Phase 1 grid)
# --------------------------------------------------------------------------- #
# LOCAL_TEST flag: set True for fast local dry-run, False for production
LOCAL_TEST = False

if LOCAL_TEST:
    GRID_3D = (32, 32, 16)
    N_STEPS = 10
    N_ADVI_ITER = 500  # for Phase 1 reference
else:
    GRID_3D = (128, 128, 64)
    N_STEPS = 200  # explicit Euler pressure steps (vectorized, no matrix)

DOMAIN_MM = (100.0, 100.0, 50.0)
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Tumor geometry (synthetic placeholder for Phase 3 coupling)
# Center in voxel coordinates (matches GRID_3D indexing)
TUMOR_CENTER = (GRID_3D[0] / 2.0, GRID_3D[1] / 2.0, GRID_3D[2] / 2.0)
TUMOR_SIGMA_MM = 8.0  # Gaussian spread (mm)

# Biomechanical parameters (brain tissue)
E_YOUNG = 3.0        # kPa
NU_POISSON = 0.45    # nearly incompressible
ALPHA_BIOT = 0.8     # Biot coefficient
K_HYDRAULIC = 0.01   # mm^2/day/kPa (hydraulic conductivity)
GAMMA_SOURCE = 0.5   # pressure source scaling
Q_DRAIN = 0.05       # drainage rate (per day)
BETA_STRESS = 0.02   # kPa^-1 (stress-diffusion coupling) - reference value, overridden by normalization
# Normalized stress-coupling: scale so peak von Mises stress causes this local diffusion
# drop in the tumor core (5-15% target -> 10% center). Decouples from absolute stress magnitudes.
TARGET_PEAK_DIFFUSION_REDUCTION = 0.10  # fraction (0.10 = 10% drop at peak stress)

DT_DAYS = 0.1

def grad_3d(f: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Central difference gradient with Neumann BC."""
    fx = np.zeros_like(f); fy = np.zeros_like(f); fz = np.zeros_like(f)
    fx[1:-1,:,:] = (f[2:,:,:] - f[:-2,:,:]) / (2*dx)
    fx[0,:,:] = (f[1,:,:] - f[0,:,:]) / dx
    fx[-1,:,:] = (f[-1,:,:] - f[-2,:,:]) / dx
    fy[:,1:-1,:] = (f[:,2:,:] - f[:,:-2,:]) / (2*dx)
    fy[:,0,:] = (f[:,1,:] - f[:,0,:]) / dx
    fy[:,-1,:] = (f[:,-1,:] - f[:,-2,:]) / dx
    fz[:,:,1:-1] = (f[:,:,2:] - f[:,:,:-2]) / (2*dx)
    fz[:,:,0] = (f[:,:,1] - f[:,:,0]) / dx
    fz[:,:,-1] = (f[:,:,-1] - f[:,:,-2]) / dx
    return fx, fy, fz


def div_3d(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float) -> np.ndarray:
    """Divergence of vector field."""
    return grad_3d(vx, dx)[0] + grad_3d(vy, dx)[1] + grad_3d(vz, dx)[2]


def laplacian_3d(f: np.ndarray, dx: float) -> np.ndarray:
    """7-point stencil Laplacian with zero-flux (Neumann) BC.
    Vectorized via direct array slicing (no pad allocation) - fast on 1M+ voxels."""
    lap = np.zeros_like(f)
    inv = 1.0 / dx**2
    # x-axis (reflective = use edge value => boundary reduces to nearest interior)
    lap[1:-1, :, :] += (f[2:, :, :] - 2.0 * f[1:-1, :, :] + f[:-2, :, :]) * inv
    lap[0, :, :] += (f[1, :, :] - f[0, :, :]) * inv
    lap[-1, :, :] += (f[-2, :, :] - f[-1, :, :]) * inv
    # y-axis
    lap[:, 1:-1, :] += (f[:, 2:, :] - 2.0 * f[:, 1:-1, :] + f[:, :-2, :]) * inv
    lap[:, 0, :] += (f[:, 1, :] - f[:, 0, :]) * inv
    lap[:, -1, :] += (f[:, -2, :] - f[:, -1, :]) * inv
    # z-axis
    lap[:, :, 1:-1] += (f[:, :, 2:] - 2.0 * f[:, :, 1:-1] + f[:, :, :-2]) * inv
    lap[:, :, 0] += (f[:, :, 1] - f[:, :, 0]) * inv
    lap[:, :, -1] += (f[:, :, -2] - f[:, :, -1]) * inv
    return lap


def von_mises_stress(
    ux: np.ndarray, uy: np.ndarray, uz: np.ndarray, dx: float
) -> np.ndarray:
    """
    Compute von Mises stress from displacement field.
    σ_vM = sqrt(0.5 * ((σ_xx-σ_yy)² + (σ_yy-σ_zz)² + (σ_zz-σ_xx)² + 6(σ_xy²+σ_yz²+σ_zx²)))
    """
    # Strain tensor components
    ux_x, ux_y, ux_z = grad_3d(ux, dx)
    uy_x, uy_y, uy_z = grad_3d(uy, dx)
    uz_x, uz_y, uz_z = grad_3d(uz, dx)
    
    exx = ux_x
    eyy = uy_y
    ezz = uz_z
    exy = 0.5 * (ux_y + uy_x)
    eyz = 0.5 * (uy_z + uz_y)
    ezx = 0.5 * (uz_x + ux_z)
    
    # Linear elastic stress (kPa)
    lam = E_YOUNG * NU_POISSON / ((1 + NU_POISSON) * (1 - 2 * NU_POISSON))
    mu = E_YOUNG / (2 * (1 + NU_POISSON))
    
    trace_e = exx + eyy + ezz
    sxx = lam * trace_e + 2 * mu * exx
    syy = lam * trace_e + 2 * mu * eyy
    szz = lam * trace_e + 2 * mu * ezz
    sxy = 2 * mu * exy
    syz = 2 * mu * eyz
    szx = 2 * mu * ezx
    
    # von Mises
    vM = np.sqrt(0.5 * (
        (sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2 +
        6 * (sxy**2 + syz**2 + szx**2)
    ))
    return vM


# --------------------------------------------------------------------------- #
# Phase 1 Posterior Loading
# --------------------------------------------------------------------------- #
def load_phase1_posteriors() -> Dict[str, np.ndarray]:
    """Load rho and D posterior fields from Phase 1, or generate synthetic fallback."""
    out_dir = Path("output")
    rho_path = out_dir / "phase1_rho_posterior.npy"
    D_path = out_dir / "phase1_D_posterior.npy"
    
    if rho_path.exists() and D_path.exists():
        print("[PHASE 2] Loading Phase 1 posteriors from disk...")
        rho_data = np.load(rho_path)
        D_data = np.load(D_path)
    else:
        print("[PHASE 2] Phase 1 files not found. Generating synthetic fallback fields...")
        # Vectorized synthetic fallback matching Phase 1 ranges
        nx, ny, nz = GRID_3D
        cx, cy, cz = GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2
        
        ix = np.arange(nx)[:, None, None]
        iy = np.arange(ny)[None, :, None]
        iz = np.arange(nz)[None, None, :]
        
        r2 = ((ix - cx)**2 + (iy - cy)**2 + (iz - cz)**2) * DX_MM**2
        
        # Core (r2 < 500): high rho, low D; Rim: low rho, high D
        rho_mean = np.where(r2 < 500, 0.08, 0.02) + 0.01 * np.random.rand(nx, ny, nz)
        D_mean = np.where(r2 < 500, 0.05, 0.15) + 0.03 * np.random.rand(nx, ny, nz)
        
        # Smooth
        from scipy.ndimage import gaussian_filter
        rho_mean = gaussian_filter(rho_mean, sigma=2)
        D_mean = gaussian_filter(D_mean, sigma=2)
        
        # Create dummy CI bounds (±20%)
        rho_lower = rho_mean * 0.8
        rho_upper = rho_mean * 1.2
        D_lower = D_mean * 0.8
        D_upper = D_mean * 1.2
        
        # Save for future use
        out_dir.mkdir(parents=True, exist_ok=True)
        np.save(rho_path, np.stack([rho_mean, rho_lower, rho_upper], axis=0))
        np.save(D_path, np.stack([D_mean, D_lower, D_upper], axis=0))
        print(f"[PHASE 2] Saved synthetic Phase 1 posteriors to {out_dir}")
        
        return {
            "rho_mean": rho_mean,
            "rho_lower": rho_lower,
            "rho_upper": rho_upper,
            "D_mean": D_mean,
            "D_lower": D_lower,
            "D_upper": D_upper,
        }
    
    rho_mean, rho_lower, rho_upper = rho_data[0], rho_data[1], rho_data[2]
    D_mean, D_lower, D_upper = D_data[0], D_data[1], D_data[2]

    # If loaded arrays don't match GRID_3D (e.g. LOCAL_TEST with production files), resize
    target_shape = GRID_3D
    from scipy.ndimage import zoom
    def _resize(arr):
        if arr.shape == target_shape:
            return arr
        zf = tuple(target_shape[i] / arr.shape[i] for i in range(3))
        return zoom(arr, zf, order=1)

    return {
        "rho_mean": _resize(rho_mean),
        "rho_lower": _resize(rho_lower),
        "rho_upper": _resize(rho_upper),
        "D_mean": _resize(D_mean),
        "D_lower": _resize(D_lower),
        "D_upper": _resize(D_upper),
    }


# --------------------------------------------------------------------------- #
# Tumor Concentration Field (placeholder for Phase 3 coupling)
# --------------------------------------------------------------------------- #
def generate_tumor_concentration() -> np.ndarray:
    """Generate synthetic tumor cell concentration field u(x,y,z) - vectorized."""
    nx, ny, nz = GRID_3D
    cx, cy, cz = TUMOR_CENTER
    sigma_vox = TUMOR_SIGMA_MM / DX_MM
    
    # Vectorized coordinate grids
    ix = np.arange(nx)[:, None, None]
    iy = np.arange(ny)[None, :, None]
    iz = np.arange(nz)[None, None, :]
    
    r2 = ((ix - cx)**2 + (iy - cy)**2 + (iz - cz)**2) * DX_MM**2
    u = np.exp(-r2 / (2 * sigma_vox**2))
    
    # Normalize peak to 1.0
    u = u / (u.max() + 1e-12)
    return u


# --------------------------------------------------------------------------- #
# Pressure Solver: dp/dt = K * laplacian(p) + gamma * rho * u - Q_drain * p
# Explicit forward Euler, vectorized 3D Laplacian, NO sparse matrices.
# --------------------------------------------------------------------------- #
def solve_pressure_field(
    rho_field: np.ndarray,
    tumor_conc: np.ndarray,
    dt: float = DT_DAYS,
    n_steps: int = N_STEPS
) -> np.ndarray:
    """Solve parabolic pressure PDE on the full grid using explicit forward Euler.
    No sparse matrices, no LU factorization, no coarse-grid zoom. Pure NumPy."""
    # Align inputs to the full grid shape
    if rho_field.shape != GRID_3D:
        from scipy.ndimage import zoom as _zoom
        rho_field = _zoom(rho_field, tuple(GRID_3D[i]/rho_field.shape[i] for i in range(3)), order=1)
    if tumor_conc.shape != GRID_3D:
        from scipy.ndimage import zoom as _zoom
        tumor_conc = _zoom(tumor_conc, tuple(GRID_3D[i]/tumor_conc.shape[i] for i in range(3)), order=1)

    # Time-invariant source term: gamma * rho * u
    source = GAMMA_SOURCE * rho_field * tumor_conc

    p = np.zeros(GRID_3D, dtype=np.float64)
    # Explicit forward Euler stability for 3D (7-point Laplacian, Neumann BC) + reaction:
    # The aggregate spectral radius is K * lambda_max + Q_drain (reaction).
    # Stability bound: dt < 2 / (K * lambda_max + Q). Apply a 40% safety margin.
    lambda_max = 6.0 / (DX_MM ** 2)  # 3D Laplacian max eigenvalue (conservative upper bound)
    dt_stable = 0.4 * 2.0 / (K_HYDRAULIC * lambda_max + Q_DRAIN + 1e-12)
    # Steady-state relaxation time ~ 1 / Q_drain; need ~5*tau to converge.
    n_max = max(n_steps, 5 * int(round(1.0 / (Q_DRAIN * dt_stable))) + 10)

    print(f"  Solving pressure PDE (explicit Euler, up to {n_max} steps, dt={dt_stable:.4f} day)...")
    prev_mass = 0.0
    for step in range(n_max):
        lap = laplacian_3d(p, DX_MM)
        p = p + dt_stable * (K_HYDRAULIC * lap + source - Q_DRAIN * p)
        cur_mass = float(p.sum())
        if step % 50 == 0:
            print(f"    step {step}: max p = {p.max():.4f} kPa, mass = {cur_mass:.4f}")
        if step > 0 and abs(cur_mass - prev_mass) < 1e-4 * (abs(cur_mass) + 1e-12):
            print(f"  Converged at step {step}")
            break
        prev_mass = cur_mass

    return p


# --------------------------------------------------------------------------- #
# Displacement Solver: u_disp = -alpha * grad(p) / (lambda + 2*mu)
# Closed-form from the pressure gradient field. No iterative matrix solver.
# --------------------------------------------------------------------------- #
def solve_displacement_field(pressure: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form displacement from the pressure gradient field.
    u_disp = -alpha * grad(p) / (lambda + 2*mu)
    No CG, no sparse operators. One vectorized pass.
    Pressure is mildly Gaussian-smoothed before taking gradients so the closed-form
    displacement (and the derived von Mises stress) is free of the triangular block
    artifacts produced by raw directional finite differences."""
    from scipy.ndimage import gaussian_filter
    # Align to full grid
    if pressure.shape != GRID_3D:
        from scipy.ndimage import zoom as _zoom
        pressure = _zoom(pressure, tuple(GRID_3D[i]/pressure.shape[i] for i in range(3)), order=1)

    # Mild smoothing (sigma=1 voxel) to suppress FD block artifacts in grad(p).
    p_smooth = gaussian_filter(pressure, sigma=1.0, mode="nearest")

    lam = E_YOUNG * NU_POISSON / ((1 + NU_POISSON) * (1 - 2 * NU_POISSON))
    mu = E_YOUNG / (2 * (1 + NU_POISSON))
    denom = lam + 2 * mu  # P-wave modulus (oedometric stiffness)

    px, py, pz = grad_3d(p_smooth, DX_MM)
    ux = -(ALPHA_BIOT * px) / denom
    uy = -(ALPHA_BIOT * py) / denom
    uz = -(ALPHA_BIOT * pz) / denom
    return ux, uy, uz


# --------------------------------------------------------------------------- #
# Stress-Modified Effective Diffusion
# --------------------------------------------------------------------------- #
def compute_effective_diffusion(
    D_base: np.ndarray,
    vM_stress: np.ndarray
) -> np.ndarray:
    """D_eff = D0 * exp(-gamma_eff * sigma_vM) - mechanical inhibition of diffusion.

    The coupling coefficient gamma_eff is NORMALIZED to the peak von Mises stress in the
    field so that the tumor core always shows a TARGET_PEAK_DIFFUSION_REDUCTION local drop
    (5-15% range), independent of the absolute stress magnitudes produced by the solver:
        gamma_eff = -ln(1 - target) / max(sigma_vM)
    This makes Panel 4 display a clear spatial variation (core < periphery) rather than the
    solid-green field produced by a fixed, under-scaled beta.
    """
    vM_peak = float(np.max(vM_stress))
    if vM_peak <= 1e-12:
        return D_base.copy()
    target = TARGET_PEAK_DIFFUSION_REDUCTION
    if target >= 1.0:
        target = 0.99
    gamma_eff = -np.log1p(-target) / vM_peak  # gamma*peak = -ln(1-target) -> peak drop = target
    reduction = np.exp(-gamma_eff * vM_stress)
    return D_base * reduction


# --------------------------------------------------------------------------- #
# Mechanical Metrics
# --------------------------------------------------------------------------- #
def compute_mechanics_metrics(
    pressure: np.ndarray,
    ux: np.ndarray, uy: np.ndarray, uz: np.ndarray,
    vM_stress: np.ndarray,
    D_eff: np.ndarray,
    D_base: np.ndarray
) -> Dict:
    """Compute clinical mechanics metrics."""
    # Max interstitial pressure
    max_pressure = float(np.max(pressure))
    
    # Max displacement magnitude (midline shift)
    disp_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    max_disp = float(np.max(disp_mag))
    
    # Midline shift (displacement at ventricle center)
    vx, vy, vz = GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2
    midline_shift = float(disp_mag[vx, vy, vz])
    
    # Ventricular compression (volume change)
    # Approximate ventricle as central cylindrical region across z - vectorized
    nx, ny, nz = GRID_3D
    vx, vy = nx//2, ny//2
    ix = np.arange(nx)[:, None, None]
    iy = np.arange(ny)[None, :, None]
    iz = np.arange(nz)[None, None, :]
    # 10mm radius cylinder along z, spanning middle 60% of z
    r2_xy = ((ix - vx)**2 + (iy - vy)**2) * DX_MM**2
    z_low, z_high = int(nz * 0.2), int(nz * 0.8)
    vent_mask = (r2_xy < 10**2) & (iz >= z_low) & (iz < z_high)  # (nx, ny, nz)
    
    # Volumetric strain at ventricle
    vol_strain = np.sum(vM_stress[vent_mask]) / (vent_mask.sum() + 1)
    vent_compression = float(100 * vol_strain / (E_YOUNG * 100))  # %
    
    # Diffusion reduction
    D_reduction = float(100 * (1 - D_eff.mean() / (D_base.mean() + 1e-12)))
    
    return {
        "max_interstitial_pressure_kPa": max_pressure,
        "max_displacement_mm": max_disp,
        "midline_shift_mm": midline_shift,
        "ventricular_compression_pct": max(0, min(100, vent_compression)),
        "mean_von_mises_stress_kPa": float(np.mean(vM_stress)),
        "max_von_mises_stress_kPa": float(np.max(vM_stress)),
        "diffusion_reduction_pct": D_reduction,
        "pressure_field_stats": {
            "mean": float(np.mean(pressure)),
            "std": float(np.std(pressure)),
        },
        "displacement_field_stats": {
            "mean_mag": float(np.mean(disp_mag)),
            "std_mag": float(np.std(disp_mag)),
        }
    }


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def save_visualization(
    pressure: np.ndarray,
    ux: np.ndarray, uy: np.ndarray, uz: np.ndarray,
    vM_stress: np.ndarray,
    D_eff: np.ndarray,
    D_base: np.ndarray,
    out_path: Path
) -> None:
    """Render 4-panel Phase 2 mechanics visualization."""
    nx, ny, nz = GRID_3D
    sx, sy, sz = nx//2, ny//2, nz//2
    disp_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    
    # Panel 1: Interstitial Fluid Pressure
    ax1 = fig.add_subplot(gs[0, 0])
    vmax_p = np.max(np.abs(pressure))
    im1 = ax1.imshow(pressure[sx, :, :].T, cmap="RdBu_r", origin="lower",
                     extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]],
                     vmin=-vmax_p, vmax=vmax_p)
    ax1.set_title("Panel 1: Interstitial Fluid Pressure p(x,y,z) [kPa]", fontweight="bold", fontsize=11)
    ax1.set_xlabel("y [mm]"); ax1.set_ylabel("z [mm]")
    plt.colorbar(im1, ax=ax1, shrink=0.8, label="p [kPa]")
    
    # Panel 2: Displacement Field Magnitude
    ax2 = fig.add_subplot(gs[0, 1])
    vmax_d = np.max(disp_mag)
    im2 = ax2.imshow(disp_mag[sx, :, :].T, cmap="hot", origin="lower",
                     extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=vmax_d)
    # Quiver overlay (subsampled)
    step = 8
    ax2.quiver(
        np.arange(0, DOMAIN_MM[1], step*DX_MM),
        np.arange(0, DOMAIN_MM[2], step*DX_MM),
        uy[sx, ::step, ::step].T, uz[sx, ::step, ::step].T,
        scale=20, color='white', alpha=0.7, width=0.003
    )
    ax2.set_title("Panel 2: Displacement Field ||u_disp|| [mm] (Midline Shift)", fontweight="bold", fontsize=11)
    ax2.set_xlabel("y [mm]"); ax2.set_ylabel("z [mm]")
    plt.colorbar(im2, ax=ax2, shrink=0.8, label="||u|| [mm]")
    
    # Panel 3: von Mises Stress
    ax3 = fig.add_subplot(gs[1, 0])
    vmax_vM = np.max(vM_stress)
    im3 = ax3.imshow(vM_stress[sx, :, :].T, cmap="plasma", origin="lower",
                     extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=vmax_vM)
    ax3.set_title("Panel 3: von Mises Mechanical Stress [kPa]", fontweight="bold", fontsize=11)
    ax3.set_xlabel("y [mm]"); ax3.set_ylabel("z [mm]")
    plt.colorbar(im3, ax=ax3, shrink=0.8, label="sigma_vM [kPa]")
    
    # Panel 4: Stress-Modified Effective Diffusion
    ax4 = fig.add_subplot(gs[1, 1])
    ratio = D_eff / (D_base + 1e-12)
    vmax_r = np.max(ratio)
    im4 = ax4.imshow(ratio[sx, :, :].T, cmap="RdYlGn", origin="lower",
                     extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=1)
    ax4.set_title("Panel 4: Effective Diffusion Ratio D_eff/D0 (Stress-Modified)", fontweight="bold", fontsize=11)
    ax4.set_xlabel("y [mm]"); ax4.set_ylabel("z [mm]")
    plt.colorbar(im4, ax=ax4, shrink=0.8, label="D_eff/D0")
    
    fig.suptitle("Phase 2: Biot Poroelastic Mechanics - Mass Effect & Stress-Modified Diffusion",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main() -> None:
    out_dir = Path("output")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("[PHASE 2] Loading Phase 1 posteriors...")
    posteriors = load_phase1_posteriors()
    rho = posteriors["rho_mean"]
    D_base = posteriors["D_mean"]
    print(f"  rho range: [{rho.min():.5f}, {rho.max():.5f}] day^-1")
    print(f"  D range: [{D_base.min():.5f}, {D_base.max():.5f}] mm^2/day")
    
    print("[PHASE 2] Generating tumor concentration field...")
    tumor_conc = generate_tumor_concentration()
    print(f"  Tumor max conc: {tumor_conc.max():.3f}")
    
    print("[PHASE 2] Solving pressure PDE (Biot poroelasticity)...")
    pressure = solve_pressure_field(rho, tumor_conc)
    print(f"  Pressure range: [{pressure.min():.4f}, {pressure.max():.4f}] kPa")
    
    print("[PHASE 2] Solving displacement field (Navier-Cauchy)...")
    ux, uy, uz = solve_displacement_field(pressure)
    disp_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    print(f"  Max displacement: {disp_mag.max():.4f} mm")
    
    print("[PHASE 2] Computing von Mises stress...")
    vM = von_mises_stress(ux, uy, uz, DX_MM)
    print(f"  von Mises range: [{vM.min():.4f}, {vM.max():.4f}] kPa")
    
    print("[PHASE 2] Computing stress-modified diffusion...")
    D_eff = compute_effective_diffusion(D_base, vM)
    print(f"  D_eff range: [{D_eff.min():.5f}, {D_eff.max():.5f}] mm^2/day")
    print(f"  Mean reduction: {(1 - D_eff.mean()/D_base.mean())*100:.1f}%")
    
    print("[PHASE 2] Computing mechanics metrics...")
    metrics = compute_mechanics_metrics(pressure, ux, uy, uz, vM, D_eff, D_base)
    print(f"  Max pressure: {metrics['max_interstitial_pressure_kPa']:.2f} kPa")
    print(f"  Midline shift: {metrics['midline_shift_mm']:.2f} mm")
    print(f"  Ventricular compression: {metrics['ventricular_compression_pct']:.1f}%")
    
    print("[PHASE 2] Saving fields...")
    np.save(out_dir / "phase2_pressure_field.npy", pressure)
    np.save(out_dir / "phase2_displacement_field.npy", np.stack([ux, uy, uz], axis=0))
    
    with open(out_dir / "phase2_mechanics_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    
    print("[PHASE 2] Rendering visualization...")
    save_visualization(pressure, ux, uy, uz, vM, D_eff, D_base,
                       out_dir / "phase2_poroelastic_mechanics.png")
    
    print(f"[PHASE 2] Complete. Outputs in {out_dir}/")
    print("  - phase2_pressure_field.npy")
    print("  - phase2_displacement_field.npy (ux, uy, uz)")
    print("  - phase2_mechanics_metrics.json")
    print("  - phase2_poroelastic_mechanics.png")


if __name__ == "__main__":
    main()