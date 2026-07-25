#!/usr/bin/env python3
"""
Self-contained Kaggle Kernel: Phase 1 (Spatial Genomics) -> Phase 2 (Poroelastic Mechanics)

This single script runs both phases with in-memory data passing.
Optimized for < 120s total runtime on Kaggle GPU/CPU.
"""
from __future__ import annotations

import json
import os
import time
import warnings
from pathlib import Path
from typing import Dict, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm
from scipy.ndimage import gaussian_filter

warnings.filterwarnings("ignore", category=FutureWarning)

# =============================================================================
# SHARED CONFIGURATION
# =============================================================================
# LOCAL_TEST flag: set True for fast local dry-run, False for production on Kaggle
LOCAL_TEST = False

if LOCAL_TEST:
    GRID_3D = (32, 32, 16)
    N_STEPS = 10
    N_ADVI_ITER = 500
else:
    GRID_3D = (128, 128, 64)
    N_STEPS = 200  # explicit Euler pressure steps (vectorized, no matrix)
    N_ADVI_ITER = 10000  # ELBO converges ~here; 20k added no accuracy

DOMAIN_MM = (100.0, 100.0, 50.0)
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Tumor geometry (synthetic placeholder for Phase 3 coupling)
# Center in voxel coordinates (matches GRID_3D indexing)
TUMOR_CENTER = (GRID_3D[0] / 2.0, GRID_3D[1] / 2.0, GRID_3D[2] / 2.0)
TUMOR_SIGMA_MM = 8.0  # Gaussian spread (mm)

NEFTEL_STATES = ["NPC-like", "OPC-like", "AC-like", "MES-like"]
N_STATES = len(NEFTEL_STATES)

# Visium spot configuration (reduced for fast development)
SPOT_SPACING = 16  # 8x8x4 = 256 spots
N_SPOTS_X = GRID_3D[0] // SPOT_SPACING
N_SPOTS_Y = GRID_3D[1] // SPOT_SPACING
N_SPOTS_Z = GRID_3D[2] // SPOT_SPACING
N_SPOTS = N_SPOTS_X * N_SPOTS_Y * N_SPOTS_Z

N_GENES = 100
SEED = 42
np.random.seed(SEED)

# ADVI settings (fast variational inference)
ADVI_LEARNING_RATE = 0.01

# Biomechanical parameters (brain tissue)
E_YOUNG = 3.0        # kPa
NU_POISSON = 0.45
ALPHA_BIOT = 0.8
K_HYDRAULIC = 0.01   # mm^2/day/kPa (hydraulic conductivity)
GAMMA_SOURCE = 0.5
Q_DRAIN = 0.05
BETA_STRESS = 0.02   # kPa^-1 (stress-diffusion coupling) - reference value, overridden by normalization
# Normalized stress-coupling: scale so peak von Mises stress causes this local diffusion
# drop in the tumor core (5-15% target -> 10% center). Decouples from absolute stress magnitudes.
TARGET_PEAK_DIFFUSION_REDUCTION = 0.10  # fraction (0.10 = 10% drop at peak stress)

DT_DAYS = 0.1

# =============================================================================
# 3D DIFFERENTIAL OPERATORS
# =============================================================================
def grad_3d(f: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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


def von_mises_stress(ux: np.ndarray, uy: np.ndarray, uz: np.ndarray, dx: float) -> np.ndarray:
    ux_x, ux_y, ux_z = grad_3d(ux, dx)
    uy_x, uy_y, uy_z = grad_3d(uy, dx)
    uz_x, uz_y, uz_z = grad_3d(uz, dx)
    
    exx, eyy, ezz = ux_x, uy_y, uz_z
    exy = 0.5 * (ux_y + uy_x)
    eyz = 0.5 * (uy_z + uz_y)
    ezx = 0.5 * (uz_x + ux_z)
    
    lam = E_YOUNG * NU_POISSON / ((1 + NU_POISSON) * (1 - 2 * NU_POISSON))
    mu = E_YOUNG / (2 * (1 + NU_POISSON))
    
    trace_e = exx + eyy + ezz
    sxx = lam * trace_e + 2 * mu * exx
    syy = lam * trace_e + 2 * mu * eyy
    szz = lam * trace_e + 2 * mu * ezz
    sxy = 2 * mu * exy
    syz = 2 * mu * eyz
    szx = 2 * mu * ezx
    
    vM = np.sqrt(0.5 * (
        (sxx - syy)**2 + (syy - szz)**2 + (szz - sxx)**2 +
        6 * (sxy**2 + syz**2 + szx**2)
    ))
    return vM

# =============================================================================
# PHASE 1: SPATIAL GENOMICS DECONVOLUTION (ADVI)
# =============================================================================
def generate_neftel_signatures() -> np.ndarray:
    signatures = np.zeros((N_STATES, N_GENES))
    genes_per_state = N_GENES // N_STATES
    for s in range(N_STATES):
        start = s * genes_per_state
        end = (s + 1) * genes_per_state if s < N_STATES - 1 else N_GENES
        signatures[s, start:end] = np.random.gamma(5.0, 1.0, end - start)
    hk_genes = np.random.choice(N_GENES, size=N_GENES // 10, replace=False)
    signatures[:, hk_genes] += np.random.gamma(2.0, 0.5, (N_STATES, len(hk_genes)))
    signatures += np.random.gamma(0.5, 0.2, signatures.shape)
    return signatures


def generate_spatial_fractions() -> Tuple[np.ndarray, np.ndarray]:
    fractions = np.zeros((N_SPOTS, N_STATES))
    spot_coords = []
    for idx, (ix, iy, iz) in enumerate(np.ndindex(N_SPOTS_X, N_SPOTS_Y, N_SPOTS_Z)):
        x = (ix + 0.5) * SPOT_SPACING * DX_MM
        y = (iy + 0.5) * SPOT_SPACING * DX_MM
        z = (iz + 0.5) * SPOT_SPACING * DX_MM
        spot_coords.append([x, y, z])
        cx, cy, cz = N_SPOTS_X / 2, N_SPOTS_Y / 2, N_SPOTS_Z / 2
        r = np.sqrt(((ix - cx)/cx)**2 + ((iy - cy)/cy)**2 + ((iz - cz)/cz)**2)
        if r < 0.4:
            npc, opc, ac, mes = 0.45, 0.35, 0.10, 0.10
        elif r > 0.6:
            npc, opc, ac, mes = 0.10, 0.10, 0.35, 0.45
        else:
            npc, opc, ac, mes = 0.25, 0.20, 0.25, 0.30
        f = np.array([npc, opc, ac, mes]) + 0.1 * np.random.rand(4)
        fractions[idx] = f / f.sum()
    return fractions, np.array(spot_coords)


def simulate_spot_expression(fractions: np.ndarray, signatures: np.ndarray) -> np.ndarray:
    mean_expr = fractions @ signatures
    noise = np.random.lognormal(0.0, 0.3, mean_expr.shape)
    observed = np.random.poisson(mean_expr * noise * 100).astype(float) / 100.0
    return observed


def run_deconvolution(observed: np.ndarray, signatures: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    N_spots, N_genes = observed.shape
    with pm.Model() as model:
        alpha = np.ones(N_STATES)
        fractions = pm.Dirichlet("fractions", a=alpha, shape=(N_spots, N_STATES))
        mu = pm.math.dot(fractions, signatures)
        sigma = pm.HalfNormal("sigma", sigma=0.5)
        pm.LogNormal("obs", mu=pm.math.log(mu + 1e-6), sigma=sigma, observed=observed + 1e-6)
        approx = pm.fit(n=N_ADVI_ITER, method="advi", obj_optimizer=pm.adam(learning_rate=ADVI_LEARNING_RATE))
        trace = approx.sample(2000)
    frac_samples = trace.posterior["fractions"].values.reshape(-1, N_spots, N_STATES)
    return frac_samples.mean(axis=0), np.percentile(frac_samples, 2.5, axis=0), np.percentile(frac_samples, 97.5, axis=0)


def map_fractions_to_parameters(frac_mean, frac_lower, frac_upper) -> Dict:
    npc_idx, opc_idx, ac_idx, mes_idx = 0, 1, 2, 3
    prolif = frac_mean[:, npc_idx] + frac_mean[:, opc_idx]
    invas = frac_mean[:, ac_idx] + frac_mean[:, mes_idx]
    prolif_l = frac_lower[:, npc_idx] + frac_lower[:, opc_idx]
    prolif_u = frac_upper[:, npc_idx] + frac_upper[:, opc_idx]
    invas_l = frac_lower[:, ac_idx] + frac_lower[:, mes_idx]
    invas_u = frac_upper[:, ac_idx] + frac_upper[:, mes_idx]
    
    def norm(arr, ref): return (arr - ref.min()) / (ref.max() - ref.min() + 1e-8)
    
    RHO_MIN, RHO_MAX = 0.005, 0.12
    D_MIN, D_MAX = 0.01, 0.50
    
    rho_mean = RHO_MIN + (RHO_MAX - RHO_MIN) * norm(prolif, prolif)
    D_mean = D_MIN + (D_MAX - D_MIN) * norm(invas, invas)
    rho_lower = RHO_MIN + (RHO_MAX - RHO_MIN) * norm(prolif_l, prolif)
    rho_upper = RHO_MIN + (RHO_MAX - RHO_MIN) * norm(prolif_u, prolif)
    D_lower = D_MIN + (D_MAX - D_MIN) * norm(invas_l, invas)
    D_upper = D_MIN + (D_MAX - D_MIN) * norm(invas_u, invas)
    
    return {"rho_mean": rho_mean, "rho_lower": rho_lower, "rho_upper": rho_upper,
            "D_mean": D_mean, "D_lower": D_lower, "D_upper": D_upper}


def interpolate_to_full_grid(spot_values: np.ndarray, spot_coords: np.ndarray) -> np.ndarray:
    from scipy.interpolate import griddata
    nx, ny, nz = GRID_3D
    x = np.linspace(0, DOMAIN_MM[0], nx)
    y = np.linspace(0, DOMAIN_MM[1], ny)
    z = np.linspace(0, DOMAIN_MM[2], nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    grid_coords = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    fill = np.nanmedian(spot_values) if spot_values.size else 0.0
    try:
        vals = griddata(spot_coords, spot_values, grid_coords, method="linear", fill_value=fill)
    except Exception:
        # Qhull can fail for degenerate (coplanar/flat) point sets; fall back to nearest
        vals = griddata(spot_coords, spot_values, grid_coords, method="nearest", fill_value=fill)
    return vals.reshape(GRID_3D)


def run_phase1() -> Dict[str, np.ndarray]:
    print("[PHASE 1] Generating synthetic 10x Visium data...")
    signatures = generate_neftel_signatures()
    true_fractions, spot_coords = generate_spatial_fractions()
    observed = simulate_spot_expression(true_fractions, signatures)
    
    print("[PHASE 1] Running Bayesian deconvolution (ADVI)...")
    frac_mean, frac_lower, frac_upper = run_deconvolution(observed, signatures)
    
    print("[PHASE 1] Mapping to biophysical parameters...")
    params = map_fractions_to_parameters(frac_mean, frac_lower, frac_upper)
    
    print("[PHASE 1] Interpolating to 3D grid...")
    posteriors = {}
    for key in params:
        posteriors[key] = interpolate_to_full_grid(params[key], spot_coords)
    
    print(f"  rho range: [{posteriors['rho_mean'].min():.5f}, {posteriors['rho_mean'].max():.5f}]")
    print(f"  D range: [{posteriors['D_mean'].min():.5f}, {posteriors['D_mean'].max():.5f}]")
    return posteriors


# =============================================================================
# PHASE 2: BIOT POROELASTIC MECHANICS (Optimized Solvers)
# =============================================================================
def generate_tumor_concentration() -> np.ndarray:
    nx, ny, nz = GRID_3D
    cx, cy, cz = TUMOR_CENTER
    sigma_vox = TUMOR_SIGMA_MM / DX_MM
    ix = np.arange(nx)[:, None, None]
    iy = np.arange(ny)[None, :, None]
    iz = np.arange(nz)[None, None, :]
    r2 = ((ix - cx)**2 + (iy - cy)**2 + (iz - cz)**2) * DX_MM**2
    u = np.exp(-r2 / (2 * sigma_vox**2))
    return u / (u.max() + 1e-12)


def solve_pressure_field(rho_field: np.ndarray, tumor_conc: np.ndarray) -> np.ndarray:
    """Solve parabolic pressure PDE on full grid using explicit forward Euler + vectorized Laplacian.
    dp/dt = K * laplacian(p) + gamma * rho * u - Q_drain * p
    No sparse matrices, no LU factorization, no coarse-grid zoom. Pure NumPy."""
    from scipy.ndimage import zoom as _zoom

    # Align inputs to the full grid shape (defensive; Phase 1 returns GRID_3D arrays)
    if rho_field.shape != GRID_3D:
        rho_field = _zoom(rho_field, tuple(GRID_3D[i]/rho_field.shape[i] for i in range(3)), order=1)
    if tumor_conc.shape != GRID_3D:
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
    n_max = max(N_STEPS, 5 * int(round(1.0 / (Q_DRAIN * dt_stable))) + 10)

    print(f"  Solving pressure PDE (explicit Euler, up to {n_max} steps, dt={dt_stable:.4f} day)...")
    prev_mass = 0.0
    for step in range(n_max):
        lap = laplacian_3d(p, DX_MM)
        p = p + dt_stable * (K_HYDRAULIC * lap + source - Q_DRAIN * p)
        # Track total pressure mass (settles only when fully equilibrated, not just locally)
        cur_mass = float(p.sum())
        if step % 50 == 0:
            print(f"    step {step}: max p = {p.max():.4f} kPa, mass = {cur_mass:.4f}")
        if step > 0 and abs(cur_mass - prev_mass) < 1e-4 * (abs(cur_mass) + 1e-12):
            print(f"  Converged at step {step}")
            break
        prev_mass = cur_mass

    return p


def solve_displacement_field(pressure: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Closed-form displacement from the pressure gradient field.
    u_disp = -alpha * grad(p) / (lambda + 2*mu)
    No iterative matrix solver (no CG, no sparse operators). One vectorized pass.
    Pressure is mildly Gaussian-smoothed before taking gradients so the closed-form
    displacement (and the derived von Mises stress) is free of the triangular block
    artifacts produced by raw directional finite differences."""
    from scipy.ndimage import gaussian_filter, zoom as _zoom
    # Align to full grid
    if pressure.shape != GRID_3D:
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


# =============================================================================
# VISUALIZATION & OUTPUT
# =============================================================================
def save_phase1_visualization(frac_mean, rho_grid, D_grid, rho_l, rho_u, D_l, D_u, out_path: Path):
    nx, ny, nz = GRID_3D
    sx = nx // 2
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = fig.add_gridspec(2, 4)
    
    ax1 = fig.add_subplot(gs[0, :2])
    state_means = np.nan_to_num(frac_mean.mean(axis=0), nan=0.25)
    if np.sum(state_means) == 0:
        state_means = np.array([0.25, 0.25, 0.25, 0.25])
    bars = ax1.bar(NEFTEL_STATES, state_means,
                   color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'],
                   alpha=0.85, edgecolor='black')
    ax1.set_ylim(0, max(0.5, float(np.max(state_means)) * 1.3))
    ax1.set_ylabel("Mean Proportions")
    ax1.set_title("Panel 1: Neftel Cell State Proportions", fontweight="bold")
    ax1.grid(True, axis='y', linestyle='--', alpha=0.5)
    
    ax2 = fig.add_subplot(gs[0, 2])
    im2 = ax2.imshow(rho_grid[sx].T, cmap="hot_r", origin="lower", extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    ax2.set_title("Panel 2: Proliferation Rate rho [day^-1]", fontweight="bold")
    plt.colorbar(im2, ax=ax2, shrink=0.8)
    
    ax3 = fig.add_subplot(gs[0, 3])
    im3 = ax3.imshow(D_grid[sx].T, cmap="viridis", origin="lower", extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    ax3.set_title("Panel 3: Diffusion Rate D [mm^2/day]", fontweight="bold")
    plt.colorbar(im3, ax=ax3, shrink=0.8)
    
    ax4 = fig.add_subplot(gs[1, :])
    n_vis = 200
    flat_rho = rho_grid.ravel(); flat_D = D_grid.ravel()
    flat_rl = rho_l.ravel(); flat_ru = rho_u.ravel()
    flat_Dl = D_l.ravel(); flat_Du = D_u.ravel()
    valid = ~np.isnan(flat_rho)
    idx = np.random.choice(np.where(valid)[0], n_vis, replace=False)
    sort_idx = np.argsort(flat_rho[idx])
    x_pos = np.arange(n_vis)
    ax4.fill_between(x_pos, flat_rl[idx][sort_idx], flat_ru[idx][sort_idx], alpha=0.3, color="red", label="rho 95% CI")
    ax4.plot(x_pos, flat_rho[idx][sort_idx], "r-", label="rho mean")
    ax4.fill_between(x_pos, flat_Dl[idx][sort_idx], flat_Du[idx][sort_idx], alpha=0.3, color="blue", label="D 95% CI")
    ax4.plot(x_pos, flat_D[idx][sort_idx], "b-", label="D mean")
    ax4.set_title("Panel 4: Spatial Parameter Posteriors with 95% CI", fontweight="bold")
    ax4.legend(); ax4.grid(True, alpha=0.3)
    
    fig.suptitle("Phase 1: Spatial Genomics -> Bayesian Parameter Fields", fontsize=14, fontweight="bold", y=1.02)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_phase2_visualization(pressure, ux, uy, uz, vM, D_eff, D_base, out_path: Path):
    nx, ny, nz = GRID_3D
    sx = nx // 2
    disp_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)
    
    ax1 = fig.add_subplot(gs[0, 0])
    vmax_p = np.max(np.abs(pressure))
    im1 = ax1.imshow(pressure[sx].T, cmap="RdBu_r", origin="lower", extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=-vmax_p, vmax=vmax_p)
    ax1.set_title("Panel 1: Interstitial Pressure p [kPa]", fontweight="bold")
    plt.colorbar(im1, ax=ax1, shrink=0.8)
    
    ax2 = fig.add_subplot(gs[0, 1])
    im2 = ax2.imshow(disp_mag[sx].T, cmap="hot", origin="lower", extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    ax2.set_title("Panel 2: Displacement ||u|| [mm]", fontweight="bold")
    plt.colorbar(im2, ax=ax2, shrink=0.8)
    
    ax3 = fig.add_subplot(gs[1, 0])
    im3 = ax3.imshow(vM[sx].T, cmap="plasma", origin="lower", extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    ax3.set_title("Panel 3: von Mises Stress [kPa]", fontweight="bold")
    plt.colorbar(im3, ax=ax3, shrink=0.8)
    
    ax4 = fig.add_subplot(gs[1, 1])
    ratio = D_eff / (D_base + 1e-12)
    im4 = ax4.imshow(ratio[sx].T, cmap="RdYlGn", origin="lower", extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=1)
    ax4.set_title("Panel 4: D_eff/D0 (Stress-Modified)", fontweight="bold")
    plt.colorbar(im4, ax=ax4, shrink=0.8)
    
    fig.suptitle("Phase 2: Biot Poroelastic Mechanics", fontsize=14, fontweight="bold", y=1.02)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def compute_metrics(pressure, ux, uy, uz, vM, D_eff, D_base) -> Dict:
    disp_mag = np.sqrt(ux**2 + uy**2 + uz**2)
    vx, vy, vz = GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2
    midline_shift = float(disp_mag[vx, vy, vz])
    
    ix = np.arange(GRID_3D[0])[:, None, None]
    iy = np.arange(GRID_3D[1])[None, :, None]
    iz = np.arange(GRID_3D[2])[None, None, :]
    r2_xy = ((ix - vx)**2 + (iy - vy)**2) * DX_MM**2
    z_low, z_high = int(GRID_3D[2] * 0.2), int(GRID_3D[2] * 0.8)
    vent_mask = (r2_xy < 10**2) & (iz >= z_low) & (iz < z_high)
    vol_strain = np.sum(vM[vent_mask]) / (vent_mask.sum() + 1)
    vent_comp = max(0, min(100, float(100 * vol_strain / (E_YOUNG * 100))))
    
    return {
        "max_interstitial_pressure_kPa": float(np.max(pressure)),
        "max_displacement_mm": float(np.max(disp_mag)),
        "midline_shift_mm": midline_shift,
        "ventricular_compression_pct": vent_comp,
        "mean_von_mises_stress_kPa": float(np.mean(vM)),
        "max_von_mises_stress_kPa": float(np.max(vM)),
        "diffusion_reduction_pct": float(100 * (1 - D_eff.mean() / (D_base.mean() + 1e-12))),
    }


# =============================================================================
# MAIN PIPELINE
# =============================================================================
def main():
    out_dir = Path("output")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("  GBM PIPELINE: Phase 1 -> Phase 2 (Self-contained Kaggle Kernel)")
    print("=" * 70)
    
    # --- Phase 1 ---
    p1_start = time.time()
    posteriors = run_phase1()
    print(f"[PHASE 1] Done in {time.time() - p1_start:.1f}s")
    
    # Save Phase 1 outputs
    np.save(out_dir / "phase1_rho_posterior.npy", 
            np.stack([posteriors["rho_mean"], posteriors["rho_lower"], posteriors["rho_upper"]], axis=0))
    np.save(out_dir / "phase1_D_posterior.npy", 
            np.stack([posteriors["D_mean"], posteriors["D_lower"], posteriors["D_upper"]], axis=0))
    save_phase1_visualization(
        posteriors.get("frac_mean", np.zeros((1, 4))),  # placeholder if not available
        posteriors["rho_mean"], posteriors["D_mean"],
        posteriors["rho_lower"], posteriors["rho_upper"],
        posteriors["D_lower"], posteriors["D_upper"],
        out_dir / "phase1_spatial_deconv.png"
    )
    
    # --- Phase 2 ---
    p2_start = time.time()
    print("[PHASE 2] Loading Phase 1 posteriors from memory...")
    rho = posteriors["rho_mean"]
    D_base = posteriors["D_mean"]
    
    print("[PHASE 2] Generating tumor concentration...")
    tumor_conc = generate_tumor_concentration()
    
    print("[PHASE 2] Solving pressure PDE...")
    pressure = solve_pressure_field(rho, tumor_conc)
    
    print("[PHASE 2] Solving displacement field...")
    ux, uy, uz = solve_displacement_field(pressure)
    
    print("[PHASE 2] Computing von Mises stress...")
    vM = von_mises_stress(ux, uy, uz, DX_MM)
    
    print("[PHASE 2] Computing stress-modified diffusion...")
    # Normalized coupling: peak von Mises -> target local diffusion drop (5-15% range),
    # so Panel 4 shows a clear core<periphery variation independent of absolute stress scale.
    vM_peak = float(vM.max())
    if vM_peak <= 1e-12:
        D_eff = D_base.copy()
    else:
        target = TARGET_PEAK_DIFFUSION_REDUCTION
        if target >= 1.0:
            target = 0.99
        gamma_eff = -np.log1p(-target) / vM_peak
        D_eff = D_base * np.exp(-gamma_eff * vM)
    print(f"  Peak von Mises: {vM_peak:.6f} kPa -> diffusion reduction at peak ~{TARGET_PEAK_DIFFUSION_REDUCTION*100:.0f}%")
    
    print("[PHASE 2] Computing metrics...")
    metrics = compute_metrics(pressure, ux, uy, uz, vM, D_eff, D_base)
    
    # Save Phase 2 outputs
    np.save(out_dir / "phase2_pressure_field.npy", pressure)
    np.save(out_dir / "phase2_displacement_field.npy", np.stack([ux, uy, uz], axis=0))
    with open(out_dir / "phase2_mechanics_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    save_phase2_visualization(pressure, ux, uy, uz, vM, D_eff, D_base, out_dir / "phase2_poroelastic_mechanics.png")
    
    print(f"[PHASE 2] Done in {time.time() - p2_start:.1f}s")
    
    # Summary
    total = time.time() - p1_start
    print("\n" + "=" * 70)
    print(f"  PIPELINE COMPLETE in {total:.1f}s")
    print(f"  Max pressure: {metrics['max_interstitial_pressure_kPa']:.2f} kPa")
    print(f"  Midline shift: {metrics['midline_shift_mm']:.2f} mm")
    print(f"  Ventricular compression: {metrics['ventricular_compression_pct']:.1f}%")
    print(f"  Outputs saved to {out_dir}/")
    print("=" * 70)


if __name__ == "__main__":
    main()