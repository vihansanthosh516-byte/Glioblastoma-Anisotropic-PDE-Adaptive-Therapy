#!/usr/bin/env python3
"""
Phase 4: Virtual Clinical Therapy Simulation (Stupp Protocol) — FIXED

Fixes applied:
1. diff/react always computed (moved outside conditional blocks)
2. Radiation applied ONCE PER DAY (not every time-step) via day-tracking
3. Anisotropic diffusion sign verified: J = -D∇u, ∂u/∂t = ∇·(D∇u) = -∇·J (correct)

References:
- Stupp et al., NEJM 2005 (standard protocol)
- Rockne et al., Cancer Research 2010 (LQ model in PDE)
- Swanson et al., Cell 2003 (FK model)
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
GRID_3D = (128, 128, 64)
DOMAIN_MM = (100.0, 100.0, 50.0)
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Biophysical parameters
RHO_MIN, RHO_MAX = 0.005, 0.12      # day^-1
K_CARRYING = 1.0
R_ANISO = 0.85
D0_BASE = 0.05                       # mm^2/day

# Therapy parameters (Stupp protocol)
SURGERY_DAY = 15
RAD_START_DAY = 20
RAD_END_DAY = 50
RAD_DOSE_GY = 2.0                    # Gy per fraction
RAD_DAYS_PER_WEEK = 5
ALPHA_RAD = 0.05                     # Gy^-1
BETA_RAD = 0.03                      # Gy^-2

CHEMO_START_DAY = 20
CHEMO_END_DAY = 80
GAMMA_CHEMO = 0.015                  # clearance rate
CHEMO_PK_HALFLIFE = 1.5              # days
CHEMO_DAILY_DOSE = 1.0               # normalized units

T_TOTAL_DAYS = 90.0


# --------------------------------------------------------------------------- #
# Phase 3 Field Generators
# --------------------------------------------------------------------------- #
def generate_dti_tensor_field() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx),
        np.linspace(-1, 1, ny),
        np.linspace(-1, 1, nz),
        indexing="ij"
    )

    cc_mask = (np.abs(Y) < 0.15) & (np.abs(Z) < 0.1)
    cst_mask = (X**2 + Y**2 > 0.15) & (X**2 + Y**2 < 0.35) & (Z > -0.6)
    slf_mask = (np.abs(X) < 0.1) & (Z > -0.3)

    e1 = np.zeros((nx, ny, nz, 3))
    FA = np.zeros((nx, ny, nz))

    e1[cc_mask] = [1, 0, 0];  FA[cc_mask] = 0.75
    e1[cst_mask] = [0, 0, 1]; FA[cst_mask] = 0.70
    e1[slf_mask] = [0, 1, 0]; FA[slf_mask] = 0.65

    iso_mask = ~(cc_mask | cst_mask | slf_mask)
    FA[iso_mask] = 0.15
    np.random.seed(42)
    rand = np.random.randn(nx, ny, nz, 3)
    rand_norm = np.linalg.norm(rand, axis=-1, keepdims=True) + 1e-8
    e1[iso_mask] = rand[iso_mask] / rand_norm[iso_mask]

    D_base = D0_BASE * np.ones(GRID_3D)
    lam1 = D_base * (1 + 3 * FA)
    lam2 = D_base * (1 - FA)
    lam3 = lam2

    tensor = np.zeros((nx, ny, nz, 3, 3))
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                v = e1[i, j, k]
                if np.linalg.norm(v) < 1e-8:
                    v = np.array([1.0, 0.0, 0.0])
                if abs(v[0]) < 0.9:
                    w = np.array([1.0, 0.0, 0.0])
                else:
                    w = np.array([0.0, 1.0, 0.0])
                e2 = w - np.dot(w, v) * v
                e2 = e2 / (np.linalg.norm(e2) + 1e-8)
                e3 = np.cross(v, e2)
                R = np.column_stack([v, e2, e3])
                Dd = np.diag([lam1[i,j,k], lam2[i,j,k], lam3[i,j,k]])
                tensor[i, j, k] = R @ Dd @ R.T

    tensor = (1 - R_ANISO) * D0_BASE * np.eye(3) + R_ANISO * tensor
    return tensor, FA, e1


def generate_rho_field() -> np.ndarray:
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx),
        np.linspace(-1, 1, ny),
        np.linspace(-1, 1, nz),
        indexing="ij"
    )
    rho = RHO_MIN + (RHO_MAX - RHO_MIN) * (
        np.exp(-((X - 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12) +
        np.exp(-((X + 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12)
    )
    rho = gaussian_filter(rho, sigma=(1.5, 1.5, 1.0))
    rho = 0.005 + 0.115 * (rho - rho.min()) / (rho.max() - rho.min() + 1e-12)
    return rho


# --------------------------------------------------------------------------- #
# 3D Differential Operators (Fick's Law: J = -D∇u, ∂u/∂t = ∇·(D∇u) = -∇·J)
# --------------------------------------------------------------------------- #
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


def div_3d(vx: np.ndarray, vy: np.ndarray, vz: np.ndarray, dx: float) -> np.ndarray:
    return grad_3d(vx, dx)[0] + grad_3d(vy, dx)[1] + grad_3d(vz, dx)[2]


def anisotropic_diffusion(u: np.ndarray, tensor: np.ndarray, dx: float) -> np.ndarray:
    """Returns ∇·(D∇u). Fick's law: J = -D∇u. Returns -∇·J = ∇·(D∇u)."""
    ux, uy, uz = grad_3d(u, dx)
    Jx = -(tensor[:,:,:,0,0]*ux + tensor[:,:,:,0,1]*uy + tensor[:,:,:,0,2]*uz)
    Jy = -(tensor[:,:,:,1,0]*ux + tensor[:,:,:,1,1]*uy + tensor[:,:,:,1,2]*uz)
    Jz = -(tensor[:,:,:,2,0]*ux + tensor[:,:,:,2,1]*uy + tensor[:,:,:,2,2]*uz)
    return -div_3d(Jx, Jy, Jz, dx)


# --------------------------------------------------------------------------- #
# Therapy Models
# --------------------------------------------------------------------------- #
def surgical_resection(
    u: np.ndarray,
    core_threshold: float = 0.05,
) -> np.ndarray:
    """Remove only the high-density tumor core (u >= core_threshold).
    Low-density infiltrating cells (u < core_threshold) are preserved to
    model microscopic residual disease and enable post-surgical recurrence.
    """
    u_resected = u.copy()
    core_mask = u >= core_threshold
    u_resected[core_mask] = 0.0
    return u_resected


def radiation_kill_fraction(alpha: float, beta: float, dose_gy: float) -> float:
    return 1.0 - np.exp(-alpha * dose_gy - beta * dose_gy**2)


def radiation_mask(
    u_pre_surg: np.ndarray,
    threshold: float = 0.02,
    margin_mm: float = 3.0,
) -> np.ndarray:
    mask = u_pre_surg >= threshold
    if not np.any(mask):
        return np.zeros_like(u_pre_surg, dtype=bool)
    from scipy.ndimage import distance_transform_edt
    dist = distance_transform_edt(~mask) * DX_MM
    return dist <= margin_mm


def chemo_concentration(day: float, start: float, end: float,
                        half_life: float = CHEMO_PK_HALFLIFE,
                        daily_dose: float = CHEMO_DAILY_DOSE) -> float:
    if day < start or day > end:
        return 0.0
    k_el = np.log(2) / half_life
    acc = 1.0 / (1.0 - np.exp(-k_el))
    return daily_dose * acc


# --------------------------------------------------------------------------- #
# CFL Time Step
# --------------------------------------------------------------------------- #
def compute_cfl_dt(tensor: np.ndarray, dx: float, safety: float = 0.2) -> Tuple[float, int]:
    max_diag = np.max(np.abs(tensor[:, :, :, [0,1,2], [0,1,2]]))
    max_eig = max_diag * 1.5
    dt = safety * dx**2 / (2 * 3 * max_eig)
    n_steps = max(1, int(np.ceil(T_TOTAL_DAYS / dt)))
    actual_dt = T_TOTAL_DAYS / n_steps
    print(f"[CFL] max_eig~{max_eig:.4f}, dx={dx:.3f} -> dt={actual_dt:.5f} days, "
          f"n_steps={n_steps} (T={T_TOTAL_DAYS}d)")
    return actual_dt, n_steps


# --------------------------------------------------------------------------- #
# Main Time-Stepping Loop (FIXED)
# --------------------------------------------------------------------------- #
def run_therapy_simulation(
    tensor: np.ndarray,
    rho: np.ndarray,
    dt: float,
    n_steps: int,
    dx: float,
) -> Tuple[np.ndarray, List[float], List[float], np.ndarray, np.ndarray]:
    
    u = np.zeros(GRID_3D)
    cx, cy, cz = GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2
    # Gaussian seed: peak ~0.5 at center, decaying to ~0.01 at ~4 voxels
    # This creates a core (u >= 0.05) of ~3mm radius with infiltrating tails
    for i in range(GRID_3D[0]):
        for j in range(GRID_3D[1]):
            for k in range(GRID_3D[2]):
                r2 = ((i - cx)**2 + (j - cy)**2 + (k - cz)**2) * dx**2
                u[i, j, k] = 0.5 * np.exp(-r2 / (2 * 4.0**2))  # sigma = 4mm

    volume_history = []
    time_points = []
    u_pre_surg = None
    u_post_surg = None

    current_day = 0.0

    rad_kill = radiation_kill_fraction(ALPHA_RAD, BETA_RAD, RAD_DOSE_GY)
    rad_mask_arr = None

    # Track radiation application: apply ONCE PER DAY, not every time-step
    last_rad_day = -1

    for step in range(n_steps):
        # --- Therapy events ---
        
        # Surgery (Day 15)
        if abs(current_day - SURGERY_DAY) < dt/2:
            u_pre_surg = u.copy()
            u = surgical_resection(u, core_threshold=0.05)
            u_post_surg = u.copy()
            rad_mask_arr = radiation_mask(u_pre_surg, threshold=0.02, margin_mm=3.0)
            print(f"[DAY {current_day:.1f}] Surgery: volume={np.sum(u)*dx**3:.1f} mm^3")

        # Radiation (Days 20-50, weekdays only) — ONCE PER DAY
        if RAD_START_DAY <= current_day <= RAD_END_DAY:
            day_int = int(np.floor(current_day))
            day_of_week = day_int % 7
            # Apply radiation once when we cross into a new weekday
            if day_of_week < 5 and rad_mask_arr is not None and day_int != last_rad_day:
                u[rad_mask_arr] *= (1.0 - rad_kill)
                last_rad_day = day_int

        # Chemotherapy (Days 20-80) — continuous sink
        if CHEMO_START_DAY <= current_day <= CHEMO_END_DAY:
            C_t = chemo_concentration(current_day, CHEMO_START_DAY, CHEMO_END_DAY)
            if C_t > 0:
                u -= GAMMA_CHEMO * C_t * u * dt

        # --- PDE Step (diff/react ALWAYS computed) ---
        diff = anisotropic_diffusion(u, tensor, dx)
        react = rho * u * (1.0 - u / K_CARRYING)
        u = u + dt * (diff + react)
        u = np.clip(u, 0.0, K_CARRYING)

        # Neumann BC
        u[0,:,:]=u[1,:,:]; u[-1,:,:]=u[-2,:,:]
        u[:,0,:]=u[:,1,:]; u[:,-1,:]=u[:,-2,:]
        u[:,:,0]=u[:,:,1]; u[:,:,-1]=u[:,:,-2]

        current_day += dt
        volume_history.append(np.sum(u) * dx**3)
        time_points.append(current_day)

        if step % max(1, n_steps//20) == 0:
            print(f"  Step {step}/{n_steps}: day={current_day:.1f}, "
                  f"vol={volume_history[-1]:.1f} mm^3, max_u={np.max(u):.3f}")

    return u, volume_history, time_points, u_pre_surg, u_post_surg


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def compute_therapy_metrics(
    volume_history: List[float],
    time_points: List[float],
    u_pre_surg: np.ndarray,
    u_final: np.ndarray,
    dx: float,
) -> Dict:
    pre_vol = np.sum(u_pre_surg) * dx**3 if u_pre_surg is not None else 0.0
    post_surg_vol = min(volume_history) if volume_history else 0.0
    reduction_pct = (pre_vol - post_surg_vol) / pre_vol * 100 if pre_vol > 0 else 0.0
    final_vol = volume_history[-1] if volume_history else 0.0

    days_to_recur = None
    for v, t in zip(volume_history, time_points):
        if v >= pre_vol * 0.95 and t > SURGERY_DAY:
            days_to_recur = t
            break

    if u_final is not None and np.any(u_final >= 0.02):
        mask = u_final >= 0.02
        coords = np.argwhere(mask)
        center = np.array([GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2])
        rel = coords - center
        dx_mm = np.max(np.abs(rel[:, 0])) * dx
        dy_mm = np.max(np.abs(rel[:, 1])) * dx
        asym = (dx_mm - dy_mm) / (dx_mm + dy_mm + 1e-8)
    else:
        asym = 0.0

    return {
        "pre_surgical_volume_mm3": float(pre_vol),
        "post_surgical_min_volume_mm3": float(post_surg_vol),
        "resection_reduction_pct": float(reduction_pct),
        "days_to_recurrence": float(days_to_recur) if days_to_recur else -1.0,
        "final_day90_volume_mm3": float(final_vol),
        "recurrence_asymmetry_index": float(asym),
        "simulation_days": T_TOTAL_DAYS,
        "grid": GRID_3D,
        "dx_mm": dx,
    }


# --------------------------------------------------------------------------- #
# Visualization (Dynamic vmax + yellow contour)
# --------------------------------------------------------------------------- #
def save_outputs(
    u_final: np.ndarray,
    u_pre_surg: np.ndarray,
    u_post_surg: np.ndarray,
    volume_history: List[float],
    time_points: List[float],
    metrics: Dict,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    nx, ny, nz = GRID_3D
    sx, sy, sz = nx//2, ny//2, nz//2

    vmax_pre = max(0.02, float(np.max(u_pre_surg)) if u_pre_surg is not None else 0.02)
    vmax_post = max(0.02, float(np.max(u_post_surg)) if u_post_surg is not None else 0.02)
    vmax_final = max(0.02, float(np.max(u_final)))

    fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)

    # Panel 1: Volume trajectory
    axes[0, 0].plot(time_points, volume_history, 'k-', linewidth=2)
    axes[0, 0].axvline(SURGERY_DAY, color='red', linestyle='--', linewidth=2, label='Surgery (Day 15)')
    axes[0, 0].axvspan(RAD_START_DAY, RAD_END_DAY, alpha=0.2, color='orange', label='Radiotherapy')
    axes[0, 0].axvspan(CHEMO_START_DAY, CHEMO_END_DAY, alpha=0.15, color='blue', label='TMZ Chemo')
    axes[0, 0].set_xlabel('Time (days)')
    axes[0, 0].set_ylabel('Tumor Volume (mm$^3$)')
    axes[0, 0].set_title('Tumor Volume Trajectory (Stupp Protocol)', fontweight='bold')
    axes[0, 0].legend(fontsize=9)
    axes[0, 0].grid(True, alpha=0.3)

    # Panel 2: Pre-surgical
    if u_pre_surg is not None:
        im2 = axes[0, 1].imshow(u_pre_surg[sx,:,:].T, cmap='hot', origin='lower',
                                 vmin=0, vmax=vmax_pre,
                                 extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        axes[0, 1].contour(u_pre_surg[sx,:,:].T, levels=[0.02], colors='yellow', linewidths=1.5)
        axes[0, 1].set_title('Pre-Surgical Lesion (Day 15, Sagittal)', fontweight='bold')
        axes[0, 1].set_xlabel('y [mm]'); axes[0, 1].set_ylabel('z [mm]')
        plt.colorbar(im2, ax=axes[0, 1], shrink=0.8, label='Density')

    # Panel 3: Post-resection + RT zone
    if u_post_surg is not None:
        im3 = axes[1, 0].imshow(u_post_surg[sx,:,:].T, cmap='hot', origin='lower',
                                 vmin=0, vmax=vmax_post,
                                 extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        axes[1, 0].contour(u_post_surg[sx,:,:].T, levels=[0.02], colors='yellow', linewidths=1.5)
        rad_mask_arr = radiation_mask(u_pre_surg, threshold=0.02, margin_mm=10.0)
        axes[1, 0].contour(rad_mask_arr[sx,:,:].T, levels=[0.5], colors='cyan', linewidths=2, linestyles='--')
        axes[1, 0].set_title('Post-Resection Cavity + RT Target (Sagittal)', fontweight='bold')
        axes[1, 0].set_xlabel('y [mm]'); axes[1, 0].set_ylabel('z [mm]')
        plt.colorbar(im3, ax=axes[1, 0], shrink=0.8, label='Density')

    # Panel 4: Day 90 recurrence
    im4 = axes[1, 1].imshow(u_final[sx,:,:].T, cmap='hot', origin='lower',
                             vmin=0, vmax=vmax_final,
                             extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
    axes[1, 1].contour(u_final[sx,:,:].T, levels=[metrics.get('threshold_used', 0.02)], 
                       colors='yellow', linewidths=1.5)
    axes[1, 1].set_title('Recurrent Lesion at Day 90 (Sagittal)', fontweight='bold')
    axes[1, 1].set_xlabel('y [mm]'); axes[1, 1].set_ylabel('z [mm]')
    plt.colorbar(im4, ax=axes[1, 1], shrink=0.8, label='Density')

    fig.suptitle('Phase 4: Virtual Stupp Protocol Therapy Simulation (T=90 days)',
                 fontsize=14, fontweight='bold', y=1.01)
    fig.savefig(out_dir / 'phase4_virtual_therapy.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

    with open(out_dir / 'phase4_therapy_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    print(f"[SAVE] Figure -> {out_dir / 'phase4_virtual_therapy.png'}")
    print(f"[SAVE] Metrics -> {out_dir / 'phase4_therapy_metrics.json'}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    os.makedirs("output", exist_ok=True)

    print("[PHASE 4] Generating 3D DTI tensor field...")
    tensor, FA, e1 = generate_dti_tensor_field()
    print(f"  Tensor shape: {tensor.shape}")
    print(f"  FA range: [{FA.min():.3f}, {FA.max():.3f}]")

    print("[PHASE 4] Generating proliferation field rho...")
    rho = generate_rho_field()
    print(f"  rho range: [{rho.min():.5f}, {rho.max():.5f}] day^-1")

    dt, n_steps = compute_cfl_dt(tensor, DX_MM)

    print(f"[PHASE 4] Running 90-day therapy simulation: {n_steps} steps, dt={dt:.5f} days...")
    u_final, vol_hist, t_pts, u_pre, u_post = run_therapy_simulation(
        tensor, rho, dt, n_steps, DX_MM
    )

    print("[PHASE 4] Computing therapy metrics...")
    metrics = compute_therapy_metrics(vol_hist, t_pts, u_pre, u_final, DX_MM)

    print(f"[METRICS] Pre-surg volume: {metrics['pre_surgical_volume_mm3']:.1f} mm^3")
    print(f"[METRICS] Post-surg min: {metrics['post_surgical_min_volume_mm3']:.1f} mm^3")
    print(f"[METRICS] Resection reduction: {metrics['resection_reduction_pct']:.1f}%")
    print(f"[METRICS] Days to recurrence: {metrics['days_to_recurrence']:.1f}")
    print(f"[METRICS] Day 90 volume: {metrics['final_day90_volume_mm3']:.1f} mm^3")
    print(f"[METRICS] Recurrence asymmetry: {metrics['recurrence_asymmetry_index']:.3f}")

    save_outputs(u_final, u_pre, u_post, vol_hist, t_pts, metrics, Path("output"))

    print("[PHASE 4] Complete.")


if __name__ == "__main__":
    main()