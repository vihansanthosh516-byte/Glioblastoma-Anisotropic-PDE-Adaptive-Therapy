#!/usr/bin/env python3
"""
Phase 3: 3D Anisotropic DTI Tensor Integration for GBM (STABLE)

Full 3D reaction-diffusion with patient-specific DTI-derived diffusion tensors.
Fixed CFL stability, density clamping, and directional metrics.

PDE: ∂u/∂t = ∇·(D_cell∇u) + ρ u (1 - u/K)  with K=1.0, u∈[0,1]

References:
- Hormuth et al., Cancer Research 2017 (anisotropic glioma PDE)
- Swanson et al., Cell 2003 (FK model)
- Tournier et al., NeuroImage 2007 (CSD)
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
# Configuration
# --------------------------------------------------------------------------- #
GRID_3D = (128, 128, 64)        # (nx, ny, nz)
DOMAIN_MM = (100.0, 100.0, 50.0)  # mm in x, y, z
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Biophysical parameters
RHO_MIN, RHO_MAX = 0.005, 0.12   # day^-1
K_CARRYING = 1.0                 # normalized carrying capacity (density cap)
R_ANISO = 0.85                   # DTI guidance strength [0,1]
D0_BASE = 0.05                   # mm^2/day baseline

# Simulation time
T_TOTAL_DAYS = 30.0


# --------------------------------------------------------------------------- #
# CFL-Safe Time Step Calculation
# --------------------------------------------------------------------------- #
def compute_cfl_dt(tensor: np.ndarray, dx: float, safety: float = 0.25) -> Tuple[float, int]:
    """
    Compute maximum stable time step for 3D explicit anisotropic diffusion.
    CFL: dt <= safety * dx^2 / (2 * ndim * max_eigenvalue)
    """
    # Maximum eigenvalue across all voxels
    # For speed, sample diagonal max as upper bound
    max_diag = np.max(np.abs(tensor[:, :, :, [0,1,2], [0,1,2]]))
    max_eig = max_diag * 1.5  # conservative factor for off-diagonals
    
    dt = safety * dx**2 / (2 * 3 * max_eig)
    n_steps = max(1, int(np.ceil(T_TOTAL_DAYS / dt)))
    actual_dt = T_TOTAL_DAYS / n_steps
    
    print(f"[CFL] max_eig~{max_eig:.4f}, dx={dx:.3f} → dt={actual_dt:.5f} days, "
          f"n_steps={n_steps} (target T={T_TOTAL_DAYS}d)")
    return actual_dt, n_steps


# --------------------------------------------------------------------------- #
# Phase 1/2 Field Generators
# --------------------------------------------------------------------------- #
def load_phase1_fields() -> Tuple[np.ndarray, np.ndarray]:
    """Generate D_base(x,y,z) and rho(x,y,z) matching Phase 1 stats."""
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx),
        np.linspace(-1, 1, ny),
        np.linspace(-1, 1, nz),
        indexing="ij"
    )

    # D: high at invasive fronts (periphery)
    D = 0.01 + 0.49 * (
        np.exp(-((X - 0.4)**2 + (Y - 0.4)**2 + Z**2) / 0.15) +
        np.exp(-((X + 0.4)**2 + (Y - 0.4)**2 + Z**2) / 0.15)
    )
    D = gaussian_filter(D, sigma=(1.5, 1.5, 1.0))
    D = 0.01 + 0.49 * (D - D.min()) / (D.max() - D.min() + 1e-12)

    # rho: high in core (proliferative)
    rho = 0.005 + 0.115 * (
        np.exp(-((X - 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12) +
        np.exp(-((X + 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12)
    )
    rho = gaussian_filter(rho, sigma=(1.5, 1.5, 1.0))
    rho = 0.005 + 0.115 * (rho - rho.min()) / (rho.max() - rho.min() + 1e-12)

    return D, rho


# --------------------------------------------------------------------------- #
# 3D DTI Tensor Construction
# --------------------------------------------------------------------------- #
def construct_dti_tensor_field(
    D_base: np.ndarray,
    r_aniso: float = R_ANISO,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Construct D_cell = D0 * ((1-r)I + r * D_water_norm) from synthetic tracts.
    Returns: tensor (nx,ny,nz,3,3), FA, e1
    """
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx),
        np.linspace(-1, 1, ny),
        np.linspace(-1, 1, nz),
        indexing="ij"
    )

    # Tract masks
    cc_mask = (np.abs(Y) < 0.15) & (np.abs(Z) < 0.1)          # corpus callosum (x-dir)
    cst_mask = (X**2 + Y**2 > 0.15) & (X**2 + Y**2 < 0.35) & (Z > -0.6)  # CST (z-dir)
    slf_mask = (np.abs(X) < 0.1) & (Z > -0.3)                  # SLF (y-dir)

    e1 = np.zeros((nx, ny, nz, 3))
    FA = np.zeros((nx, ny, nz))

    # Assign principal directions
    e1[cc_mask] = [1, 0, 0]
    FA[cc_mask] = 0.75
    e1[cst_mask] = [0, 0, 1]
    FA[cst_mask] = 0.70
    e1[slf_mask] = [0, 1, 0]
    FA[slf_mask] = 0.65

    # Isotropic remainder
    iso_mask = ~(cc_mask | cst_mask | slf_mask)
    FA[iso_mask] = 0.15
    np.random.seed(42)
    rand = np.random.randn(nx, ny, nz, 3)
    rand_norm = np.linalg.norm(rand, axis=-1, keepdims=True) + 1e-8
    e1[iso_mask] = rand[iso_mask] / rand_norm[iso_mask]

    # Eigenvalues
    lam1 = D_base * (1 + 3 * FA)
    lam2 = D_base * (1 - FA)
    lam3 = lam2

    # Build tensor: R @ diag(lam1,lam2,lam3) @ R.T
    tensor = np.zeros((nx, ny, nz, 3, 3))
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                v = e1[i, j, k]
                if np.linalg.norm(v) < 1e-8:
                    v = np.array([1.0, 0.0, 0.0])
                # Gram-Schmidt
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

    # Blend with isotropic: D_cell = (1-r)*D0*I + r*tensor
    I3 = np.eye(3)
    D0 = D0_BASE
    tensor = (1 - r_aniso) * D0 * I3 + r_aniso * tensor

    return tensor, FA, e1


# --------------------------------------------------------------------------- #
# 3D Gradient / Divergence / Flux
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


def anisotropic_flux(u: np.ndarray, tensor: np.ndarray, dx: float
                    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    ux, uy, uz = grad_3d(u, dx)
    Jx = -(tensor[:,:,:,0,0]*ux + tensor[:,:,:,0,1]*uy + tensor[:,:,:,0,2]*uz)
    Jy = -(tensor[:,:,:,1,0]*ux + tensor[:,:,:,1,1]*uy + tensor[:,:,:,1,2]*uz)
    Jz = -(tensor[:,:,:,2,0]*ux + tensor[:,:,:,2,1]*uy + tensor[:,:,:,2,2]*uz)
    return Jx, Jy, Jz


# --------------------------------------------------------------------------- #
# Single Time Step (Explicit Euler with Clamping)
# --------------------------------------------------------------------------- #
def step_anisotropic_3d(
    u: np.ndarray,
    tensor: np.ndarray,
    rho: np.ndarray,
    dx: float,
    dt: float,
    K: float = K_CARRYING,
) -> np.ndarray:
    # 1. Anisotropic diffusion
    Jx, Jy, Jz = anisotropic_flux(u, tensor, dx)
    diff = -div_3d(Jx, Jy, Jz, dx)

    # 2. Logistic reaction
    react = rho * u * (1.0 - u / K)

    # 3. Explicit Euler
    u_new = u + dt * (diff + react)

    # 4. Neumann BC (zero flux)
    u_new[0,:,:] = u_new[1,:,:]
    u_new[-1,:,:] = u_new[-2,:,:]
    u_new[:,0,:] = u_new[:,1,:]
    u_new[:,-1,:] = u_new[:,-2,:]
    u_new[:,:,0] = u_new[:,:,1]
    u_new[:,:,-1] = u_new[:,:,-2]

    # 5. HARD CLAMP to [0, K]
    u_new = np.clip(u_new, 0.0, K)

    return u_new


# --------------------------------------------------------------------------- #
# Isotropic Baseline (for comparison)
# --------------------------------------------------------------------------- #
def run_isotropic_baseline(
    D_iso: np.ndarray,
    rho: np.ndarray,
    dt: float,
    n_steps: int,
    dx: float,
    K: float = K_CARRYING,
) -> np.ndarray:
    u = np.zeros(GRID_3D)
    cx, cy, cz = GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2
    u[cx-2:cx+2, cy-2:cy+2, cz-2:cz+2] = 0.1

    for _ in range(n_steps):
        ux, uy, uz = grad_3d(u, dx)
        Jx, Jy, Jz = -D_iso * ux, -D_iso * uy, -D_iso * uz
        diff = -div_3d(Jx, Jy, Jz, dx)
        react = rho * u * (1.0 - u / K)
        u = u + dt * (diff + react)
        u = np.clip(u, 0.0, K)
        # BC
        u[0,:,:]=u[1,:,:]; u[-1,:,:]=u[-2,:,:]
        u[:,0,:]=u[:,1,:]; u[:,-1,:]=u[:,-2,:]
        u[:,:,0]=u[:,:,1]; u[:,:,-1]=u[:,:,-2]

    return u


# --------------------------------------------------------------------------- #
# Metrics (Directional Invasion)
# --------------------------------------------------------------------------- #
def compute_anisotropy_metrics(
    u: np.ndarray,
    tensor: np.ndarray,
    dx: float,
) -> Dict:
    """Compute directional invasion metrics with dynamic density threshold."""
    nx, ny, nz = u.shape
    center = np.array([nx//2, ny//2, nz//2])

    # Global peak density
    u_max = float(np.max(u))

    # Dynamic threshold: 20% of peak density, minimum 0.001
    threshold = max(0.20 * u_max, 0.001)
    mask = u >= threshold
    if not np.any(mask):
        return {
            "mean_tensor_anisotropy_ratio": 1.0,
            "max_invasion_x_mm": 0.0,
            "max_invasion_y_mm": 0.0,
            "max_invasion_z_mm": 0.0,
            "asymmetry_index": 0.0,
            "preferred_migration_angle_deg": 0.0,
            "simulation_days": T_TOTAL_DAYS,
            "grid": GRID_3D,
            "dx_mm": dx,
            "threshold_used": threshold,
            "peak_density": u_max,
        }

    # Anisotropy ratio in tumor region (where u >= threshold)
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[:, :, :, ::-1]  # descending
    aniso_ratio = evals[:, :, :, 0] / (evals[:, :, :, 2] + 1e-12)
    mean_aniso = float(np.mean(aniso_ratio[mask]))

    # Max invasion distance along each axis from center
    coords = np.argwhere(mask)
    rel = coords - center
    dist_x = np.max(np.abs(rel[:, 0])) * dx
    dist_y = np.max(np.abs(rel[:, 1])) * dx
    dist_z = np.max(np.abs(rel[:, 2])) * dx

    # Asymmetry: X vs Y (corpus callosum is X-direction)
    asym = (dist_x - dist_y) / (dist_x + dist_y + 1e-8)

    # Preferred angle from principal eigenvector at boundary
    # Boundary: mask voxels with at least one neighbor below threshold
    boundary = mask & (
        (np.roll(u, 1, axis=0) < threshold) | (np.roll(u, -1, axis=0) < threshold) |
        (np.roll(u, 1, axis=1) < threshold) | (np.roll(u, -1, axis=1) < threshold) |
        (np.roll(u, 1, axis=2) < threshold) | (np.roll(u, -1, axis=2) < threshold)
    )
    # Remove wrap-around artifacts at domain edges
    boundary[0,:,:] = False; boundary[-1,:,:] = False
    boundary[:,0,:] = False; boundary[:,-1,:] = False
    boundary[:,:,0] = False; boundary[:,:,-1] = False

    if np.any(boundary):
        evecs = np.linalg.eigh(tensor[boundary])[1]  # (N,3,3) ascending
        e1_boundary = evecs[:, :, -1]  # last column = largest eigenvalue
        angles = np.arctan2(e1_boundary[:, 1], e1_boundary[:, 0]) * 180 / np.pi
        pref_angle = float(np.mean(angles))
    else:
        pref_angle = 0.0

    return {
        "mean_tensor_anisotropy_ratio": mean_aniso,
        "max_invasion_x_mm": float(dist_x),
        "max_invasion_y_mm": float(dist_y),
        "max_invasion_z_mm": float(dist_z),
        "asymmetry_index": float(asym),
        "preferred_migration_angle_deg": pref_angle,
        "simulation_days": T_TOTAL_DAYS,
        "grid": GRID_3D,
        "dx_mm": dx,
        "threshold_used": threshold,
        "peak_density": u_max,
    }


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def save_outputs_3d(
    u_aniso: np.ndarray,
    u_iso: np.ndarray,
    tensor: np.ndarray,
    metrics: Dict,
    out_dir: Path,
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    nx, ny, nz = GRID_3D
    sx, sy, sz = nx//2, ny//2, nz//2

    # FA map
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[:, :, :, ::-1]
    fa = np.sqrt(0.5 * ((evals[:,:,:,0]-evals[:,:,:,1])**2 +
                         (evals[:,:,:,1]-evals[:,:,:,2])**2 +
                         (evals[:,:,:,2]-evals[:,:,:,0])**2) /
                  (evals[:,:,:,0]**2 + evals[:,:,:,1]**2 + evals[:,:,:,2]**2 + 1e-12))

    # Dynamic vmax for better visualization
    vmax_aniso = max(0.02, float(np.max(u_aniso)))
    vmax_iso = max(0.02, float(np.max(u_iso)))

    # Invasion threshold from metrics
    threshold = metrics.get("threshold_used", 0.10)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    # Row 1: Anisotropic tumor slices
    aniso_slices = [u_aniso[sx,:,:], u_aniso[:,sy,:], u_aniso[:,:,sz]]
    aniso_titles = ["Sagittal (x=mid)", "Coronal (y=mid)", "Axial (z=mid)"]
    for ax, slc, title in zip(axes[0], aniso_slices, aniso_titles):
        im = ax.imshow(slc.T, cmap="hot", origin="lower", vmin=0, vmax=vmax_aniso,
                       extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        # Contour at invasion threshold
        ax.contour(slc.T, levels=[threshold], colors="yellow", linewidths=2,
                   extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.set_title(f"Anisotropic u: {title}")
        ax.set_xlabel("mm"); ax.set_ylabel("mm")
        plt.colorbar(im, ax=ax, shrink=0.8, label="Density")

    # Row 2: Isotropic baseline
    iso_slices = [u_iso[sx,:,:], u_iso[:,sy,:], u_iso[:,:,sz]]
    iso_titles = ["Sagittal", "Coronal", "Axial"]
    for ax, slc, title in zip(axes[1], iso_slices, iso_titles):
        im = ax.imshow(slc.T, cmap="hot", origin="lower", vmin=0, vmax=vmax_iso,
                       extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.contour(slc.T, levels=[threshold], colors="yellow", linewidths=2,
                   extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.set_title(f"Isotropic u: {title}")
        ax.set_xlabel("mm"); ax.set_ylabel("mm")
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle("Phase 3: 3D Anisotropic vs Isotropic Tumor Growth (T=30 days)",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.savefig(out_dir / "phase3_3d_anisotropic_dti.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    with open(out_dir / "phase3_3d_dti_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)

    print(f"[SAVE] Figure -> {out_dir / 'phase3_3d_anisotropic_dti.png'}")
    print(f"[SAVE] Metrics -> {out_dir / 'phase3_3d_dti_metrics.json'}")


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    os.makedirs("output", exist_ok=True)

    print("[PHASE 3] Generating 3D DTI tensor field...")
    D_base, rho = load_phase1_fields()
    tensor, FA, e1 = construct_dti_tensor_field(D_base, R_ANISO)
    print(f"  Tensor shape: {tensor.shape}")
    print(f"  FA range: [{FA.min():.3f}, {FA.max():.3f}]")
    print(f"  D_base range: [{D_base.min():.4f}, {D_base.max():.4f}] mm²/day")

    # CFL time step
    dt, n_steps = compute_cfl_dt(tensor, DX_MM)

    # Isotropic baseline (trace/3)
    D_iso = np.trace(tensor, axis1=3, axis2=4) / 3.0
    print("[PHASE 3] Running isotropic baseline...")
    u_iso = run_isotropic_baseline(D_iso, rho, dt, n_steps, DX_MM)

    # Anisotropic simulation
    u = np.zeros(GRID_3D)
    cx, cy, cz = GRID_3D[0]//2, GRID_3D[1]//2, GRID_3D[2]//2
    u[cx-2:cx+2, cy-2:cy+2, cz-2:cz+2] = 0.1

    print(f"[PHASE 3] Running 3D FULL TENSOR simulation: {n_steps} steps, dt={dt:.5f} days...")
    for step in range(n_steps):
        u = step_anisotropic_3d(u, tensor, rho, DX_MM, dt)
        if step % max(1, n_steps//10) == 0:
            vol = np.sum(u) * DX_MM**3
            max_u = np.max(u)
            print(f"  Step {step}/{n_steps}: volume={vol:.2f} mm³, max_density={max_u:.3f}")

    print("[PHASE 3] Computing anisotropy metrics...")
    metrics = compute_anisotropy_metrics(u, tensor, DX_MM)

    print(f"[METRICS] Mean aniso ratio: {metrics['mean_tensor_anisotropy_ratio']:.2f}")
    print(f"[METRICS] Invasion X/Y/Z: {metrics['max_invasion_x_mm']:.2f} / "
          f"{metrics['max_invasion_y_mm']:.2f} / {metrics['max_invasion_z_mm']:.2f} mm")
    print(f"[METRICS] Asymmetry index (X-Y): {metrics['asymmetry_index']:.3f}")
    print(f"[METRICS] Preferred angle: {metrics['preferred_migration_angle_deg']:.1f} deg")

    save_outputs_3d(u, u_iso, tensor, metrics, Path("output"))

    print("[PHASE 3] Complete.")


if __name__ == "__main__":
    main()