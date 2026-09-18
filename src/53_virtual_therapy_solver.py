#!/usr/bin/env python3
"""
Phase 4: Virtual Therapy Solver

Ingests Phase 3 anisotropic concentration field and applies a 60-day virtual
therapy protocol to predict treatment response and recurrence dynamics.

PDE: dc/dt = nabla.(D nabla c) + rho*c*(1-c) - K_drug(t,x)*c

Therapy protocol: 60-day timeline with 2 cycles of cytotoxic chemotherapy
- Cycle 1: Days 1-5 (K_max), Days 6-21 (exponential washout, half-life=3d)
- Cycle 2: Days 22-26 (K_max), Days 27-43 (exponential washout, half-life=3d)
- Days 44-60: Post-therapy surveillance

Drug penetration modulated by anisotropic diffusion tensor from Phase 3.
Spatially heterogeneous K_drug(x) = K_max * (c/c_max)^alpha * exp(-dist/lambda)

Outputs:
- output/phase4_virtual_therapy.png (4-panel figure)
- output/phase4_therapy_metrics.json
- output/phase4_treated_concentration.npy
- output/phase4_untreated_concentration.npy
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter, zoom as _zoom

warnings.filterwarnings("ignore", category=FutureWarning)


# =============================================================================
# Configuration
# =============================================================================
LOCAL_TEST = False  # set True for fast local dry-run

if LOCAL_TEST:
    GRID_3D = (32, 32, 16)
    T_THERAPY_DAYS = 10.0
else:
    GRID_3D = (128, 128, 64)
    T_THERAPY_DAYS = 60.0

DOMAIN_MM = (100.0, 100.0, 50.0)
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Biophysical parameters (from Phase 1/3)
K_CARRYING = 1.0
CFL_SAFETY = 0.25

# Therapy parameters
K_MAX = 0.8          # max daily kill rate (day^-1)
HALF_LIFE_DAYS = 3.0 # drug washout half-life (days)
CYCLE_ON = 5         # days of active dosing per cycle
CYCLE_OFF = 16       # days between cycles (washout)
N_CYCLES = 2         # number of cycles
ALPHA_PENETRATION = 0.5  # drug penetration exponent (c/c_max)^alpha
LAMBDA_PENETRATION = 5.0 # mm, spatial decay length

# Tumor seed
TUMOR_CENTER_VOX = (GRID_3D[0] // 2, GRID_3D[1] // 2, GRID_3D[2] // 2)


# =============================================================================
# 3D Differential Operators (same pattern as Phase 2/3)
# =============================================================================
def grad_3d(f: np.ndarray, dx: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
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
    return (grad_3d(vx, dx)[0] + grad_3d(vy, dx)[1] + grad_3d(vz, dx)[2])


def _align_to_grid(arr: np.ndarray) -> np.ndarray:
    if arr.shape == GRID_3D:
        return arr
    from scipy.ndimage import zoom as _zoom
    return _zoom(arr, tuple(GRID_3D[i] / arr.shape[i] for i in range(3)), order=1)


# =============================================================================
# Phase 3 Ingestion
# =============================================================================
def _load_phase3_aniso(out_dir: Path) -> np.ndarray:
    """Load Phase 3 anisotropic concentration field."""
    path = out_dir / "phase3_aniso_concentration.npy"
    if path.exists():
        print("[PHASE 4] Loading Phase 3 anisotropic concentration...")
        c_aniso = np.load(path)
        return _align_to_grid(c_aniso)
    print("[PHASE 4] Phase 3 aniso file not found; generating synthetic fallback...")
    return _synthetic_aniso_field()


def _synthetic_aniso_field() -> np.ndarray:
    """Fallback synthetic anisotropic field matching Phase 3 stats."""
    nx, ny, nz = GRID_3D
    cx, cy, cz = TUMOR_CENTER_VOX
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), np.linspace(-1, 1, nz),
        indexing="ij",
    )
    # Elongated along x-axis (corpus callosum direction)
    r2 = ((X - 0.3)**2 / 0.04 + (Y + 0.2)**2 / 0.09 + Z**2 / 0.16)
    c = np.exp(-r2)
    c = gaussian_filter(c, sigma=(2, 1, 1))
    return c / (c.max() + 1e-12)


def _load_phase3_tensor(out_dir: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load or reconstruct Phase 3 DTI tensor, FA, and e1 fields."""
    tensor_path = out_dir / "phase3_tensor.npy"
    if tensor_path.exists():
        print("[PHASE 4] Loading Phase 3 tensor field...")
        tensor = np.load(tensor_path)
        tensor = _align_to_grid(tensor)
        # Reconstruct FA and e1 from tensor
        FA, e1 = _decompose_tensor(tensor)
        return tensor, FA, e1
    print("[PHASE 4] Phase 3 tensor not found; reconstructing from Phase 3 outputs...")
    return _reconstruct_tensor_from_phase3_outputs()


def _decompose_tensor(tensor: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Extract FA and principal eigenvector from symmetric 3x3 tensor."""
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[..., ::-1]  # descending
    lam1, lam2, lam3 = evals[..., 0], evals[..., 1], evals[..., 2]
    fa = np.sqrt(0.5 * ((lam1 - lam2)**2 + (lam2 - lam3)**2 + (lam3 - lam1)**2))
    fa /= np.sqrt(lam1**2 + lam2**2 + lam3**2 + 1e-12)
    # Get principal eigenvector
    evecs = np.linalg.eigh(tensor)[1]  # ascending order
    e1 = evecs[..., -1]  # last column = largest eigenvalue
    return fa, e1


def _reconstruct_tensor_from_phase3_outputs() -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reconstruct tensor from Phase 3 D_base and FA/e1 logic (matches Phase 3)."""
    # Fallback: recreate the same synthetic tracts as Phase 3
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), np.linspace(-1, 1, nz),
        indexing="ij",
    )
    cc_mask = (np.abs(Y) < 0.15) & (np.abs(Z) < 0.1)
    cst_mask = (X**2 + Y**2 > 0.15) & (X**2 + Y**2 < 0.35) & (Z > -0.6)
    slf_mask = (np.abs(X) < 0.1) & (Z > -0.3)

    e1 = np.zeros((nx, ny, nz, 3))
    FA = np.full((nx, ny, nz), 0.15, dtype=np.float64)

    e1[cc_mask] = [1, 0, 0]; FA[cc_mask] = 0.75
    e1[cst_mask] = [0, 0, 1]; FA[cst_mask] = 0.70
    e1[slf_mask] = [0, 1, 0]; FA[slf_mask] = 0.65

    iso_mask = ~(cc_mask | cst_mask | slf_mask)
    np.random.seed(42)
    rand = np.random.randn(nx, ny, nz, 3)
    rand_norm = np.linalg.norm(rand, axis=-1, keepdims=True) + 1e-8
    e1[iso_mask] = rand[iso_mask] / rand_norm[iso_mask]

    D0_eff = np.full(GRID_3D, 0.05, dtype=np.float64)  # baseline
    d0 = (1.0 - FA) * D0_eff
    d_d = 3.0 * FA * D0_eff

    ex = e1[..., 0]; ey = e1[..., 1]; ez = e1[..., 2]
    tensor = np.zeros((nx, ny, nz, 3, 3), dtype=np.float64)
    tensor[..., 0, 0] = d0 + d_d * ex * ex
    tensor[..., 1, 1] = d0 + d_d * ey * ey
    tensor[..., 2, 2] = d0 + d_d * ez * ez
    tensor[..., 0, 1] = d_d * ex * ey
    tensor[..., 1, 0] = d_d * ex * ey
    tensor[..., 0, 2] = d_d * ex * ez
    tensor[..., 2, 0] = d_d * ex * ez
    tensor[..., 1, 2] = d_d * ey * ez
    tensor[..., 2, 1] = d_d * ey * ez

    return tensor, FA, e1


def _align_to_grid(arr: np.ndarray) -> np.ndarray:
    if arr.shape == GRID_3D:
        return arr
    from scipy.ndimage import zoom as _zoom
    return _zoom(arr, tuple(GRID_3D[i] / arr.shape[i] for i in range(3)), order=1)


# =============================================================================
# Phase 4: Therapy PDE Solver
# =============================================================================
def _drug_kill_rate(t: float, c: np.ndarray, c_max: float, tensor: np.ndarray,
                    e1: np.ndarray) -> np.ndarray:
    """
    Spatially heterogeneous drug kill rate K_drug(t, x).
    K(t) = K_max * drug_schedule(t) * (c/c_max)^alpha * exp(-distance_to_tumor/lambda)
    """
    # Temporal schedule: 2 cycles of 5 days on, 16 days off
    cycle_len = CYCLE_ON + CYCLE_OFF
    day_in_cycle = t % cycle_len
    if day_in_cycle < CYCLE_ON:
        k_t = 1.0
    else:
        # Exponential washout
        k_t = np.exp(-(day_in_cycle - CYCLE_ON) * np.log(2) / HALF_LIFE_DAYS)

    # Spatial modulation: drug penetrates better where tumor is denser
    spatial_mod = (c / (c_max + 1e-12)) ** ALPHA_PENETRATION

    # Distance-based decay from tumor center
    nx, ny, nz = GRID_3D
    cx, cy, cz = TUMOR_CENTER_VOX
    ix = np.arange(nx)[:, None, None]
    iy = np.arange(ny)[None, :, None]
    iz = np.arange(nz)[None, None, :]
    dist = np.sqrt((ix - cx)**2 + (iy - cy)**2 + (iz - cz)**2) * DX_MM
    dist_mod = np.exp(-dist / LAMBDA_PENETRATION)

    return K_MAX * k_t * spatial_mod * dist_mod


def compute_cfl_dt(tensor: np.ndarray, dx: float, safety: float = CFL_SAFETY) -> Tuple[float, int]:
    diag_max = float(np.max(np.abs(tensor[..., [0,1,2], [0,1,2]])))
    max_eig = max(diag_max * 1.5, 1e-9)
    dt_cap = safety * dx ** 2 / (2 * 3 * max_eig)
    n_steps = max(1, int(np.ceil(T_THERAPY_DAYS / dt_cap)))
    actual_dt = T_THERAPY_DAYS / n_steps
    print(f"[CFL] max_eig~{max_eig:.4f}, dx={dx:.3f} -> dt={actual_dt:.5f} days, "
          f"n_steps={n_steps} (target T={T_THERAPY_DAYS}d)")
    return actual_dt, n_steps


def _make_initial_seed() -> np.ndarray:
    c = np.zeros(GRID_3D, dtype=np.float64)
    cx, cy, cz = TUMOR_CENTER_VOX
    c[cx-2:cx+2, cy-2:cy+2, cz-2:cz+2] = 0.1
    return c


def step_anisotropic_3d(
    c: np.ndarray, tensor: np.ndarray, rho: np.ndarray,
    dx: float, dt: float, K_drug: np.ndarray,
    K: float = K_CARRYING,
) -> np.ndarray:
    """One explicit Euler step: dc/dt = nabla.(D nabla c) + rho*c*(1-c/K) - K_drug*c."""
    # Diffusion term
    cx, cy, cz = grad_3d(c, dx)
    D00 = tensor[..., 0, 0]; D01 = tensor[..., 0, 1]; D02 = tensor[..., 0, 2]
    D10 = tensor[..., 1, 0]; D11 = tensor[..., 1, 1]; D12 = tensor[..., 1, 2]
    D20 = tensor[..., 2, 0]; D21 = tensor[..., 2, 1]; D22 = tensor[..., 2, 2]

    Jx = -(D00 * cx + D01 * cy + D02 * cz)
    Jy = -(D10 * cx + D11 * cy + D12 * cz)
    Jz = -(D20 * cx + D21 * cy + D22 * cz)

    diff = -div_3d(Jx, Jy, Jz, dx)

    # Reaction + drug kill
    react = rho * c * (1.0 - c / K_CARRYING)
    kill = -K_drug * c

    c_new = c + dt * (diff + react + kill)

    # Neumann BC
    c_new[0, :, :] = c_new[1, :, :]; c_new[-1, :, :] = c_new[-2, :, :]
    c_new[:, 0, :] = c_new[:, 1, :]; c_new[:, -1, :] = c_new[:, -2, :]
    c_new[:, :, 0] = c_new[:, :, 1]; c_new[:, :, -1] = c_new[:, :, -1]
    return np.clip(c_new, 0.0, K_CARRYING)


def run_untreated_baseline(
    c_init: np.ndarray, tensor: np.ndarray, rho: np.ndarray,
    dt: float, n_steps: int, dx: float,
) -> np.ndarray:
    """Run simulation without drug for comparison."""
    c = c_init.copy()
    for step in range(n_steps):
        t = step * dt
        K_drug = np.zeros_like(c)  # no drug
        c = step_anisotropic_3d(c, tensor, rho, dx, dt, K_drug)
        if step % 20 == 0:
            vol = float(np.sum(c) * DX_MM ** 3)
            print(f"  [Untreated] step {step}/{n_steps}: volume={vol:.2f} mm^3")
    return c


def compute_therapy_metrics(
    c_treated: np.ndarray, c_untreated: np.ndarray,
    c_init: np.ndarray, tensor: np.ndarray, dx: float,
) -> Dict:
    """Compute therapy response metrics."""
    nx, ny, nz = GRID_3D
    center = np.array([nx // 2, ny // 2, nz // 2])

    # Threshold at 20% of peak
    u_max = float(np.max(c_treated))
    threshold = max(0.20 * u_max, 1e-3)
    mask_treated = c_treated >= threshold
    mask_untreated = c_untreated >= threshold

    # Volume
    vol_treated = float(np.sum(mask_treated) * dx ** 3)
    vol_untreated = float(np.sum(mask_untreated) * dx ** 3)
    vol_reduction_pct = 100.0 * (vol_untreated - vol_treated) / (vol_untreated + 1e-12)

    # Log-kill (total cell kill)
    total_init = float(np.sum(c_init))
    total_treated = float(np.sum(c_treated))
    log_kill = np.log10(total_init / (total_treated + 1e-12))

    # Time to recurrence (first day after cycle 2 where volume > 1.5x min post-therapy volume)
    # For simplicity, track volume at each step - this is computed during simulation
    # We'll return the final metrics; time-to-recurrence requires tracking during sim

    # Spatial selectivity: invasion directionality
    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[..., ::-1]
    aniso_ratio = evals[..., 0] / (evals[..., 2] + 1e-12)
    mean_aniso = float(np.mean(aniso_ratio[mask_treated])) if np.any(mask_treated) else 1.0

    # Max invasion distance per axis
    if np.any(mask_treated):
        coords = np.argwhere(mask_treated)
        rel = coords - center
        dist_x = float(np.max(np.abs(rel[:, 0])) * dx)
        dist_y = float(np.max(np.abs(rel[:, 1])) * dx)
        dist_z = float(np.max(np.abs(rel[:, 2])) * dx)
        asym = (dist_x - dist_y) / (dist_x + dist_y + 1e-8)
    else:
        dist_x = dist_y = dist_z = 0.0
        asym = 0.0

    return {
        "phase": 4,
        "volume_reduction_pct": vol_reduction_pct,
        "treated_volume_mm3": vol_treated,
        "untreated_volume_mm3": vol_untreated,
        "log_kill": float(log_kill),
        "mean_tensor_anisotropy_ratio": float(mean_aniso),
        "max_invasion_x_mm": dist_x,
        "max_invasion_y_mm": dist_y,
        "max_invasion_z_mm": dist_z,
        "asymmetry_index": asym,
        "simulation_days": T_THERAPY_DAYS,
        "grid": list(GRID_3D),
        "dx_mm": dx,
        "threshold_used": threshold,
        "peak_density": u_max,
    }


def save_therapy_visualization(
    c_treated: np.ndarray, c_untreated: np.ndarray,
    tensor: np.ndarray, metrics: Dict, out_path: Path,
) -> None:
    nx, ny, nz = GRID_3D
    sx, sy, sz = nx // 2, ny // 2, nz // 2

    evals = np.linalg.eigvalsh(tensor)
    evals = np.sort(evals, axis=-1)[..., ::-1]
    lam1, lam2, lam3 = evals[..., 0], evals[..., 1], evals[..., 2]
    fa_map = np.sqrt(0.5 * ((lam1 - lam2)**2 + (lam2 - lam3)**2 + (lam3 - lam1)**2))
    fa_map /= np.sqrt(lam1**2 + lam2**2 + lam3**2 + 1e-12)

    vmax_t = max(0.02, float(np.max(c_treated)))
    vmax_u = max(0.02, float(np.max(c_untreated)))
    threshold = metrics.get("threshold_used", 0.1)

    fig, axes = plt.subplots(2, 3, figsize=(15, 10), constrained_layout=True)

    # Row 1: Treated
    treated_slices = [c_treated[sx, :, :], c_treated[:, sy, :], c_treated[:, :, sz]]
    titles = ["Sagittal (x=mid)", "Coronal (y=mid)", "Axial (z=mid)"]
    for ax, slc, title in zip(axes[0], treated_slices, titles):
        im = ax.imshow(slc.T, cmap="hot", origin="lower", vmin=0, vmax=vmax_t,
                       extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.contour(slc.T, levels=[threshold], colors="lime", linewidths=2,
                   extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.set_title(f"Treated c: {title}")
        ax.set_xlabel("mm"); ax.set_ylabel("mm")
        plt.colorbar(im, ax=ax, shrink=0.8, label="Density")

    # Row 2: Untreated
    untreated_slices = [c_untreated[sx, :, :], c_untreated[:, sy, :], c_untreated[:, :, sz]]
    for ax, slc, title in zip(axes[1], untreated_slices, titles):
        im = ax.imshow(slc.T, cmap="hot", origin="lower", vmin=0, vmax=vmax_u,
                       extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.contour(slc.T, levels=[threshold], colors="lime", linewidths=2,
                   extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]])
        ax.set_title(f"Untreated c: {title}")
        ax.set_xlabel("mm"); ax.set_ylabel("mm")
        plt.colorbar(im, ax=ax, shrink=0.8)

    fig.suptitle(f"Phase 4: Virtual Therapy Response (T={T_THERAPY_DAYS} days)",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[SAVE] Figure -> {out_path}")


def save_therapy_metrics(metrics: Dict, out_path: Path) -> None:
    with open(out_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[SAVE] Metrics -> {out_path}")


# =============================================================================
# Main Phase 4 Entry
# =============================================================================
def run_phase4(
    out_dir: Path,
    rho: Optional[np.ndarray] = None,
    D_eff: Optional[np.ndarray] = None,
) -> Dict:
    """Run Phase 4 virtual therapy. Accepts Phase 1/3 fields in-memory or loads from disk."""
    print("=" * 70)
    print("  PHASE 4: VIRTUAL THERAPY SOLVER")
    print("=" * 70)

    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- Load Phase 1 rho ----
    if rho is None:
        rho_path = out_dir / "phase1_rho_posterior.npy"
        if rho_path.exists():
            rho_data = np.load(rho_path)
            rho = _align_to_grid(rho_data[0])
        else:
            print("[PHASE 4] Phase 1 rho not found; using synthetic fallback...")
            rho = _synthetic_rho_field()
    else:
        rho = _align_to_grid(rho)
    print(f"  rho range: [{rho.min():.5f}, {rho.max():.5f}] day^-1")

    # ---- Load Phase 3 aniso concentration (initial condition) ----
    c_init = _load_phase3_aniso(out_dir)
    c_init = np.clip(c_init, 0.0, K_CARRYING)
    print(f"  c_init range: [{c_init.min():.4f}, {c_init.max():.4f}]")

    # ---- Load Phase 3 tensor ----
    tensor, FA, e1 = _load_phase3_tensor(out_dir)
    print(f"  Tensor shape: {tensor.shape}, FA range: [{FA.min():.3f}, {FA.max():.3f}]")

    # ---- CFL time step ----
    dt, n_steps = compute_cfl_dt(tensor, DX_MM)
    print(f"  Therapy horizon: {T_THERAPY_DAYS} days, dt={dt:.5f} days, n_steps={n_steps}")

    # ---- Untreated baseline ----
    print("[PHASE 4] Running untreated baseline...")
    c_untreated = run_untreated_baseline(c_init.copy(), tensor, rho, dt, n_steps, DX_MM)

    # ---- Treated simulation ----
    print("[PHASE 4] Running treated simulation...")
    c_treated = c_init.copy()
    c_max = float(np.max(c_init))
    report_every = max(1, n_steps // 10)

    for step in range(n_steps):
        t = step * dt
        K_drug = _drug_kill_rate(t, c_treated, c_max, tensor, e1)
        c_treated = step_anisotropic_3d(c_treated, tensor, rho, DX_MM, dt, K_drug)

        if step % report_every == 0 or step == n_steps - 1:
            vol = float(np.sum(c_treated) * DX_MM ** 3)
            k_today = float(np.mean(K_drug))
            print(f"  step {step}/{n_steps} (t={t:.1f}d): volume={vol:.2f} mm^3, "
                  f"max_c={c_treated.max():.3f}, K_drug_mean={k_today:.4f}")

    # ---- Metrics + Outputs ----
    print("[PHASE 4] Computing therapy metrics...")
    metrics = compute_therapy_metrics(c_treated, c_init, c_init, tensor, DX_MM)
    print(f"  Volume reduction: {metrics['volume_reduction_pct']:.1f}%")
    print(f"  Log-kill: {metrics['log_kill']:.2f}")
    print(f"  Invasion X/Y/Z: {metrics['max_invasion_x_mm']:.1f} / "
          f"{metrics['max_invasion_y_mm']:.1f} / {metrics['max_invasion_z_mm']:.1f} mm")

    np.save(out_dir / "phase4_treated_concentration.npy", c_treated)
    np.save(out_dir / "phase4_untreated_concentration.npy", c_init)  # untreated final
    save_therapy_visualization(c_treated, c_init, tensor, metrics,
                               out_dir / "phase4_virtual_therapy.png")
    save_therapy_metrics(metrics, out_dir / "phase4_therapy_metrics.json")
    print("[PHASE 4] Complete.")
    return metrics


def _synthetic_rho_field() -> np.ndarray:
    """Fallback rho field if Phase 1 files missing."""
    nx, ny, nz = GRID_3D
    X, Y, Z = np.meshgrid(
        np.linspace(-1, 1, nx), np.linspace(-1, 1, ny), np.linspace(-1, 1, nz),
        indexing="ij",
    )
    rho = 0.005 + 0.115 * (
        np.exp(-((X - 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12) +
        np.exp(-((X + 0.2)**2 + (Y + 0.3)**2 + (Z + 0.1)**2) / 0.12)
    )
    rho = gaussian_filter(rho, sigma=(1.5, 1.5, 1.0))
    rho = RHO_MIN + (RHO_MAX - RHO_MIN) * (rho - rho.min()) / (rho.max() - rho.min() + 1e-12)
    return rho


def main() -> None:
    print("=" * 70)
    print("  PHASE 4: VIRTUAL THERAPY SOLVER (standalone)")
    print("=" * 70)
    run_phase4(Path("output"))


if __name__ == "__main__":
    main()