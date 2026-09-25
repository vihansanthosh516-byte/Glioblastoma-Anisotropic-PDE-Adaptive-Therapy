#!/usr/bin/env python3
"""Phase 3 Extension: Optimized 3D Volumetric Anisotropic Tumor Growth (8-Patient Cohort).

Optimizations for speed:
- Vectorized time stepping with pre-allocated arrays
- Spot-check positive-definite only (not per-voxel)
- Reduced I/O: save only final volumes, not full tensor fields
- Pre-compute face diffusivities once
- Fused divergence + step operations
- Single NPZ with compressed outputs
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path as _Path
from typing import Dict, Tuple

import numpy as np

warnings.filterwarnings("ignore")

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# 3D Grid and Physical Constants (OPTIMIZED DEFAULTS)
# --------------------------------------------------------------------------- #
GRID_SIZE = 50          # 50³ = 125K voxels (reduce to 40³ = 64K for faster runs)
DX = 1.0                # mm voxel spacing
DIM = 3

# Diffusion coefficients
D_WHITE = 0.013         # mm²/day (along tracts)
D_GRAY = 0.0013         # mm²/day (isotropic baseline)

from mu_glioma_loader import load_mu_glioma_params, real_cohort_stats
_params = load_mu_glioma_params()
_cohort = real_cohort_stats()

# Proliferation (real MU-Glioma median)
RHO = _cohort["rho_median"]
K = 1.0                 # normalized density

# TMZ PK parameters
TMZ_HALF_LIFE = 0.075
K_EL = np.log(2) / TMZ_HALF_LIFE
C_PEAK = 10.0
EC50 = 5.0
HILL_COEFF = 2.0
E_MAX = RHO * 1000.0

print(f"[Params] RHO = {RHO:.6f} /day, E_MAX = {E_MAX:.4f} "
      f"(cohort n={_cohort['n_total']}, "
      f"rho_median={_cohort['rho_median']:.6f})")

# Dosing: 5-on / 23-off
DOSE_DAYS_ON = 5
CYCLE_DAYS = 28

# Simulation (OPTIMIZED: fewer steps)
DT = 0.05               # days (CFL: dx²/(2*3*0.013) ≈ 12.8, 0.05 is safe)
SIM_DAYS = 180
N_STEPS = int(SIM_DAYS / DT)  # 3600 steps (was 4500)
SAVE_INTERVAL = 100     # save less frequently

# Initial tumor (legacy defaults; real seed uses V0 via initial_tumor_seed)
TUMOR_CENTER = (GRID_SIZE // 2 - 5, GRID_SIZE // 2 - 5, GRID_SIZE // 2)
TUMOR_RADIUS = 2

# --------------------------------------------------------------------------- #
# Cohort Definition (8 real MU-Glioma patients)
# --------------------------------------------------------------------------- #
# Select 8 real MU-Glioma patients: positive rho, R² > 0.7,
# spread across the rho range
_positive = {pid: p for pid, p in _params.items()
             if p["rho_per_day"] > 0 and p["r_squared"] > 0.7}
_sorted = sorted(_positive.items(), key=lambda x: x[1]["rho_per_day"])
if len(_sorted) < 8:
    raise RuntimeError(f"Only {len(_sorted)} eligible patients; need 8")
_indices = np.linspace(0, len(_sorted) - 1, 8).astype(int)
COHORT_PATIENTS = [_sorted[i][0] for i in _indices]

print(f"[Cohort] Selected {len(COHORT_PATIENTS)} real patients:")
for _pid in COHORT_PATIENTS:
    _p = _params[_pid]
    print(f"  {_pid}: rho={_p['rho_per_day']:.6f}, "
          f"V0={_p['V0_mm3']:.0f}, R²={_p['r_squared']:.3f}")

# Preserve the 4 orientation groups; map onto the 8 real patient IDs
_TRACT_ORIENTATIONS = [
    np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),  # Diagonal xy-plane
    np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
    np.array([0.0, 0.0, 1.0]),                  # Vertical (z-axis)
    np.array([0.0, 0.0, 1.0]),
    np.array([1.0, 0.0, 0.0]),                  # Horizontal (x-axis)
    np.array([1.0, 0.0, 0.0]),
    np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),  # Oblique
    np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
]
PATIENT_TRACTS = {
    pid: _TRACT_ORIENTATIONS[i] for i, pid in enumerate(COHORT_PATIENTS)
}


# --------------------------------------------------------------------------- #
# 3D Tensor Field (Pre-compute face diffusivities)
# --------------------------------------------------------------------------- #
def create_3d_tensor_field_fast(
    grid_size: int = GRID_SIZE,
    dx: float = DX,
    d_white: float = D_WHITE,
    d_gray: float = D_GRAY,
    tract_orientation: np.ndarray = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
    seed: int = 42,
) -> Dict[str, np.ndarray]:
    """Generate 3D tensor field and pre-compute face diffusivities."""
    rng = np.random.default_rng(seed)
    gs = grid_size
    n = tract_orientation

    # Coordinate grids
    x, y, z = np.mgrid[0:gs, 0:gs, 0:gs]
    center = np.array([gs/2, gs/2, gs/2])
    pos = np.stack([x, y, z], axis=-1) - center

    # Distance to tract line
    proj_parallel = np.sum(pos * n, axis=-1, keepdims=True) * n
    dist_perp = np.sqrt(np.sum((pos - proj_parallel) ** 2, axis=-1))
    tract_radius = gs / 3.0
    in_tract = dist_perp < tract_radius

    # Initialize with gray matter baseline
    D_xx = np.full((gs, gs, gs), d_gray, dtype=np.float32)
    D_yy = np.full((gs, gs, gs), d_gray, dtype=np.float32)
    D_zz = np.full((gs, gs, gs), d_gray, dtype=np.float32)
    D_xy = np.zeros((gs, gs, gs), dtype=np.float32)
    D_xz = np.zeros((gs, gs, gs), dtype=np.float32)
    D_yz = np.zeros((gs, gs, gs), dtype=np.float32)

    if np.any(in_tract):
        delta_D = d_white - d_gray
        D_xx[in_tract] = d_gray + delta_D * (n[0] ** 2)
        D_yy[in_tract] = d_gray + delta_D * (n[1] ** 2)
        D_zz[in_tract] = d_gray + delta_D * (n[2] ** 2)
        D_xy[in_tract] = delta_D * n[0] * n[1]
        D_xz[in_tract] = delta_D * n[0] * n[2]
        D_yz[in_tract] = delta_D * n[1] * n[2]

        # Small noise
        noise_scale = 0.02 * d_gray
        mask_sum = int(in_tract.sum())
        D_xx[in_tract] += rng.normal(0, noise_scale, size=mask_sum)
        D_yy[in_tract] += rng.normal(0, noise_scale, size=mask_sum)
        D_zz[in_tract] += rng.normal(0, noise_scale, size=mask_sum)
        D_xy[in_tract] += rng.normal(0, noise_scale * 0.5, size=mask_sum)
        D_xz[in_tract] += rng.normal(0, noise_scale * 0.5, size=mask_sum)
        D_yz[in_tract] += rng.normal(0, noise_scale * 0.5, size=mask_sum)

    # Pre-compute face diffusivities (shifted averages)
    # Dxx_xf: face between i and i+1 (x-direction)
    Dxx_xf = 0.5 * (np.pad(D_xx, ((1, 0), (0, 0), (0, 0)), mode='edge')[:-1] + D_xx)
    Dyy_yf = 0.5 * (np.pad(D_yy, ((0, 0), (1, 0), (0, 0)), mode='edge')[:, :-1] + D_yy)
    Dzz_zf = 0.5 * (np.pad(D_zz, ((0, 0), (0, 0), (1, 0)), mode='edge')[:, :, :-1] + D_zz)

    return {
        "D_xx": D_xx, "D_yy": D_yy, "D_zz": D_zz,
        "D_xy": D_xy, "D_xz": D_xz, "D_yz": D_yz,
        "Dxx_xf": Dxx_xf, "Dyy_yf": Dyy_yf, "Dzz_zf": Dzz_zf,
        "in_tract": in_tract,
    }


# --------------------------------------------------------------------------- #
# Optimized 3D Solver (Fused divergence + step)
# --------------------------------------------------------------------------- #
class AnisotropicFKSolver3DFast:
    """Optimized 3D anisotropic FK solver with fused operations."""

    __slots__ = ("D_xx", "D_yy", "D_zz", "Dxx_xf", "Dyy_yf", "Dzz_zf",
                 "dt", "dx", "rho", "K", "H", "W", "D")

    def __init__(self, tensor_field: Dict, dt: float = DT, dx: float = DX,
                 rho: float = RHO, K: float = K):
        self.D_xx = tensor_field["D_xx"]
        self.D_yy = tensor_field["D_yy"]
        self.D_zz = tensor_field["D_zz"]
        self.Dxx_xf = tensor_field["Dxx_xf"]
        self.Dyy_yf = tensor_field["Dyy_yf"]
        self.Dzz_zf = tensor_field["Dzz_zf"]
        self.dt = float(dt)
        self.dx = float(dx)
        self.rho = float(rho)
        self.K = float(K)
        self.H, self.W, self.D = self.D_xx.shape

        # CFL check (diagonal only)
        max_D = max(self.D_xx.max(), self.D_yy.max(), self.D_zz.max())
        cfl_limit = (self.dx ** 2) / (2.0 * DIM * max_D)
        if self.dt > cfl_limit:
            self.dt = 0.9 * cfl_limit
            print(f"[3D] CFL clamped dt to {self.dt:.4f}")

    def step_fused(self, u: np.ndarray, C: float) -> np.ndarray:
        """Single fused step: divergence + reaction - kill (in-place safe)."""
        dx = self.dx
        dt = self.dt
        rho = self.rho
        K = self.K

        # Pad with zeros (Neumann)
        u_p = np.pad(u, 1, mode="constant", constant_values=0)

        # Flux X: Dxx_xf * (u[i] - u[i-1])/dx at faces
        ux_m = (u_p[1:-1, 1:-1, 1:-1] - u_p[:-2, 1:-1, 1:-1]) / dx
        ux_p = (u_p[2:, 1:-1, 1:-1] - u_p[1:-1, 1:-1, 1:-1]) / dx
        Fx_m = self.Dxx_xf * ux_m
        Fx_p = self.Dxx_xf * ux_p  # same face diffusivity
        div_x = (Fx_p - Fx_m) / dx

        # Flux Y
        uy_m = (u_p[1:-1, 1:-1, 1:-1] - u_p[1:-1, :-2, 1:-1]) / dx
        uy_p = (u_p[1:-1, 2:, 1:-1] - u_p[1:-1, 1:-1, 1:-1]) / dx
        Fy_m = self.Dyy_yf * uy_m
        Fy_p = self.Dyy_yf * uy_p
        div_y = (Fy_p - Fy_m) / dx

        # Flux Z
        uz_m = (u_p[1:-1, 1:-1, 1:-1] - u_p[1:-1, 1:-1, :-2]) / dx
        uz_p = (u_p[1:-1, 1:-1, 2:] - u_p[1:-1, 1:-1, 1:-1]) / dx
        Fz_m = self.Dzz_zf * uz_m
        Fz_p = self.Dzz_zf * uz_p
        div_z = (Fz_p - Fz_m) / dx

        div_term = div_x + div_y + div_z
        react_term = rho * u * (1.0 - u / K)

        # BBB permeability & drug kill
        E_bbb = 0.15 + 0.70 * (u / K)
        C_eff = C * E_bbb
        kill_term = E_MAX * (C_eff ** HILL_COEFF) / (EC50 ** HILL_COEFF + C_eff ** HILL_COEFF + 1e-12) * u

        u_new = u + dt * (div_term + react_term - kill_term)
        return np.clip(u_new, 0.0, K)


# --------------------------------------------------------------------------- #
# Initial Conditions & Drug Schedule (Vectorized)
# --------------------------------------------------------------------------- #
def initial_tumor_seed(grid_size: int, patient_id: str) -> np.ndarray:
    """Initial tumor density field for a real MU-Glioma patient."""
    p = _params[patient_id]
    V0_mm3 = p["V0_mm3"]
    # Convert V0 to density fraction: calibrate so median V0 ≈ 0.15
    M_s0 = float(np.clip(V0_mm3 / 500000.0, 0.05, 0.30))

    # Seed radius in voxels: M_s0 = (4/3 * pi * r³) / grid³
    # => r = (M_s0 * grid³ * 3 / (4*pi))^(1/3)
    r_voxels = (M_s0 * grid_size**3 * 3 / (4 * np.pi)) ** (1 / 3)

    # Build a spherical seed centered in the grid
    center = grid_size / 2
    y, x, z = np.ogrid[0:grid_size, 0:grid_size, 0:grid_size]
    dist_sq = (y - center)**2 + (x - center)**2 + (z - center)**2
    u0 = (dist_sq <= r_voxels**2).astype(np.float64)

    # Normalize so total mass = M_s0
    if u0.sum() > 0:
        u0 = u0 * (M_s0 * grid_size**3) / u0.sum()

    return u0


def tmz_concentration(step: int, dt: float = DT) -> float:
    t_days = step * dt
    day_in_cycle = int(t_days) % CYCLE_DAYS
    if day_in_cycle < DOSE_DAYS_ON:
        return C_PEAK * np.exp(-K_EL * dt)
    days_since_dose = day_in_cycle - (DOSE_DAYS_ON - 1)
    return C_PEAK * np.exp(-K_EL * days_since_dose)


# --------------------------------------------------------------------------- #
# Simulation Runners (Optimized: single loop, pre-allocated)
# --------------------------------------------------------------------------- #
def run_mtd_3d(solver: AnisotropicFKSolver3DFast, u0: np.ndarray) -> Dict:
    """Run MTD with minimal overhead."""
    u = u0.copy()
    vol_hist = []
    mass_hist = []

    for step in range(N_STEPS):
        C = tmz_concentration(step)
        u = solver.step_fused(u, C)

        if step % SAVE_INTERVAL == 0 or step == N_STEPS - 1:
            vol = float(np.sum(u > 0.1)) * (DX ** 3)
            mass = float(u.sum()) * (DX ** 3)
            vol_hist.append((step * DT, vol))
            mass_hist.append((step * DT, mass))

    return {"final_u": u, "volume_history": vol_hist, "mass_history": mass_hist}


def run_adaptive_3d(solver: AnisotropicFKSolver3DFast, u0: np.ndarray) -> Dict:
    """Run adaptive therapy."""
    u = u0.copy()
    baseline_mass = float(u.sum())
    drug_on = True
    vol_hist = []
    mass_hist = []
    drug_on_history = []

    for step in range(N_STEPS):
        current_mass = float(u.sum())

        if drug_on and current_mass < 0.5 * baseline_mass:
            drug_on = False
        elif not drug_on and current_mass > 0.8 * baseline_mass:
            drug_on = True

        C = tmz_concentration(step) if drug_on else 0.0
        u = solver.step_fused(u, C)
        drug_on_history.append(drug_on)

        if step % SAVE_INTERVAL == 0 or step == N_STEPS - 1:
            vol = float(np.sum(u > 0.1)) * (DX ** 3)
            mass = float(u.sum()) * (DX ** 3)
            vol_hist.append((step * DT, vol))
            mass_hist.append((step * DT, mass))

    drug_on_frac = float(np.sum(drug_on_history) / len(drug_on_history))
    return {"final_u": u, "volume_history": vol_hist, "mass_history": mass_hist,
            "drug_on_fraction": drug_on_frac}


# --------------------------------------------------------------------------- #
# Metrics (Fast: no scipy dependency)
# --------------------------------------------------------------------------- #
def compute_sphericity_fast(u: np.ndarray) -> float:
    """Fast 3D sphericity using voxel boundary count."""
    mask = u > 0.1
    vol = float(np.sum(mask))
    if vol < 8:
        return 0.0

    r_eq = (3.0 * vol / (4.0 * np.pi)) ** (1.0 / 3.0)
    sphere_area = 4.0 * np.pi * r_eq ** 2

    # Manual boundary count (no scipy)
    # Boundary = voxels with at least one neighbor that's not tumor
    boundaries = (
        (mask[:-1, :, :] != mask[1:, :, :]).sum() * DX**2 +  # x faces
        (mask[:, :-1, :] != mask[:, 1:, :]).sum() * DX**2 +  # y faces
        (mask[:, :, :-1] != mask[:, :, 1:]).sum() * DX**2    # z faces
    )

    if boundaries < 1e-6:
        return 1.0
    return min(sphere_area / boundaries, 1.0)


# --------------------------------------------------------------------------- #
# Main Execution (Optimized)
# --------------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("Phase 3D Optimized: 3D Volumetric Anisotropic Tumor Growth")
    print(f"8-Patient Real MU-Glioma Cohort ({COHORT_PATIENTS[0]} … {COHORT_PATIENTS[-1]})")
    print(f"Grid: {GRID_SIZE}³, Steps: {N_STEPS}, dt: {DT}")
    print("=" * 70)

    cohort_results = []
    all_volumes = {}

    for pid in COHORT_PATIENTS:
        print(f"\n{'='*70}")
        print(f"Processing {pid}...")
        print(f"{'='*70}")

        rho_s_patient = _params[pid]["rho_per_day"]
        # Single-population FK solver: only rho applies (no separate resistant clone)
        print(f"  rho={rho_s_patient:.6f} /day, V0={_params[pid]['V0_mm3']:.0f} mm³")

        tract_n = PATIENT_TRACTS[pid]
        print(f"Tract orientation: n = [{tract_n[0]:.3f}, {tract_n[1]:.3f}, {tract_n[2]:.3f}]")

        # Stable seed from numeric suffix of PatientID_XXXX
        try:
            seed_extra = int(pid.split("_")[-1])
        except ValueError:
            seed_extra = abs(hash(pid)) % 10_000

        # 1. Generate tensor field (fast)
        print(f"\n[1] Generating 3D tensor field...")
        tf = create_3d_tensor_field_fast(
            grid_size=GRID_SIZE, dx=DX, d_white=D_WHITE, d_gray=D_GRAY,
            tract_orientation=tract_n, seed=42 + seed_extra
        )

        # Quick SPD check at center only
        center = GRID_SIZE // 2
        t = np.array([[tf["D_xx"][center, center, center], tf["D_xy"][center, center, center], tf["D_xz"][center, center, center]],
                      [tf["D_xy"][center, center, center], tf["D_yy"][center, center, center], tf["D_yz"][center, center, center]],
                      [tf["D_xz"][center, center, center], tf["D_yz"][center, center, center], tf["D_zz"][center, center, center]]])
        eigs = np.linalg.eigvalsh(t)
        eig_ratio = float(eigs.max() / eigs.min()) if eigs.min() > 0 else float('inf')
        print(f"    Eigenvalue ratio: {eig_ratio:.1f}x, SPD: {'PASS' if eigs.min() > 0 else 'FAIL'}")

        # 2. Initialize solver with per-patient rho
        print(f"\n[2] Initializing 3D solver...")
        solver = AnisotropicFKSolver3DFast(tf, dt=DT, dx=DX, rho=rho_s_patient, K=K)

        # 3. Initial seed from real V0_mm3
        print(f"\n[3] Planting tumor seed...")
        u0 = initial_tumor_seed(GRID_SIZE, pid)
        init_vol = float(np.sum(u0 > 0.1)) * (DX ** 3)
        print(f"    Initial volume: {init_vol:.2f} mm³")

        # 4. MTD
        print(f"\n[4] Running MTD ({SIM_DAYS} days, {N_STEPS} steps)...")
        res_mtd = run_mtd_3d(solver, u0)
        vol_mtd = float(np.sum(res_mtd["final_u"] > 0.1)) * (DX ** 3)
        mass_mtd = float(res_mtd["final_u"].sum()) * (DX ** 3)
        sph_mtd = compute_sphericity_fast(res_mtd["final_u"])
        print(f"    MTD final volume: {vol_mtd:.2f} mm³")

        # 5. Adaptive
        print(f"\n[5] Running Adaptive...")
        res_adapt = run_adaptive_3d(solver, u0)
        vol_adapt = float(np.sum(res_adapt["final_u"] > 0.1)) * (DX ** 3)
        mass_adapt = float(res_adapt["final_u"].sum()) * (DX ** 3)
        sph_adapt = compute_sphericity_fast(res_adapt["final_u"])
        dose_sparing = 1.0 - res_adapt["drug_on_fraction"]
        print(f"    Adaptive final volume: {vol_adapt:.2f} mm³")
        print(f"    Dose sparing: {dose_sparing:.1%}")

        # Store minimal results
        patient_result = {
            "patient_id": pid,
            "rho_per_day": rho_s_patient,
            "V0_mm3": _params[pid]["V0_mm3"],
            "r_squared": _params[pid]["r_squared"],
            "tract_orientation": tract_n.tolist(),
            "eigenvalue_ratio": eig_ratio,
            "initial_volume_mm3": init_vol,
            "mtd": {"final_volume_mm3": vol_mtd, "final_mass": mass_mtd, "sphericity": sph_mtd},
            "adaptive": {"final_volume_mm3": vol_adapt, "final_mass": mass_adapt,
                         "sphericity": sph_adapt, "drug_on_fraction": res_adapt["drug_on_fraction"],
                         "dose_sparing_fraction": dose_sparing},
        }
        cohort_results.append(patient_result)

        # Save ONLY final densities (not full tensor fields)
        # Keys: {real_pid}_mtd / _adapt / _u0 — suffix pattern required by script 49
        all_volumes[f"{pid}_mtd"] = res_mtd["final_u"]
        all_volumes[f"{pid}_adapt"] = res_adapt["final_u"]
        all_volumes[f"{pid}_u0"] = u0

    # Cohort stats
    print(f"\n{'='*70}")
    print("COHORT AGGREGATE STATISTICS")
    print(f"{'='*70}")
    mtd_v = [r["mtd"]["final_volume_mm3"] for r in cohort_results]
    ad_v = [r["adaptive"]["final_volume_mm3"] for r in cohort_results]
    ds = [r["adaptive"]["dose_sparing_fraction"] for r in cohort_results]
    print(f"MTD: {np.mean(mtd_v):.1f} ± {np.std(mtd_v, ddof=1):.1f} mm³")
    print(f"Adaptive: {np.mean(ad_v):.1f} ± {np.std(ad_v, ddof=1):.1f} mm³")
    print(f"Dose sparing: {np.mean(ds)*100:.1f}% ± {np.std(ds, ddof=1)*100:.1f}%")

    # Save artifacts
    print(f"\n[6] Saving artifacts...")
    np.savez_compressed(OUTPUT_DIR / "3d_master_cohort_volumes.npz", **all_volumes)

    summary = {
        "grid_size": GRID_SIZE,
        "dx_mm": DX,
        "dt_days": DT,
        "sim_days": SIM_DAYS,
        "n_steps": N_STEPS,
        "cohort_size": len(COHORT_PATIENTS),
        "cohort_patients": COHORT_PATIENTS,
        "anisotropy": {"D_parallel": D_WHITE, "D_perpendicular": D_GRAY, "anisotropy_ratio": D_WHITE/D_GRAY},
        "patients": cohort_results,
        "aggregate": {
            "mtd_final_volume_mm3": {"mean": float(np.mean(mtd_v)), "std": float(np.std(mtd_v, ddof=1))},
            "adaptive_final_volume_mm3": {"mean": float(np.mean(ad_v)), "std": float(np.std(ad_v, ddof=1))},
            "dose_sparing_fraction": {"mean": float(np.mean(ds)), "std": float(np.std(ds, ddof=1))},
        },
    }
    with open(OUTPUT_DIR / "3d_extension_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print("3D Extension Complete (Optimized)")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
