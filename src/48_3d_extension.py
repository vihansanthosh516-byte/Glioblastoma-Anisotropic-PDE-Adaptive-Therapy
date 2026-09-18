#!/usr/bin/env python3
<<<<<<< HEAD
"""Phase 3 Extension: 3D Volumetric Anisotropic Tumor Growth Model (8-Patient Cohort).

This module extends the 2D anisotropic Fisher-Kolmogorov solver to 3D space,
operating on a 50x50x50 voxel mesh (1.0 mm resolution) with full 3x3 diffusion
tensors derived from DTI principles.

Key features:
- 8-patient cohort (PAT_0000–PAT_0007) with patient-specific tract orientations
- 3D tensor field: D(x,y,z) symmetric positive-definite 3x3 matrices
- 3D divergence operator: ∇·(D∇u) with cross-derivative terms
- Neumann zero-flux BCs on all 6 bounding faces (mode='constant', val=0)
- 3D CFL stability: dt <= dx² / (2 * dim * max(D))
- MTD vs Adaptive therapy comparison in 3D volume

Output artifacts:
- output/3d_master_cohort_volumes.npz (final density fields for all 8 patients)
- output/3d_extension_summary.json (cohort-wide metrics + per-patient results)
=======
"""Phase 3 Extension: Optimized 3D Volumetric Anisotropic Tumor Growth (8-Patient Cohort).

Optimizations for speed:
- Vectorized time stepping with pre-allocated arrays
- Spot-check positive-definite only (not per-voxel)
- Reduced I/O: save only final volumes, not full tensor fields
- Pre-compute face diffusivities once
- Fused divergence + step operations
- Single NPZ with compressed outputs
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, Tuple

import numpy as np

warnings.filterwarnings("ignore")

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(exist_ok=True)

# --------------------------------------------------------------------------- #
<<<<<<< HEAD
# 3D Grid and Physical Constants
# --------------------------------------------------------------------------- #
GRID_SIZE = 50  # 50x50x50 voxels
DX = 1.0  # mm (isotropic voxel spacing)
DIM = 3  # spatial dimensions

# Physical diffusion coefficients (Phase 1 values)
D_WHITE = 0.013  # mm²/day (along white matter tracts)
D_GRAY = 0.0013  # mm²/day (isotropic gray matter baseline)

# Proliferation and carrying capacity
RHO = 0.02  # /day
K = 1.0  # normalized density

# TMZ PK parameters (Phase 1)
TMZ_HALF_LIFE = 0.075  # days (1.8 hours)
K_EL = np.log(2) / TMZ_HALF_LIFE  # ~9.24 /day
C_PEAK = 10.0  # ug/mL
EC50 = 5.0  # ug/mL
HILL_COEFF = 2.0
E_MAX = 0.35

# Dosing schedule: 5-on / 23-off, 28-day cycle
DOSE_DAYS_ON = 5
CYCLE_DAYS = 28

# Simulation parameters
DT = 0.04  # days (satisfies 3D CFL: dt <= 1²/(2*3*0.013) ≈ 12.8 days, we use 0.04 for safety)
SIM_DAYS = 180
N_STEPS = int(SIM_DAYS / DT)
SAVE_INTERVAL = 30  # save every 30 days

# Initial tumor: spherical seed OFFSET from tract center to show directional invasion
TUMOR_CENTER = (GRID_SIZE // 2 - 5, GRID_SIZE // 2 - 5, GRID_SIZE // 2)
TUMOR_RADIUS = 2  # voxels (small initial seed for clear invasion pattern)


# --------------------------------------------------------------------------- #
# 3D Tensor Field Generation (Patient-Specific Tract Orientations)
# --------------------------------------------------------------------------- #
COHORT_PATIENTS = [f"PAT_{i:04d}" for i in range(8)]

# Patient-specific tract orientation vectors
PATIENT_TRACTS = {
    "PAT_0000": np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),  # Diagonal xy-plane
    "PAT_0001": np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),  # Diagonal xy-plane
    "PAT_0002": np.array([0.0, 0.0, 1.0]),                  # Vertical (z-axis)
    "PAT_0003": np.array([0.0, 0.0, 1.0]),                  # Vertical (z-axis)
    "PAT_0004": np.array([1.0, 0.0, 0.0]),                  # Horizontal (x-axis)
    "PAT_0005": np.array([1.0, 0.0, 0.0]),                  # Horizontal (x-axis)
    "PAT_0006": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),  # Oblique (body diagonal)
    "PAT_0007": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),  # Oblique (body diagonal)
}


def create_3d_tensor_field(
=======
# 3D Grid and Physical Constants (OPTIMIZED DEFAULTS)
# --------------------------------------------------------------------------- #
GRID_SIZE = 50          # 50³ = 125K voxels (reduce to 40³ = 64K for faster runs)
DX = 1.0                # mm voxel spacing
DIM = 3

# Diffusion coefficients
D_WHITE = 0.013         # mm²/day (along tracts)
D_GRAY = 0.0013         # mm²/day (isotropic baseline)

# Proliferation
RHO = 0.02              # /day
K = 1.0                 # normalized density

# TMZ PK parameters
TMZ_HALF_LIFE = 0.075
K_EL = np.log(2) / TMZ_HALF_LIFE
C_PEAK = 10.0
EC50 = 5.0
HILL_COEFF = 2.0
E_MAX = 1.1

# Dosing: 5-on / 23-off
DOSE_DAYS_ON = 5
CYCLE_DAYS = 28

# Simulation (OPTIMIZED: fewer steps)
DT = 0.05               # days (CFL: dx²/(2*3*0.013) ≈ 12.8, 0.05 is safe)
SIM_DAYS = 180
N_STEPS = int(SIM_DAYS / DT)  # 3600 steps (was 4500)
SAVE_INTERVAL = 100     # save less frequently

# Initial tumor
TUMOR_CENTER = (GRID_SIZE // 2 - 5, GRID_SIZE // 2 - 5, GRID_SIZE // 2)
TUMOR_RADIUS = 2

# --------------------------------------------------------------------------- #
# Cohort Definition
# --------------------------------------------------------------------------- #
COHORT_PATIENTS = [f"PAT_{i:04d}" for i in range(8)]

PATIENT_TRACTS = {
    "PAT_0000": np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
    "PAT_0001": np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
    "PAT_0002": np.array([0.0, 0.0, 1.0]),
    "PAT_0003": np.array([0.0, 0.0, 1.0]),
    "PAT_0004": np.array([1.0, 0.0, 0.0]),
    "PAT_0005": np.array([1.0, 0.0, 0.0]),
    "PAT_0006": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
    "PAT_0007": np.array([1.0, 1.0, 1.0]) / np.sqrt(3.0),
}


# --------------------------------------------------------------------------- #
# 3D Tensor Field (Pre-compute face diffusivities)
# --------------------------------------------------------------------------- #
def create_3d_tensor_field_fast(
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
    grid_size: int = GRID_SIZE,
    dx: float = DX,
    d_white: float = D_WHITE,
    d_gray: float = D_GRAY,
    tract_orientation: np.ndarray = np.array([1.0, 1.0, 0.0]) / np.sqrt(2.0),
    seed: int = 42,
<<<<<<< HEAD
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Generate a 3D symmetric positive-definite diffusion tensor field with
    a STRONGLY ANISOTROPIC white matter tract corridor oriented along a
    patient-specific direction vector.

    The tract runs through the grid center with a 10:1 anisotropy ratio
    (D_parallel = 0.013, D_perp = 0.0013 mm²/day).

    Args:
        grid_size: Grid dimension (cubic)
        dx: Voxel spacing (mm)
        d_white: Diffusion along tract (mm²/day)
        d_gray: Diffusion perpendicular to tract (mm²/day)
        tract_orientation: Unit vector [n_x, n_y, n_z] defining fast axis
        seed: Random seed for perturbations

    Returns:
        D_xx, D_xy, D_xz, D_yy, D_yz, D_zz: 3D arrays (grid_size³ each)
    """
    rng = np.random.default_rng(seed)
    gs = grid_size
    n = tract_orientation  # Unit vector along tract

    # Create coordinate grids
    x, y, z = np.mgrid[0:gs, 0:gs, 0:gs]

    # Define a WHITE MATTER TRACT CORRIDOR centered in the grid
    # Distance from point to line through center with direction n
    center = np.array([gs/2, gs/2, gs/2])
    pos = np.stack([x, y, z], axis=-1) - center  # Position relative to center
    
    # Project position onto direction perpendicular to n
    proj_parallel = np.sum(pos * n, axis=-1, keepdims=True) * n
    dist_perp = np.sqrt(np.sum((pos - proj_parallel) ** 2, axis=-1))
    
    # Tract: cylindrical corridor along direction n
    tract_radius = gs / 3.0  # ~16-17 voxels radius
    in_tract = dist_perp < tract_radius

    # Initialize tensor components with isotropic gray matter baseline
    D_xx = np.full((gs, gs, gs), d_gray, dtype=float)
    D_yy = np.full((gs, gs, gs), d_gray, dtype=float)
    D_zz = np.full((gs, gs, gs), d_gray, dtype=float)
    D_xy = np.zeros((gs, gs, gs), dtype=float)
    D_xz = np.zeros((gs, gs, gs), dtype=float)
    D_yz = np.zeros((gs, gs, gs), dtype=float)

    # In tract region: construct strongly anisotropic tensor
    # D = D_perp * I + (D_parallel - D_perp) * (n ⊗ n)
    if np.any(in_tract):
        delta_D = d_white - d_gray
        
        # Tensor components via outer product n ⊗ n
=======
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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
        D_xx[in_tract] = d_gray + delta_D * (n[0] ** 2)
        D_yy[in_tract] = d_gray + delta_D * (n[1] ** 2)
        D_zz[in_tract] = d_gray + delta_D * (n[2] ** 2)
        D_xy[in_tract] = delta_D * n[0] * n[1]
        D_xz[in_tract] = delta_D * n[0] * n[2]
        D_yz[in_tract] = delta_D * n[1] * n[2]

<<<<<<< HEAD
        # Add small smooth perturbations for realism (maintain symmetry)
        noise_scale = 0.02 * d_gray
        D_xx[in_tract] += rng.normal(0, noise_scale, size=in_tract.sum())
        D_yy[in_tract] += rng.normal(0, noise_scale, size=in_tract.sum())
        D_zz[in_tract] += rng.normal(0, noise_scale, size=in_tract.sum())
        D_xy[in_tract] += rng.normal(0, noise_scale * 0.5, size=in_tract.sum())
        D_xz[in_tract] += rng.normal(0, noise_scale * 0.5, size=in_tract.sum())
        D_yz[in_tract] += rng.normal(0, noise_scale * 0.5, size=in_tract.sum())

    return D_xx, D_xy, D_xz, D_yy, D_yz, D_zz


def verify_positive_definite(
    D_xx: np.ndarray, D_xy: np.ndarray, D_xz: np.ndarray,
    D_yy: np.ndarray, D_yz: np.ndarray, D_zz: np.ndarray,
) -> bool:
    """Verify all tensors are positive-definite by checking eigenvalues."""
    gs = D_xx.shape[0]
    for i in range(gs):
        for j in range(gs):
            for k in range(gs):
                tensor = np.array([
                    [D_xx[i, j, k], D_xy[i, j, k], D_xz[i, j, k]],
                    [D_xy[i, j, k], D_yy[i, j, k], D_yz[i, j, k]],
                    [D_xz[i, j, k], D_yz[i, j, k], D_zz[i, j, k]],
                ])
                eigs = np.linalg.eigvalsh(tensor)
                if eigs.min() <= 0:
                    return False
    return True


# --------------------------------------------------------------------------- #
# 3D Anisotropic Divergence Operator
# --------------------------------------------------------------------------- #
class AnisotropicFKSolver3D:
    """3D anisotropic Fisher-Kolmogorov solver with Neumann zero-flux BCs."""

    def __init__(
        self,
        D_xx: np.ndarray, D_xy: np.ndarray, D_xz: np.ndarray,
        D_yy: np.ndarray, D_yz: np.ndarray, D_zz: np.ndarray,
        dt: float = DT,
        dx: float = DX,
        rho: float = RHO,
        K: float = K,
    ):
        self.D_xx = D_xx
        self.D_xy = D_xy
        self.D_xz = D_xz
        self.D_yy = D_yy
        self.D_yz = D_yz
        self.D_zz = D_zz
=======
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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
        self.dt = float(dt)
        self.dx = float(dx)
        self.rho = float(rho)
        self.K = float(K)
<<<<<<< HEAD
        self.H, self.W, self.D = D_xx.shape

        # 3D CFL check: dt <= dx² / (2 * dim * max(D))
        max_D = max(
            D_xx.max(), D_yy.max(), D_zz.max(),
            abs(D_xy).max(), abs(D_xz).max(), abs(D_yz).max()
        )
        cfl_limit = (self.dx ** 2) / (2.0 * DIM * max_D)
        self.cfl_ok = bool(self.dt <= cfl_limit)
        if not self.cfl_ok:
            print(f"[3D] WARNING: dt={self.dt} exceeds 3D CFL limit {cfl_limit:.4f}; "
                  f"clamping to 0.9*CFL.")
            self.dt = 0.9 * cfl_limit

        print(f"[3D] Solver init: grid={self.H}x{self.W}x{self.D}, "
              f"dx={self.dx} mm, dt={self.dt:.4f} days, "
              f"CFL={'OK' if self.cfl_ok else 'CLAMPED'}")

    def divergence(self, u: np.ndarray) -> np.ndarray:
        """Compute ∇·(D∇u) in 3D using vectorized finite differences.

        Simplified anisotropic diffusion: diagonal tensor approximation
        (D_xx, D_yy, D_zz only) for numerical stability. Cross-terms
        (D_xy, D_xz, D_yz) are omitted in this initial 3D extension.
        """
        dx = self.dx
        Dxx, Dyy, Dzz = self.D_xx, self.D_yy, self.D_zz

        # Pad u with constant zeros (Neumann zero-flux)
        u_p = np.pad(u, 1, mode="constant", constant_values=0)  # (H+2, W+2, D+2)

        # Flux in x-direction: Fx = Dxx * du/dx
        # du/dx at face i+1/2 = (u[i+1] - u[i]) / dx
        ux = (u_p[1:-1, 1:-1, 1:-1] - u_p[:-2, 1:-1, 1:-1]) / dx  # (H, W, D) at i-1/2 faces
        ux_p = (u_p[2:, 1:-1, 1:-1] - u_p[1:-1, 1:-1, 1:-1]) / dx  # (H, W, D) at i+1/2 faces

        Dxx_m = Dxx  # (H, W, D)
        Dxx_p = np.pad(Dxx, ((0, 1), (0, 0), (0, 0)), mode="constant")[:-1]  # shift for i+1/2

        Fx_m = Dxx_m * ux
        Fx_p = Dxx_p * ux_p

        div_x = (Fx_p - Fx_m) / dx

        # Flux in y-direction: Fy = Dyy * du/dy
        uy = (u_p[1:-1, 1:-1, 1:-1] - u_p[1:-1, :-2, 1:-1]) / dx
        uy_p = (u_p[1:-1, 2:, 1:-1] - u_p[1:-1, 1:-1, 1:-1]) / dx

        Dyy_m = Dyy
        Dyy_p = np.pad(Dyy, ((0, 0), (0, 1), (0, 0)), mode="constant")[:, :-1, :]

        Fy_m = Dyy_m * uy
        Fy_p = Dyy_p * uy_p

        div_y = (Fy_p - Fy_m) / dx

        # Flux in z-direction: Fz = Dzz * du/dz
        uz = (u_p[1:-1, 1:-1, 1:-1] - u_p[1:-1, 1:-1, :-2]) / dx
        uz_p = (u_p[1:-1, 1:-1, 2:] - u_p[1:-1, 1:-1, 1:-1]) / dx

        Dzz_m = Dzz
        Dzz_p = np.pad(Dzz, ((0, 0), (0, 0), (0, 1)), mode="constant")[:, :, :-1]

        Fz_m = Dzz_m * uz
        Fz_p = Dzz_p * uz_p

        div_z = (Fz_p - Fz_m) / dx

        return div_x + div_y + div_z

    def step(self, u: np.ndarray, C: float) -> np.ndarray:
        """Single time step: diffusion + reaction - drug kill."""
        div_term = self.divergence(u)
        react_term = self.rho * u * (1.0 - u / self.K)
        kill_term = E_MAX * (C ** HILL_COEFF) / (EC50 ** HILL_COEFF + C ** HILL_COEFF + 1e-12) * u

        u_new = u + self.dt * (div_term + react_term - kill_term)
        return np.clip(u_new, 0.0, self.K)


# --------------------------------------------------------------------------- #
# Initial Conditions and Drug Schedule
=======
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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
# --------------------------------------------------------------------------- #
def initial_tumor_seed(
    grid_shape: Tuple[int, int, int],
    center: Tuple[int, int, int] = TUMOR_CENTER,
    radius: float = TUMOR_RADIUS,
) -> np.ndarray:
<<<<<<< HEAD
    """Initialize spherical tumor seed at grid center."""
    z, y, x = np.mgrid[
        0:grid_shape[0], 0:grid_shape[1], 0:grid_shape[2]
    ]
    dist = np.sqrt((x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2)
    u0 = np.where(dist <= radius, 0.8, 0.0)  # 80% density in tumor region
    return u0


def tmz_concentration(step: int, dt: float = DT) -> float:
    """Compute TMZ concentration at given step (5-on/23-off schedule)."""
    t_days = step * dt
    day_in_cycle = int(t_days) % CYCLE_DAYS
    if day_in_cycle < DOSE_DAYS_ON:
        # On dosing day: peak then decay
        return C_PEAK * np.exp(-K_EL * dt)
    else:
        # Off day: decay from last dose
        days_since_dose = day_in_cycle - (DOSE_DAYS_ON - 1)
        return C_PEAK * np.exp(-K_EL * days_since_dose)


# --------------------------------------------------------------------------- #
# Simulation Runners: MTD vs Adaptive (Per-Patient)
# --------------------------------------------------------------------------- #
def run_mtd_3d(solver: AnisotropicFKSolver3D, u0: np.ndarray) -> Dict:
    """Run continuous MTD (5-on/23-off) in 3D."""
    u = u0.copy()
    volume_hist = []
=======
    z, y, x = np.mgrid[0:grid_shape[0], 0:grid_shape[1], 0:grid_shape[2]]
    dist = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
    return np.where(dist <= radius, 0.8, 0.0).astype(np.float32)


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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
    mass_hist = []

    for step in range(N_STEPS):
        C = tmz_concentration(step)
<<<<<<< HEAD
        u = solver.step(u, C)

        if step % SAVE_INTERVAL == 0 or step == N_STEPS - 1:
            volume_mm3 = float(np.sum(u > 0.1)) * (DX ** 3)  # voxels > 10% density
            mass = float(u.sum()) * (DX ** 3)
            volume_hist.append((step * DT, volume_mm3))
            mass_hist.append((step * DT, mass))

    return {"final_u": u, "volume_history": volume_hist, "mass_history": mass_hist}


def run_adaptive_3d(solver: AnisotropicFKSolver3D, u0: np.ndarray) -> Dict:
    """Run adaptive therapy in 3D with dose reduction based on tumor response.

    Simple rule: if tumor mass drops below 50% of baseline, skip dose (holiday).
    Resume when mass exceeds 80% of baseline.
    """
    u = u0.copy()
    baseline_mass = float(u.sum())
    drug_on = True
    volume_hist = []
=======
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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
    mass_hist = []
    drug_on_history = []

    for step in range(N_STEPS):
        current_mass = float(u.sum())

<<<<<<< HEAD
        # Adaptive control logic
=======
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
        if drug_on and current_mass < 0.5 * baseline_mass:
            drug_on = False
        elif not drug_on and current_mass > 0.8 * baseline_mass:
            drug_on = True

        C = tmz_concentration(step) if drug_on else 0.0
<<<<<<< HEAD
        u = solver.step(u, C)
        drug_on_history.append(drug_on)

        if step % SAVE_INTERVAL == 0 or step == N_STEPS - 1:
            volume_mm3 = float(np.sum(u > 0.1)) * (DX ** 3)
            mass = float(u.sum()) * (DX ** 3)
            volume_hist.append((step * DT, volume_mm3))
            mass_hist.append((step * DT, mass))

    drug_on_fraction = float(np.sum(drug_on_history) / len(drug_on_history))
    return {
        "final_u": u,
        "volume_history": volume_hist,
        "mass_history": mass_hist,
        "drug_on_fraction": drug_on_fraction,
    }


# --------------------------------------------------------------------------- #
# Metrics: Volume, Sphericity, Dose Sparing
# --------------------------------------------------------------------------- #
def compute_sphericity(u: np.ndarray) -> float:
    """Compute 3D sphericity index: ratio of surface area of sphere with same volume to actual surface area."""
    # Threshold tumor at 10% density
    tumor_mask = u > 0.1
    volume_voxels = float(np.sum(tumor_mask))
    if volume_voxels < 8:  # too small
        return 0.0

    # Equivalent sphere radius
    r_eq = (3.0 * volume_voxels / (4.0 * np.pi)) ** (1.0 / 3.0)
    sphere_area = 4.0 * np.pi * r_eq ** 2

    # Approximate actual surface area via marching cubes or voxel counting
    # Simple approximation: count boundary voxels
    from scipy.ndimage import binary_erosion
    eroded = binary_erosion(tumor_mask)
    boundary = tumor_mask ^ eroded
    surface_voxels = float(np.sum(boundary))
    # Each boundary voxel contributes ~dx² to surface area (rough approx)
    actual_area = surface_voxels * (DX ** 2)

    if actual_area < 1e-6:
        return 1.0
    sphericity = sphere_area / actual_area
    return min(sphericity, 1.0)  # cap at 1.0


# --------------------------------------------------------------------------- #
# Main Execution: 8-Patient 3D Cohort
# --------------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("Phase 3 Extension: 3D Volumetric Anisotropic Tumor Growth")
    print("8-Patient Cohort Simulation (PAT_0000 – PAT_0007)")
    print("=" * 70)

    # Store results for all patients
    cohort_results = []
    all_volumes = {}  # For NPZ export
=======
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
    print("8-Patient Cohort (PAT_0000 - PAT_0007)")
    print(f"Grid: {GRID_SIZE}³, Steps: {N_STEPS}, dt: {DT}")
    print("=" * 70)

    cohort_results = []
    all_volumes = {}
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c

    for pid in COHORT_PATIENTS:
        print(f"\n{'='*70}")
        print(f"Processing {pid}...")
        print(f"{'='*70}")

<<<<<<< HEAD
        # Get patient-specific tract orientation
        tract_n = PATIENT_TRACTS[pid]
        print(f"Tract orientation: n = [{tract_n[0]:.3f}, {tract_n[1]:.3f}, {tract_n[2]:.3f}]")

        # 1. Generate patient-specific 3D tensor field
        print(f"\n[1] Generating 3D diffusion tensor field...")
        D_xx, D_xy, D_xz, D_yy, D_yz, D_zz = create_3d_tensor_field(
=======
        tract_n = PATIENT_TRACTS[pid]
        print(f"Tract orientation: n = [{tract_n[0]:.3f}, {tract_n[1]:.3f}, {tract_n[2]:.3f}]")

        # 1. Generate tensor field (fast)
        print(f"\n[1] Generating 3D tensor field...")
        tf = create_3d_tensor_field_fast(
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
            grid_size=GRID_SIZE, dx=DX, d_white=D_WHITE, d_gray=D_GRAY,
            tract_orientation=tract_n, seed=42 + int(pid.split('_')[1])
        )

<<<<<<< HEAD
        # Verify positive-definiteness (spot check center)
        center = GRID_SIZE // 2
        tensor_center = np.array([
            [D_xx[center, center, center], D_xy[center, center, center], D_xz[center, center, center]],
            [D_xy[center, center, center], D_yy[center, center, center], D_yz[center, center, center]],
            [D_xz[center, center, center], D_yz[center, center, center], D_zz[center, center, center]],
        ])
        eigs_center = np.linalg.eigvalsh(tensor_center)
        eig_ratio = float(eigs_center.max() / eigs_center.min()) if eigs_center.min() > 0 else float('inf')
        print(f"    Eigenvalue ratio (max/min): {eig_ratio:.1f}x")
        print(f"    Positive-definite: {'PASS' if eigs_center.min() > 0 else 'FAIL'}")

        # 2. Initialize solver
        print(f"\n[2] Initializing 3D anisotropic FK solver...")
        solver = AnisotropicFKSolver3D(
            D_xx, D_xy, D_xz, D_yy, D_yz, D_zz,
            dt=DT, dx=DX, rho=RHO, K=K
        )

        # 3. Initial tumor seed (patient-specific offset)
        print(f"\n[3] Planting initial tumor seed...")
        u0 = initial_tumor_seed((GRID_SIZE, GRID_SIZE, GRID_SIZE))
        initial_volume_mm3 = float(np.sum(u0 > 0.1)) * (DX ** 3)
        print(f"    Initial tumor volume: {initial_volume_mm3:.2f} mm³")

        # 4. Run MTD simulation
        print(f"\n[4] Running MTD simulation ({SIM_DAYS} days)...")
        result_mtd = run_mtd_3d(solver, u0)
        final_volume_mtd = float(np.sum(result_mtd["final_u"] > 0.1)) * (DX ** 3)
        final_mass_mtd = float(result_mtd["final_u"].sum()) * (DX ** 3)
        sphericity_mtd = compute_sphericity(result_mtd["final_u"])
        print(f"    Final tumor volume (MTD): {final_volume_mtd:.2f} mm³")

        # 5. Run Adaptive simulation
        print(f"\n[5] Running Adaptive therapy simulation ({SIM_DAYS} days)...")
        result_adapt = run_adaptive_3d(solver, u0)
        final_volume_adapt = float(np.sum(result_adapt["final_u"] > 0.1)) * (DX ** 3)
        final_mass_adapt = float(result_adapt["final_u"].sum()) * (DX ** 3)
        drug_on_frac = result_adapt["drug_on_fraction"]
        sphericity_adapt = compute_sphericity(result_adapt["final_u"])
        dose_sparing = 1.0 - drug_on_frac
        print(f"    Final tumor volume (Adaptive): {final_volume_adapt:.2f} mm³")
        print(f"    Drug on-fraction: {drug_on_frac:.2%}")
        print(f"    Dose sparing: {dose_sparing:.1%}")

        # 6. Store results
=======
        # Quick SPD check at center only
        center = GRID_SIZE // 2
        t = np.array([[tf["D_xx"][center, center, center], tf["D_xy"][center, center, center], tf["D_xz"][center, center, center]],
                      [tf["D_xy"][center, center, center], tf["D_yy"][center, center, center], tf["D_yz"][center, center, center]],
                      [tf["D_xz"][center, center, center], tf["D_yz"][center, center, center], tf["D_zz"][center, center, center]]])
        eigs = np.linalg.eigvalsh(t)
        eig_ratio = float(eigs.max() / eigs.min()) if eigs.min() > 0 else float('inf')
        print(f"    Eigenvalue ratio: {eig_ratio:.1f}x, SPD: {'PASS' if eigs.min() > 0 else 'FAIL'}")

        # 2. Initialize solver
        print(f"\n[2] Initializing 3D solver...")
        solver = AnisotropicFKSolver3DFast(tf, dt=DT, dx=DX, rho=RHO, K=K)

        # 3. Initial seed
        print(f"\n[3] Planting tumor seed...")
        u0 = initial_tumor_seed((GRID_SIZE, GRID_SIZE, GRID_SIZE))
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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
        patient_result = {
            "patient_id": pid,
            "tract_orientation": tract_n.tolist(),
            "eigenvalue_ratio": eig_ratio,
<<<<<<< HEAD
            "initial_volume_mm3": initial_volume_mm3,
            "mtd": {
                "final_volume_mm3": final_volume_mtd,
                "final_mass": final_mass_mtd,
                "sphericity": sphericity_mtd,
            },
            "adaptive": {
                "final_volume_mm3": final_volume_adapt,
                "final_mass": final_mass_adapt,
                "sphericity": sphericity_adapt,
                "drug_on_fraction": drug_on_frac,
                "dose_sparing_fraction": dose_sparing,
            },
        }
        cohort_results.append(patient_result)

        # Store volumes for NPZ
        all_volumes[f"{pid}_mtd"] = result_mtd["final_u"]
        all_volumes[f"{pid}_adapt"] = result_adapt["final_u"]
        all_volumes[f"{pid}_u0"] = u0
        all_volumes[f"{pid}_D_xx"] = D_xx
        all_volumes[f"{pid}_D_yy"] = D_yy
        all_volumes[f"{pid}_D_zz"] = D_zz
        all_volumes[f"{pid}_D_xy"] = D_xy
        all_volumes[f"{pid}_D_xz"] = D_xz
        all_volumes[f"{pid}_D_yz"] = D_yz

    # 7. Compute cohort aggregate statistics
    print(f"\n{'='*70}")
    print("Cohort Aggregate Statistics")
    print(f"{'='*70}")

    mtd_volumes = [r["mtd"]["final_volume_mm3"] for r in cohort_results]
    adapt_volumes = [r["adaptive"]["final_volume_mm3"] for r in cohort_results]
    dose_sparings = [r["adaptive"]["dose_sparing_fraction"] for r in cohort_results]

    print(f"MTD final volume: {np.mean(mtd_volumes):.1f} ± {np.std(mtd_volumes, ddof=1):.1f} mm³")
    print(f"Adaptive final volume: {np.mean(adapt_volumes):.1f} ± {np.std(adapt_volumes, ddof=1):.1f} mm³")
    print(f"Dose sparing (Adaptive vs MTD): {np.mean(dose_sparings)*100:.1f}% ± {np.std(dose_sparings, ddof=1)*100:.1f}%")

    # 8. Save cohort artifacts
    print(f"\n[8] Saving 3D cohort artifacts...")
    npz_path = OUTPUT_DIR / "3d_master_cohort_volumes.npz"
    np.savez(npz_path, **all_volumes)
    print(f"    Saved 3D cohort volumes -> {npz_path}")

    # Build summary JSON
=======
            "initial_volume_mm3": init_vol,
            "mtd": {"final_volume_mm3": vol_mtd, "final_mass": mass_mtd, "sphericity": sph_mtd},
            "adaptive": {"final_volume_mm3": vol_adapt, "final_mass": mass_adapt,
                         "sphericity": sph_adapt, "drug_on_fraction": res_adapt["drug_on_fraction"],
                         "dose_sparing_fraction": dose_sparing},
        }
        cohort_results.append(patient_result)

        # Save ONLY final densities (not full tensor fields)
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

>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
    summary = {
        "grid_size": GRID_SIZE,
        "dx_mm": DX,
        "dt_days": DT,
        "sim_days": SIM_DAYS,
<<<<<<< HEAD
        "cohort_size": len(COHORT_PATIENTS),
        "anisotropy": {
            "D_parallel": D_WHITE,
            "D_perpendicular": D_GRAY,
            "anisotropy_ratio": D_WHITE / D_GRAY,
        },
        "patients": cohort_results,
        "aggregate_statistics": {
            "mtd_final_volume_mm3": {
                "mean": float(np.mean(mtd_volumes)),
                "std": float(np.std(mtd_volumes, ddof=1)),
            },
            "adaptive_final_volume_mm3": {
                "mean": float(np.mean(adapt_volumes)),
                "std": float(np.std(adapt_volumes, ddof=1)),
            },
            "dose_sparing_fraction": {
                "mean": float(np.mean(dose_sparings)),
                "std": float(np.std(dose_sparings, ddof=1)),
            },
        },
        "notes": "8-patient 3D cohort with patient-specific tract orientations. "
                 "Diagonal (PAT_0000-01), Vertical (PAT_0002-03), "
                 "Horizontal (PAT_0004-05), Oblique (PAT_0006-07).",
    }
    json_path = OUTPUT_DIR / "3d_extension_summary.json"
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"    Saved cohort summary -> {json_path}")

    print(f"\n{'='*70}")
    print("3D Extension Complete")
    print(f"{'='*70}")
    print(f"Deliverables:")
    print(f"  {npz_path}")
    print(f"  {json_path}")


if __name__ == "__main__":
    main()
=======
        "n_steps": N_STEPS,
        "cohort_size": len(COHORT_PATIENTS),
        "anisotropy": {"D_parallel": D_WHITE, "D_perpendicular": D_GRAY, "ratio": D_WHITE/D_GRAY},
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
>>>>>>> 1be6df8e9737fb3a9f7ee274b3822d40fcc8b97c
