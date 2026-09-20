#!/usr/bin/env python3
"""
Month 3, Week 2: Fisher-Kolmogorov PDE Solver — Finite Difference (Non-Periodic)
Morphogen diffusion with reaction term: ∂ρ/∂t = D∇²ρ + rρ(1-ρ)

Calibrated for clinical glioblastoma invasion velocity: 10-50 µm/hr
Uses explicit finite differences with Dirichlet BCs (ρ=0 at boundaries)
to eliminate periodic boundary artifacts.

Key fixes from spectral version:
- Non-periodic boundaries (Dirichlet ρ=0) — no ghost mass at corners
- Radial-average front detector — immune to 2D filling artifacts
- Longer run time (t=5.0) — front moves 20 px, measurable above noise
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Tuple

import numpy as np
import torch


# ---------------------------------------------------------------------------
# Project paths (absolute; independent of cwd)
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------
def initialize_morphogen_field(
    grid_size: Tuple[int, int],
    device: torch.device,
    initial_conditions: str = "tumor_core",
) -> torch.Tensor:
    """Initialize morphogen concentration field."""
    H, W = grid_size
    rho = torch.zeros(H, W, device=device, dtype=torch.float32)

    if initial_conditions == "tumor_core":
        center_h, center_w = H // 2, W // 2
        y = torch.arange(H, device=device).float() - center_h
        x = torch.arange(W, device=device).float() - center_w
        Y, X = torch.meshgrid(y, x, indexing="ij")
        dist_sq = X**2 + Y**2
        rho = torch.exp(-dist_sq / (2 * 15**2)) * 0.8  # Smooth Gaussian, no Gibbs

    return rho


# ---------------------------------------------------------------------------
# Finite Difference Solver (Explicit, Dirichlet BCs)
# ---------------------------------------------------------------------------
class FDSolverFK:
    """
    Explicit finite difference solver for Fisher-Kolmogorov PDE:
        ∂ρ/∂t = D∇²ρ + rρ(1-ρ)

    Spatial discretization: 5-point stencil Laplacian
        ∇²ρ ≈ (ρ_{i+1,j} + ρ_{i-1,j} + ρ_{i,j+1} + ρ_{i,j-1} - 4ρ_{i,j}) / dx²

    Time stepping: Forward Euler (simple, stable for dt ≤ dx²/(4D))
    Boundary conditions: Dirichlet ρ=0 at all edges (tumor cannot escape domain)
    """

    def __init__(
        self,
        grid_size: Tuple[int, int],
        D: float = 10.0,
        r: float = 0.4,
        dt: float = 0.001,
        n_steps: int = 5000,
        save_interval: int = 500,
        device: torch.device = None,
    ):
        self.H, self.W = grid_size
        self.grid_size = grid_size
        self.D = D
        self.r = r
        self.dt = dt
        self.n_steps = n_steps
        self.save_interval = save_interval
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # Spatial step (dx = 1 pixel = 5 µm)
        self.dx = 1.0
        self.dx2 = self.dx * self.dx

        # Diffusion coefficient for finite difference stencil
        self.diff_coeff = D * dt / self.dx2

        # Stability check
        max_stable_dt = self.dx2 / (4.0 * D)
        if dt > max_stable_dt:
            print(f"[WARNING] dt={dt} exceeds stability limit {max_stable_dt:.6f}")
            print(f"          Consider reducing dt for stability")

        # Precompute Laplacian stencil indices for vectorized update
        # We'll use slicing for the 5-point stencil
        self._setup_indices()

        print(f"[FK-FD] Finite Difference solver: {grid_size[0]}x{grid_size[1]}, D={D}, r={r}, dt={dt}")
        print(f"[FK-FD] Analytical c = 2*sqrt(D*r) = {2 * np.sqrt(D * r):.4f} px/unit-time")
        print(f"[FK-FD] Clinical velocity: {2 * np.sqrt(D * r) * 5:.1f} µm/hr")
        print(f"[FK-FD] Diffusion CFL: dt_max = {max_stable_dt:.6f}, using dt={dt} ({'OK' if dt <= max_stable_dt else 'UNSTABLE'})")
        print(f"[FK-FD] Total simulated time: {n_steps * dt:.2f} unit-time")
        print(f"[FK-FD] Expected front travel: {2 * np.sqrt(D * r) * n_steps * dt:.1f} px")

    def _setup_indices(self):
        """Precompute slicing indices for vectorized 5-point Laplacian."""
        H, W = self.H, self.W
        # Interior indices (excluding boundaries where Dirichlet BC applies)
        self.interior = (slice(1, H-1), slice(1, W-1))
        # Neighbors for 5-point stencil
        self.up = (slice(0, H-2), slice(1, W-1))
        self.down = (slice(2, H), slice(1, W-1))
        self.left = (slice(1, H-1), slice(0, W-2))
        self.right = (slice(1, H-1), slice(2, W))

    def reaction(self, rho: torch.Tensor) -> torch.Tensor:
        """N(ρ) = r * ρ * (1 - ρ)"""
        return self.r * rho * (1.0 - rho)

    def laplacian(self, rho: torch.Tensor) -> torch.Tensor:
        """5-point stencil Laplacian with Dirichlet BC (ρ=0 at boundaries)."""
        # Result array (same shape as rho)
        lap = torch.zeros_like(rho)

        # Interior points only (boundaries stay zero = Dirichlet)
        lap[self.interior] = (
            rho[self.up] + rho[self.down] + rho[self.left] + rho[self.right] - 4.0 * rho[self.interior]
        ) / self.dx2

        return lap

    def rhs(self, rho: torch.Tensor) -> torch.Tensor:
        """Right-hand side: D∇²ρ + rρ(1-ρ)"""
        diff_term = self.D * self.laplacian(rho)
        react_term = self.reaction(rho)
        return diff_term + react_term

    def step(self, rho: torch.Tensor) -> torch.Tensor:
        """Single RK4 time step for better accuracy (less numerical diffusion)."""
        # RK4 stages
        k1 = self.rhs(rho)
        k2 = self.rhs(rho + 0.5 * self.dt * k1)
        k3 = self.rhs(rho + 0.5 * self.dt * k2)
        k4 = self.rhs(rho + self.dt * k3)

        rho_new = rho + (self.dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # Enforce Dirichlet BC: ρ=0 at boundaries
        rho_new[0, :] = 0.0
        rho_new[-1, :] = 0.0
        rho_new[:, 0] = 0.0
        rho_new[:, -1] = 0.0

        # Clamp to [0, 1] for physical consistency
        rho_new = torch.clamp(rho_new, 0.0, 1.0)

        # NaN/inf guard
        if not torch.isfinite(rho_new).all():
            n_bad = int((~torch.isfinite(rho_new)).sum().item())
            raise RuntimeError(
                f"Non-finite values ({n_bad}) detected in rho — aborting. Reduce dt."
            )

        return rho_new

    def detect_front_radial(self, field: torch.Tensor, threshold: float = 0.1) -> float:
        """
        Robust front detector: radial average profile.

        Tracks the LEADING EDGE (threshold=0.1) where analytical speed c=2√(Dr) applies.
        The 50% level (threshold=0.5) lags behind and underestimates speed.
        """
        H, W = self.H, self.W
        y = torch.arange(H, device=self.device).float() - H // 2
        x = torch.arange(W, device=self.device).float() - W // 2
        Y, X = torch.meshgrid(y, x, indexing='ij')
        dist = torch.sqrt(X**2 + Y**2).flatten()
        vals = field.flatten()

        # Bin by radius and average
        r_max = int(dist.max().item())
        bins = torch.arange(0, r_max + 1, device=self.device).float()
        profile = torch.zeros(len(bins) - 1, device=self.device)

        for i in range(len(bins) - 1):
            mask = (dist >= bins[i]) & (dist < bins[i + 1])
            if mask.any():
                profile[i] = vals[mask].mean()

        # Find largest radius where profile > threshold
        above = (profile > threshold).nonzero(as_tuple=True)[0]
        if len(above) == 0:
            return 0.0
        return float(bins[above[-1] + 1].item())

    def detect_front_mass(self, field: torch.Tensor, threshold: float = 0.1) -> float:
        """
        A2: Front position = radius of disk with same mass above threshold.

        Robust integral measurement: total mass above threshold ≈ π * r_front².
        Use low threshold (0.1) to track leading edge of wave.
        """
        mass = (field > threshold).float().sum().item()
        if mass == 0:
            return 0.0
        return np.sqrt(mass / np.pi)

    def detect_front_axis(self, field: torch.Tensor, threshold: float = 0.5) -> float:
        """
        A3: Track front along central x-axis (1D slice).

        Find rightmost x where field crosses threshold on y=center line.
        Robust to 2D filling artifacts.
        """
        mid = self.H // 2
        line = field[mid, :].cpu().numpy()
        crossings = np.where(line > threshold)[0]
        if len(crossings) == 0:
            return 0.0
        return float(crossings.max() - self.W // 2)

    def run(
        self,
        n_steps: int = None,
        save_interval: int = None,
        initial_conditions: str = "tumor_core",
    ) -> Tuple[np.ndarray, dict]:
        if n_steps is None:
            n_steps = self.n_steps
        if save_interval is None:
            save_interval = self.save_interval

        rho = initialize_morphogen_field(self.grid_size, self.device, initial_conditions)

        snapshots = [rho.cpu().numpy()]
        front_positions = [self.detect_front_radial(rho, threshold=0.1)]  # Track leading edge
        front_positions_mass = [self.detect_front_mass(rho, threshold=0.1)]
        front_positions_axis = [self.detect_front_axis(rho, threshold=0.1)]  # Lower threshold
        times = [0.0]

        t0 = time.perf_counter()
        aborted_at = None

        for step in range(1, n_steps + 1):
            try:
                rho = self.step(rho)
            except RuntimeError as e:
                print(f"[FK-FD] ABORT at step {step}: {e}")
                aborted_at = step
                break

            # Record front position every 10 steps (for speed vs accuracy tradeoff)
            if step % 10 == 0:
                front_positions.append(self.detect_front_radial(rho, threshold=0.1))
                front_positions_mass.append(self.detect_front_mass(rho, threshold=0.1))
                front_positions_axis.append(self.detect_front_axis(rho, threshold=0.1))
                times.append(step * self.dt)

            if step % save_interval == 0:
                snapshots.append(rho.cpu().numpy())

            if step % 1000 == 0:
                elapsed = time.perf_counter() - t0
                print(f"[FK-FD] Step {step}/{n_steps} (t={step*self.dt:.2f}) front={front_positions[-1]:.2f} ({elapsed:.1f}s)")

        elapsed = time.perf_counter() - t0
        print(f"[FK-FD] Completed in {elapsed:.2f}s")
        if aborted_at is not None:
            print(f"[FK-FD] Run was aborted at step {aborted_at}; metrics reflect partial run.")

        # Wave speed from linear fit of front position vs time
        # Use the radial detector (most robust)
        if len(front_positions) > 1:
            times_np = np.array(times)
            fronts_np = np.array(front_positions)
            n = len(times_np)
            lo, hi = max(1, int(0.10 * n)), max(2, int(0.90 * n))
            coeffs = np.polyfit(times_np[lo:hi], fronts_np[lo:hi], 1)
            numerical_wave_speed = float(coeffs[0])  # px/unit-time
            analytical = 2 * np.sqrt(self.D * self.r)
            speed_error = abs(numerical_wave_speed - analytical) / analytical * 100
        else:
            numerical_wave_speed = 0.0
            analytical = 2 * np.sqrt(self.D * self.r)
            speed_error = 0.0

        # Also compute speeds from other detectors for comparison
        def compute_speed(times_np, fronts_np):
            if len(fronts_np) < 2:
                return 0.0, 0.0
            n = len(times_np)
            lo, hi = max(1, int(0.10 * n)), max(2, int(0.90 * n))
            coeffs = np.polyfit(times_np[lo:hi], fronts_np[lo:hi], 1)
            return float(coeffs[0]), abs(float(coeffs[0]) - analytical) / analytical * 100

        times_np = np.array(times)
        speed_mass, err_mass = compute_speed(times_np, np.array(front_positions_mass))
        speed_axis, err_axis = compute_speed(times_np, np.array(front_positions_axis))

        metrics = {
            "grid_size": list(self.grid_size),
            "n_steps": n_steps,
            "steps_completed": len(times) * 10 if len(times) > 1 else 0,
            "aborted": aborted_at is not None,
            "aborted_at_step": aborted_at,
            "D": self.D,
            "r": self.r,
            "dt": self.dt,
            "dx": self.dx,
            "analytical_wave_speed": float(analytical),
            "numerical_wave_speed": float(numerical_wave_speed),
            "speed_error_percent": float(speed_error),
            "numerical_wave_speed_mass": float(speed_mass),
            "speed_error_percent_mass": float(err_mass),
            "numerical_wave_speed_axis": float(speed_axis),
            "speed_error_percent_axis": float(err_axis),
            "front_positions_radial": [float(x) for x in front_positions],
            "front_positions_mass": [float(x) for x in front_positions_mass],
            "front_positions_axis": [float(x) for x in front_positions_axis],
            "times": [float(t) for t in times],
            "runtime_seconds": float(elapsed),
            "wave_speed_um_per_hr": float(numerical_wave_speed * 5),
        }

        return np.array(snapshots, dtype=np.float32), metrics


# ---------------------------------------------------------------------------
# Export
# ---------------------------------------------------------------------------
def export_results(field_history: np.ndarray, metrics: dict) -> None:
    print("[EXPORT] Saving Fisher-Kolmogorov results...")

    fk_path = OUTPUT_DIR / "fk_field_history.npy"
    metrics_path = OUTPUT_DIR / "fk_metrics.json"

    np.save(str(fk_path), field_history)

    with open(str(metrics_path), "w") as f:
        json.dump(metrics, f, indent=2)

    assert fk_path.exists(), f"FAILED TO WRITE: {fk_path}"
    assert metrics_path.exists(), f"FAILED TO WRITE: {metrics_path}"

    print(f"[SAVE] {fk_path}  shape={field_history.shape}  size={fk_path.stat().st_size}")
    print(f"[SAVE] {metrics_path}  size={metrics_path.stat().st_size}")
    print(
        f"  Wave speed: analytical={metrics['analytical_wave_speed']:.4f}, "
        f"numerical={metrics['numerical_wave_speed']:.4f}, "
        f"error={metrics['speed_error_percent']:.1f}%"
    )
    print(f"  Clinical velocity: {metrics['wave_speed_um_per_hr']:.2f} µm/hr")
    print(
        f"  Cross-check: mass-detector={metrics['numerical_wave_speed_mass']:.4f} "
        f"(err={metrics['speed_error_percent_mass']:.1f}%), "
        f"axis-detector={metrics['numerical_wave_speed_axis']:.4f} "
        f"(err={metrics['speed_error_percent_axis']:.1f}%)"
    )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    # Calibrated for c = 2*sqrt(D*r) = 4.0 px/unit-time = 20 µm/hr
    # at 5 µm/px, 1 hr/unit-time
    D = 10.0
    r = 0.4
    dt = 0.001       # Stable for explicit FD (CFL: dt ≤ dx²/(4D) = 0.025)
    n_steps = 20000  # t_total = 20.0, front travels ~80 px (asymptotic regime)
    save_interval = 2000

    print(f"[FK-FD] Calibrated params: D={D}, r={r}, dt={dt}")
    print(f"[FK-FD] Target wave speed: {2 * np.sqrt(D * r) * 5:.1f} µm/hr")

    solver = FDSolverFK(
        grid_size=(512, 512),
        D=D,
        r=r,
        dt=dt,
        n_steps=n_steps,
        save_interval=save_interval,
    )

    field_history, metrics = solver.run(
        n_steps=n_steps,
        save_interval=save_interval,
        initial_conditions="tumor_core",
    )

    export_results(field_history, metrics)

    print("\n[SUCCESS] Month 3 Week 2 Complete: Fisher-Kolmogorov PDE Solver (Finite Difference)")
    print(f"  Analytical wave speed: {metrics['analytical_wave_speed']:.4f} px/unit-time")
    print(f"  Numerical (radial):    {metrics['numerical_wave_speed']:.4f} px/unit-time (err={metrics['speed_error_percent']:.1f}%)")
    print(f"  Numerical (mass):      {metrics['numerical_wave_speed_mass']:.4f} px/unit-time (err={metrics['speed_error_percent_mass']:.1f}%)")
    print(f"  Numerical (axis):      {metrics['numerical_wave_speed_axis']:.4f} px/unit-time (err={metrics['speed_error_percent_axis']:.1f}%)")
    print(f"  Clinical velocity:     {metrics['wave_speed_um_per_hr']:.1f} µm/hr")


if __name__ == "__main__":
    main()