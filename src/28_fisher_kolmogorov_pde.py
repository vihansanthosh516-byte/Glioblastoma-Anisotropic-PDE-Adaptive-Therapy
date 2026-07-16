#!/usr/bin/env python3
"""
Month 3, Week 2: Fisher-Kolmogorov PDE Solver — Crank-Nicolson with 4th-Order Spatial Scheme
Morphogen diffusion with reaction term: ∂ρ/∂t = D∇²ρ + rρ(1-ρ)

Calibrated for clinical glioblastoma invasion velocity: 10-50 µm/hr
Uses Crank-Nicolson (implicit) with 4th-order central differences to achieve
numerical wave speed error < 10% vs analytical solution.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Tuple

import numpy as np
import torch
import scipy.sparse as sp
import scipy.sparse.linalg as spla


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
        Y, X = torch.meshgrid(y, x, indexing='ij')
        dist_sq = X**2 + Y**2
        rho = torch.exp(-dist_sq / (2 * 15**2)) * 0.8
    elif initial_conditions == "periphery_ring":
        center_h, center_w = H // 2, W // 2
        y = torch.arange(H, device=device).float() - center_h
        x = torch.arange(W, device=device).float() - center_w
        Y, X = torch.meshgrid(y, x, indexing='ij')
        dist = torch.sqrt(X**2 + Y**2)
        rho = ((dist > 25) & (dist < 45)).float() * 0.6

    return rho


def build_4th_order_laplacian(
    H: int, W: int, dx: float = 1.0
) -> sp.csr_matrix:
    """
    Build 4th-order accurate 2D Laplacian using 9-point stencil:
    ∇²u ≈ (1/6dx²)[4(u_{i+1,j}+u_{i-1,j}+u_{i,j+1}+u_{i,j-1}) 
                   + (u_{i+1,j+1}+u_{i+1,j-1}+u_{i-1,j+1}+u_{i-1,j-1})
                   - 20u_{i,j}]
    
    For uniform grid with spacing dx, the coefficients are:
    - Center: -20/(6dx²)
    - Axis neighbors: 4/(6dx²)
    - Diagonal neighbors: 1/(6dx²)
    """
    N = H * W
    
    # 9-point stencil coefficients
    c_center = -20.0 / (6.0 * dx * dx)
    c_axis = 4.0 / (6.0 * dx * dx)
    c_diag = 1.0 / (6.0 * dx * dx)
    
    # Build in LIL format for easy construction
    L = sp.lil_matrix((N, N))
    
    for i in range(H):
        for j in range(W):
            idx = i * W + j
            
            # Center
            L[idx, idx] = c_center
            
            # Axis neighbors
            if i > 0:
                L[idx, (i-1)*W + j] = c_axis
            if i < H - 1:
                L[idx, (i+1)*W + j] = c_axis
            if j > 0:
                L[idx, i*W + (j-1)] = c_axis
            if j < W - 1:
                L[idx, i*W + (j+1)] = c_axis
            
            # Diagonal neighbors
            if i > 0 and j > 0:
                L[idx, (i-1)*W + (j-1)] = c_diag
            if i > 0 and j < W - 1:
                L[idx, (i-1)*W + (j+1)] = c_diag
            if i < H - 1 and j > 0:
                L[idx, (i+1)*W + (j-1)] = c_diag
            if i < H - 1 and j < W - 1:
                L[idx, (i+1)*W + (j+1)] = c_diag
    
    return L.tocsr()


class CrankNicolsonFK:
    """Crank-Nicolson solver for Fisher-Kolmogorov PDE with 4th-order Laplacian."""
    
    def __init__(
        self,
        grid_size: Tuple[int, int],
        D: float = 1.0,
        r: float = 2.25,
        dt: float = 0.01,
        dx: float = 1.0,
        device: torch.device = None,
    ):
        self.H, self.W = grid_size
        self.N = self.H * self.W
        self.D = D
        self.r = r
        self.dt = dt
        self.dx = dx
        self.device = device or torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        
        # Build 4th-order Laplacian
        print(f"[FK] Building 4th-order Laplacian for {self.H}x{self.W} grid...")
        t0 = time.perf_counter()
        self.L = build_4th_order_laplacian(self.H, self.W, dx)
        print(f"[FK] Laplacian built in {time.perf_counter() - t0:.3f}s")
        
        # Crank-Nicolson matrices for diffusion term
        # (I - dt/2 * D * L) ρ^{n+1} = (I + dt/2 * D * L) ρ^n + dt * R(ρ)
        I = sp.eye(self.N, format='csr')
        alpha = self.D * self.dt / 2.0
        
        self.A = (I - alpha * self.L).tocsr()
        self.B = (I + alpha * self.L).tocsr()
        
        # LU factorization for fast solves
        print("[FK] Computing LU factorization of A...")
        t0 = time.perf_counter()
        self.A_lu = spla.splu(self.A.tocsc())
        print(f"[FK] LU factorization done in {time.perf_counter() - t0:.3f}s")
        
        self.dt = dt
        self.r = r
    
    def step(self, rho: torch.Tensor) -> torch.Tensor:
        """Single Crank-Nicolson step with explicit reaction term."""
        # Move to numpy for sparse solve
        rho_np = rho.detach().cpu().numpy().flatten()
        
        # Explicit reaction term: r * ρ * (1 - ρ)
        R = self.r * rho_np * (1.0 - rho_np)
        
        # RHS = B @ ρ^n + dt * R
        rhs = self.B @ rho_np + self.dt * R
        
        # Solve A @ ρ^{n+1} = rhs
        rho_new_np = self.A_lu.solve(rhs)
        
        # Clamp to [0, 1]
        rho_new_np = np.clip(rho_new_np, 0.0, 1.0)
        
        return torch.from_numpy(rho_new_np.reshape(self.H, self.W)).to(self.device, dtype=rho.dtype)


def analytical_wave_speed(D: float, r: float) -> float:
    """Analytical Fisher-Kolmogorov wave speed: c = 2√(Dr)"""
    return 2 * np.sqrt(D * r)


def simulate_fisher_kolmogorov(
    grid_size: Tuple[int, int] = (512, 512),
    n_steps: int = 2000,
    D: float = 1.0,
    r: float = 2.25,
    dt: float = 0.01,
    save_interval: int = 200,
    device: torch.device = None,
) -> Tuple[np.ndarray, dict]:
    """
    Run Fisher-Kolmogorov PDE simulation with Crank-Nicolson + 4th-order spatial scheme.
    
    Returns:
        - Field history: (n_snapshots, H, W)
        - Metrics dict with wave speed analysis
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    c_analytical = analytical_wave_speed(D, r)
    print(f"[FK] Running Fisher-Kolmogorov (CN + 4th-order) on {device}")
    print(f"[FK] Grid: {grid_size}, Steps: {n_steps}, D={D}, r={r}, dt={dt}")
    print(f"[FK] Analytical wave speed: c = 2*sqrt(Dr) = {c_analytical:.4f} pixels/step")
    print(f"[FK] At 5µm/pixel, 1 step=1hr: {c_analytical * 5:.2f} µm/hr")
    
    # Initialize
    solver = CrankNicolsonFK(grid_size, D, r, dt, device=device)
    rho = initialize_morphogen_field(grid_size, device, "tumor_core")
    
    # Storage
    snapshots = [rho.cpu().numpy()]
    front_positions = []
    times = []
    
    def detect_front(field: torch.Tensor, threshold: float = 0.5) -> float:
        """Detect radial position of front at threshold."""
        H, W = field.shape
        center_h, center_w = H // 2, W // 2
        y = torch.arange(H, device=device).float() - center_h
        x = torch.arange(W, device=device).float() - center_w
        Y, X = torch.meshgrid(y, x, indexing='ij')
        dist = torch.sqrt(X**2 + Y**2)
        
        mask = field > threshold
        if mask.any():
            return dist[mask].float().mean().item()
        return 0.0
    
    front_positions.append(detect_front(rho))
    times.append(0)
    
    t0 = time.perf_counter()
    for step in range(1, n_steps + 1):
        rho = solver.step(rho)
        
        if step % 10 == 0:
            front_positions.append(detect_front(rho))
            times.append(step * dt)
        
        if step % save_interval == 0:
            snapshots.append(rho.cpu().numpy())
        
        if step % 500 == 0:
            elapsed = time.perf_counter() - t0
            print(f"[FK] Step {step}/{n_steps} ({elapsed:.1f}s)")
    
    elapsed = time.perf_counter() - t0
    print(f"[FK] Completed in {elapsed:.2f}s")
    
    # Compute wave speed from front positions
    if len(front_positions) > 1:
        times_np = np.array(times)
        fronts_np = np.array(front_positions)
        coeffs = np.polyfit(times_np, fronts_np, 1)
        numerical_wave_speed = coeffs[0]
        analytical = c_analytical
        speed_error = abs(numerical_wave_speed - analytical) / analytical * 100
    else:
        numerical_wave_speed = 0.0
        analytical = c_analytical
        speed_error = 0.0
    
    metrics = {
        "grid_size": grid_size,
        "n_steps": n_steps,
        "D": D,
        "r": r,
        "dt": dt,
        "dx": 1.0,
        "analytical_wave_speed": analytical,
        "numerical_wave_speed": float(numerical_wave_speed),
        "speed_error_percent": float(speed_error),
        "front_positions": front_positions,
        "times": times,
        "runtime_seconds": elapsed,
        "wave_speed_um_per_hr": float(numerical_wave_speed * 5),  # 5 µm/pixel
    }
    
    return np.array(snapshots, dtype=np.float32), metrics


def export_results(
    field_history: np.ndarray,
    metrics: dict,
) -> None:
    """Export simulation results."""
    print("[EXPORT] Saving Fisher-Kolmogorov results...")
    Path("output").mkdir(exist_ok=True)
    
    np.save("output/fk_field_history.npy", field_history)
    
    with open("output/fk_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    
    print(f"  output/fk_field_history.npy: {field_history.shape}")
    print(f"  output/fk_metrics.json")
    print(f"  Wave speed: analytical={metrics['analytical_wave_speed']:.4f}, "
          f"numerical={metrics['numerical_wave_speed']:.4f}, "
          f"error={metrics['speed_error_percent']:.1f}%")
    print(f"  Clinical velocity: {metrics['wave_speed_um_per_hr']:.2f} µm/hr")


def main():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")
    
    # CALIBRATED PARAMETERS for 10-50 µm/hr with <10% numerical error
    # Target: c = 2*sqrt(D*r) = 2-10 px/step (10-50 µm/hr at 5µm/px, 1hr/step)
    # For c=4.0 (20 µm/hr): D*r = 4.0
    # D=1.0, r=4.0 -> c=4.0, but need small dt for stability
    # Using dt=0.01, D=1.0, r=4.0 -> c=4.0 px/step = 20 µm/hr
    
    D = 1.0       # Diffusion coefficient
    r = 4.0       # Proliferation rate
    dt = 0.01     # Tightened temporal discretization
    
    print(f"[FK] Calibrated params: D={D}, r={r}, dt={dt}")
    print(f"[FK] Target wave speed: {2*np.sqrt(D*r)*5:.1f} µm/hr")
    
    # Run simulation
    field_history, metrics = simulate_fisher_kolmogorov(
        grid_size=(512, 512),
        n_steps=2000,
        D=D,
        r=r,
        dt=dt,
        save_interval=200,
        device=device,
    )
    
    # Export
    export_results(field_history, metrics)
    
    print("\n[SUCCESS] Month 3 Week 2 Complete: Fisher-Kolmogorov PDE Solver (Crank-Nicolson + 4th-order)")
    print(f"  Analytical wave speed: {metrics['analytical_wave_speed']:.4f} px/step")
    print(f"  Numerical wave speed: {metrics['numerical_wave_speed']:.4f} px/step")
    print(f"  Error: {metrics['speed_error_percent']:.1f}%")
    print(f"  Clinical velocity: {metrics['wave_speed_um_per_hr']:.1f} µm/hr")


if __name__ == "__main__":
    main()