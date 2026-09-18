#!/usr/bin/env python3
"""
Phase 10: FNO Benchmark - Compare FNO vs Finite Difference PDE Solver
"""
import torch
import time
import numpy as np
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent))
from fno_solver import FNO3d


def benchmark_fno(
    model_path="output/fno_model.pth",
    grid_size=64,
    n_runs=10,
    device=None
):
    device = device or ("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[Benchmark] Device: {device}")

    # Load trained FNO
    checkpoint = torch.load(model_path, map_location=device)
    state_dict = checkpoint["model_state"] if "model_state" in checkpoint else checkpoint
    
    model = FNO3d(in_channels=2, out_channels=1, modes=8, width=20).to(device)
    model.load_state_dict(state_dict)
    model.eval()

    # Create test input
    batch_size = 1
    x = torch.randn(batch_size, 2, grid_size, grid_size, grid_size, device=device)

    # Warmup
    with torch.no_grad():
        for _ in range(3):
            _ = model(x)

    # Benchmark FNO
    torch.cuda.synchronize() if device == "cuda" else None
    start = time.perf_counter()
    for _ in range(n_runs):
        with torch.no_grad():
            _ = model(x)
    torch.cuda.synchronize() if device == "cuda" else None
    fno_time = (time.perf_counter() - start) / n_runs

    print(f"[Benchmark] FNO inference: {fno_time*1000:.2f} ms per forward pass ({grid_size}^3 grid)")

    # Estimate FD time (theoretical)
    # 3D FD with 5-point stencil: ~O(grid^3) operations per step
    # For 1500 steps on 64^3: ~1500 * 64^3 * 27 ops ≈ 3.2e9 ops
    # CPU ~10^9 ops/sec → ~3 seconds per trajectory
    # FNO: single forward pass ~milliseconds
    fd_estimate = 3.0  # seconds for 1500 steps on 64^3
    speedup = fd_estimate / (fno_time / 1000)
    
    print(f"[Benchmark] Estimated FD time (1500 steps): {fd_estimate:.1f} s")
    print(f"[Benchmark] FNO speedup: ~{speedup:.0f}x faster")
    print(f"[Benchmark] FNO enables real-time control: {fno_time*1000:.1f} ms per prediction")

    return {
        "fno_time_ms": fno_time * 1000,
        "fd_estimate_s": fd_estimate,
        "speedup": speedup
    }


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Benchmark FNO vs FD")
    parser.add_argument("--model", type=str, default="output/fno_model.pth")
    parser.add_argument("--grid", type=int, default=64)
    parser.add_argument("--runs", type=int, default=10)
    args = parser.parse_args()

    benchmark_fno(model_path=args.model, grid_size=args.grid, n_runs=args.runs)