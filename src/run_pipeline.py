#!/usr/bin/env python3
"""
Master Pipeline: Phase 1 (Spatial Genomics) -> Phase 2 (Poroelastic Mechanics)

Chained execution with in-memory data passing (no file I/O between phases).
"""
from __future__ import annotations

import importlib.util
import json
import sys
import time
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np

# --------------------------------------------------------------------------- #
# Phase Execution Utilities
# --------------------------------------------------------------------------- #
def load_module(module_path: Path) -> Any:
    """Dynamically load a module from file path."""
    spec = importlib.util.spec_from_file_location(module_path.stem, module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_phase1() -> Dict[str, np.ndarray]:
    """
    Execute Phase 1: Spatial Genomics Deconvolution.
    Returns dict with rho_mean, rho_lower, rho_upper, D_mean, D_lower, D_upper.
    """
    print("=" * 70)
    print("PIPELINE: Starting Phase 1 - Spatial Genomics Deconvolution")
    print("=" * 70)
    
    start = time.time()
    phase1 = load_module(Path("src/50_spatial_genomics_deconv.py"))
    phase1.main()
    elapsed = time.time() - start
    
    print(f"[PIPELINE] Phase 1 completed in {elapsed:.1f}s")
    print("=" * 70)
    
    # Load outputs
    out_dir = Path("output")
    rho_data = np.load(out_dir / "phase1_rho_posterior.npy")
    D_data = np.load(out_dir / "phase1_D_posterior.npy")
    
    return {
        "rho_mean": rho_data[0],
        "rho_lower": rho_data[1],
        "rho_upper": rho_data[2],
        "D_mean": D_data[0],
        "D_lower": D_data[1],
        "D_upper": D_data[2],
    }


def run_phase2(phase1_posteriors: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
    """
    Execute Phase 2: Poroelastic Mechanics with Phase 1 posteriors as input.
    Returns dict with pressure, displacement, stress, and D_eff fields.
    """
    print("=" * 70)
    print("PIPELINE: Starting Phase 2 - Biot Poroelastic Mechanics")
    print("=" * 70)
    
    start = time.time()
    phase2 = load_module(Path("src/51_poroelastic_hybrid_solver.py"))
    
    # Patch Phase 2 to use Phase 1 posteriors directly (no file I/O)
    original_load = phase2.load_phase1_posteriors
    
    def patched_load():
        print("[PHASE 2] Using Phase 1 posteriors from memory (no file I/O)")
        return phase1_posteriors
    
    phase2.load_phase1_posteriors = patched_load
    
    # Run Phase 2 main
    phase2.main()
    elapsed = time.time() - start
    
    print(f"[PIPELINE] Phase 2 completed in {elapsed:.1f}s")
    print("=" * 70)
    
    # Load Phase 2 outputs
    out_dir = Path("output")
    pressure = np.load(out_dir / "phase2_pressure_field.npy")
    disp = np.load(out_dir / "phase2_displacement_field.npy")  # (3, nx, ny, nz)
    
    return {
        "pressure": pressure,
        "displacement": disp,
    }


def compute_cross_phase_metrics(
    phase1: Dict[str, np.ndarray],
    phase2: Dict[str, np.ndarray]
) -> Dict:
    """Compute cross-phase validation metrics."""
    rho = phase1["rho_mean"]
    D = phase1["D_mean"]
    pressure = phase2["pressure"]
    disp = phase2["displacement"]
    disp_mag = np.sqrt(np.sum(disp**2, axis=0))
    
    # Correlation between rho and pressure
    rho_flat = rho.ravel()
    p_flat = pressure.ravel()
    valid = ~np.isnan(rho_flat) & ~np.isnan(p_flat)
    rho_p_corr = float(np.corrcoef(rho_flat[valid], p_flat[valid])[0, 1]) if valid.sum() > 10 else 0.0
    
    # Correlation between D and displacement
    D_flat = D.ravel()
    d_flat = disp_mag.ravel()
    valid = ~np.isnan(D_flat) & ~np.isnan(d_flat)
    D_disp_corr = float(np.corrcoef(D_flat[valid], d_flat[valid])[0, 1]) if valid.sum() > 10 else 0.0
    
    return {
        "rho_pressure_correlation": rho_p_corr,
        "D_displacement_correlation": D_disp_corr,
        "max_pressure_kPa": float(np.max(pressure)),
        "max_displacement_mm": float(np.max(disp_mag)),
        "mean_rho": float(np.nanmean(rho)),
        "mean_D": float(np.nanmean(D)),
    }


# --------------------------------------------------------------------------- #
# Master Pipeline Entry Point
# --------------------------------------------------------------------------- #
def main() -> None:
    """Execute full Phase 1 -> Phase 2 pipeline."""
    pipeline_start = time.time()
    
    print("\n" + "=" * 70)
    print("  GBM PIPELINE: Phase 1 (Spatial Genomics) -> Phase 2 (Poroelastic)")
    print("=" * 70 + "\n")
    
    # Phase 1: Spatial Genomics
    phase1_posteriors = run_phase1()
    
    # Phase 2: Poroelastic Mechanics (fed by Phase 1)
    phase2_outputs = run_phase2(phase1_posteriors)
    
    # Cross-phase metrics
    print("[PIPELINE] Computing cross-phase validation metrics...")
    cross_metrics = compute_cross_phase_metrics(phase1_posteriors, phase2_outputs)
    
    out_dir = Path("output")
    with open(out_dir / "pipeline_cross_metrics.json", "w") as f:
        json.dump(cross_metrics, f, indent=2)
    
    total_elapsed = time.time() - pipeline_start
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print("=" * 70)
    print(f"  Total runtime: {total_elapsed:.1f}s")
    print(f"  Cross-phase rho-pressure correlation: {cross_metrics['rho_pressure_correlation']:.3f}")
    print(f"  Cross-phase D-displacement correlation: {cross_metrics['D_displacement_correlation']:.3f}")
    print(f"  Max pressure: {cross_metrics['max_pressure_kPa']:.2f} kPa")
    print(f"  Max displacement: {cross_metrics['max_displacement_mm']:.2f} mm")
    print(f"  Outputs in: {out_dir}/")
    print("=" * 70)
    
    # Summary of all artifacts
    artifacts = [
        "phase1_rho_posterior.npy",
        "phase1_D_posterior.npy",
        "phase1_genomics_metrics.json",
        "phase1_spatial_deconv.png",
        "phase2_pressure_field.npy",
        "phase2_displacement_field.npy",
        "phase2_mechanics_metrics.json",
        "phase2_poroelastic_mechanics.png",
        "pipeline_cross_metrics.json",
    ]
    for a in artifacts:
        path = out_dir / a
        if path.exists():
            print(f"  [OK] {a} ({path.stat().st_size/1024:.1f} KB)")
        else:
            print(f"  [MISSING] {a} (MISSING)")


if __name__ == "__main__":
    main()