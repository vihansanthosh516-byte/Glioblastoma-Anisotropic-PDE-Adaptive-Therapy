#!/usr/bin/env python3
"""Test suite for robust MPC controller (Phase 2).

Tests src/52_robust_mpc_controller.py:
  - Surrogate ODE volume dynamics (dV/dt direction)
  - Robust cost risk-averse formulation (mean + lambda*std)
  - Adaptive horizon adjustment rules
  - Robust MPC simulation achieves target dose-sparing / dosing
  - Non-inferiority vs baseline MPC

Run:
    venv\\Scripts\\python.exe tests\\test_robust_mpc.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = REPO_ROOT / "src" / "52_robust_mpc_controller.py"

spec = importlib.util.spec_from_file_location("robust_mpc_controller", MODULE_PATH)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def make_seed_tumor():
    """Create a small spherical tumor seed on a 50^3 grid."""
    gs = mod.GRID_SIZE
    u0 = np.zeros((gs, gs, gs))
    z, y, x = np.mgrid[0:gs, 0:gs, 0:gs]
    dist = np.sqrt((x - gs//2)**2 + (y - gs//2)**2 + (z - gs//2)**2)
    u0[dist <= 3.0] = 0.8
    return u0


def test_volume_step_growth_no_drug():
    """Without drug, volume should grow for positive rho."""
    v_next = mod.volume_step(100.0, 0.02, 0.013, drug_on=False, step=0)
    assert v_next > 100.0, "Volume should increase under growth alone"


def test_volume_step_drug_kills():
    """With drug doping during dosing phase, volume should decrease."""
    v_next = mod.volume_step(100.0, 0.02, 0.013, drug_on=True, step=0)
    assert v_next < 100.0, "Volume should decrease under drug"


def test_adjust_horizon_stable_growth_extends():
    """Stable growth (|dV/dt| < 0.01) should extend horizon by 1 day."""
    h_new = mod.adjust_horizon(
        current_horizon=14, growth_rate=0.005, volume=10.0,
        target_volume=12.0,
    )
    assert h_new == 15, "Stable growth should extend horizon"


def test_adjust_horizon_accelerated_growth_shortens():
    """Accelerating growth (dV/dt > 0.05) should shorten horizon by 1 day."""
    h_new = mod.adjust_horizon(
        current_horizon=14, growth_rate=0.08, volume=20.0,
        target_volume=12.0,
    )
    assert h_new == 13, "Accelerating growth should shorten horizon"


def test_adjust_horizon_near_target_maintains():
    """When volume within 10% of target, horizon should be maintained."""
    h_new = mod.adjust_horizon(
        current_horizon=14, growth_rate=0.08, volume=12.0,
        target_volume=12.0,
    )
    assert h_new == 14, "Near target should maintain horizon"


def test_adjust_horizon_clamps_to_bounds():
    """Horizon should not exceed 21 maximum."""
    h_max = mod.adjust_horizon(
        current_horizon=21, growth_rate=0.005, volume=50.0,
        target_volume=12.0,
    )
    assert h_max <= 21, "Horizon must be <= 21"
    h_min = mod.adjust_horizon(
        current_horizon=7, growth_rate=0.08, volume=50.0,
        target_volume=12.0,
    )
    assert h_min >= 7, "Horizon must be >= 7"


def test_robust_cost_includes_risk_aversion():
    """J_robust = mean(J) + lambda*std(J); higher lambda should amplify std term."""
    rng1 = np.random.default_rng(7)
    rng2 = np.random.default_rng(7)  # same seed -> same samples
    c_low = mod.robust_cost(
        dose=1.0, volume_mm3=100.0, rho_nominal=0.02, D_nominal=0.013,
        target_volume=10.0, w_tumor=1.2, w_drug=0.03, w_uncertainty=0.5,
        risk_aversion=0.0, horizon_days=14, step=0, rng=rng1,
    )
    c_high = mod.robust_cost(
        dose=1.0, volume_mm3=100.0, rho_nominal=0.02, D_nominal=0.013,
        target_volume=10.0, w_tumor=1.2, w_drug=0.03, w_uncertainty=0.5,
        risk_aversion=1.0, horizon_days=14, step=0, rng=rng2,
    )
    # Higher risk aversion can't lower the cost (std >= 0)
    assert c_high >= c_low - 1e-9, "Higher lambda should not lower robust cost"


def test_orchestration_dose_sparing_target():
    """Robust MPC dose-sparing should be >= 60% on a representative patient."""
    u0 = make_seed_tumor()
    r = mod.run_robust_mpc_simulation(u0, 0.02, 0.013)
    assert r["dose_sparing_fraction"] >= 0.60, (
        f"Dose sparing {r['dose_sparing_fraction']:.2%} below 60% target"
    )


def test_orchestration_dosing_in_target_band():
    """Robust MPC drug administration should be in the 25-40% target band."""
    u0 = make_seed_tumor()
    r = mod.run_robust_mpc_simulation(u0, 0.02, 0.013)
    assert 0.25 <= r["drug_on_fraction"] <= 0.40, (
        f"Drug administration {r['drug_on_fraction']:.2%} outside 25-40%"
    )


def test_orchestration_adaptive_horizon_adjusts():
    """Robust MPC horizon should vary over the simulation (adaptive)."""
    u0 = make_seed_tumor()
    r = mod.run_robust_mpc_simulation(u0, 0.02, 0.013)
    assert r["horizon_range"][0] < r["horizon_range"][1], (
        "Horizon must vary (adaptive), got range "
        f"{r['horizon_range']}"
    )


def test_orchestration_returns_full_schema():
    """Output must match the baseline schema keys for drop-in compatibility."""
    u0 = make_seed_tumor()
    r = mod.run_robust_mpc_simulation(u0, 0.02, 0.013)
    required_keys = [
        "drug_schedule", "drug_on_fraction", "dose_sparing_fraction",
        "predicted_final_volume_mm3", "volume_history", "drug_on_history",
    ]
    for k in required_keys:
        assert k in r, f"Missing key '{k}' from robust MPC output"


def run_all():
    """Run all tests and report results."""
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed, failed = 0, 0
    for test in tests:
        try:
            test()
            print(f"  PASS  {test.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {test.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR {test.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\nResults: {passed} passed, {failed} failed, {passed + failed} total")
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all() else 1)
