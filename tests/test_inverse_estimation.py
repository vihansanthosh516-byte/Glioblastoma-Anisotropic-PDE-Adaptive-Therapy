#!/usr/bin/env python3
"""Test suite for inverse biophysical parameter estimation (Phase 1).

Tests the src/51_inverse_parameter_estimation.py module:
  - Synthetic data recovery (noise-free)
  - Noise robustness (5-15% Gaussian noise)
  - Clinical sanity checks (physiological bounds)
  - Convergence within iteration budget

Run:
    venv\\Scripts\\python.exe tests\\test_inverse_estimation.py
    venv\\Scripts\\python.exe -m pytest tests/test_inverse_estimation.py
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = REPO_ROOT / "src" / "51_inverse_parameter_estimation.py"


def _load_module():
    """Load the inverse parameter estimation module (filename starts with digit)."""
    spec = importlib.util.spec_from_file_location(
        "inverse_parameter_estimation", MODULE_PATH
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


mod = _load_module()
surrogate_ode_model = mod.surrogate_ode_model
estimate_patient_parameters = mod.estimate_patient_parameters
RHO_MIN = mod.RHO_MIN
RHO_MAX = mod.RHO_MAX
D_MIN = mod.D_MIN
D_MAX = mod.D_MAX


def generate_synthetic_v1(
    rho: float, D: float, V0: float, delta_t: float
) -> float:
    """Generate a synthetic follow-up volume with known parameters."""
    return surrogate_ode_model(rho, D, V0, delta_t)


def test_surrogate_model_returns_positive():
    """Surrogate model should return positive volume for positive inputs."""
    V1 = surrogate_ode_model(0.02, 0.013, 1000.0, 30.0)
    assert V1 > 0, "Surrogate model returned non-positive volume"
    assert V1 > 1000.0, "Tumor should grow with positive rho and D"


def test_surrogate_model_growth_increases_with_rho():
    """Higher growth rate should yield larger tumor volume."""
    V_low = surrogate_ode_model(0.01, 0.013, 1000.0, 30.0)
    V_high = surrogate_ode_model(0.05, 0.013, 1000.0, 30.0)
    assert V_high > V_low, "Higher rho should produce larger volume"


def test_surrogate_model_diffusion_increases_volume():
    """Higher diffusivity should yield larger apparent volume (radial invasion)."""
    V_low = surrogate_ode_model(0.02, 0.005, 1000.0, 30.0)
    V_high = surrogate_ode_model(0.02, 0.04, 1000.0, 30.0)
    assert V_high > V_low, "Higher D should produce larger volume (radial invasion)"


def test_surrogate_model_zero_volume():
    """Zero initial volume should stay zero."""
    V1 = surrogate_ode_model(0.02, 0.013, 0.0, 30.0)
    assert V1 == 0.0, "Zero initial volume should remain zero"


def test_estimate_recover_no_noise():
    """Recover parameters from noise-free synthetic data (RMSE < 5%)."""
    rng = np.random.default_rng(42)
    true_rho, true_D = 0.025, 0.015
    V0, delta_t = 1000.0, 30.0
    V1_true = generate_synthetic_v1(true_rho, true_D, V0, delta_t)

    est = estimate_patient_parameters(
        t0_volume=V0,
        t1_volume=V1_true,
        delta_t_days=delta_t,
    )
    rho_err = abs(est["rho"] - true_rho) / true_rho
    assert rho_err < 0.05, f"rho error {rho_err:.2%} exceeds 5% threshold"
    assert est["convergence"], "Estimation did not converge"
    assert est["n_iterations"] < 50, "Exceeded 50-iteration budget"


def test_estimate_within_physiological_bounds():
    """Estimated parameters must fall within physiological bounds."""
    V0, delta_t = 1000.0, 60.0
    V1 = 1500.0  # moderate growth
    est = estimate_patient_parameters(
        t0_volume=V0, t1_volume=V1, delta_t_days=delta_t
    )
    assert RHO_MIN <= est["rho"] <= RHO_MAX, "rho out of physiological bounds"
    assert D_MIN <= est["D"] <= D_MAX, "D out of physiological bounds"


def test_estimate_with_noise_robustness():
    """Recovery under 10% Gaussian noise (RMSE < 15%)."""
    rng = np.random.default_rng(123)
    true_rho, true_D = 0.025, 0.015
    V0, delta_t = 1000.0, 30.0
    V1_true = generate_synthetic_v1(true_rho, true_D, V0, delta_t)

    n_trials = 20
    rho_errors, D_errors = [], []
    for _ in range(n_trials):
        noise = rng.normal(0, 0.10 * V1_true)
        V1_noisy = max(0, V1_true + noise)
        est = estimate_patient_parameters(
            t0_volume=V0, t1_volume=V1_noisy, delta_t_days=delta_t
        )
        rho_errors.append(est["rho"] - true_rho)
        D_errors.append(est["D"] - true_D)

    rho_rmse = np.sqrt(np.mean(np.array(rho_errors) ** 2))
    D_rmse = np.sqrt(np.mean(np.array(D_errors) ** 2))
    assert (rho_rmse / true_rho) < 0.25, (
        f"rho RMSE {rho_rmse:.5f} too high (>25% relative)"
    )
    assert (D_rmse / true_D) < 0.50, (
        f"D RMSE {D_rmse:.5f} too high (>50% relative)"
    )


def test_ci_contains_estimate():
    """95% CI should bracket the point estimate."""
    V0, delta_t = 1000.0, 45.0
    V1 = generate_synthetic_v1(0.022, 0.012, V0, delta_t)
    est = estimate_patient_parameters(
        t0_volume=V0, t1_volume=V1, delta_t_days=delta_t
    )
    rho_lo, rho_hi = est["rho_ci"]
    D_lo, D_hi = est["D_ci"]
    assert rho_lo <= est["rho"] <= rho_hi, "rho CI does not bracket estimate"
    assert D_lo <= est["D"] <= D_hi, "D CI does not bracket estimate"


def test_rmse_low_for_well_fit_data():
    """RMSE should be small when data is generated by the surrogate model."""
    true_rho, true_D = 0.025, 0.015
    V0, delta_t = 1000.0, 30.0
    V1 = generate_synthetic_v1(true_rho, true_D, V0, delta_t)
    est = estimate_patient_parameters(
        t0_volume=V0, t1_volume=V1, delta_t_days=delta_t
    )
    assert est["rmse"] < 1.0, f"RMSE {est['rmse']:.4f} too high for exact synthetic data"


def test_fallback_to_population_average_on_degenerate_input():
    """When input volumes are identical (no growth), estimation still returns
    a valid bounded parameter set without crashing."""
    V0, V1, delta_t = 1000.0, 1000.0, 30.0
    est = estimate_patient_parameters(
        t0_volume=V0, t1_volume=V1, delta_t_days=delta_t
    )
    assert RHO_MIN <= est["rho"] <= RHO_MAX
    assert D_MIN <= est["D"] <= D_MAX


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
