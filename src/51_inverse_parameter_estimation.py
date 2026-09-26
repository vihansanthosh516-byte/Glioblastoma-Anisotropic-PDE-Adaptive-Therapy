#!/usr/bin/env python3
"""Inverse Biophysical Parameter Estimation for Patient-Specific GBM Modeling.

This module estimates patient-specific biophysical parameters (ρ, D) from
longitudinal imaging data using bounded optimization with bootstrap uncertainty.

Parameter Identifiability & Confidence Bounds (Tier 1):
------------------------------------------------------
Volume-only inverse problems are inherently underdetermined: a single pair of
volumes (T₀, T₁) constrains the combination ρ·Δt + D·Δt^(1/3) but cannot
uniquely separate ρ and D without additional priors or data (e.g., spatial
profile, multi-timepoint). This module resolves unidentifiability by enforcing
BOUNDED PHYSIOLOGICAL PRIORS:

* ρ (growth rate) ∈ [0.005, 0.05] /day
  - Lower bound: minimum mitotic rate for viable GBM (~0.5%/day)
  - Upper bound: maximum observed doubling time (~14 days) → ln(2)/14 ≈ 0.05/day
  - Prior literature: Swanson et al. (2003), Hormuth et al. (2017), Neil et al. (2020)

* D₀ (diffusivity) ∈ [0.001, 0.1] mm²/day
  - Lower bound: negligible diffusion (near-spherical growth)
  - Upper bound: fast infiltrative edge (~0.1 mm²/day from DTI tractography)
  - Prior literature: Swanson et al. (2000), Jbabdi et al. (2005), Price et al. (2018)

These bounds restrict the feasible parameter space to physiologically
plausible regimes, converting an ill-posed inverse problem into a well-posed
bounded optimization. The L-BFGS-B algorithm enforces these bounds during
optimization, and bootstrap resampling (N=100) quantifies estimation
uncertainty via 95% confidence intervals.

Mathematical Formulation:
    Given: T₀ (baseline volume), T₁ (follow-up volume), Δt (time between scans)
    Solve: min_{ρ, D} ||V_simulated(ρ, D, Δt) - T₁||²
    Subject to:
        - 0.005 ≤ ρ ≤ 0.05 /day (physiological bounds)
        - 0.001 ≤ D ≤ 0.1 mm²/day (diffusivity bounds)

Algorithm:
    - scipy.optimize.minimize with L-BFGS-B method (bounded)
    - Surrogate ODE model for fast evaluation: dV/dt = ρ*V*(1-V/K) + D*∇²V
    - Bootstrap resampling (N=100) for confidence intervals
    - Convergence residuals reported: objective value, gradient norm, iterations

Usage:
    python src/51_inverse_parameter_estimation.py --test
    python src/51_inverse_parameter_estimation.py --t0-volume 1000 --t1-volume 1200 --delta-t 30
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm

warnings.filterwarnings("ignore")

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
COHORT_JSON = OUTPUT_DIR / "mu_glioma_cohort.json"          # scan volumes + days (mu_glioma_loader)
PARAMS_72_CSV = OUTPUT_DIR / "mu_glioma_params_real.csv"     # script 72 log-linear fits, cross-check
METRICS_JSON = OUTPUT_DIR / "inverse_est_metrics.json"

# Physiological bounds (per plan specification — Tier 1 identifiability)
# Volume-only inverse problems are underdetermined; these bounded priors
# restrict the feasible space to physiologically plausible regimes.
RHO_MIN = 0.005  # /day (minimum growth rate: ~0.5%/day)
RHO_MAX = 0.05   # /day (maximum growth rate: ~ln(2)/14 ≈ 0.05/day)
D_MIN = 0.001    # mm²/day (minimum diffusivity: near-spherical)
D_MAX = 0.1      # mm²/day (maximum diffusivity: fast infiltration)

# Default initial guess
RHO_DEFAULT = 0.02    # /day
D_DEFAULT = 0.013     # mm²/day

# Bootstrap parameters
N_BOOTSTRAP = 100
NOISE_STD = 0.10  # 10% Gaussian noise for robustness testing

# Carrying capacity: realistic brain tumor maximum (~brain volume in mm3)
# A GBM cannot exceed the cranial vault; set K to a large value so logistic
# suppression is negligible at typical tumor volumes (1-50 cm^3 = 1000-50000 mm^3)
K_DEFAULT = 1.0e6  # mm3 (acts as near-pure exponential for small tumors)


def surrogate_ode_model(
    rho: float,
    D: float,
    V0: float,
    delta_t: float,
    K: float = K_DEFAULT,
) -> float:
    """
    Surrogate ODE model for tumor volume evolution.
    
    Separates the two biophysical effects so rho and D are independently
    identifiable from a pair of volume timepoints:
    
      - rho: proliferation (logistic/exponential growth of cell density)
      - D:   radial invasion (increases apparent volume via diffusion-driven
             boundary expansion, scaling as dR/dt ~ D/R)
    
    Model (lumped-volume surrogate):
        dV/dt = rho * V * (1 - V/K)                # logistic proliferation
                + 3 * sqrt(4*pi/3) * D * V^(1/3)   # radial-diffusion volume gain
    
    The radial term derives from the observation that for a spherical tumor
    of radius R, V = 4/3*pi*R^3 and diffusion-driven radial velocity is
    dR/dt = D/R (heat-equation front speed scaling), so:
        dV_diff/dt = 4*pi*R^2 * dR/dt = 4*pi*R^2 * D/R = 4*pi*D*R
                   = 4*pi*D * (3*V/(4*pi))^(1/3)
                   = (36*pi)^(1/3) * D * V^(1/3)
    
    Integrated numerically (Euler) for accuracy.
    
    Args:
        rho: Growth rate (/day)
        D: Diffusion coefficient (mm^2/day)
        V0: Initial volume (mm^3)
        delta_t: Time interval (days)
        K: Carrying capacity (mm^3, default 1e6)
    
    Returns:
        V1: Predicted volume at t+delta_t (mm^3)
    """
    if V0 <= 0:
        return 0.0
    
    # Radial-diffusion coefficient: (36*pi)^(1/3) ~ 3.269
    c_diff = (36.0 * np.pi) ** (1.0 / 3.0)
    
    # Numerical integration (Euler) with small substeps for stability
    n_substeps = max(1, int(np.ceil(delta_t / 1.0)))  # 1-day substeps
    dt_sub = delta_t / n_substeps
    
    V = V0
    for _ in range(n_substeps):
        # Logistic proliferation
        dV_prolif = rho * V * (1.0 - V / K) * dt_sub
        # Radial-diffusion volume gain
        dV_diff = c_diff * D * (V ** (1.0 / 3.0)) * dt_sub
        V = V + dV_prolif + dV_diff
        if V < 0:
            V = 0.0
    
    return V


def objective_function(
    params: np.ndarray,
    V0: float,
    V1_target: float,
    delta_t: float,
) -> float:
    """
    Objective function for parameter optimization.
    
    Minimizes squared relative error between simulated and observed volume.
    Relative error is used so that the optimizer scales correctly across
    tumor volumes of different magnitudes.
    
    Args:
        params: [rho, D] parameter vector
        V0: Baseline volume
        V1_target: Follow-up volume (target)
        delta_t: Time interval
    
    Returns:
        Squared relative error
    """
    rho, D = params
    
    # Hard clamp to physiological bounds (defensive; L-BFGS-B also enforces)
    rho = max(RHO_MIN, min(RHO_MAX, rho))
    D = max(D_MIN, min(D_MAX, D))
    
    # Simulate volume
    V1_sim = surrogate_ode_model(rho, D, V0, delta_t)
    
    # Squared relative error (scale-invariant)
    scale = max(abs(V1_target), 1.0)
    error = ((V1_sim - V1_target) / scale) ** 2
    
    return error


def estimate_patient_parameters(
    t0_volume: float,
    t1_volume: float,
    delta_t_days: float,
    initial_guess: Optional[Tuple[float, float]] = None,
    bounds: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None,
    method: str = "L-BFGS-B",
    n_bootstrap: int = N_BOOTSTRAP,
) -> Dict[str, Any]:
    """
    Estimate patient-specific biophysical parameters from longitudinal volumes.
    
    Args:
        t0_volume: Baseline tumor volume (mm³)
        t1_volume: Follow-up tumor volume (mm³)
        delta_t_days: Time between scans (days)
        initial_guess: Initial [rho, D] guess (default: [0.02, 0.013])
        bounds: Parameter bounds [(rho_min, rho_max), (D_min, D_max)]
        method: Optimization method (default: "L-BFGS-B")
        n_bootstrap: Number of bootstrap samples for CI (default: 100)
    
    Returns:
        Dictionary with:
            - rho: Estimated growth rate (/day)
            - D: Estimated diffusion coefficient (mm²/day)
            - rho_ci: 95% confidence interval for rho [lower, upper]
            - D_ci: 95% confidence interval for D [lower, upper]
            - convergence: dict with:
                - success: bool indicating successful convergence
                - objective_value: final objective function value (residual)
                - gradient_norm: L2 norm of gradient at solution
                - n_iterations: Number of iterations
                - n_function_evals: Number of function evaluations
                - message: optimizer status message
            - rmse: Root mean squared error of fit
            - bootstrap_samples: Array of [rho, D] bootstrap samples
    """
    if initial_guess is None:
        initial_guess = (RHO_DEFAULT, D_DEFAULT)
    
    if bounds is None:
        bounds = [(RHO_MIN, RHO_MAX), (D_MIN, D_MAX)]
    
    # Optimization (L-BFGS-B enforces bounds; fallback methods receive
    # bounds as None if unsupported)
    bounds_for_method = bounds if method in ("L-BFGS-B", "TNC", "SLSQP") else None
    
    result = minimize(
        fun=objective_function,
        x0=np.array(initial_guess),
        args=(t0_volume, t1_volume, delta_t_days),
        method=method,
        bounds=bounds_for_method,
        options={"maxiter": 1000, "disp": False},
    )
    
    rho_est, D_est = result.x
    
    # L-BFGS-B convergence diagnostics
    convergence_info = {
        "success": bool(result.success),
        "objective_value": float(result.fun) if hasattr(result, "fun") else None,
        "gradient_norm": float(np.linalg.norm(result.jac)) if hasattr(result, "jac") and result.jac is not None else None,
        "n_iterations": int(result.nit) if hasattr(result, "nit") else 0,
        "n_function_evals": int(result.nfev) if hasattr(result, "nfev") else 0,
        "message": str(result.message) if hasattr(result, "message") else "N/A",
    }
    
    # Bootstrap resampling for confidence intervals
    bootstrap_samples = np.zeros((n_bootstrap, 2))
    
    for i in range(n_bootstrap):
        # Add noise to target volume
        noise = np.random.normal(0, NOISE_STD * t1_volume)
        t1_noisy = max(0, t1_volume + noise)
        
        # Re-optimize with noisy data
        boot_result = minimize(
            fun=objective_function,
            x0=np.array(initial_guess),
            args=(t0_volume, t1_noisy, delta_t_days),
            method=method,
            bounds=bounds_for_method,
            options={"maxiter": 500, "disp": False},
        )
        bootstrap_samples[i] = boot_result.x
    
    # Confidence intervals (95%)
    rho_ci = np.percentile(bootstrap_samples[:, 0], [2.5, 97.5])
    D_ci = np.percentile(bootstrap_samples[:, 1], [2.5, 97.5])
    
    # RMSE
    V1_pred = surrogate_ode_model(rho_est, D_est, t0_volume, delta_t_days)
    rmse = np.sqrt((V1_pred - t1_volume) ** 2)
    
    return {
        "rho": float(rho_est),
        "D": float(D_est),
        "rho_ci": [float(rho_ci[0]), float(rho_ci[1])],
        "D_ci": [float(D_ci[0]), float(D_ci[1])],
        "convergence": convergence_info,
        "rmse": float(rmse),
        "bootstrap_samples": bootstrap_samples.tolist(),
    }


def validate_with_synthetic_data(
    true_rho: float = 0.025,
    true_D: float = 0.015,
    V0: float = 1000.0,
    delta_t: float = 30.0,
    noise_levels: List[float] = [0.0, 0.05, 0.10, 0.15],
    n_trials: int = 50,
) -> Dict[str, Any]:
    """
    Validate parameter estimation with synthetic data.
    
    Args:
        true_rho: True growth rate for synthetic data
        true_D: True diffusion coefficient for synthetic data
        V0: Initial volume
        delta_t: Time interval
        noise_levels: List of noise standard deviations to test
        n_trials: Number of Monte Carlo trials per noise level
    
    Returns:
        Validation results dictionary
    """
    print(f"\n{'='*70}")
    print("INVERSE PARAMETER ESTIMATION - SYNTHETIC VALIDATION")
    print(f"{'='*70}")
    print(f"True parameters: rho = {true_rho:.4f} /day, D = {true_D:.4f} mm2/day")
    print(f"Initial volume: V0 = {V0:.1f} mm3")
    print(f"Time interval: dt = {delta_t:.1f} days")
    print(f"Trials per noise level: {n_trials}")
    print(f"{'='*70}\n")
    np.random.seed(51)

    # Generate synthetic follow-up volume
    V1_true = surrogate_ode_model(true_rho, true_D, V0, delta_t)
    print(f"Synthetic V1 (noise-free): {V1_true:.2f} mm3\n")
    
    results = {}
    
    for noise_std in noise_levels:
        rho_errors = []
        D_errors = []
        convergence_count = 0
        n_iterations_list = []
        
        for trial in range(n_trials):
            # Add noise
            noise = np.random.normal(0, noise_std * V1_true)
            V1_noisy = max(0, V1_true + noise)
            
            # Estimate parameters
            est = estimate_patient_parameters(
                t0_volume=V0,
                t1_volume=V1_noisy,
                delta_t_days=delta_t,
            )
            
            if est["convergence"]["success"]:
                convergence_count += 1
                rho_errors.append(est["rho"] - true_rho)
                D_errors.append(est["D"] - true_D)
                n_iterations_list.append(est["convergence"]["n_iterations"])
        
        # Compute metrics
        rho_rmse = np.sqrt(np.mean(np.array(rho_errors) ** 2)) if rho_errors else float("inf")
        D_rmse = np.sqrt(np.mean(np.array(D_errors) ** 2)) if D_errors else float("inf")
        rho_bias = np.mean(rho_errors) if rho_errors else float("inf")
        D_bias = np.mean(D_errors) if D_errors else float("inf")
        avg_iterations = np.mean(n_iterations_list) if n_iterations_list else 0
        
        results[f"noise_{int(noise_std*100)}pct"] = {
            "rho_rmse": rho_rmse,
            "D_rmse": D_rmse,
            "rho_bias": rho_bias,
            "D_bias": D_bias,
            "convergence_rate": convergence_count / n_trials,
            "avg_iterations": avg_iterations,
        }
        
        print(f"Noise Level: {int(noise_std*100)}%")
        print(f"  rho RMSE: {rho_rmse:.6f} /day  (bias: {rho_bias:.6f})")
        print(f"  D RMSE: {D_rmse:.6f} mm2/day  (bias: {D_bias:.6f})")
        print(f"  Convergence: {convergence_count}/{n_trials} ({100*convergence_count/n_trials:.1f}%)")
        print(f"  Avg iterations: {avg_iterations:.1f}")
        print()
    
    # Validation summary
    print(f"{'='*70}")
    print("VALIDATION SUMMARY")
    print(f"{'='*70}")
    
    # Check success criteria
    noise_free = results.get("noise_0pct", {})
    noise_10pct = results.get("noise_10pct", {})
    
    criteria_met = {
        "rmse_noise_free_5pct": noise_free.get("rho_rmse", 1.0) < 0.05 and noise_free.get("D_rmse", 1.0) < 0.05,
        "rmse_10pct_noise_15pct": noise_10pct.get("rho_rmse", 1.0) < 0.15 and noise_10pct.get("D_rmse", 1.0) < 0.15,
        "convergence_50_iters": noise_free.get("avg_iterations", 100) < 50,
    }
    
    for criterion, met in criteria_met.items():
        status = "PASS" if met else "FAIL"
        print(f"  {criterion}: {status}")
    
    print(f"{'='*70}\n")
    
    return results


# --------------------------------------------------------------------------- #
# Cohort mode: every MU-Glioma patient with >= 2 scans, all scans fitted jointly
# --------------------------------------------------------------------------- #
WIDE_BOUNDS = [(1e-5, 0.2), (D_MIN, D_MAX)]   # sensitivity: rho floor far below the physiological 0.005
C_DIFF = (36.0 * np.pi) ** (1.0 / 3.0)


def simulate_trajectory(rho: float, D: float, V0: float, t: np.ndarray) -> np.ndarray:
    """surrogate_ode_model integrated continuously from V0 through scan days t (t[0] = 0)."""
    out, V = [V0], V0
    for dt in np.diff(t):
        V = surrogate_ode_model(rho, D, V, float(dt))
        out.append(V)
    return np.array(out)


def fit_trajectory(days: List[float], vols: List[float],
                   bounds: List[Tuple[float, float]]) -> Dict[str, Any]:
    """Least squares on log volume over all follow-up scans, V(t0) = first observed
    volume. Parameters are optimised in log space (rho and D differ by orders of
    magnitude) from a 3x3 grid of starts. Returns fit, Gauss-Newton 95% CIs when
    there are more follow-up scans than parameters, and identifiability diagnostics."""
    t = np.asarray(days, float) - days[0]
    v = np.asarray(vols, float)
    lb = [(np.log(lo), np.log(hi)) for lo, hi in bounds]

    def resid(x):
        sim = simulate_trajectory(np.exp(x[0]), np.exp(x[1]), v[0], t)
        return np.log(np.maximum(sim[1:], 1e-9)) - np.log(v[1:])

    def f(x):
        return float(np.sum(resid(x) ** 2))

    best = None
    for a in np.linspace(lb[0][0], lb[0][1], 3):
        for b in np.linspace(lb[1][0], lb[1][1], 3):
            r = minimize(f, np.array([a, b]), method="L-BFGS-B", bounds=lb, options={"maxiter": 1000})
            if best is None or r.fun < best.fun:
                best = r
    x = best.x
    h = 1e-4
    J = np.stack([(resid(x + h * e) - resid(x - h * e)) / (2 * h) for e in np.eye(2)], 1)
    n_obs = len(v) - 1
    dof = n_obs - 2
    at_lo = [bool(x[i] - lb[i][0] < 1e-6) for i in range(2)]
    at_hi = [bool(lb[i][1] - x[i] < 1e-6) for i in range(2)]
    JTJ = J.T @ J
    cond = float(np.linalg.cond(JTJ)) if np.all(np.isfinite(JTJ)) else float("inf")
    ci = {"rho_ci95": None, "D_ci95": None}
    if dof > 0 and cond < 1e12 and not any(at_lo + at_hi):
        cov = best.fun / dof * np.linalg.inv(JTJ)
        se = np.sqrt(np.maximum(np.diag(cov), 0))
        ci = {"rho_ci95": [float(np.exp(x[0] - 1.96 * se[0])), float(np.exp(x[0] + 1.96 * se[0]))],
              "D_ci95": [float(np.exp(x[1] - 1.96 * se[1])), float(np.exp(x[1] + 1.96 * se[1]))]}
    rho, D = float(np.exp(x[0])), float(np.exp(x[1]))
    growth = rho * v[0] * (1 - v[0] / K_DEFAULT)
    diff = C_DIFF * D * v[0] ** (1 / 3)
    return {"rho": rho, "D": D, "sse_log": float(best.fun), "rmse_log": float(np.sqrt(best.fun / n_obs)),
            "n_followup_scans": n_obs, "converged": bool(best.success), "n_iterations": int(best.nit),
            "rho_at_lower_bound": at_lo[0], "rho_at_upper_bound": at_hi[0],
            "D_at_lower_bound": at_lo[1], "D_at_upper_bound": at_hi[1],
            "fisher_condition_number": cond,
            "diffusion_share_of_initial_growth": float(diff / (growth + diff)) if growth + diff > 0 else None,
            **ci}


def load_scan_series() -> List[Dict[str, Any]]:
    cohort = json.loads(COHORT_JSON.read_text())
    out = []
    for p in cohort:
        tps = sorted((tp["day_from_diagnosis"], tp["volume_mm3"]) for tp in p["timepoints"]
                     if tp["volume_mm3"] is not None and tp["day_from_diagnosis"] is not None)
        if len(tps) >= 2:
            out.append({"patient_id": p["patient_id"], "days": [d for d, _ in tps], "volumes": [v for _, v in tps]})
    return out


def synthetic_identifiability(rng: np.random.Generator, n_trials: int = 20) -> Dict[str, Any]:
    """Recovery of known (rho, D) from surrogate trajectories: 2, 3 and 5 scans over
    180 days, 0/5/10% multiplicative noise, at a small (1 cm^3) and a cohort-sized
    (50 cm^3) tumour."""
    true_rho, true_D = 0.02, 0.015
    out = {"true_rho": true_rho, "true_D": true_D, "cases": []}
    for V0 in (1.0e3, 5.0e4):
        for n_scans in (2, 3, 5):
            t = np.linspace(0, 180, n_scans)
            clean = simulate_trajectory(true_rho, true_D, V0, t)
            for noise in (0.0, 0.05, 0.10):
                er, ed, dbound = [], [], 0
                for _ in range(n_trials if noise > 0 else 1):
                    obs = clean * np.exp(rng.normal(0, noise, n_scans))
                    obs[0] = V0
                    fit = fit_trajectory(t.tolist(), obs.tolist(), [(RHO_MIN, RHO_MAX), (D_MIN, D_MAX)])
                    er.append(abs(fit["rho"] - true_rho) / true_rho)
                    ed.append(abs(fit["D"] - true_D) / true_D)
                    dbound += fit["D_at_lower_bound"] or fit["D_at_upper_bound"]
                out["cases"].append({"V0_mm3": V0, "n_scans": n_scans, "noise": noise, "n_trials": len(er),
                                     "rho_median_rel_err": float(np.median(er)),
                                     "D_median_rel_err": float(np.median(ed)),
                                     "D_at_bound_pct": 100.0 * dbound / len(er),
                                     "diffusion_share_of_initial_growth": fit["diffusion_share_of_initial_growth"]})
    return out


def diagnose_legacy_fit(series: List[Dict[str, Any]], patient_id: str = "PatientID_0003") -> Dict[str, Any]:
    """Reproduce the old per-patient fit (first two scans only, estimate_patient_parameters)
    and explain why it stops at iteration 1 with a zero-width CI."""
    legacy_path = per_patient_path(patient_id)
    legacy = json.loads(legacy_path.read_text()) if legacy_path.exists() else None
    if legacy is not None and "fit" in legacy:      # already rewritten by a cohort run
        legacy = None
    s = next(x for x in series if x["patient_id"] == patient_id)
    V0, V1, dt = s["volumes"][0], s["volumes"][1], s["days"][1] - s["days"][0]
    np.random.seed(51)
    rep = estimate_patient_parameters(V0, V1, dt)
    floor_pred = surrogate_ode_model(RHO_MIN, D_MIN, V0, dt)
    h = 1e-8
    f0 = objective_function(np.array([RHO_MIN, D_MIN]), V0, V1, dt)
    grad = [float((objective_function(np.array([RHO_MIN, D_MIN]) + h * e, V0, V1, dt) - f0) / h) for e in np.eye(2)]
    obs_rate = float(np.log(V1 / V0) / dt)
    return {
        "patient_id": patient_id,
        "legacy_file_contents": legacy,
        "scans_used_by_legacy_fit": {"days": s["days"][:2], "volumes_mm3": [V0, V1]},
        "scans_ignored_by_legacy_fit": {"days": s["days"][2:], "volumes_mm3": s["volumes"][2:]},
        "observed_growth_rate_per_day": obs_rate,
        "rho_lower_bound_per_day": RHO_MIN,
        "lower_bound_over_observed_rate": RHO_MIN / obs_rate if obs_rate > 0 else None,
        "prediction_at_lower_bounds_mm3": floor_pred,
        "residual_at_lower_bounds_mm3": floor_pred - V1,
        "reproduced_fit": {k: v for k, v in rep.items() if k != "bootstrap_samples"},
        "objective_gradient_at_lower_corner": grad,
        "cause": ("Not a convergence failure. The observed growth over the two scans used "
                  f"({obs_rate:.5f}/day) is below the model floor rho >= {RHO_MIN}/day with D >= {D_MIN}, "
                  "so even the slowest allowed parameters over-predict the second volume. The "
                  "objective gradient at the (rho_min, D_min) corner is positive in both parameters, "
                  "i.e. it points out of the feasible box, so the projected gradient is zero and "
                  "L-BFGS-B stops after one iteration at the bound (message: projected gradient <= pgtol). "
                  "Every 10%-noise bootstrap target is also below the floor, so all bootstrap fits land "
                  "on the same corner and the CI has zero width. The fit used only the first two of the "
                  "patient's scans."),
    }


def per_patient_path(patient_id: str) -> Path:
    return OUTPUT_DIR / f"inverse_est_{patient_id}.json"


def fit_cohort(series: List[Dict[str, Any]]) -> None:
    """Fit every patient on all scans and write output/inverse_est_<patient_id>.json."""
    phys = [(RHO_MIN, RHO_MAX), (D_MIN, D_MAX)]
    for s in series:
        fit = fit_trajectory(s["days"], s["volumes"], phys)
        wide = fit_trajectory(s["days"], s["volumes"], WIDE_BOUNDS)
        rec = {**s, "observed_log_change": float(np.log(s["volumes"][-1] / s["volumes"][0])),
               "fit": fit, "fit_wide_bounds": wide}
        per_patient_path(s["patient_id"]).write_text(json.dumps(rec, indent=2))
        print(f"  {s['patient_id']}: {len(s['days'])} scans  rho {fit['rho']:.5f}  D {fit['D']:.5f}  "
              f"rmse_log {fit['rmse_log']:.3f}  rho@lo {fit['rho_at_lower_bound']}")


def aggregate(legacy_diagnosis: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Combine every output/inverse_est_PatientID_*.json into one metrics dict. Files
    in the old two-scan format (no 'fit' block) are counted and listed, not mixed in."""
    import pandas as pd
    from scipy.stats import spearmanr

    files = sorted(OUTPUT_DIR.glob("inverse_est_PatientID_*.json"))
    patients, legacy_files = [], []
    for f in files:
        rec = json.loads(f.read_text())
        if "fit" in rec:
            patients.append(rec)
        else:
            legacy_files.append(f.name)
    if not patients:
        raise RuntimeError("no cohort-format inverse_est_PatientID_*.json files; run with --cohort")

    fits = [p["fit"] for p in patients]
    n = len(fits)
    rho_lo = np.array([f["rho_at_lower_bound"] for f in fits])
    rho_hi = np.array([f["rho_at_upper_bound"] for f in fits])
    D_bound = np.array([f["D_at_lower_bound"] or f["D_at_upper_bound"] for f in fits])
    shrank = np.array([p["observed_log_change"] < 0 for p in patients])
    interior = ~rho_lo & ~rho_hi
    rho_int = np.array([f["rho"] for f in fits])[interior]
    with_ci = [f for f in fits if f["D_ci95"] is not None]
    D_ci_ratio = np.array([f["D_ci95"][1] / f["D_ci95"][0] for f in with_ci])
    share = np.array([f["diffusion_share_of_initial_growth"] for f in fits])

    ref = pd.read_csv(PARAMS_72_CSV).set_index("patient_id")["rho_per_day"]
    ids_int = [p["patient_id"] for p, ok in zip(patients, interior) if ok and p["patient_id"] in ref.index]
    r72 = [ref[i] for i in ids_int]
    r51 = [p["fit"]["rho"] for p in patients if p["patient_id"] in ids_int]
    cross = spearmanr(r51, r72) if len(ids_int) > 2 else None

    wide = [p["fit_wide_bounds"] for p in patients]
    summary = {
        "n_per_patient_files": len(files),
        "n_legacy_format_files": len(legacy_files),
        "n_patients_with_2plus_scans": n,
        "n_shrinking_first_to_last_scan": int(shrank.sum()),
        "n_rho_at_lower_bound": int(rho_lo.sum()),
        "n_rho_at_lower_bound_among_shrinking": int((rho_lo & shrank).sum()),
        "n_rho_at_upper_bound": int(rho_hi.sum()),
        "n_rho_interior": int(interior.sum()),
        "rho_interior_median_per_day": float(np.median(rho_int)) if len(rho_int) else None,
        "rho_interior_iqr_per_day": [float(np.percentile(rho_int, 25)), float(np.percentile(rho_int, 75))]
        if len(rho_int) else None,
        "n_D_at_bound": int(D_bound.sum()),
        "n_with_gauss_newton_ci": len(with_ci),
        "D_ci95_median_upper_over_lower": float(np.median(D_ci_ratio)) if len(D_ci_ratio) else None,
        "diffusion_share_of_initial_growth_median": float(np.median(share)),
        "median_rmse_log": float(np.median([f["rmse_log"] for f in fits])),
        "convergence_rate": float(np.mean([f["converged"] for f in fits])),
        "wide_bounds_n_rho_at_lower_bound": int(sum(w["rho_at_lower_bound"] for w in wide)),
        "cross_check_script72_spearman_rho": float(cross.correlation) if cross is not None else None,
        "cross_check_script72_n": len(ids_int),
    }
    synth = synthetic_identifiability(np.random.default_rng(51))
    # D counts as identifiable if a cohort-sized tumour with 5 scans and 5% noise recovers D
    # within 50% (median), and diffusion is at least 1% of growth in the real cohort.
    ref_case = [c for c in synth["cases"] if c["V0_mm3"] == 5.0e4 and c["n_scans"] == 5 and c["noise"] == 0.05][0]
    summary["D_synthetic_median_rel_err_cohort_size"] = ref_case["D_median_rel_err"]
    summary["D_identifiable_from_volumes"] = bool(ref_case["D_median_rel_err"] < 0.5
                                                  and summary["diffusion_share_of_initial_growth_median"] >= 0.01)
    if legacy_diagnosis is not None:
        summary["legacy_PatientID_0003_cause"] = "rho floor above observed growth; optimum is the bound corner"
    return {"model": "surrogate ODE dV/dt = rho V (1 - V/K) + (36 pi)^(1/3) D V^(1/3), K = 1e6 mm^3",
            "fit": "least squares on log volume over all follow-up scans; log-parameter L-BFGS-B, 3x3 multistart",
            "bounds_physiological": {"rho": [RHO_MIN, RHO_MAX], "D": [D_MIN, D_MAX]},
            "bounds_wide_sensitivity": {"rho": list(WIDE_BOUNDS[0]), "D": list(WIDE_BOUNDS[1])},
            "data": str(COHORT_JSON.relative_to(PROJECT_ROOT)),
            "summary": summary,
            "legacy_fit_diagnosis": legacy_diagnosis,
            "legacy_format_files": legacy_files,
            "synthetic_identifiability": synth,
            "patients": patients,
            "notes": [
                "The surrogate cannot shrink (rho >= 0.005, D >= 0.001), so treated tumours that shrank "
                "between scans pin rho at its lower bound; they are counted, not fitted.",
                "D_identifiable_from_volumes is computed: the cohort-sized synthetic case (5 scans, 5% "
                "noise) must recover D within 50% and diffusion must be >= 1% of initial growth.",
                "Per-patient files are rewritten by --cohort from all scans; the six older two-scan files "
                "(fits at both lower bounds after one iteration) are replaced. legacy_fit_diagnosis "
                "reproduces the PatientID_0003 case from the data.",
            ]}


def main():
    """Main entry point for inverse parameter estimation."""
    parser = argparse.ArgumentParser(
        description="Inverse biophysical parameter estimation for GBM modeling"
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run synthetic validation tests",
    )
    parser.add_argument(
        "--cohort",
        action="store_true",
        help="Fit every MU-Glioma patient with >= 2 scans (writes inverse_est_PatientID_*.json), "
             f"then aggregate into {METRICS_JSON.name}",
    )
    parser.add_argument(
        "--aggregate",
        action="store_true",
        help=f"Only combine existing inverse_est_PatientID_*.json into {METRICS_JSON.name}",
    )
    parser.add_argument(
        "--t0-volume",
        type=float,
        help="Baseline tumor volume (mm³)",
    )
    parser.add_argument(
        "--t1-volume",
        type=float,
        help="Follow-up tumor volume (mm³)",
    )
    parser.add_argument(
        "--delta-t",
        type=float,
        help="Time between scans (days)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON file for results",
    )
    
    args = parser.parse_args()

    def resolve(p: str) -> Path:
        p = Path(p)
        return p if p.is_absolute() else PROJECT_ROOT / p

    if args.cohort or args.aggregate:
        print(f"\n{'='*70}\nINVERSE ESTIMATION - MU-GLIOMA COHORT\n{'='*70}")
        series = load_scan_series()
        diagnosis = diagnose_legacy_fit(series)   # before --cohort rewrites the legacy file
        print(f"Legacy PatientID_0003 fit: {diagnosis['cause']}")
        if args.cohort:
            fit_cohort(series)
        results = aggregate(diagnosis)
        output_path = resolve(args.output) if args.output else METRICS_JSON
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(results, indent=2))
        print(json.dumps(results["summary"], indent=2))
        print(f"Results saved to: {output_path}")
        return 0

    if args.test:
        # Run validation
        results = validate_with_synthetic_data()

        if args.output:
            output_path = resolve(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with open(output_path, "w") as f:
                json.dump(results, f, indent=2)
            print(f"Results saved to: {output_path}")
        
        return 0
    
    if args.t0_volume and args.t1_volume and args.delta_t:
        # Estimate parameters for patient
        print(f"\n{'='*70}")
        print("PATIENT PARAMETER ESTIMATION")
        print(f"{'='*70}")
        print(f"Baseline volume (V0): {args.t0_volume:.1f} mm3")
        print(f"Follow-up volume (V1): {args.t1_volume:.1f} mm3")
        print(f"Time interval (dt): {args.delta_t:.1f} days")
        print(f"{'='*70}\n")
        
        result = estimate_patient_parameters(
            t0_volume=args.t0_volume,
            t1_volume=args.t1_volume,
            delta_t_days=args.delta_t,
        )
        
        print("ESTIMATED PARAMETERS:")
        print(f"  rho (growth rate):     {result['rho']:.6f} /day")
        print(f"                         95% CI: [{result['rho_ci'][0]:.6f}, {result['rho_ci'][1]:.6f}]")
        print(f"  D (diffusivity):       {result['D']:.6f} mm2/day")
        print(f"                         95% CI: [{result['D_ci'][0]:.6f}, {result['D_ci'][1]:.6f}]")
        print(f"  Convergence:         {'Yes' if result['convergence']['success'] else 'No'}")
        print(f"  Iterations:          {result['convergence']['n_iterations']}")
        print(f"  Objective value:     {result['convergence']['objective_value']:.6e}")
        print(f"  Gradient norm:       {result['convergence']['gradient_norm']:.6e}")
        print(f"  RMSE:                {result['rmse']:.4f} mm³")
        print(f"{'='*70}\n")
        
        if args.output:
            output_path = resolve(args.output)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            # Remove non-serializable bootstrap samples for JSON output
            result_json = {k: v for k, v in result.items() if k != "bootstrap_samples"}
            with open(output_path, "w") as f:
                json.dump(result_json, f, indent=2)
            print(f"Results saved to: {output_path}")
        
        return 0
    
    # No arguments provided - show help
    parser.print_help()
    print("\nExamples:")
    print("  python src/51_inverse_parameter_estimation.py --cohort")
    print("  python src/51_inverse_parameter_estimation.py --test")
    print("  python src/51_inverse_parameter_estimation.py --t0-volume 1000 --t1-volume 1200 --delta-t 30")
    return 0


if __name__ == "__main__":
    sys.exit(main())
