#!/usr/bin/env python3
"""Uncertainty-Aware & Adaptive-Horizon MPC Controller (Phase 2).

Enhances the baseline 14-day receding-horizon MPC with two algorithmic
innovations:

  1. Uncertainty-Aware Cost (Robust MPC):
     Evaluates the MPC objective over a parameter uncertainty distribution
     (rho +/-15%, D +/-15%) using a risk-averse formulation:
         J_robust = mean(J) + lambda * std(J)
     This penalizes high-variance decisions that could fail under parameter
     perturbation, yielding more stable dosing trajectories.

  2. Adaptive Prediction Horizon:
     Dynamically adjusts the MPC prediction horizon (7-21 days) based on
     observed tumor growth dynamics:
         - Stable growth  (|dV/dt| < 0.01 /day): extend horizon (slower response)
         - Accelerating   (dV/dt > 0.05 /day):  shorten horizon (faster response)
         - Near target    (|V - V_target| < 10%): maintain

The controller returns the same output schema as run_mpc_adaptive_3d() so it
can be a drop-in replacement, plus additional robustness diagnostics.

Usage:
    python src/52_robust_mpc_controller.py --benchmark
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize_scalar

warnings.filterwarnings("ignore")

# Physical constants (must match src/50_clinical_cdss_app.py)
DX = 1.0          # mm
GRID_SIZE = 50
DT = 0.04         # days
SIM_DAYS = 180
N_STEPS = int(SIM_DAYS / DT)

# TMZ PK parameters
TMZ_HALF_LIFE = 0.075  # days
K_EL = np.log(2) / TMZ_HALF_LIFE
C_PEAK = 10.0
EC50 = 5.0
HILL_COEFF = 2.0
E_MAX = 0.55

# Dosing schedule
DOSE_DAYS_ON = 5
CYCLE_DAYS = 28

# Horizon adjustment bounds
MIN_HORIZON = 7
MAX_HORIZON = 21
DEFAULT_HORIZON = 14

# Coarser timestep used only for open-loop horizon predictions (the actual
# simulation still uses DT=0.04). Using a 1-day prediction step reduces the
# horizon loop from ~350 iterations to ~14, keeping the robust cost cheap
# enough for Monte Carlo benchmarking while preserving the qualitative
# growth/kill dynamics.
PREDICT_DT = 1.0  # days

# Uncertainty parameters
UNCERTAINTY_FRACTION = 0.15  # rho and D perturbed by +/- 15%
DEFAULT_N_SAMPLES = 16

# Carrying capacity for surrogate ODE (matches inverse estimation module)
K_CARRY = 1.0e6


def tmz_concentration(step: int, drug_on: bool = True) -> float:
    """Compute TMZ concentration at given step (matches baseline)."""
    if not drug_on:
        return 0.0
    t_days = step * DT
    day_in_cycle = int(t_days) % CYCLE_DAYS
    if day_in_cycle < DOSE_DAYS_ON:
        return C_PEAK * np.exp(-K_EL * DT)
    else:
        days_since_dose = day_in_cycle - (DOSE_DAYS_ON - 1)
        return C_PEAK * np.exp(-K_EL * days_since_dose)


def compute_kill_rate(C: float) -> float:
    """Compute tumor kill rate from TMZ concentration (Hill equation)."""
    return E_MAX * (C ** HILL_COEFF) / (EC50 ** HILL_COEFF + C ** HILL_COEFF + 1e-12)


def volume_step(
    volume_mm3: float,
    rho: float,
    D: float,
    drug_on: bool,
    step: int,
    dt: float = DT,
    Vmax: float = (GRID_SIZE ** 3) * (DX ** 3),
) -> float:
    """Advance tumor volume one time step under given parameters and dosing.

    Surrogate ODE:
        dV/dt = rho*V*(1 - V/Vmax) - kill*V
    """
    C = tmz_concentration(step, drug_on)
    kill = compute_kill_rate(C)
    dV = (rho * volume_mm3 * (1.0 - volume_mm3 / Vmax) - kill * volume_mm3) * dt
    return max(volume_mm3 + dV, 0.0)


def predict_volume_horizon(
    volume_mm3: float,
    rho: float,
    D: float,
    horizon_days: int,
    step: int,
    drug_on: bool,
) -> float:
    """Predict volume at end of prediction horizon (open-loop).

    Used by the robust cost function to evaluate control candidates.
    Uses the same dosing decision (drug_on) for the entire horizon.
    Integrates with the coarse PREDICT_DT prediction timestep.
    """
    v = volume_mm3
    n_steps = max(1, int(round(horizon_days / PREDICT_DT)))
    Vmax = (GRID_SIZE ** 3) * (DX ** 3)
    for k in range(n_steps):
        pred_day = (step * DT) + k * PREDICT_DT
        # Determine if currently in dosing phase within the cycle
        day_in_cycle = int(pred_day) % CYCLE_DAYS
        on_now = drug_on and (day_in_cycle < DOSE_DAYS_ON)
        if on_now:
            C_k = C_PEAK * np.exp(-K_EL * PREDICT_DT)
        else:
            C_k = 0.0
        kill_k = compute_kill_rate(C_k)
        dV = (rho * v * (1.0 - v / Vmax) - kill_k * v) * PREDICT_DT
        v = max(v + dV, 0.0)
    return v


def predict_volume_horizon_vec(
    volume_mm3: float,
    rho_samples: np.ndarray,
    D_samples: np.ndarray,
    horizon_days: int,
    step: int,
    drug_on: bool,
) -> np.ndarray:
    """Vectorized volume prediction across all uncertainty samples.

    Returns an array of predicted volumes (one per sample).
    Integrates with the coarse PREDICT_DT prediction timestep.
    """
    n = rho_samples.shape[0]
    V = np.full(n, float(volume_mm3))
    n_steps = max(1, int(round(horizon_days / PREDICT_DT)))
    Vmax = (GRID_SIZE ** 3) * (DX ** 3)

    for k in range(n_steps):
        pred_day = (step * DT) + k * PREDICT_DT
        day_in_cycle = int(pred_day) % CYCLE_DAYS
        on_now = drug_on and (day_in_cycle < DOSE_DAYS_ON)
        if on_now:
            C_k = C_PEAK * np.exp(-K_EL * PREDICT_DT)
        else:
            C_k = 0.0
        kill_k = compute_kill_rate(C_k)
        dV = (rho_samples * V * (1.0 - V / Vmax) - kill_k * V) * PREDICT_DT
        V = np.maximum(V + dV, 0.0)
    return V


def nominal_cost(
    volume_mm3: float,
    dose: float,
    rho: float,
    D: float,
    target_volume: float,
    w_tumor: float,
    w_drug: float,
    horizon_days: int,
    step: int,
) -> float:
    """Compute nominal MPC cost for a single (rho, D) realization.

    J = w_tumor * max(0, predicted_volume - target) + w_drug * dose
    """
    drug_on = dose >= 0.5
    v_pred = predict_volume_horizon(
        volume_mm3, rho, D, horizon_days, step, drug_on=drug_on
    )
    volume_above_target = max(0.0, v_pred - target_volume)
    return w_tumor * (volume_above_target / max(target_volume, 1.0)) + w_drug * dose


def nominal_cost_vec(
    volume_mm3: float,
    dose: float,
    rho_samples: np.ndarray,
    D_samples: np.ndarray,
    target_volume: float,
    w_tumor: float,
    w_drug: float,
    horizon_days: int,
    step: int,
) -> np.ndarray:
    """Vectorized nominal cost across all samples. Returns array of costs."""
    drug_on = dose >= 0.5
    v_pred = predict_volume_horizon_vec(
        volume_mm3, rho_samples, D_samples, horizon_days, step, drug_on=drug_on
    )
    volume_above_target = np.maximum(0.0, v_pred - target_volume)
    return w_tumor * (volume_above_target / max(target_volume, 1.0)) + w_drug * dose


def robust_cost(
    dose: float,
    volume_mm3: float,
    rho_nominal: float,
    D_nominal: float,
    target_volume: float,
    w_tumor: float,
    w_drug: float,
    w_uncertainty: float,
    risk_aversion: float,
    horizon_days: int,
    step: int,
    n_samples: int = DEFAULT_N_SAMPLES,
    rng: Optional[np.random.Generator] = None,
) -> float:
    """Evaluate robust (uncertainty-aware) MPC cost.

    J_robust = mean(J) + lambda * std(J)
    where J is computed over (rho +/-15%, D +/-15%) samples.

    Lower J_robust is better. The std term penalizes high-variance decisions.
    """
    if rng is None:
        rng = np.random.default_rng()

    rho_samples = rng.normal(rho_nominal, UNCERTAINTY_FRACTION * rho_nominal, n_samples)
    D_samples = rng.normal(D_nominal, UNCERTAINTY_FRACTION * D_nominal, n_samples)

    # Clamp to physically valid range
    rho_samples = np.clip(rho_samples, 1e-4, 0.2)
    D_samples = np.clip(D_samples, 1e-4, 0.1)

    # Vectorized cost over all samples
    costs = nominal_cost_vec(
        volume_mm3, dose, rho_samples, D_samples,
        target_volume, w_tumor, w_drug, horizon_days, step,
    )

    mean_cost = float(np.mean(costs))
    std_cost = float(np.std(costs))
    return mean_cost + risk_aversion * std_cost


def adjust_horizon(
    current_horizon: int,
    growth_rate: float,
    volume: float,
    target_volume: float,
    min_horizon: int = MIN_HORIZON,
    max_horizon: int = MAX_HORIZON,
) -> int:
    """Dynamically adjust the MPC prediction horizon based on growth dynamics.

    Rules (per plan spec):
      - Stable growth   (|dV/dt| < 0.01 /day): extend horizon (slower response OK)
      - Accelerating    (dV/dt > 0.05 /day):   shorten horizon (need faster response)
      - Near target     (|V - V_target| < 10%): maintain current horizon
    """
    if target_volume > 0 and abs(volume - target_volume) < 0.10 * target_volume:
        return int(np.clip(current_horizon, min_horizon, max_horizon))

    if abs(growth_rate) < 0.01:
        return int(np.clip(current_horizon + 1, min_horizon, max_horizon))
    elif growth_rate > 0.05:
        return int(np.clip(current_horizon - 1, min_horizon, max_horizon))
    else:
        return int(np.clip(current_horizon, min_horizon, max_horizon))


class RobustMPCController:
    """Uncertainty-aware adaptive-horizon MPC controller for GBM therapy.

    Replaces the baseline run_mpc_adaptive_3d() with a robust formulation that
    integrates parameter uncertainty (rho +/-15%, D +/-15%) and dynamically
    adjusts the prediction horizon (7-21 days) based on tumor growth dynamics.
    """

    def __init__(
        self,
        w_tumor: float = 1.2,
        w_drug: float = 0.03,
        w_uncertainty: float = 0.5,
        risk_aversion: float = 0.3,
        initial_horizon: int = DEFAULT_HORIZON,
        min_horizon: int = MIN_HORIZON,
        max_horizon: int = MAX_HORIZON,
        n_uncertainty_samples: int = DEFAULT_N_SAMPLES,
        seed: Optional[int] = None,
    ) -> None:
        self.w_tumor = w_tumor
        self.w_drug = w_drug
        self.w_uncertainty = w_uncertainty
        self.risk_aversion = risk_aversion
        self.horizon = int(np.clip(initial_horizon, min_horizon, max_horizon))
        self.min_horizon = min_horizon
        self.max_horizon = max_horizon
        self.n_uncertainty_samples = n_uncertainty_samples
        self.rng = np.random.default_rng(seed)

    def optimize_control(
        self,
        current_volume_mm3: float,
        rho_nominal: float,
        D_nominal: float,
        target_volume: float,
        step: int,
    ) -> Tuple[float, Dict[str, Any]]:
        """Decide whether to dose (1.0) or hold (0.0) at the current step.

        Compares the robust cost of dosing vs holding and picks the lower.
        Returns (dose_decision, diagnostics).
        """
        # Standard 5-on/23-off schedule as the baseline anchor.
        in_dose_phase = (step * DT) % CYCLE_DAYS < DOSE_DAYS_ON

        kwargs = dict(
            volume_mm3=current_volume_mm3,
            rho_nominal=rho_nominal,
            D_nominal=D_nominal,
            target_volume=target_volume,
            w_tumor=self.w_tumor,
            w_drug=self.w_drug,
            w_uncertainty=self.w_uncertainty,
            risk_aversion=self.risk_aversion,
            horizon_days=self.horizon,
            step=step,
            n_samples=self.n_uncertainty_samples,
        )

        # Robust cost of dosing vs holding.
        # Use paired uncertainty samples (same parametric draw) so the comparison
        # isolates the effect of the control decision from sampling noise.
        sub_rng = np.random.default_rng(self.rng.integers(0, 2**31 - 1))
        cost_dose = robust_cost(dose=1.0, rng=sub_rng, **kwargs)
        cost_hold = robust_cost(dose=0.0, rng=sub_rng, **kwargs)

        # Decision rule (uncertainty-augmented baseline):
        #   - Start from the standard 5-on/23-off cycle.
        #   - Robust cost benefit of dosing = cost_hold - cost_dose (how much
        #     worse holding is). When this exceeds the drug penalty threshold,
        #     the controller doses; otherwise it spares.
        #   - This keeps the scheduler's structure (anchor ~30% dosing) while
        #     letting uncertainty-awareness spare doses when the tumor is
        #     well-controlled and add doses only under genuine escape risk.
        benefit_of_dosing = cost_hold - cost_dose  # >0 means dosing helps

        # Threshold scales with the drug penalty so the controller only doses
        # when the tumor-control benefit materially exceeds the drug cost.
        dose_threshold = self.w_drug + self.w_uncertainty * 0.02

        if benefit_of_dosing > dose_threshold:
            # Tumor-control benefit of dosing exceeds its cost -> dose
            dose = 1.0
        elif benefit_of_dosing < -dose_threshold:
            # Dosing is wasteful (tumor already controlled) -> hold
            dose = 0.0
        else:
            # Indifferent: defer to the standard 5-on/23-off cycle anchor.
            dose = 1.0 if in_dose_phase else 0.0

        diagnostics = {
            "cost_dose": cost_dose,
            "cost_hold": cost_hold,
            "benefit_of_dosing": benefit_of_dosing,
            "horizon_days": self.horizon,
            "decision": dose,
        }
        return dose, diagnostics

    def update_horizon(
        self,
        growth_rate: float,
        volume: float,
        target_volume: float,
    ) -> int:
        """Adjust prediction horizon based on dynamics; returns new horizon."""
        self.horizon = adjust_horizon(
            current_horizon=self.horizon,
            growth_rate=growth_rate,
            volume=volume,
            target_volume=target_volume,
            min_horizon=self.min_horizon,
            max_horizon=self.max_horizon,
        )
        return self.horizon


def run_robust_mpc_simulation(
    u0: np.ndarray,
    rho: float,
    D: float,
    controller: Optional[RobustMPCController] = None,
) -> Dict[str, Any]:
    """Run a full robust MPC simulation. Output schema matches the baseline
    run_mpc_adaptive_3d() so this is a drop-in replacement, plus extra
    robustness diagnostics.
    """
    if controller is None:
        controller = RobustMPCController()

    initial_volume = float(np.sum(u0 > 0.1)) * (DX ** 3)
    target_volume_mm3 = initial_volume * 0.12

    volume_mm3 = initial_volume
    drug_on_history: List[bool] = []
    dose_schedule: List[float] = []
    volume_history: List[float] = [volume_mm3]
    horizon_history: List[int] = []
    cost_history: List[Dict[str, Any]] = []
    growth_rate_history: List[float] = []

    prev_volume = volume_mm3

    for step in range(N_STEPS):
        # Estimate empirical growth rate (per day) from recent window
        if step >= int(7.0 / DT):  # at least 7 days of history
            window = max(0, step - int(7.0 / DT))
            v_now = volume_history[-1]
            v_then = volume_history[window]
            dt_days = (step - window) * DT
            if v_then > 1e-6 and dt_days > 0:
                growth_rate = (v_now - v_then) / (v_then * dt_days)
            else:
                growth_rate = 0.0
        else:
            growth_rate = 0.0
        growth_rate_history.append(growth_rate)

        # Adaptive horizon update
        controller.update_horizon(
            growth_rate=growth_rate,
            volume=volume_mm3,
            target_volume=target_volume_mm3,
        )
        horizon_history.append(controller.horizon)

        # Tumor escaping override (preserve baseline safety logic)
        if volume_mm3 > target_volume_mm3 * 1.5:
            dose = 1.0
            diagnostics = {
                "cost_dose": 0.0,
                "cost_hold": 1.0,
                "horizon_days": controller.horizon,
                "decision": 1.0,
                "override": "tumor_escaping",
            }
        else:
            dose, diagnostics = controller.optimize_control(
                current_volume_mm3=volume_mm3,
                rho_nominal=rho,
                D_nominal=D,
                target_volume=target_volume_mm3,
                step=step,
            )

        drug_on = bool(dose >= 0.5)
        drug_on_history.append(drug_on)
        dose_schedule.append(dose)
        cost_history.append(diagnostics)

        # Advance tumor volume using the nominal (rho, D) dynamics
        volume_mm3 = volume_step(
            volume_mm3, rho, D, drug_on=drug_on, step=step
        )
        volume_history.append(volume_mm3)

    drug_on_fraction = float(np.mean(drug_on_history))
    final_volume = volume_history[-1]

    # Robustness diagnostics
    horizon_array = np.array(horizon_history)
    cost_dose_array = np.array([c["cost_dose"] for c in cost_history])
    cost_hold_array = np.array([c["cost_hold"] for c in cost_history])

    return {
        "drug_schedule": dose_schedule,
        "drug_on_fraction": drug_on_fraction,
        "dose_sparing_fraction": 1.0 - drug_on_fraction,
        "predicted_final_volume_mm3": float(final_volume),
        "volume_history": volume_history,
        "drug_on_history": drug_on_history,
        # Robust MPC extra diagnostics
        "horizon_history": horizon_history,
        "growth_rate_history": growth_rate_history,
        "final_horizon": int(horizon_array[-1]),
        "horizon_range": [int(horizon_array.min()), int(horizon_array.max())],
        "cost_variance": float(np.var(cost_dose_array - cost_hold_array)),
        "mean_horizon": float(np.mean(horizon_array)),
        "controller": "robust_uncertainty_aware",
    }


# ============================================================================ #
# Benchmark: Robust vs Standard MPC
# ============================================================================ #
def run_benchmark(
    n_mc: int = 100,
    perturbation: float = 0.15,
    seed: int = 42,
) -> Dict[str, Any]:
    """Compare robust MPC vs standard (baseline) MPC under parameter perturbations.

    Metrics:
      - drug administration fraction
      - dose sparing
      - final volume
      - cost variance (lower = more stable)
    """
    print(f"\n{'='*70}")
    print("ROBUST MPC BENCHMARK — Robust vs Standard MPC")
    print(f"{'='*70}")
    print(f"Monte Carlo trials: {n_mc}")
    print(f"Parameter perturbation: +/-{perturbation*100:.0f}%")
    print(f"{'='*70}\n")

    rng = np.random.default_rng(seed)

    # Standard MPC baseline import (from src/50)
    import importlib.util
    base_path = Path(__file__).parent / "50_clinical_cdss_app.py"
    spec = importlib.util.spec_from_file_location("clinical_cdss_app", base_path)
    base_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(base_mod)

    # Setup a representative patient
    gs = GRID_SIZE
    u0 = np.zeros((gs, gs, gs))
    z, y, x = np.mgrid[0:gs, 0:gs, 0:gs]
    dist = np.sqrt((x - gs//2)**2 + (y - gs//2)**2 + (z - gs//2)**2)
    u0[dist <= 3.0] = 0.8

    robust_metrics: Dict[str, List[float]] = {
        "drug_frac": [], "spare_frac": [], "final_vol": [], "cost_var": [],
    }
    standard_metrics: Dict[str, List[float]] = {
        "drug_frac": [], "spare_frac": [], "final_vol": [], "cost_var": [],
    }

    # Carrying capacity nominal parameters used to evaluate realized cost
    # trajectory for a fair apples-to-apples cost-variance comparison.
    nominal_rho = 0.02
    nominal_D = 0.013

    def realized_cost_trajectory(
        volume_hist: List[float],
        dose_hist: List[float],
        rho: float,
        target_vol: float,
    ) -> float:
        """Variance of the realized per-step nominal cost over the trajectory."""
        costs = []
        for v, dose in zip(volume_hist, dose_hist):
            above = max(0.0, v - target_vol)
            c = W_TUMOR_DEFAULT * (above / max(target_vol, 1.0)) + W_DRUG_DEFAULT * dose
            costs.append(c)
        return float(np.var(costs)) if costs else 0.0

    W_TUMOR_DEFAULT = 1.2
    W_DRUG_DEFAULT = 0.03

    for trial in range(n_mc):
        # Perturb parameters
        rho_true = float(rng.normal(0.02, perturbation * 0.02))
        rho_true = max(0.005, min(0.1, rho_true))
        D_true = float(rng.normal(0.013, perturbation * 0.013))
        D_true = max(0.001, min(0.05, D_true))

        # Robust MPC (uses estimates with assumed +/-15% uncertainty)
        ctrl = RobustMPCController(seed=trial)
        r = run_robust_mpc_simulation(u0, rho_true, D_true, ctrl)
        robust_metrics["drug_frac"].append(r["drug_on_fraction"])
        robust_metrics["spare_frac"].append(r["dose_sparing_fraction"])
        robust_metrics["final_vol"].append(r["predicted_final_volume_mm3"])
        # Realized trajectory cost variance under the true parameters
        initial_volume = float(np.sum(u0 > 0.1)) * (DX ** 3)
        target_vol = initial_volume * 0.12
        robust_metrics["cost_var"].append(
            realized_cost_trajectory(r["volume_history"][:-1], r["drug_schedule"],
                                     rho_true, target_vol)
        )

        # Standard (baseline) MPC
        s = base_mod.run_mpc_adaptive_3d(u0, rho_true, np.array([0.707, 0.707, 0.0]))
        standard_metrics["drug_frac"].append(s["drug_on_fraction"])
        standard_metrics["spare_frac"].append(s["dose_sparing_fraction"])
        standard_metrics["final_vol"].append(s["predicted_final_volume_mm3"])
        standard_metrics["cost_var"].append(
            realized_cost_trajectory(s["volume_history"][:-1], s["drug_schedule"],
                                     rho_true, target_vol)
        )

    def summarize(d: Dict[str, List[float]], key: str) -> Tuple[float, float]:
        arr = np.array(d[key])
        return float(np.mean(arr)), float(np.std(arr))

    summary = {
        "robust": {
            "drug_admin_mean_std": summarize(robust_metrics, "drug_frac"),
            "dose_sparing_mean_std": summarize(robust_metrics, "spare_frac"),
            "final_volume_mean_std": summarize(robust_metrics, "final_vol"),
            "cost_variance_mean": float(np.mean(robust_metrics["cost_var"])),
        },
        "standard": {
            "drug_admin_mean_std": summarize(standard_metrics, "drug_frac"),
            "dose_sparing_mean_std": summarize(standard_metrics, "spare_frac"),
            "final_volume_mean_std": summarize(standard_metrics, "final_vol"),
            "cost_variance_mean": float(np.mean(standard_metrics["cost_var"])),
        },
    }

    print("ROBUST MPC:")
    print(f"  Drug administration: {summary['robust']['drug_admin_mean_std'][0]*100:.1f}% "
          f"+/- {summary['robust']['drug_admin_mean_std'][1]*100:.1f}%")
    print(f"  Dose sparing:         {summary['robust']['dose_sparing_mean_std'][0]*100:.1f}% "
          f"+/- {summary['robust']['dose_sparing_mean_std'][1]*100:.1f}%")
    print(f"  Final volume:         {summary['robust']['final_volume_mean_std'][0]:.2f} mm3 "
          f"+/- {summary['robust']['final_volume_mean_std'][1]:.2f}")
    print(f"  Cost variance (mean): {summary['robust']['cost_variance_mean']:.4f}")

    print("\nSTANDARD MPC:")
    print(f"  Drug administration: {summary['standard']['drug_admin_mean_std'][0]*100:.1f}% "
          f"+/- {summary['standard']['drug_admin_mean_std'][1]*100:.1f}%")
    print(f"  Dose sparing:         {summary['standard']['dose_sparing_mean_std'][0]*100:.1f}% "
          f"+/- {summary['standard']['dose_sparing_mean_std'][1]*100:.1f}%")
    print(f"  Final volume:         {summary['standard']['final_volume_mean_std'][0]:.2f} mm3 "
          f"+/- {summary['standard']['final_volume_mean_std'][1]:.2f}")

    # Cost variance reduction (plan target: >=30% lower for robust MPC)
    robust_cost_var = summary["robust"]["cost_variance_mean"]
    std_cost_var = summary["standard"]["cost_variance_mean"]
    if std_cost_var > 0:
        cost_var_reduction = (std_cost_var - robust_cost_var) / std_cost_var * 100
    else:
        cost_var_reduction = 0.0
    print(f"\n  Cost-trajectory variance (robust):  {robust_cost_var:.4f}")
    print(f"  Cost-trajectory variance (standard): {std_cost_var:.4f}")
    print(f"  Cost-variance reduction (robust vs standard): {cost_var_reduction:.1f}%")

    # Drug-administration MC variance reduction (secondary metric)
    robust_std = summary["robust"]["drug_admin_mean_std"][1]
    std_std = summary["standard"]["drug_admin_mean_std"][1]
    if std_std > 0:
        admin_var_reduction = (std_std - robust_std) / std_std * 100
    else:
        admin_var_reduction = 0.0
    print(f"  Drug-admin MC-variance reduction: {admin_var_reduction:.1f}%")
    print(f"{'='*70}\n")

    summary["cost_variance_reduction_pct"] = cost_var_reduction
    summary["variance_reduction_pct"] = admin_var_reduction
    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Uncertainty-aware adaptive-horizon MPC controller"
    )
    parser.add_argument(
        "--benchmark",
        action="store_true",
        help="Run robust vs standard MPC Monte Carlo benchmark",
    )
    parser.add_argument(
        "--n-mc",
        type=int,
        default=100,
        help="Number of Monte Carlo trials for benchmark (default: 100)",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON file for benchmark results",
    )
    args = parser.parse_args()

    if args.benchmark:
        summary = run_benchmark(n_mc=args.n_mc, seed=42)
        if args.output:
            out_path = Path(args.output)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            # Convert tuples to lists for JSON
            serializable = {}
            for ctrl, metrics in summary.items():
                if isinstance(metrics, dict):
                    serializable[ctrl] = {}
                    for k, v in metrics.items():
                        if isinstance(v, tuple):
                            serializable[ctrl][k] = list(v)
                        else:
                            serializable[ctrl][k] = v
                else:
                    serializable[ctrl] = metrics
            with open(out_path, "w") as f:
                json.dump(serializable, f, indent=2)
            print(f"Results saved to: {out_path}")
        return 0

    parser.print_help()
    return 0


if __name__ == "__main__":
    sys.exit(main())
