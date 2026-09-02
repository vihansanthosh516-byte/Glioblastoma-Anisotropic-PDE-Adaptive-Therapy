#!/usr/bin/env python3
"""Trust signal extraction from the robust MPC controller.

The robust MPC evaluates each control decision over a +/-15% parameter
perturbation ensemble. The margin `benefit_of_dosing = cost_hold - cost_dose`
is a scalar; its *distribution* across the ensemble is an uncertainty signal.
This module exposes that distribution and converts it into a normalized
trust score in [0, 1]:

    trust = |benefit_mean| / (|benefit_mean| + lambda * benefit_std)

- trust -> 1 : the mechanistic model is confident (large margin, small spread)
- trust -> 0 : the model is uncertain (margin comparable to / smaller than noise)

A high trust means "defer to MPC"; a low trust means "delegate to RL".

IMPORTANT — consistent PK model: the base MPC horizon predictor collapses the
TMZ half-life (0.075 d) to ~zero concentration over a 1-day prediction step,
so dosing appears to do nothing. That makes its uncertainty signal degenerate.
This module uses a corrected effective-daily-kill predictor that matches the
simulator's PK (time-averaged kill over a dosing day), so the benefit margin
is meaningful and the trust score is informative.
"""
from __future__ import annotations

import importlib.util
import math
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
_mpc_spec = importlib.util.spec_from_file_location(
    "robust_mpc", ROOT / "src" / "52_robust_mpc_controller.py"
)
_mpc = importlib.util.module_from_spec(_mpc_spec)
assert _mpc_spec.loader is not None
_mpc_spec.loader.exec_module(_mpc)

RobustMPCController = _mpc.RobustMPCController
UNCERTAINTY_FRACTION = _mpc.UNCERTAINTY_FRACTION
DEFAULT_N_SAMPLES = _mpc.DEFAULT_N_SAMPLES
C_PEAK = _mpc.C_PEAK
K_EL = _mpc.K_EL
CYCLE_DAYS = _mpc.CYCLE_DAYS
DOSE_DAYS_ON = _mpc.DOSE_DAYS_ON


# --------------------------------------------------------------------------- #
# Corrected predictor (consistent with the simulator's PK)
# --------------------------------------------------------------------------- #
def _effective_daily_kill(
    rho: float,
    D: float,
    on_dose_day: bool,
) -> float:
    """Time-averaged kill on one day under the simulator's PK semantics.

    The instantaneous Hill kill varies within a dosing day (peak right after a
    bolus, then ~0.075 d half-life decay). Averaging over the day with the
    simulator's DT=0.04 substeps yields the effective daily kill used here.
    """
    dt = 0.04
    n = int(round(1.0 / dt))
    if not on_dose_day:
        return 0.0
    kills = []
    for i in range(n):
        t_days = i * dt
        if t_days < DOSE_DAYS_ON:
            C = C_PEAK * math.exp(-K_EL * t_days)
        else:
            C = C_PEAK * math.exp(-K_EL * (t_days - (DOSE_DAYS_ON - 1)))
        kills.append(_mpc.compute_kill_rate(C))
    return float(np.mean(kills))


def _predict_volume_corrected(
    volume_mm3: float,
    rho: float,
    D: float,
    horizon_days: int,
    step: int,
    drug_on: bool,
    target_volume: float,
    w_tumor: float,
    w_drug: float,
) -> float:
    """Horizon volume prediction with corrected daily kill.

    Walks the horizon in 1-day increments, applying the effective daily kill on
    dosing days (5-on/28-off cycle anchored at step) and logistic growth
    otherwise. Returns the predicted volume at the end of the horizon.
    """
    v = volume_mm3
    Vmax = (_mpc.GRID_SIZE ** 3) * (_mpc.DX ** 3)
    start_day = step * _mpc.DT
    for k in range(horizon_days):
        day = start_day + k
        in_dose_phase = drug_on and (int(day) % CYCLE_DAYS < DOSE_DAYS_ON)
        kill = _effective_daily_kill(rho, D, in_dose_phase)
        dV = (rho * v * (1.0 - v / Vmax) - kill * v)
        v = max(v + dV, 0.0)
    return v


def _benefit_margin_corrected(
    volume_mm3: float,
    rho: float,
    D: float,
    horizon_days: int,
    step: int,
    target_volume: float,
    w_tumor: float,
    w_drug: float,
) -> float:
    """benefit_of_dosing = cost_hold - cost_dose (corrected PK)."""
    v_hold = _predict_volume_corrected(
        volume_mm3, rho, D, horizon_days, step, drug_on=False,
        target_volume=target_volume, w_tumor=w_tumor, w_drug=w_drug,
    )
    v_dose = _predict_volume_corrected(
        volume_mm3, rho, D, horizon_days, step, drug_on=True,
        target_volume=target_volume, w_tumor=w_tumor, w_drug=w_drug,
    )
    cost_hold = w_tumor * (max(0.0, v_hold - target_volume) / max(target_volume, 1.0))
    cost_dose = w_tumor * (max(0.0, v_dose - target_volume) / max(target_volume, 1.0)) + w_drug
    return cost_hold - cost_dose


def benefit_distribution(
    controller: RobustMPCController,
    current_volume_mm3: float,
    rho_nominal: float,
    D_nominal: float,
    target_volume: float,
    step: int,
    n_samples: int = DEFAULT_N_SAMPLES,
    seed: Optional[int] = None,
) -> np.ndarray:
    """Return per-sample benefit-of-dosing across the parameter ensemble.

    Uses the same +/-15% parameter perturbation scheme as the robust cost, but
    with the corrected PK predictor so the margin is meaningful.
    """
    rng = np.random.default_rng(seed)
    rho_samples = rng.normal(rho_nominal, UNCERTAINTY_FRACTION * rho_nominal, n_samples)
    D_samples = rng.normal(D_nominal, UNCERTAINTY_FRACTION * D_nominal, n_samples)
    rho_samples = np.clip(rho_samples, 1e-4, 0.2)
    D_samples = np.clip(D_samples, 1e-4, 0.1)

    margins = np.empty(n_samples)
    for i in range(n_samples):
        margins[i] = _benefit_margin_corrected(
            current_volume_mm3, rho_samples[i], D_samples[i],
            controller.horizon, step, target_volume,
            controller.w_tumor, controller.w_drug,
        )
    return margins


def compute_trust_signal(
    controller: RobustMPCController,
    current_volume_mm3: float,
    rho_nominal: float,
    D_nominal: float,
    target_volume: float,
    step: int,
    n_samples: int = DEFAULT_N_SAMPLES,
    lambda_scale: float = 1.0,
    seed: Optional[int] = None,
) -> Dict[str, Any]:
    """Compute the normalized trust signal plus diagnostics.

    Returns:
        dict with keys:
          - trust: float in [0, 1]; 1 = confident MPC, 0 = uncertain
          - benefit_mean: mean benefit-of-dosing over the ensemble
          - benefit_std: std of benefit-of-dosing over the ensemble
          - benefit_cv: coefficient of variation (std/|mean|), inf-safe
          - mpc_decision: the corrected controller decision (dose if benefit>0)
          - decision_confident: bool, trust > 0.5
    """
    benefits = benefit_distribution(
        controller, current_volume_mm3, rho_nominal, D_nominal,
        target_volume, step, n_samples=n_samples, seed=seed,
    )
    benefit_mean = float(np.mean(benefits))
    benefit_std = float(np.std(benefits))
    denom = abs(benefit_mean) + lambda_scale * benefit_std
    trust = abs(benefit_mean) / denom if denom > 0 else 0.0
    trust = float(np.clip(trust, 0.0, 1.0))

    benefit_cv = (
        benefit_std / abs(benefit_mean) if abs(benefit_mean) > 1e-12 else float("inf")
    )

    mpc_decision = 1.0 if benefit_mean > 0 else 0.0

    return {
        "trust": trust,
        "benefit_mean": benefit_mean,
        "benefit_std": benefit_std,
        "benefit_cv": benefit_cv,
        "mpc_decision": mpc_decision,
        "decision_confident": trust > 0.5,
    }


# --------------------------------------------------------------------------- #
# Calibration-based trust (model-vs-reality)
# --------------------------------------------------------------------------- #
class CalibrationTrust:
    """Online trust based on MPC's predictive calibration.

    Tracks the running error between MPC's horizon prediction and the observed
    tumor volume. When the model's predictions match reality, trust is high and
    control stays with MPC. When the model diverges (e.g. because the surrogate
    ignores resistance, PK mismatch, parameter drift), trust drops and control
    should be handed to the data-driven RL policy.

    trust = 1 / (1 + lambda * relative_mae)

    where relative_mae is the normalized mean absolute error of the model's
    horizon predictions over a trailing window.
    """

    def __init__(self, window: int = 14, lambda_cal: float = 5.0):
        self.window = window
        self.lambda_cal = lambda_cal
        self._preds: list[float] = []
        self._actuals: list[float] = []

    def record(self, predicted: float, actual: float) -> None:
        self._preds.append(float(predicted))
        self._actuals.append(float(actual))
        if len(self._preds) > self.window:
            self._preds.pop(0)
            self._actuals.pop(0)

    def trust(self) -> float:
        if len(self._preds) < 3:
            return 1.0  # not enough evidence -> trust model initially
        preds = np.asarray(self._preds)
        actuals = np.asarray(self._actuals)
        rel_err = np.mean(np.abs(preds - actuals) / np.maximum(actuals, 1e-6))
        return float(np.clip(1.0 / (1.0 + self.lambda_cal * rel_err), 0.0, 1.0))

    @property
    def n_observations(self) -> int:
        return len(self._preds)


if __name__ == "__main__":
    ctrl = RobustMPCController(seed=0)
    for rho in (0.005, 0.02, 0.05):
        sig = compute_trust_signal(ctrl, 5000.0, rho, 0.013, 600.0, step=100)
        print(f"rho={rho}: trust={sig['trust']:.3f} "
              f"benefit={sig['benefit_mean']:+.4f}±{sig['benefit_std']:.4f} "
              f"decision={sig['mpc_decision']}")