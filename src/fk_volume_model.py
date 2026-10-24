"""Shrinkage and predictive log-volume components for Track B.

The caller must construct validation cases from outer-training patients only. These helpers do not inspect
any future target unless that target is explicitly supplied as a training-only validation case.
"""
from __future__ import annotations

import math
from collections import Counter

import numpy as np
from scipy.optimize import nnls


def _finite_scalar(value, name: str) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def shrink_log_growth_slope(patient_slope: float, population_slope: float, *,
                            observation_gap_days: float, time_constant_days: float) -> float:
    """Shrink the patient log-volume slope toward the training-population slope.

    The patient weight is 1-exp(-gap/tau), so short observation gaps borrow more from the population.
    Negative patient slopes remain valid and can forecast shrinkage.
    """
    patient_slope = _finite_scalar(patient_slope, "patient_slope")
    population_slope = _finite_scalar(population_slope, "population_slope")
    gap = _finite_scalar(observation_gap_days, "observation_gap_days")
    tau = _finite_scalar(time_constant_days, "time_constant_days")
    if gap <= 0 or tau <= 0:
        raise ValueError("observation gap and time constant must be positive")
    patient_weight = -math.expm1(-gap / tau)
    return float(population_slope + patient_weight * (patient_slope - population_slope))


def select_shrinkage_time_constant(validation_cases, *, candidates_days=(30, 90, 180, 365)) -> dict:
    """Select tau by log-volume MSE on caller-supplied outer-training validation forecasts.

    Each case must contain patient_slope, population_slope, observation_gap_days, horizon_days, last_volume,
    and target_volume. Optional patient_id makes origins within a patient share one total weight. Cases without
    patient_id are treated as independent patients. The caller must keep all cases inside outer training.
    """
    cases = list(validation_cases)
    candidates = tuple(float(value) for value in candidates_days)
    if not cases:
        raise ValueError("validation_cases must be nonempty")
    if not candidates or any(not np.isfinite(value) or value <= 0 for value in candidates):
        raise ValueError("candidate time constants must be positive and finite")
    if len(set(candidates)) != len(candidates):
        raise ValueError("candidate time constants must be unique")

    required = ("patient_slope", "population_slope", "observation_gap_days", "horizon_days",
                "last_volume", "target_volume")
    normalized = []
    for index, case in enumerate(cases):
        values = {key: _finite_scalar(case[key], key) for key in required}
        if values["observation_gap_days"] <= 0 or values["horizon_days"] <= 0:
            raise ValueError("validation observation gaps and horizons must be positive")
        if values["last_volume"] <= 0 or values["target_volume"] <= 0:
            raise ValueError("validation volumes must be positive")
        patient_id = str(case.get("patient_id", f"case-{index}")).strip()
        if not patient_id:
            raise ValueError("validation patient IDs must be nonempty")
        values["patient_id"] = patient_id
        normalized.append(values)

    errors_by_candidate = {}
    for tau in candidates:
        squared_errors_by_patient = {}
        for case in normalized:
            slope = shrink_log_growth_slope(
                case["patient_slope"], case["population_slope"],
                observation_gap_days=case["observation_gap_days"], time_constant_days=tau,
            )
            predicted_log_volume = math.log(case["last_volume"]) + slope * case["horizon_days"]
            squared_errors_by_patient.setdefault(case["patient_id"], []).append(
                (predicted_log_volume - math.log(case["target_volume"])) ** 2
            )
        patient_mean_errors = [float(np.mean(values)) for values in squared_errors_by_patient.values()]
        errors_by_candidate[tau] = float(np.mean(patient_mean_errors))
    selected = min(candidates, key=lambda tau: (errors_by_candidate[tau], tau))
    return {"time_constant_days": float(selected), "candidate_mse_log_volume": errors_by_candidate,
            "n_validation_forecasts": len(normalized),
            "n_validation_patients": len({case["patient_id"] for case in normalized}),
            "selection_weighting": "equal patient weight; forecast origins averaged within patient",
            "selection_scope": "caller-supplied outer-training validation cases"}


def fit_log_volume_variance(horizons_days, residuals_log_volume, *, patient_ids=None) -> dict:
    """Fit nonnegative linear predictive variance a+b*h with equal total weight per training patient."""
    horizons = np.asarray(horizons_days, dtype=float)
    residuals = np.asarray(residuals_log_volume, dtype=float)
    if (horizons.ndim != 1 or residuals.ndim != 1 or len(horizons) != len(residuals)
            or len(horizons) < 2 or not np.isfinite(horizons).all() or not np.isfinite(residuals).all()):
        raise ValueError("horizons and residuals must be matching finite 1D arrays with at least two rows")
    if np.any(horizons < 0) or len(np.unique(horizons)) < 2:
        raise ValueError("variance fitting needs at least two distinct nonnegative horizons")
    if patient_ids is None:
        ids = [f"case-{index}" for index in range(len(horizons))]
    else:
        if isinstance(patient_ids, (str, bytes)) or len(patient_ids) != len(horizons):
            raise ValueError("patient_ids must match the number of training residuals")
        ids = [str(value).strip() for value in patient_ids]
        if any(not value for value in ids):
            raise ValueError("patient IDs must be nonempty")
    counts = Counter(ids)
    weights = np.asarray([1.0 / math.sqrt(counts[patient_id]) for patient_id in ids])
    design = np.column_stack((np.ones(len(horizons)), horizons))
    coefficients, residual_norm = nnls(design * weights[:, None], residuals**2 * weights)
    return {"intercept": float(coefficients[0]), "slope_per_day": float(coefficients[1]),
            "training_sse_variance": float(residual_norm**2), "n_training_residuals": int(len(horizons)),
            "n_training_patients": len(counts),
            "variance_weighting": "equal patient weight; forecast origins averaged within patient"}


def sample_positive_volumes(last_volume: float, *, horizon_days: float, log_growth_slope: float,
                            variance_fit: dict, n_draws: int = 200, seed: int = 0) -> np.ndarray:
    """Draw positive future volumes from a lognormal model with variance a+b*h."""
    last_volume = _finite_scalar(last_volume, "last_volume")
    horizon = _finite_scalar(horizon_days, "horizon_days")
    slope = _finite_scalar(log_growth_slope, "log_growth_slope")
    if last_volume <= 0 or horizon < 0:
        raise ValueError("last_volume must be positive and horizon_days nonnegative")
    if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 1:
        raise ValueError("n_draws must be a positive integer")
    intercept = _finite_scalar(variance_fit["intercept"], "variance intercept")
    variance_slope = _finite_scalar(variance_fit["slope_per_day"], "variance slope")
    if intercept < 0 or variance_slope < 0:
        raise ValueError("predictive variance coefficients must be nonnegative")
    mean_log_volume = math.log(last_volume) + slope * horizon
    variance = intercept + variance_slope * horizon
    draws = np.exp(np.random.default_rng(seed).normal(mean_log_volume, math.sqrt(variance), int(n_draws)))
    if not np.isfinite(draws).all() or np.any(draws <= 0):
        raise ValueError("predictive log-volume draws overflowed or were nonpositive")
    return draws
