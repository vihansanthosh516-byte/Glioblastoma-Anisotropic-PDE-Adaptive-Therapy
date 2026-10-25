#!/usr/bin/env python3
"""Residual hybrid, transition-count forecast, and mass calibration.

The model couples an image-guided screened-diffusion ranking prior to a
ridge-regularized voxel residual.  A stability gate controls the residual
weight.  A separate Gamma-Poisson rate shrinker predicts the total number of
annotated additions from the previous interval; every spatial comparator can
then use the same expected-count constraint.

This file contains model primitives only.  It does not load patient data,
choose folds, tune parameters, or run a cohort analysis.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


@dataclass(frozen=True)
class VoxelSample:
    """One forecast origin, with future labels isolated in ``target``."""

    patient_id: str
    features: np.ndarray
    log_pde: np.ndarray
    gate: np.ndarray
    target: np.ndarray
    prior_addition_count: int
    prior_candidate_voxels: int
    prior_interval_days: float
    forecast_interval_days: float

    def validate(self) -> None:
        x = np.asarray(self.features)
        pde = np.asarray(self.log_pde).reshape(-1)
        gate = np.asarray(self.gate).reshape(-1)
        target = np.asarray(self.target).reshape(-1)
        if x.ndim != 2 or x.shape[0] != pde.size or pde.size != gate.size or gate.size != target.size:
            raise ValueError("sample features, PDE score, gate, and target have inconsistent shapes")
        if x.shape[0] == 0 or x.shape[1] == 0:
            raise ValueError("sample has no candidate voxels or features")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(pde)) or not np.all(np.isfinite(gate)):
            raise ValueError("sample contains non-finite predictors")
        if np.any((gate < 0) | (gate > 1)):
            raise ValueError("stability gate values must lie in [0, 1]")
        if np.any((target != 0) & (target != 1)):
            raise ValueError("target must be binary")
        if self.prior_addition_count < 0:
            raise ValueError("prior addition count must be nonnegative")
        if self.prior_candidate_voxels <= 0:
            raise ValueError("prior candidate region must contain voxels")
        if self.prior_addition_count > self.prior_candidate_voxels:
            raise ValueError("prior additions cannot exceed the prior candidate region")
        if self.prior_interval_days <= 0 or self.forecast_interval_days <= 0:
            raise ValueError("both scan intervals must be positive")


@dataclass(frozen=True)
class ResidualModel:
    """Fitted score ``beta * PDE + gate * residual(features)``."""

    feature_mean: np.ndarray
    feature_scale: np.ndarray
    parameters: np.ndarray
    ridge: float
    pde_penalty: float
    use_pde: bool
    use_gate: bool

    def predict_score(
        self,
        features: np.ndarray,
        log_pde: np.ndarray,
        gate: np.ndarray,
        *,
        use_pde: bool | None = None,
        use_gate: bool | None = None,
    ) -> np.ndarray:
        x = np.asarray(features, dtype=np.float64)
        pde = np.asarray(log_pde, dtype=np.float64).reshape(-1)
        gate_values = np.asarray(gate, dtype=np.float64).reshape(-1)
        if x.ndim != 2 or x.shape[1] != self.feature_mean.size or x.shape[0] != pde.size:
            raise ValueError("prediction feature shape does not match fitted model")
        if gate_values.size != pde.size:
            raise ValueError("prediction gate shape does not match PDE score")
        x_scaled = np.clip((x - self.feature_mean) / self.feature_scale, -10.0, 10.0)
        residual = self.parameters[1] + x_scaled @ self.parameters[2:]
        use_pde = self.use_pde if use_pde is None else use_pde
        use_gate = self.use_gate if use_gate is None else use_gate
        applied_gate = gate_values if use_gate else np.ones_like(gate_values)
        pde_component = self.parameters[0] * pde if use_pde else np.zeros_like(pde)
        return pde_component + applied_gate * residual


@dataclass(frozen=True)
class CountPrior:
    """Cohort addition intensity and pseudo-exposure in voxel-days."""

    cohort_rate_per_voxel_day: float
    prior_exposure_voxel_days: float


@dataclass(frozen=True)
class TreeBaseline:
    """Non-PDE nonlinear benchmark fitted to the same voxel features."""

    classifier: Any
    max_leaf_nodes: int

    def predict_score(self, features: np.ndarray) -> np.ndarray:
        probability = np.asarray(self.classifier.predict_proba(features), dtype=np.float64)[:, 1]
        probability = np.clip(probability, 1e-8, 1.0 - 1e-8)
        return np.log(probability) - np.log1p(-probability)


def _balanced_rows(
    samples: Sequence[VoxelSample],
    max_per_class: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Deterministically sample each patient/interval with equal class weight."""
    if max_per_class < 1:
        raise ValueError("max_per_class must be positive")
    group_counts: dict[str, int] = {}
    for sample in samples:
        sample.validate()
        group_counts[sample.patient_id] = group_counts.get(sample.patient_id, 0) + 1
    if not group_counts:
        raise ValueError("no training samples supplied")

    x_parts: list[np.ndarray] = []
    pde_parts: list[np.ndarray] = []
    gate_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    weight_parts: list[np.ndarray] = []
    for index, sample in enumerate(samples):
        y = np.asarray(sample.target, dtype=np.uint8).reshape(-1)
        positive = np.flatnonzero(y == 1)
        negative = np.flatnonzero(y == 0)
        if positive.size == 0 or negative.size == 0:
            continue
        rng = np.random.default_rng(seed + index)
        n_pos = min(max_per_class, positive.size)
        n_neg = min(max_per_class, negative.size)
        pos_idx = positive if positive.size <= max_per_class else rng.choice(positive, n_pos, replace=False)
        neg_idx = negative if negative.size <= max_per_class else rng.choice(negative, n_neg, replace=False)
        chosen = np.concatenate((pos_idx, neg_idx))
        y_chosen = y[chosen].astype(np.float64)
        patient_interval_weight = 1.0 / group_counts[sample.patient_id]
        row_weights = np.empty(chosen.size, dtype=np.float64)
        row_weights[:n_pos] = patient_interval_weight * 0.5 / n_pos
        row_weights[n_pos:] = patient_interval_weight * 0.5 / n_neg
        x_parts.append(np.asarray(sample.features[chosen], dtype=np.float64))
        pde_parts.append(np.asarray(sample.log_pde[chosen], dtype=np.float64))
        gate_parts.append(np.asarray(sample.gate[chosen], dtype=np.float64))
        y_parts.append(y_chosen)
        weight_parts.append(row_weights)
    if not x_parts:
        raise ValueError("training fold has no samples with both positive and negative candidates")
    x = np.vstack(x_parts)
    pde = np.concatenate(pde_parts)
    gate = np.concatenate(gate_parts)
    y = np.concatenate(y_parts)
    weights = np.concatenate(weight_parts)
    weights /= weights.sum()
    return x, pde, gate, y, weights


def fit_residual_model(
    samples: Sequence[VoxelSample],
    *,
    ridge: float = 1.0,
    pde_penalty: float = 1e-3,
    max_per_class_per_interval: int = 4000,
    seed: int = 20261024,
    use_pde: bool = True,
    use_gate: bool = True,
) -> ResidualModel:
    """Fit a low-capacity logistic residual with a nonnegative PDE weight.

    Training rows are sampled within patient and interval.  The intercept is
    later mass-calibrated, so class-balanced sampling targets spatial ranking
    and does not supply the predicted positive count.
    """
    if ridge < 0 or pde_penalty < 0:
        raise ValueError("penalties must be nonnegative")
    x, pde, gate, y, weights = _balanced_rows(samples, max_per_class_per_interval, seed)
    pde_fit = pde if use_pde else np.zeros_like(pde)
    gate_fit = gate if use_gate else np.ones_like(gate)
    mean = np.average(x, axis=0, weights=weights)
    variance = np.average((x - mean) ** 2, axis=0, weights=weights)
    scale = np.sqrt(variance)
    scale[scale < 1e-8] = 1.0
    x = np.clip((x - mean) / scale, -10.0, 10.0)
    dimensions = x.shape[1]

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        beta_pde = parameters[0]
        theta0 = parameters[1]
        theta = parameters[2:]
        residual = theta0 + x @ theta
        score = beta_pde * pde_fit + gate_fit * residual
        probability = expit(score)
        loss = np.logaddexp(0.0, score) - y * score
        penalty = 0.5 * ridge * float(np.dot(theta, theta))
        penalty += 0.5 * pde_penalty * beta_pde * beta_pde
        weighted_error = weights * (probability - y)
        gradient = np.empty_like(parameters)
        gradient[0] = float(np.dot(weighted_error, pde_fit)) + pde_penalty * beta_pde
        gradient[1] = float(np.dot(weighted_error, gate_fit))
        gradient[2:] = x.T @ (weighted_error * gate_fit) + ridge * theta
        return float(np.dot(weights, loss) + penalty), gradient

    initial = np.zeros(dimensions + 2, dtype=np.float64)
    fit = minimize(
        objective,
        initial,
        jac=True,
        method="L-BFGS-B",
        bounds=[(0.0, None) if use_pde else (0.0, 0.0)] + [(None, None)] * (dimensions + 1),
        options={"maxiter": 500, "ftol": 1e-10, "gtol": 1e-6},
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"residual model optimization failed: {fit.message}")
    return ResidualModel(
        feature_mean=mean,
        feature_scale=scale,
        parameters=np.asarray(fit.x, dtype=np.float64),
        ridge=float(ridge),
        pde_penalty=float(pde_penalty),
        use_pde=bool(use_pde),
        use_gate=bool(use_gate),
    )


def fit_tree_baseline(
    samples: Sequence[VoxelSample],
    *,
    max_leaf_nodes: int = 15,
    max_iter: int = 100,
    max_per_class_per_interval: int = 1000,
    seed: int = 20261024,
) -> TreeBaseline:
    """Fit a compact histogram-gradient-boosting non-PDE baseline."""
    from sklearn.ensemble import HistGradientBoostingClassifier

    if max_leaf_nodes < 2 or max_iter < 1:
        raise ValueError("invalid gradient-boosting settings")
    x, _, _, y, weights = _balanced_rows(samples, max_per_class_per_interval, seed)
    classifier = HistGradientBoostingClassifier(
        max_iter=max_iter,
        learning_rate=0.08,
        max_leaf_nodes=max_leaf_nodes,
        l2_regularization=1.0,
        early_stopping=False,
        random_state=seed,
    )
    classifier.fit(x, y.astype(np.uint8), sample_weight=weights)
    return TreeBaseline(classifier=classifier, max_leaf_nodes=max_leaf_nodes)


def fit_count_prior(samples: Sequence[VoxelSample]) -> CountPrior:
    """Fit equal-patient additions per candidate-voxel-day from training histories."""
    by_patient: dict[str, list[VoxelSample]] = {}
    for sample in samples:
        sample.validate()
        by_patient.setdefault(sample.patient_id, []).append(sample)
    if not by_patient:
        raise ValueError("no training samples supplied")
    patient_rates = []
    patient_exposures = []
    for patient_samples in by_patient.values():
        counts = sum(sample.prior_addition_count for sample in patient_samples)
        exposure = sum(
            sample.prior_candidate_voxels * sample.prior_interval_days
            for sample in patient_samples
        )
        if exposure > 0:
            patient_rates.append(counts / exposure)
            patient_exposures.append(
                float(
                    np.median(
                        [
                            sample.prior_candidate_voxels * sample.prior_interval_days
                            for sample in patient_samples
                        ]
                    )
                )
            )
    if not patient_rates:
        raise ValueError("training histories contain no positive exposure")
    # The cohort prior gives each patient equal influence.  Its pseudo-exposure
    # equals one median training history exposure, avoiding a hand-picked scale.
    rate = float(np.mean(patient_rates))
    prior_exposure = float(np.median(patient_exposures))
    return CountPrior(max(rate, 0.0), max(prior_exposure, 1e-6))


def predict_addition_count(
    prior: CountPrior,
    prior_addition_count: int,
    prior_candidate_voxels: int,
    prior_interval_days: float,
    forecast_interval_days: float,
    candidate_count: int,
) -> float:
    """Predict expected next-interval additions from past count and durations."""
    if prior_addition_count < 0 or candidate_count < 0:
        raise ValueError("addition and candidate counts must be nonnegative")
    if prior_candidate_voxels <= 0:
        raise ValueError("prior candidate region must contain voxels")
    if prior_interval_days <= 0 or forecast_interval_days <= 0:
        raise ValueError("both scan intervals must be positive")
    posterior_rate = (
        prior_addition_count
        + prior.cohort_rate_per_voxel_day * prior.prior_exposure_voxel_days
    ) / (prior_candidate_voxels * prior_interval_days + prior.prior_exposure_voxel_days)
    predicted = posterior_rate * candidate_count * forecast_interval_days
    return float(np.clip(predicted, 0.0, float(candidate_count)))


def predict_cohort_rate_count(
    prior: CountPrior,
    forecast_interval_days: float,
    candidate_count: int,
) -> float:
    """Cohort-rate-only count baseline, without patient history additions."""
    if forecast_interval_days <= 0 or candidate_count < 0:
        raise ValueError("forecast interval must be positive and candidate count nonnegative")
    return float(
        np.clip(
            prior.cohort_rate_per_voxel_day * candidate_count * forecast_interval_days,
            0.0,
            candidate_count,
        )
    )


def calibrate_probability_mass(scores: np.ndarray, expected_count: float) -> np.ndarray:
    """Apply a shared logistic intercept so probabilities sum to expected count."""
    values = np.asarray(scores, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(values)):
        raise ValueError("scores must be finite")
    if not np.isfinite(expected_count) or expected_count < 0:
        raise ValueError("expected_count must be finite and nonnegative")
    n = values.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    if expected_count <= 1e-12:
        return np.zeros(n, dtype=np.float64)
    if expected_count >= n - 1e-12:
        return np.ones(n, dtype=np.float64)
    centered = values - float(np.median(values))
    bound = max(40.0, float(np.max(np.abs(centered))) + 40.0)
    low, high = -bound, bound
    for _ in range(100):
        intercept = 0.5 * (low + high)
        mass = float(np.sum(expit(intercept + centered)))
        if mass < expected_count:
            low = intercept
        else:
            high = intercept
    return expit(0.5 * (low + high) + centered)


def patient_average_precision(target: np.ndarray, probability: np.ndarray) -> float:
    """Compute per-sample average precision; return NaN for undefined cases."""
    from sklearn.metrics import average_precision_score

    y = np.asarray(target, dtype=np.uint8).reshape(-1)
    p = np.asarray(probability, dtype=np.float64).reshape(-1)
    if y.size != p.size or y.size == 0 or np.unique(y).size < 2:
        return float("nan")
    return float(average_precision_score(y, p))
