#!/usr/bin/env python3
"""Train direction-specific spatial residual heads for signed mask changes.

The addition head combines an externally supplied outward spatial prior with
MRI and prior-transition features. The removal head uses an inward boundary
or retreat prior with the same residual form. Both heads learn spatial ranks;
their probabilities are set by the bounded count law in
``163_signed_transition.py``. Targets remain annotation changes, not verified
biological growth or response.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


@dataclass(frozen=True)
class SignedVoxelSample:
    """One patient interval for one transition direction."""

    patient_id: str
    features: np.ndarray
    spatial_prior_score: np.ndarray
    stability_gate: np.ndarray
    target: np.ndarray

    def validate(self) -> None:
        x = np.asarray(self.features)
        prior = np.asarray(self.spatial_prior_score).reshape(-1)
        gate = np.asarray(self.stability_gate).reshape(-1)
        target = np.asarray(self.target).reshape(-1)
        if not self.patient_id:
            raise ValueError("patient_id must be nonempty")
        if x.ndim != 2 or x.shape[0] != prior.size or prior.size != gate.size or gate.size != target.size:
            raise ValueError("sample features, prior, gate, and target have inconsistent shapes")
        if x.shape[0] == 0 or x.shape[1] == 0:
            raise ValueError("sample has no voxels or features")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(prior)) or not np.all(np.isfinite(gate)):
            raise ValueError("sample predictors must be finite")
        if np.any((gate < 0) | (gate > 1)):
            raise ValueError("stability gate must lie in [0, 1]")
        if np.any((target != 0) & (target != 1)):
            raise ValueError("target must be binary")


def signed_samples_from_maps(
    patient_id: str,
    feature_maps: Mapping[str, np.ndarray],
    feature_names: Sequence[str],
    targets: object,
    addition_prior_map: np.ndarray,
    removal_prior_map: np.ndarray,
    addition_gate_map: np.ndarray,
    removal_gate_map: np.ndarray,
) -> tuple[SignedVoxelSample, SignedVoxelSample]:
    """Extract separate risk-set samples from forecast-origin feature maps.

    ``targets`` must expose addition/removal risk and event masks. Feature,
    prior, and gate maps must come from the forecast origin and earlier scans.
    This function never accepts future labels; callers pass only target masks
    produced after feature construction.
    """
    if not feature_names or any(name not in feature_maps for name in feature_names):
        raise ValueError("feature_names must select existing feature maps")
    targets.validate()
    add_risk = np.asarray(targets.addition_risk, dtype=bool)
    remove_risk = np.asarray(targets.removal_risk, dtype=bool)
    add_target = np.asarray(targets.additions, dtype=bool)
    remove_target = np.asarray(targets.removals, dtype=bool)
    add_prior = np.asarray(addition_prior_map, dtype=np.float64)
    remove_prior = np.asarray(removal_prior_map, dtype=np.float64)
    add_gate = np.asarray(addition_gate_map, dtype=np.float64)
    remove_gate = np.asarray(removal_gate_map, dtype=np.float64)
    shape = add_risk.shape
    maps = [np.asarray(feature_maps[name], dtype=np.float64) for name in feature_names]
    if any(array.shape != shape for array in maps + [
        remove_risk, add_target, remove_target, add_prior, remove_prior, add_gate, remove_gate
    ]):
        raise ValueError("signed target and feature maps must share one grid")
    features = np.stack(maps, axis=-1)
    addition_sample = SignedVoxelSample(
        patient_id=patient_id,
        features=features[add_risk],
        spatial_prior_score=add_prior[add_risk],
        stability_gate=add_gate[add_risk],
        target=add_target[add_risk].astype(np.uint8),
    )
    removal_sample = SignedVoxelSample(
        patient_id=patient_id,
        features=features[remove_risk],
        spatial_prior_score=remove_prior[remove_risk],
        stability_gate=remove_gate[remove_risk],
        target=remove_target[remove_risk].astype(np.uint8),
    )
    addition_sample.validate()
    removal_sample.validate()
    return addition_sample, removal_sample


@dataclass(frozen=True)
class SignedResidualHead:
    """Spatial ranking score: beta*prior + gate*(intercept + scaled residual)."""

    feature_mean: np.ndarray
    feature_scale: np.ndarray
    parameters: np.ndarray
    ridge: float
    prior_penalty: float
    n_training_patients: int

    def predict_score(
        self,
        features: np.ndarray,
        spatial_prior_score: np.ndarray,
        stability_gate: np.ndarray,
    ) -> np.ndarray:
        x = np.asarray(features, dtype=np.float64)
        prior = np.asarray(spatial_prior_score, dtype=np.float64).reshape(-1)
        gate = np.asarray(stability_gate, dtype=np.float64).reshape(-1)
        if x.ndim != 2 or x.shape[1] != self.feature_mean.size or x.shape[0] != prior.size or prior.size != gate.size:
            raise ValueError("prediction features, prior, gate, or fitted dimensions differ")
        if not np.all(np.isfinite(x)) or not np.all(np.isfinite(prior)) or not np.all(np.isfinite(gate)):
            raise ValueError("prediction inputs must be finite")
        if np.any((gate < 0) | (gate > 1)):
            raise ValueError("stability gate must lie in [0, 1]")
        scaled = np.clip((x - self.feature_mean) / self.feature_scale, -10.0, 10.0)
        residual = self.parameters[1] + scaled @ self.parameters[2:]
        return self.parameters[0] * prior + gate * residual


def _balanced_rows(
    samples: Sequence[SignedVoxelSample],
    *,
    max_per_class: int,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    if isinstance(max_per_class, (bool, np.bool_)) or int(max_per_class) != max_per_class or max_per_class < 1:
        raise ValueError("max_per_class must be positive")
    usable_samples: list[SignedVoxelSample] = []
    group_counts: dict[str, int] = {}
    for sample in samples:
        sample.validate()
        target = np.asarray(sample.target, dtype=np.uint8).reshape(-1)
        if np.any(target == 1) and np.any(target == 0):
            usable_samples.append(sample)
            group_counts[sample.patient_id] = group_counts.get(sample.patient_id, 0) + 1
    x_parts: list[np.ndarray] = []
    prior_parts: list[np.ndarray] = []
    gate_parts: list[np.ndarray] = []
    target_parts: list[np.ndarray] = []
    weight_parts: list[np.ndarray] = []
    used_patients: set[str] = set()
    for index, sample in enumerate(usable_samples):
        target = np.asarray(sample.target, dtype=np.uint8).reshape(-1)
        positive = np.flatnonzero(target == 1)
        negative = np.flatnonzero(target == 0)
        rng = np.random.default_rng(seed + index)
        n_positive = min(max_per_class, positive.size)
        n_negative = min(max_per_class, negative.size)
        selected_positive = positive if positive.size <= max_per_class else rng.choice(positive, n_positive, replace=False)
        selected_negative = negative if negative.size <= max_per_class else rng.choice(negative, n_negative, replace=False)
        selected = np.concatenate((selected_positive, selected_negative))
        n_positive_selected = selected_positive.size
        interval_weight = 1.0 / group_counts[sample.patient_id]
        weights = np.concatenate((
            np.full(n_positive_selected, interval_weight * 0.5 / n_positive_selected),
            np.full(selected_negative.size, interval_weight * 0.5 / selected_negative.size),
        ))
        x_parts.append(np.asarray(sample.features[selected], dtype=np.float64))
        prior_parts.append(np.asarray(sample.spatial_prior_score, dtype=np.float64)[selected])
        gate_parts.append(np.asarray(sample.stability_gate, dtype=np.float64)[selected])
        target_parts.append(target[selected].astype(np.float64))
        weight_parts.append(weights)
        used_patients.add(sample.patient_id)
    if len(used_patients) < 2:
        raise ValueError("fewer than two training patients have both transition classes")
    x = np.vstack(x_parts)
    prior = np.concatenate(prior_parts)
    gate = np.concatenate(gate_parts)
    target = np.concatenate(target_parts)
    weights = np.concatenate(weight_parts)
    weights /= weights.sum()
    return x, prior, gate, target, weights, len(used_patients)


def fit_signed_residual_head(
    samples: Sequence[SignedVoxelSample],
    *,
    ridge: float = 1.0,
    prior_penalty: float = 1e-3,
    max_per_class_per_interval: int = 4000,
    seed: int = 20261025,
) -> SignedResidualHead:
    """Fit one direction's patient-weighted logistic residual on training data."""
    if (
        not np.isfinite(ridge)
        or not np.isfinite(prior_penalty)
        or ridge < 0
        or prior_penalty < 0
    ):
        raise ValueError("penalties must be nonnegative")
    x, prior, gate, target, weights, n_patients = _balanced_rows(
        samples, max_per_class=max_per_class_per_interval, seed=seed
    )
    mean = np.average(x, axis=0, weights=weights)
    variance = np.average((x - mean) ** 2, axis=0, weights=weights)
    scale = np.sqrt(variance)
    scale[scale < 1e-8] = 1.0
    scaled = np.clip((x - mean) / scale, -10.0, 10.0)

    def objective(parameters: np.ndarray) -> tuple[float, np.ndarray]:
        prior_weight = parameters[0]
        residual = parameters[1] + scaled @ parameters[2:]
        score = prior_weight * prior + gate * residual
        probability = expit(score)
        loss = np.logaddexp(0.0, score) - target * score
        penalty = 0.5 * ridge * float(np.dot(parameters[2:], parameters[2:]))
        penalty += 0.5 * prior_penalty * prior_weight * prior_weight
        error = weights * (probability - target)
        gradient = np.empty_like(parameters)
        gradient[0] = float(np.dot(error, prior)) + prior_penalty * prior_weight
        gradient[1] = float(np.dot(error, gate))
        gradient[2:] = scaled.T @ (error * gate) + ridge * parameters[2:]
        return float(np.dot(weights, loss) + penalty), gradient

    initial = np.zeros(x.shape[1] + 2, dtype=np.float64)
    fit = minimize(
        objective,
        initial,
        jac=True,
        method="L-BFGS-B",
        bounds=[(0.0, None)] + [(None, None)] * (x.shape[1] + 1),
        options={"maxiter": 500, "ftol": 1e-10, "gtol": 1e-6},
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"signed spatial head fit failed: {fit.message}")
    return SignedResidualHead(
        feature_mean=mean,
        feature_scale=scale,
        parameters=np.asarray(fit.x, dtype=np.float64),
        ridge=float(ridge),
        prior_penalty=float(prior_penalty),
        n_training_patients=n_patients,
    )


@dataclass(frozen=True)
class SignedSpatialHybrid:
    """Separate addition and removal heads; count law is applied downstream."""

    addition_head: SignedResidualHead
    removal_head: SignedResidualHead

    @classmethod
    def fit(
        cls,
        addition_samples: Sequence[SignedVoxelSample],
        removal_samples: Sequence[SignedVoxelSample],
        *,
        ridge: float = 1.0,
        prior_penalty: float = 1e-3,
        max_per_class_per_interval: int = 4000,
        seed: int = 20261025,
    ) -> "SignedSpatialHybrid":
        """Fit both direction-specific heads from training patients only."""
        add_head = fit_signed_residual_head(
            addition_samples,
            ridge=ridge,
            prior_penalty=prior_penalty,
            max_per_class_per_interval=max_per_class_per_interval,
            seed=seed,
        )
        remove_head = fit_signed_residual_head(
            removal_samples,
            ridge=ridge,
            prior_penalty=prior_penalty,
            max_per_class_per_interval=max_per_class_per_interval,
            seed=seed + 1,
        )
        return cls(addition_head=add_head, removal_head=remove_head)


def predict_signed_score_maps(
    model: SignedSpatialHybrid,
    feature_maps: Mapping[str, np.ndarray],
    addition_feature_names: Sequence[str],
    removal_feature_names: Sequence[str],
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    addition_prior_map: np.ndarray,
    removal_prior_map: np.ndarray,
    addition_gate_map: np.ndarray,
    removal_gate_map: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply both fitted spatial heads only inside their forecast-origin risks."""
    add_risk = np.asarray(addition_risk, dtype=bool)
    remove_risk = np.asarray(removal_risk, dtype=bool)
    add_prior = np.asarray(addition_prior_map, dtype=np.float64)
    remove_prior = np.asarray(removal_prior_map, dtype=np.float64)
    add_gate = np.asarray(addition_gate_map, dtype=np.float64)
    remove_gate = np.asarray(removal_gate_map, dtype=np.float64)
    if add_risk.ndim != 3 or any(
        array.shape != add_risk.shape
        for array in (remove_risk, add_prior, remove_prior, add_gate, remove_gate)
    ):
        raise ValueError("risk, prior, and gate maps must share one 3-D grid")
    if np.any(add_risk & remove_risk):
        raise ValueError("addition and removal risk sets must be disjoint")

    def head_score(
        head: SignedResidualHead,
        names: Sequence[str],
        risk: np.ndarray,
        prior: np.ndarray,
        gate: np.ndarray,
    ) -> np.ndarray:
        if not names or any(name not in feature_maps for name in names):
            raise ValueError("feature names must select existing feature maps")
        selected_maps = [np.asarray(feature_maps[name], dtype=np.float64) for name in names]
        if any(array.shape != risk.shape for array in selected_maps):
            raise ValueError("feature maps must share the forecast grid")
        full_score = np.zeros(risk.shape, dtype=np.float64)
        if np.any(risk):
            design = np.column_stack([array[risk] for array in selected_maps])
            full_score[risk] = head.predict_score(design, prior[risk], gate[risk])
        return full_score

    return (
        head_score(model.addition_head, addition_feature_names, add_risk, add_prior, add_gate),
        head_score(model.removal_head, removal_feature_names, remove_risk, remove_prior, remove_gate),
    )
