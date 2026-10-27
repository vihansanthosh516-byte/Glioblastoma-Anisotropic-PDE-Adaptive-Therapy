#!/usr/bin/env python3
"""Exploratory patient-grouped evaluation of signed Track B transitions.

This runner compares a two-head spatial-prior residual against prior-only and
MRI-only arms. A bounded shared-frailty count law constrains all probability
maps. MU-Glioma-Post has been examined in earlier project work, so results are
exploratory even with patient-grouped folds. Endpoints are supplied mask
transitions, not verified biological tumor growth or response.

Count predictions condition on the realized next-scan interval from the
retrospective manifest. They are not pre-schedule forecasts unless the
interval is supplied at forecast origin or replaced by a prespecified horizon.
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Sequence

import numpy as np
from scipy.optimize import brentq
from scipy.ndimage import distance_transform_edt, gaussian_filter
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
from sklearn.model_selection import GroupKFold


ROOT = Path(__file__).resolve().parent.parent
FEATURES_PATH = ROOT / "src/160_growth_zone_features.py"
BASE_RUNNER_PATH = ROOT / "src/162_growth_zone_cv.py"
TRANSITIONS_PATH = ROOT / "src/163_signed_transition.py"
HEADS_PATH = ROOT / "src/164_signed_spatial_heads.py"
DEFAULT_OUTPUT_DIR = ROOT / "output"
FORECAST_HORIZON = {
    "estimand": "retrospective prediction conditional on realized next-scan interval",
    "interval_source": "actual future interval in the chained scan manifest",
    "known_at_forecast_origin": False,
    "deployment_interpretation": "not a pre-schedule forecast; supply the interval at origin or prespecify a fixed horizon",
}


def _logical_input_path(path: Path, repo_root: Path = ROOT) -> str:
    """Describe an input by repo-relative name without exposing external paths."""
    try:
        resolved = Path(path).expanduser().resolve()
        relative = resolved.relative_to(Path(repo_root).resolve())
    except (OSError, RuntimeError, ValueError):
        return "external-input"
    return relative.as_posix() or "."


def _load_module(alias: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(alias, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


BASE = _load_module("signed_cv_base_runner", BASE_RUNNER_PATH)
FEATURES = _load_module("signed_cv_features", FEATURES_PATH)
TRANSITIONS = _load_module("signed_cv_transition_model", TRANSITIONS_PATH)
HEADS = _load_module("signed_cv_spatial_heads", HEADS_PATH)

ADDITION_FEATURES = tuple(FEATURES.FEATURE_NAMES) + (
    "inside_distance_mm",
    "prior_removal_density",
)
REMOVAL_FEATURES = (
    "inside_distance_mm",
    "distance_to_cavity_mm",
    "t1n_z",
    "t1c_z",
    "t2w_z",
    "t2f_z",
    "label2_indicator",
    "prior_addition_density",
    "prior_removal_density",
)


@dataclass(frozen=True)
class RawSignedInterval:
    """Origin-only predictors, future targets, and realized next-scan gap."""

    sample_id: str
    target_interval_id: str
    patient_id: str
    addition_features: np.ndarray
    addition_prior: np.ndarray
    addition_target: np.ndarray
    removal_features: np.ndarray
    removal_prior: np.ndarray
    removal_target: np.ndarray
    history_count: Any
    target_count: Any
    boundary_counts: dict[str, int]


def _interval_features(
    patient_id: str,
    previous: Any,
    current: Any,
    future_labels: np.ndarray,
    future_affine: np.ndarray,
    parameters: Any,
    target_interval_id: str,
    previous_interval_days: float,
    forecast_interval_days: float,
    perturbation_radius_mm: float,
) -> RawSignedInterval:
    spacing = FEATURES.validate_aligned_visits(previous, current)
    current_images, current_labels, current_affine, _ = FEATURES.validate_visit(current, "current")
    if future_labels.shape != current_labels.shape or not np.allclose(
        future_affine, current_affine, rtol=0.0, atol=1e-3
    ):
        raise ValueError(f"{patient_id}: future mask is not aligned to forecast origin")

    visit = FEATURES.build_visit_features(previous, current, parameters)
    previous_labels = np.asarray(previous.labels)
    previous_tumor = FEATURES.annotated_core_mask(previous_labels)
    current_tumor = visit.source
    prior_removals = previous_tumor & ~current_tumor
    sigma = tuple(float(parameters.prior_smoothing_mm / value) for value in spacing)
    prior_removal_density = gaussian_filter(prior_removals.astype(np.float32), sigma=sigma)
    peak = float(np.max(prior_removal_density))
    if peak > 0:
        prior_removal_density /= peak
    inside_distance = distance_transform_edt(current_tumor, sampling=spacing).astype(np.float32)

    maps = dict(visit.maps)
    maps["inside_distance_mm"] = inside_distance
    maps["prior_removal_density"] = prior_removal_density.astype(np.float32, copy=False)
    add_target = TRANSITIONS.make_signed_transition_targets(
        current_labels,
        future_labels,
        visit.candidate,
        removal_domain=visit.source,
    )
    count_target = TRANSITIONS.bounded_observation_from_targets(
        patient_id, add_target, forecast_interval_days
    )
    boundary = TRANSITIONS.make_boundary_sensitivity_bounds(
        current_labels,
        future_labels,
        visit.candidate,
        spacing,
        perturbation_radius_mm=perturbation_radius_mm,
        removal_domain=visit.source,
    )

    add_design = np.column_stack([maps[name][add_target.addition_risk] for name in ADDITION_FEATURES])
    remove_design = np.column_stack([maps[name][add_target.removal_risk] for name in REMOVAL_FEATURES])
    if add_design.shape[0] == 0 or remove_design.shape[0] == 0:
        raise ValueError(f"{patient_id}: empty addition or removal risk set")

    previous_images, previous_labels, _, previous_spacing = FEATURES.validate_visit(previous, "previous")
    previous_support = FEATURES.brain_support(previous_images, previous_labels)
    previous_source, previous_candidate, _ = FEATURES._candidate_region(
        previous_labels,
        previous_support,
        previous_spacing,
        parameters.candidate_radius_mm,
        parameters.block_cavity,
    )
    history_targets = TRANSITIONS.make_signed_transition_targets(
        previous_labels,
        current_labels,
        previous_candidate,
        removal_domain=previous_source,
    )
    history_count = TRANSITIONS.bounded_observation_from_targets(
        patient_id, history_targets, previous_interval_days
    )
    return RawSignedInterval(
        sample_id=target_interval_id,
        target_interval_id=target_interval_id,
        patient_id=patient_id,
        addition_features=add_design.astype(np.float32, copy=False),
        addition_prior=visit.log_pde_field[add_target.addition_risk].astype(np.float32, copy=False),
        addition_target=add_target.additions[add_target.addition_risk].astype(np.uint8),
        removal_features=remove_design.astype(np.float32, copy=False),
        removal_prior=(-inside_distance[add_target.removal_risk] / 10.0).astype(np.float32, copy=False),
        removal_target=add_target.removals[add_target.removal_risk].astype(np.uint8),
        history_count=history_count,
        target_count=count_target,
        boundary_counts={
            "addition_certain": int(boundary.addition_certain.sum()),
            "addition_nominal": int(boundary.nominal.additions.sum()),
            "addition_possible": int(boundary.addition_possible.sum()),
            "removal_certain": int(boundary.removal_certain.sum()),
            "removal_nominal": int(boundary.nominal.removals.sum()),
            "removal_possible": int(boundary.removal_possible.sum()),
        },
    )


def _load_signed_intervals(
    *,
    data_root: Path,
    patient_manifest: Path,
    pair_manifest: Path,
    parameters: Any,
    max_patients: int | None,
    perturbation_radius_mm: float,
) -> list[RawSignedInterval]:
    triplets = BASE._manifest_triplets(patient_manifest, pair_manifest)
    patient_ids = sorted({patient_id for patient_id, _, _ in triplets})
    if max_patients is not None:
        if max_patients < 2:
            raise ValueError("max_patients must be at least 2")
        keep = set(patient_ids[:max_patients])
        triplets = [row for row in triplets if row[0] in keep]
    intervals: list[RawSignedInterval] = []
    for index, (patient_id, first, second) in enumerate(triplets, start=1):
        tp1, tp2, tp3 = int(first["tp_in"]), int(first["tp_out"]), int(second["tp_out"])
        previous = BASE._load_visit(data_root, patient_id, tp1)
        current = BASE._load_visit(data_root, patient_id, tp2)
        future_labels, future_affine = BASE._load_mask(data_root, patient_id, tp3)
        interval = _interval_features(
            patient_id,
            previous,
            current,
            future_labels,
            future_affine,
            parameters,
            f"{patient_id}:{tp2}-{tp3}",
            float(first["dt_days"]),
            float(second["dt_days"]),
            perturbation_radius_mm,
        )
        intervals.append(interval)
        print(f"loaded signed interval {index}/{len(triplets)}: {patient_id}", flush=True)
    return intervals


def _head_samples(
    intervals: Sequence[RawSignedInterval],
    *,
    arm: str,
) -> tuple[list[Any], list[Any]]:
    additions = []
    removals = []
    for row in intervals:
        add_prior = row.addition_prior if arm == "hybrid" else np.zeros_like(row.addition_prior)
        remove_prior = row.removal_prior if arm == "hybrid" else np.zeros_like(row.removal_prior)
        additions.append(
            HEADS.SignedVoxelSample(
                patient_id=row.patient_id,
                features=row.addition_features,
                spatial_prior_score=add_prior,
                stability_gate=np.ones(row.addition_target.size, dtype=np.float32),
                target=row.addition_target,
            )
        )
        removals.append(
            HEADS.SignedVoxelSample(
                patient_id=row.patient_id,
                features=row.removal_features,
                spatial_prior_score=remove_prior,
                stability_gate=np.ones(row.removal_target.size, dtype=np.float32),
                target=row.removal_target,
            )
        )
    return additions, removals


def _unique_count_observations(intervals: Sequence[RawSignedInterval]) -> list[Any]:
    by_interval: dict[str, Any] = {}
    for row in intervals:
        previous = by_interval.get(row.target_interval_id)
        if previous is not None and previous != row.target_count:
            raise ValueError("duplicate patient transition has conflicting count observations")
        by_interval[row.target_interval_id] = row.target_count
    return list(by_interval.values())


def _forecast_count_pair(
    count_model: Any,
    row: RawSignedInterval,
    *,
    addition_rate_multiplier: float,
    removal_rate_multiplier: float,
) -> tuple[Any, Any]:
    """Forecast counts from origin history, risk-set sizes, and realized horizon.

    The realized future interval length is used as a conditional horizon.
    Future addition and removal counts are not passed to the count model.
    """
    arguments = {
        "prior_observations": [row.history_count],
        "future_addition_risk_voxels": row.addition_target.size,
        "future_removal_risk_voxels": row.removal_target.size,
        "future_interval_days": row.target_count.interval_days,
    }
    uncalibrated = count_model.forecast(**arguments)
    calibrated = count_model.forecast(
        **arguments,
        future_addition_rate_multiplier=addition_rate_multiplier,
        future_removal_rate_multiplier=removal_rate_multiplier,
    )
    return uncalibrated, calibrated


def _crossfit_count_rate_calibration(
    training_intervals: Sequence[RawSignedInterval],
    *,
    seed: int,
    max_inner_splits: int = 5,
) -> dict[str, Any]:
    """Fit future-rate multipliers from inner-fold predictions only.

    The count law is cross-fitted by patient inside the outer training set.
    Separate monotone multipliers then match the total out-of-fold expected
    additions and removals to their observed training totals. This calibrates
    aggregate count bias; it does not guarantee patient-level or interval-level
    calibration.
    """
    patients = np.asarray([row.patient_id for row in training_intervals])
    unique_patients = np.unique(patients)
    if unique_patients.size < 2:
        raise ValueError("count calibration needs at least two training patients")
    n_inner_splits = min(max_inner_splits, int(unique_patients.size))
    splitter = GroupKFold(n_splits=n_inner_splits)
    oof_forecasts: list[Any] = []
    oof_targets: list[Any] = []
    for inner_fold, (inner_train_idx, inner_valid_idx) in enumerate(
        splitter.split(np.zeros(len(training_intervals)), groups=patients), start=1
    ):
        inner_train = [training_intervals[index] for index in inner_train_idx]
        inner_valid = [training_intervals[index] for index in inner_valid_idx]
        train_ids = {row.patient_id for row in inner_train}
        valid_ids = {row.patient_id for row in inner_valid}
        if train_ids & valid_ids:
            raise RuntimeError("patient leakage detected in inner count calibration")
        count_model = TRANSITIONS.SharedFrailtyBinomialPrior.fit(
            _unique_count_observations(inner_train)
        )
        for row in inner_valid:
            forecast = count_model.forecast(
                prior_observations=[row.history_count],
                future_addition_risk_voxels=row.addition_target.size,
                future_removal_risk_voxels=row.removal_target.size,
                future_interval_days=row.target_count.interval_days,
            )
            oof_forecasts.append(forecast)
            oof_targets.append(row.target_count)

    if len(oof_forecasts) != len(training_intervals):
        raise RuntimeError("inner count calibration did not produce one prediction per interval")

    def calibrate_direction(direction: str) -> dict[str, Any]:
        if direction == "addition":
            target_total = float(sum(row.addition_count for row in oof_targets))
            rate_name = "addition_rate_per_voxel_day"
            risk_name = "future_addition_risk_voxels"
        else:
            target_total = float(sum(row.removal_count for row in oof_targets))
            rate_name = "removal_rate_per_voxel_day"
            risk_name = "future_removal_risk_voxels"

        def expected_total(multiplier: float) -> float:
            total = 0.0
            for forecast, observation in zip(oof_forecasts, oof_targets):
                rate = float(getattr(forecast, rate_name))
                risk = int(getattr(forecast, risk_name))
                interval_days = float(forecast.future_interval_days)
                nodes = np.asarray(forecast.frailty_nodes, dtype=np.float64)
                weights = np.asarray(forecast.posterior_weights, dtype=np.float64)
                probabilities = -np.expm1(-nodes * rate * multiplier * interval_days)
                total += risk * float(np.dot(weights, probabilities))
            return total

        uncalibrated_total = expected_total(1.0)
        if target_total <= 0:
            multiplier = 0.0
            boundary = False
        elif uncalibrated_total <= 0:
            multiplier = 1.0
            boundary = True
        else:
            upper = 1.0
            while expected_total(upper) < target_total and upper < 1e8:
                upper *= 2.0
            if expected_total(upper) < target_total:
                multiplier = upper
                boundary = True
            else:
                multiplier = float(
                    brentq(
                        lambda value: expected_total(value) - target_total,
                        0.0,
                        upper,
                        xtol=1e-10,
                        rtol=1e-10,
                        maxiter=200,
                    )
                )
                boundary = False
        calibrated_total = expected_total(multiplier)
        return {
            "rate_multiplier": multiplier,
            "observed_oof_total": int(target_total),
            "uncalibrated_oof_expected_total": uncalibrated_total,
            "calibrated_oof_expected_total": calibrated_total,
            "at_boundary": boundary,
        }

    return {
        "method": "direction-specific aggregate moment calibration on inner patient-grouped out-of-fold predictions",
        "n_inner_splits": n_inner_splits,
        "n_training_patients": int(unique_patients.size),
        "n_oof_intervals": len(oof_forecasts),
        "addition": calibrate_direction("addition"),
        "removal": calibrate_direction("removal"),
    }


def _rank_probabilities(
    addition_score: np.ndarray,
    removal_score: np.ndarray,
    count_forecast: Any,
) -> tuple[np.ndarray, np.ndarray]:
    n_add = int(count_forecast.future_addition_risk_voxels)
    n_remove = int(count_forecast.future_removal_risk_voxels)
    total = n_add + n_remove
    add_risk = np.zeros(total, dtype=bool)
    remove_risk = np.zeros(total, dtype=bool)
    add_risk[:n_add] = True
    remove_risk[n_add:] = True
    add_scores = np.zeros(total, dtype=np.float64)
    remove_scores = np.zeros(total, dtype=np.float64)
    add_scores[:n_add] = np.asarray(addition_score, dtype=np.float64)
    remove_scores[n_add:] = np.asarray(removal_score, dtype=np.float64)
    forecast = TRANSITIONS.rank_constrained_signed_forecast(
        add_scores.reshape(total, 1, 1),
        remove_scores.reshape(total, 1, 1),
        add_risk.reshape(total, 1, 1),
        remove_risk.reshape(total, 1, 1),
        count_forecast,
    )
    return (
        forecast.addition_probability[:n_add, 0, 0],
        forecast.removal_probability[n_add:, 0, 0],
    )


def _equal_weight_probability_mixture(
    hybrid_probability: np.ndarray,
    prior_probability: np.ndarray,
) -> np.ndarray:
    """Return a fixed equal-weight linear pool of two spatial marginals.

    Both component maps use the same count law, so their mixture preserves the
    expected event count while hedging between learned and prior-only location
    forecasts. The 1/2 weight is fixed, not fitted on the evaluation patients.
    """
    hybrid = np.asarray(hybrid_probability, dtype=np.float64)
    prior = np.asarray(prior_probability, dtype=np.float64)
    if hybrid.shape != prior.shape or hybrid.ndim != 1:
        raise ValueError("mixture inputs must be one-dimensional maps with matching shape")
    if (
        not np.all(np.isfinite(hybrid))
        or not np.all(np.isfinite(prior))
        or np.any((hybrid < 0) | (hybrid > 1))
        or np.any((prior < 0) | (prior > 1))
    ):
        raise ValueError("mixture inputs must contain finite probabilities in [0, 1]")
    return 0.5 * hybrid + 0.5 * prior


def _metric_row(target: np.ndarray, score: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    y = np.asarray(target, dtype=np.uint8)
    if y.size == 0 or y.min() == y.max():
        return {"average_precision": float("nan"), "auroc": float("nan"), "brier": float("nan")}
    return {
        "average_precision": float(average_precision_score(y, score)),
        "auroc": float(roc_auc_score(y, score)),
        "brier": float(brier_score_loss(y, probability)),
    }


def _aggregate_by_patient(records: Sequence[dict[str, Any]], arm: str, outcome: str, metric: str) -> dict[str, float]:
    grouped: dict[str, list[float]] = {}
    key = f"{arm}_{outcome}_{metric}"
    for record in records:
        value = float(record[key])
        if np.isfinite(value):
            grouped.setdefault(record["patient_id"], []).append(value)
    return {patient_id: float(np.mean(values)) for patient_id, values in grouped.items()}


def _bootstrap_difference(
    left: dict[str, float],
    right: dict[str, float],
    *,
    seed: int,
    draws: int = 2000,
) -> dict[str, float]:
    patients = sorted(set(left) & set(right))
    if len(patients) < 2:
        return {"n_patients": len(patients), "mean_difference": float("nan"), "ci95_low": float("nan"), "ci95_high": float("nan")}
    differences = np.asarray([left[patient] - right[patient] for patient in patients], dtype=np.float64)
    rng = np.random.default_rng(seed)
    sampled = rng.integers(0, len(patients), size=(draws, len(patients)))
    means = np.mean(differences[sampled], axis=1)
    return {
        "n_patients": len(patients),
        "mean_difference": float(np.mean(differences)),
        "ci95_low": float(np.quantile(means, 0.025)),
        "ci95_high": float(np.quantile(means, 0.975)),
    }


def evaluate(
    intervals: Sequence[RawSignedInterval],
    *,
    n_splits: int = 5,
    seed: int = 20261025,
) -> dict[str, Any]:
    patients = np.asarray([row.patient_id for row in intervals])
    unique_patients = np.unique(patients)
    if n_splits < 2 or n_splits > unique_patients.size:
        raise ValueError("n_splits must be between 2 and the number of patients")
    splitter = GroupKFold(n_splits=n_splits)
    records: list[dict[str, Any]] = []
    fold_summaries = []
    for fold, (train_idx, test_idx) in enumerate(splitter.split(np.zeros(len(intervals)), groups=patients), start=1):
        train = [intervals[index] for index in train_idx]
        test = [intervals[index] for index in test_idx]
        train_ids = {row.patient_id for row in train}
        test_ids = {row.patient_id for row in test}
        if train_ids & test_ids:
            raise RuntimeError("patient leakage detected across outer fold")
        count_calibration = _crossfit_count_rate_calibration(
            train,
            seed=seed + 1000 * fold,
        )
        train_counts = _unique_count_observations(train)
        count_model = TRANSITIONS.SharedFrailtyBinomialPrior.fit(train_counts)
        hybrid_add, hybrid_remove = _head_samples(train, arm="hybrid")
        mri_add, mri_remove = _head_samples(train, arm="mri")
        hybrid = HEADS.SignedSpatialHybrid.fit(
            hybrid_add,
            hybrid_remove,
            ridge=1.0,
            prior_penalty=1e-3,
            max_per_class_per_interval=2000,
            seed=seed + fold,
        )
        mri_model = HEADS.SignedSpatialHybrid.fit(
            mri_add,
            mri_remove,
            ridge=1.0,
            prior_penalty=1e-3,
            max_per_class_per_interval=2000,
            seed=seed + 100 + fold,
        )
        for row in test:
            uncalibrated_count_forecast, count_forecast = _forecast_count_pair(
                count_model,
                row,
                addition_rate_multiplier=count_calibration["addition"]["rate_multiplier"],
                removal_rate_multiplier=count_calibration["removal"]["rate_multiplier"],
            )
            hybrid_add_score = hybrid.addition_head.predict_score(
                row.addition_features,
                row.addition_prior,
                np.ones(row.addition_target.size, dtype=np.float32),
            )
            hybrid_remove_score = hybrid.removal_head.predict_score(
                row.removal_features,
                row.removal_prior,
                np.ones(row.removal_target.size, dtype=np.float32),
            )
            mri_add_score = mri_model.addition_head.predict_score(
                row.addition_features,
                np.zeros(row.addition_target.size, dtype=np.float32),
                np.ones(row.addition_target.size, dtype=np.float32),
            )
            mri_remove_score = mri_model.removal_head.predict_score(
                row.removal_features,
                np.zeros(row.removal_target.size, dtype=np.float32),
                np.ones(row.removal_target.size, dtype=np.float32),
            )
            score_arms = {
                "hybrid": (hybrid_add_score, hybrid_remove_score),
                "mri_only": (mri_add_score, mri_remove_score),
                "prior_only": (row.addition_prior.astype(np.float64), row.removal_prior.astype(np.float64)),
            }
            record: dict[str, Any] = {
                "patient_id": row.patient_id,
                "sample_id": row.sample_id,
                "target_interval_id": row.target_interval_id,
                "fold": fold,
                "addition_count_observed": row.target_count.addition_count,
                "removal_count_observed": row.target_count.removal_count,
                "addition_risk_voxels": int(row.addition_target.size),
                "removal_risk_voxels": int(row.removal_target.size),
                "future_interval_days": float(row.target_count.interval_days),
                "history_addition_count": int(row.history_count.addition_count),
                "history_removal_count": int(row.history_count.removal_count),
                "addition_count_predicted": count_forecast.expected_additions,
                "removal_count_predicted": count_forecast.expected_removals,
                "addition_count_predicted_uncalibrated": uncalibrated_count_forecast.expected_additions,
                "removal_count_predicted_uncalibrated": uncalibrated_count_forecast.expected_removals,
                "count_joint_log_score": count_forecast.log_probability(
                    row.target_count.addition_count, row.target_count.removal_count
                ),
                "count_joint_log_score_uncalibrated": uncalibrated_count_forecast.log_probability(
                    row.target_count.addition_count, row.target_count.removal_count
                ),
                "calibrated_joint_count_log_score": count_forecast.log_probability(
                    row.target_count.addition_count, row.target_count.removal_count
                ),
                "uncalibrated_joint_count_log_score": uncalibrated_count_forecast.log_probability(
                    row.target_count.addition_count, row.target_count.removal_count
                ),
                "calibrated_addition_count_absolute_error": abs(
                    count_forecast.expected_additions - row.target_count.addition_count
                ),
                "uncalibrated_addition_count_absolute_error": abs(
                    uncalibrated_count_forecast.expected_additions - row.target_count.addition_count
                ),
                "calibrated_removal_count_absolute_error": abs(
                    count_forecast.expected_removals - row.target_count.removal_count
                ),
                "uncalibrated_removal_count_absolute_error": abs(
                    uncalibrated_count_forecast.expected_removals - row.target_count.removal_count
                ),
                "addition_future_rate_multiplier": count_calibration["addition"]["rate_multiplier"],
                "removal_future_rate_multiplier": count_calibration["removal"]["rate_multiplier"],
                "boundary_counts": row.boundary_counts,
            }
            probability_arms: dict[str, tuple[np.ndarray, np.ndarray]] = {}
            for arm, (add_score, remove_score) in score_arms.items():
                add_probability, remove_probability = _rank_probabilities(
                    add_score, remove_score, count_forecast
                )
                probability_arms[arm] = (add_probability, remove_probability)
                for outcome, target, score, probability in (
                    ("addition", row.addition_target, add_score, add_probability),
                    ("removal", row.removal_target, remove_score, remove_probability),
                ):
                    for metric, value in _metric_row(target, score, probability).items():
                        record[f"{arm}_{outcome}_{metric}"] = value
            mixture_add = _equal_weight_probability_mixture(
                probability_arms["hybrid"][0], probability_arms["prior_only"][0]
            )
            mixture_remove = _equal_weight_probability_mixture(
                probability_arms["hybrid"][1], probability_arms["prior_only"][1]
            )
            for outcome, target, probability in (
                ("addition", row.addition_target, mixture_add),
                ("removal", row.removal_target, mixture_remove),
            ):
                for metric, value in _metric_row(target, probability, probability).items():
                    record[f"prior_anchored_mixture_{outcome}_{metric}"] = value
            records.append(record)
        fold_summaries.append(
            {
                "fold": fold,
                "n_train_patients": len(train_ids),
                "n_test_patients": len(test_ids),
                "train_patients": sorted(train_ids),
                "test_patients": sorted(test_ids),
                "count_model_training_patients": count_model.n_patients,
                "frailty_shape": count_model.frailty_shape,
                "frailty_shape_at_boundary": count_model.shape_at_fit_boundary,
                "count_rate_calibration": count_calibration,
            }
        )
        print(f"evaluated patient fold {fold}/{n_splits}", flush=True)

    contrasts = {}
    for outcome in ("addition", "removal"):
        for metric in ("average_precision", "brier"):
            hybrid_values = _aggregate_by_patient(records, "hybrid", outcome, metric)
            for comparator in ("prior_only", "mri_only"):
                comparator_values = _aggregate_by_patient(records, comparator, outcome, metric)
                contrast = _bootstrap_difference(
                    hybrid_values,
                    comparator_values,
                    seed=seed + len(contrasts),
                )
                if metric == "brier":
                    contrast = {**contrast, "direction": "negative favors hybrid"}
                else:
                    contrast = {**contrast, "direction": "positive favors hybrid"}
                contrasts[f"hybrid_minus_{comparator}_{outcome}_{metric}"] = contrast
            mixture_values = _aggregate_by_patient(
                records, "prior_anchored_mixture", outcome, metric
            )
            for comparator in ("hybrid", "prior_only", "mri_only"):
                comparator_values = _aggregate_by_patient(records, comparator, outcome, metric)
                contrast = _bootstrap_difference(
                    mixture_values,
                    comparator_values,
                    seed=seed + len(contrasts),
                )
                direction = "negative favors mixture" if metric == "brier" else "positive favors mixture"
                contrasts[f"prior_anchored_mixture_minus_{comparator}_{outcome}_{metric}"] = {
                    **contrast,
                    "direction": direction,
                }
        for metric in ("count_absolute_error",):
            calibrated_values = _aggregate_by_patient(records, "calibrated", outcome, metric)
            uncalibrated_values = _aggregate_by_patient(records, "uncalibrated", outcome, metric)
            contrast = _bootstrap_difference(
                calibrated_values,
                uncalibrated_values,
                seed=seed + len(contrasts),
            )
            contrasts[f"calibrated_minus_uncalibrated_{outcome}_{metric}"] = {
                **contrast,
                "direction": "negative favors calibrated",
            }
    calibrated_log_scores = _aggregate_by_patient(records, "calibrated", "joint", "count_log_score")
    uncalibrated_log_scores = _aggregate_by_patient(records, "uncalibrated", "joint", "count_log_score")
    contrasts["calibrated_minus_uncalibrated_joint_count_log_score"] = {
        **_bootstrap_difference(
            calibrated_log_scores,
            uncalibrated_log_scores,
            seed=seed + len(contrasts),
        ),
        "direction": "positive favors calibrated",
    }
    return {
        "status": "exploratory; MU-Glioma-Post has prior project use",
        "endpoint": "future annotated label-1/3 additions and removals inside forecast-origin risk sets",
        "forecast_horizon": dict(FORECAST_HORIZON),
        "method": "direction-specific spatial-prior residual heads plus bounded shared-frailty binomial counts, rank-constrained probabilities, and a fixed equal-weight hybrid/prior probability mixture",
        "n_patients": int(unique_patients.size),
        "n_intervals": len(intervals),
        "n_splits": n_splits,
        "fixed_hyperparameters": {
            "ridge": 1.0,
            "prior_penalty": 1e-3,
            "max_voxels_per_class_per_interval": 2000,
            "correlated_score_field_used_for_primary_scores": False,
            "count_rate_calibration": "direction-specific moment matching using inner patient-grouped out-of-fold forecasts",
            "prior_anchored_mixture_weight": 0.5,
        },
        "folds": fold_summaries,
        "patient_bootstrap_contrasts": contrasts,
        "interval_records": records,
        "limitations": [
            "Observed segmentation transitions are not verified biological growth or treatment response.",
            "Count predictions use the realized next-scan interval from retrospective data; they are conditional-horizon estimates, not pre-schedule forecasts.",
            "Boundary perturbation radius is a sensitivity assumption, not measured segmentation noise.",
            "MU results are exploratory because the cohort has prior project use.",
            "Spatial heads use fixed hyperparameters; no inner-fold tuning was performed.",
            "Count rate correction matches aggregate inner-fold expected counts and does not guarantee calibration for individual patients or intervals.",
            "The equal-weight hybrid/prior probability mixture is a fixed model-averaging sensitivity arm, not a calibrated or learned weight.",
            "The calibration layer was added after exploratory results on this same cohort; nested calibration limits within-run leakage but does not create an independent confirmation set.",
            "One shared frailty induces count dependence; residual local dependence is represented only in separate sensitivity draws.",
        ],
    }


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=BASE.DATA_ROOT)
    parser.add_argument("--patient-manifest", type=Path, default=BASE.PATIENT_MANIFEST)
    parser.add_argument("--pair-manifest", type=Path, default=BASE.PAIR_MANIFEST)
    parser.add_argument("--n-splits", type=int, default=5)
    parser.add_argument("--max-patients", type=int)
    parser.add_argument("--perturbation-radius-mm", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--acknowledge-exploratory-mu", action="store_true")
    args = parser.parse_args()
    if not args.acknowledge_exploratory_mu:
        parser.error("pass --acknowledge-exploratory-mu to run exploratory MU data analysis")
    output = args.output or DEFAULT_OUTPUT_DIR / f"track_b_signed_cv_{datetime.now(timezone.utc):%Y%m%dT%H%M%SZ}.json"
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")
    start = time.perf_counter()
    parameters = FEATURES.DiffusionParameters()
    intervals = _load_signed_intervals(
        data_root=args.data_root,
        patient_manifest=args.patient_manifest,
        pair_manifest=args.pair_manifest,
        parameters=parameters,
        max_patients=args.max_patients,
        perturbation_radius_mm=args.perturbation_radius_mm,
    )
    result = evaluate(intervals, n_splits=args.n_splits, seed=args.seed)
    result["runtime_seconds"] = float(time.perf_counter() - start)
    result["source_paths"] = {
        "data_root": _logical_input_path(args.data_root),
        "patient_manifest": _logical_input_path(args.patient_manifest),
        "pair_manifest": _logical_input_path(args.pair_manifest),
    }
    result["run_configuration"] = {
        "n_splits": args.n_splits,
        "max_patients": args.max_patients,
        "seed": args.seed,
        "perturbation_radius_mm": args.perturbation_radius_mm,
        "diffusion_parameters": {
            "edge_sensitivity": parameters.edge_sensitivity,
            "screening_length_mm": parameters.screening_length_mm,
            "label2_log_permeability": parameters.label2_log_permeability,
            "candidate_radius_mm": parameters.candidate_radius_mm,
            "prior_smoothing_mm": parameters.prior_smoothing_mm,
            "block_cavity": parameters.block_cavity,
            "solver_tolerance": parameters.solver_tolerance,
            "solver_max_iterations": parameters.solver_max_iterations,
        },
    }
    result["run_provenance"] = BASE.RUN_PROVENANCE.build_run_provenance(
        BASE.ROOT,
        manifests={
            "patient_manifest": args.patient_manifest,
            "pair_manifest": args.pair_manifest,
        },
        source_files={
            "signed_transition_runner": Path(__file__),
            "growth_zone_runner": Path(BASE.__file__),
            "growth_zone_features": Path(FEATURES.__file__),
            "signed_transition_model": Path(TRANSITIONS.__file__),
            "signed_spatial_heads": Path(HEADS.__file__),
            "run_provenance_helper": Path(BASE.RUN_PROVENANCE.__file__),
        },
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    result = _json_safe(result)
    output.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
    print(f"wrote exploratory signed-transition result: {output}", flush=True)
    print(json.dumps({key: result[key] for key in ("status", "n_patients", "n_intervals", "n_splits", "runtime_seconds", "patient_bootstrap_contrasts")}, indent=2, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
