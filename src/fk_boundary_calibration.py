"""Estimate Track B effective boundary covariance from patient-held-out residuals.

This is a method-of-moments variogram fit for the planned exponential-plus-nugget covariance. Each case must
record which patients trained its forecast; the held-out patient must be absent from that list.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
from scipy import ndimage
from scipy.optimize import nnls

from src.fk_spatial_uncertainty import deterministic_surface_anchors


def make_cross_fitted_boundary_residual_case(*, patient_id, origin_id, training_patient_ids,
                                             target_mask, pre_boundary_probability,
                                             spacing_mm=1.0, n_anchors: int = 256) -> dict:
    """Create one calibration residual vector after forecasting a held-out patient from history only.

    The target mask is used only to define target-surface anchors and residual values. The caller must provide
    the patient IDs used to fit every forecast nuisance parameter; this function rejects self-inclusion.
    """
    if patient_id is None or origin_id is None:
        raise ValueError("patient_id and origin_id must be present")
    patient_id = str(patient_id).strip()
    origin_id = str(origin_id).strip()
    if not patient_id or not origin_id:
        raise ValueError("patient_id and origin_id must be nonempty")
    if isinstance(training_patient_ids, (str, bytes)):
        raise ValueError("training_patient_ids must be a collection of patient identifiers")
    training_ids = [str(value).strip() for value in training_patient_ids]
    if not training_ids or any(not value for value in training_ids):
        raise ValueError("training_patient_ids must contain nonempty identifiers")
    if patient_id in set(training_ids):
        raise ValueError("held-out patient must not appear in forecast training IDs")

    target = np.asarray(target_mask, dtype=bool)
    probability = np.asarray(pre_boundary_probability, dtype=float)
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if (target.ndim != 3 or probability.shape != target.shape or not target.any() or target.all()
            or not np.isfinite(probability).all() or np.any((probability < 0) | (probability > 1))):
        raise ValueError("target and forecast probability must be valid matching 3D arrays")
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    forecast_mask = probability >= 0.5
    if not forecast_mask.any() or forecast_mask.all():
        raise ValueError("median pre-boundary forecast must be nonempty and not fill the domain")

    anchors = deterministic_surface_anchors(target, spacing_mm=spacing, n_anchors=n_anchors)
    anchor_indices = (anchors / spacing).T
    target_distance = (ndimage.distance_transform_edt(target, sampling=spacing)
                       - ndimage.distance_transform_edt(~target, sampling=spacing))
    forecast_distance = (ndimage.distance_transform_edt(forecast_mask, sampling=spacing)
                         - ndimage.distance_transform_edt(~forecast_mask, sampling=spacing))
    observed = ndimage.map_coordinates(target_distance, anchor_indices, order=1, mode="nearest")
    predicted = ndimage.map_coordinates(forecast_distance, anchor_indices, order=1, mode="nearest")
    return {
        "patient_id": patient_id,
        "origin_id": origin_id,
        "training_patient_ids": training_ids,
        "residuals_mm": np.asarray(predicted - observed, dtype=float),
        "anchor_points_mm": anchors,
        "target_used_for": "held-out residual measurement only",
    }


def fit_cross_fitted_boundary_covariance(calibration_cases, *, ell_candidates_mm=(5.0, 10.0, 20.0)) -> dict:
    """Fit sigma_b, tau_b, and ell from patient-held-out boundary residual vectors.

    Each record contains patient_id, training_patient_ids, residuals_mm, and anchor_points_mm. For distinct
    anchors, half the squared residual difference has expectation
    sigma_b^2 * (1 - exp(-distance / ell)) + tau_b^2. Origins receive equal total weight within patient.
    """
    cases = list(calibration_cases)
    ell_values = tuple(float(value) for value in ell_candidates_mm)
    if len(cases) < 2:
        raise ValueError("need at least two patient-held-out calibration cases")
    if len(ell_values) < 2 or any(not np.isfinite(value) or value <= 0 for value in ell_values):
        raise ValueError("need at least two positive finite ell candidates")
    if len(set(ell_values)) != len(ell_values):
        raise ValueError("ell candidates must be unique")

    normalized = []
    origin_counts = {}
    seen_origins = set()
    for case in cases:
        patient_value = case.get("patient_id")
        origin_value = case.get("origin_id", len(normalized))
        if patient_value is None or origin_value is None:
            raise ValueError("patient_id and origin_id must be present")
        patient_id = str(patient_value).strip()
        origin_id = str(origin_value).strip()
        raw_training_ids = case["training_patient_ids"]
        if isinstance(raw_training_ids, (str, bytes)):
            raise ValueError("training_patient_ids must be a collection of patient identifiers")
        training_values = [str(value).strip() for value in raw_training_ids]
        if not training_values or any(not value for value in training_values):
            raise ValueError("training_patient_ids must contain nonempty identifiers")
        training_ids = set(training_values)
        if not patient_id or patient_id in training_ids:
            raise ValueError("each calibration patient's identity must be excluded from its forecast training set")
        if not origin_id or (patient_id, origin_id) in seen_origins:
            raise ValueError("calibration origin identifiers must be nonempty and unique within patient")
        seen_origins.add((patient_id, origin_id))
        residuals = np.asarray(case["residuals_mm"], dtype=float)
        points = np.asarray(case["anchor_points_mm"], dtype=float)
        if (residuals.ndim != 1 or len(residuals) < 3 or not np.isfinite(residuals).all()
                or points.shape != (len(residuals), 3) or not np.isfinite(points).all()):
            raise ValueError("each case needs at least three finite residuals and matching finite 3D anchors")
        pair_i, pair_j = np.triu_indices(len(residuals), k=1)
        distances = np.linalg.norm(points[pair_i] - points[pair_j], axis=1)
        keep = distances > 0
        if not keep.any():
            raise ValueError("each case needs at least one pair of distinct anchor locations")
        semivariogram = 0.5 * (residuals[pair_i[keep]] - residuals[pair_j[keep]]) ** 2
        normalized.append((patient_id, distances[keep], semivariogram))
        origin_counts[patient_id] = origin_counts.get(patient_id, 0) + 1

    if len(origin_counts) < 2:
        raise ValueError("need calibration origins from at least two distinct patients")

    weights = []
    for patient_id, distances, semivariogram in normalized:
        pair_weight = 1.0 / np.sqrt(len(distances) * origin_counts[patient_id])
        weights.append(np.full(len(distances), pair_weight, dtype=float))
    total_weight_squared = float(sum(weight @ weight for weight in weights))

    candidate_fits = {}
    for ell in ell_values:
        design_parts, response_parts = [], []
        for (_, distances, semivariogram), case_weights in zip(normalized, weights):
            design = np.column_stack((1.0 - np.exp(-distances / ell), np.ones(len(distances))))
            design_parts.append(design * case_weights[:, None])
            response_parts.append(semivariogram * case_weights)
        coefficients, residual_norm = nnls(np.vstack(design_parts), np.concatenate(response_parts))
        candidate_fits[ell] = {
            "sigma_b_mm": float(np.sqrt(coefficients[0])),
            "tau_b_mm": float(np.sqrt(coefficients[1])),
            "weighted_rmse_mm2": float(residual_norm / np.sqrt(total_weight_squared)),
        }
    selected_ell = min(ell_values, key=lambda value: candidate_fits[value]["weighted_rmse_mm2"])
    selected = candidate_fits[selected_ell]
    return {
        "sigma_b_mm": selected["sigma_b_mm"],
        "tau_b_mm": selected["tau_b_mm"],
        "ell_mm": float(selected_ell),
        "candidate_fits": candidate_fits,
        "n_calibration_origins": int(len(cases)),
        "n_calibration_patients": int(len(origin_counts)),
        "n_anchor_pairs": int(sum(len(distances) for _, distances, _ in normalized)),
        "selection_boundary": bool(selected_ell == min(ell_values) or selected_ell == max(ell_values)),
        "estimator": "patient-balanced nonnegative method-of-moments semivariogram least squares",
        "does_not_model_mean_boundary_bias": True,
        "scope": "effective held-out mask, registration, treatment, and model discrepancy; not scanner noise",
    }


def assign_patient_folds(patient_ids, *, n_folds: int = 5, seed: int = 0) -> dict[str, int]:
    """Make deterministic, balanced folds after normalizing patient identifiers."""
    if isinstance(patient_ids, (str, bytes)):
        raise ValueError("patient_ids must be a collection of identifiers")
    ids = [str(value).strip() for value in patient_ids]
    if not ids or any(not value for value in ids):
        raise ValueError("patient_ids must be nonempty")
    if len(set(ids)) != len(ids):
        raise ValueError("patient_ids must be unique after normalization")
    if isinstance(n_folds, (bool, np.bool_)) or int(n_folds) != n_folds or n_folds < 2:
        raise ValueError("n_folds must be an integer of at least two")
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise ValueError("seed must be an integer")
    n_folds = int(n_folds)
    if n_folds > len(ids):
        raise ValueError("n_folds cannot exceed the number of patients")
    shuffled = np.asarray(sorted(ids), dtype=object)
    shuffled = shuffled[np.random.default_rng(int(seed)).permutation(len(shuffled))]
    return {str(patient_id): int(index % n_folds) for index, patient_id in enumerate(shuffled)}


def run_cross_fitted_boundary_calibration(records: Sequence, fit_nuisance, forecast, *,
                                          n_folds: int = 5, seed: int = 0,
                                          minimum_history: int = 2, spacing_mm=1.0,
                                          n_anchors: int = 256,
                                          ell_candidates_mm=(5.0, 10.0, 20.0)) -> dict:
    """Generate inner patient-held-out boundary residuals and fit their effective covariance.

    Call this only on an outer-training cohort. ``fit_nuisance`` sees inner-training patients only, and
    ``forecast`` sees one held-out patient's history plus the nuisance fit, never that patient's target mask.
    The target enters afterward only to measure a cross-fitted residual. The returned covariance includes
    model, treatment, registration, and mask-reference error; it is not a scanner-noise estimate.
    """
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)) or not records:
        raise ValueError("records must be a nonempty sequence of patient series")
    folds = assign_patient_folds([record.patient_id for record in records], n_folds=n_folds, seed=seed)
    residual_cases = []

    def fit_with_training_ids(training_records):
        training_ids = tuple(sorted(record.patient_id for record in training_records))
        return {"nuisance": fit_nuisance(training_records), "training_patient_ids": training_ids}

    def forecast_with_training_ids(history, context):
        raw_prediction = forecast(history, context["nuisance"])
        if not isinstance(raw_prediction, Mapping):
            raise ValueError("cross-fitted forecast must return a mapping with pre_boundary_probability")
        prediction = dict(raw_prediction)
        prediction["_training_patient_ids"] = context["training_patient_ids"]
        return prediction

    def score_residual(prediction, target):
        if "pre_boundary_probability" not in prediction:
            raise ValueError("cross-fitted forecast must return pre_boundary_probability")
        case = make_cross_fitted_boundary_residual_case(
            patient_id=target.patient_id,
            origin_id=f"target-{target.target_index}-day-{target.day:.12g}",
            training_patient_ids=prediction["_training_patient_ids"],
            target_mask=target.mask,
            pre_boundary_probability=prediction["pre_boundary_probability"],
            spacing_mm=spacing_mm,
            n_anchors=n_anchors,
        )
        residual_cases.append(case)
        return {"residual_case_index": len(residual_cases) - 1}

    from src.fk_patient_cv import run_outer_patient_cv
    evaluation = run_outer_patient_cv(
        records, folds, fit_with_training_ids, forecast_with_training_ids, score_residual,
        n_folds=n_folds, minimum_history=minimum_history,
    )
    covariance_fit = fit_cross_fitted_boundary_covariance(
        residual_cases, ell_candidates_mm=ell_candidates_mm,
    )
    return {
        "covariance_fit": covariance_fit,
        "residual_cases": residual_cases,
        "fold_audit": evaluation["fold_audit"],
        "n_residual_cases": len(residual_cases),
        "scope": "inner patient-held-out calibration on outer-training records only",
        "noise_interpretation": "effective forecast-residual covariance; not scanner noise",
    }
