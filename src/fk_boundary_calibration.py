"""Estimate Track B effective boundary covariance from patient-held-out residuals.

This is a method-of-moments variogram fit for the planned exponential-plus-nugget covariance. Each case must
record which patients trained its forecast; the held-out patient must be absent from that list.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import nnls


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
