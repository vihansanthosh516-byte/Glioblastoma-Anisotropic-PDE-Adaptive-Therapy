from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_boundary_calibration import (  # noqa: E402
    assign_patient_folds,
    fit_cross_fitted_boundary_covariance,
    make_cross_fitted_boundary_residual_case,
    run_cross_fitted_boundary_calibration,
)
from src.fk_patient_cv import PatientSeries  # noqa: E402
from src.fk_spatial_uncertainty import exponential_covariance_matrix  # noqa: E402


def test_cross_fitted_variogram_estimator_recovers_simulated_covariance_scale():
    points = np.asarray([(x * 4.0, y * 4.0, 0.0) for x in range(4) for y in range(4)])
    covariance = exponential_covariance_matrix(
        points, sigma_b_mm=2.0, tau_b_mm=0.5, ell_mm=10.0,
    )
    factor = np.linalg.cholesky(covariance)
    rng = np.random.default_rng(782)
    cases = []
    for index in range(200):
        patient_id = f"p{index}"
        residuals = factor @ rng.normal(size=len(points))
        cases.append({
            "patient_id": patient_id,
            "origin_id": "origin0",
            "training_patient_ids": [f"p{other}" for other in range(200) if other != index],
            "residuals_mm": residuals,
            "anchor_points_mm": points,
        })

    result = fit_cross_fitted_boundary_covariance(cases)

    assert result["ell_mm"] == 10.0
    assert abs(result["sigma_b_mm"] - 2.0) < 0.5
    assert abs(result["tau_b_mm"] - 0.5) < 0.35
    assert result["n_calibration_patients"] == 200
    assert result["n_anchor_pairs"] == 200 * 120
    assert "not scanner noise" in result["scope"]


def test_covariance_fit_rejects_calibration_patient_in_its_training_set():
    points = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    case = {
        "patient_id": "p1",
        "training_patient_ids": ["p1", "p2"],
        "residuals_mm": np.asarray([0.0, 1.0, -1.0]),
        "anchor_points_mm": points,
    }

    with pytest.raises(ValueError, match="excluded from its forecast training set"):
        fit_cross_fitted_boundary_covariance([case, {**case, "patient_id": "p2"}])


def test_residual_case_uses_held_out_target_only_after_forecast():
    shape = (16, 16, 16)
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    target = sum((axis - value) ** 2 for axis, value in zip(coords, center)) <= 3.0**2
    shifted = np.roll(target, 1, axis=0)

    case = make_cross_fitted_boundary_residual_case(
        patient_id="heldout", origin_id="visit2", training_patient_ids=["train1", "train2"],
        target_mask=target, pre_boundary_probability=shifted.astype(float),
        spacing_mm=1.0, n_anchors=24,
    )

    assert case["residuals_mm"].shape == (24,)
    assert np.any(np.abs(case["residuals_mm"]) > 0)
    assert case["patient_id"] not in case["training_patient_ids"]
    assert case["target_used_for"] == "held-out residual measurement only"


def test_residual_case_rejects_target_patient_in_training_ids():
    target = np.zeros((8, 8, 8), dtype=bool)
    target[2:5, 2:5, 2:5] = True

    with pytest.raises(ValueError, match="held-out patient must not appear"):
        make_cross_fitted_boundary_residual_case(
            patient_id="p1", origin_id="o1", training_patient_ids=["p1", "p2"],
            target_mask=target, pre_boundary_probability=target.astype(float), n_anchors=10,
        )


def _series(patient_id):
    shape = (16, 16, 16)
    base = np.zeros(shape, dtype=bool)
    base[5:11, 5:11, 5:11] = True
    masks = [base.copy(), base.copy(), base.copy(), np.roll(base, 1, axis=0)]
    return PatientSeries(patient_id, np.asarray([0.0, 20.0, 40.0, 60.0]), tuple(masks))


def test_inner_patient_folds_are_balanced_deterministic_and_cover_each_id_once():
    patient_ids = [f"p{i}" for i in range(12)]
    first = assign_patient_folds(patient_ids, n_folds=4, seed=101)
    second = assign_patient_folds(patient_ids, n_folds=4, seed=101)

    assert first == second
    assert set(first) == set(patient_ids)
    assert {fold: list(first.values()).count(fold) for fold in range(4)} == {0: 3, 1: 3, 2: 3, 3: 3}


def test_cross_fitted_boundary_calibration_excludes_each_patient_from_its_forecast_fit():
    records = [_series(f"p{i}") for i in range(6)]
    forecasted = []

    def fit_nuisance(training_records):
        return {"training_ids": {record.patient_id for record in training_records}}

    def forecast(history, nuisance):
        assert not hasattr(history, "target")
        assert history.patient_id not in nuisance["training_ids"]
        forecasted.append(history.patient_id)
        return {"pre_boundary_probability": history.masks[-1].astype(float)}

    result = run_cross_fitted_boundary_calibration(
        records, fit_nuisance, forecast, n_folds=3, seed=33, n_anchors=24,
        ell_candidates_mm=(2.0, 5.0, 10.0),
    )

    assert len(forecasted) == len(set(forecasted)) == 6
    assert result["covariance_fit"]["n_calibration_patients"] == 6
    assert len(result["residual_cases"]) == 6
    assert all(case["patient_id"] not in case["training_patient_ids"]
               for case in result["residual_cases"])
    assert all(fold["n_test_patients"] == 2 for fold in result["fold_audit"])


def test_cross_fitted_boundary_calibration_requires_pre_boundary_probability():
    records = [_series(f"p{i}") for i in range(4)]
    with pytest.raises(ValueError, match="pre_boundary_probability"):
        run_cross_fitted_boundary_calibration(
            records, lambda training: {}, lambda history, nuisance: {},
            n_folds=2, seed=5, n_anchors=12, ell_candidates_mm=(2.0, 5.0),
        )
