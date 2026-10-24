from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_boundary_calibration import fit_cross_fitted_boundary_covariance  # noqa: E402
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
