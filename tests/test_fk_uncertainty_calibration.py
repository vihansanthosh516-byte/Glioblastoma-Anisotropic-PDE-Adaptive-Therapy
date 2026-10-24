from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_mask_inverse as fi  # noqa: E402
from src import fk_synthetic_validation as fv  # noqa: E402


def ball(shape=(16, 16, 16), radius=3.0):
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    return sum((axis - c) ** 2 for axis, c in zip(coords, center)) <= radius**2


def test_boundary_disagreement_is_zero_for_identical_and_positive_for_different_forecasts():
    mask = ball()
    moved = np.roll(mask, 2, axis=0)

    assert fv.ensemble_boundary_disagreement([mask, mask.copy()]) == 0.0
    score = fv.ensemble_boundary_disagreement([mask, moved])
    assert 0.0 < score <= 1.0


def test_boundary_disagreement_rejects_invalid_ensembles():
    mask = ball()
    with np.testing.assert_raises(ValueError):
        fv.ensemble_boundary_disagreement([mask])
    with np.testing.assert_raises(ValueError):
        fv.ensemble_boundary_disagreement([mask, np.zeros((4, 4, 4), bool)])


def test_conformal_error_bound_uses_calibration_cases_only_and_is_bounded():
    scores = np.array([0.1, 0.2, 0.3, 0.4])
    errors = np.array([0.1, 0.4, 0.2, 0.3])
    result = fv.conformal_error_upper_bounds(
        scores, errors, np.array([0.2, 0.5]), alpha=0.25, score_floor=1e-6
    )

    ratios = errors / (scores + 1e-6)
    assert result["quantile"] == np.max(ratios)
    assert result["calibration_n"] == 4
    assert result["bounds"].shape == (2,)
    assert np.all((result["bounds"] >= 0) & (result["bounds"] <= 1))
    assert result["bounds"][1] >= result["bounds"][0]


def test_conformal_error_bound_is_uninformative_when_calibration_set_is_too_small():
    result = fv.conformal_error_upper_bounds(
        [0.1, 0.2], [0.1, 0.2], [0.0, 0.5], alpha=0.1
    )
    assert np.all(result["bounds"] == 1.0)
    assert np.isinf(result["quantile"])


def test_conformal_error_bound_rejects_bad_inputs():
    with np.testing.assert_raises(ValueError):
        fv.conformal_error_upper_bounds([0.1], [1.2], [0.1])
    with np.testing.assert_raises(ValueError):
        fv.conformal_error_upper_bounds([0.1], [0.2], [0.1], alpha=1.0)
    with np.testing.assert_raises(ValueError):
        fv.conformal_error_upper_bounds([0.1], [0.2, 0.3], [0.1])


def test_bootstrap_forecast_ensemble_refits_only_perturbed_training_masks(monkeypatch):
    shape = (16, 16, 16)
    training = [ball(shape, radius=r) for r in (2.0, 2.5, 3.0)]
    days = [0.0, 10.0, 20.0]
    domain = ball(shape, radius=3.1)
    fit_inputs = []

    def fake_fit(masks, fit_days, h, domain=None):
        assert all(not np.any(mask & ~domain) for mask in masks)
        fit_inputs.append(([int(mask.sum()) for mask in masks], tuple(fit_days)))
        return {"rho": 0.04, "D": 0.05, "sse_log_volume": 0.25, "dice": [0.95, 0.8]}

    def fake_forward(mask0, domain, rho, D, forecast_days, h):
        assert tuple(forecast_days) == (0.0, 30.0)
        return [np.roll(mask0, 1, axis=0)]

    monkeypatch.setattr(fi, "fit_fk_masks", fake_fit)
    monkeypatch.setattr(fi, "forward_masks", fake_forward)
    first = fv.bootstrap_forecast_ensemble(
        training, days, 30.0, h=1.0, sigma_mm=0.5, correlation_mm=1.5,
        n_draws=3, seed=11, domain=domain,
    )
    second = fv.bootstrap_forecast_ensemble(
        training, days, 30.0, h=1.0, sigma_mm=0.5, correlation_mm=1.5,
        n_draws=3, seed=11, domain=domain,
    )

    assert len(fit_inputs) == 6
    assert all(days_seen == tuple(days) for _, days_seen in fit_inputs)
    assert len(first["forecasts"]) == 3
    assert all(np.array_equal(a, b) for a, b in zip(first["forecasts"], second["forecasts"]))
    assert first["disagreement"] == second["disagreement"]
    assert 0.0 <= first["disagreement"] <= 1.0


def test_next_scan_comparison_reports_training_spatial_residual(monkeypatch):
    shape = (16, 16, 16)
    training = [ball(shape, radius=r) for r in (2.0, 2.5, 3.0)]
    days = [0.0, 10.0, 20.0]

    monkeypatch.setattr(fi, "fit_fk_masks", lambda *args, **kwargs: {
        "rho": 0.04, "D": 0.05, "sse_log_volume": 0.25, "dice": [0.95, 0.8]
    })

    def fake_forward(mask0, domain, rho, D, forecast_days, h):
        return [training[1].copy(), np.roll(training[2], 1, axis=0)]

    monkeypatch.setattr(fi, "forward_masks", fake_forward)
    result = fv.compare_next_scan(training, days, 30.0, ball(shape, 4.0), h=1.0)

    assert result["training_spatial_residual"] > 0.0
    assert result["training_spatial_residual"] <= 1.0
    assert result["training_volume_residual_score"] == 0.25 / 1.25
