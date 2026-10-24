from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_volume_model import (  # noqa: E402
    fit_log_volume_variance,
    sample_positive_volumes,
    select_shrinkage_time_constant,
    shrink_log_growth_slope,
)


def test_shrinkage_respects_observation_gap_and_can_keep_negative_growth():
    short_gap = shrink_log_growth_slope(-0.01, 0.02, observation_gap_days=5, time_constant_days=90)
    long_gap = shrink_log_growth_slope(-0.01, 0.02, observation_gap_days=500, time_constant_days=90)

    assert short_gap > -0.01 and short_gap < 0.02
    assert long_gap < short_gap
    assert long_gap < 0.0


def test_time_constant_selection_uses_only_supplied_training_validation_cases():
    cases = [
        {"patient_slope": 0.02, "population_slope": 0.0, "observation_gap_days": 20,
         "horizon_days": 60, "last_volume": 100.0, "target_volume": 100.0},
        {"patient_slope": -0.01, "population_slope": 0.0, "observation_gap_days": 30,
         "horizon_days": 60, "last_volume": 100.0, "target_volume": 100.0},
    ]
    selected = select_shrinkage_time_constant(cases, candidates_days=(30, 90, 180, 365))

    assert selected["time_constant_days"] in (30, 90, 180, 365)
    assert set(selected["candidate_mse_log_volume"]) == {30, 90, 180, 365}
    altered = [dict(case) for case in cases]
    altered[0]["target_volume"] = 1000.0
    changed = select_shrinkage_time_constant(altered, candidates_days=(30, 90, 180, 365))
    assert changed["candidate_mse_log_volume"] != selected["candidate_mse_log_volume"]


def test_log_volume_variance_fit_is_nonnegative_and_predictive_draws_are_positive_reproducible():
    fit = fit_log_volume_variance([30, 60, 90, 120], [0.10, -0.15, 0.20, -0.25])
    first = sample_positive_volumes(
        100.0, horizon_days=60, log_growth_slope=-0.002, variance_fit=fit, n_draws=100, seed=9
    )
    second = sample_positive_volumes(
        100.0, horizon_days=60, log_growth_slope=-0.002, variance_fit=fit, n_draws=100, seed=9
    )

    assert fit["intercept"] >= 0 and fit["slope_per_day"] >= 0
    assert np.all(first > 0)
    assert np.array_equal(first, second)
    assert np.isclose(np.log(first).mean(), np.log(100.0) - 0.002 * 60, atol=0.08)


def test_zero_residual_variance_produces_deterministic_positive_volume():
    fit = fit_log_volume_variance([30, 60], [0.0, 0.0])
    draws = sample_positive_volumes(
        80.0, horizon_days=45, log_growth_slope=-0.01, variance_fit=fit, n_draws=4, seed=3
    )
    assert fit["intercept"] == 0.0 and fit["slope_per_day"] == 0.0
    assert np.allclose(draws, 80.0 * np.exp(-0.45))
