import itertools
import math

import numpy as np
import pytest
from scipy.special import logsumexp

from src.growth_zone_joint_law import (
    _summarize_cell,
    fit_nb_regression,
    log_elementary_symmetric,
    negative_binomial_logpmf,
    sample_conditional_bernoulli,
    solve_mass_intercept,
)


def test_log_elementary_symmetric_matches_brute_force():
    scores = np.array([-0.7, 0.2, 0.9, -0.1, 0.4])
    k = 2
    exact = logsumexp(
        [sum(scores[index] for index in subset) for subset in itertools.combinations(range(len(scores)), k)]
    )
    assert log_elementary_symmetric(scores, k) == pytest.approx(exact, abs=1e-12)


def test_conditional_bernoulli_sampler_has_requested_cardinality():
    rng = np.random.default_rng(42)
    scores = np.array([-1.0, 0.2, 0.5, 1.1, -0.4, 0.8])
    for k in range(len(scores) + 1):
        mask = sample_conditional_bernoulli(scores, k, rng)
        assert mask.dtype == bool
        assert int(mask.sum()) == k


def test_mass_intercept_hits_expected_count():
    scores = np.array([-2.0, -0.4, 0.1, 0.5, 1.3, 2.1])
    probabilities = solve_mass_intercept(scores, expected_count=2.4)
    assert probabilities.sum() == pytest.approx(2.4, abs=1e-10)
    assert np.array_equal(np.argsort(scores), np.argsort(probabilities))


def test_negative_binomial_logpmf_matches_scipy():
    from scipy.stats import nbinom

    k = np.arange(8)
    mean = 4.2
    alpha = 0.7
    size = 1.0 / alpha
    probability = size / (size + mean)
    expected = nbinom.logpmf(k, size, probability)
    assert np.allclose(negative_binomial_logpmf(k, mean, alpha), expected, atol=1e-12)


def test_nb_regression_recovers_constant_mean_and_returns_finite_fit():
    rng = np.random.default_rng(7)
    n = 500
    previous = rng.integers(0, 20, size=n)
    ratio = rng.uniform(0.75, 1.25, size=n)
    features = np.column_stack([np.ones(n), np.log1p(previous), np.log(ratio)])
    observed = np.full(n, 8, dtype=int)
    fit = fit_nb_regression(features[:, :1], observed)
    assert np.isfinite(fit["dispersion"])
    assert fit["dispersion"] > 0
    assert math.exp(fit["coefficients"][0]) == pytest.approx(8.0, rel=0.03)


def test_summary_ignores_undefined_metrics_and_reports_empty_as_nan():
    summary = _summarize_cell([
        {"ap": 0.7, "other": 2.0},
        {"ap": float("nan"), "other": 4.0},
    ])
    assert summary["ap"] == pytest.approx(0.7)
    assert summary["other"] == pytest.approx(3.0)
    assert math.isnan(_summarize_cell([{"ap": float("nan")}])["ap"])
