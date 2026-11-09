from __future__ import annotations

from itertools import product

import numpy as np
import pytest

from trackc_ge_betting import (
    _accepts_mean,
    _validated_delta,
    _validated_unit_values,
    _variance_path,
    ge_betting_interval,
    ge_betting_interval_unit,
)


def test_unit_interval_validates_domain_and_delta():
    for values in ([], [[0.2, 0.3]], [np.nan], [-0.01], [1.01]):
        with pytest.raises(ValueError):
            _validated_unit_values(values)
    for delta in (0.0, 1.0, np.nan, True, "0.05"):
        with pytest.raises(ValueError):
            _validated_delta(delta)


def test_variance_path_matches_predictable_formula():
    values = np.asarray([0.2, 0.8, 0.4])
    path = _variance_path(values)
    means = [0.5]
    variances = [0.25]
    for t, _ in enumerate(values[:-1], start=1):
        variances.append(
            (0.25 + sum((values[j - 1] - means[j - 1]) ** 2 for j in range(1, t + 1)))
            / (t + 1)
        )
        means.append((0.5 + float(values[:t].sum())) / (t + 1))
    assert path == pytest.approx(variances)


@pytest.mark.parametrize("values,truth", [([0.0] * 7, 0.0), ([1.0] * 7, 1.0), ([0.5] * 7, 0.5)])
def test_deterministic_interval_contains_extreme_and_constant_means(values, truth):
    interval = ge_betting_interval_unit(values, delta=0.025)
    assert interval is not None
    assert 0.0 <= interval[0] <= truth <= interval[1] <= 1.0


def test_accepted_candidate_means_form_one_interval_for_small_binary_sequences():
    grid = np.linspace(0.0, 1.0, 201)
    for values in product((0.0, 1.0), repeat=4):
        x = _validated_unit_values(values)
        variance = _variance_path(x)
        accepted = [
            _accepts_mean(x, float(mean), 0.025, variance)
            for mean in grid
        ]
        positions = np.flatnonzero(accepted)
        if positions.size:
            assert positions[-1] - positions[0] + 1 == positions.size


def test_exact_bernoulli_coverage_for_n7_grid():
    n = 7
    delta = 0.025
    sequences = list(product((0.0, 1.0), repeat=n))
    intervals = [ge_betting_interval_unit(sequence, delta=delta) for sequence in sequences]
    for probability in np.linspace(0.0, 1.0, 41):
        coverage = 0.0
        for sequence, interval in zip(sequences, intervals):
            successes = int(sum(sequence))
            mass = probability**successes * (1.0 - probability) ** (n - successes)
            if interval is not None and interval[0] - 1e-12 <= probability <= interval[1] + 1e-12:
                coverage += mass
        assert coverage >= 1.0 - delta - 1e-12


def test_project_all_zero_sample_obeys_documented_rare_event_lower_bound():
    # D=0 maps to X=0.5. With n=7 and per-contrast delta=.025, a deterministic
    # uniformly valid interval must include means near +/-323.6 days here.
    intervals = ge_betting_interval([0.0] * 7, support=(-790.0, 790.0), delta=0.025)
    assert intervals is not None
    assert intervals[0] <= -323.5
    assert intervals[1] >= 323.5


def test_exact_coverage_at_seven_worlds_for_rare_event_bound_distributions():
    n = 7
    delta = 0.025
    probability = 0.409
    for extreme, truth in (
        (1.0, 0.5 + probability / 2.0),
        (0.0, 0.5 - probability / 2.0),
    ):
        coverage = 0.0
        for events in product((False, True), repeat=n):
            values = [extreme if event else 0.5 for event in events]
            event_count = sum(events)
            mass = probability**event_count * (1.0 - probability) ** (n - event_count)
            interval = ge_betting_interval_unit(values, delta=delta)
            if interval is not None and interval[0] - 1e-12 <= truth <= interval[1] + 1e-12:
                coverage += mass
        assert coverage >= 1.0 - delta - 1e-12


def test_project_support_mapping_and_rejection():
    interval = ge_betting_interval([-10.0, 0.0, 10.0], support=(-790.0, 790.0))
    assert interval is not None
    assert -790.0 <= interval[0] <= interval[1] <= 790.0
    with pytest.raises(ValueError, match="outside the declared support"):
        ge_betting_interval([791.0], support=(-790.0, 790.0))
    with pytest.raises(ValueError, match="increasing"):
        ge_betting_interval([0.0], support=(1.0, 1.0))
