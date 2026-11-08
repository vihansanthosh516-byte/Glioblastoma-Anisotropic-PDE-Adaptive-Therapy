import numpy as np
import pytest

from trackc_max_t_coverage_audit import _draw_standardized, _wilson_interval, run_coverage_audit


def test_standardized_distribution_draws_are_seeded_and_finite():
    for name in ("normal", "gamma2", "student_t3", "uniform"):
        first = np.random.Generator(np.random.PCG64DXSM(17))
        second = np.random.Generator(np.random.PCG64DXSM(17))
        np.testing.assert_array_equal(_draw_standardized(name, first, 25), _draw_standardized(name, second, 25))
    with pytest.raises(ValueError, match="unknown distribution"):
        _draw_standardized("unknown", np.random.default_rng(1), 3)


def test_wilson_interval_bounds_zero_and_all_success_counts():
    low_zero, high_zero = _wilson_interval(0, 100, 0.95)
    low_all, high_all = _wilson_interval(100, 100, 0.95)
    assert 0.0 <= low_zero < high_zero < 0.1
    assert 0.9 < low_all < high_all <= 1.0


def test_fast_coverage_audit_has_all_predeclared_pair_cells_and_no_patient_data():
    first = run_coverage_audit(repetitions=8, bootstrap_resamples=20, world_counts=(5,))
    second = run_coverage_audit(repetitions=8, bootstrap_resamples=20, world_counts=(5,))
    assert first == second
    assert len(first["cells"]) == len(first["distribution_pairs"]) == 6
    assert first["patient_data_used"] is False
    assert first["trackc_policy_outcomes_used"] is False
    assert all(cell["monte_carlo_repetitions"] == 8 for cell in first["cells"])


@pytest.mark.parametrize(
    "kwargs",
    [
        {"repetitions": 0},
        {"bootstrap_resamples": 0},
        {"world_counts": ()},
        {"world_counts": (1,)},
        {"world_counts": (True,)},
    ],
)
def test_coverage_audit_rejects_invalid_settings(kwargs):
    with pytest.raises(ValueError):
        run_coverage_audit(**kwargs)
