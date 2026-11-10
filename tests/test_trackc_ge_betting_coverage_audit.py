"""Mechanics checks for the outcome-free GE-betting coverage screen."""
import math

import numpy as np
import pytest

import trackc_ge_betting as ge
import trackc_ge_betting_coverage_audit as audit


def test_task_grid_matches_protocol():
    tasks = audit.build_tasks()
    assert len(tasks) == 30
    assert [task[0] for task in tasks] == list(range(30))
    assert {task[1] for task in tasks} == {7, 14, 30, 100, 207}
    assert audit.DELTA_PER_FAMILY == 0.025
    assert audit.JOINT_COVERAGE_LOWER_BOUND == 0.93
    assert audit.PRECISION_TARGET_DAYS == 14.0
    assert audit.SUPPORT == (-790.0, 790.0)


def test_rare_event_probability_matches_all_zero_lower_bound():
    assert 790.0 * audit.rare_event_probability(7) == pytest.approx(323.60, abs=0.01)
    assert 790.0 * audit.rare_event_probability(200) == pytest.approx(14.44, abs=0.01)
    assert 790.0 * audit.rare_event_probability(207) < 14.0 < 790.0 * audit.rare_event_probability(206)


@pytest.mark.parametrize("name", ["rare_pos", "rare_neg", "two_point", "uniform", "narrow_normal", "skew_beta"])
def test_draws_stay_in_support_and_match_declared_mean(name):
    rng = np.random.Generator(np.random.PCG64DXSM(7))
    days, mean = audit.draw_days(name, rng, (4000, 30))
    assert days.min() >= -790.0 and days.max() <= 790.0
    standard_error = days.std() / math.sqrt(days.size)
    assert abs(days.mean() - mean) < 5.0 * standard_error + 1e-9


def test_closed_form_half_widths():
    assert audit.hoeffding_half_width_days(14) == pytest.approx(1580.0 * math.sqrt(math.log(80.0) / 28.0))
    constant = np.zeros((1, 10))
    assert audit.student_t_half_width_days(constant)[0] == 0.0
    expected = 1580.0 * 7.0 * math.log(160.0) / 27.0
    assert audit.empirical_bernstein_half_width_days(constant)[0] == pytest.approx(expected)


def test_wilson_interval_bounds():
    lower, upper = audit.wilson_interval(1950, 2000)
    assert 0.0 < lower < 0.975 < upper < 1.0
    with pytest.raises(ValueError):
        audit.wilson_interval(3, 2)


def test_bisection_interval_matches_dense_grid_acceptance():
    """The reported endpoints must bound the accepted set found by a grid scan."""
    rng = np.random.Generator(np.random.PCG64DXSM(11))
    grid = np.linspace(0.0, 1.0, 401)
    for n in (5, 7, 14, 30):
        for _ in range(6):
            x = rng.beta(0.7, 1.6, size=n)
            interval = ge.ge_betting_interval_unit(x, delta=0.025)
            assert interval is not None
            path = ge._variance_path(x)
            accepted = np.array([ge._accepts_mean(x, float(m), 0.025, path) for m in grid])
            inside = (grid >= interval[0] - 1e-9) & (grid <= interval[1] + 1e-9)
            assert not (accepted & ~inside).any()
            interior = (grid >= interval[0] + 1e-9) & (grid <= interval[1] - 1e-9)
            assert accepted[interior].all()


def test_small_cell_runs_and_reports_all_methods():
    cell = audit.run_cell((0, 7, "narrow_normal_narrow_normal", "narrow_normal", "narrow_normal", 40))
    assert set(cell["methods"]) == set(audit.METHODS)
    ge_row = cell["methods"]["ge_betting"]
    assert ge_row["inference_failures"] == 0
    assert ge_row["inversion_inconsistencies"] == 0
    assert 0.0 <= ge_row["joint_coverage"] <= 1.0
    # At n=7 the Hoeffding half-width exceeds the support, so the clipped
    # interval can be no wider than the support itself.
    assert audit.hoeffding_half_width_days(7) > 790.0
    assert 0.0 < cell["methods"]["hoeffding"]["median_larger_half_width_days"] <= 790.0
