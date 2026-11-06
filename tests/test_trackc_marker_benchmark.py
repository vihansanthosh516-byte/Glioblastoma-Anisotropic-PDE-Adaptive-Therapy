"""Pure checks for the hypothetical marker policy and benchmark definitions.

These tests use hand-built arrays. They do not run Track C simulations.
"""
from __future__ import annotations

import sys
import importlib.util
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
import trackc_regime as rg  # noqa: E402
import trackc_sim as tc  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "trackc_marker_benchmark", REPO / "src" / "206_trackc_marker_benchmark.py"
)
marker_benchmark = importlib.util.module_from_spec(_spec)
sys.modules[_spec.name] = marker_benchmark
_spec.loader.exec_module(marker_benchmark)
MARKER_LEVELS = marker_benchmark.MARKER_LEVELS
apply_marker_schedule_policy = marker_benchmark.apply_marker_schedule_policy
choose_marker_schedules = marker_benchmark.choose_marker_schedules
comparison_policy_specs = marker_benchmark.comparison_policy_specs
fit_marker_schedule_policy = marker_benchmark.fit_marker_schedule_policy
frontier_at = marker_benchmark.frontier_at
marker_from_uniforms = marker_benchmark.marker_from_uniforms
marker_policy_outcomes = marker_benchmark.marker_policy_outcomes
robust_training_scores = marker_benchmark.robust_training_scores
select_global_fixed_schedule = marker_benchmark.select_global_fixed_schedule
upper_hull = marker_benchmark.upper_hull
protocol_record = marker_benchmark._protocol
SELECTION_RULE_DESCRIPTION = marker_benchmark.SELECTION_RULE_DESCRIPTION


def _selection_outcomes(ttp_rano, courses_rano, ttp_baseline=None, courses_baseline=None):
    ttp_baseline = ttp_rano if ttp_baseline is None else ttp_baseline
    courses_baseline = courses_rano if courses_baseline is None else courses_baseline
    return {
        "3b_rano": (np.asarray(ttp_rano, float), np.asarray(courses_rano, float)),
        "3a_baseline": (np.asarray(ttp_baseline, float), np.asarray(courses_baseline, float)),
    }


def test_marker_from_uniforms_hits_requested_sensitivity_and_specificity():
    responder = np.array([True] * 10 + [False] * 10)
    uniforms = np.tile(np.arange(10) / 10.0, 2)
    calls = marker_from_uniforms(responder, uniforms, sensitivity=0.70, specificity=0.80)

    assert calls[responder].mean() == pytest.approx(0.70)
    assert (~calls[~responder]).mean() == pytest.approx(0.80)
    assert len(MARKER_LEVELS) == 6 and MARKER_LEVELS[0] == 0.50 and MARKER_LEVELS[-1] == 1.00


def test_robust_score_uses_worse_endpoint_and_mean_courses_from_both_endpoints():
    ttp_rano = np.array([[500.0, 500.0], [400.0, 400.0]])
    ttp_baseline = np.array([[200.0, 200.0], [450.0, 450.0]])
    courses_rano = np.array([[2.0, 2.0], [0.0, 0.0]])
    courses_baseline = np.array([[4.0, 4.0], [0.0, 0.0]])

    score, cost = robust_training_scores(
        ttp_rano, courses_rano, ttp_baseline, courses_baseline, np.array([0, 1]), lam=5.0
    )

    np.testing.assert_allclose(cost, [3.0, 0.0])
    np.testing.assert_allclose(score, [180.0, 400.0])


def test_benchmark_selector_uses_patientwise_worst_endpoint_utility():
    # Aggregating endpoints first would favor schedule 0, which wins one endpoint
    # per patient. The plan's patient-wise worst-endpoint utility must select 1.
    outcomes = _selection_outcomes(
        [[730.0, 100.0], [300.0, 300.0]],
        [[0.0, 0.0], [0.0, 0.0]],
        ttp_baseline=[[100.0, 730.0], [300.0, 300.0]],
        courses_baseline=[[0.0, 0.0], [0.0, 0.0]],
    )

    result = select_global_fixed_schedule(outcomes, ["split", "balanced"], lam=0.0)

    assert result["schedule"] == "balanced"
    assert result["training_robust_utility"] == pytest.approx(300.0)


def test_protocol_metadata_matches_patientwise_selection_rule():
    protocol = protocol_record("mixture")

    assert protocol["selection_rule"] == SELECTION_RULE_DESCRIPTION
    assert "mean across selection patients of min" in protocol["selection_rule"]
    assert "lambda*courses_RANO" in protocol["selection_rule"]
    assert "lambda*courses_baseline" in protocol["selection_rule"]


def test_schedule_selector_does_not_turn_near_ties_into_course_ties():
    outcomes = _selection_outcomes(
        [[500.0000005], [500.0]],
        [[1.0], [0.0]],
    )

    result = select_global_fixed_schedule(outcomes, ["slightly_better", "lower_courses"], lam=0.0)

    assert result["schedule"] == "slightly_better"


def test_schedule_selector_breaks_exact_utility_ties_by_fewer_courses():
    outcomes = _selection_outcomes(
        [[500.0], [500.0]],
        [[2.0], [0.0]],
    )

    result = select_global_fixed_schedule(outcomes, ["more_courses", "fewer_courses"], lam=0.0)

    assert result["schedule"] == "fewer_courses"


def test_fit_selects_different_fixed_schedules_for_observed_marker_groups():
    marker = np.array([False, False, True, True])
    outcomes = _selection_outcomes(
        [[500, 500, 100, 100], [100, 100, 500, 500]],
        [[1, 1, 1, 1], [1, 1, 1, 1]],
    )

    fitted = fit_marker_schedule_policy(marker, outcomes, ["negative_best", "positive_best"], lam=0.0,
                                        min_group_n=2)

    assert fitted["negative"]["schedule"] == "negative_best"
    assert fitted["positive"]["schedule"] == "positive_best"
    assert not fitted["negative"]["used_global_fallback"]
    assert not fitted["positive"]["used_global_fallback"]


def test_small_observed_group_falls_back_to_global_selection():
    marker = np.array([False, False, False, True])
    outcomes = _selection_outcomes(
        [[400, 400, 400, 200], [100, 100, 100, 500]],
        [[1, 1, 1, 1], [1, 1, 1, 1]],
    )

    fitted = choose_marker_schedules(marker, outcomes, ["global_best", "minority_best"], lam=0.0,
                                     min_group_n=2)

    assert fitted["positive"]["used_global_fallback"]
    assert fitted["positive"]["schedule"] == fitted["negative"]["global_schedule"]
    assert fitted["negative"]["selection_n"] == 3


def test_course_penalty_changes_selected_schedule_and_tie_prefers_fewer_courses():
    marker = np.zeros(4, bool)
    outcomes = _selection_outcomes(
        [[500, 500, 500, 500], [480, 480, 480, 480]],
        [[2, 2, 2, 2], [0, 0, 0, 0]],
    )

    time_only = fit_marker_schedule_policy(marker, outcomes, ["more_courses", "fewer_courses"], lam=0.0)
    cost_aware = fit_marker_schedule_policy(marker, outcomes, ["more_courses", "fewer_courses"], lam=10.0)

    assert time_only["negative"]["schedule"] == "more_courses"
    assert cost_aware["negative"]["schedule"] == "fewer_courses"


def test_global_fixed_comparator_uses_separate_selection_outcomes_and_course_cost():
    outcomes = _selection_outcomes(
        [[500, 500, 500, 500], [480, 480, 480, 480]],
        [[2, 2, 2, 2], [0, 0, 0, 0]],
    )

    time_only = select_global_fixed_schedule(outcomes, ["more_courses", "fewer_courses"], lam=0.0)
    cost_aware = select_global_fixed_schedule(outcomes, ["more_courses", "fewer_courses"], lam=10.0)

    assert time_only["schedule"] == "more_courses"
    assert cost_aware["schedule"] == "fewer_courses"


def test_pure_application_uses_only_marker_calls_and_selected_schedule_indices():
    marker = np.array([True, False, True])
    ttp = np.array([[10.0, 11.0, 12.0], [20.0, 21.0, 22.0]])
    courses = np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]])
    fitted = {
        "positive": {"schedule_index": 1},
        "negative": {"schedule_index": 0},
    }

    ttp_result, courses_result, schedule_index = apply_marker_schedule_policy(marker, fitted, ttp, courses)
    ttp_direct, courses_direct, direct_index = marker_policy_outcomes(marker, 1, 0, ttp, courses)

    np.testing.assert_array_equal(schedule_index, [1, 0, 1])
    np.testing.assert_array_equal(ttp_result, [20.0, 11.0, 22.0])
    np.testing.assert_array_equal(courses_result, [2.0, 1.0, 2.0])
    np.testing.assert_array_equal(ttp_result, ttp_direct)
    np.testing.assert_array_equal(courses_result, courses_direct)
    np.testing.assert_array_equal(schedule_index, direct_index)


def test_fixed_frontier_interpolates_only_inside_supported_course_range():
    hull = upper_hull([(0, 100), (1, 200), (2, 250), (3, 240)])

    assert frontier_at(hull, 1.5) == pytest.approx(225.0)
    assert np.isnan(frontier_at(hull, -0.1))
    assert np.isnan(frontier_at(hull, 3.1))


def test_comparator_catalog_covers_fixed_family_scan_rules_and_current_pf_rules():
    schedules, names = rg.candidate_schedules()
    arms = comparison_policy_specs()

    assert schedules.shape == (56, tc.N_CYCLES)
    assert len(names) == 56
    assert {"soc6", "raw_zhang_1.0_0.5", "raw_0.6_0.6", "nadir_k3_on1.25", "resp_k3_th0.75"} <= set(arms)
    assert {"pf_wide_containment", "pf_tail_trigger_candidate", "pf_holiday_containment"} <= set(arms)
    assert any(name.startswith("every2_") for name in names)
