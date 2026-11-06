"""Tests for selection-only fixed-schedule dose matching."""
from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

from trackc_equal_course_frontier import randomize_schedule_mix, select_schedule_mix  # noqa: E402


def _metric(index, name, courses, days):
    return {
        "schedule_index": index,
        "schedule_name": name,
        "mean_courses": courses,
        "mean_worst_endpoint_days": days,
    }


def test_selects_adjacent_upper_frontier_schedules_and_interpolates_without_extrapolation():
    metrics = [
        _metric(0, "none", 0.0, 100.0),
        _metric(1, "three", 2.0, 160.0),
        _metric(2, "six", 4.0, 170.0),  # below the upper hull
        _metric(3, "twelve", 6.0, 250.0),
    ]

    match = select_schedule_mix(metrics, 3.0)

    assert match["status"] == "matched"
    assert match["lower_schedule_name"] == "three"
    assert match["upper_schedule_name"] == "twelve"
    assert match["upper_mix_probability"] == pytest.approx(0.25)
    assert match["selection_frontier_mean_worst_endpoint_days"] == pytest.approx(182.5)
    assert [row["schedule_name"] for row in match["frontier_points"]] == ["none", "three", "twelve"]


def test_exact_course_point_uses_best_schedule_at_that_course_level():
    match = select_schedule_mix(
        [
            _metric(0, "a", 1.0, 120.0),
            _metric(1, "b", 1.0, 130.0),
            _metric(2, "c", 3.0, 180.0),
        ],
        1.0,
    )

    assert match["status"] == "matched"
    assert match["lower_schedule_name"] == match["upper_schedule_name"] == "b"
    assert match["upper_mix_probability"] == 0.0
    assert match["selection_frontier_mean_worst_endpoint_days"] == 130.0


@pytest.mark.parametrize("target", [-0.1, 7.0])
def test_refuses_negative_or_out_of_frontier_target_without_extrapolating(target):
    if target < 0:
        with pytest.raises(ValueError, match="non-negative"):
            select_schedule_mix([_metric(0, "a", 0, 10)], target)
    else:
        result = select_schedule_mix(
            [_metric(0, "a", 0, 10), _metric(1, "b", 6, 50)], target
        )
        assert result["status"] == "outside_fixed_schedule_frontier"
        assert result["observed_course_range"] == [0.0, 6.0]
        assert "lower_schedule_index" not in result


def test_randomized_mix_is_seeded_and_uses_no_patient_outcomes():
    low = np.asarray([1, 0, 0, 0], dtype=np.uint8)
    high = np.asarray([1, 1, 1, 1], dtype=np.uint8)

    first, first_meta = randomize_schedule_mix(low, high, 0.4, 20, 177001)
    repeated, repeated_meta = randomize_schedule_mix(low, high, 0.4, 20, 177001)
    changed, _ = randomize_schedule_mix(low, high, 0.4, 20, 177002)

    np.testing.assert_array_equal(first, repeated)
    assert first_meta == repeated_meta
    assert not np.array_equal(first, changed)
    assert first.shape == (20, 4)
    assert first_meta["nominal_upper_schedule_probability"] == 0.4
    assert first_meta["n_assigned_upper_schedule"] == int(first_meta["realized_upper_schedule_fraction"] * 20)
    assert len(first_meta["assignment_sha256"]) == 64


def test_schedule_mixture_rejects_invalid_inputs():
    with pytest.raises(ValueError, match="matching nonempty"):
        randomize_schedule_mix(np.asarray([1]), np.asarray([1, 0]), 0.5, 4, 1)
    with pytest.raises(ValueError, match="probability"):
        randomize_schedule_mix(np.asarray([1]), np.asarray([0]), 1.5, 4, 1)
    with pytest.raises(ValueError, match="positive integer"):
        randomize_schedule_mix(np.asarray([1]), np.asarray([0]), 0.5, 0, 1)
