from __future__ import annotations

import numpy as np
import pytest

from src.trackc_integrated_policy_utils import (
    patient_robust_utility,
    select_schedule_by_plan_utility,
)


def test_selector_uses_mean_patientwise_worst_endpoint_utility():
    # Schedule 0 wins each endpoint on a different patient. Aggregating endpoints
    # before taking the minimum would select it (415 days), but the frozen plan's
    # patient-wise robust score is only 100 days, so the balanced schedule must win.
    ttp_baseline = np.array([[730.0, 100.0], [300.0, 300.0]])
    ttp_rano = np.array([[100.0, 730.0], [300.0, 300.0]])
    courses = np.zeros((2, 2))

    result = select_schedule_by_plan_utility(
        ttp_baseline,
        ttp_rano,
        courses,
        courses,
        ["split_endpoint_wins", "balanced"],
        tau=730.0,
        course_penalty=5.0,
    )

    assert result["schedule_index"] == 1
    assert result["schedule"] == "balanced"
    np.testing.assert_allclose(result["mean_patient_utility"], [100.0, 300.0])


def test_patient_utility_applies_each_endpoint_course_cost_before_worst_case():
    utility = patient_robust_utility(
        ttp_baseline=np.array([500.0, 400.0]),
        ttp_rano=np.array([300.0, 600.0]),
        courses_baseline=np.array([1.0, 0.0]),
        courses_rano=np.array([0.0, 2.0]),
        tau=730.0,
        course_penalty=5.0,
    )

    np.testing.assert_allclose(utility, [300.0, 400.0])


def test_selector_breaks_utility_ties_by_fewer_courses_then_catalog_order():
    times = np.full((3, 2), 500.0)
    courses_baseline = np.array([[3.0, 3.0], [2.0, 2.0], [2.0, 2.0]])
    courses_rano = courses_baseline.copy()

    result = select_schedule_by_plan_utility(
        times,
        times,
        courses_baseline,
        courses_rano,
        ["more_courses", "first_minimum", "second_minimum"],
        tau=730.0,
        course_penalty=0.0,
    )

    assert result["schedule_index"] == 1
    assert result["schedule"] == "first_minimum"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tau": 0.0},
        {"tau": float("nan")},
        {"course_penalty": -1.0},
    ],
)
def test_utility_rejects_invalid_horizon_or_penalty(kwargs):
    args = dict(
        ttp_baseline=[100.0],
        ttp_rano=[100.0],
        courses_baseline=[0.0],
        courses_rano=[0.0],
        tau=730.0,
        course_penalty=5.0,
    )
    args.update(kwargs)
    with pytest.raises(ValueError):
        patient_robust_utility(**args)
