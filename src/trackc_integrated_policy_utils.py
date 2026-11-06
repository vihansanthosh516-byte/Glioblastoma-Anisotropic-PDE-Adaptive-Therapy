"""Selection utilities frozen in the integrated Track A/B/C research plan.

These functions operate on development-world outcomes only. They do not fit a
therapy model, access patient data, or evaluate held-out worlds.
"""
from __future__ import annotations

import numpy as np


def patient_robust_utility(
    ttp_baseline,
    ttp_rano,
    courses_baseline,
    courses_rano,
    *,
    tau=730.0,
    course_penalty=5.0,
):
    """Return each patient's worse endpoint-specific restricted utility.

    Times are restricted at ``tau`` before the endpoint-specific course costs
    are applied. Inputs may be scalars or broadcastable arrays.
    """
    tau = _finite_nonnegative(tau, "tau", strictly_positive=True)
    course_penalty = _finite_nonnegative(course_penalty, "course_penalty")
    tb, tr, cb, cr = np.broadcast_arrays(
        np.asarray(ttp_baseline, dtype=float),
        np.asarray(ttp_rano, dtype=float),
        np.asarray(courses_baseline, dtype=float),
        np.asarray(courses_rano, dtype=float),
    )
    if not all(np.isfinite(values).all() for values in (tb, tr, cb, cr)):
        raise ValueError("times and course counts must be finite")
    if any(np.any(values < 0) for values in (tb, tr, cb, cr)):
        raise ValueError("times and course counts must be non-negative")

    utility_baseline = np.minimum(tb, tau) - course_penalty * cb
    utility_rano = np.minimum(tr, tau) - course_penalty * cr
    return np.minimum(utility_baseline, utility_rano)


def select_schedule_by_plan_utility(
    ttp_baseline,
    ttp_rano,
    courses_baseline,
    courses_rano,
    schedule_names,
    *,
    tau=730.0,
    course_penalty=5.0,
):
    """Select on training patients by mean patient-wise worst-endpoint utility.

    Arrays must be schedule-major with shape ``(n_schedules, n_patients)``.
    Exact utility ties go to fewer mean courses averaged across endpoints, then
    to the earlier schedule in ``schedule_names``.
    """
    arrays = tuple(
        np.asarray(values, dtype=float)
        for values in (ttp_baseline, ttp_rano, courses_baseline, courses_rano)
    )
    if any(values.ndim != 2 for values in arrays):
        raise ValueError("schedule outcomes must be two-dimensional")
    if any(values.shape != arrays[0].shape for values in arrays[1:]):
        raise ValueError("all schedule outcome arrays must have the same shape")
    n_schedules, n_patients = arrays[0].shape
    if n_schedules < 1 or n_patients < 1:
        raise ValueError("schedule outcomes must include schedules and patients")
    names = tuple(str(name) for name in schedule_names)
    if len(names) != n_schedules or len(set(names)) != len(names):
        raise ValueError("schedule_names must be unique and match the schedule count")

    utility = patient_robust_utility(
        *arrays, tau=tau, course_penalty=course_penalty
    )
    mean_patient_utility = utility.mean(axis=1)
    mean_courses = 0.5 * (arrays[2].mean(axis=1) + arrays[3].mean(axis=1))
    order = np.lexsort((np.arange(n_schedules), mean_courses, -mean_patient_utility))
    best = int(order[0])
    return {
        "schedule_index": best,
        "schedule": names[best],
        "mean_patient_utility": mean_patient_utility,
        "mean_courses": mean_courses,
        "n_training_patients": int(n_patients),
        "selection_rule": (
            "mean patient-wise minimum of endpoint-specific restricted time minus "
            "endpoint-specific course cost; ties prefer fewer mean courses, then catalog order"
        ),
    }


def _finite_nonnegative(value, name, *, strictly_positive=False):
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite and non-negative") from exc
    if not np.isfinite(number) or number < 0 or (strictly_positive and number == 0):
        qualifier = "positive" if strictly_positive else "non-negative"
        raise ValueError(f"{name} must be finite and {qualifier}")
    return number


__all__ = ["patient_robust_utility", "select_schedule_by_plan_utility"]
