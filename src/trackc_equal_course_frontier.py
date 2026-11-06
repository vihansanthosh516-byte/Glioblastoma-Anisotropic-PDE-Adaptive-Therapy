"""Selection-only fixed-schedule frontier matching for synthetic policy studies.

This module picks a randomized mixture of two fixed schedules at a target mean
course use. The schedule frontier and target dose must come from non-evaluation
data; evaluation outcomes are never inputs to these helpers.
"""
from __future__ import annotations

import math
import hashlib

import numpy as np


_TOL = 1e-9


def select_schedule_mix(schedule_metrics, target_mean_courses: float) -> dict:
    """Select the upper fixed-schedule envelope at a calibration dose level.

    ``schedule_metrics`` rows need ``schedule_index``, ``schedule_name``,
    ``mean_courses`` and ``mean_worst_endpoint_days``. The last quantity is the
    selection-cohort mean of the patient's earlier baseline/RANO-like endpoint,
    each restricted at the same horizon. Mixtures between adjacent upper-hull
    schedules interpolate expected mean course use and expected restricted time.
    Extrapolation beyond the observed schedule range is refused.
    """
    try:
        target = float(target_mean_courses)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("target_mean_courses must be finite and non-negative") from exc
    if not math.isfinite(target) or target < 0:
        raise ValueError("target_mean_courses must be finite and non-negative")

    rows = list(schedule_metrics)
    if not rows:
        raise ValueError("schedule_metrics must contain at least one schedule")

    best_at_x = {}
    seen_indices = set()
    seen_names = set()
    for row in rows:
        try:
            index = int(row["schedule_index"])
            name = row["schedule_name"]
            courses = float(row["mean_courses"])
            days = float(row["mean_worst_endpoint_days"])
        except (KeyError, TypeError, ValueError, OverflowError) as exc:
            raise ValueError("each schedule metric needs index, name, course use, and endpoint days") from exc
        if isinstance(row["schedule_index"], (bool, np.bool_)) or index != row["schedule_index"] or index < 0:
            raise ValueError("schedule_index must be a non-negative integer")
        if not isinstance(name, str) or not name.strip() or name in seen_names:
            raise ValueError("schedule names must be nonempty and unique")
        if index in seen_indices:
            raise ValueError("schedule indices must be unique")
        if not math.isfinite(courses) or courses < 0 or not math.isfinite(days):
            raise ValueError("schedule course use and endpoint days must be finite; course use non-negative")
        seen_indices.add(index)
        seen_names.add(name)
        x_key = courses
        candidate = {"x": x_key, "y": days, "index": index, "name": name}
        incumbent = best_at_x.get(x_key)
        if incumbent is None or days > incumbent["y"] or (
            days == incumbent["y"] and index < incumbent["index"]
        ):
            best_at_x[x_key] = candidate

    ordered = [best_at_x[x] for x in sorted(best_at_x)]
    hull = []
    for point in ordered:
        while len(hull) >= 2:
            left, middle = hull[-2], hull[-1]
            left_slope_num = (middle["y"] - left["y"]) * (point["x"] - middle["x"])
            right_slope_num = (point["y"] - middle["y"]) * (middle["x"] - left["x"])
            if left_slope_num <= right_slope_num:
                hull.pop()
            else:
                break
        hull.append(point)

    serialized_hull = [
        {
            "schedule_index": int(point["index"]),
            "schedule_name": point["name"],
            "mean_courses": float(point["x"]),
            "mean_worst_endpoint_days": float(point["y"]),
        }
        for point in hull
    ]
    result = {
        "target_mean_courses": target,
        "frontier_points": serialized_hull,
        "selection_rule": "upper concave hull of fixed schedules; linear randomized mixture; no extrapolation",
    }
    if target < hull[0]["x"] - _TOL or target > hull[-1]["x"] + _TOL:
        result.update({
            "status": "outside_fixed_schedule_frontier",
            "observed_course_range": [float(hull[0]["x"]), float(hull[-1]["x"])],
        })
        return result

    for point in hull:
        if abs(target - point["x"]) <= _TOL:
            result.update({
                "status": "matched",
                "lower_schedule_index": int(point["index"]),
                "lower_schedule_name": point["name"],
                "upper_schedule_index": int(point["index"]),
                "upper_schedule_name": point["name"],
                "upper_mix_probability": 0.0,
                "selection_frontier_mean_worst_endpoint_days": float(point["y"]),
            })
            return result

    for lower, upper in zip(hull, hull[1:]):
        if lower["x"] < target < upper["x"]:
            weight = (target - lower["x"]) / (upper["x"] - lower["x"])
            frontier_days = lower["y"] + weight * (upper["y"] - lower["y"])
            result.update({
                "status": "matched",
                "lower_schedule_index": int(lower["index"]),
                "lower_schedule_name": lower["name"],
                "upper_schedule_index": int(upper["index"]),
                "upper_schedule_name": upper["name"],
                "upper_mix_probability": float(weight),
                "selection_frontier_mean_worst_endpoint_days": float(frontier_days),
            })
            return result

    raise ArithmeticError("target was within the hull range but no matching segment was found")


def randomize_schedule_mix(lower_schedule, upper_schedule, upper_probability: float, n_patients: int, seed: int):
    """Assign a seeded independent schedule draw to each patient, without outcome input."""
    lower = np.asarray(lower_schedule)
    upper = np.asarray(upper_schedule)
    if lower.ndim != 1 or upper.shape != lower.shape or lower.size == 0:
        raise ValueError("lower and upper schedules must be matching nonempty one-dimensional arrays")
    if not np.isin(lower, (False, True, 0, 1)).all() or not np.isin(upper, (False, True, 0, 1)).all():
        raise ValueError("schedule entries must be binary")
    try:
        probability = float(upper_probability)
        n = int(n_patients)
        seed_value = int(seed)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("probability, patient count, and seed must be valid numeric values") from exc
    if not math.isfinite(probability) or not 0 <= probability <= 1:
        raise ValueError("upper_probability must be finite and in [0, 1]")
    if isinstance(n_patients, (bool, np.bool_)) or n != n_patients or n < 1:
        raise ValueError("n_patients must be a positive integer")
    if isinstance(seed, (bool, np.bool_)) or seed_value != seed or seed_value < 0:
        raise ValueError("seed must be a non-negative integer")

    assigned_upper = np.random.default_rng(seed_value).random(n) < probability
    schedules = np.broadcast_to(lower.astype(bool), (n, lower.size)).copy()
    schedules[assigned_upper] = upper.astype(bool)
    return schedules, {
        "n_patients": n,
        "n_assigned_upper_schedule": int(assigned_upper.sum()),
        "realized_upper_schedule_fraction": float(assigned_upper.mean()),
        "nominal_upper_schedule_probability": probability,
        "assignment_seed": seed_value,
        "assignment_sha256": hashlib.sha256(
            assigned_upper.astype(np.uint8).tobytes()
        ).hexdigest(),
    }


__all__ = ["randomize_schedule_mix", "select_schedule_mix"]
