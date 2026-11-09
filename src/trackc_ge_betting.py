"""Fixed-horizon Gaussian-efficient betting interval for bounded means.

This is an independent implementation of the deterministic construction in
Martinez-Taboada and Ramdas (2026), arXiv:2608.21694v2. It is a development
candidate only; it is not wired into the Track C evaluation runner.
"""
from __future__ import annotations

from math import exp, isfinite, log, pi, sqrt
from numbers import Real
from typing import Sequence

import numpy as np
from scipy.special import ndtri


_LOG_SQRT_2PI = 0.5 * log(2.0 * pi)


def _validated_unit_values(values: Sequence[float]) -> np.ndarray:
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size < 1:
        raise ValueError("values must be a nonempty one-dimensional sequence")
    if not np.isfinite(x).all() or (x < 0.0).any() or (x > 1.0).any():
        raise ValueError("values must be finite and within [0, 1]")
    return x


def _validated_delta(delta: Real) -> float:
    if isinstance(delta, (bool, np.bool_)) or not isinstance(delta, Real):
        raise ValueError("delta must be finite and strictly between 0 and 1")
    value = float(delta)
    if not isfinite(value) or not 0.0 < value < 1.0:
        raise ValueError("delta must be finite and strictly between 0 and 1")
    return value


def _variance_path(x: np.ndarray) -> np.ndarray:
    """Return vhat_0,...,vhat_(n-1), shared across candidate means."""
    path = np.empty(x.size, dtype=float)
    running_sum = 0.25
    running_x = 0.0
    previous_mean = 0.5
    path[0] = 0.25
    for index in range(1, x.size):
        observation = float(x[index - 1])
        running_sum += (observation - previous_mean) ** 2
        running_x += observation
        t = index
        path[index] = running_sum / (t + 1)
        previous_mean = (0.5 + running_x) / (t + 1)
    return path


def _psi(p: float) -> float:
    if not 0.0 < p <= 1.0:
        raise ArithmeticError("betting probability left (0, 1]")
    if p == 1.0:
        return 0.0
    z = float(ndtri(p))
    return exp(-0.5 * z * z - _LOG_SQRT_2PI - log(p))


def _tail_wealth(
    x: np.ndarray,
    candidate_mean: float,
    delta: float,
    variance_path: np.ndarray,
    *,
    upper_tail: bool,
) -> float:
    """Compute one terminal wealth, capped at its fixed rejection threshold."""
    threshold = 2.0 / delta
    wealth = 1.0
    n = x.size
    for index, observation in enumerate(x):
        if wealth <= 0.0 or wealth >= threshold:
            break
        p = (delta / 2.0) * wealth
        if p <= 0.0:
            break
        raw_bet = _psi(p) / sqrt((n - index) * variance_path[index])
        if upper_tail:
            cap = float("inf") if candidate_mean == 0.0 else 1.0 / candidate_mean
            increment = float(observation) - candidate_mean
        else:
            cap = float("inf") if candidate_mean == 1.0 else 1.0 / (1.0 - candidate_mean)
            increment = candidate_mean - float(observation)
        bet = min(raw_bet, cap)
        factor = max(0.0, 1.0 + bet * increment)
        wealth = min(threshold, wealth * factor)
    return wealth


def _accepts_mean(x: np.ndarray, candidate_mean: float, delta: float, variance_path: np.ndarray) -> bool:
    threshold = 2.0 / delta
    return (
        _tail_wealth(x, candidate_mean, delta, variance_path, upper_tail=True) < threshold
        and _tail_wealth(x, candidate_mean, delta, variance_path, upper_tail=False) < threshold
    )


def ge_betting_interval_unit(
    values: Sequence[float],
    *,
    delta: float = 0.025,
    tolerance: float = 1e-10,
    max_iterations: int = 80,
) -> tuple[float, float] | None:
    """Return the deterministic GE-betting interval for a mean in [0, 1].

    `delta` is the two-sided noncoverage level for one mean. For two primary
    Track C contrasts, `delta=0.025` per contrast gives a 95% familywise level
    by Bonferroni. Return ``None`` if numerical inversion finds an empty set.
    The deterministic procedure uses the paper's default betting cap c=1 and
    does not use external randomization.
    """
    x = _validated_unit_values(values)
    delta_value = _validated_delta(delta)
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if isinstance(tolerance, bool) or not isinstance(tolerance, Real):
        raise ValueError("tolerance must be finite and positive")
    tolerance_value = float(tolerance)
    if not isfinite(tolerance_value) or tolerance_value <= 0.0:
        raise ValueError("tolerance must be finite and positive")

    threshold = 2.0 / delta_value
    variance_path = _variance_path(x)

    def upper_wealth(mean: float) -> float:
        return _tail_wealth(x, mean, delta_value, variance_path, upper_tail=True)

    def lower_wealth(mean: float) -> float:
        return _tail_wealth(x, mean, delta_value, variance_path, upper_tail=False)

    if upper_wealth(1.0) >= threshold or lower_wealth(0.0) >= threshold:
        return None

    if upper_wealth(0.0) < threshold:
        lower = 0.0
    else:
        lo, hi = 0.0, 1.0
        for _ in range(max_iterations):
            mid = (lo + hi) / 2.0
            if mid == lo or mid == hi or hi - lo <= tolerance_value:
                break
            if upper_wealth(mid) < threshold:
                hi = mid
            else:
                lo = mid
        lower = hi

    if lower_wealth(1.0) < threshold:
        upper = 1.0
    else:
        lo, hi = 0.0, 1.0
        for _ in range(max_iterations):
            mid = (lo + hi) / 2.0
            if mid == lo or mid == hi or hi - lo <= tolerance_value:
                break
            if lower_wealth(mid) < threshold:
                lo = mid
            else:
                hi = mid
        upper = lo

    if lower > upper:
        return None
    midpoint = (lower + upper) / 2.0
    if not _accepts_mean(x, midpoint, delta_value, variance_path):
        return None
    return lower, upper


def ge_betting_interval(
    values: Sequence[float],
    *,
    support: tuple[float, float] = (-790.0, 790.0),
    delta: float = 0.025,
) -> tuple[float, float] | None:
    """Return a GE-betting interval on a finite real support.

    For Track C, the current utility-derived paired-difference support is
    [-790, 790] days. This function only performs the affine transformation;
    it does not validate the utility bound's scientific applicability.
    """
    if len(support) != 2:
        raise ValueError("support must contain (lower, upper)")
    lower_support, upper_support = map(float, support)
    if not isfinite(lower_support) or not isfinite(upper_support) or lower_support >= upper_support:
        raise ValueError("support bounds must be finite and increasing")
    d = np.asarray(values, dtype=float)
    if d.ndim != 1 or d.size < 1 or not np.isfinite(d).all():
        raise ValueError("values must be a nonempty finite one-dimensional sequence")
    if (d < lower_support).any() or (d > upper_support).any():
        raise ValueError("values fall outside the declared support")
    scaled = (d - lower_support) / (upper_support - lower_support)
    interval = ge_betting_interval_unit(scaled, delta=delta)
    if interval is None:
        return None
    width = upper_support - lower_support
    return lower_support + width * interval[0], lower_support + width * interval[1]
