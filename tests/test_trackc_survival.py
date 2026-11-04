"""Survival statistics (src/trackc_survival.py) on synthetic data with a known answer."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_survival as sv  # noqa: E402


def exp_data(n, lam0, lam1, cap, seed, quantum=None):
    rng = np.random.default_rng(seed)
    g = np.r_[np.zeros(n, int), np.ones(n, int)]
    t = rng.exponential(1.0 / np.where(g == 0, lam0, lam1))
    if quantum:
        t = np.ceil(t / quantum) * quantum                 # events seen only at scan visits: many ties
    e = t < cap
    return np.minimum(t, cap), e, g


def test_km_median_matches_empirical_median_without_censoring():
    rng = np.random.default_rng(1)
    t = rng.exponential(100.0, 2001)
    assert sv.km_median(t, np.ones_like(t, bool)) == pytest.approx(float(np.median(t)), rel=1e-3)


def test_interpolated_median_on_a_scan_grid():
    # 100 patients; 30 found at the day-84 scan, 30 at day 168, 40 never: S(84) = 0.7, S(168) = 0.4
    t = np.r_[np.full(30, 84.0), np.full(30, 168.0), np.full(40, 1456.0)]
    e = np.r_[np.ones(60, bool), np.zeros(40, bool)]
    assert sv.km_median(t, e) == 168.0                                   # step median sits on the scan day
    assert sv.km_median_interp(t, e) == pytest.approx(84.0 + 84.0 * (0.7 - 0.5) / (0.7 - 0.4))   # 140 d
    # first gap: S(84) = 0.4 -> crossing inside [0, 84]
    t2 = np.r_[np.full(60, 84.0), np.full(40, 1456.0)]
    e2 = np.r_[np.ones(60, bool), np.zeros(40, bool)]
    assert sv.km_median_interp(t2, e2) == pytest.approx(84.0 * 0.5 / 0.6)
    assert sv.km_median_interp(np.full(10, 5.0), np.zeros(10, bool)) == float("inf")


def test_interpolated_median_agrees_with_step_median_for_fine_data():
    rng = np.random.default_rng(11)
    t = rng.exponential(100.0, 3001)
    ev = np.ones_like(t, bool)
    assert sv.km_median_interp(t, ev) == pytest.approx(sv.km_median(t, ev), rel=2e-3)


def test_km_median_is_infinite_when_curve_stays_above_half():
    t = np.r_[np.full(30, 400.0), np.full(70, 1000.0)]
    e = np.r_[np.ones(30, bool), np.zeros(70, bool)]
    assert sv.km_median(t, e) == float("inf")
    assert sv.km_survival(t, e, 500.0) == pytest.approx(0.7)


def test_km_handles_censoring():
    t = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    e = np.array([1, 0, 1, 0, 1], bool)
    # S(1) = 4/5; at t = 3 two subjects are at risk beyond the censored one: 3 at risk -> S = 0.8 * 2/3
    assert sv.km_survival(t, e, 3.0) == pytest.approx(0.8 * 2.0 / 3.0)


def test_cox_recovers_hazard_ratio_without_ties():
    t, e, g = exp_data(3000, 0.004, 0.008, cap=600.0, seed=2)
    hr, (lo, hi) = sv.cox_hr(t, e, g)
    assert lo < 2.0 < hi and hr == pytest.approx(2.0, rel=0.1)


def test_cox_with_heavy_ties_still_close_with_efron():
    t, e, g = exp_data(3000, 0.004, 0.008, cap=600.0, seed=3, quantum=84.0)
    hr, (lo, hi) = sv.cox_hr(t, e, g)
    assert hr == pytest.approx(2.0, rel=0.15)
    assert np.unique(t[e]).size <= 8                        # ties really are heavy here


def test_cox_equal_groups_gives_hazard_ratio_near_one_and_swap_inverts():
    t, e, g = exp_data(2000, 0.005, 0.005, cap=500.0, seed=4)
    hr, (lo, hi) = sv.cox_hr(t, e, g)
    assert lo < 1.0 < hi
    hr_swapped, _ = sv.cox_hr(t, e, 1 - g)
    assert hr * hr_swapped == pytest.approx(1.0, rel=1e-3)


def test_cox_needs_events_in_both_groups():
    with pytest.raises(ValueError):
        sv.cox_hr(np.array([1.0, 2.0, 3.0, 4.0]), np.array([1, 1, 0, 0], bool), np.array([0, 0, 1, 1]))
