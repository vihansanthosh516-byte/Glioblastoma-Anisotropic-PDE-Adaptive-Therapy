"""Tests for observation-based treatment stopping in synthetic Track C simulations."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_sim as tc  # noqa: E402


class AlwaysTreat:
    uses_filter = False

    def reset(self, context, rng):
        self.reset_context = context

    def act(self, cycle, readings):
        return np.ones(readings.shape, dtype=bool)

    def advance(self, cycle, treatment):
        pass


def _world(scan_every=1, *, prog_factor=10.0, missed=False):
    scen = tc.draw_scenarios(1, seed=8102, sigmas=(1e-10,))
    scen.update(
        r=np.array([0.03]),
        m=np.array([10.0]),
        f0=np.array([0.01]),
        c=np.array([0.0]),
        eps=np.array([0.0]),
        B0=np.array([0.1]),
        sigma=np.array([1e-10]),
        scan_every=np.array([scan_every]),
        dose=np.array([0.0]),
        prog_factor=np.array([prog_factor]),
    )
    if missed:
        scen["miss_p"] = np.array([1.0])
        scen["miss_visible"] = np.array([True])
    return scen


def test_observed_stop_waits_for_two_positive_scheduled_scans():
    for scan_every in (1, 3):
        out = tc.simulate_two_endpoints(
            _world(scan_every), AlwaysTreat(), seed=[19, 0],
            stop_endpoint="observed_two_scan", budget=12,
        )
        interval = scan_every * tc.CYCLE_DAYS
        assert np.isfinite(out["observed_stop_day"]).all()
        assert out["observed_stop_day"][0] == out["ttp_rano"][0] + interval
        assert out["courses"][0] == out["observed_stop_day"][0] / tc.CYCLE_DAYS


def test_missing_scheduled_scans_prevent_confirmation():
    out = tc.simulate_two_endpoints(
        _world(missed=True), AlwaysTreat(), seed=[20, 0],
        stop_endpoint="observed_two_scan", budget=4,
    )
    assert np.isnan(out["observed_stop_day"]).all()
    assert out["courses"].tolist() == [4]


def test_hidden_progression_threshold_does_not_control_observed_stop():
    low_threshold = tc.simulate_two_endpoints(
        _world(prog_factor=1.2), AlwaysTreat(), seed=[21, 0],
        stop_endpoint="observed_two_scan", budget=12,
    )
    high_threshold = tc.simulate_two_endpoints(
        _world(prog_factor=10.0), AlwaysTreat(), seed=[21, 0],
        stop_endpoint="observed_two_scan", budget=12,
    )
    assert np.array_equal(low_threshold["courses"], high_threshold["courses"])
    assert np.array_equal(low_threshold["observed_stop_day"], high_threshold["observed_stop_day"])
    assert low_threshold["ttp_baseline"][0] < high_threshold["ttp_baseline"][0]


def test_shared_observed_stop_state_resets_only_on_scheduled_missing_scans():
    baseline, nadir, count, stopped, newly = tc.update_observed_two_scan_state(
        readings=np.array([0.0]),
        scheduled=True,
        observed_baseline_log=np.array([np.nan]),
        observed_nadir_log=np.array([np.nan]),
        positive_scan_count=np.array([0]),
        stopped=np.array([False]),
    )
    assert baseline[0] == nadir[0] == 0.0
    assert count.tolist() == [0]
    assert not stopped[0] and not newly[0]

    first_positive = np.log(1.5)
    baseline, nadir, count, stopped, newly = tc.update_observed_two_scan_state(
        first_positive, True, baseline, nadir, count, stopped
    )
    assert count.tolist() == [1]
    assert not stopped[0] and not newly[0]

    # A NaN at an unscheduled cycle does not change the confirmation count.
    _, _, count, stopped, _ = tc.update_observed_two_scan_state(
        np.nan, False, baseline, nadir, count, stopped
    )
    assert count.tolist() == [1]

    # A missing scheduled scan resets the two-scan confirmation sequence.
    _, _, count, stopped, _ = tc.update_observed_two_scan_state(
        np.nan, True, baseline, nadir, count, stopped
    )
    assert count.tolist() == [0]
    assert not stopped[0]
