"""Stage 2 world family and standard arm (src/trackc_calibration.py). Synthetic parameters only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_sim as tc  # noqa: E402

OLD_WORLD = {"world": 0, "log_r_shift": 0.0, "log_m_scale": 0.0, "log_f0_shift": 0.0, "B0_center": 0.275,
             "c": 0.1, "mu": 0.0, "nu": 0.0}


def test_target_is_the_published_interval_in_days():
    assert cal.TARGET_PFS_DAYS[0] == pytest.approx(143.07, abs=0.01)     # 4.7 months
    assert cal.TARGET_PFS_DAYS[1] == pytest.approx(185.68, abs=0.01)     # 6.1 months
    assert cal.in_target(168.0) and not cal.in_target(252.0) and not cal.in_target(84.0)


def test_draw_worlds_is_deterministic_stratified_and_inside_the_axes():
    a, b = cal.draw_worlds(40, 5), cal.draw_worlds(40, 5)
    assert a == b and a != cal.draw_worlds(40, 6)
    for k, (name, lo, hi) in enumerate(cal.AXES):
        vals = np.array([w[name] for w in a])
        assert (vals >= lo).all() and (vals <= hi).all()
        bins = np.floor((vals - lo) / (hi - lo) * 40).astype(int).clip(0, 39)
        assert sorted(bins.tolist()) == list(range(40)), name             # Latin hypercube: one point per stratum
    assert {w["mu"] for w in a} <= set(cal.MU_NU_LEVELS) and {w["nu"] for w in a} <= set(cal.MU_NU_LEVELS)
    counts = [sum(w["mu"] == lvl for w in a) for lvl in cal.MU_NU_LEVELS]
    assert max(counts) - min(counts) <= 2                                  # levels used with near-equal weight


def test_world_ranges_at_the_old_setting_are_the_old_ranges():
    rr = cal.world_ranges(OLD_WORLD)
    for name in ("r", "m", "f0"):
        assert rr[name][0] == pytest.approx(tc.RANGES[name][0]) and rr[name][1] == pytest.approx(tc.RANGES[name][1])
    assert rr["B0"][0] == pytest.approx(0.175) and rr["B0"][1] == pytest.approx(0.375)
    assert rr["c"][0] == rr["c"][1] == 0.1
    assert rr["eps"] == tc.RANGES["eps"]


def test_shifts_move_the_ranges_as_documented():
    w = dict(OLD_WORLD, log_r_shift=np.log(2.0), log_m_scale=np.log(0.1), log_f0_shift=np.log(10.0))
    rr = cal.world_ranges(w)
    assert rr["r"][0] == pytest.approx(0.002) and rr["r"][1] == pytest.approx(0.1)
    assert rr["m"][0] == pytest.approx(0.6) and rr["m"][1] == pytest.approx(2.0)
    assert rr["f0"][0] == pytest.approx(0.01) and rr["f0"][1] == pytest.approx(cal.F0_CAP)    # 1.0 capped at 0.9


def test_patients_stay_inside_their_world():
    w = dict(OLD_WORLD, log_r_shift=0.5, B0_center=0.4, c=0.2, mu=0.01, nu=0.04)
    scen = cal.draw_patients(w, 300, 7)
    rr = cal.world_ranges(w)
    for name in ("r", "m", "f0", "B0"):
        assert (scen[name] >= rr[name][0] - 1e-12).all() and (scen[name] <= rr[name][1] + 1e-12).all(), name
    assert (scen["c"] == 0.2).all() and (scen["mu"] == 0.01).all() and (scen["nu"] == 0.04).all()
    assert len(scen["r"]) == 300 and (scen["f0"] <= cal.F0_CAP).all()


def test_standard_arm_detections_fall_on_scan_days_and_courses_are_capped():
    scen = cal.draw_patients(OLD_WORLD, 200, 8)
    ttp, event, courses = cal.run_arm(scen)
    assert courses.max() <= cal.STD_COURSES
    det = ttp[event]
    assert event.any() and np.allclose(det % (tc.CYCLE_DAYS * cal.SCAN_EVERY), 0.0)
    assert (ttp[~event] == tc.HORIZON).all()
    assert (courses[event & (ttp <= 84)] <= 3).all()                      # called at the first scan: at most 3 courses


def test_batching_worlds_does_not_change_any_patient():
    wa, wb = OLD_WORLD, dict(OLD_WORLD, world=1, log_r_shift=0.6, log_m_scale=-1.0, mu=0.01, nu=0.01, c=0.2)
    a, b = cal.draw_patients(wa, 60, 21), cal.draw_patients(wb, 60, 22)
    alone = cal.run_arm(a)[0], cal.run_arm(b)[0]
    both = cal.run_arm(cal.concat([a, b]))[0]
    np.testing.assert_allclose(both[:60], alone[0], atol=1e-9)
    np.testing.assert_allclose(both[60:], alone[1], atol=1e-9)


def test_zero_switching_arrays_leave_the_dynamics_unchanged():
    scen = tc.draw_scenarios(30, 3)
    plain = cal.run_arm(scen)[0]
    with_zero = cal.run_arm(dict(scen, mu=np.zeros(30), nu=np.zeros(30)))[0]
    np.testing.assert_allclose(with_zero, plain, atol=1e-12)


def test_control_flags_use_the_published_intervals():
    assert cal.control_flags(0.87, 0.99) == (True, True)
    assert cal.control_flags(0.60, 1.50) == (False, False)
    assert cal.control_flags(0.75, 1.38) == (True, True) and cal.control_flags(1.01, 0.70) == (False, False)
    assert not any(np.isnan(x) for x in (cal.DOSE_DENSE_HR + cal.TWELVE_VS_SIX_HR))
    assert cal.control_flags(float("nan"), float("nan")) == (False, False)


def test_dose_dense_control_on_the_old_ranges_shows_a_benefit():
    scen = tc.draw_scenarios(600, 31)                                      # old model (a pilot gave HR about 0.8 at x2)
    out = cal.dose_dense_control(scen, mult=2.0)
    assert out["hr"] < 0.95                                                # doubling the drug effect lengthens PFS here
    assert out["median_dose_dense"] > out["median_standard"]
    same = cal.dose_dense_control(scen, mult=1.0)
    assert same["hr"] == pytest.approx(1.0, abs=1e-6)                      # multiplier 1 is the standard arm


def test_landmark_control_cohort_and_timing():
    w = dict(OLD_WORLD, log_r_shift=0.0, log_m_scale=np.log(0.15), log_f0_shift=0.0)
    scen = cal.draw_patients(w, 600, 32)
    out = cal.landmark_control(scen)
    ttp0, ev0, _ = cal.run_arm(scen)
    assert out["n_landmark"] == int((ttp0 > cal.LANDMARK_DAYS).sum()) and out["n_total"] == 600
    assert 0 < out["n_landmark"] <= 600 and np.isfinite(out["hr"])
    # no extra drug can be given to a patient who already had a progression call: HR is finite and positive
    assert out["hr"] > 0 and out["ci"][0] < out["hr"] < out["ci"][1]


def test_hazard_ratio_returns_nan_for_an_arm_without_events():
    t = np.array([100.0, 200.0, 300.0])
    hr, ci = cal.hazard_ratio(t, np.array([1, 1, 1], bool), t, np.zeros(3, bool))
    assert np.isnan(hr) and np.isnan(ci[0])


def test_pfs_stats_reports_both_medians_and_survival_after_each_scan_day():
    t = np.r_[np.full(30, 84.0), np.full(30, 168.0), np.full(40, 1456.0)]
    e = np.r_[np.ones(60, bool), np.zeros(40, bool)]
    st = cal.pfs_stats(t, e)
    assert st["median_step"] == 168.0 and st["median_interp"] == pytest.approx(140.0)
    assert st["S84"] == pytest.approx(0.7) and st["S168"] == pytest.approx(0.4) and st["event_fraction"] == pytest.approx(0.6)
    assert cal.in_target(st["median_interp"]) is False                    # 140 d is just under 143 d
