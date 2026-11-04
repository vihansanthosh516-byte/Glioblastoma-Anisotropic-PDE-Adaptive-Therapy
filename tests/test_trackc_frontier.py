"""Script 202 (feedback rules against the fixed-schedule frontier) and the ResponseGuided rule with its regime-map arms.
Synthetic rows and synthetic patients only; no result file is read."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_policies as pol  # noqa: E402
import trackc_regime as rg  # noqa: E402


def _load():
    spec = importlib.util.spec_from_file_location("stage2_202", REPO / "src" / "202_trackc_frontier_analysis.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def _row(world, arm, courses, rmst, ep="3b_rano", tau=730):
    return {"world": world, "endpoint": ep, "arm": arm, "tau": tau, "rmst": rmst, "mean_courses": courses}


# ---------------------------------------------------------------- script 202
def test_upper_hull_keeps_the_concave_envelope():
    s = _load()
    hull = s.upper_hull([(0, 100), (2, 110), (3, 108), (4, 114), (6, 115), (2, 105)])
    assert hull == [(0.0, 100.0), (2.0, 110.0), (4.0, 114.0), (6.0, 115.0)]
    assert s.frontier_at(hull, 3.0) == pytest.approx(112.0)
    assert np.isnan(s.frontier_at(hull, 7.0)) and np.isnan(s.frontier_at(hull, -0.5))
    assert np.isnan(s.frontier_at([(1.0, 5.0)], 1.0))


def test_family_names():
    s = _load()
    assert [s.family(a) for a in ("soc6", "std_capmix", "front12", "best_fixed_lam5", "paced3", "raw_0.6_0.6",
                                  "nadir_k3_on1.25", "resp_k3_th0.75", "informed_lam10", "informed_family_max")] == \
        ["fixed"] * 5 + ["feedback"] * 3 + ["informed"] * 2


def _synthetic():
    rows = []
    for w in range(6):
        for arm, c, r in (("front3", 2.0, 110.0), ("front9", 4.0, 114.0), ("front12", 6.0, 115.0), ("soc6", 5.0, 100.0)):
            rows.append(_row(w, arm, c, r))
        rows.append(_row(w, "resp_k3_th0.75", 3.0, 113.0))        # frontier at 3 courses is 112: +1 d
        rows.append(_row(w, "nadir_k3_on1.25", 3.0, 111.0))       # -1 d
        rows.append(_row(w, "raw_0.6_0.6", 4.0, 114.0))           # on the frontier
        rows.append(_row(w, "informed_family_max", 3.0, 120.0))   # +8 d headroom
        rows.append(_row(w, "nadir_k6_on1.25", 9.0, 130.0))       # outside the fixed range
    return rows


def test_gaps_and_summary_on_synthetic_rows():
    s = _load()
    gaps = {(g["world"], g["arm"]): g for g in s.gaps_for(_synthetic()) if g["tau"] == 730}
    assert gaps[(0, "resp_k3_th0.75")]["gap"] == pytest.approx(1.0) and gaps[(0, "nadir_k3_on1.25")]["gap"] == pytest.approx(-1.0)
    assert gaps[(0, "raw_0.6_0.6")]["gap"] == pytest.approx(0.0) and gaps[(0, "informed_family_max")]["gap"] == pytest.approx(8.0)
    assert np.isnan(gaps[(0, "nadir_k6_on1.25")]["gap"])
    assert not any(g["arm"] == "front12" for g in s.gaps_for(_synthetic()))            # fixed arms are not scored
    res = s.analyse(_synthetic())
    arms = res["endpoints"]["3b_rano"]["730"]["arms"]
    assert res["n_worlds"] == 6
    assert arms["resp_k3_th0.75"]["n_above_frontier"] == 6 and arms["resp_k3_th0.75"]["median_gap_days"] == pytest.approx(1.0)
    assert arms["resp_k3_th0.75"]["sign_test_p_two_sided"] == pytest.approx(2 * 0.5 ** 6)
    assert arms["nadir_k3_on1.25"]["n_below_frontier"] == 6 and arms["raw_0.6_0.6"]["n_on_frontier"] == 6
    assert arms["raw_0.6_0.6"]["sign_test_p_two_sided"] is None
    assert arms["nadir_k6_on1.25"]["n_worlds"] == 0 and arms["nadir_k6_on1.25"]["n_outside_fixed_range"] == 6
    best = res["endpoints"]["3b_rano"]["730"]["best_feedback_arm_per_world"]
    assert best["most_often_best_arm"] == "resp_k3_th0.75" and best["median_gap_days"] == pytest.approx(1.0)
    assert arms["informed_family_max"]["family"] == "informed" and arms["informed_family_max"]["median_gap_days"] == pytest.approx(8.0)


def test_write_gaps_creates_a_csv(tmp_path):
    s = _load()
    s.write_gaps(_synthetic(), tmp_path)
    lines = (tmp_path / "gaps_by_world.csv").read_text(encoding="utf-8").splitlines()
    assert lines[0].split(",") == ["world", "endpoint", "tau", "arm", "family", "mean_courses", "rmst", "frontier", "gap"]
    assert len(lines) == 1 + 6 * 5                                   # five scored arms per world, one endpoint, one tau


# ---------------------------------------------------------------- ResponseGuided
def _fresh(induction=3, theta=0.75):
    p = pol.ResponseGuided(induction, theta)
    p.reset({}, None)
    return p


NAN2 = np.array([np.nan, np.nan])


def test_induction_courses_are_given_to_everyone():
    p = _fresh()
    assert p.act(0, np.log([1.0, 1.0])).all() and p.act(1, NAN2).all() and p.act(2, NAN2).all()


def test_a_failed_scan_stops_the_drug_for_good_and_between_scans_the_decision_stands():
    p = _fresh()
    p.act(0, np.log([1.0, 1.0]))
    p.act(1, NAN2)
    p.act(2, NAN2)
    assert p.act(3, np.log([0.5, 0.9])).tolist() == [True, False]        # patient 1 is above 0.75 x baseline
    assert p.act(4, NAN2).tolist() == [True, False] and p.act(5, NAN2).tolist() == [True, False]
    assert p.act(6, np.log([0.8, 0.4])).tolist() == [False, False]       # patient 0 fails now; patient 1 never comes back


def test_a_scan_exactly_at_theta_counts_as_a_response_and_a_missing_baseline_keeps_going():
    p = _fresh()
    p.act(0, np.array([0.0, np.nan]))                                    # patient 1 has no first reading
    p.act(1, NAN2)
    p.act(2, NAN2)
    out = p.act(3, np.array([np.log(0.75), np.log(5.0)]))
    assert out.tolist() == [True, True]


def test_theta_must_be_positive_and_the_arms_are_registered():
    with pytest.raises(ValueError):
        pol.ResponseGuided(3, 0.0)
    arms = rg.policy_arms()
    for name in ("resp_k3_th0.5", "resp_k3_th0.75", "resp_k3_th1.0", "resp_k6_th0.75"):
        assert name in arms and arms[name][1] == 12 and arms[name][0]().uses_filter is False


def test_weak_kill_patients_stop_early_and_strong_kill_patients_continue():
    weak = {"world": 0, "log_r_shift": 0.0, "log_m_scale": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
            "c": 0.05, "mu": 0.0, "nu": 0.0}
    strong = dict(weak, log_m_scale=0.0)
    factory, budget = rg.policy_arms()["resp_k3_th0.75"]
    used = {}
    for name, w in (("weak", weak), ("strong", strong)):
        scen = cal.draw_patients(w, 300, [11, 1])
        _, courses = rg.run_policy(scen, factory, budget, "rano", [3, 4])
        used[name] = float(np.mean(courses))
    assert used["weak"] <= 3.5                                           # nobody responds: about the three induction courses
    assert used["strong"] > used["weak"] + 1.0                           # responders keep being treated
