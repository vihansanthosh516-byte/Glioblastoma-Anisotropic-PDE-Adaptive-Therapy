"""Regime-map helpers (src/trackc_regime.py). Synthetic checks only."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_regime as rg  # noqa: E402
import trackc_sim as tc  # noqa: E402


def test_candidate_schedules_are_unique_short_and_include_the_standard_arm():
    sched, names = rg.candidate_schedules()
    assert sched.shape == (len(names), tc.N_CYCLES) and len(set(names)) == len(names)
    assert sched.sum(1).max() <= tc.BUDGET and sched.sum(1).min() == 0          # front0 is no treatment
    assert len({tuple(r) for r in sched.tolist()}) == len(sched)                 # no duplicate rows
    soc = np.zeros(tc.N_CYCLES, bool)
    soc[:6] = True
    assert any((r == soc).all() for r in sched) and "front6" in names
    assert 40 <= len(names) <= 120


def test_rmst_and_paired_stats():
    ttp = np.array([100.0, 500.0, 900.0, 1456.0])
    assert rg.rmst(ttp, 365.0) == pytest.approx((100 + 365 + 365 + 365) / 4)
    ref = np.array([100.0, 400.0, 800.0, 1000.0])
    st = rg.paired_stats(ttp, ref, 730.0, alpha=0.25)
    d = np.array([0.0, 100.0, 0.0, 0.0])                                         # capped at 730: 100,500,730,730 vs 100,400,730,730
    assert st["diff"] == pytest.approx(d.mean()) and st["frac_worse"] == 0.0
    assert st["regret_cvar"] == pytest.approx(0.0) and st["ci_lo"] < st["diff"] < st["ci_hi"]
    worse = rg.paired_stats(ref, ttp, 730.0, alpha=0.25)
    assert worse["diff"] == pytest.approx(-25.0) and worse["regret_cvar"] == pytest.approx(-100.0)


def test_evaluate_schedules_matches_one_at_a_time_runs():
    scen = tc.draw_scenarios(15, 51)
    sched, names = rg.candidate_schedules()
    pick = [names.index(n) for n in ("front0", "front6", "every3_from0", "block2_gap4_block4")]
    for endpoint in ("rano", "baseline"):
        ttp, courses = rg.evaluate_schedules(scen, sched[pick], endpoint)
        for row, i in enumerate(pick):
            one = dict(scen, scan_every=np.full(15, rg.SCAN_EVERY))
            from trackc_policies import Schedule
            ref = tc.simulate(one, Schedule(sched[i]), [0, 0], endpoint=endpoint)
            np.testing.assert_allclose(ttp[row], ref["ttp"], atol=1e-9)
            assert courses[row].tolist() == ref["courses"].tolist()
    assert (ttp[0] <= ttp[1]).mean() > 0.9                                        # treating beats not treating


def test_untreated_arm_has_the_closed_form_rano_time():
    scen = tc.draw_scenarios(10, 52)
    scen["c"] = np.zeros(10)
    sched, names = rg.candidate_schedules()
    ttp, _ = rg.evaluate_schedules(scen, sched[[names.index("front0")]], "rano")
    t_cross = np.log(tc.RANO_FACTOR * (1 - scen["B0"]) / (1 - tc.RANO_FACTOR * scen["B0"])) / scen["r"]
    step = tc.CYCLE_DAYS * rg.SCAN_EVERY
    np.testing.assert_allclose(ttp[0], np.minimum(np.ceil(t_cross / step) * step, tc.HORIZON), atol=1e-6)


def test_policy_arms_run_and_respect_their_budgets():
    scen = tc.draw_scenarios(10, 53)
    arms = rg.policy_arms()
    assert {"soc6", "front12", "nadir_k3_on1.25"} <= set(arms)
    for name in ("soc6", "front3", "raw_0.6_0.6", "nadir_k3_on1.25"):
        factory, budget = arms[name]
        ttp, courses = rg.run_policy(scen, factory, budget, "rano", [6, 0])
        assert courses.max() <= budget and np.isfinite(ttp).all()
