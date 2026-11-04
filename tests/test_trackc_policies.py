"""Cheap dosing rules (src/trackc_policies.py). Synthetic readings; no claim about any patient."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_policies as pol  # noqa: E402
import trackc_sim as tc  # noqa: E402


def test_schedule_shared_mask_gives_the_same_action_to_everyone():
    mask = np.zeros(tc.N_CYCLES, bool)
    mask[[0, 1, 4]] = True
    p = pol.Schedule(mask)
    y = np.zeros(5)
    assert [bool(p.act(j, y)[0]) for j in range(6)] == [True, True, False, False, True, False]
    assert p.act(0, y).shape == (5,) and p.act(0, y).all()


def test_schedule_per_patient_mask_follows_each_row():
    mask = np.zeros((3, tc.N_CYCLES), bool)
    mask[0, 0] = mask[1, 1] = mask[2, 0] = mask[2, 1] = True
    p = pol.Schedule(mask)
    assert p.act(0, np.zeros(3)).tolist() == [True, False, True]
    assert p.act(1, np.zeros(3)).tolist() == [False, True, True]


def test_schedule_matches_periodic_in_the_simulator():
    scen = tc.draw_scenarios(12, 41)
    mask = np.zeros(tc.N_CYCLES, bool)
    mask[::2] = True
    a = tc.simulate(scen, pol.Schedule(mask), [3, 0])
    b = tc.simulate(scen, tc.Periodic(2), [3, 0])
    np.testing.assert_allclose(a["ttp"], b["ttp"], atol=1e-12)
    assert a["courses"].tolist() == b["courses"].tolist()


def test_nadir_containment_decisions_on_a_hand_made_reading_sequence():
    p = pol.NadirContainment(induction=1, a_on=1.25, a_off=1.0)
    p.reset(None, None)
    logs = np.log([1.0, 0.5, 0.45, 0.50, 0.57, 0.50, 0.44])           # one patient; burden relative to baseline
    acts = [bool(p.act(j, np.array([v]))[0]) for j, v in enumerate(logs)]
    # j0 induction: on. j1 gap 0 -> off. j2 new nadir -> off. j3 gap ln(0.50/0.45) = 11% -> off.
    # j4 gap ln(0.57/0.45) = 27% >= 25% -> on. j5 gap 11%, still above a_off = 1.0 -> stays on. j6 new nadir, gap 0 -> off
    assert acts == [True, False, False, False, True, True, False]


def test_nadir_containment_holds_its_decision_when_the_scan_is_missing():
    p = pol.NadirContainment(induction=0, a_on=1.2, a_off=1.0)
    p.reset(None, None)
    p.act(0, np.log(np.array([1.0, 1.0])))
    on = p.act(1, np.log(np.array([1.5, 1.0])))                         # first patient is 50% above the nadir
    assert on.tolist() == [True, False]
    held = p.act(2, np.array([np.nan, np.nan]))
    assert held.tolist() == [True, False]
    held = p.act(3, np.array([np.nan, np.nan]))
    assert held.tolist() == [True, False]


def test_nadir_containment_needs_a_off_not_above_a_on():
    with pytest.raises(ValueError):
        pol.NadirContainment(0, 1.0, 1.2)


def test_nadir_containment_runs_in_the_simulator_with_scans_every_third_cycle():
    scen = tc.draw_scenarios(20, 42)
    scen["scan_every"] = np.full(20, 3)
    out = tc.simulate(scen, pol.NadirContainment(3, 1.25, 1.0), [4, 0], endpoint="rano")
    assert (out["courses"] >= 0).all() and out["courses"].max() <= tc.BUDGET
    assert (out["courses"][out["ttp"] > 3 * tc.CYCLE_DAYS] >= 1).all()      # induction courses were given
    assert np.isfinite(out["ttp"]).all()
