"""Script 204 (what the RANO-like headroom consists of). Synthetic arrays and tiny simulated worlds; no result file is read."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_regime as rg  # noqa: E402


def _load(file="204_trackc_endpoint_mechanism.py", name="stage2_204"):
    spec = importlib.util.spec_from_file_location(name, REPO / "src" / file)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


WEAK = {"world": 0, "log_r_shift": 0.0, "log_m_scale": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
        "c": 0.05, "mu": 0.0, "nu": 0.0}
MIX = {"world": 1, "log_r_shift": 0.0, "log_m_scale_non": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
       "c": 0.05, "q": 0.4, "mu": 0.0, "nu": 0.0}


# ---------------------------------------------------------------- the burden run reproduces trackc_sim.simulate
def test_burden_run_reproduces_the_progression_times_and_courses_of_simulate():
    s = _load()
    sched, names = rg.candidate_schedules()
    for world, seed in ((WEAK, [5, 1]), (MIX, [6, 1])):
        scen = cal.draw_patients(world, 24, seed)
        pick = [names.index(n) for n in ("front0", "front3", "front6", "front12", "every2_from0", "every3_from1", "every6_from2")]
        for scan_every in (rg.SCAN_EVERY, 1):
            ttp_ref, courses_ref = rg.evaluate_schedules(scen, sched[pick], "rano", scan_every=scan_every)
            got = s.evaluate_schedules_burden(scen, sched[pick], scan_every=scan_every)
            assert np.allclose(ttp_ref, got["ttp"]) and (courses_ref == got["courses"]).all()
        assert got["burden"].shape == (len(pick), 24, s.N_REC)
        assert np.allclose(got["burden"][:, :, 0], scen["B0"][None, :]) and (got["burden"] > 0).all()


def test_the_reference_used_at_a_call_is_the_scanned_nadir_with_the_floor():
    s = _load()
    sched, names = rg.candidate_schedules()
    scen = cal.draw_patients(WEAK, 30, [7, 1])
    got = s.evaluate_schedules_burden(scen, sched[[names.index("front12")]])
    called = got["ttp"][0] <= 728.0
    assert called.any()
    ref = got["ref_at_call"][0]
    assert np.isnan(ref[~called]).all() and np.isfinite(ref[called]).all()
    assert (ref[called] <= scen["B0"][called] + 1e-12).all()                       # the reference never exceeds the baseline
    assert (ref[called] >= cal.tc.RANO_FLOOR * scen["B0"][called] - 1e-12).all()   # and never goes below the floor


# ---------------------------------------------------------------- measures on a synthetic burden path
def _paths():
    k = 27
    a = np.ones(k)                                                    # never moves: never called, never above a level
    b = np.r_[np.ones(5), np.full(k - 5, 1.3)]                        # 1.3 from boundary 5 (day 140): above 1.2, below 1.4; called at 6
    c = np.ones(k)
    c[3:6] = 0.5                                                      # a scanned dip to 0.5, regrowth to 0.8 at boundary 6
    c[6:] = 0.8
    d = np.r_[np.ones(4), np.full(4, 1.25), np.full(k - 8, 1.5)]      # 1.25 from boundary 4, 1.5 from boundary 8
    burden = np.vstack([a, b, c, d]) * 2.0                            # baseline 2.0 for all four
    ttp = np.array([1456.0, 168.0, 168.0, 252.0])
    return burden, ttp, np.array([3, 6, 6, 4])


def test_patient_metrics_on_known_paths():
    s = _load()
    burden, ttp, courses = _paths()
    m = s.patient_metrics(burden, ttp, courses)
    assert m["called"].tolist() == [False, True, True, True]
    assert np.isnan(m["ratio_call"][0]) and m["ratio_call"][1] == pytest.approx(1.3)
    assert m["ratio_call"][2] == pytest.approx(0.8) and m["ratio_call"][3] == pytest.approx(1.5)
    assert m["base"].tolist() == [730.0, 140.0, 730.0, 112.0]                          # first boundary at or above 1.2 x baseline
    assert m["c14"].tolist() == [730.0, 730.0, 730.0, 224.0]                           # first boundary at or above 1.4 x baseline
    assert m["days_above"].tolist() == [0.0, 22 * 28.0, 0.0, 23 * 28.0]                # boundaries 5..26 and 4..26
    assert m["log_ratio"][0] == pytest.approx(0.0) and m["log_ratio"][1] > 0 > m["log_ratio"][2]
    assert m["scan_dip"].tolist() == [False, False, True, False]
    assert m["rano"].tolist() == [730.0, 168.0, 168.0, 252.0]                          # a call beyond 730 d is capped
    summ = s.summarise_arm(m)
    assert summ["frac_called"] == pytest.approx(0.75) and summ["frac_called_below_B0"] == pytest.approx(1 / 3)
    assert summ["frac_called_below_1p2"] == pytest.approx(1 / 3) and summ["frac_scan_dip"] == pytest.approx(0.25)
    assert summ["median_ratio_call"] == pytest.approx(1.3)
    assert summ["rmst730_c14"] == pytest.approx((730 + 730 + 730 + 224) / 4)


def test_summary_without_any_call_has_undefined_ratios():
    s = _load()
    m = s.patient_metrics(np.ones((4, s.N_REC)), np.full(4, 1456.0), np.zeros(4))
    summ = s.summarise_arm(m)
    assert summ["frac_called"] == 0.0 and np.isnan(summ["median_ratio_call"]) and np.isnan(summ["frac_called_below_1p2"])
    assert summ["rmst730_rano"] == 730.0 and summ["rmst730_base"] == 730.0 and summ["rmst730_c14"] == 730.0


# ---------------------------------------------------------------- frontier and sweep helpers
def test_upper_envelope_keeps_only_the_concave_vertices_and_interpolates_inside_the_range():
    s = _load()
    hull = s.upper_envelope([0, 1, 2, 3, 2.5, 1.0], [0, 3, 4, 4.5, 3.0, 1.0])       # (2.5, 3) lies below the segment; (1, 1) below (1, 3)
    assert hull == [(0, 0), (1, 3), (2, 4), (3, 4.5)]
    assert s.frontier_at(hull, 1.5) == pytest.approx(3.5) and s.frontier_at(hull, 3.0) == pytest.approx(4.5)
    assert np.isnan(s.frontier_at(hull, 3.2)) and np.isnan(s.frontier_at(hull, -0.1))
    flat = s.upper_envelope([0, 1, 2], [0, 1, 4])                                    # (1, 1) is below the chord from (0, 0) to (2, 4)
    assert flat == [(0, 0), (2, 4)] and s.frontier_at(flat, 1) == pytest.approx(2.0)


def test_sweep_gaps_on_a_hand_made_case():
    s = _load()
    courses = np.array([[0, 0, 0, 0], [2, 2, 2, 2], [4, 4, 4, 4]], float)
    ttp = np.array([[100, 100, 100, 100], [300, 100, 100, 100], [300, 200, 100, 100]], float)
    res, hull = s.sweep_gaps(ttp, courses)
    assert hull == [(0, 100.0), (2, 150.0), (4, 175.0)]
    g = res["informed"]
    assert g["courses"] == pytest.approx(1.5) and g["rmst"] == pytest.approx(175.0)      # patient 0 ties go to fewer courses
    assert g["frontier"] == pytest.approx(137.5) and g["gap"] == pytest.approx(37.5)
    only_ones = s.sweep_gaps(ttp[:, :1], courses[:, :1] + 1.0)[0]["informed"]             # fixed courses 1..5; informed uses 3 (inside range)
    assert np.isfinite(only_ones["gap"])
    # two fixed schedules that each use 2 courses on average; the informed arm gives every patient the better one with no course
    ttp2 = np.array([[100.0, 300.0], [300.0, 100.0]])
    crs2 = np.array([[4.0, 0.0], [0.0, 4.0]])
    res2, _ = s.sweep_gaps(ttp2, crs2)
    assert res2["informed"]["courses"] == pytest.approx(0.0) and res2["informed"]["rmst"] == pytest.approx(300.0)
    assert np.isnan(res2["informed"]["gap"])                                              # fewer courses than any fixed schedule: undefined


# ---------------------------------------------------------------- arms and worlds
def test_informed_is_at_least_as_good_as_every_fixed_arm_and_the_costly_version_uses_fewer_courses():
    s = _load()
    sched, names = rg.candidate_schedules()
    sel = cal.draw_patients(WEAK, 30, [11, 1])
    ev = cal.draw_patients(WEAK, 40, [12, 1])
    arms, best = s.pick_arms(s.evaluate_schedules_burden(sel, sched), s.evaluate_schedules_burden(ev, sched), names)
    assert best in names and set(arms) == set(s.ARMS)
    mean = {a: float(np.minimum(v["ttp"], 730.0).mean()) for a, v in arms.items()}
    for a in ("front6", "front3", "front12", "best_fixed"):
        assert mean["informed"] >= mean[a] - 1e-9
    assert float(arms["informed_lam10"]["courses"].mean()) <= float(arms["informed"]["courses"].mean()) + 1e-9
    assert float(arms["front3"]["courses"].max()) <= 3 and float(arms["front12"]["courses"].max()) <= 12


def test_a_tiny_world_gives_finite_rows_for_both_families():
    s = _load()
    for w in (WEAK, MIX):
        out = s.evaluate_world(w, 24, 24, (11, 12))
        assert [r["arm"] for r in out["arms"]] == list(s.ARMS)
        for r in out["arms"]:
            for m in ("mean_courses", "rmst730_rano", "rmst730_c14", "rmst730_base", "mean_days_above", "mean_log_ratio", "frac_called", "frac_scan_dip", "median_ratio_728"):
                assert np.isfinite(r[m]), (r["arm"], m)
            assert 0.0 <= r["frac_called"] <= 1.0 and r["rmst730_rano"] <= 730.0 and r["rmst730_base"] <= 730.0
        assert len(out["sweep"]) == 2 * len(s.SCAN_SWEEP)
        assert {r["scan_days"] for r in out["sweep"]} == {28 * x for x in s.SCAN_SWEEP}
        for r in out["sweep"]:
            assert np.isfinite(r["rmst730"]) and r["rmst730"] <= 730.0 and np.isfinite(r["courses"])
        # at the primary interval the sweep and the arm table use the same simulation, so the informed restricted time agrees
        inf_arm = next(r for r in out["arms"] if r["arm"] == "informed")["rmst730_rano"]
        inf_sw = next(r for r in out["sweep"] if r["arm"] == "informed" and r["scan_days"] == 84)["rmst730"]
        assert inf_arm == pytest.approx(inf_sw)


# ---------------------------------------------------------------- summary over worlds and the labels
def _rows(n_worlds, front3, informed):
    s = _load()
    base = {m: 0.0 for m in s.METRICS}
    rows = []
    for w in range(n_worlds):
        for arm in s.ARMS:
            r = dict(base, world=w, arm=arm, note="")
            if arm == "front3":
                r.update(front3)
            if arm == "informed":
                r.update(informed)
            rows.append(r)
    return rows


def test_verdict_follows_the_rule_fixed_in_the_docstring():
    s = _load()
    f3 = {"rmst730_rano": 300.0, "rmst730_c14": 400.0, "rmst730_base": 400.0, "frac_called_below_1p2": 0.6, "mean_log_ratio": 0.5}

    def label(**kw):
        inf = {"rmst730_rano": 320.0, "rmst730_c14": 401.0, "rmst730_base": 401.0, "frac_called_below_1p2": 0.6, "mean_log_ratio": 0.5}
        inf.update(kw)
        return s.analyse(_rows(5, f3, inf))["contrasts"]["informed_minus_front3"]["verdict"]

    assert label() == "mostly the scan grid"
    assert label(frac_called_below_1p2=0.3) == "mostly the nadir reference"
    assert label(rmst730_c14=418.0, rmst730_base=418.0, mean_log_ratio=0.3) == "mostly disease control"
    assert label(rmst730_c14=410.0, rmst730_base=410.0) == "mixed"
    assert label(rmst730_rano=299.0) == "no gain"
    assert s.analyse(_rows(5, f3, f3))["contrasts"]["informed_minus_best_fixed"]["verdict"] is None


def test_sweep_verdict_compares_the_gap_at_28_days_with_the_gap_at_84_days():
    s = _load()

    def rows(g28, g84):
        out = []
        for arm in ("informed", "informed_lam10"):
            for days, g in ((28, g28), (56, 20.0), (84, g84), (112, 30.0)):
                for w in range(3):
                    out.append({"world": w, "scan_days": days, "arm": arm, "courses": 2.5, "rmst730": 300.0, "frontier": 300.0 - g, "gap": g,
                                "front12_rmst730": 290.0, "front6_rmst730": 280.0})
        return out

    dep = s.analyse_sweep(rows(5.0, 30.0))
    assert dep["verdict_informed"] == "scan-grid dependent" and dep["informed"]["28"]["n_above"] == 3
    assert dep["informed"]["84"]["gap_q25_q50_q75"] == [30.0, 30.0, 30.0]
    ind = s.analyse_sweep(rows(28.0, 30.0))
    assert ind["verdict_informed"] == "not scan-grid dependent"


def test_seed_blocks_are_distinct_and_inside_the_development_block():
    s = _load()
    h = _load("203_trackc_value_of_information.py", "stage2_203_for_204")
    used = []
    for cfg in cal.CONFIGS.values():
        used += [v for k, v in cfg.items() if k.endswith("_seed") and isinstance(v, int)]
    used += [x for pair in h.H_SEEDS.values() for x in pair]
    m = [x for pair in s.M_SEEDS.values() for x in pair]
    assert len(set(m)) == len(m) and not set(m) & set(used)
    assert all(20261800 <= x <= 20262099 for x in m)
    assert set(s.M_SEEDS) <= set(cal.CONFIGS)
