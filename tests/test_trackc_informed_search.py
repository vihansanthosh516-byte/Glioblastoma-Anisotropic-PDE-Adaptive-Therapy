"""Script 205 (how far below the best open-loop schedule the informed arm and the fixed frontier lie). Synthetic arrays and tiny
simulated worlds; no result file is read."""
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


def _load(file="205_trackc_informed_search.py", name="stage2_205"):
    spec = importlib.util.spec_from_file_location(name, REPO / "src" / file)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


WEAK = {"world": 0, "log_r_shift": 0.0, "log_m_scale": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
        "c": 0.05, "mu": 0.0, "nu": 0.0}
MIX = {"world": 1, "log_r_shift": 0.0, "log_m_scale_non": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
       "c": 0.05, "q": 0.4, "mu": 0.0, "nu": 0.0}


def _mask(*cycles, width=52):
    m = np.zeros(width, bool)
    m[list(cycles)] = True
    return m


# ---------------------------------------------------------------- schedules and moves
def test_truncate_cuts_late_courses_and_removes_duplicates():
    s = _load()
    sched, _ = rg.candidate_schedules()
    t = s.truncate(sched)
    assert not t[:, s.NC:].any()
    assert len(t) <= len(sched) and len({row.tobytes() for row in t}) == len(t)
    assert (sched[:, s.NC:].any(axis=1)).any()                    # some schedules of the family do have courses after cycle 26
    cut = sched.copy()
    cut[:, s.NC:] = False
    assert {row.tobytes() for row in t} == {row.tobytes() for row in cut}


def test_neighbor_masks_toggle_and_swap_exactly_one_move():
    s = _load()
    base = np.stack([_mask(0, 3, 7), _mask(), _mask(25, 26)])
    nb = s.neighbor_masks(base)
    assert nb.shape == (2 * s.NC - 1, 3, 52)
    for c in range(s.NC):                                          # toggles
        expect = base.copy()
        expect[:, c] = ~base[:, c]
        assert (nb[c] == expect).all()
    for c in range(s.NC - 1):                                      # swaps of neighbouring cycles
        expect = base.copy()
        expect[:, c], expect[:, c + 1] = base[:, c + 1], base[:, c]
        assert (nb[s.NC + c] == expect).all()
    assert (nb[:, :, s.NC:] == base[None, :, s.NC:]).all()         # nothing changes after the first NC cycles
    changed = (nb != base[None]).sum(axis=2)
    assert (changed[:s.NC] == 1).all()                              # a toggle changes one cell
    assert (changed[s.NC:] <= 2).all() and (changed[s.NC:] % 2 == 0).all()   # a swap changes two cells, or none when both are equal


def test_neighbor_masks_respect_the_course_budget():
    s = _load()
    full = _mask(*range(12))
    nb = s.neighbor_masks(full[None])
    assert (nb.sum(axis=2) <= 12).all()
    assert (nb[12:s.NC, 0] == full).all()                          # toggling an empty cycle on would give 13 courses: unchanged
    assert not (nb[0, 0] == full).all()                            # toggling a course off is allowed
    small = s.neighbor_masks(_mask(0)[None], budget=1)
    assert (small.sum(axis=2) <= 1).all()


# ---------------------------------------------------------------- the searches on synthetic values
def test_climb_patients_reaches_a_target_and_never_decreases():
    s = _load()
    targets = np.stack([_mask(2, 5, 9), _mask(0), _mask(3, 4, 5, 6)])
    start = np.zeros_like(targets)

    def value_fn(nb, ids):
        return -(nb != targets[ids][None]).sum(axis=2).astype(float)

    v0 = -(start != targets).sum(axis=1).astype(float)
    cur, val, used, conv = s.climb_patients(start, v0, value_fn, rounds=10)
    assert (cur == targets).all() and (val == 0).all() and conv.all()
    assert (val >= v0).all()
    assert (used <= (start != targets).sum(axis=1)).all()          # a swap can fix two cells at once, a toggle one
    cur1, val1, used1, conv1 = s.climb_patients(start, v0, value_fn, rounds=1)
    assert (val1 >= v0).all() and (used1 <= 1).all() and not conv1[2]   # the third patient needs more than one move


def test_climb_patients_leaves_a_patient_whose_start_is_already_best():
    s = _load()
    targets = np.stack([_mask(1, 2), _mask(4)])
    start = targets.copy()

    def value_fn(nb, ids):
        return -(nb != targets[ids][None]).sum(axis=2).astype(float)

    cur, val, used, conv = s.climb_patients(start, np.zeros(2), value_fn, rounds=5)
    assert (cur == start).all() and (used == 0).all() and conv.all()


def test_climb_fixed_follows_an_increasing_path_to_a_local_optimum():
    s = _load()
    target = _mask(1, 4, 8)

    def objective(masks):
        return -(masks != target[None]).sum(axis=1).astype(float)

    path = s.climb_fixed(np.zeros(52, bool), objective, rounds=10)
    vals = [v for _, v in path]
    assert all(b > a for a, b in zip(vals, vals[1:])) and (path[-1][0] == target).all() and vals[-1] == 0.0
    short = s.climb_fixed(np.zeros(52, bool), objective, rounds=1)
    assert len(short) == 2


# ---------------------------------------------------------------- envelope and labels
def test_upper_envelope_and_frontier_at():
    s = _load()
    hull = s.upper_envelope([0, 1, 2, 3, 1.5], [0, 10, 12, 13, 9])
    assert [round(h[0], 6) for h in hull] == [0.0, 1.0, 2.0, 3.0]
    assert s.frontier_at(hull, 1.5) == pytest.approx(11.0)
    assert np.isnan(s.frontier_at(hull, 3.5)) and np.isnan(s.frontier_at(hull, -0.1))


def test_labels_follow_the_rule_fixed_in_the_docstring():
    s = _load()
    assert s.label_g1(1.99) == "best of 56 is near the optimum"
    assert s.label_g1(2.0) == "intermediate" and s.label_g1(4.99) == "intermediate"
    assert s.label_g1(5.0) == "best of 56 understates the informed optimum"
    assert s.label_room(9.99) == "room above fixed dosing stays under 10 d"
    assert s.label_room(10.0) == "room above fixed dosing is 10 d or more"


# ---------------------------------------------------------------- simulation
def test_per_patient_evaluation_equals_evaluate_schedules():
    s = _load()
    sched, names = rg.candidate_schedules()
    t = s.truncate(sched)
    scen = cal.draw_patients(WEAK, 10, [21, 1])
    pick = [0, 3, 7, len(t) - 1]
    for endpoint in ("rano", "baseline"):
        ttp_ref, c_ref = rg.evaluate_schedules(scen, t[pick], endpoint)
        nb = np.repeat(t[pick][:, None, :], 10, axis=1)            # (moves, patients, width): the same schedule for every patient
        ttp, c = s.eval_per_patient(scen, nb, endpoint)
        assert np.allclose(ttp, ttp_ref) and (c == c_ref).all()


def test_a_tiny_world_gives_finite_rows_and_the_search_never_loses_to_the_family():
    s = _load()
    out = s.evaluate_world((MIX, 6, 8, (31, 32), 1, 1))             # one mixture world keeps the run short
    assert len(out["search"]) == len(out["front"]) == len(s.ENDPOINTS) * len(s.LAMS)
    for r in out["search"]:
        assert all(np.isfinite(v) for k, v in r.items() if k != "endpoint")
        assert r["mean_gain"] >= -1e-9 and 0.0 <= r["frac_improved"] <= 1.0 and 0.0 <= r["frac_not_converged"] <= 1.0
        assert r["mean_value_search"] >= r["mean_value_family"] - 1e-9
    for r in out["front"]:
        assert r["gap_family_vs_family_frontier"] >= -1e-6                      # the family maximum lies on or above its own frontier
        assert r["frontier_gain_at_family_courses"] >= -1e-9 and r["frontier_gain_at_search_courses"] >= -1e-9
        assert r["gap_search_vs_searched_frontier"] <= r["gap_search_vs_family_frontier"] + 1e-9
        assert r["gap_family_vs_searched_frontier"] <= r["gap_family_vs_family_frontier"] + 1e-9


def test_analysis_summarises_rows_and_attaches_labels():
    s = _load()
    search = [{"world": w, "endpoint": ep, "lam": lam, "mean_gain": g, "frac_improved": 0.5, "frac_not_converged": 0.1,
               "courses_family": 2.0, "courses_search": 2.5}
              for w, g in ((0, 1.0), (1, 3.0), (2, 9.0)) for ep in s.ENDPOINTS for lam in s.LAMS]
    keys = ("gap_family_vs_family_frontier", "gap_search_vs_family_frontier", "gap_family_vs_searched_frontier",
            "gap_search_vs_searched_frontier", "frontier_gain_at_family_courses", "frontier_gain_at_search_courses")
    front = [dict({k: float(i) for k in keys}, world=w, endpoint=ep, lam=lam)
             for i, w in enumerate((0, 1, 2)) for ep in s.ENDPOINTS for lam in s.LAMS]
    res = s.analyse(search, front)
    assert res["n_worlds"] == 3
    b = res["by_endpoint"]["rano"]
    assert b["lam0"]["mean_gain_q25_q50_q75"][1] == pytest.approx(3.0)
    assert b["label_g1"] == "intermediate" and b["label_room"] == "room above fixed dosing stays under 10 d"
    assert b["lam0"]["gap_search_vs_searched_frontier"]["n_above"] == 2        # the values 0, 1, 2: two above 0.05


def test_seed_blocks_are_distinct_and_inside_the_development_block():
    s = _load()
    h = _load("203_trackc_value_of_information.py", "stage2_203_for_205")
    m = _load("204_trackc_endpoint_mechanism.py", "stage2_204_for_205")
    used = []
    for cfg in cal.CONFIGS.values():
        used += [v for k, v in cfg.items() if k.endswith("_seed") and isinstance(v, int)]
    used += [x for pair in h.H_SEEDS.values() for x in pair]
    used += [x for pair in m.M_SEEDS.values() for x in pair]
    g = [x for pair in s.G_SEEDS.values() for x in pair]
    assert len(set(g)) == len(g) and not set(g) & set(used)
    assert all(20261800 <= x <= 20262099 for x in g)
    assert set(s.G_SEEDS) <= set(cal.CONFIGS)
    assert all(pair[1] == pair[0] + 1 for pair in s.G_SEEDS.values())          # the registry takes (low, high) ranges
