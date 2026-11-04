"""Script 203 (value of knowing a patient-level quantity). Synthetic score matrices and tiny simulated worlds; no result file is read."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_calibration as cal  # noqa: E402


def _load():
    spec = importlib.util.spec_from_file_location("stage2_203", REPO / "src" / "203_trackc_value_of_information.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- helpers on synthetic scores
def test_quantile_labels_make_equal_groups_and_share_the_cut_points():
    s = _load()
    x = np.arange(100.0)
    lab = s.quantile_labels(x, 4)
    assert [int((lab == g).sum()) for g in range(4)] == [25, 25, 25, 25]
    # evaluation values use the cut points of the pooled reference, so a shifted set lands in the upper groups
    assert s.quantile_labels(np.array([-5.0, 10.0, 60.0, 500.0]), 4, reference=x).tolist() == [0, 0, 2, 3]
    assert set(s.quantile_labels(np.zeros(10), 4).tolist()) <= {0, 1, 2, 3}          # ties never fail


def _two_types(n=40):
    """Type A patients (even index) gain from schedule 0, type B (odd index) from schedule 1."""
    typ = np.arange(n) % 2
    score = np.zeros((2, n))
    score[0, typ == 0], score[1, typ == 1] = 10.0, 10.0
    return score, typ


def test_group_value_is_perfect_when_the_label_separates_the_types_and_equal_to_fixed_when_it_does_not():
    s = _load()
    score, typ = _two_types()
    fixed = s.group_value(score, score, np.zeros(40, int), np.zeros(40, int))
    perfect = s.perfect_value(score)
    informative = s.group_value(score, score, typ, typ)
    useless = s.group_value(score, score, (np.arange(40) // 2) % 2, (np.arange(40) // 2) % 2)   # both groups mix the types equally
    assert fixed == pytest.approx(5.0) and perfect == pytest.approx(10.0)
    assert informative == pytest.approx(10.0) and useless == pytest.approx(5.0)
    assert s.share(informative, fixed, perfect) == pytest.approx(1.0) and s.share(useless, fixed, perfect) == pytest.approx(0.0)


def test_group_value_uses_selection_patients_for_the_choice_and_evaluation_patients_for_the_score():
    s = _load()
    sel = np.array([[5.0, 5.0], [0.0, 0.0]])             # on the selection set schedule 0 is best
    ev = np.array([[1.0, 1.0], [9.0, 9.0]])              # on the evaluation set schedule 1 would be better
    assert s.group_value(sel, ev, np.zeros(2, int), np.zeros(2, int)) == pytest.approx(1.0)   # the choice is not tuned on ev
    assert s.perfect_value(ev) == pytest.approx(9.0)


def test_group_without_selection_patients_falls_back_to_the_overall_best():
    s = _load()
    sel = np.array([[3.0, 3.0, 3.0], [1.0, 1.0, 1.0]])
    ev = np.array([[3.0, 3.0], [1.0, 1.0]])
    assert s.group_value(sel, ev, np.array([0, 0, 0]), np.array([0, 1])) == pytest.approx(3.0)


def test_share_is_undefined_when_nothing_can_be_gained():
    s = _load()
    assert np.isnan(s.share(1.0, 1.0, 1.0)) and np.isnan(s.share(1.0, 1.0, 1.0 + 1e-9))


def test_score_matrix_caps_the_time_and_charges_the_courses():
    s = _load()
    ttp = np.array([[100.0, 900.0], [50.0, 800.0]])
    courses = np.array([[2.0, 4.0], [0.0, 1.0]])
    out = s.score_matrix(ttp, courses, 10.0)
    assert out[0, 0] == pytest.approx(100.0 - 20.0 - 2e-6) and out[0, 1] == pytest.approx(730.0 - 40.0 - 4e-6)
    assert out[1, 0] == pytest.approx(50.0) and out[1, 1] == pytest.approx(730.0 - 10.0 - 1e-6)
    assert out[0, 0] < s.score_matrix(ttp, np.zeros_like(courses), 10.0)[0, 0]


def test_analyse_summarises_shares_over_worlds_and_counts_the_best_single_quantity():
    s = _load()
    rows = []
    for w in range(4):
        for lam in (0.0,):
            base = {"world": w, "lam": lam, "v_fixed": 10.0, "v_perfect": 30.0}
            rows.append(dict(base, quantity="a", v_group=25.0, gain_days=15.0, share=0.75))
            rows.append(dict(base, quantity="b", v_group=15.0, gain_days=5.0, share=0.25))
            rows.append(dict(base, quantity="a+b", v_group=28.0, gain_days=18.0, share=0.9))
    res = s.analyse(rows)
    block = res["by_lambda"]["0"]
    assert res["n_worlds"] == 4 and block["perfect_minus_fixed_days"] == [20.0, 20.0, 20.0]
    assert block["quantities"]["a"]["share_q25_q50_q75"] == [0.75, 0.75, 0.75] and block["quantities"]["a"]["n_worlds_share_ge_half"] == 4
    assert block["quantities"]["b"]["n_worlds_share_ge_half"] == 0
    assert block["best_single_quantity_counts"] == {"a": 4}                      # the pair 'a+b' is not a single quantity


# ---------------------------------------------------------------- groups and tiny simulated worlds
WEAK = {"world": 0, "log_r_shift": 0.0, "log_m_scale": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
        "c": 0.05, "mu": 0.0, "nu": 0.0}
MIX = {"world": 1, "log_r_shift": 0.0, "log_m_scale_non": float(np.log(0.05)), "log_f0_shift": 0.0, "B0_center": 0.3,
       "c": 0.05, "q": 0.4, "mu": 0.0, "nu": 0.0}


def test_quantity_groups_cover_singles_pairs_and_the_responder_class_only_for_the_mixture_family():
    s = _load()
    for w, with_class in ((WEAK, False), (MIX, True)):
        sel, ev = cal.draw_patients(w, 80, [5, 1]), cal.draw_patients(w, 80, [6, 1])
        groups = s.quantity_groups(sel, ev, w)
        assert ("class" in groups) is with_class
        assert len([g for g in groups if "+" not in g and g != "class"]) == 5 and len([g for g in groups if "+" in g]) == 10
        for name, (a, b) in groups.items():
            assert a.shape == (80,) and b.shape == (80,)
            top = 1 if name == "class" else (s.K_SINGLE - 1 if "+" not in name else s.K_PAIR * s.K_PAIR - 1)
            assert a.min() >= 0 and a.max() <= top and b.min() >= 0 and b.max() <= top
    assert groups["class"][0].dtype.kind in "iu"


def test_a_tiny_world_gives_ordered_finite_values():
    s = _load()
    rows = s.evaluate_world(WEAK, 40, 40, (11, 12))
    assert len(rows) == len(s.LAMBDAS) * 15
    for r in rows:
        assert np.isfinite(r["v_fixed"]) and np.isfinite(r["v_perfect"]) and np.isfinite(r["v_group"])
        assert r["v_perfect"] >= r["v_group"] - 1e-9 and r["v_perfect"] >= r["v_fixed"] - 1e-9     # a per-patient maximum is an upper bound
        assert r["share"] != r["share"] or r["share"] <= 1.0 + 1e-9                                  # NaN or at most 1
    assert len(s.evaluate_world(MIX, 40, 40, (11, 12))) == len(s.LAMBDAS) * 16


def test_seed_blocks_are_distinct_and_inside_the_development_block():
    s = _load()
    used = []
    for cfg in cal.CONFIGS.values():
        used += [v for k, v in cfg.items() if k.endswith("_seed") and isinstance(v, int)]
    h = [x for pair in s.H_SEEDS.values() for x in pair]
    assert len(set(h)) == len(h) and not set(h) & set(used)
    assert all(20261800 <= x <= 20262099 for x in h)
    assert set(s.H_SEEDS) == set(cal.CONFIGS)
