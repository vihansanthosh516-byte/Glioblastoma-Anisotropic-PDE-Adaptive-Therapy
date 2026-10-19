"""Tests for src/rano_calc.py, fixed in Amendment 17 (A17.3). Run once after the code is committed."""
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("rano_calc", ROOT / "src" / "rano_calc.py")
rc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rc)

PAIRS = [(0.8, 0.1), (1.0, 0.1), (1.4, 0.1), (1.4, 0.3), (1.2, 0.2), (1.8, 0.25),
         (0.35, 0.2), (0.5, 0.5), (0.2, 0.3), (2.5, 0.7), (1.0, 1.0), (1.4, 1.2)]


def test1_closed_form_vs_monte_carlo():
    rng = np.random.default_rng(20261018)
    n = 200000
    for ratio, s in PAIRS:
        m = math.log(ratio)
        x = m + s * rng.standard_normal(n)
        mc = {"P_PD": float((x >= rc.LN_PD).mean()), "P_PR": float((x <= rc.LN_PR).mean())}
        mc["P_SD"] = 1.0 - mc["P_PD"] - mc["P_PR"]
        fw = rc.forward(ratio, s)
        bw = rc.backward(1000.0, 1000.0 * ratio, s)  # posterior r | obs ~ N(m, s^2): same normal
        for k in mc:
            assert abs(fw[k] - mc[k]) < 0.005, (ratio, s, k)
            assert abs(bw[k] - mc[k]) < 0.005, (ratio, s, k)


def test2_probabilities_sum_to_one():
    rng = np.random.default_rng(1)
    for _ in range(50):
        v0 = float(rng.uniform(500, 50000)); v1 = float(rng.uniform(100, 80000)); s = float(rng.uniform(0.01, 1.5))
        for d in (rc.backward(v0, v1, s), rc.forward(v1 / v0, s)):
            assert abs(sum(d.values()) - 1.0) < 1e-9
            assert all(-1e-12 <= p <= 1 + 1e-12 for p in d.values())


def test3_monotone_and_half_at_line():
    s = 0.3
    ps = [rc.backward(1000.0, 1000.0 * r, s)["P_PD"] for r in np.linspace(0.5, 3.0, 40)]
    assert all(b >= a for a, b in zip(ps, ps[1:]))
    assert abs(rc.backward(1000.0, 1400.0, s)["P_PD"] - 0.5) < 1e-9


def test4_symmetry_forward():
    assert abs(rc.forward(1.40, 0.4)["P_PD"] - 0.5) < 1e-9


def test5_s_zero_matches_rule():
    rng = np.random.default_rng(2)
    for _ in range(1000):
        v0 = float(rng.uniform(500, 50000)); v1 = float(rng.uniform(100, 100000))
        r = v1 / v0
        want = "PD" if r >= 1.40 else ("PR" if r <= 0.35 else "SD")
        got = rc.backward(v0, v1, 0.0)
        assert got["P_" + want] == 1.0 and rc.classify(v0, v1) == want


def test6_floor_flag():
    assert rc.rano_calc(499.9, 2000.0, 0.2)["below_measurability_floor"] is True
    assert rc.rano_calc(2000.0, 499.9, 0.2)["below_measurability_floor"] is True
    assert rc.rano_calc(500.0, 500.0, 0.2)["below_measurability_floor"] is False


def test7_bracket_equals_a15_file():
    t2 = json.loads((ROOT / "output" / "rano_noise" / "rano_tool_dependence.json").read_text())["TEST2_mdc_vs_progression_line"]
    for b in ("0", "1", "2"):
        lo, hi = rc.bracket(int(b))
        want = sorted((t2["bands"][b]["s_change_upper_bound"], t2["bands"][b]["s_change_lower_bound"]))
        assert (lo, hi) == (want[0], want[1])
