"""Stage 2 simulator options (src/trackc_sim.py, added 2026-11-04): regular scan interval ("scan_every")
and the RANO-like endpoint (simulate(..., endpoint="rano")). Synthetic parameters only.
Both options are off by default; the golden-value test pins the default outputs recorded from the code
before the options were added (record: scratch run on 2026-11-04, draw_scenarios(8, 3), noise key [7, 0])."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_sim as tc  # noqa: E402

GOLDEN = {
    "periodic1": ([1363.480165, 956.295921, 518.992476, 814.849302, 214.440321, 928.527043, 499.492699, 398.920121],
                  [12, 12, 12, 12, 8, 12, 12, 12]),
    "periodic3": ([1431.187373, 1067.119889, 151.730168, 665.843746, 74.710055, 755.748058, 67.260896, 162.087427],
                  [12, 12, 2, 8, 1, 9, 1, 2]),
    "raw0.6": ([1407.392387, 1034.471517, 561.150753, 814.849302, 162.961173, 943.147738, 507.579893, 440.256302],
               [12, 12, 12, 12, 3, 12, 12, 9]),
}


def _policy(name):
    return {"periodic1": tc.Periodic(1), "periodic3": tc.Periodic(3), "raw0.6": tc.RawThreshold(0.6, 0.6)}[name]


@pytest.mark.parametrize("name", list(GOLDEN))
def test_defaults_reproduce_recorded_outputs(name):
    out = tc.simulate(tc.draw_scenarios(8, 3), _policy(name), [7, 0])
    ttp, courses = GOLDEN[name]
    np.testing.assert_allclose(out["ttp"], ttp, atol=2e-6)
    assert out["courses"].tolist() == courses


class Spy:
    """Wraps a policy, records what it was shown and what it decided (before the course budget)."""
    def __init__(self, inner):
        self.inner, self.uses_filter = inner, getattr(inner, "uses_filter", False)
        self.seen, self.acts = [], []

    def reset(self, scen, rng):
        self.seen, self.acts = [], []
        self.inner.reset(scen, rng)

    def act(self, j, y):
        self.seen.append(y.copy())
        a = self.inner.act(j, y)
        self.acts.append(np.array(a, bool))
        return a

    def advance(self, j, treat):
        self.inner.advance(j, treat)


def test_scan_every_hides_readings_between_scans():
    scen = tc.draw_scenarios(6, 4)
    scen["scan_every"] = np.array([1, 2, 3, 3, 4, 6])
    spy = Spy(tc.RawThreshold(0.6, 0.6))
    tc.simulate(scen, spy, [1, 0])
    for j, y in enumerate(spy.seen):
        should_see = (j % scen["scan_every"]) == 0
        assert (~np.isnan(y) == should_see).all(), j
    assert not np.isnan(spy.seen[0]).any()               # the baseline scan always exists


def test_without_scan_every_nothing_is_hidden():
    scen = tc.draw_scenarios(4, 4)
    spy = Spy(tc.RawThreshold(0.6, 0.6))
    tc.simulate(scen, spy, [1, 0])
    assert not any(np.isnan(y).any() for y in spy.seen)


def test_raw_rule_holds_its_decision_between_scans():
    scen = tc.draw_scenarios(6, 5)
    scen["scan_every"] = np.full(6, 3)
    spy = Spy(tc.RawThreshold(0.6, 0.6))
    tc.simulate(scen, spy, [2, 0])
    for j in range(1, len(spy.acts)):
        if j % 3 != 0:
            assert (spy.acts[j] == spy.acts[j - 1]).all(), j


class NoTreat:
    uses_filter = False

    def reset(self, scen, rng):
        pass

    def act(self, j, y):
        return np.zeros(y.shape, bool)

    def advance(self, j, treat):
        pass


def _untreated_world():
    scen = tc.draw_scenarios(4, 6)
    scen["r"] = np.array([0.01, 0.02, 0.005, 0.03])
    scen["B0"] = np.array([0.10, 0.20, 0.30, 0.05])
    scen["c"] = np.zeros(4)
    scen["f0"] = np.full(4, 0.01)
    return scen


@pytest.mark.parametrize("scan_every", [1, 3])
def test_rano_endpoint_untreated_matches_closed_form(scan_every):
    scen = _untreated_world()
    scen["scan_every"] = np.full(4, scan_every)
    r, B0 = scen["r"], scen["B0"]
    t_cross = np.log(tc.RANO_FACTOR * (1 - B0) / (1 - tc.RANO_FACTOR * B0)) / r
    step = tc.CYCLE_DAYS * scan_every
    expected = np.ceil(t_cross / step) * step                 # first scan at or after the crossing
    out = tc.simulate(scen, NoTreat(), [3, 0], endpoint="rano")
    np.testing.assert_allclose(out["ttp"], np.minimum(expected, tc.HORIZON), atol=1e-9)


def test_baseline_endpoint_untreated_matches_closed_form_and_precedes_rano():
    scen = _untreated_world()
    r, B0 = scen["r"], scen["B0"]
    t_base = np.log(tc.PROG_FACTOR * (1 - B0) / (1 - tc.PROG_FACTOR * B0)) / r
    base = tc.simulate(scen, NoTreat(), [3, 0])
    np.testing.assert_allclose(base["ttp"], np.minimum(t_base, tc.HORIZON), atol=0.01)
    rano = tc.simulate(scen, NoTreat(), [3, 0], endpoint="rano")
    assert (rano["ttp"] > base["ttp"]).all()                  # 1.4 x nadir is later than 1.2 x baseline here


def _rano_reference(scen, policy_fn, scan_every):
    """Independent re-implementation of the RANO-like rule for a policy that ignores readings."""
    n = len(scen["r"])
    S, R = (1 - scen["f0"]) * scen["B0"], scen["f0"] * scen["B0"]
    base, nadir = S + R, S + R
    det, used = np.full(n, np.nan), np.zeros(n, int)
    for j in range(tc.N_CYCLES):
        scan = (j % scan_every) == 0
        B = S + R
        if j > 0:
            hit = scan & np.isnan(det) & (B >= tc.RANO_FACTOR * np.maximum(nadir, tc.RANO_FLOOR * base))
            det = np.where(hit, j * tc.CYCLE_DAYS, det)
            nadir = np.where(scan, np.minimum(nadir, B), nadir)
        treat = policy_fn(j, n) & (used < tc.BUDGET) & np.isnan(det)
        used = used + treat
        S, R = tc.run_cycle(S, R, scen, treat.astype(float), j * tc.CYCLE_DAYS)
        if not np.isnan(det).any():
            break
    det = np.where(np.isnan(det), tc.HORIZON, det)
    return det, used


@pytest.mark.parametrize("scan_every", [1, 3])
@pytest.mark.parametrize("q", [1, 3])
def test_rano_endpoint_with_treatment_matches_reference(scan_every, q):
    scen = tc.draw_scenarios(10, 8)
    scen["scan_every"] = np.full(10, scan_every)
    out = tc.simulate(scen, tc.Periodic(q), [4, 0], endpoint="rano")
    det, used = _rano_reference(scen, lambda j, n: np.full(n, j % q == 0), scan_every)
    np.testing.assert_allclose(out["ttp"], det, atol=1e-9)
    assert out["courses"].tolist() == used.tolist()


def test_budget_option_caps_courses_and_matches_oracle_rollout():
    import trackc_oracle as to
    scen = tc.draw_scenarios(6, 12)
    out = tc.simulate(scen, tc.Periodic(1), [5, 0], budget=6)
    assert out["courses"].max() <= 6 and out["courses"].max() >= 2
    sched = np.zeros(tc.N_CYCLES, bool)
    sched[:6] = True
    for i in range(6):
        p = to._scalar(scen, i)
        assert to.rollout(p, sched[None])[0] == pytest.approx(out["ttp"][i], abs=1e-6)
    assert tc.simulate(scen, tc.Periodic(1), [5, 0])["courses"].max() <= tc.BUDGET   # default unchanged


def test_dose_multiplier_scales_the_drug_effect_on_resistance():
    p = {key: np.array([val], float) for key, val in dict(r=0.01, m=10.0, c=0.0, eps=0.1, a_rs=1.0, a_sr=1.0).items()}
    shifts = []
    for dose in (1.0, 2.0):
        S, R = np.array([0.19]), np.array([0.01])
        start = np.log(R / S)
        for j in range(3):
            S, R = tc.run_cycle(S, R, p, np.array([dose]), j * tc.CYCLE_DAYS)
        shifts.append(float((np.log(R / S) - start)[0]))
    assert shifts[1] == pytest.approx(2.0 * shifts[0], rel=1e-6)     # Lemma 2: shift = k (1 - eps) * dose * drug-days


def test_dose_option_is_off_by_default_and_changes_outcomes_when_set():
    scen = tc.draw_scenarios(6, 13)
    plain = tc.simulate(scen, tc.Periodic(1), [5, 0])
    same = tc.simulate(dict(scen, dose=np.ones(6)), tc.Periodic(1), [5, 0])
    np.testing.assert_allclose(same["ttp"], plain["ttp"], atol=1e-12)
    double = tc.simulate(dict(scen, dose=np.full(6, 2.0)), tc.Periodic(1), [5, 0])
    assert not np.allclose(double["ttp"], plain["ttp"])


def test_rano_floor_delays_progression_after_a_deep_response(monkeypatch):
    """After a deep response the nadir is tiny. Without a floor, a 40% regrowth of a nearly invisible tumour
    counts as progression at the next scan. The floor (10% of baseline) stands in for the minimum measurable
    lesion and delays the call until the tumour is measurable again."""
    scen = tc.draw_scenarios(1, 9)
    scen.update(r=np.array([0.03]), m=np.array([20.0]), f0=np.array([1e-4]), c=np.zeros(1), eps=np.zeros(1),
                B0=np.array([0.2]), sigma=np.array([0.05]), scan_every=np.array([1]))
    monkeypatch.setattr(tc, "RANO_FLOOR", 0.0)
    no_floor = tc.simulate(scen, tc.Periodic(1), [5, 0], endpoint="rano")["ttp"][0]
    monkeypatch.setattr(tc, "RANO_FLOOR", 0.1)
    with_floor = tc.simulate(scen, tc.Periodic(1), [5, 0], endpoint="rano")["ttp"][0]
    assert no_floor < 400                                      # sensible: the resistant clone regrows within a year
    assert with_floor >= no_floor + 2 * tc.CYCLE_DAYS          # the floor delays the call by at least two scans
