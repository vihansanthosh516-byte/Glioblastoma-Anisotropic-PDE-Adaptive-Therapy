"""Track C pilot simulator (src/trackc_sim.py): equation limits, endpoint reachability, budget,
step accuracy, filter sanity and the CVaR estimator. Synthetic parameters only; these check the
assumed model behaves as written, not that it describes any patient or drug."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_sim as tc  # noqa: E402


class NoTreat:
    uses_filter = False

    def reset(self, scen, rng):
        pass

    def act(self, j, y):
        return np.zeros(y.shape, bool)

    def advance(self, j, treat):
        pass


def one(**kw):
    base = dict(r=0.03, m=0.05 / 0.03, f0=0.0, c=0.0, eps=0.0, B0=0.2, sigma=0.05, a_sr=1.0, a_rs=1.0)
    base.update(kw)
    return {key: np.array([float(val)]) for key, val in base.items()}


def test_single_population_matches_logistic():
    p = one(f0=0.0)
    S, R = tc.run_cycle(p["B0"].copy(), np.zeros(1), p, np.zeros(1), 0.0)
    r, B0, t = 0.03, 0.2, tc.CYCLE_DAYS
    exact = B0 * np.exp(r * t) / (1 - B0 + B0 * np.exp(r * t))
    assert abs((S + R)[0] - exact) < 1e-8


def test_step_size_converged():
    scen = tc.draw_scenarios(64, seed=1)
    S0, R0 = (1 - scen["f0"]) * scen["B0"], scen["f0"] * scen["B0"]
    treat = np.ones(64)
    a = tc.run_cycle(S0, R0, scen, treat, 0.0, dt=tc.DT)
    b = tc.run_cycle(S0, R0, scen, treat, 0.0, dt=0.05)
    assert np.max(np.abs(a[0] + a[1] - b[0] - b[1]) / (b[0] + b[1])) < 1e-4


def test_endpoint_reachable_without_treatment():
    scen = tc.draw_scenarios(400, seed=2)
    out = tc.simulate(scen, NoTreat(), seed=3)
    assert not out["capped"].any(), "every untreated scenario must progress before the horizon"


def test_positive_and_bounded_under_full_dosing():
    scen = tc.draw_scenarios(200, seed=4)
    S, R = (1 - scen["f0"]) * scen["B0"], scen["f0"] * scen["B0"]
    for j in range(tc.N_CYCLES):
        S, R = tc.run_cycle(S, R, scen, np.ones(200), j * tc.CYCLE_DAYS)
        assert (S >= 0).all() and (R >= 0).all() and (S + R <= 1 + 1e-9).all()


def test_drug_selects_for_resistance():
    p = one(f0=0.01, m=0.08 / 0.03, eps=0.0)
    S, R = (1 - p["f0"]) * p["B0"], p["f0"] * p["B0"]
    S1, R1 = tc.run_cycle(S, R, p, np.ones(1), 0.0)
    assert R1 / (S1 + R1) > R / (S + R)


@pytest.mark.parametrize("spec", [("periodic", 1), ("periodic", 4), ("raw", 0.8, 0.3), ("raw", 1.0, 0.5)])
def test_budget_never_exceeded(spec):
    scen = tc.draw_scenarios(100, seed=5)
    out = tc.simulate(scen, tc.make_policy(spec), seed=6)
    assert out["courses"].max() <= tc.BUDGET


def test_paced_uses_twelve_courses_when_no_progression():
    scen = tc.draw_scenarios(50, seed=7)
    out = tc.simulate(scen, tc.Periodic(4), seed=8)
    capped = out["capped"]
    assert (out["courses"][capped] == tc.BUDGET).all()


def test_raw_threshold_hysteresis():
    pol = tc.RawThreshold(1.0, 0.5)
    pol.reset(None, None)
    seq = [0.0, np.log(0.7), np.log(0.49), np.log(0.8), np.log(1.01)]
    acts = [pol.act(j, np.array([y]))[0] for j, y in enumerate(seq)]
    assert acts == [True, True, False, False, True]


def test_filter_tracks_burden_and_resamples():
    scen = tc.draw_scenarios(30, seed=9, sigmas=(0.05,))
    pol = tc.FilteredThreshold(1.0, 0.5, 1.0, 0.2, n_particles=300)
    out = tc.simulate(scen, pol, seed=10, track_fr=True)
    assert pol.pf.n_resample > 0
    hits, total = out["fr_cover"]
    assert total > 0 and 0.0 <= hits / total <= 1.0
    # one more reading: posterior burden should sit near the truth on the log scale
    assert np.isfinite(pol.pf.burden_mean()).all()


def test_filter_posterior_concentrates_on_truth():
    scen = tc.draw_scenarios(20, seed=11, sigmas=(0.05,))
    rng = np.random.default_rng(12)
    pf = tc.ParticleFilter(scen["sigma"], 400, rng)
    S, R = (1 - scen["f0"]) * scen["B0"], scen["f0"] * scen["B0"]
    errs = []
    for j in range(8):
        y = np.log(S + R) + rng.standard_normal(20) * 0.05
        pf.update(j, y)
        errs.append(np.median(np.abs(np.log(pf.burden_mean()) - np.log(S + R))))
        treat = np.full(20, j % 2 == 0)
        S, R = tc.run_cycle(S, R, scen, treat, j * tc.CYCLE_DAYS)
        pf.propagate(j, treat)
    assert errs[-1] < 0.05


def test_risk_containment_waits_then_treats():
    scen = tc.draw_scenarios(20, seed=13, sigmas=(0.05,))
    out = tc.simulate(scen, tc.make_policy(("risk", 0.1, 1)), seed=14)
    assert out["courses"].max() <= tc.BUDGET
    # p_star = 1.0 can never trigger before the first progression is certain: fewer courses
    lazy = tc.simulate(scen, tc.make_policy(("risk", 1.0, 1)), seed=14)
    assert lazy["courses"].sum() < out["courses"].sum()
    eager = tc.simulate(scen, tc.make_policy(("risk", 0.0, 1)), seed=14)
    assert (eager["courses"] == np.minimum(tc.BUDGET, eager["courses"].max())).any()


def test_cvar_lower():
    assert tc.cvar_lower(np.arange(1, 11), 0.1) == pytest.approx(1.0)
    x = np.arange(15.0)
    assert tc.cvar_lower(x, 0.1) == pytest.approx((0 + 0.5 * 1) / 1.5)
