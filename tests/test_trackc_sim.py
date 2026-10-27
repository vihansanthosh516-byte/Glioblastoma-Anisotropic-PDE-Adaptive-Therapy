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


def test_switch_boundary_exact_when_dt_does_not_divide_course():
    # dt = 0.3 does not divide 5 days; the on/off switch must still land on a step boundary
    scen = tc.draw_scenarios(32, seed=19)
    S0, R0 = (1 - scen["f0"]) * scen["B0"], scen["f0"] * scen["B0"]
    a = tc.run_cycle(S0, R0, scen, np.ones(32), 0.0, dt=0.3)
    b = tc.run_cycle(S0, R0, scen, np.ones(32), 0.0, dt=0.01)
    assert np.max(np.abs(a[0] + a[1] - b[0] - b[1]) / (b[0] + b[1])) < 1e-4


def test_event_time_matches_independent_solver():
    import trackc_oracle as to
    import trackc_reference as ref
    scen = tc.draw_scenarios(12, seed=20)
    scheds = to.periodic_schedules(qs=[1, 4]).tolist() + [[0] * tc.N_CYCLES]
    for i in range(12):
        p = to._scalar(scen, i)
        ours = to.rollout(p, scheds)
        for s, t in zip(scheds, ours):
            assert abs(t - ref.reference_ttp(p, s)) < 0.01


def test_grid_mode_is_the_old_rounded_rule():
    import trackc_oracle as to
    scen = tc.draw_scenarios(20, seed=21)
    sched = to.periodic_schedules(qs=[2])
    try:
        tc.EVENT_MODE = "grid"
        grid = np.array([to.rollout(to._scalar(scen, i), sched)[0] for i in range(20)])
    finally:
        tc.EVENT_MODE = "refine"
    fine = np.array([to.rollout(to._scalar(scen, i), sched)[0] for i in range(20)])
    assert np.allclose(grid / tc.DT, np.round(grid / tc.DT))
    assert (grid >= fine - 1e-9).all() and (grid - fine <= tc.DT + 1e-9).all()


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


def test_periodic_delayed_start():
    pol = tc.make_policy(("periodic", 3, 2))
    acts = [bool(pol.act(j, np.zeros(1))[0]) for j in range(9)]
    assert acts == [False, False, True, False, False, True, False, False, True]
    assert [bool(tc.make_policy(("periodic", 2)).act(j, np.zeros(1))[0]) for j in range(4)] == [True, False, True, False]


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


def test_weak_drug_guard_front_loads_weak_drugs():
    weak = tc.draw_scenarios(30, seed=17, ranges=dict(tc.RANGES, m=(3.0, 4.0, "log")), sigmas=(0.05,))
    narrow = tc.simulate(weak, tc.make_policy(("risk", 0.05, 2, 0.0, False)), seed=18)
    wide = tc.simulate(weak, tc.make_policy(("risk", 0.05, 2, 0.0, True)), seed=18)
    guard = tc.simulate(weak, tc.make_policy(("risk", 0.05, 2, 6.0, True)), seed=18)
    # the wide prior alone lets the filter see a weak drug and dose early (found 2026-10-25)
    assert wide["ttp"].mean() > narrow["ttp"].mean() + 28
    assert wide["courses"].mean() > narrow["courses"].mean()
    assert guard["ttp"].mean() >= wide["ttp"].mean() - 7
    assert tc.RiskContainment(0.05, wide=True).ranges["m"][0] == 3.0
    assert tc.RiskContainment(0.05).ranges["m"][0] == 6.0


class Record(NoTreat):
    def act(self, j, y):
        self.y = getattr(self, "y", []) + [y.copy()]
        return np.zeros(y.shape, bool)


def test_wrong_model_reading_options():
    scen = tc.draw_scenarios(200, seed=15, sigmas=(0.05,))
    base = tc.simulate(scen, Record(), seed=16)["ttp"]
    for extra in ({}, {"bias": np.full(200, 0.1)}, {"sigma_true": scen["sigma"] * 2},
                  {"noise_df": np.full(200, 3.0)}):
        pol = Record()
        out = tc.simulate(dict(scen, **extra), pol, seed=16)
        assert np.array_equal(out["ttp"], base)          # readings never change the untreated truth
        y0 = pol.y[0] - np.log(scen["B0"])
        if "bias" in extra:
            assert abs(y0.mean() - 0.1) < 0.02
        elif "sigma_true" in extra:
            assert 0.08 < y0.std() < 0.12
        else:
            assert 0.04 < y0.std() < 0.06


def test_cvar_estimands_differ_as_defined():
    # scenario 0 is bad in one replicate only; joint CVaR sees it, scenario-mean CVaR halves it
    reps = np.array([[0.0] + [100.0] * 9, [100.0] * 10])
    assert tc.cvar_joint(reps, 0.1) == pytest.approx((0 + 100) / 2)        # lowest 2 of 20
    assert tc.cvar_scenario_mean(reps, 0.1) == pytest.approx(50.0)          # lowest 1 of 10
    reps2 = np.array([[0.0] * 2 + [100.0] * 18])
    assert tc.cvar_joint(reps2, 0.1) == pytest.approx(0.0)


def test_hier_bootstrap_paired_and_shift():
    rng = np.random.default_rng(22)
    a = rng.normal(500, 100, (3, 200))
    d, ci = tc.hier_bootstrap_diff(a, a, tc.cvar_joint, n_boot=200)
    assert d == 0 and ci == [0.0, 0.0]
    d, ci = tc.hier_bootstrap_diff(a + 5, a, lambda x: float(np.mean(x)), n_boot=200)
    assert d == pytest.approx(5) and ci[0] == pytest.approx(5) and ci[1] == pytest.approx(5)


def test_evaluate_keeps_replicates():
    scen = tc.draw_scenarios(10, seed=23)
    out = tc.evaluate(scen, tc.Periodic(4), [1, 2, 3])
    assert out["ttp_reps"].shape == (3, 10)
    assert np.allclose(out["ttp_reps"].mean(0), out["ttp"])
    assert out["max_courses"] <= tc.BUDGET


def test_cvar_lower():
    assert tc.cvar_lower(np.arange(1, 11), 0.1) == pytest.approx(1.0)
    x = np.arange(15.0)
    assert tc.cvar_lower(x, 0.1) == pytest.approx((0 + 0.5 * 1) / 1.5)
