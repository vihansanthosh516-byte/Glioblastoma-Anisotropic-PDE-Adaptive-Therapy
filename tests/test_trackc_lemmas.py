"""Exact identities of the Track C model (src/trackc_sim.py), plus brute-force checks of what they imply.
Synthetic parameters only; these check algebra and model behaviour, not any patient or drug.

Setting: a_rs = a_sr = 1 (the default), B = S + R, f = R / B, k = m r, u in {0, 1} drug on/off.

Lemma 1 (exact, pointwise):
    df/dt = f (1 - f) [ k (1 - eps) u - c r (1 - B) ] + mu u (1 - f) - nu (1 - u) f
  so  d logit f / dt = k (1 - eps) u - c r (1 - B) + mu u / f - nu (1 - u) / (1 - f).
Lemma 2 (c = mu = nu = 0):  logit f(t) = logit f0 + k (1 - eps) D(t),  D = drug-days so far.
  The resistant fraction depends on the schedule only through the total drug-days, not their timing.
Lemma 3 (no drug, c = 0): B is logistic and  TTP = (1 / r) ln[ 1.2 (1 - B0) / (1 - 1.2 B0) ].
  With c > 0 the growth is slower, so this is a lower bound.
Consequence (small burden, B << 1, c = mu = nu = 0): ln B(t) = ln B0 + r t - K(D(t)) with K increasing in D, and
  front-loaded dosing has the largest D(t) at every t, so it maximises the time to progression. Feedback can
  gain only through crowding (B near the carrying capacity K = 1), through the cost c, or through switching.
  The brute-force tests below check this behaviour on short horizons: the gain over front-loaded dosing is
  small when B0 is small and large when B0 is large. They are behaviour checks of the model, not theorems.
"""
from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_oracle as to  # noqa: E402
import trackc_sim as tc  # noqa: E402


def logit(x):
    return np.log(x / (1.0 - x))


def test_lemma1_derivative_identity_at_random_states():
    rng = np.random.default_rng(11)
    n = 4000
    p = {"r": rng.uniform(0.005, 0.05, n), "m": rng.uniform(3.0, 20.0, n), "c": rng.uniform(0.0, 0.25, n),
         "eps": rng.uniform(0.0, 0.2, n), "a_rs": np.ones(n), "a_sr": np.ones(n),
         "mu": rng.uniform(0.0, 0.05, n), "nu": rng.uniform(0.0, 0.05, n)}
    S, R = rng.uniform(0.001, 0.6, n), rng.uniform(0.001, 0.4, n)
    u = rng.integers(0, 2, n).astype(float)
    dS, dR = tc._deriv(S, R, p, u)
    B = S + R
    f = R / B
    fdot = (dR * B - R * (dS + dR)) / B ** 2
    k = p["m"] * p["r"]
    formula = (f * (1 - f) * (k * (1 - p["eps"]) * u - p["c"] * p["r"] * (1 - B))
               + p["mu"] * u * (1 - f) - p["nu"] * (1 - u) * f)
    assert np.allclose(fdot, formula, rtol=1e-10, atol=1e-13)


def _scalar_params(**kw):
    base = {"r": 0.01, "m": 10.0, "c": 0.0, "eps": 0.1, "a_rs": 1.0, "a_sr": 1.0}
    base.update(kw)
    return {key: np.array([val], float) for key, val in base.items()}


def test_lemma2_resistant_fraction_depends_only_on_drug_days():
    p = _scalar_params()
    f0, B0 = 0.02, 0.2
    k = float(p["m"][0] * p["r"][0])
    schedules = {"front": [1] * 6 + [0] * 8, "alternate": [1, 0] * 6 + [0, 0], "late": [0] * 7 + [1] * 6 + [0]}
    finals = {}
    for name, sched in schedules.items():
        S, R = np.array([(1 - f0) * B0]), np.array([f0 * B0])
        for j, a in enumerate(sched):
            S, R = tc.run_cycle(S, R, p, np.array([float(a)]), j * tc.CYCLE_DAYS)
        finals[name] = float(np.log(R[0] / S[0]))
    drug_days = 6 * tc.COURSE_DAYS
    expected = float(logit(f0)) + k * (1 - p["eps"][0]) * drug_days
    for name, val in finals.items():
        assert val == pytest.approx(expected, abs=1e-6), name
    assert max(finals.values()) - min(finals.values()) < 1e-6


def _rk4_aug(p, state, u, h):
    """RK4 for (S, R, I1, I2, I3, D) with I1' = 1 - B, I2' = u / f, I3' = (1 - u) / (1 - f), D' = u."""
    def deriv(x):
        S, R = x[0], x[1]
        dS, dR = tc._deriv(S, R, p, u)
        B = S + R
        f = R / B
        return np.array([dS, dR, 1.0 - B, u / f, (1.0 - u) / (1.0 - f), u], float)
    k1 = deriv(state)
    k2 = deriv(state + 0.5 * h * k1)
    k3 = deriv(state + 0.5 * h * k2)
    k4 = deriv(state + h * k3)
    return state + h / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


@pytest.mark.parametrize("c,mu,nu", [(0.0, 0.0, 0.0), (0.2, 0.0, 0.0), (0.1, 0.004, 0.0), (0.15, 0.003, 0.01)])
def test_lemma1_integrated_form(c, mu, nu):
    kw = {"c": c}
    if mu or nu:
        kw.update(mu=mu, nu=nu)
    p = {key: float(val[0]) for key, val in _scalar_params(**kw).items()}   # plain floats: the state is a vector
    f0, B0 = 0.05, 0.15
    state = np.array([(1 - f0) * B0, f0 * B0, 0.0, 0.0, 0.0, 0.0])
    h = 0.05
    for day in range(0, 230):                          # 230 days; drug on for days 0-4 of cycles 0, 1, 3, 5
        cyc, in_cycle = divmod(day, tc.CYCLE_DAYS)
        u = 1.0 if (cyc in (0, 1, 3, 5) and in_cycle < tc.COURSE_DAYS) else 0.0
        for _ in range(int(round(1.0 / h))):
            state = _rk4_aug(p, state, u, h)
    S, R, I1, I2, I3, D = state
    k = p["m"] * p["r"]
    pred = logit(f0) + k * (1 - p["eps"]) * D - c * p["r"] * I1 + mu * I2 - nu * I3
    assert float(np.log(R / S)) == pytest.approx(float(pred), abs=1e-7)
    assert D == pytest.approx(4 * tc.COURSE_DAYS)


def _untreated_ttp(p, f0, B0, level_factor=1.2):
    S, R = np.array([(1 - f0) * B0]), np.array([f0 * B0])
    ttp = np.full(1, np.nan)
    for j in range(tc.N_CYCLES):
        S, R = tc.run_cycle(S, R, p, np.zeros(1), j * tc.CYCLE_DAYS, prog_level=np.array([level_factor * B0]), ttp=ttp)
        if not np.isnan(ttp).any():
            break
    return float(ttp[0])


@pytest.mark.parametrize("r,B0", [(0.01, 0.1), (0.002, 0.05), (0.05, 0.4), (0.001, 0.3), (0.02, 0.45)])
def test_lemma3_untreated_time_to_progression(r, B0):
    exact = np.log(1.2 * (1 - B0) / (1 - 1.2 * B0)) / r
    p = _scalar_params(r=r, c=0.0)
    for f0 in (0.0001, 0.3):                            # with c = 0 the resistant share does not matter
        assert _untreated_ttp(p, f0, B0) == pytest.approx(exact, abs=0.01)
    slower = _scalar_params(r=r, c=0.25)
    assert _untreated_ttp(slower, 0.5, B0) >= exact - 0.01   # cost only slows growth


def _gain_over_front_loaded(c_rng, b0_rng, n, seed, n_cycles=10, budget=4):
    rng = dict(tc.RANGES, r=(0.02, 0.05, "log"), c=c_rng, B0=b0_rng)
    scen = tc.draw_scenarios(n, seed, ranges=rng)
    rows = np.array([np.isin(np.arange(n_cycles), cc) for m in range(budget + 1)
                     for cc in itertools.combinations(range(n_cycles), m)])
    front = np.zeros(n_cycles, bool)
    front[:budget] = True
    gains = []
    for i in range(n):
        p = to._scalar(scen, i)
        gains.append(to.rollout(p, rows, n_cycles).max() - to.rollout(p, front[None], n_cycles)[0])
    return np.array(gains)


def test_front_loaded_is_nearly_optimal_without_cost_when_burden_is_small():
    small = _gain_over_front_loaded((0.0, 0.0, "lin"), (0.05, 0.10, "lin"), n=8, seed=5)
    assert (small >= -1e-6).all()                       # brute force includes front-loaded, so never negative
    assert small.max() < 10.0 and small.mean() < 4.0    # observed on 16 scenarios: mean 1.3 d, max 4.2 d


def test_gain_from_perfect_timing_grows_with_starting_burden():
    small = _gain_over_front_loaded((0.0, 0.0, "lin"), (0.05, 0.10, "lin"), n=8, seed=5)
    large = _gain_over_front_loaded((0.0, 0.0, "lin"), (0.40, 0.50, "lin"), n=8, seed=5)
    assert large.mean() > 8.0                           # observed on 16 scenarios: mean 21.9 d
    assert large.mean() > 3.0 * small.mean()
