from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_sim as tc

_SPEC = importlib.util.spec_from_file_location(
    "stage2_208", REPO / "src" / "208_trackc_robust_belief_mpc.py"
)
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)
RobustBeliefMPC = _MODULE.RobustBeliefMPC
robust_endpoint_utility = _MODULE.robust_endpoint_utility


def _policy(*, budget=4, scan_every=2, forecast_cycles=4):
    return RobustBeliefMPC(
        scan_every=scan_every,
        budget=budget,
        course_penalty=5.0,
        n_particles=8,
        prediction_particles=4,
        forecast_cycles=forecast_cycles,
    )


def _scenarios(n=2):
    return tc.draw_scenarios(n, 1234)


def test_robust_utility_uses_worse_endpoint_and_course_penalty():
    got = robust_endpoint_utility(
        [100.0, 500.0],
        [160.0, 200.0],
        [1.0, 2.0],
        [2.0, 1.0],
        tau=365.0,
        course_penalty=5.0,
    )
    np.testing.assert_allclose(got, [95.0, 195.0])


def test_robust_utility_caps_time_at_tau():
    got = robust_endpoint_utility(
        500.0, 400.0, 2.0, 1.0, tau=365.0, course_penalty=5.0
    )
    assert float(got) == pytest.approx(355.0)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"tau": 0},
        {"tau": float("nan")},
        {"course_penalty": -1},
        {"ttp_baseline": -1},
    ],
)
def test_robust_utility_rejects_invalid_inputs(kwargs):
    args = dict(
        ttp_baseline=100.0,
        ttp_rano=120.0,
        courses_baseline=1.0,
    )
    args.update(kwargs)
    with pytest.raises(ValueError):
        robust_endpoint_utility(**args)


def _actions_for_seed(seed):
    scen = _scenarios()
    policy = _policy()
    policy.reset(scen, np.random.default_rng(seed))
    y0 = np.array([-2.5, -2.1])
    actions = [policy.act(0, y0)]
    policy.advance(0, actions[-1])
    actions.append(policy.act(1, np.full(2, np.nan)))
    policy.advance(1, actions[-1])
    actions.append(policy.act(2, np.array([-2.4, -2.0])))
    return np.asarray(actions)


def test_policy_actions_are_deterministic_for_fixed_seed_and_readings():
    np.testing.assert_array_equal(_actions_for_seed(77), _actions_for_seed(77))


def test_policy_does_not_read_true_parameter_arrays_from_scenarios():
    first = _scenarios()
    second = {k: v.copy() for k, v in first.items()}
    for key in ("r", "m", "f0", "c", "eps", "B0"):
        second[key] = second[key] * 0.5
    p1, p2 = _policy(), _policy()
    p1.reset(first, np.random.default_rng(91))
    p2.reset(second, np.random.default_rng(91))
    y = np.array([-2.5, -2.1])
    a1 = p1.act(0, y)
    a2 = p2.act(0, y)
    np.testing.assert_array_equal(a1, a2)
    p1.advance(0, a1)
    p2.advance(0, a2)
    a1b = p1.act(1, np.full(2, np.nan))
    a2b = p2.act(1, np.full(2, np.nan))
    p1.advance(1, a1b)
    p2.advance(1, a2b)
    np.testing.assert_array_equal(
        p1.act(2, np.array([-2.4, -2.0])),
        p2.act(2, np.array([-2.4, -2.0])),
    )


def test_missing_readings_do_not_update_particle_weights():
    scen = _scenarios()
    policy = _policy(scan_every=2, forecast_cycles=4)
    policy.reset(scen, np.random.default_rng(19))
    first = policy.act(0, np.array([-2.5, -2.1]))
    policy.advance(0, first)
    weights_before_missing_cycle = policy.pf.logw.copy()
    missing_action = policy.act(1, np.full(2, np.nan))
    np.testing.assert_array_equal(policy.pf.logw, weights_before_missing_cycle)
    policy.advance(1, missing_action)

    missed_scan_action = policy.act(2, np.array([np.nan, np.nan]))
    np.testing.assert_array_equal(missed_scan_action, [False, False])
    np.testing.assert_array_equal(
        policy.pf.logw[0], weights_before_missing_cycle[0]
    )


def test_policy_never_exceeds_its_course_budget():
    scen = _scenarios()
    policy = _policy(budget=1, scan_every=1, forecast_cycles=2)
    policy.reset(scen, np.random.default_rng(53))
    total = np.zeros(2, dtype=int)
    for j in range(8):
        action = policy.act(j, np.array([-2.5 + j * 0.02, -2.1 + j * 0.02]))
        assert action.dtype == bool
        total += action
        policy.advance(j, action)
    assert np.all(total <= 1)
    np.testing.assert_array_equal(policy.act(8, np.array([-2.0, -1.8])), [False, False])


def test_invalid_scan_or_forecast_dimensions_are_rejected():
    with pytest.raises(ValueError):
        RobustBeliefMPC(scan_every=5)
    with pytest.raises(ValueError):
        RobustBeliefMPC(scan_every=3, forecast_cycles=4)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"scan_every": 2.9},
        {"budget": 3.5},
        {"n_particles": 8.1},
        {"prediction_particles": 2.5},
        {"forecast_cycles": 4.2},
        {"scan_every": True},
    ],
)
def test_fractional_or_boolean_integer_settings_are_rejected(kwargs):
    with pytest.raises(ValueError):
        RobustBeliefMPC(**kwargs)


def _baseline_nadir_policy():
    return RobustBeliefMPC(
        scan_every=3,
        budget=4,
        course_penalty=5.0,
        n_particles=256,
        prediction_particles=8,
        forecast_cycles=6,
    )


def _capture_rano_forecast(monkeypatch, policy):
    observed = []

    def capture(state, rows, j, pattern, nadir_start, t_base_now, t_rano_now):
        observed.append(t_rano_now.copy())
        return np.zeros(rows.size)

    monkeypatch.setattr(policy, "_forecast_pattern", capture)
    return observed


def test_baseline_scan_replaces_prior_nadir_and_stable_rescan_is_not_progression(monkeypatch):
    policy = _baseline_nadir_policy()
    policy.reset(_scenarios(n=1), np.random.default_rng(20261105))
    prior_mean = policy.pf.burden_mean()[0]
    baseline = np.array([np.log(0.45)])

    policy.act(0, baseline)
    posterior_mean = policy.pf.burden_mean()[0]

    assert posterior_mean > prior_mean + 0.10
    assert policy.nadir_est[0] == pytest.approx(posterior_mean)

    # Hold the synthetic patient at the same burden so only nadir bookkeeping is tested.
    monkeypatch.setattr(policy.pf, "propagate", lambda j, treat: None)
    for j in (0, 1, 2):
        policy.advance(j, np.array([False]))
    rano_calls = _capture_rano_forecast(monkeypatch, policy)
    policy.act(3, baseline)

    assert len(rano_calls) == policy.patterns.shape[0]
    assert all(np.isnan(call).all() for call in rano_calls)


def test_first_observed_scan_initializes_nadir_without_inventing_prior_progression(monkeypatch):
    policy = _baseline_nadir_policy()
    policy.reset(_scenarios(n=1), np.random.default_rng(20261106))
    assert np.isnan(policy.nadir_est[0])

    policy.act(0, np.array([np.nan]))
    assert np.isnan(policy.nadir_est[0])

    monkeypatch.setattr(policy.pf, "propagate", lambda j, treat: None)
    for j in (0, 1, 2):
        policy.advance(j, np.array([False]))
    first_reading = np.array([np.log(0.45)])
    rano_calls = _capture_rano_forecast(monkeypatch, policy)
    policy.act(3, first_reading)

    assert np.isfinite(policy.nadir_est[0])
    assert policy.nadir_est[0] == pytest.approx(policy.pf.burden_mean()[0])
    assert all(np.isnan(call).all() for call in rano_calls)
