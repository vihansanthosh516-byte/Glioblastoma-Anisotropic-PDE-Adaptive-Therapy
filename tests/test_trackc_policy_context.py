from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import trackc_sim as tc

_SPEC = importlib.util.spec_from_file_location(
    "trackc_policy_benchmark_context_test",
    SRC_DIR / "211_trackc_policy_benchmark.py",
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("could not load Track C policy benchmark")
benchmark = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = benchmark
_SPEC.loader.exec_module(benchmark)


class ResetKeyRecorder:
    uses_filter = False

    def reset(self, scen, rng):
        del rng
        self.reset_keys = set(scen)
        self.sigma_received = np.asarray(scen["sigma"]).copy()

    def act(self, cycle, readings):
        del cycle
        return np.zeros(np.asarray(readings).shape, dtype=bool)

    def advance(self, cycle, treatment):
        del cycle, treatment


def test_two_endpoint_simulator_passes_only_sigma_to_policy_by_default():
    scen = tc.draw_scenarios(2, seed=881001)
    policy = ResetKeyRecorder()

    tc.simulate_two_endpoints(scen, policy, seed=881002, budget=2)

    assert policy.reset_keys == {"sigma"}
    np.testing.assert_array_equal(policy.sigma_received, scen["sigma"])


def test_two_endpoint_simulator_full_context_is_an_explicit_opt_in():
    scen = tc.draw_scenarios(1, seed=881003)
    policy = ResetKeyRecorder()

    tc.simulate_two_endpoints(
        scen, policy, seed=881004, budget=2, policy_context="full"
    )

    assert policy.reset_keys == set(scen)


def test_two_endpoint_simulator_rejects_unknown_policy_context_before_reset():
    scen = tc.draw_scenarios(1, seed=881005)
    policy = ResetKeyRecorder()

    with pytest.raises(ValueError, match="policy_context"):
        tc.simulate_two_endpoints(
            scen, policy, seed=881006, budget=2, policy_context="latent"
        )

    assert not hasattr(policy, "reset_keys")


def test_fixed_schedule_and_course_cap_policies_accept_sigma_only_context():
    context = {"sigma": np.asarray([0.1, 0.2])}
    rng = np.random.default_rng(881007)
    schedule = np.zeros(tc.N_CYCLES, dtype=bool)

    fixed = benchmark.FixedSchedulePolicy(schedule)
    fixed.reset(context, rng)
    assert fixed._active.shape == (2, tc.N_CYCLES)

    holiday = tc.HolidayContainment(
        p_star=0.1, horizon=1, max_run=2, n_particles=8
    )
    capped = benchmark.CourseCappedPolicy(holiday, max_courses=2)
    capped.reset(context, np.random.default_rng(881008))
    assert capped._used.shape == (2,)
    assert holiday.run.shape == (2,)


def test_course_cap_strips_extra_scenario_fields_before_reset():
    scen = tc.draw_scenarios(2, seed=881011)
    policy = ResetKeyRecorder()
    capped = benchmark.CourseCappedPolicy(policy, max_courses=2)

    capped.reset(scen, np.random.default_rng(881012))

    assert policy.reset_keys == {"sigma"}
    np.testing.assert_array_equal(policy.sigma_received, scen["sigma"])


def test_benchmark_evaluation_requests_observation_only_context(monkeypatch):
    captured = {}

    def fake_simulate_two_endpoints(scen, policy, *, seed, stop_endpoint, budget,
                                    course_cutoff_days, policy_context):
        del policy, seed, stop_endpoint, budget, course_cutoff_days
        captured["policy_context"] = policy_context
        n = len(scen["r"])
        return {
            "ttp_baseline": np.full(n, 50.0),
            "capped_baseline": np.zeros(n, dtype=bool),
            "ttp_rano": np.full(n, 40.0),
            "capped_rano": np.zeros(n, dtype=bool),
            "courses": np.ones(n),
            "courses_through_cutoff": np.ones(n),
            "stop_endpoint": "observed_two_scan",  # integrated v2 refuses any other rule
            "observed_stop_day": np.full(n, np.nan),  # no observed stop in this fake
        }

    monkeypatch.setattr(benchmark.tc, "simulate_two_endpoints", fake_simulate_two_endpoints)
    scen = tc.draw_scenarios(2, seed=881009)
    settings = benchmark.PolicyBenchmarkSettings(
        tau_days=730, mpc_particles=8, prediction_particles=2, q_episodes=1
    )

    benchmark._evaluate_policy(
        ResetKeyRecorder(), scen, (881010,), settings
    )

    assert captured["policy_context"] == "observable"
