"""Focused contract tests for the Track C tabular Q-learning baseline."""
from __future__ import annotations

from collections.abc import Iterator, Mapping
from pathlib import Path
import importlib.util
import sys

import numpy as np


REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
SPEC = importlib.util.spec_from_file_location(
    "trackc_qlearning_baseline", REPO / "src" / "207_trackc_qlearning_baseline.py"
)
qlearn = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = qlearn
SPEC.loader.exec_module(qlearn)


def test_state_binning_uses_observed_burden_trend_time_and_course_count():
    state = qlearn.bin_state(
        cycle=3,
        log_relative_burden=np.log(0.6),
        trend_delta=-0.2,
        courses_used=7,
        course_cap=2,
    )
    assert state == (3, 1, 0, 2)
    assert qlearn.bin_state(0, 0.0, 0.0, 0, course_cap=2) == (0, 3, 1, 0)


def test_bellman_update_applies_one_step_target():
    q = np.zeros(qlearn.q_shape(course_cap=2))
    state = qlearn.bin_state(1, 0.0, 0.0, 0, course_cap=2)
    next_state = qlearn.bin_state(2, 0.1, 0.2, 1, course_cap=2)
    q[state + (0,)] = 2.0
    q[next_state + (0,)] = 1.0
    q[next_state + (1,)] = 4.0

    updated = qlearn.bellman_update(
        q, state, 0, reward=3.0, next_state=next_state,
        learning_rate=0.5, discount=0.9,
    )

    assert np.isclose(updated, 4.3)
    assert q[state + (0,)] == updated


def test_policy_is_deterministic_and_enforces_its_course_cap():
    q = np.zeros(qlearn.q_shape(course_cap=1))
    initial = qlearn.bin_state(0, 0.0, 0.0, 0, course_cap=1)
    q[initial + (1,)] = 5.0
    policy = qlearn.QTablePolicy(q, course_cap=1, scan_every=3)

    policy.reset({}, np.random.default_rng(4))
    first = policy.act(0, np.array([0.0]))
    policy.reset({}, np.random.default_rng(99))
    first_again = policy.act(0, np.array([0.0]))
    assert first_again.tolist() == first.tolist()

    policy.reset({}, np.random.default_rng(4))
    first = policy.act(0, np.array([0.0]))
    policy.advance(0, first)
    second = policy.act(1, np.array([np.nan]))

    assert first.tolist() == [True]
    assert second.tolist() == [False]
    assert policy.scan_every == 3


class _NoTruthAccess(Mapping):
    forbidden = {"r", "m", "f0", "c", "eps", "B0", "mu", "nu"}

    def __getitem__(self, key):
        if key in self.forbidden:
            raise AssertionError(f"policy accessed true patient parameter {key}")
        raise KeyError(key)

    def __iter__(self) -> Iterator:
        return iter(())

    def __len__(self) -> int:
        return 0


def test_policy_reset_never_reads_true_patient_parameters():
    policy = qlearn.QTablePolicy(
        np.zeros(qlearn.q_shape(course_cap=2)), course_cap=2, scan_every=3
    )
    policy.reset(_NoTruthAccess(), np.random.default_rng(9))
    chosen = policy.act(0, np.array([0.0, 0.4]))
    policy.advance(0, chosen)
    policy.act(1, np.array([np.nan, np.nan]))

    assert chosen.shape == (2,)


def test_training_uses_terminal_ttp_reward_and_structured_seed(monkeypatch):
    captured = []

    def fake_simulate(scenarios, policy, seed, *, endpoint, budget):
        captured.append((list(seed), endpoint, budget))
        policy.reset(scenarios, np.random.default_rng(seed))
        for cycle, reading in ((0, np.array([0.0])), (1, np.array([np.nan]))):
            action = policy.act(cycle, reading)
            policy.advance(cycle, action)
        return {"ttp": np.array([56.0])}

    monkeypatch.setattr(qlearn.tc, "simulate", fake_simulate)
    fitted = qlearn.train_q_learning(
        {"r": np.array([0.01])},
        episodes=1,
        course_cap=1,
        course_penalty=0.0,
        scan_every=1,
        learning_rate=1.0,
        epsilon=0.0,
        min_epsilon=0.0,
        seed=[5, 7],
    )

    terminal_state = qlearn.bin_state(1, 0.0, 0.0, 0, course_cap=1)
    initial_state = qlearn.bin_state(0, 0.0, 0.0, 0, course_cap=1)
    assert fitted.q_values[terminal_state + (0,)] == 56.0
    assert fitted.q_values[initial_state + (0,)] == 56.0
    assert captured == [([5, 7, 0], "rano", 1)]


def test_training_requires_undiscounted_objective():
    with np.testing.assert_raises(ValueError):
        qlearn.train_q_learning({"r": np.array([0.01])}, discount=0.9)


def test_evaluation_forwards_fitted_scan_interval_and_course_cap(monkeypatch):
    policy = qlearn.QTablePolicy(
        np.zeros(qlearn.q_shape(course_cap=2)), course_cap=2, scan_every=4
    )
    captured = {}

    def fake_simulate(scenarios, actual_policy, seed, *, endpoint, budget):
        captured.update(scenarios=scenarios, policy=actual_policy, endpoint=endpoint,
                        budget=budget, seed=seed)
        return {"ttp": np.array([28.0]), "courses": np.array([1])}

    monkeypatch.setattr(qlearn.tc, "simulate", fake_simulate)
    result = qlearn.evaluate_q_policy({"r": np.array([0.01])}, policy, [12, 13], endpoint="baseline")

    assert captured["scenarios"]["scan_every"] == 4
    assert captured["budget"] == 2
    assert captured["endpoint"] == "baseline"
    assert result["ttp"].tolist() == [28.0]


def test_shared_training_uses_earliest_endpoint_on_one_path(monkeypatch):
    captured = {}

    def fake_simulate_two_endpoints(scenarios, policy, seed, *, stop_endpoint, budget):
        captured.update(scenarios=scenarios, stop_endpoint=stop_endpoint,
                        budget=budget, seed=list(seed))
        policy.reset(scenarios, np.random.default_rng(seed))
        for cycle in range(9):
            action = policy.act(cycle, np.array([0.0]))
            policy.advance(cycle, action)
        return {
            "ttp_baseline": np.array([500.0]),
            "ttp_rano": np.array([200.0]),
            "courses": np.array([0]),
        }

    monkeypatch.setattr(qlearn.tc, "simulate_two_endpoints", fake_simulate_two_endpoints)
    fitted = qlearn.train_q_learning_shared_endpoints(
        {"r": np.array([0.01])}, episodes=1, course_cap=1, course_penalty=0.0,
        scan_every=1, tau_days=365.0, learning_rate=1.0, epsilon=0.0,
        min_epsilon=0.0, seed=[5, 7],
    )

    event_state = qlearn.bin_state(7, 0.0, 0.0, 0, course_cap=1)
    after_event_state = qlearn.bin_state(8, 0.0, 0.0, 0, course_cap=1)
    assert fitted.q_values[event_state + (0,)] == 200.0
    assert fitted.q_values[after_event_state + (0,)] == 0.0
    assert captured["stop_endpoint"] == "either"
    assert captured["budget"] == 1
    assert captured["seed"] == [5, 7, 0]
    assert captured["scenarios"]["scan_every"] == 1


def test_shared_evaluator_forwards_scan_and_budget(monkeypatch):
    policy = qlearn.QTablePolicy(
        np.zeros(qlearn.q_shape(course_cap=2)), course_cap=2, scan_every=4
    )
    captured = {}

    def fake_simulate_two_endpoints(scenarios, actual_policy, seed, *, stop_endpoint, budget):
        captured.update(scenarios=scenarios, policy=actual_policy, stop_endpoint=stop_endpoint,
                        budget=budget, seed=seed)
        return {"ttp_baseline": np.array([56.0]), "ttp_rano": np.array([84.0]),
                "courses": np.array([1])}

    monkeypatch.setattr(qlearn.tc, "simulate_two_endpoints", fake_simulate_two_endpoints)
    result = qlearn.evaluate_q_policy_shared_endpoints({"r": np.array([0.01])}, policy, [12, 13])

    assert captured["scenarios"]["scan_every"] == 4
    assert captured["budget"] == 2
    assert captured["stop_endpoint"] == "either"
    assert captured["policy"] is policy
    assert result["ttp_baseline"].tolist() == [56.0]
    assert result["ttp_rano"].tolist() == [84.0]


def test_shared_endpoint_trainer_and_evaluator_run_on_synthetic_scenarios():
    scenarios = qlearn.tc.draw_scenarios(3, seed=20261106, sigmas=(0.05,))
    policy = qlearn.train_q_learning_shared_endpoints(
        scenarios, episodes=1, course_cap=2, course_penalty=5.0,
        scan_every=4, tau_days=365.0, learning_rate=0.1,
        epsilon=0.0, min_epsilon=0.0, seed=[20261106, 1],
    )
    outcome = qlearn.evaluate_q_policy_shared_endpoints(
        scenarios, policy, seed=[20261106, 2]
    )

    assert outcome["stop_endpoint"] == "either"
    assert outcome["ttp_baseline"].shape == (3,)
    assert outcome["ttp_rano"].shape == (3,)
    assert outcome["courses"].shape == (3,)
    assert np.isfinite(outcome["ttp_baseline"]).all()
    assert np.isfinite(outcome["ttp_rano"]).all()
    assert np.all((outcome["courses"] >= 0) & (outcome["courses"] <= 2))
