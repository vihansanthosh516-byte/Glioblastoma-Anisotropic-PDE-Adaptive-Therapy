"""Contract tests for paired Track C signal-policy benchmarking; synthetic only."""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import importlib.util
import sys

import numpy as np
import pytest


REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))
import trackc_sim as tc  # noqa: E402
from trackc_integrated_signal_model import (  # noqa: E402
    SignalModelParameters,
    SignalObservations,
    generate_integrated_signals,
)

SPEC = importlib.util.spec_from_file_location(
    "trackc_paired_signal_benchmark",
    REPO / "src" / "209_trackc_paired_signal_benchmark.py",
)
benchmark = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = benchmark
SPEC.loader.exec_module(benchmark)


def _parameters():
    return SignalModelParameters(q=0.8, a=0.5, b=0.4, d=0.6)


def test_matched_arms_change_only_declared_signal_channel_or_dependence():
    arms = benchmark.matched_signal_arms(_parameters())

    assert set(arms) == {
        "joint",
        "conditional_independence",
        "molecular_only",
        "mri_only",
        "no_signal",
    }
    joint = arms["joint"]
    independent = arms["conditional_independence"]
    assert joint.channels == independent.channels == ("molecular", "mri")
    assert joint.parameters.q == independent.parameters.q == 0.8
    assert joint.parameters.a == independent.parameters.a == 0.5
    assert joint.parameters.b == independent.parameters.b == 0.4
    assert joint.parameters.d == 0.6
    assert independent.parameters.d == 0.0
    assert arms["no_signal"].channels == ()


def test_matched_arms_reject_zero_dependence_primary_model():
    with pytest.raises(ValueError, match="nonzero residual dependence"):
        benchmark.matched_signal_arms(
            SignalModelParameters(q=0.8, a=0.5, b=0.4, d=0.0)
        )


def test_runner_pairs_worlds_signals_seeds_and_hides_latent_truth(monkeypatch):
    scenarios = tc.draw_scenarios(2, seed=20261120, sigmas=(0.05,))
    observations = SignalObservations(
        molecular_report=np.array([True, False]),
        mri_score=np.array([0.2, -0.3]),
    )
    arms = benchmark.matched_signal_arms(_parameters())
    calls = []

    def fake_simulate_two_endpoints(scen, policy, *, seed, stop_endpoint, budget):
        policy.reset(scen, np.random.default_rng(seed))
        arm_index = list(arms).index(policy.arm_name)
        used_observations = policy.policy.signal_observations
        calls.append(
            {
                "arm": policy.arm_name,
                "stop_endpoint": stop_endpoint,
                "seed": seed,
                "budget": budget,
                "scan_every": scen["scan_every"],
                "reset_keys": policy.reset_context_keys,
                "reports": None if used_observations is None else used_observations.molecular_report.copy(),
                "scores": None if used_observations is None else used_observations.mri_score.copy(),
            }
        )
        return {
            "ttp_baseline": np.array([100.0 + arm_index, tc.HORIZON]),
            "ttp_rano": np.array([80.0 + arm_index, tc.HORIZON]),
            "capped_baseline": np.array([False, True]),
            "capped_rano": np.array([False, True]),
            "courses": np.ones(2, dtype=int),
            "stop_endpoint": stop_endpoint,
        }

    monkeypatch.setattr(benchmark.tc, "simulate_two_endpoints", fake_simulate_two_endpoints)
    result = benchmark.run_paired_signal_world(
        world_id="synthetic-world-01",
        scenarios=scenarios,
        observations=observations,
        arms=arms,
        generating_parameters=_parameters(),
        seed=20261121,
        scan_every=2,
        budget=3,
        tau=730.0,
        course_penalty=5.0,
        policy_kwargs={
            "n_particles": 32,
            "prediction_particles": 4,
            "forecast_cycles": 2,
        },
    )

    assert len(calls) == 5
    assert {call["seed"] for call in calls} == {20261121}
    assert {call["stop_endpoint"] for call in calls} == {"either"}
    assert {call["budget"] for call in calls} == {3}
    assert {call["scan_every"] for call in calls} == {2}
    assert {call["reset_keys"] for call in calls} == {("sigma",)}
    for call in calls:
        if call["arm"] in {"joint", "conditional_independence"}:
            np.testing.assert_array_equal(call["reports"], observations.molecular_report)
            np.testing.assert_array_equal(call["scores"], observations.mri_score)
        elif call["arm"] == "molecular_only":
            np.testing.assert_array_equal(call["reports"], observations.molecular_report)
            np.testing.assert_array_equal(call["scores"], np.zeros(2))
        elif call["arm"] == "mri_only":
            np.testing.assert_array_equal(call["reports"], np.zeros(2, dtype=bool))
            np.testing.assert_array_equal(call["scores"], observations.mri_score)
        else:
            assert call["arm"] == "no_signal"
            assert call["reports"] is None and call["scores"] is None

    assert result["world_id"] == "synthetic-world-01"
    assert result["seed"] == 20261121
    assert result["tau"] == 730.0
    assert result["generating_signal_parameters"] == {"q": 0.8, "a": 0.5, "b": 0.4, "d": 0.6}
    assert result["joint_likelihood_matches_generator"] is True
    assert result["policy_settings"]["n_particles"] == 32
    assert result["policy_settings"]["prediction_particles"] == 4
    assert result["policy_settings"]["forecast_cycles"] == 2
    assert result["policy_settings"]["stop_endpoint"] == "either"
    assert result["stop_endpoint"] == "either"
    assert "one shared treatment path" in result["endpoint_treatment_stop_rule"].lower()
    assert set(result["courses"]) == set(arms)
    assert set(result["endpoints"]) == {"baseline", "rano"}
    for arm_index, arm_name in enumerate(arms):
        expected = np.array([75.0 + arm_index, 725.0])
        np.testing.assert_array_equal(result["robust_utility"][arm_name], expected)
        np.testing.assert_array_equal(
            result["endpoints"]["baseline"][arm_name]["ttp"],
            np.array([100.0 + arm_index, tc.HORIZON]),
        )
        np.testing.assert_array_equal(
            result["endpoints"]["baseline"][arm_name]["capped"], [False, True]
        )
        np.testing.assert_array_equal(
            result["endpoints"]["rano"][arm_name]["ttp"],
            np.array([80.0 + arm_index, tc.HORIZON]),
        )
        np.testing.assert_array_equal(
            result["endpoints"]["rano"][arm_name]["capped"], [False, True]
        )
    assert "scan_every" not in scenarios


def test_runner_rejects_primary_pair_with_different_marginal_likelihoods():
    arms = benchmark.matched_signal_arms(_parameters())
    arms["conditional_independence"] = replace(
        arms["conditional_independence"],
        parameters=SignalModelParameters(q=0.7, a=0.5, b=0.4, d=0.0),
    )

    with pytest.raises(ValueError, match="q, a, and b must match"):
        benchmark.run_paired_signal_world(
            world_id="synthetic-world-02",
            scenarios=tc.draw_scenarios(1, seed=20261122),
            observations=SignalObservations(np.array([True]), np.array([0.0])),
            arms=arms,
            generating_parameters=_parameters(),
            seed=20261123,
            scan_every=1,
            budget=2,
        )


def test_runner_rejects_fractional_course_counts():
    with pytest.raises(ValueError, match="course counts must be integers"):
        benchmark._validate_outcome(
            {
                "ttp_baseline": np.array([10.0]),
                "ttp_rano": np.array([8.0]),
                "capped_baseline": np.array([False]),
                "capped_rano": np.array([False]),
                "courses": np.array([0.5]),
                "stop_endpoint": "either",
            },
            n_patients=1,
            budget=2,
        )


def test_runner_smoke_executes_all_arms_on_one_tiny_synthetic_world():
    scenarios = tc.draw_scenarios(2, seed=20261124, sigmas=(0.05,))
    observations, _ = generate_integrated_signals(
        2,
        q=0.8,
        a=0.5,
        b=0.4,
        d=0.6,
        seed=20261125,
        f0=scenarios["f0"],
    )
    result = benchmark.run_paired_signal_world(
        world_id="synthetic-smoke-only",
        scenarios=scenarios,
        observations=observations,
        arms=benchmark.matched_signal_arms(_parameters()),
        generating_parameters=_parameters(),
        seed=20261126,
        scan_every=1,
        budget=2,
        policy_kwargs={
            "n_particles": 32,
            "prediction_particles": 4,
            "forecast_cycles": 1,
        },
    )

    assert result["analysis_type"] == "synthetic_virtual_world"
    assert set(result["endpoints"]) == {"baseline", "rano"}
    assert set(result["endpoints"]["baseline"]) == set(benchmark.matched_signal_arms(_parameters()))
    assert set(result["endpoints"]["rano"]) == set(benchmark.matched_signal_arms(_parameters()))
    for arm_name in result["robust_utility"]:
        assert result["robust_utility"][arm_name].shape == (2,)
        assert np.isfinite(result["robust_utility"][arm_name]).all()


def test_runner_rejects_world_with_mismatched_patient_vector_shapes():
    scenarios = tc.draw_scenarios(2, seed=20261127, sigmas=(0.05,))
    scenarios["r"] = scenarios["r"][:1]

    with pytest.raises(ValueError, match="r must have shape"):
        benchmark.run_paired_signal_world(
            world_id="synthetic-world-bad-shape",
            scenarios=scenarios,
            observations=SignalObservations(np.array([True, False]), np.array([0.1, -0.1])),
            arms=benchmark.matched_signal_arms(_parameters()),
            generating_parameters=_parameters(),
            seed=20261129,
            scan_every=1,
            budget=2,
        )
