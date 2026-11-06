"""Contract tests for the paired synthetic policy benchmark."""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import trackc_sim as tc  # noqa: E402
from trackc_integrated_signal_model import (  # noqa: E402
    SignalModelParameters,
    generate_integrated_signals,
)

SOURCE = REPO / "src" / "211_trackc_policy_benchmark.py"
SPEC = importlib.util.spec_from_file_location("trackc_policy_benchmark", SOURCE)
assert SPEC is not None and SPEC.loader is not None
benchmark = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = benchmark
SPEC.loader.exec_module(benchmark)


def _quiet_scenario():
    """One slow synthetic patient whose late schedule falls after 730 days."""
    return {
        "r": np.asarray([0.0001]),
        "m": np.asarray([6.0]),
        "f0": np.asarray([0.02]),
        "c": np.asarray([0.0]),
        "eps": np.asarray([0.0]),
        "B0": np.asarray([0.05]),
        "sigma": np.asarray([0.1]),
        "a_sr": np.asarray([1.0]),
        "a_rs": np.asarray([1.0]),
        "scan_every": np.asarray([3]),
    }


def test_restricted_course_count_excludes_courses_starting_at_or_after_cutoff():
    scenario = _quiet_scenario()
    schedule = np.zeros(tc.N_CYCLES, dtype=bool)
    schedule[26:32] = True
    outcome = tc.simulate_two_endpoints(
        scenario,
        benchmark.FixedSchedulePolicy(schedule),
        seed=123,
        stop_endpoint="either",
        budget=6,
        course_cutoff_days=730.0,
    )

    assert outcome["courses"].tolist() == [6]
    assert outcome["courses_through_cutoff"].tolist() == [1]


def test_restricted_utility_uses_courses_through_tau_not_later_courses():
    outcome = {
        "ttp_baseline": np.asarray([500.0, 900.0]),
        "ttp_rano": np.asarray([450.0, 1_000.0]),
        "courses": np.asarray([9, 10]),
        "courses_through_cutoff": np.asarray([2, 3]),
    }

    utility = benchmark.robust_utility(
        outcome, tau_days=730.0, course_penalty=5.0
    )

    np.testing.assert_allclose(utility, [440.0, 715.0])


def test_conditional_independence_control_preserves_candidate_qab_support():
    candidates = (
        SignalModelParameters(q=0.6, a=0.2, b=0.4, d=-0.3),
        SignalModelParameters(q=0.6, a=0.2, b=0.4, d=0.0),
        SignalModelParameters(q=0.8, a=-0.1, b=0.7, d=0.3),
        SignalModelParameters(q=0.8, a=-0.1, b=0.7, d=0.0),
    )

    control = benchmark.conditional_independence_grid(candidates)

    assert {(row.q, row.a, row.b) for row in control} == {
        (row.q, row.a, row.b) for row in candidates
    }
    assert all(row.d == 0.0 for row in control)


def test_full_policy_harness_runs_all_arms_on_shared_endpoint_path(monkeypatch):
    received_f0_bounds = []
    original_loader = benchmark._load_source_module

    def recording_loader(module_name, filename):
        module = original_loader(module_name, filename)
        if filename == "214_trackc_loo_model_average.py":
            original_controller = module.LOOModelAveragedBeliefMPC

            def recording_controller(*args, **kwargs):
                received_f0_bounds.append(tuple(kwargs["f0_bounds"]))
                return original_controller(*args, **kwargs)

            monkeypatch.setattr(module, "LOOModelAveragedBeliefMPC", recording_controller)
        return module

    monkeypatch.setattr(benchmark, "_load_source_module", recording_loader)
    selection = tc.draw_scenarios(2, seed=2026110601, sigmas=(0.05,))
    evaluation = tc.draw_scenarios(2, seed=2026110602, sigmas=(0.05,))
    observations, _truth = generate_integrated_signals(
        2,
        q=0.6,
        a=0.2,
        b=0.4,
        d=0.3,
        seed=2026110603,
        f0=evaluation["f0"],
    )
    models = (
        SignalModelParameters(q=0.6, a=0.2, b=0.4, d=0.3),
    )
    settings = benchmark.PolicyBenchmarkSettings(
        course_budget=6,
        tau_days=84.0,
        n_f0_nodes=8,
        mpc_particles=2,
        prediction_particles=1,
        forecast_cycles=3,
        q_episodes=1,
    )

    result = benchmark.run_shared_policy_world(
        world_id="tiny-mechanics-test",
        selection_scenarios=selection,
        evaluation_scenarios=evaluation,
        evaluation_observations=observations,
        model_parameters=models,
        selection_seed=[2026110604],
        q_training_seed=[2026110605],
        evaluation_seed=[2026110606],
        settings=settings,
    )

    expected = {
        "loo_dependence_aware_mpc",
        "loo_conditional_independence_mpc",
        "no_signal_mpc",
        "q_learning",
        "selection_best_fixed",
        "standard_six_course",
    }
    assert set(result["arm_outcomes"]) == expected
    assert result["primary_paired_contrast"]["robust_utility_difference"].shape == (2,)
    for arm in result["arm_outcomes"].values():
        assert arm["ttp_baseline"].shape == (2,)
        assert arm["ttp_rano"].shape == (2,)
        assert arm["courses"].shape == (2,)
        assert arm["courses_through_cutoff"].shape == (2,)
        assert (arm["courses_through_cutoff"] <= arm["courses"]).all()
        assert (arm["courses"] <= settings.course_budget).all()
    assert result["analysis_type"] == "synthetic_virtual_world_development"
    assert result["settings"]["signal_f0_bounds"] == tuple(map(float, tc.RANGES["f0"][:2]))
    assert received_f0_bounds == [settings.signal_f0_bounds, settings.signal_f0_bounds]


def test_harness_requires_budget_that_can_include_six_course_comparator():
    selection = tc.draw_scenarios(1, seed=2026110607, sigmas=(0.05,))
    evaluation = tc.draw_scenarios(1, seed=2026110608, sigmas=(0.05,))
    observations, _truth = generate_integrated_signals(
        1, q=0.6, a=0.2, b=0.4, d=0.3, seed=2026110609, f0=evaluation["f0"]
    )
    settings = benchmark.PolicyBenchmarkSettings(
        course_budget=5,
        tau_days=84.0,
        n_f0_nodes=8,
        mpc_particles=2,
        prediction_particles=1,
        forecast_cycles=3,
        q_episodes=1,
    )

    with pytest.raises(ValueError, match="at least 6"):
        benchmark.run_shared_policy_world(
            world_id="invalid-budget",
            selection_scenarios=selection,
            evaluation_scenarios=evaluation,
            evaluation_observations=observations,
            model_parameters=(SignalModelParameters(0.6, 0.2, 0.4, 0.3),),
            selection_seed=[2026110610],
            q_training_seed=[2026110611],
            evaluation_seed=[2026110612],
            settings=settings,
        )


def test_harness_rejects_f0_values_outside_declared_signal_prior_before_fitting():
    selection = tc.draw_scenarios(2, seed=2026110613, sigmas=(0.05,))
    evaluation = tc.draw_scenarios(2, seed=2026110614, sigmas=(0.05,))
    evaluation["f0"][0] = 0.12
    observations, _truth = generate_integrated_signals(
        2,
        q=0.6,
        a=0.2,
        b=0.4,
        d=0.3,
        seed=2026110615,
        f0=evaluation["f0"],
    )

    with pytest.raises(ValueError, match="outside the declared signal_f0_bounds"):
        benchmark.run_shared_policy_world(
            world_id="f0-prior-mismatch",
            selection_scenarios=selection,
            evaluation_scenarios=evaluation,
            evaluation_observations=observations,
            model_parameters=(SignalModelParameters(0.6, 0.2, 0.4, 0.3),),
            selection_seed=[2026110616],
            q_training_seed=[2026110617],
            evaluation_seed=[2026110618],
        )
