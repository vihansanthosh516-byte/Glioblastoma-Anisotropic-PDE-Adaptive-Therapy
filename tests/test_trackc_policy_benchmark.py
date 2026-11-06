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


def test_baseline_catalog_covers_scan_rules_and_state_estimation_with_provenance():
    settings = benchmark.PolicyBenchmarkSettings(
        mpc_particles=7,
        prediction_particles=1,
        n_f0_nodes=8,
        forecast_cycles=3,
        q_episodes=1,
    )

    catalog = benchmark.baseline_policy_catalog(settings)

    assert "raw_zhang_1.0_0.5" in catalog
    assert "nadir_k3_on1.25" in catalog
    assert "resp_k3_th0.75" in catalog
    assert "filtered_threshold_pilot2" in catalog
    assert "pf_wide_containment" in catalog
    assert "pf_tail_trigger_candidate" in catalog
    assert "pf_holiday_containment" in catalog
    assert catalog["raw_zhang_1.0_0.5"]["family"] == "scan_rule"
    assert catalog["pf_wide_containment"]["family"] == "state_estimation"
    filter_grid = [name for name in catalog if name.startswith("filtered_threshold_grid_")]
    assert len(filter_grid) == 128
    assert all(catalog[name]["family"] == "filtered_threshold_grid" for name in filter_grid)
    assert all(not catalog[name]["evaluate_all"] for name in filter_grid)
    filtered_settings = {
        (catalog[name]["factory"]().a_on, catalog[name]["factory"]().a_off,
         catalog[name]["factory"]().a_R, catalog[name]["factory"]().b)
        for name in filter_grid
    }
    assert len(filtered_settings) == 128
    raw_grid = [name for name in catalog if name.startswith("raw_threshold_grid_")]
    assert len(raw_grid) == 24
    assert all(catalog[name]["family"] == "raw_threshold_grid" for name in raw_grid)
    assert all(not catalog[name]["evaluate_all"] for name in raw_grid)
    raw_settings = {
        (catalog[name]["factory"]().a_on, catalog[name]["factory"]().a_off)
        for name in raw_grid
    }
    assert len(raw_settings) == 24
    assert "pilot_results.json" in catalog["filtered_threshold_pilot2"]["provenance"]
    assert all(spec["configuration"] and spec["provenance"] for spec in catalog.values())
    assert catalog["pf_wide_containment"]["factory"]().n_particles == 7


def test_course_cap_wrapper_stops_policy_after_arm_specific_limit():
    class AlwaysTreat:
        uses_filter = False

        def reset(self, scenarios, rng):
            del scenarios, rng

        def act(self, cycle, observations):
            del cycle
            return np.ones(observations.shape, dtype=bool)

        def advance(self, cycle, treat):
            del cycle, treat

    policy = benchmark.CourseCappedPolicy(AlwaysTreat(), max_courses=2)
    policy.reset({"r": np.ones(1), "sigma": np.ones(1)}, np.random.default_rng(1))
    decisions = []
    for cycle in range(4):
        action = policy.act(cycle, np.zeros(1))
        decisions.append(bool(action[0]))
        policy.advance(cycle, action)

    assert decisions == [True, True, False, False]


def test_full_policy_harness_runs_all_arms_on_shared_endpoint_path(monkeypatch):
    received_f0_bounds = []
    received_scan_intervals = []
    original_loader = benchmark._load_source_module
    original_simulator = tc.simulate_two_endpoints

    def recording_simulator(scenarios, *args, **kwargs):
        received_scan_intervals.append(np.asarray(scenarios.get("scan_every")).copy())
        return original_simulator(scenarios, *args, **kwargs)

    monkeypatch.setattr(tc, "simulate_two_endpoints", recording_simulator)

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
        SignalModelParameters(q=0.6, a=0.2, b=0.4, d=0.1),
    )
    generation_parameters = SignalModelParameters(q=0.6, a=0.2, b=0.4, d=0.3)
    settings = benchmark.PolicyBenchmarkSettings(
        course_budget=6,
        tau_days=84.0,
        n_f0_nodes=8,
        mpc_particles=2,
        prediction_particles=1,
        forecast_cycles=3,
        q_episodes=1,
    )
    full_catalog = benchmark.baseline_policy_catalog(settings)
    test_catalog = {
        name: spec for name, spec in full_catalog.items()
        if not name.startswith(("filtered_threshold_grid_", "raw_threshold_grid_"))
    }
    filter_grid_rows = [
        (name, spec) for name, spec in full_catalog.items()
        if name.startswith("filtered_threshold_grid_")
    ][:2]
    raw_grid_rows = [
        (name, spec) for name, spec in full_catalog.items()
        if name.startswith("raw_threshold_grid_")
    ][:2]
    test_catalog.update(filter_grid_rows + raw_grid_rows)
    monkeypatch.setattr(benchmark, "baseline_policy_catalog", lambda _settings: test_catalog)

    result = benchmark.run_shared_policy_world(
        world_id="tiny-mechanics-test",
        selection_scenarios=selection,
        evaluation_scenarios=evaluation,
        evaluation_observations=observations,
        model_parameters=models,
        signal_generation_parameters=generation_parameters,
        signal_generation_seed=2026110603,
        signal_generation_mri_template_sha256=None,
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
    expected.update(
        name for name, spec in test_catalog.items()
        if spec["evaluate_all"] or name in {
            result["selected_baseline_by_family"]["filtered_threshold_grid"]["name"],
            result["selected_baseline_by_family"]["raw_threshold_grid"]["name"],
        }
    )
    assert set(result["arm_outcomes"]) == expected
    assert set(result["selected_baseline_by_family"]) == {
        "scan_rule", "state_estimation", "filtered_threshold_grid", "raw_threshold_grid"
    }
    assert result["baseline_selection"]["tie_rule"] == "first catalog entry wins exact selection-score ties"
    assert set(result["baseline_selection"]["scores"]) == set(test_catalog)
    assert sum(name.startswith("filtered_threshold_grid_") for name in result["arm_outcomes"]) == 1
    assert sum(name.startswith("raw_threshold_grid_") for name in result["arm_outcomes"]) == 1
    for family in ("filtered_threshold_grid", "raw_threshold_grid"):
        assert result["baseline_arm_metadata"][
            result["selected_baseline_by_family"][family]["name"]
        ]["evaluated_on_holdout"]
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
    assert result["signal_generation"] == {
        "generator": "trackc_integrated_signal_model.generate_integrated_signals",
        "parameters": {"q": 0.6, "a": 0.2, "b": 0.4, "d": 0.3},
        "seed": 2026110603,
        "f0_source": "evaluation_scenarios['f0']",
        "mri_template_sha256": None,
    }
    assert result["inference_model_parameters"] == [
        {"q": 0.6, "a": 0.2, "b": 0.4, "d": 0.1}
    ]
    assert result["signal_generation_in_inference_grid"] is False
    assert set(result["input_fingerprints"]) == {
        "selection_scenarios_sha256",
        "evaluation_scenarios_sha256",
        "evaluation_observations_sha256",
    }
    assert all(len(value) == 64 for value in result["input_fingerprints"].values())
    assert received_f0_bounds == [settings.signal_f0_bounds, settings.signal_f0_bounds]
    assert len(received_scan_intervals) == (
        settings.q_episodes + 1 + len(test_catalog) + len(result["arm_outcomes"])
    )
    assert all(np.all(interval == settings.scan_every) for interval in received_scan_intervals)


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
            signal_generation_parameters=SignalModelParameters(0.6, 0.2, 0.4, 0.3),
            signal_generation_seed=2026110600,
            signal_generation_mri_template_sha256=None,
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
            signal_generation_parameters=SignalModelParameters(0.6, 0.2, 0.4, 0.3),
            signal_generation_seed=2026110615,
            signal_generation_mri_template_sha256=None,
            selection_seed=[2026110616],
            q_training_seed=[2026110617],
            evaluation_seed=[2026110618],
        )


def test_input_fingerprint_is_key_order_independent_and_value_sensitive():
    left = {"b": np.asarray([2.0]), "a": np.asarray([1, 3], dtype=np.int64)}
    reversed_order = {"a": left["a"], "b": left["b"]}
    changed_value = {"a": left["a"], "b": np.asarray([2.5])}

    digest = benchmark._fingerprint_array_mapping(left)

    assert digest == benchmark._fingerprint_array_mapping(reversed_order)
    assert digest != benchmark._fingerprint_array_mapping(changed_value)


def test_signal_generation_metadata_validates_seed_and_optional_template_hash():
    parameters = SignalModelParameters(q=0.6, a=0.2, b=0.4, d=0.3)

    assert benchmark._validate_signal_generation_metadata(parameters, 42, None) == (42, None)
    assert benchmark._validate_signal_generation_metadata(parameters, 42, "A" * 64) == (
        42, "a" * 64
    )
    with pytest.raises(ValueError, match="one integer"):
        benchmark._validate_signal_generation_metadata(parameters, (42, 43), None)
    with pytest.raises(ValueError, match="64-character hex digest"):
        benchmark._validate_signal_generation_metadata(parameters, 42, "not-a-hash")


def test_harness_rejects_reused_signal_generation_seed_before_fitting():
    selection = tc.draw_scenarios(1, seed=2026110619, sigmas=(0.05,))
    evaluation = tc.draw_scenarios(1, seed=2026110620, sigmas=(0.05,))
    parameters = SignalModelParameters(0.6, 0.2, 0.4, 0.3)
    observations, _truth = generate_integrated_signals(
        1, q=0.6, a=0.2, b=0.4, d=0.3, seed=2026110621, f0=evaluation["f0"]
    )

    with pytest.raises(ValueError, match="signal-generation.*seeds must be distinct"):
        benchmark.run_shared_policy_world(
            world_id="reused-signal-seed",
            selection_scenarios=selection,
            evaluation_scenarios=evaluation,
            evaluation_observations=observations,
            model_parameters=(parameters,),
            signal_generation_parameters=parameters,
            signal_generation_seed=2026110622,
            signal_generation_mri_template_sha256=None,
            selection_seed=[2026110622],
            q_training_seed=[2026110623],
            evaluation_seed=[2026110624],
        )
