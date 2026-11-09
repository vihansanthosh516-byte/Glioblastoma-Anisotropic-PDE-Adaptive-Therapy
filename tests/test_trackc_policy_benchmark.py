"""Contract tests for the paired synthetic policy benchmark."""
from __future__ import annotations

import importlib.util
import inspect
from copy import deepcopy
from pathlib import Path
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

import trackc_sim as tc  # noqa: E402
from trackc_world_summary import summarize_benchmark_worlds  # noqa: E402
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


def test_ambiguity_controller_is_opt_in_by_default():
    parameter = inspect.signature(benchmark.run_shared_policy_world).parameters[
        "include_ambiguity_controller"
    ]
    assert parameter.default is False


def test_mpc_assumed_future_scan_loss_requires_observed_stop():
    with pytest.raises(ValueError, match="requires observed_two_scan"):
        benchmark.PolicyBenchmarkSettings(mpc_assumed_future_scan_miss_probability=0.1)
    settings = benchmark.PolicyBenchmarkSettings(
        stop_endpoint="observed_two_scan",
        mpc_assumed_future_scan_miss_probability=0.25,
    )
    assert settings.mpc_assumed_future_scan_miss_probability == pytest.approx(0.25)


def test_scan_loss_scenarios_reject_invalid_or_misleading_inputs():
    scenario = tc.draw_scenarios(2, seed=1234)
    for bad in ([-0.1, 0.2], [0.1, 1.1], [0.1, np.nan], [[0.1], [0.2]]):
        with pytest.raises(ValueError, match="miss_p"):
            benchmark._validate_scenarios(dict(scenario, miss_p=np.asarray(bad)), "test")
    with pytest.raises(ValueError, match="requires miss_p"):
        benchmark._validate_scenarios(dict(scenario, miss_p_hi=np.zeros(2)), "test")
    with pytest.raises(ValueError, match="requires visible missing scans"):
        benchmark._validate_scenarios(dict(scenario, miss_p=np.array([0.0, 0.2])), "test")
    with pytest.raises(ValueError, match="requires visible missing scans"):
        benchmark._validate_scenarios(
            dict(scenario, miss_p=np.zeros(2), miss_p_hi=np.array([0.0, 0.2])), "test"
        )
    assert benchmark._validate_scenarios(dict(scenario, miss_p=np.zeros(2)), "test") == 2
    for visibility in ([False, False], [True, False], [1, 1]):
        with pytest.raises(ValueError, match="miss_visible"):
            benchmark._validate_scenarios(
                dict(scenario, miss_p=np.zeros(2), miss_visible=np.asarray(visibility)), "test"
            )
    assert benchmark._validate_scenarios(
        dict(scenario, miss_p=np.array([0.0, 1.0]), miss_visible=np.ones(2, dtype=bool)), "test"
    ) == 2


def test_benchmark_blocks_legacy_repeated_reading_from_unseen_cycle():
    class RecordPolicy:
        uses_filter = False

        def reset(self, scen, rng):
            self.readings = []

        def act(self, cycle, y):
            self.readings.append(y.copy())
            return np.zeros(1, dtype=bool)

        def advance(self, cycle, treat):
            pass

    scenario = tc.draw_scenarios(1, seed=1234)
    scenario.update(scan_every=np.array([3]), sigma_true=np.zeros(1), miss_p=np.array([0.5]))
    seed = next(
        value for value in range(1000)
        if (lambda u: u[0, 1] >= 0.5 and u[0, 2] >= 0.5 and u[0, 3] < 0.5)(
            np.random.default_rng([value, 1]).random((1, tc.N_CYCLES))
        )
    )
    policy = RecordPolicy()
    tc.simulate_two_endpoints(scenario, policy, seed, stop_endpoint="observed_two_scan")
    assert np.isnan(policy.readings[1][0]) and np.isnan(policy.readings[2][0])
    assert np.isfinite(policy.readings[3][0])
    assert policy.readings[3][0] != policy.readings[0][0]
    with pytest.raises(ValueError, match="requires visible missing scans"):
        benchmark._validate_scenarios(scenario, "test")


def test_policy_rng_consumption_does_not_change_shared_scan_noise():
    class RecordPolicy:
        uses_filter = False

        def __init__(self, reset_draws):
            self.reset_draws = reset_draws

        def reset(self, scen, rng):
            del scen
            self.readings = []
            rng.random(self.reset_draws)

        def act(self, cycle, y):
            self.readings.append(y.copy())
            return np.zeros(1, dtype=bool)

        def advance(self, cycle, treat):
            pass

    scenario = tc.draw_scenarios(1, seed=2026110901, sigmas=(0.1,))
    settings = benchmark.PolicyBenchmarkSettings(
        stop_endpoint="observed_two_scan", scan_every=3,
    )
    first = RecordPolicy(reset_draws=0)
    second = RecordPolicy(reset_draws=1000)

    benchmark._evaluate_policy(first, scenario, (2026110902,), settings)
    benchmark._evaluate_policy(second, scenario, (2026110902,), settings)

    assert len(first.readings) == len(second.readings) == tc.N_CYCLES
    assert all(np.array_equal(left, right, equal_nan=True)
               for left, right in zip(first.readings, second.readings))


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


def _calibration_cohort(n, scenario_seed, signal_seed):
    scenarios = tc.draw_scenarios(n, seed=scenario_seed, sigmas=(0.05,))
    observations, _truth = generate_integrated_signals(
        n, q=0.6, a=0.2, b=0.4, d=0.3,
        seed=signal_seed, f0=scenarios["f0"],
    )
    return scenarios, observations


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


def test_robust_utility_takes_each_patients_worse_endpoint_before_averaging():
    outcome = {
        "ttp_baseline": np.asarray([100.0, 700.0]),
        "capped_baseline": np.asarray([False, False]),
        "ttp_rano": np.asarray([700.0, 100.0]),
        "capped_rano": np.asarray([False, False]),
        "courses": np.asarray([0, 0]),
        "courses_through_cutoff": np.asarray([0, 0]),
        "stop_endpoint": "either",
    }

    utility = benchmark.robust_utility(
        outcome, tau_days=730.0, course_penalty=5.0
    )

    # The estimand is mean_i min(T_baseline_i, T_RANO_i, tau), which is 100.
    # Taking min(mean(T_baseline), mean(T_RANO)) would instead give 400.
    np.testing.assert_allclose(utility, [100.0, 100.0])
    assert float(utility.mean()) == 100.0


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


def test_fixed_schedule_controls_expose_six_and_twelve_course_references():
    full_budget = benchmark.PolicyBenchmarkSettings(course_budget=12)
    controls = benchmark._fixed_schedule_control_catalog(full_budget)

    assert set(controls) == {
        "standard_six_course",
        "standard_twelve_course_front_loaded",
        "paced_twelve_course_every_other_cycle",
    }
    assert controls["standard_six_course"]["cycle_indices"] == list(range(6))
    assert controls["standard_twelve_course_front_loaded"]["cycle_indices"] == list(range(12))
    assert controls["paced_twelve_course_every_other_cycle"]["cycle_indices"] == list(range(0, 24, 2))
    assert controls["paced_twelve_course_every_other_cycle"]["nominal_interval_days"] == 56
    assert all(spec["schedule"].sum() == spec["planned_courses"] for spec in controls.values())

    reduced_budget = benchmark.PolicyBenchmarkSettings(course_budget=6)
    reduced_controls = benchmark._fixed_schedule_control_catalog(reduced_budget)
    assert set(reduced_controls) == {"standard_six_course"}


def test_full_policy_harness_runs_all_arms_on_shared_endpoint_path(monkeypatch):
    received_f0_bounds = []
    received_scan_intervals = []
    received_stop_endpoints = []
    received_future_scan_miss_probabilities = []
    original_loader = benchmark._load_source_module
    original_simulator = tc.simulate_two_endpoints

    def recording_simulator(scenarios, *args, **kwargs):
        received_scan_intervals.append(np.asarray(scenarios.get("scan_every")).copy())
        received_stop_endpoints.append(kwargs.get("stop_endpoint"))
        return original_simulator(scenarios, *args, **kwargs)

    monkeypatch.setattr(tc, "simulate_two_endpoints", recording_simulator)

    def recording_loader(module_name, filename):
        module = original_loader(module_name, filename)
        if filename == "214_trackc_loo_model_average.py":
            original_controller = module.LOOModelAveragedBeliefMPC

            def recording_controller(*args, **kwargs):
                received_f0_bounds.append(tuple(kwargs["f0_bounds"]))
                received_future_scan_miss_probabilities.append(
                    kwargs["assumed_future_scan_miss_probability"]
                )
                return original_controller(*args, **kwargs)

            monkeypatch.setattr(module, "LOOModelAveragedBeliefMPC", recording_controller)
        return module

    monkeypatch.setattr(benchmark, "_load_source_module", recording_loader)
    selection = tc.draw_scenarios(2, seed=2026110601, sigmas=(0.05,))
    calibration, calibration_observations = _calibration_cohort(
        2, scenario_seed=2026110611, signal_seed=2026110618
    )
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
        SignalModelParameters(q=0.8, a=-0.1, b=0.7, d=-0.2),
    )
    generation_parameters = SignalModelParameters(q=0.6, a=0.2, b=0.4, d=0.3)
    settings = benchmark.PolicyBenchmarkSettings(
        course_budget=12,
        stop_endpoint="observed_two_scan",
        mpc_assumed_future_scan_miss_probability=0.25,
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
        calibration_scenarios=calibration,
        evaluation_scenarios=evaluation,
        selection_patient_ids=("selection-0", "selection-1"),
        calibration_patient_ids=("calibration-0", "calibration-1"),
        evaluation_patient_ids=("evaluation-0", "evaluation-1"),
        calibration_observations=calibration_observations,
        evaluation_observations=observations,
        model_parameters=models,
        signal_generation_parameters=generation_parameters,
        signal_generation_seed=2026110603,
        calibration_signal_generation_seed=2026110618,
        signal_generation_mri_template_sha256=None,
        selection_scenario_seed=2026110610,
        calibration_scenario_seed=2026110611,
        evaluation_scenario_seed=2026110612,
        selection_seed=[2026110604],
        q_training_seed=[2026110605],
        calibration_seed=[2026110619],
        evaluation_seed=[2026110606],
        schedule_mixture_seed=[2026110620],
        settings=settings,
        include_ambiguity_controller=True,
        primary_contrast="minimax",
    )

    # Exercise the adapter with the actual result schema. This copied record is
    # synthetic test data, not evidence of two independent simulated worlds.
    result["world_provenance"] = {
        "family": "single",
        "scan_loss": None,
        "world_generation_seed": 900001,
        "world_parameters_sha256": "c" * 64,
    }
    copied = deepcopy(result)
    copied["world_id"] = "copied-schema-test"
    copied["world_provenance"] = {
        "family": "single",
        "scan_loss": None,
        "world_generation_seed": 900002,
        "world_parameters_sha256": "d" * 64,
    }
    for index, key in enumerate(("selection_scenarios_sha256", "calibration_scenarios_sha256", "evaluation_scenarios_sha256")):
        copied["input_fingerprints"][key] = format(100 + index, "064x")
    copied["seeds"] = {key: [int(x) + 1000000 for x in value] for key, value in copied["seeds"].items()}
    copied["signal_generation"]["seed"] += 1000000
    summary = summarize_benchmark_worlds([result, copied], resamples=100, seed=42)
    assert summary["n_worlds"] == 2
    assert summary["primary_contrast_name"] == "minimax"
    assert summary["equal_world_mean_difference"] == pytest.approx(
        result["paired_contrasts"]["minimax"]["robust_utility_difference"].mean()
    )

    expected = {
        "loo_dependence_aware_mpc",
        "loo_conditional_independence_mpc",
        "no_signal_mpc",
        "q_learning",
        "selection_best_fixed",
        "standard_six_course",
        "standard_twelve_course_front_loaded",
        "paced_twelve_course_every_other_cycle",
        "finite_set_minimax_mpc",
        "finite_set_d0_minimax_mpc",
    }
    expected.update(
        name for name, spec in test_catalog.items()
        if spec["evaluate_all"] or name in {
            result["selected_baseline_by_family"]["filtered_threshold_grid"]["name"],
            result["selected_baseline_by_family"]["raw_threshold_grid"]["name"],
        }
    )
    matched_arm_names = {
        row["evaluation_arm"]
        for row in result["matched_fixed_schedule_controls"].values()
        if row["status"] == "evaluated"
    }
    assert set(result["arm_outcomes"]) == expected | matched_arm_names
    assert received_stop_endpoints and set(received_stop_endpoints) == {"observed_two_scan"}
    assert received_future_scan_miss_probabilities
    assert set(received_future_scan_miss_probabilities) == {0.25}
    assert set(result["calibration_course_use"]) == set(result["matched_fixed_schedule_controls"])
    for policy_name, control in result["matched_fixed_schedule_controls"].items():
        assert control["calibration_mean_courses"] == pytest.approx(
            result["calibration_course_use"][policy_name]
        )
        assert control["selection_mix"]["target_mean_courses"] == pytest.approx(
            control["calibration_mean_courses"]
        )
        if control["status"] == "evaluated":
            assert control["evaluation_arm"] in result["arm_outcomes"]
            assert control["evaluation_course_use_difference"] == pytest.approx(
                control["evaluation_mean_courses"] - control["calibration_mean_courses"]
            )
            assert control["paired_worst_endpoint_days_difference"].shape == (2,)
            assert control["paired_robust_utility_difference"].shape == (2,)
            assert len(control["assignment"]["assignment_sha256"]) == 64
    assert set(result["selected_baseline_by_family"]) == {
        "scan_rule", "state_estimation", "filtered_threshold_grid", "raw_threshold_grid"
    }
    assert result["baseline_selection"]["tie_rule"] == "first catalog entry wins exact selection-score ties"
    fixed_controls = result["fixed_schedule_control_catalog"]
    assert fixed_controls["standard_twelve_course_front_loaded"]["cycle_indices"] == list(range(12))
    assert fixed_controls["paced_twelve_course_every_other_cycle"]["nominal_interval_days"] == 56
    assert set(result["baseline_selection"]["scores"]) == set(test_catalog)
    assert sum(name.startswith("filtered_threshold_grid_") for name in result["arm_outcomes"]) == 1
    assert sum(name.startswith("raw_threshold_grid_") for name in result["arm_outcomes"]) == 1
    for family in ("filtered_threshold_grid", "raw_threshold_grid"):
        assert result["baseline_arm_metadata"][
            result["selected_baseline_by_family"][family]["name"]
        ]["evaluated_on_holdout"]
    assert result["primary_contrast_name"] == "minimax"
    assert result["primary_paired_contrast"]["robust_utility_difference"].shape == (2,)
    assert result["primary_paired_contrast"]["contrast"] == (
        "finite_set_minimax_mpc minus finite_set_d0_minimax_mpc"
    )
    assert result["paired_contrasts"]["loo"]["contrast"] == (
        "loo_dependence_aware_mpc minus loo_conditional_independence_mpc"
    )
    assert result["paired_contrasts"]["minimax"]["robust_utility_difference"].shape == (2,)
    assert result["exploratory_controller"]["included"] is True
    assert result["exploratory_controller"]["model_set_size"] == 2
    assert result["exploratory_controller"]["role"] == "primary synthetic development candidate"
    assert result["exploratory_d0_controller"]["included"] is True
    assert result["exploratory_d0_controller"]["model_set_size"] == len(
        benchmark.conditional_independence_grid(models)
    )
    assert result["exploratory_d0_controller"]["dependence_rule"] == (
        "same q/a/b support as the full set, with d fixed to 0"
    )
    assert result["exploratory_d0_controller"]["role"] == "matched primary control"
    assert set(result["arm_outcomes"]["finite_set_minimax_mpc"]) >= {
        "ttp_baseline", "ttp_rano", "courses", "robust_utility"
    }
    assert set(result["arm_outcomes"]["finite_set_d0_minimax_mpc"]) >= {
        "ttp_baseline", "ttp_rano", "courses", "robust_utility"
    }
    for arm in result["arm_outcomes"].values():
        assert arm["ttp_baseline"].shape == (2,)
        assert arm["ttp_rano"].shape == (2,)
        assert arm["courses"].shape == (2,)
        assert arm["courses_through_cutoff"].shape == (2,)
        assert arm["capped_baseline"].dtype == np.bool_
        assert arm["capped_rano"].dtype == np.bool_
        assert (arm["courses_through_cutoff"] <= arm["courses"]).all()
        assert (arm["courses"] <= settings.course_budget).all()
    assert result["analysis_type"] == "synthetic_virtual_world_development"
    assert result["settings"]["stop_endpoint"] == "observed_two_scan"
    assert result["settings"]["mpc_assumed_future_scan_miss_probability"] == pytest.approx(0.25)
    assert result["settings"]["signal_f0_bounds"] == tuple(map(float, tc.RANGES["f0"][:2]))
    assert result["signal_generation"] == {
        "generator": "trackc_integrated_signal_model.generate_integrated_signals",
        "parameters": {"q": 0.6, "a": 0.2, "b": 0.4, "d": 0.3},
        "seed": 2026110603,
        "f0_source": "evaluation_scenarios['f0']",
        "mri_template_sha256": None,
    }
    assert result["inference_model_parameters"] == [
        {"q": 0.6, "a": 0.2, "b": 0.4, "d": 0.1},
        {"q": 0.8, "a": -0.1, "b": 0.7, "d": -0.2},
    ]
    assert result["signal_generation_in_inference_grid"] is False
    assert set(result["input_fingerprints"]) == {
        "selection_scenarios_sha256",
        "calibration_scenarios_sha256",
        "evaluation_scenarios_sha256",
        "selection_patient_ids_sha256",
        "calibration_patient_ids_sha256",
        "evaluation_patient_ids_sha256",
        "calibration_observations_sha256",
        "evaluation_observations_sha256",
    }
    assert all(len(value) == 64 for value in result["input_fingerprints"].values())
    assert result["cohorts"]["selection"]["n_patients"] == 2
    assert result["cohorts"]["calibration"]["n_patients"] == 2
    assert result["cohorts"]["evaluation"]["n_patients"] == 2
    assert result["seeds"]["selection_scenario_generation"] == [2026110610]
    assert result["seeds"]["calibration_scenario_generation"] == [2026110611]
    assert result["seeds"]["evaluation_scenario_generation"] == [2026110612]
    assert result["scenario_seed_provenance"].startswith("caller-declared;")
    assert "each arm has its own action-dependent path" in result["path_rule"]
    assert "share one treatment path" not in result["path_rule"]
    assert received_f0_bounds == [settings.signal_f0_bounds] * 4
    assert len(received_scan_intervals) == (
        settings.q_episodes + 1 + len(test_catalog)
        + len(result["calibration_course_use"]) + len(result["arm_outcomes"])
    )
    assert all(np.all(interval == settings.scan_every) for interval in received_scan_intervals)


def test_minimax_primary_requires_including_finite_set_controllers():
    required = {
        name: None
        for name, parameter in inspect.signature(benchmark.run_shared_policy_world).parameters.items()
        if parameter.default is inspect.Parameter.empty
    }
    required.update({
        "world_id": "missing-minimax-arms",
        "include_ambiguity_controller": False,
        "primary_contrast": "minimax",
    })
    with pytest.raises(ValueError, match="requires the finite-set controllers"):
        benchmark.run_shared_policy_world(**required)


def test_harness_rejects_overlapping_synthetic_patient_ids_before_policy_fit(monkeypatch):
    selection = tc.draw_scenarios(1, seed=2026110631, sigmas=(0.05,))
    calibration, calibration_observations = _calibration_cohort(
        1, scenario_seed=2026110637, signal_seed=2026110638
    )
    evaluation = tc.draw_scenarios(1, seed=2026110632, sigmas=(0.05,))
    parameters = SignalModelParameters(0.6, 0.2, 0.4, 0.3)
    observations, _truth = generate_integrated_signals(
        1, q=0.6, a=0.2, b=0.4, d=0.3, seed=2026110633, f0=evaluation["f0"]
    )

    def unexpected_fit(*_args, **_kwargs):
        pytest.fail("policy fitting must not start for overlapping cohorts")

    monkeypatch.setattr(benchmark, "select_fixed_schedule", unexpected_fit)
    with pytest.raises(ValueError, match="pairwise disjoint"):
        benchmark.run_shared_policy_world(
            world_id="overlapping-cohorts",
            selection_scenarios=selection,
            calibration_scenarios=calibration,
            evaluation_scenarios=evaluation,
            selection_patient_ids=("same-virtual-patient",),
            calibration_patient_ids=("calibration-0",),
            evaluation_patient_ids=("same-virtual-patient",),
            calibration_observations=calibration_observations,
            evaluation_observations=observations,
            model_parameters=(parameters,),
            signal_generation_parameters=parameters,
            signal_generation_seed=2026110633,
            calibration_signal_generation_seed=2026110638,
            signal_generation_mri_template_sha256=None,
            selection_scenario_seed=2026110631,
            calibration_scenario_seed=2026110637,
            evaluation_scenario_seed=2026110632,
            selection_seed=[2026110634],
            q_training_seed=[2026110635],
            calibration_seed=[2026110639],
            evaluation_seed=[2026110636],
            schedule_mixture_seed=[2026110640],
        )


def test_harness_requires_budget_that_can_include_six_course_comparator():
    selection = tc.draw_scenarios(1, seed=2026110607, sigmas=(0.05,))
    calibration, calibration_observations = _calibration_cohort(
        1, scenario_seed=2026110613, signal_seed=2026110614
    )
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
            calibration_scenarios=calibration,
            evaluation_scenarios=evaluation,
            selection_patient_ids=("selection-0",),
            calibration_patient_ids=("calibration-0",),
            evaluation_patient_ids=("evaluation-0",),
            calibration_observations=calibration_observations,
            evaluation_observations=observations,
            model_parameters=(SignalModelParameters(0.6, 0.2, 0.4, 0.3),),
            signal_generation_parameters=SignalModelParameters(0.6, 0.2, 0.4, 0.3),
            signal_generation_seed=2026110609,
            calibration_signal_generation_seed=2026110614,
            signal_generation_mri_template_sha256=None,
            selection_scenario_seed=2026110607,
            calibration_scenario_seed=2026110613,
            evaluation_scenario_seed=2026110608,
            selection_seed=[2026110610],
            q_training_seed=[2026110611],
            calibration_seed=[2026110615],
            evaluation_seed=[2026110612],
            schedule_mixture_seed=[2026110616],
            settings=settings,
        )


def test_harness_rejects_f0_values_outside_declared_signal_prior_before_fitting():
    selection = tc.draw_scenarios(2, seed=2026110613, sigmas=(0.05,))
    calibration, calibration_observations = _calibration_cohort(
        2, scenario_seed=2026110619, signal_seed=2026110620
    )
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
            calibration_scenarios=calibration,
            evaluation_scenarios=evaluation,
            selection_patient_ids=("selection-0", "selection-1"),
            calibration_patient_ids=("calibration-0", "calibration-1"),
            evaluation_patient_ids=("evaluation-0", "evaluation-1"),
            calibration_observations=calibration_observations,
            evaluation_observations=observations,
            model_parameters=(SignalModelParameters(0.6, 0.2, 0.4, 0.3),),
            signal_generation_parameters=SignalModelParameters(0.6, 0.2, 0.4, 0.3),
            signal_generation_seed=2026110615,
            calibration_signal_generation_seed=2026110620,
            signal_generation_mri_template_sha256=None,
            selection_scenario_seed=2026110613,
            calibration_scenario_seed=2026110619,
            evaluation_scenario_seed=2026110614,
            selection_seed=[2026110616],
            q_training_seed=[2026110617],
            calibration_seed=[2026110621],
            evaluation_seed=[2026110618],
            schedule_mixture_seed=[2026110622],
        )


def test_input_fingerprint_is_key_order_independent_and_value_sensitive():
    left = {"b": np.asarray([2.0]), "a": np.asarray([1, 3], dtype=np.int64)}
    reversed_order = {"a": left["a"], "b": left["b"]}
    changed_value = {"a": left["a"], "b": np.asarray([2.5])}

    digest = benchmark._fingerprint_array_mapping(left)

    assert digest == benchmark._fingerprint_array_mapping(reversed_order)
    assert digest != benchmark._fingerprint_array_mapping(changed_value)


def test_patient_ids_must_be_unique_and_canonical():
    with pytest.raises(ValueError, match="unique IDs"):
        benchmark._validate_patient_ids(("virtual-1", "virtual-1"), 2, "ids")
    with pytest.raises(ValueError, match="leading or trailing whitespace"):
        benchmark._validate_patient_ids(("virtual-1 ",), 1, "ids")


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
    calibration, calibration_observations = _calibration_cohort(
        1, scenario_seed=2026110625, signal_seed=2026110628
    )
    parameters = SignalModelParameters(0.6, 0.2, 0.4, 0.3)
    observations, _truth = generate_integrated_signals(
        1, q=0.6, a=0.2, b=0.4, d=0.3, seed=2026110622, f0=evaluation["f0"]
    )

    with pytest.raises(ValueError, match="seeds must be distinct"):
        benchmark.run_shared_policy_world(
            world_id="reused-signal-seed",
            selection_scenarios=selection,
            calibration_scenarios=calibration,
            evaluation_scenarios=evaluation,
            selection_patient_ids=("selection-0",),
            calibration_patient_ids=("calibration-0",),
            evaluation_patient_ids=("evaluation-0",),
            calibration_observations=calibration_observations,
            evaluation_observations=observations,
            model_parameters=(parameters,),
            signal_generation_parameters=parameters,
            signal_generation_seed=2026110622,
            calibration_signal_generation_seed=2026110628,
            signal_generation_mri_template_sha256=None,
            selection_scenario_seed=2026110619,
            calibration_scenario_seed=2026110625,
            evaluation_scenario_seed=2026110620,
            selection_seed=[2026110622],
            q_training_seed=[2026110623],
            calibration_seed=[2026110626],
            evaluation_seed=[2026110624],
            schedule_mixture_seed=[2026110627],
        )
