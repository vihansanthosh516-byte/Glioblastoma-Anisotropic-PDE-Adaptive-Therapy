"""Config/code parity for integrated-v2 settings not covered by efc4ea6.

The sampler axes, q-signal levels, patient axes, sigma levels and the stop
endpoint code are already compared in tests/test_trackc_integrated_sampler.py
and tests/test_trackc_policy_benchmark.py. This file covers the remaining
numeric settings that are hard-coded in the simulator, controller, benchmark
and summary modules. It is a drift alarm, not a runtime loader: no module
loads the config, and a passing test does not make the config an execution
seal.
"""
import importlib.util
import inspect
import json
import math
import sys
from pathlib import Path

import pytest

import trackc_integrated_signal_model as signal_model
import trackc_sim as tc
import trackc_world_summary as world_summary

ROOT = Path(__file__).resolve().parent.parent
CONFIG = json.loads((ROOT / "configs" / "trackABC_integrated_v2.json").read_text(encoding="utf-8"))


def _load(filename):
    name = "parity_" + filename.replace(".py", "")
    spec = importlib.util.spec_from_file_location(name, ROOT / "src" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def benchmark():
    return _load("211_trackc_policy_benchmark.py")


@pytest.fixture(scope="module")
def mpc():
    return _load("208_trackc_robust_belief_mpc.py")


def _f0_low(text):
    assert text == "0.001 * exp(-0.7)"
    return 0.001 * math.exp(-0.7)


def test_simulation_clock_matches_simulator_constants():
    clock = CONFIG["simulation_clock_and_treatment"]
    assert clock["primary_step_days"] == tc.DT
    assert clock["cycle_days"] == tc.CYCLE_DAYS
    assert clock["course_treatment_days"] == tc.COURSE_DAYS
    assert clock["simulation_horizon_days"] == tc.HORIZON
    assert clock["max_courses"] == tc.BUDGET
    assert f">= {tc.RANO_FACTOR} *" in clock["stop_rule"]
    assert f"floor {tc.RANO_FLOOR:.2f} *" in clock["stop_rule"]


def test_scan_interval_matches_benchmark_default(benchmark):
    clock = CONFIG["simulation_clock_and_treatment"]
    settings = benchmark.PolicyBenchmarkSettings()
    assert clock["scan_interval_primary_days"] == settings.scan_every * tc.CYCLE_DAYS
    for days in clock["scan_interval_sensitivities_days"]:
        assert days % tc.CYCLE_DAYS == 0
        assert 1 <= days // tc.CYCLE_DAYS <= 4


@pytest.mark.xfail(strict=True, reason="known gap: the 112-day scan sensitivity (scan_every=4) is rejected with the configured 6 forecast cycles; its forecast length is undeclared")
def test_every_scan_sensitivity_runs_with_configured_forecast(benchmark):
    for days in CONFIG["simulation_clock_and_treatment"]["scan_interval_sensitivities_days"]:
        benchmark.PolicyBenchmarkSettings(scan_every=days // tc.CYCLE_DAYS)


def test_particle_filter_prior_matches_belief_ranges(mpc):
    prior = CONFIG["particle_filter_prior"]
    ranges = mpc.INTEGRATED_BELIEF_RANGES
    assert tuple(prior["r_day_minus1_log_uniform"]) == ranges["r"][:2]
    assert tuple(prior["m_log_uniform"]) == ranges["m"][:2]
    assert ranges["f0"][0] == pytest.approx(_f0_low(prior["f0_log_uniform"][0]), rel=1e-15)
    assert ranges["f0"][1] == prior["f0_log_uniform"][1]
    assert tuple(prior["c_uniform"]) == ranges["c"][:2]
    assert tuple(prior["epsilon_uniform"]) == ranges["eps"][:2]
    assert tuple(prior["B0_uniform"]) == ranges["B0"][:2]


def test_resistance_score_reference_matches_signal_model():
    reference = CONFIG["signal_model"]["resistance_score_reference"]
    assert signal_model.F0_REFERENCE_BOUNDS[0] == pytest.approx(_f0_low(reference["low"]), rel=1e-15)
    assert signal_model.F0_REFERENCE_BOUNDS[1] == reference["high"]
    assert reference["clip_probability"] == [1e-8, 1.0 - 1e-8]


def test_controller_settings_match_benchmark_defaults(benchmark):
    controller = CONFIG["controller_comparison"]
    settings = benchmark.PolicyBenchmarkSettings()
    assert controller["particle_filters"] == settings.mpc_particles
    assert controller["common_predictive_draws_per_model"] == settings.prediction_particles
    assert controller["forecast_cycles"] == settings.forecast_cycles
    assert controller["forecast_horizon_days"] == settings.forecast_cycles * tc.CYCLE_DAYS
    assert controller["action_block_cycles"] == settings.scan_every
    assert controller["candidate_patterns"] == 2 ** settings.scan_every


def test_outcome_settings_match_benchmark_and_utility_defaults(benchmark, mpc):
    outcomes = CONFIG["outcomes_and_inference"]
    settings = benchmark.PolicyBenchmarkSettings()
    assert outcomes["progression_horizon_days"] == settings.tau_days
    assert f"minus {settings.course_penalty:g} simulated days per course" in outcomes["utility"]
    utility = inspect.signature(mpc.robust_endpoint_utility).parameters
    assert utility["tau"].default == outcomes["progression_horizon_days"]
    assert utility["course_penalty"].default == settings.course_penalty
    assert settings.course_penalty not in outcomes["course_penalty_sensitivities_days"]


def test_precision_target_and_level_are_unchanged():
    """The 14-day target and 95% level are user-fixed; this pins both copies."""
    outcomes = CONFIG["outcomes_and_inference"]
    assert outcomes["precision_target_half_width_days"] == 14
    assert outcomes["confidence_level"] == 0.95
    assert "14 simulated days" in outcomes["interval_method_proposal"]["precision_rule"]
    defaults = inspect.signature(world_summary.simultaneous_max_t_world_intervals).parameters
    assert defaults["precision_target_half_width"].default == outcomes["precision_target_half_width_days"]
    assert defaults["confidence_level"].default == outcomes["confidence_level"]


def test_paired_difference_support_follows_from_utility_settings(benchmark):
    """[-790, 790] days is the support the bounded-mean intervals assume."""
    settings = benchmark.PolicyBenchmarkSettings()
    lowest = -settings.course_penalty * settings.course_budget
    highest = settings.tau_days
    assert (lowest, highest) == (-60.0, 730.0)
    assert highest - lowest == 790.0


def test_rejected_max_t_method_stays_marked_rejected():
    proposal = CONFIG["outcomes_and_inference"]["interval_method_proposal"]
    assert proposal["status"].startswith("rejected_by_predeclared_synthetic_coverage_screen")
    assert CONFIG["outcomes_and_inference"]["inference_status"].startswith("blocked_until")


@pytest.mark.xfail(strict=True, reason="known drift: controller model grid uses 3 q levels (81 models); config declares 5 (135)")
def test_controller_q_grid_matches_config_grid():
    grid = _load("212_trackc_signal_fisher_grid.py")
    assert list(grid.Q_GRID) == CONFIG["signal_model"]["synthetic_q_signal_grid"]


@pytest.mark.xfail(strict=True, reason="known drift: RobustBeliefMPC default stop endpoint is the legacy latent rule; the benchmark overrides it")
def test_direct_mpc_default_stop_endpoint_matches_config(mpc):
    default = inspect.signature(mpc.RobustBeliefMPC.__init__).parameters["stop_endpoint"].default
    assert default == CONFIG["simulation_clock_and_treatment"]["stop_endpoint_code"]
