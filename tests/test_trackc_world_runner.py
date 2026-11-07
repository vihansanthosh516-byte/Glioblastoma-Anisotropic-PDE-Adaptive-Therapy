import pytest
import numpy as np

import trackc_seed_registry as registry
import trackc_world_schedule as schedule
from trackc_integrated_signal_model import SignalModelParameters
from trackc_world_runner import (WorldRunError, run_development_worlds, prepare_synthetic_world,
                                 execute_prepared_world, make_synthetic_executor)


def world(name, start):
    values = {role: start + i for i, role in enumerate(schedule.SCENARIO_ROLES)}
    values.update({role: 1001 + start + i for i, role in enumerate(schedule.POLICY_ROLES)})
    return {"world_id": name, **values}


def fake_result(w):
    mapping = {
        "selection_scenario": "selection_scenario_generation", "calibration_scenario": "calibration_scenario_generation",
        "evaluation_scenario": "evaluation_scenario_generation", "calibration_signal": "calibration_signal_generation",
        "selection_policy": "selection", "q_training": "q_training", "calibration_policy": "calibration_policy_simulation",
        "evaluation_policy": "evaluation", "schedule_mixture": "schedule_mixture_master",
    }
    return {"world_id": w["world_id"], "analysis_type": "synthetic_virtual_world_development",
            "seeds": {key: [w[role]] for role, key in mapping.items()},
            "signal_generation": {"seed": w["evaluation_signal"]}}


def entries():
    return [registry._entry("reserved", "dev", "reserved", [(1, 1000)], [], purpose="stage2_dev")]


def test_runs_all_preflighted_worlds_in_order():
    worlds = [world("a", 1), world("b", 20)]
    assert [r["world_id"] for r in run_development_worlds(worlds, fake_result, entries=entries())] == ["a", "b"]


def test_preflight_failure_prevents_any_execution():
    worlds = [world("a", 1), world("b", 20)]
    worlds[1]["world"] = worlds[0]["world"]
    calls = []
    with pytest.raises(ValueError, match="reused"):
        run_development_worlds(worlds, lambda w: calls.append(w), entries=entries())
    assert not calls


def test_bad_record_stops_later_worlds():
    worlds = [world("a", 1), world("b", 20), world("c", 40)]
    calls = []

    def execute(w):
        calls.append(w["world_id"])
        result = fake_result(w)
        if w["world_id"] == "b":
            result["seeds"]["evaluation"] = [0]
        return result

    with pytest.raises(WorldRunError, match="world b failed after 1") as error:
        run_development_worlds(worlds, execute, entries=entries())
    assert isinstance(error.value.__cause__, ValueError)
    assert calls == ["a", "b"]


def test_prepares_reproducible_disjoint_synthetic_inputs_without_scoring():
    w = world("a", 1)
    options = dict(family="single", patients_per_cohort=2,
                   signal_parameters=SignalModelParameters(0.7, 0.2, 0.3, 0.1),
                   model_parameters=(SignalModelParameters(0.7, 0.2, 0.3, 0.1),),
                   settings="synthetic-settings", entries=entries())
    first = prepare_synthetic_world(w, **options)
    second = prepare_synthetic_world(w, **options)
    assert first["world_generation_seed"] == w["world"]
    assert len(first["world_parameters_sha256"]) == 64
    assert first["world_parameters_sha256"] == second["world_parameters_sha256"]
    a, b = first["benchmark_kwargs"], second["benchmark_kwargs"]
    np.testing.assert_array_equal(a["evaluation_scenarios"]["f0"], b["evaluation_scenarios"]["f0"])
    assert not set(a["selection_patient_ids"]) & set(a["evaluation_patient_ids"])
    assert a["signal_generation_seed"] == w["evaluation_signal"]
    assert a["calibration_signal_generation_seed"] == w["calibration_signal"]
    assert a["include_ambiguity_controller"] is True


def test_world_fingerprint_changes_when_world_parameters_change(monkeypatch):
    import trackc_world_runner as runner
    original = runner.calibration.draw_worlds

    def changed(*args, **kwargs):
        generated = original(*args, **kwargs)
        generated[0]["c"] += 0.001
        return generated

    options = dict(family="single", patients_per_cohort=2,
                   signal_parameters=SignalModelParameters(0.7, 0.2, 0.3, 0.1),
                   model_parameters=(SignalModelParameters(0.7, 0.2, 0.3, 0.1),),
                   settings="synthetic-settings", entries=entries())
    baseline = prepare_synthetic_world(world("a", 1), **options)
    monkeypatch.setattr(runner.calibration, "draw_worlds", changed)
    modified = prepare_synthetic_world(world("a", 1), **options)
    assert baseline["world_parameters_sha256"] != modified["world_parameters_sha256"]


def test_prepared_provenance_is_attached_to_matching_fake_benchmark():
    prepared = {"world_generation_seed": 7, "family": "single", "scan_loss": None, "world_parameters_sha256": "a" * 64,
                "benchmark_kwargs": {"world_id": "sample"}}
    result = execute_prepared_world(prepared, lambda **kwargs: {"world_id": kwargs["world_id"]})
    assert result["world_provenance"]["world_parameters_sha256"] == "a" * 64
    assert result["world_provenance"]["world_generation_seed"] == 7
    with pytest.raises(ValueError, match="world ID"):
        execute_prepared_world(prepared, lambda **kwargs: {"world_id": "different"})
    with pytest.raises(ValueError, match="already contains"):
        execute_prepared_world(prepared, lambda **kwargs: {"world_id": "sample", "world_provenance": {}})


def test_composed_runner_builds_and_checks_two_synthetic_worlds():
    def benchmark(**kwargs):
        return {"world_id": kwargs["world_id"], "analysis_type": "synthetic_virtual_world_development",
                "seeds": {
                    "selection_scenario_generation": [kwargs["selection_scenario_seed"]],
                    "calibration_scenario_generation": [kwargs["calibration_scenario_seed"]],
                    "evaluation_scenario_generation": [kwargs["evaluation_scenario_seed"]],
                    "calibration_signal_generation": [kwargs["calibration_signal_generation_seed"]],
                    "selection": kwargs["selection_seed"], "q_training": kwargs["q_training_seed"],
                    "calibration_policy_simulation": kwargs["calibration_seed"],
                    "evaluation": kwargs["evaluation_seed"], "schedule_mixture_master": kwargs["schedule_mixture_seed"],
                }, "signal_generation": {"seed": kwargs["signal_generation_seed"]}}

    params = SignalModelParameters(0.7, 0.2, 0.3, 0.1)
    executor = make_synthetic_executor(family="single", patients_per_cohort=2,
        signal_parameters=params, model_parameters=(params,), settings="test", benchmark=benchmark, entries=entries())
    results = run_development_worlds([world("a", 1), world("b", 20)], executor, entries=entries())
    assert [row["world_provenance"]["world_generation_seed"] for row in results] == [1, 20]
    assert results[0]["world_provenance"]["world_parameters_sha256"] != results[1]["world_provenance"]["world_parameters_sha256"]


def test_scan_loss_is_visible_and_applied_to_every_cohort():
    params = SignalModelParameters(0.7, 0.2, 0.3, 0.1)
    options = dict(family="single", patients_per_cohort=2, signal_parameters=params,
                   model_parameters=(params,), settings="test", entries=entries())
    prepared = prepare_synthetic_world(world("a", 1), scan_loss={"miss_p": 0.2, "miss_p_hi": 0.7}, **options)
    assert prepared["scan_loss"] == {"miss_p": 0.2, "miss_p_hi": 0.7}
    for role in ("selection", "calibration", "evaluation"):
        scenario = prepared["benchmark_kwargs"][f"{role}_scenarios"]
        np.testing.assert_array_equal(scenario["miss_p"], [0.2, 0.2])
        np.testing.assert_array_equal(scenario["miss_p_hi"], [0.7, 0.7])
        assert scenario["miss_visible"].all()
    with pytest.raises(ValueError, match="finite probability"):
        prepare_synthetic_world(world("a", 1), scan_loss={"miss_p": 1.1}, **options)
