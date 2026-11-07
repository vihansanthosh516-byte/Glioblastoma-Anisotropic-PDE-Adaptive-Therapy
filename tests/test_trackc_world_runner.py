import pytest
import numpy as np

import trackc_seed_registry as registry
import trackc_world_schedule as schedule
from trackc_integrated_signal_model import SignalModelParameters
from trackc_world_runner import WorldRunError, run_development_worlds, prepare_synthetic_world


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
