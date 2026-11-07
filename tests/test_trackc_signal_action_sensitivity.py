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
SCRIPT_PATH = ROOT / "src" / "220_trackc_signal_action_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("trackc_signal_action_sensitivity_tested", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise ImportError(f"cannot load diagnostic module: {SCRIPT_PATH}")
DIAGNOSTIC = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = DIAGNOSTIC
SPEC.loader.exec_module(DIAGNOSTIC)


def _prepared_manifest():
    current = DIAGNOSTIC.build_manifest_unchecked()
    return {
        "record_type": "pre_run_synthetic_action_sensitivity_manifest",
        "protocol_status": "prepared",
        "created_utc": "2026-11-06T00:00:00+00:00",
        "script_id": DIAGNOSTIC.SCRIPT_ID,
        "project_title": DIAGNOSTIC.PROJECT_TITLE,
        **current,
        "seed_allocation_check": {
            "purpose": "stage2_dev",
            "scenario_range_free_at_prepare": True,
            "policy_range_free_at_prepare": True,
        },
        "interpretation_rule": DIAGNOSTIC.INTERPRETATION_RULE,
        "outcome_metrics_included": False,
    }


def _rows():
    return [
        {"state_id": "s1", "signal_replicate": 0, "model_set": "full_81", "pattern": "100"},
        {"state_id": "s1", "signal_replicate": 0, "model_set": "matched_d0_27", "pattern": "100"},
        {"state_id": "s1", "signal_replicate": 1, "model_set": "full_81", "pattern": "110"},
        {"state_id": "s1", "signal_replicate": 1, "model_set": "matched_d0_27", "pattern": "100"},
        {"state_id": "s2", "signal_replicate": 0, "model_set": "full_81", "pattern": "001"},
        {"state_id": "s2", "signal_replicate": 0, "model_set": "matched_d0_27", "pattern": "001"},
    ]


def test_summary_reports_signal_plan_diversity_and_paired_model_agreement():
    result = DIAGNOSTIC.summarize_action_rows(_rows())

    assert result["n_rows"] == 6
    assert result["n_paired_signal_cases"] == 3
    assert result["full_vs_d0_agreement"] == {"agreements": 2, "n_pairs": 3, "fraction": pytest.approx(2 / 3)}
    state = {row["state_id"]: row for row in result["by_state"]}["s1"]
    assert state["n_signal_replicates"] == 2
    assert state["pattern_counts"]["full_81"] == {"100": 1, "110": 1}
    assert state["pattern_counts"]["matched_d0_27"] == {"100": 2}


def test_summary_refuses_duplicate_model_state_signal_rows():
    rows = _rows()
    rows.append(dict(rows[0]))

    with pytest.raises(ValueError, match="duplicate"):
        DIAGNOSTIC.summarize_action_rows(rows)


def test_summary_refuses_unpaired_model_sets():
    rows = [row for row in _rows() if not (row["state_id"] == "s1" and row["signal_replicate"] == 1 and row["model_set"] == "matched_d0_27")]

    with pytest.raises(ValueError, match="paired"):
        DIAGNOSTIC.summarize_action_rows(rows)


def test_summary_refuses_nonbinary_action_patterns():
    rows = _rows()
    rows[0]["pattern"] = "10x"

    with pytest.raises(ValueError, match="binary"):
        DIAGNOSTIC.summarize_action_rows(rows)


def test_json_safe_converts_numpy_values_and_nonfinite_numbers():
    converted = DIAGNOSTIC.json_safe({"array": np.array([1.0, np.nan]), "scalar": np.float64(2.5)})

    assert converted == {"array": [1.0, None], "scalar": 2.5}


def test_environment_record_includes_python_numpy_and_scipy_versions():
    environment = DIAGNOSTIC.environment_record()

    assert environment["python"]
    assert environment["numpy"] == np.__version__
    assert environment["scipy"]
    assert environment["platform"]


def test_design_uses_unique_frozen_development_seeds_and_expected_size():
    design = DIAGNOSTIC.build_design()

    assert len(design) == 64
    assert {(row["state_index"], row["signal_replicate"]) for row in design} == {
        (state, replicate) for state in range(4) for replicate in range(16)
    }
    assert len({row["signal_seed"] for row in design}) == 64
    assert min(row["signal_seed"] for row in design) == 20262004
    assert max(row["signal_seed"] for row in design) == 20262067


def test_manifest_freezes_design_and_source_hashes():
    manifest = _prepared_manifest()

    assert DIAGNOSTIC.OUT_DIR.name == "trackc_signal_action_sensitivity_20261106_v2"
    assert manifest["protocol_status"] == "prepared"
    assert manifest["settings"]["synthetic_state_count"] == 4
    assert manifest["settings"]["full_model_count"] == 81
    assert manifest["settings"]["matched_d0_model_count"] == 27
    assert manifest["outcome_metrics_included"] is False
    assert "src/216_trackc_robust_signal_ambiguity.py" in manifest["source_hashes"]


def test_manifest_refuses_a_seed_collision(monkeypatch):
    class FakeRegistry:
        ENTRIES = []

        @staticmethod
        def check_free(**kwargs):
            return ["collision"]

    monkeypatch.setattr(DIAGNOSTIC, "_seed_registry_module", lambda: FakeRegistry)

    with pytest.raises(ValueError, match="seed allocation failed"):
        DIAGNOSTIC.build_manifest()


def test_prepared_manifest_refuses_changed_project_title():
    manifest = _prepared_manifest()
    manifest["project_title"] = "different project"

    with pytest.raises(ValueError, match="project_title"):
        DIAGNOSTIC.validate_prepared_manifest(manifest)


def test_prepared_manifest_refuses_changed_outcome_declaration():
    manifest = _prepared_manifest()
    manifest["outcome_metrics_included"] = True

    with pytest.raises(ValueError, match="outcome_metrics_included"):
        DIAGNOSTIC.validate_prepared_manifest(manifest)


def test_capture_plan_runs_controller_on_synthetic_state(monkeypatch):
    monkeypatch.setattr(DIAGNOSTIC, "N_PARTICLES", 32)
    monkeypatch.setattr(DIAGNOSTIC, "PREDICTION_PARTICLES", 8)
    scenario = DIAGNOSTIC.trackc_sim.draw_scenarios(1, seed=987654)
    observations, _truth = DIAGNOSTIC.generate_integrated_signals(
        1,
        q=0.7,
        a=0.75,
        b=0.75,
        d=0.75,
        seed=987655,
        f0=scenario["f0"],
    )
    models = (
        DIAGNOSTIC.SignalModelParameters(0.7, 0.5, 0.5, 0.0),
        DIAGNOSTIC.SignalModelParameters(0.7, 0.5, 0.5, 0.5),
    )

    result = DIAGNOSTIC._capture_plan(models, observations, scenario)

    assert len(result["pattern"]) == DIAGNOSTIC.FORECAST_CYCLES
    assert set(result["pattern"]) <= {"0", "1"}
    assert np.isfinite(result["worst_case_value_margin"])
    assert result["elapsed_seconds"] >= 0
