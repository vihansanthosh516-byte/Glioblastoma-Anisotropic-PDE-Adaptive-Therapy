import importlib.util
import math
import sys
from pathlib import Path

import pytest


REPO = Path(__file__).parents[1]
sys.path.insert(0, str(REPO / "src"))
SPEC = importlib.util.spec_from_file_location(
    "trackc_minimax_action_stability",
    REPO / "src" / "219_trackc_minimax_action_stability.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_model_sets_use_full_grid_and_matched_d_zero_projection():
    full, d0 = MODULE.build_model_sets()

    assert len(full) == 81
    assert len(d0) == 27
    assert set(d0) == {model for model in full if model.d == 0.0}


def test_capture_plan_returns_one_three_cycle_plan_and_all_action_values():
    scenario = MODULE._SIM.draw_scenarios(1, seed=20261990)
    observations = MODULE._AMBIGUITY.SignalObservations(
        molecular_report=[True],
        mri_score=[0.75],
    )
    models = (
        MODULE._AMBIGUITY.SignalModelParameters(0.8, 0.3, 0.6, -0.5),
        MODULE._AMBIGUITY.SignalModelParameters(0.8, 0.3, 0.6, 0.0),
    )

    result = MODULE._capture_plan(models, observations, scenario, 2, 930000001)

    assert len(result["pattern"]) == 3
    assert len(result["worst_case_values_by_pattern"]) == 8
    assert all(math.isfinite(value) for value in result["worst_case_values_by_pattern"])
    assert result["worst_case_value_margin"] >= 0.0


def test_summary_reports_adjacent_particle_agreement_without_outcome_claims():
    rows = []
    for model_set, by_seed in (
        ("full_81", {
            1: ("100", "100", "010"),
            2: ("010", "010", "010"),
            3: ("001", "111", "001"),
        }),
        ("matched_d0_27", {
            1: ("100", "010", "010"),
            2: ("100", "100", "010"),
            3: ("010", "100", "100"),
        }),
    ):
        for prediction_seed, patterns in by_seed.items():
            for particle_count, pattern in zip((16, 32, 48), patterns):
                rows.append({
                    "model_set": model_set,
                    "anchor": "discordant",
                    "prediction_seed": prediction_seed,
                    "prediction_particles": particle_count,
                    "pattern": pattern,
                })

    summary = MODULE.summarize_action_rows(rows)

    by_model = {row["model_set"]: row for row in summary["adjacent_agreement"]}
    assert by_model["full_81"]["agreements"] == [2, 1]
    assert by_model["full_81"]["fraction_by_step"] == [pytest.approx(2 / 3), pytest.approx(1 / 3)]
    assert by_model["matched_d0_27"]["agreements"] == [1, 2]
    assert summary["outcome_metrics_included"] is False


def test_summary_rejects_duplicate_or_missing_particle_levels():
    rows = [
        {"model_set": "full_81", "anchor": "a", "prediction_seed": 1,
         "prediction_particles": 16, "pattern": "100"},
        {"model_set": "full_81", "anchor": "a", "prediction_seed": 1,
         "prediction_particles": 16, "pattern": "010"},
    ]

    with pytest.raises(ValueError, match="duplicate"):
        MODULE.summarize_action_rows(rows)
