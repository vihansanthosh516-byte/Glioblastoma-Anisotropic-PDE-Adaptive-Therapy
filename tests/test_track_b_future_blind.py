from __future__ import annotations

import importlib.util
import json
import sys
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "165_signed_transition_cv.py"
SPEC = importlib.util.spec_from_file_location("track_b_signed_cv_future_blind_test", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
SIGNED = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = SIGNED
SPEC.loader.exec_module(SIGNED)


def _synthetic_visits():
    shape = (12, 12, 12)
    rng = np.random.default_rng(1204)
    base = rng.normal(100.0, 20.0, size=shape).astype(np.float32)
    images = {
        "t1n": base,
        "t1c": base * 1.1 + 3.0,
        "t2w": base * 0.8 + 7.0,
        "t2f": base * 1.3 - 2.0,
    }
    previous_labels = np.zeros(shape, dtype=np.int16)
    previous_labels[5:7, 5:7, 5:7] = 1
    current_labels = previous_labels.copy()
    current_labels[7, 6, 6] = 1
    affine = np.eye(4, dtype=np.float64)
    previous = SIGNED.FEATURES.VisitVolume(images, previous_labels, affine)
    current = SIGNED.FEATURES.VisitVolume(images, current_labels, affine)
    future_without_transitions = current_labels.copy()
    future_with_transitions = current_labels.copy()
    future_with_transitions[5, 5, 5] = 0
    future_with_transitions[8, 6, 6] = 1
    parameters = SIGNED.FEATURES.DiffusionParameters(
        screening_length_mm=3.0,
        candidate_radius_mm=3.0,
        prior_smoothing_mm=1.0,
        solver_tolerance=1e-4,
        solver_max_iterations=500,
    )
    return previous, current, affine, future_without_transitions, future_with_transitions, parameters


def _interval(future_labels):
    previous, current, affine, _, _, parameters = _synthetic_visits()
    return SIGNED._interval_features(
        "SYNTHETIC-PATIENT",
        previous,
        current,
        future_labels,
        affine,
        parameters,
        "SYNTHETIC-PATIENT:2-3",
        previous_interval_days=30.0,
        forecast_interval_days=45.0,
        perturbation_radius_mm=1.0,
    )


def test_future_mask_changes_targets_but_not_forecast_origin_predictors():
    _, _, _, no_change_mask, changed_mask, _ = _synthetic_visits()
    no_change = _interval(no_change_mask)
    changed = _interval(changed_mask)

    for field in (
        "addition_features",
        "addition_prior",
        "removal_features",
        "removal_prior",
    ):
        assert np.array_equal(getattr(no_change, field), getattr(changed, field))
    assert no_change.history_count == changed.history_count
    assert no_change.target_count.addition_count == 0
    assert no_change.target_count.removal_count == 0
    assert changed.target_count.addition_count == 1
    assert changed.target_count.removal_count == 1


def test_future_mask_alone_does_not_change_count_forecasts():
    _, _, _, no_change_mask, changed_mask, _ = _synthetic_visits()
    no_change = _interval(no_change_mask)
    changed = _interval(changed_mask)

    class RecordingCountModel:
        def __init__(self):
            self.calls = []

        def forecast(self, **kwargs):
            self.calls.append(kwargs)
            scale = kwargs["future_interval_days"]
            add_multiplier = kwargs.get("future_addition_rate_multiplier", 1.0)
            remove_multiplier = kwargs.get("future_removal_rate_multiplier", 1.0)
            return SimpleNamespace(
                expected_additions=scale * add_multiplier,
                expected_removals=scale * remove_multiplier,
            )

    model = RecordingCountModel()
    no_change_forecasts = SIGNED._forecast_count_pair(
        model,
        no_change,
        addition_rate_multiplier=0.8,
        removal_rate_multiplier=0.9,
    )
    changed_forecasts = SIGNED._forecast_count_pair(
        model,
        changed,
        addition_rate_multiplier=0.8,
        removal_rate_multiplier=0.9,
    )

    assert no_change.target_count != changed.target_count
    assert no_change_forecasts == changed_forecasts
    assert model.calls[0] == model.calls[2]
    assert model.calls[1] == model.calls[3]


def test_future_event_counts_do_not_enter_count_forecast():
    _, _, _, _, changed_mask, _ = _synthetic_visits()
    row = _interval(changed_mask)
    changed_counts = replace(
        row.target_count,
        addition_count=0,
        removal_count=0,
    )
    alternative_row = replace(row, target_count=changed_counts)

    class RecordingCountModel:
        def __init__(self):
            self.calls = []

        def forecast(self, **kwargs):
            self.calls.append(kwargs)
            scale = kwargs["future_interval_days"]
            add_multiplier = kwargs.get("future_addition_rate_multiplier", 1.0)
            remove_multiplier = kwargs.get("future_removal_rate_multiplier", 1.0)
            return SimpleNamespace(
                expected_additions=scale * add_multiplier,
                expected_removals=scale * remove_multiplier,
            )

    model = RecordingCountModel()
    predictions = SIGNED._forecast_count_pair(
        model,
        row,
        addition_rate_multiplier=0.8,
        removal_rate_multiplier=0.9,
    )
    alternative_predictions = SIGNED._forecast_count_pair(
        model,
        alternative_row,
        addition_rate_multiplier=0.8,
        removal_rate_multiplier=0.9,
    )

    assert predictions == alternative_predictions
    assert model.calls[0] == model.calls[2]
    assert model.calls[1] == model.calls[3]
    assert "addition_count" not in model.calls[0]
    assert "removal_count" not in model.calls[0]
    assert model.calls[0]["future_interval_days"] == 45.0


def test_signed_runner_labels_realized_interval_as_conditional_horizon():
    for horizon in (SIGNED.FORECAST_HORIZON, SIGNED.BASE.FORECAST_HORIZON):
        assert horizon["known_at_forecast_origin"] is False
        assert "not a pre-schedule forecast" in horizon["deployment_interpretation"]


def test_signed_evaluation_completes_on_synthetic_patient_intervals():
    intervals = []
    for index in range(6):
        previous, current, affine, _, future, parameters = _synthetic_visits()
        intervals.append(SIGNED._interval_features(
            f"SYNTHETIC-{index}",
            previous,
            current,
            future,
            affine,
            parameters,
            f"SYNTHETIC-{index}:2-3",
            previous_interval_days=30.0,
            forecast_interval_days=45.0,
            perturbation_radius_mm=1.0,
        ))

    result = SIGNED.evaluate(intervals, n_splits=3, seed=260)

    assert result["n_patients"] == 6
    assert result["n_intervals"] == 6
    assert len(result["folds"]) == 3
    assert len(result["patient_bootstrap_contrasts"]) == 23


def test_output_input_paths_do_not_expose_external_user_paths(tmp_path: Path):
    repo = tmp_path / "synthetic-repo"
    repo.mkdir()
    external_home = tmp_path / "home" / "research-user" / "private-cohort"
    source_paths = {
        "data_root": SIGNED._logical_input_path(external_home, repo),
        "patient_manifest": SIGNED._logical_input_path(repo / "data" / "patients.csv", repo),
        "pair_manifest": SIGNED._logical_input_path(external_home / "pairs.csv", repo),
    }
    serialized_output = json.dumps({"source_paths": source_paths})

    assert source_paths == {
        "data_root": "external-input",
        "patient_manifest": "data/patients.csv",
        "pair_manifest": "external-input",
    }
    assert str(tmp_path) not in serialized_output
    assert "research-user" not in serialized_output
    assert str(external_home) not in serialized_output
