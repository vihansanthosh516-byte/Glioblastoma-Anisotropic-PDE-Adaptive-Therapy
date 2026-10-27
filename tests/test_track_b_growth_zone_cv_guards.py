from __future__ import annotations

import importlib.util
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "162_growth_zone_cv.py"
SPEC = importlib.util.spec_from_file_location("track_b_growth_zone_cv_guards_test", MODULE_PATH)
CV = importlib.util.module_from_spec(SPEC)
assert SPEC is not None and SPEC.loader is not None
sys.modules[SPEC.name] = CV
SPEC.loader.exec_module(CV)


def _raw_sample(*, prior_instability: float = 10.0, prior_addition_count: int = 3) -> CV.RawSample:
    return CV.RawSample(
        sample_id="SYNTHETIC:1-2-3",
        patient_id="SYNTHETIC",
        features=np.asarray([[1.0, 0.6], [2.0, 0.8]], dtype=np.float32),
        log_pde=np.asarray([-1.0, -0.5], dtype=np.float32),
        target=np.asarray([0, 1], dtype=np.uint8),
        log_instability=np.asarray([0.0, 2.0], dtype=np.float32),
        prior_instability=np.asarray([0.0, prior_instability], dtype=np.float32),
        prior_addition_count=prior_addition_count,
        prior_candidate_voxels=20,
        prior_interval_days=30.0,
        forecast_interval_days=45.0,
        candidate_count=2,
        all_addition_count=1,
        outside_roi_addition_count=0,
    )


def test_inner_ap_selection_fails_when_every_fold_is_undefined() -> None:
    with pytest.raises(ValueError, match="undefined for all ridge candidates"):
        CV._select_best_inner_ap({0.1: [float("nan")], 1.0: [float("nan")]}, "ridge")


def test_inner_ap_selection_ignores_undefined_folds_and_keeps_grid_tie_order() -> None:
    selected = CV._select_best_inner_ap(
        {0.1: [0.7, float("nan")], 1.0: [0.7, float("nan")]},
        "ridge",
    )
    assert selected == 0.1


def test_no_history_materialization_removes_feature_gate_and_count_inputs() -> None:
    sample = _raw_sample()
    no_history = sample.materialize((2.0, 10.0), include_history=False)
    changed_history = replace(
        sample,
        prior_instability=np.asarray([0.0, 1000.0], dtype=np.float32),
        prior_addition_count=19,
        prior_candidate_voxels=50,
        prior_interval_days=120.0,
    ).materialize((2.0, 10.0), include_history=False)

    assert np.array_equal(no_history.features[:, -1], np.zeros(2, dtype=np.float32))
    assert no_history.prior_addition_count == 0
    assert np.array_equal(no_history.features, changed_history.features)
    assert np.array_equal(no_history.gate, changed_history.gate)
    assert np.array_equal(no_history.log_pde, changed_history.log_pde)
    assert np.array_equal(no_history.target, changed_history.target)
    assert no_history.prior_candidate_voxels == changed_history.prior_candidate_voxels == 1
    assert no_history.prior_interval_days == changed_history.prior_interval_days == 1.0
    assert no_history.forecast_interval_days == changed_history.forecast_interval_days == sample.forecast_interval_days
    assert not np.array_equal(sample.materialize((2.0, 10.0)).gate, no_history.gate)


def test_no_history_inner_selection_uses_only_inner_training_cohort_rate(monkeypatch) -> None:
    samples = [
        replace(
            _raw_sample(prior_instability=10.0 + index, prior_addition_count=3 + index),
            patient_id=f"SYNTHETIC-{index}",
            sample_id=f"SYNTHETIC-{index}:1-2-3",
        )
        for index in range(3)
    ]
    fitted_history_patient_ids: list[frozenset[str]] = []
    validation_patient_to_training_ids: list[tuple[frozenset[str], str]] = []
    materialized_training_samples: list[CV.MODELS.VoxelSample] = []
    active_training_ids: frozenset[str] = frozenset()

    class FakeModel:
        def predict_score(self, features, log_pde, gate):
            return np.asarray([1.0, -1.0], dtype=np.float64)

    def fit_count_prior(training_samples):
        ids = frozenset(sample.patient_id for sample in training_samples)
        fitted_history_patient_ids.append(ids)
        return ids

    def fit_residual_model(training_samples, **_kwargs):
        materialized_training_samples.extend(training_samples)
        return FakeModel()

    def predict_cohort_rate_count(training_patient_ids, _interval_days, _candidate_count):
        nonlocal active_training_ids
        active_training_ids = training_patient_ids
        return 1.0

    def forbidden_patient_history_count(*_args, **_kwargs):
        pytest.fail("no-history tuning must not use patient-specific history counts")

    def sample_ap(sample, _score, _count):
        validation_patient_to_training_ids.append((active_training_ids, sample.patient_id))
        return np.asarray([0.25, 0.75]), {"average_precision": 0.5}

    monkeypatch.setattr(CV.MODELS, "fit_count_prior", fit_count_prior)
    monkeypatch.setattr(CV.MODELS, "fit_residual_model", fit_residual_model)
    monkeypatch.setattr(CV.MODELS, "predict_cohort_rate_count", predict_cohort_rate_count)
    monkeypatch.setattr(CV.MODELS, "predict_addition_count", forbidden_patient_history_count)
    monkeypatch.setattr(CV, "_sample_ap", sample_ap)

    selected = CV._choose_ridge(
        samples,
        ridge_grid=[0.1, 1.0],
        inner_folds=2,
        seed=17,
        include_history=False,
    )

    assert selected == 0.1
    assert len(fitted_history_patient_ids) == 2
    assert all(len(ids) < len(samples) for ids in fitted_history_patient_ids)
    assert all(train_ids.isdisjoint({validation_id}) for train_ids, validation_id in validation_patient_to_training_ids)
    assert all(sample.prior_addition_count == 0 for sample in materialized_training_samples)
    assert all(np.all(sample.features[:, -1] == 0.0) for sample in materialized_training_samples)


def test_outer_fold_fits_no_history_model_with_its_own_selected_ridge(monkeypatch) -> None:
    train_raw = [
        replace(_raw_sample(), patient_id=f"TRAIN-{index}", sample_id=f"TRAIN-{index}:1-2-3")
        for index in range(3)
    ]
    test_raw = [replace(_raw_sample(), patient_id="TEST", sample_id="TEST:1-2-3")]
    selections: list[bool] = []
    fits: list[tuple[list[CV.MODELS.VoxelSample], float]] = []

    class FakeModel:
        def predict_score(self, features, *_args, **_kwargs):
            return np.asarray([1.0, -1.0], dtype=np.float64)

    def choose_ridge(_samples, _grid, _folds, _seed, *, include_history=True):
        selections.append(include_history)
        return 0.1 if include_history else 2.0

    def fit_residual_model(samples, *, ridge, **_kwargs):
        fits.append((list(samples), float(ridge)))
        return FakeModel()

    def fit_count_prior(_samples):
        return object()

    def count_prediction(*_args):
        return 1.0

    def sample_ap(_sample, _score, count):
        return np.asarray([0.25, 0.75]), {
            "average_precision": 0.5,
            "auc": 0.5,
            "brier": 0.25,
            "predicted_mass": float(count),
        }

    monkeypatch.setattr(CV, "_choose_ridge", choose_ridge)
    monkeypatch.setattr(CV, "_choose_tree_leaves", lambda *_args, **_kwargs: 7)
    monkeypatch.setattr(CV.MODELS, "fit_count_prior", fit_count_prior)
    monkeypatch.setattr(CV.MODELS, "fit_residual_model", fit_residual_model)
    monkeypatch.setattr(CV.MODELS, "fit_tree_baseline", lambda *_args, **_kwargs: FakeModel())
    monkeypatch.setattr(CV.MODELS, "predict_addition_count", count_prediction)
    monkeypatch.setattr(CV.MODELS, "predict_cohort_rate_count", count_prediction)
    monkeypatch.setattr(CV, "_sample_ap", sample_ap)

    _rows, fold_summary = CV._evaluate_outer_fold(
        train_raw,
        test_raw,
        ridge_grid=[0.1, 2.0],
        inner_folds=2,
        seed=23,
        fold=1,
    )

    assert selections == [True, False]
    assert fold_summary["selected_no_history_ridge"] == 2.0
    no_history_training_samples, no_history_ridge = fits[-1]
    assert no_history_ridge == 2.0
    assert all(sample.prior_addition_count == 0 for sample in no_history_training_samples)
    assert all(np.all(sample.features[:, -1] == 0.0) for sample in no_history_training_samples)
