from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_full_hybrid_gate as gate  # noqa: E402


def test_manifest_balances_growth_shrinkage_anisotropy_and_outer_folds():
    rows = gate.build_synthetic_manifest(replicates_per_condition=2, n_folds=2, seed=91)

    assert len(rows) == 8
    assert len({row["case_id"] for row in rows}) == 8
    assert {(row["diffusion_truth"], row["volume_trend"]) for row in rows} == {
        ("isotropic", "growth"), ("isotropic", "shrinkage"),
        ("anisotropic", "growth"), ("anisotropic", "shrinkage"),
    }
    for condition in {(row["diffusion_truth"], row["volume_trend"]) for row in rows}:
        assert {row["fold"] for row in rows
                if (row["diffusion_truth"], row["volume_trend"]) == condition} == {0, 1}


def test_synthetic_patient_generation_is_reproducible_and_spans_shrinkage():
    row = next(row for row in gate.build_synthetic_manifest(2, 2, 19)
               if row["volume_trend"] == "shrinkage")
    first, first_truth = gate.generate_synthetic_series(row)
    second, second_truth = gate.generate_synthetic_series(row)

    assert first.patient_id == second.patient_id
    assert np.array_equal(first.days, second.days)
    assert all(np.array_equal(a, b) for a, b in zip(first.masks, second.masks))
    assert first.masks[-1].sum() < first.masks[0].sum()
    assert first_truth == second_truth


def test_training_nuisance_records_only_fit_cohort_and_uses_training_series():
    rows = gate.build_synthetic_manifest(3, 3, 23)
    records = [gate.generate_synthetic_series(row)[0] for row in rows]
    training = records[:8]
    nuisance = gate.fit_training_nuisance(training)

    assert nuisance["training_patient_ids"] == tuple(sorted(record.patient_id for record in training))
    assert np.isfinite(nuisance["population_log_volume_slope_per_day"])
    assert nuisance["shrinkage_time_constant_days"] in (30.0, 90.0, 180.0, 365.0)
    assert nuisance["volume_variance_fit"]["intercept"] >= 0
    assert nuisance["volume_variance_fit"]["slope_per_day"] >= 0


def test_outer_gate_scores_each_synthetic_patient_once_and_withholds_target_mask(monkeypatch):
    observed = []

    def cheap_forecast(masks, days, target_day, **kwargs):
        last = np.asarray(masks[-1], dtype=bool)
        observed.append((len(masks), float(target_day), kwargs["population_log_volume_slope_per_day"]))
        draws = np.stack([last, last], axis=0)
        return {"mask_draws": draws, "no_pde_ranker_mask_draws": draws.copy(),
                "requested_volume_draws_voxels": np.asarray([last.sum(), last.sum()])}

    result = gate.run_synthetic_gate(
        replicates_per_condition=2, n_folds=2, seed=37, forecast_function=cheap_forecast,
    )

    assert len(result["rows"]) == 8
    assert len({row["patient_id"] for row in result["rows"]}) == 8
    assert len(observed) == 8
    assert all(history_count == 3 and target_day > 0 for history_count, target_day, _ in observed)
    assert result["summary"]["n_patients"] == 8
    assert set(result["summary"]["by_condition"]) == {
        "isotropic|growth", "isotropic|shrinkage",
        "anisotropic|growth", "anisotropic|shrinkage",
    }
