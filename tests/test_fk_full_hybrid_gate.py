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


def test_summary_weights_patients_equally_when_visit_counts_differ():
    score_template = {
        "energy_score_difference_hybrid_minus_no_pde": 0.0,
        "hybrid_consensus_dice": 0.8,
        "no_pde_consensus_dice": 0.7,
        "hybrid_variogram_score": 1.0,
        "no_pde_variogram_score": 1.2,
        "hybrid_volume_crps_voxels": 2.0,
        "no_pde_volume_crps_voxels": 3.0,
        "diffusion_truth": "isotropic",
        "volume_trend": "growth",
    }
    rows = []
    for patient_id, differences in (("p1", (0.0, 0.0)), ("p2", (-3.0,))):
        for difference in differences:
            scores = dict(score_template)
            scores["energy_score_difference_hybrid_minus_no_pde"] = difference
            rows.append({"patient_id": patient_id, "scores": scores})

    summary = gate._summarize(rows)

    assert summary["n_patients"] == 2
    assert summary["n_forecasts"] == 3
    assert summary["overall"]["n_patients"] == 2
    assert summary["overall"]["n_forecasts"] == 3
    assert summary["overall"]["mean_energy_score_difference_hybrid_minus_no_pde"] == -1.5
    assert summary["overall"]["hybrid_lower_energy_count"] == 1


def test_outer_gate_scores_planned_origins_and_withholds_target_mask(monkeypatch):
    observed = []
    mismatched_covariance = {"sigma_b_mm": 0.2, "tau_b_mm": 0.05, "ell_mm": 2.0}

    def cheap_forecast(masks, days, target_day, **kwargs):
        last = np.asarray(masks[-1], dtype=bool)
        observed.append((len(masks), float(target_day), kwargs["population_log_volume_slope_per_day"],
                         kwargs["boundary_covariance"]))
        draws = np.stack([last, last], axis=0)
        return {"mask_draws": draws, "no_pde_ranker_mask_draws": draws.copy(),
                "requested_volume_draws_voxels": np.asarray([last.sum(), last.sum()])}

    result = gate.run_synthetic_gate(
        replicates_per_condition=2, n_folds=2, seed=37, forecast_function=cheap_forecast,
        forecast_nuisance_overrides={"boundary_covariance": mismatched_covariance},
        boundary_calibration_mode="oracle",
    )

    assert len(result["rows"]) == 16
    assert len({row["patient_id"] for row in result["rows"]}) == 8
    assert len(observed) == 16
    assert all(history_count in (2, 3) and target_day > 0 and covariance == mismatched_covariance
               for history_count, target_day, _, covariance in observed)
    assert result["metadata"]["boundary_covariance_truth"] != \
        result["metadata"]["boundary_covariance_model_input"]
    assert not result["metadata"]["boundary_covariance_is_oracle_input"]
    assert result["metadata"]["boundary_covariance_calibration_status"] == "not_run_oracle"
    assert result["summary"]["n_patients"] == 8
    assert result["summary"]["n_forecasts"] == 16
    assert result["summary"]["overall"]["n_patients"] == 8
    assert result["summary"]["overall"]["n_forecasts"] == 16
    assert set(result["summary"]["by_condition"]) == {
        "isotropic|growth", "isotropic|shrinkage",
        "anisotropic|growth", "anisotropic|shrinkage",
    }


def test_outer_gate_can_estimate_boundary_covariance_from_inner_heldout_patients(monkeypatch):
    observed_covariances = []

    def cheap_series(row):
        coordinates = np.indices((16, 16, 16), dtype=float)
        volume_sign = 1.0 if row["volume_trend"] == "growth" else -0.5
        masks = []
        for scan_index in range(4):
            center_x = 7.0 + 0.2 * scan_index
            radius = 3.8 + volume_sign * 0.15 * scan_index
            mask = ((coordinates[0] - center_x) ** 2
                    + (coordinates[1] - 7.0) ** 2
                    + (coordinates[2] - 8.0) ** 2) <= radius ** 2
            masks.append(mask)
        truth = {"diffusion_truth": row["diffusion_truth"], "volume_trend": row["volume_trend"]}
        return gate.PatientSeries(row["case_id"], np.asarray([0.0, 15.0, 35.0, 60.0]), tuple(masks)), truth

    monkeypatch.setattr(gate, "generate_synthetic_series", cheap_series)

    def cheap_forecast(masks, days, target_day, **kwargs):
        last = np.asarray(masks[-1], dtype=bool)
        observed_covariances.append(dict(kwargs["boundary_covariance"]))
        draws = np.stack([last, last], axis=0)
        return {"mask_draws": draws, "no_pde_ranker_mask_draws": draws.copy(),
                "requested_volume_draws_voxels": np.asarray([last.sum(), last.sum()]),
                "pre_boundary_probability": last.astype(float)}

    result = gate.run_synthetic_gate(
        replicates_per_condition=3, n_folds=3, seed=47, forecast_function=cheap_forecast,
        forecast_config={"n_anchors": 8}, boundary_calibration_mode="cross_fitted", inner_n_folds=2,
    )

    metadata = result["metadata"]
    assert len(result["rows"]) == 24
    assert len(observed_covariances) == 72  # 48 inner calibration plus 24 outer forecasts.
    assert metadata["boundary_calibration_mode"] == "cross_fitted"
    assert metadata["boundary_covariance_model_input"] is None
    assert not metadata["boundary_covariance_is_oracle_input"]
    assert metadata["boundary_covariance_crossfit_start_input"] == gate.BOUNDARY_CALIBRATION_START
    calibration_calls = (observed_covariances[0:16] + observed_covariances[24:40]
                          + observed_covariances[48:64])
    assert calibration_calls == [gate.BOUNDARY_CALIBRATION_START] * 48
    fits = metadata["boundary_covariance_estimates_by_outer_fold"]
    assert len(fits) == 3
    boundary_folds = sum(fit["selection_boundary"] for fit in fits)
    assert metadata["boundary_covariance_search_boundary_fold_count"] == boundary_folds
    assert metadata["boundary_covariance_calibration_status"] == (
        "inconclusive_search_boundary" if boundary_folds else "all_folds_interior"
    )
    assert "diagnostic only" in metadata["boundary_covariance_status_rule"]
    assert all(fit["n_outer_training_patients"] == 8 for fit in fits)
    assert all(fit["n_calibration_patients"] == 8 for fit in fits)
    assert all(fit["n_calibration_origins"] == 16 for fit in fits)
    assert all(set(fit["estimate"]) == {"sigma_b_mm", "tau_b_mm", "ell_mm"} for fit in fits)
    assert all(np.isfinite(list(fit["estimate"].values())).all() for fit in fits)
    assert all(value >= 0 for fit in fits for value in fit["estimate"].values())


def test_cross_fitted_mode_rejects_a_fixed_boundary_covariance_override():
    with np.testing.assert_raises_regex(ValueError, "estimates it from training patients"):
        gate.run_synthetic_gate(
            forecast_nuisance_overrides={"boundary_covariance": gate.BOUNDARY_TRUTH},
            boundary_calibration_mode="cross_fitted",
        )
