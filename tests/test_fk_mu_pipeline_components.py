from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_mu_patient_runner import (
    _select_analysis_records,
    _summarize_posterior_parameters,
    run_mu_first_two_to_third_development_cv,
)
from src.fk_mu_scoring import dice_with_empty_rule, score_first_third_target
from src.fk_mu_training_nuisance import fit_mu_training_volume_nuisance
from src.fk_patient_cv import ForecastFailure, PatientSeries, run_first_third_90day_patient_cv
from src.fk_mu_development_data import MUDevelopmentCohort, is_loader_issued_cohort


def _mask(edge: int, shape=(16, 16, 16)) -> np.ndarray:
    mask = np.zeros(shape, dtype=bool)
    mask[:edge, :edge, :edge] = True
    return mask


def test_target_support_is_used_and_empty_dice_rule_is_explicit():
    support = np.zeros((4, 4, 4), dtype=bool)
    support[:2, :2, :2] = True
    prediction = np.zeros_like(support)
    prediction[0, 0, 0] = True
    target = np.zeros_like(support)
    target[0, 0, 0] = True

    result = score_first_third_target(
        {"fallback_masks": {"hybrid": prediction, "no_pde": prediction, "persistence": prediction}},
        {"patient_id": "patient-1", "mask": target},
        valid_support=support,
    )

    assert result["valid_support_voxels"] == 8
    assert result["dice_hybrid"] == 1.0
    assert result["paired_delta_dice_hybrid_minus_no_pde"] == 0.0
    assert result["surface_distance_status_hybrid"] == "undefined_partial_support"
    assert result["absolute_volume_error_mm3_hybrid"] == 0.0
    assert dice_with_empty_rule(np.zeros_like(target), np.zeros_like(target)) == 1.0
    assert dice_with_empty_rule(np.ones_like(target), np.zeros_like(target)) == 0.0

    outside = np.zeros_like(support)
    outside[3, 3, 3] = True
    with pytest.raises(ValueError, match="outside validated target support"):
        score_first_third_target(
            {"fallback_masks": {"hybrid": outside, "no_pde": prediction, "persistence": prediction}},
            {"patient_id": "patient-1", "mask": target},
            valid_support=support,
        )


def test_secondary_surface_and_volume_metrics_use_physical_spacing():
    target = np.zeros((7, 7, 7), dtype=bool)
    target[3, 3, 3] = True
    shifted = np.zeros_like(target)
    shifted[4, 3, 3] = True
    empty = np.zeros_like(target)

    result = score_first_third_target(
        {"fallback_masks": {"hybrid": shifted, "no_pde": target, "persistence": target}},
        {"patient_id": "patient-1", "mask": target},
        voxel_volume_mm3=1.0,
        voxel_spacing_mm=(2.0, 1.0, 0.5),
        valid_support=np.ones_like(target),
    )

    assert result["mean_surface_distance_mm_hybrid"] == 2.0
    assert result["hausdorff95_surface_distance_mm_hybrid"] == 2.0
    assert result["mean_surface_distance_mm_no_pde"] == 0.0
    assert result["absolute_volume_error_mm3_hybrid"] == 0.0
    assert result["relative_volume_error_hybrid"] == 0.0
    assert result["volume_crps_mm3_hybrid"] == 0.0
    assert result["volume_interval_90_covered_hybrid"] is True

    empty_result = score_first_third_target(
        {"fallback_masks": {"hybrid": target, "no_pde": empty, "persistence": empty}},
        {"patient_id": "patient-2", "mask": empty},
        voxel_volume_mm3=8.0,
        valid_support=np.ones_like(target),
    )
    assert empty_result["surface_distance_status_hybrid"] == "undefined_one_empty"
    assert empty_result["volume_error_mm3_hybrid"] == 8.0
    assert empty_result["relative_volume_error_hybrid"] is None


def test_probabilistic_brier_and_volume_coverage_use_frozen_draws():
    target = np.zeros((7, 7, 7), dtype=bool)
    target[3, 3, 3] = True
    shifted = np.zeros_like(target)
    shifted[4, 3, 3] = True
    draws = np.stack((target, shifted))

    result = score_first_third_target(
        {
            "forecast": {
                "mask_draws": draws,
                "no_pde_ranker_mask_draws": draws,
                "geometric_train_rate_mask": target,
            },
            "history_mask": target,
        },
        {"patient_id": "patient-1", "mask": target},
        valid_support=np.ones_like(target),
    )

    assert np.isclose(result["marginal_voxel_brier_hybrid"], 0.5 / target.size)
    assert result["volume_crps_mm3_hybrid"] == 0.0
    assert result["volume_interval_90_sample_count_hybrid"] == 2
    assert result["volume_interval_90_covered_hybrid"] is True


def test_posterior_summary_keeps_correlation_quantiles_and_bound_mass_json_safe():
    points = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    posterior = {
        "grid": {"points": points},
        "posterior_probabilities": np.full(4, 0.25),
        "parameter_draws_log": np.asarray([[0.0, 0.0], [1.0, 1.0]]),
        "posterior_mass_at_parameter_bounds": {
            "rho_lower": 0.1, "rho_upper": 0.2, "D_lower": 0.3, "D_upper": 0.4,
        },
        "coarse_grid_points": 9,
        "refined_grid_points": 16,
        "n_likelihood_evaluations": 16,
    }

    summary = _summarize_posterior_parameters(posterior)

    assert summary["parameter_order"] == ["log_rho", "log_D"]
    assert np.isclose(summary["posterior_correlation_log_rho_log_D"], 0.0)
    assert np.allclose(summary["posterior_mean_log_parameters"], [0.5, 0.5])
    assert summary["posterior_mass_at_parameter_bounds"]["D_upper"] == 0.4
    assert summary["posterior_draw_count"] == 2
    json.dumps(summary, allow_nan=False)


def test_training_volume_nuisance_uses_only_supplied_three_visit_patients():
    records = []
    gaps = ((40, 70), (45, 80), (50, 90), (55, 100), (60, 110))
    edges = ((2, 3, 4), (2, 4, 5), (3, 4, 5), (2, 3, 5), (3, 5, 6))
    for index, ((gap1, gap2), volume_edges) in enumerate(zip(gaps, edges)):
        records.append(PatientSeries(
            f"training-{index}",
            np.asarray([0.0, gap1, gap1 + gap2]),
            tuple(_mask(edge) for edge in volume_edges),
        ))

    fit = fit_mu_training_volume_nuisance(records)

    assert set(fit["training_patient_ids"]) == {record.patient_id for record in records}
    assert fit["audit_counts"]["n_training_patients"] == 5
    assert fit["audit_counts"]["n_positive_volume_patients"] == 5
    assert np.isfinite(fit["population_log_volume_slope_per_day"])
    assert fit["shrinkage_time_constant_days"] in (30.0, 90.0, 180.0, 365.0)


def test_fixed_day90_cv_keeps_each_target_out_of_forecast_and_uses_persistence_fallback():
    records = []
    for index in range(5):
        history1 = _mask(2)
        history2 = _mask(3)
        target = _mask(4)
        if index == 0:
            history1 = np.zeros_like(history1)
        records.append(PatientSeries(
            f"patient-{index}", np.asarray([0.0, 45.0, 130.0]), (history1, history2, target),
        ))
    folds = {record.patient_id: index for index, record in enumerate(records)}
    support = {record.patient_id: np.ones((16, 16, 16), dtype=bool) for record in records}
    forecasted = []

    def fit_nuisance(training_records):
        return {"training_ids": {record.patient_id for record in training_records}}

    def forecast(history, nuisance):
        if any(not mask.any() for mask in history.masks):
            raise ForecastFailure("synthetic empty history")
        assert history.target_index == 2
        assert len(history.masks) == 2
        assert history.forecast_day == 135.0
        assert history.patient_id not in nuisance["training_ids"]
        forecasted.append(history.patient_id)
        last = history.masks[-1]
        draws = np.stack([last, last])
        return {"forecast": {
                    "mask_draws": draws,
                    "no_pde_ranker_mask_draws": draws,
                    "geometric_train_rate_mask_draws": draws,
                },
                "history_mask": last}

    def fallback(history, nuisance):
        last = history.masks[-1]
        return {"fallback_masks": {
                    "hybrid": last,
                    "no_pde": last,
                    "persistence": last,
                    "geometric_train_rate": last,
                }}

    result = run_first_third_90day_patient_cv(
        records,
        folds,
        fit_nuisance,
        forecast,
        lambda prediction, target: score_first_third_target(
            prediction, target, valid_support=support[target.patient_id],
        ),
        fallback_forecast=fallback,
        n_folds=5,
    )

    assert len(result["rows"]) == 5
    assert len(forecasted) == 4
    assert {row["patient_id"] for row in result["rows"]} == {record.patient_id for record in records}
    assert all(row["scheduled_forecast_day"] == 135.0 for row in result["rows"])
    fallback_rows = [row for row in result["rows"] if row["forecast_status"] == "fallback"]
    assert len(fallback_rows) == 1 and fallback_rows[0]["patient_id"] == "patient-0"
    assert fallback_rows[0]["scores"]["paired_delta_dice_hybrid_minus_no_pde"] == 0.0
    assert all(row["scores"]["geometric_train_rate_status"] == "forecast"
               for row in result["rows"] if row["forecast_status"] == "ok")
    assert fallback_rows[0]["scores"]["geometric_train_rate_status"] == "persistence_fallback"


def test_gap_sensitivity_excludes_edge_days_before_refitting():
    records = []
    for index, gap in enumerate((76, 78, 79, 90, 101, 102, 104)):
        mask = _mask(2)
        records.append(PatientSeries(
            f"patient-{index}", np.asarray([0.0, 1.0, 1.0 + gap]), (mask, mask, mask),
        ))

    selected = _select_analysis_records(records, gap_sensitivity_only=True)

    assert [record.patient_id for record in selected] == ["patient-2", "patient-3", "patient-4"]


def test_public_runner_rejects_fabricated_cohort_even_with_true_audit_flags():
    fake = MUDevelopmentCohort(
        records=(), fold_by_patient={}, origin_domain_by_patient={}, target_valid_support_by_patient={},
        audit={"access": {"source_backed_registration_gate_passed": True,
                           "source_backed_target_support_gate_passed": True}},
    )
    assert not is_loader_issued_cohort(fake)

    with pytest.raises(ValueError, match="exact in-memory object"):
        run_mu_first_two_to_third_development_cv(
            fake,
            initial_boundary_covariance={"sigma_b_mm": 1.0, "tau_b_mm": 1.0, "ell_mm": 10.0},
        )
