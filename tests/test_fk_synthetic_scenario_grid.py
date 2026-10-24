from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_synthetic_scenario_grid import build_split_manifest, summarize_rows  # noqa: E402


def test_split_manifest_is_balanced_seeded_and_disjoint():
    first = build_split_manifest(replicates_per_split=2, seed=7)
    second = build_split_manifest(replicates_per_split=2, seed=7)

    assert first == second
    assert len(first) == 16
    assert len({row["case_id"] for row in first}) == len(first)
    assert len({row["seed"] for row in first}) == len(first)
    for truth in ("isotropic", "anisotropic"):
        for sigma in (0.0, 0.8):
            for split in ("calibration", "evaluation"):
                matches = [row for row in first if row["truth"] == truth and row["sigma_mm"] == sigma
                           and row["split"] == split]
                assert len(matches) == 2


def test_summary_calibrates_on_calibration_rows_only():
    rows = [
        {"case_id": "cal-0", "split": "calibration", "pde_error": 0.10, "pde_dice": 0.90,
         "truth": "isotropic", "sigma_mm": 0.0,
         "pde_hd95_mm": 1.0, "pde_symmetric_volume_error_pct": 10.0,
         "persistence_error": 0.20, "persistence_dice": 0.80,
         "persistence_hd95_mm": 2.0, "persistence_symmetric_volume_error_pct": 20.0,
         "volume_trend_error": 0.15, "volume_trend_dice": 0.85,
         "volume_trend_hd95_mm": 1.5, "volume_trend_symmetric_volume_error_pct": 15.0,
         "ensemble_disagreement": 0.20, "training_spatial_residual": 0.10,
         "training_volume_residual_score": 0.30, "composite_uncertainty": 0.20},
        {"case_id": "cal-1", "split": "calibration", "pde_error": 0.30, "pde_dice": 0.70,
         "truth": "isotropic", "sigma_mm": 0.0,
         "pde_hd95_mm": 3.0, "pde_symmetric_volume_error_pct": 30.0,
         "persistence_error": 0.10, "persistence_dice": 0.90,
         "persistence_hd95_mm": 1.0, "persistence_symmetric_volume_error_pct": 10.0,
         "volume_trend_error": 0.20, "volume_trend_dice": 0.80,
         "volume_trend_hd95_mm": 2.0, "volume_trend_symmetric_volume_error_pct": 20.0,
         "ensemble_disagreement": 0.40, "training_spatial_residual": 0.30,
         "training_volume_residual_score": 0.10, "composite_uncertainty": 0.30},
        {"case_id": "eval-0", "split": "evaluation", "pde_error": 0.10, "pde_dice": 0.90,
         "truth": "isotropic", "sigma_mm": 0.0,
         "pde_hd95_mm": 2.0, "pde_symmetric_volume_error_pct": 20.0,
         "persistence_error": 0.25, "persistence_dice": 0.75,
         "persistence_hd95_mm": 2.5, "persistence_symmetric_volume_error_pct": 25.0,
         "volume_trend_error": 0.15, "volume_trend_dice": 0.85,
         "volume_trend_hd95_mm": 1.5, "volume_trend_symmetric_volume_error_pct": 15.0,
         "ensemble_disagreement": 0.25, "training_spatial_residual": 0.20,
         "training_volume_residual_score": 0.25, "composite_uncertainty": 0.25},
    ]
    for row in rows:
        row["pde_shape_volume_trend_error"] = row["pde_error"]
        row["pde_shape_volume_trend_dice"] = row["pde_dice"]
        row["pde_shape_volume_trend_hd95_mm"] = row["pde_hd95_mm"]
        row["pde_shape_volume_trend_symmetric_volume_error_pct"] = row["pde_symmetric_volume_error_pct"]
        row["hybrid_error"] = row["pde_shape_volume_trend_error"]
    first = summarize_rows(rows, alpha=0.7)
    changed_test = [dict(row) for row in rows]
    changed_test[-1]["hybrid_error"] = 0.95
    changed_test[-1]["pde_shape_volume_trend_error"] = 0.95
    changed_test[-1]["pde_shape_volume_trend_dice"] = 0.05
    second = summarize_rows(changed_test, alpha=0.7)

    for name in first["uncertainty_by_condition"]["isotropic|sigma_mm=0"]:
        assert first["uncertainty_by_condition"]["isotropic|sigma_mm=0"][name]["calibration_quantile"] == \
            second["uncertainty_by_condition"]["isotropic|sigma_mm=0"][name]["calibration_quantile"]
    assert first["uncertainty"]["training_spatial_residual"]["evaluation_coverage"] != \
        second["uncertainty"]["training_spatial_residual"]["evaluation_coverage"]
    assert set(first["forecast_comparison"]) == {"pde", "persistence", "volume_trend", "pde_shape_volume_trend"}


def test_uncertainty_calibration_is_separate_for_each_generating_condition():
    def row(case_id, split, truth, error, score):
        return {
            "case_id": case_id,
            "split": split,
            "truth": truth,
            "sigma_mm": 0.0,
            "pde_error": error,
            "pde_dice": 1.0 - error,
            "pde_hd95_mm": 1.0,
            "pde_symmetric_volume_error_pct": 1.0,
            "persistence_error": error,
            "persistence_dice": 1.0 - error,
            "persistence_hd95_mm": 1.0,
            "persistence_symmetric_volume_error_pct": 1.0,
            "volume_trend_error": error,
            "volume_trend_dice": 1.0 - error,
            "volume_trend_hd95_mm": 1.0,
            "volume_trend_symmetric_volume_error_pct": 1.0,
            "pde_shape_volume_trend_error": error,
            "pde_shape_volume_trend_dice": 1.0 - error,
            "pde_shape_volume_trend_hd95_mm": 1.0,
            "pde_shape_volume_trend_symmetric_volume_error_pct": 1.0,
            "hybrid_error": error,
            "ensemble_disagreement": score,
            "training_spatial_residual": score,
            "training_volume_residual_score": score,
            "composite_uncertainty": score,
        }

    rows = [
        row("iso-cal-0", "calibration", "isotropic", 0.10, 0.10),
        row("iso-cal-1", "calibration", "isotropic", 0.20, 0.20),
        row("iso-eval", "evaluation", "isotropic", 0.15, 0.10),
        row("aniso-cal-0", "calibration", "anisotropic", 0.80, 0.10),
        row("aniso-cal-1", "calibration", "anisotropic", 0.90, 0.20),
        row("aniso-eval", "evaluation", "anisotropic", 0.85, 0.10),
    ]

    result = summarize_rows(rows, alpha=0.5)
    uncertainty = result["uncertainty_by_condition"]

    assert np.isclose(
        uncertainty["isotropic|sigma_mm=0"]["training_spatial_residual"]["calibration_quantile"], 1.0
    )
    assert np.isclose(
        uncertainty["anisotropic|sigma_mm=0"]["training_spatial_residual"]["calibration_quantile"], 8.0
    )
    assert uncertainty["isotropic|sigma_mm=0"]["training_spatial_residual"]["n_calibration"] == 2
    assert uncertainty["anisotropic|sigma_mm=0"]["training_spatial_residual"]["n_calibration"] == 2
