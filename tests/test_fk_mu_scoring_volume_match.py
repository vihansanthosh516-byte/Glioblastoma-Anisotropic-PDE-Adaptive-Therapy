from __future__ import annotations

import numpy as np
import pytest

from src.fk_mu_scoring import score_first_third_target


def test_no_pde_point_mask_matches_hybrid_thresholded_volume():
    hybrid_draws = np.zeros((2, 2, 2, 2), dtype=bool)
    hybrid_draws[0].ravel()[[0, 1]] = True
    hybrid_draws[1].ravel()[[0, 2]] = True
    no_pde_draws = np.zeros_like(hybrid_draws)
    no_pde_draws[0].ravel()[[0, 1]] = True
    no_pde_draws[1].ravel()[[2, 3]] = True
    target = np.zeros((2, 2, 2), dtype=bool)
    target.ravel()[[0, 1, 2]] = True

    assert np.all(hybrid_draws.reshape(2, -1).sum(axis=1) == 2)
    assert np.all(no_pde_draws.reshape(2, -1).sum(axis=1) == 2)
    result = score_first_third_target(
        {
            "forecast": {
                "mask_draws": hybrid_draws,
                "no_pde_ranker_mask_draws": no_pde_draws,
            },
            "history_mask": target,
        },
        {"patient_id": "synthetic-patient", "mask": target},
        valid_support=np.ones_like(target),
    )

    assert result["hybrid_prediction_voxels"] == 3
    assert result["no_pde_prediction_voxels"] == 3
    assert result["dice_hybrid"] == result["dice_no_pde_matched_volume"] == 1.0
    assert result["paired_delta_dice_hybrid_minus_no_pde"] == 0.0
    assert result["no_pde_matched_volume_point_forecast_rule"] == (
        "rank no-PDE inclusion probabilities to hybrid point-mask voxel count; stable flat-index ties"
    )


def test_fallback_masks_reject_unequal_hybrid_and_no_pde_volumes():
    hybrid = np.zeros((4, 4, 4), dtype=bool)
    hybrid[:2, :2, :2] = True
    no_pde = np.zeros_like(hybrid)
    no_pde[0, 0, 0] = True

    with pytest.raises(ValueError, match="fallback masks must have equal volumes"):
        score_first_third_target(
            {"fallback_masks": {"hybrid": hybrid, "no_pde": no_pde, "persistence": hybrid}},
            {"patient_id": "synthetic-patient", "mask": hybrid},
        )
