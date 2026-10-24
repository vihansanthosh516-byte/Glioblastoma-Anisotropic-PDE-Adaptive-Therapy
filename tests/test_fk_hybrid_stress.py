from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_hybrid_forecast as hybrid  # noqa: E402
from src import fk_mask_inverse as fi  # noqa: E402


def ball(shape=(16, 16, 16), radius=3.0):
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    return sum((axis - c) ** 2 for axis, c in zip(coords, center)) <= radius**2


def forecast_kwargs():
    return {
        "prior": {"mean": np.log([0.04, 0.05]), "covariance": np.eye(2) * 0.5},
        "h": 2.0,
        "boundary_covariance": {"sigma_b_mm": 0.3, "tau_b_mm": 0.1, "ell_mm": 5.0},
        "population_log_volume_slope_per_day": 0.0,
        "shrinkage_time_constant_days": 90.0,
        "volume_variance_fit": {"intercept": 0.0, "slope_per_day": 0.0},
        "rho_bounds": (0.01, 0.1),
        "D_bounds": (0.01, 0.1),
        "n_coarse": 3,
        "n_refine": 3,
        "n_draws": 3,
        "n_anchors": 24,
        "seed": 41,
    }


def test_hybrid_forecast_runs_on_anisotropic_growth_history():
    initial = ball(radius=2.5)
    domain = np.ones(initial.shape, dtype=bool)
    tensor = np.zeros(initial.shape + (6,), dtype=np.float32)
    tensor[..., :3] = (0.08, 0.04, 0.02)
    days = np.asarray([0.0, 10.0, 20.0])
    masks = [initial, *fi.forward_masks_tensor(initial, domain, 0.04, tensor, days, 2.0)]

    result = hybrid.hybrid_posterior_predictive_forecast(
        masks, days, 30.0, domain=domain, **forecast_kwargs()
    )

    assert result["mask_draws"].shape == (3, *initial.shape)
    assert result["no_pde_ranker_mask_draws"].shape == result["mask_draws"].shape
    assert np.isfinite(result["parameter_draws_log"]).all()
    assert "no target mask used" in result["scope"]


def test_hybrid_size_component_can_forecast_shrinking_volume():
    domain = np.ones((16, 16, 16), dtype=bool)
    masks = [ball(radius=4.0), ball(radius=3.0), ball(radius=2.0)]
    days = np.asarray([0.0, 20.0, 40.0])

    result = hybrid.hybrid_posterior_predictive_forecast(
        masks, days, 60.0, domain=domain, **forecast_kwargs()
    )

    requested = result["requested_volume_draws_voxels"]
    assert result["patient_log_volume_slope_per_day"] < 0
    assert result["shrunk_log_volume_slope_per_day"] < 0
    assert np.all(requested < int(masks[-1].sum()))
    assert np.array_equal(result["final_volume_draws_voxels"], requested)
