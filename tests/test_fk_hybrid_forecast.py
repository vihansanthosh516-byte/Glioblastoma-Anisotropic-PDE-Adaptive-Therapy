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
    return sum((axis - center_i) ** 2 for axis, center_i in zip(coords, center)) <= radius**2


def install_fixed_posterior(monkeypatch, n_draws):
    logs = np.tile(np.log([0.03, 0.05]), (n_draws, 1))
    monkeypatch.setattr(
        hybrid, "infer_joint_log_parameter_posterior",
        lambda *args, **kwargs: {"parameter_draws_log": logs, "marker": "posterior-called"},
    )
    monkeypatch.setattr(
        hybrid.fi, "forward_concentrations",
        lambda mask, domain, rho, diffusion, days, h: [mask.astype(float)],
    )


def forecast(monkeypatch, *, seed, n_draws=4, uncertain=False, population_slope=0.0):
    mask = ball()
    install_fixed_posterior(monkeypatch, n_draws)
    boundary_covariance = {
        "sigma_b_mm": 0.4 if uncertain else 0.0,
        "tau_b_mm": 0.1 if uncertain else 0.0,
        "ell_mm": 5.0,
    }
    variance_fit = {"intercept": 0.02, "slope_per_day": 0.0001} if uncertain else {
        "intercept": 0.0, "slope_per_day": 0.0,
    }
    return hybrid.hybrid_posterior_predictive_forecast(
        [mask.copy(), mask.copy(), mask.copy()], [0.0, 20.0, 40.0], 60.0,
        prior={"mean": np.log([0.03, 0.05]), "covariance": np.eye(2)}, h=1.0,
        boundary_covariance=boundary_covariance,
        population_log_volume_slope_per_day=population_slope, shrinkage_time_constant_days=90.0,
        volume_variance_fit=variance_fit,
        n_draws=n_draws, n_anchors=32, seed=seed,
    )


def test_history_only_hybrid_draws_preserve_sampled_volume_without_boundary_noise(monkeypatch):
    result = forecast(monkeypatch, seed=13, n_draws=4, population_slope=0.01)
    counts = result["mask_draws"].reshape(4, -1).sum(axis=1)
    baseline_counts = result["no_pde_ranker_mask_draws"].reshape(4, -1).sum(axis=1)
    geometric_count = int(result["geometric_train_rate_mask"].sum())

    assert result["mask_draws"].shape == (4, 16, 16, 16)
    assert result["no_pde_ranker_mask_draws"].shape == result["mask_draws"].shape
    assert result["geometric_train_rate_mask"].shape == result["mask_draws"].shape[1:]
    assert np.array_equal(counts, result["requested_volume_draws_voxels"])
    assert np.array_equal(baseline_counts, result["requested_volume_draws_voxels"])
    expected_geometric_volume = int(np.rint(ball().sum() * np.exp(0.01 * 20.0)))
    assert result["geometric_train_rate_requested_volume_voxels"] == expected_geometric_volume
    assert geometric_count == expected_geometric_volume
    assert np.array_equal(counts, result["final_volume_draws_voxels"])
    assert result["posterior"]["marker"] == "posterior-called"
    assert "no target mask used" in result["scope"]


def test_hybrid_forecast_seed_reproduces_volume_and_boundary_draws(monkeypatch):
    first = forecast(monkeypatch, seed=29, n_draws=4, uncertain=True)
    second = forecast(monkeypatch, seed=29, n_draws=4, uncertain=True)

    assert np.array_equal(first["mask_draws"], second["mask_draws"])
    assert np.array_equal(first["geometric_train_rate_mask"], second["geometric_train_rate_mask"])
    assert np.array_equal(first["requested_volume_draws_voxels"], second["requested_volume_draws_voxels"])


def test_log_volume_slope_recovers_linear_log_growth():
    days = np.asarray([0.0, 20.0, 70.0])
    volumes = np.exp(4.0 + 0.012 * days)

    assert np.isclose(hybrid.estimate_log_volume_slope(volumes, days), 0.012)


def test_small_synthetic_case_runs_full_pde_posterior_forecast():
    initial = ball(shape=(16, 16, 16), radius=2.5)
    domain = np.ones(initial.shape, dtype=bool)
    days = np.asarray([0.0, 20.0, 40.0])
    masks = [initial, *fi.forward_masks(initial, domain, 0.04, 0.05, days, 2.0)]
    result = hybrid.hybrid_posterior_predictive_forecast(
        masks, days, 60.0,
        prior={"mean": np.log([0.04, 0.05]), "covariance": np.eye(2) * 0.5}, h=2.0,
        boundary_covariance={"sigma_b_mm": 0.5, "tau_b_mm": 0.2, "ell_mm": 5.0},
        population_log_volume_slope_per_day=0.0, shrinkage_time_constant_days=90.0,
        volume_variance_fit={"intercept": 0.0, "slope_per_day": 0.0}, domain=domain,
        rho_bounds=(0.01, 0.1), D_bounds=(0.01, 0.1), n_coarse=3, n_refine=3,
        n_draws=3, n_anchors=24, seed=5,
    )

    assert result["mask_draws"].shape == (3, *initial.shape)
    assert result["no_pde_ranker_mask_draws"].shape == result["mask_draws"].shape
    assert result["geometric_train_rate_mask"].shape == result["mask_draws"].shape[1:]
    assert result["posterior"]["refined_grid_points"] > result["posterior"]["coarse_grid_points"]
    assert np.isclose(result["posterior"]["posterior_probabilities"].sum(), 1.0)
    assert np.all(result["mask_draws"].reshape(3, -1).sum(axis=1) > 0)
