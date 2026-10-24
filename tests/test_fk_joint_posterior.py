from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_joint_posterior import (  # noqa: E402
    fit_log_parameter_prior,
    log_gaussian_prior_density,
    make_log_parameter_grid,
    normalize_grid_posterior,
    refine_log_parameter_grid,
    sample_joint_log_parameters,
)


def test_log_grid_uses_requested_limits_and_quadrature_area():
    grid = make_log_parameter_grid(
        rho_bounds=(1e-4, 1e-2), D_bounds=(1e-3, 1.0), n_coarse=5
    )
    expected_area = np.log(1e-2 / 1e-4) * np.log(1.0 / 1e-3)

    assert grid["points"].shape == (25, 2)
    assert np.isclose(np.exp(grid["points"][:, 0]).min(), 1e-4)
    assert np.isclose(np.exp(grid["points"][:, 0]).max(), 1e-2)
    assert np.isclose(np.exp(grid["points"][:, 1]).min(), 1e-3)
    assert np.isclose(np.exp(grid["points"][:, 1]).max(), 1.0)
    assert np.isclose(grid["quadrature_weights"].sum(), expected_area)


def test_empirical_prior_retains_parameter_correlation():
    pairs = np.asarray([[-4.0, -3.0], [-3.5, -2.5], [-3.0, -2.0], [-2.5, -1.5]])
    prior = fit_log_parameter_prior(pairs)

    assert prior["mean"].shape == (2,)
    assert prior["covariance"].shape == (2, 2)
    assert prior["covariance"][0, 1] > 0
    assert np.linalg.eigvalsh(prior["covariance"]).min() > 0
    assert prior["n_training_pairs"] == 4


def test_grid_posterior_normalizes_and_samples_from_joint_points():
    grid = make_log_parameter_grid(
        rho_bounds=(1e-4, 1e-2), D_bounds=(1e-3, 1.0), n_coarse=21
    )
    prior = {"mean": np.asarray([-6.0, -3.0]), "covariance": np.asarray([[1.0, 0.7], [0.7, 1.0]])}
    log_prior = log_gaussian_prior_density(grid["points"], prior)
    center = np.asarray([-5.7, -2.6])
    log_likelihood = -0.5 * np.sum((grid["points"] - center) ** 2 / 0.15, axis=1)
    posterior = normalize_grid_posterior(log_likelihood, log_prior, grid["quadrature_weights"])
    samples = sample_joint_log_parameters(grid["points"], posterior["probabilities"], n_draws=10000, seed=11)

    assert np.isclose(posterior["probabilities"].sum(), 1.0)
    assert np.all(posterior["probabilities"] >= 0)
    assert np.allclose(samples, sample_joint_log_parameters(
        grid["points"], posterior["probabilities"], n_draws=10000, seed=11
    ))
    assert np.linalg.norm(samples.mean(axis=0) - center) < 0.4
    assert np.corrcoef(samples.T)[0, 1] > 0.1


def test_local_refinement_increases_resolution_around_high_mass_region():
    coarse = make_log_parameter_grid(
        rho_bounds=(1e-4, 1e-2), D_bounds=(1e-3, 1.0), n_coarse=21
    )
    center = np.asarray([-4.8, -2.8])
    log_density = -0.5 * np.sum((coarse["points"] - center) ** 2 / 0.05, axis=1)
    refined = refine_log_parameter_grid(coarse, log_density, mass_fraction=0.9, n_refine=11)

    assert len(refined["rho_log_grid"]) > len(coarse["rho_log_grid"])
    assert len(refined["D_log_grid"]) > len(coarse["D_log_grid"])
    assert refined["points"].shape[0] == len(refined["rho_log_grid"]) * len(refined["D_log_grid"])
    assert np.all(refined["quadrature_weights"] > 0)
    assert np.isclose(np.exp(refined["rho_log_grid"]).min(), 1e-4)
    assert np.isclose(np.exp(refined["D_log_grid"]).max(), 1.0)
