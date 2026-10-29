from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_mask_inverse as fi  # noqa: E402
from src.fk_grid_likelihood import (  # noqa: E402
    evaluate_parameter_grid_log_likelihood,
    gaussian_residual_log_likelihood,
    infer_joint_log_parameter_posterior,
    shape_observation_log_likelihood,
    volume_match_concentration,
)
from src import fk_grid_likelihood as fgl  # noqa: E402


def ball(shape=(16, 16, 16), radius=3.0):
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    return sum((axis - center_i) ** 2 for axis, center_i in zip(coords, center)) <= radius**2


def ellipsoid(shape=(24, 24, 24), axes=(6.0, 4.0, 2.5)):
    coords = np.ogrid[tuple(slice(0, n) for n in shape)]
    center = tuple((n - 1) / 2 for n in shape)
    return sum(((axis - center_i) / axis_i) ** 2 for axis, center_i, axis_i in zip(coords, center, axes)) <= 1.0


def test_continuous_forward_fields_threshold_to_the_existing_binary_forecasts():
    mask = ball()
    domain = np.ones(mask.shape, bool)
    days = [0.0, 10.0, 20.0]
    fields = fi.forward_concentrations(mask, domain, 0.04, 0.05, days, 2.0)
    masks = fi.forward_masks(mask, domain, 0.04, 0.05, days, 2.0)

    assert len(fields) == len(masks) == 2
    assert all(field.shape == mask.shape and np.isfinite(field).all() for field in fields)
    assert all(np.array_equal(field > fi.ria.U_VISIBLE, binary) for field, binary in zip(fields, masks))


def test_volume_matcher_selects_exact_top_concentration_count_with_stable_ties():
    concentration = np.ones((3, 3, 3), dtype=float)
    concentration[0, 0, 0] = 3.0
    concentration[0, 0, 1] = 2.0
    domain = np.ones(concentration.shape, bool)
    first = volume_match_concentration(concentration, domain, target_voxels=3)
    second = volume_match_concentration(concentration, domain, target_voxels=3)

    assert int(first.sum()) == 3
    assert np.array_equal(first, second)
    assert first[0, 0, 0] and first[0, 0, 1] and first[0, 0, 2]


def test_gaussian_residual_likelihood_matches_closed_form_for_correlated_noise():
    residual = np.asarray([0.0, 0.0])
    covariance = np.asarray([[2.0, 0.5], [0.5, 1.0]])
    expected = -0.5 * (2 * np.log(2 * np.pi) + np.linalg.slogdet(covariance)[1])
    assert np.isclose(gaussian_residual_log_likelihood(residual, covariance), expected)


def test_volume_matched_shape_likelihood_prefers_the_generating_pde_parameters():
    initial = ellipsoid()
    domain = np.ones(initial.shape, bool)
    days = [0.0, 60.0, 120.0]
    true_rho, true_D = 0.04, 0.05
    target = fi.forward_masks(initial, domain, true_rho, true_D, days, 2.0)[-1]
    true_field = fi.forward_concentrations(initial, domain, true_rho, true_D, days, 2.0)[-1]
    wrong_field = fi.forward_concentrations(initial, domain, 0.04, 0.15, days, 2.0)[-1]
    noise = {"sigma_b_mm": 1.0, "tau_b_mm": 0.2, "ell_mm": 5.0}

    true_fit = shape_observation_log_likelihood(target, true_field, domain, spacing_mm=2.0,
                                                covariance_parameters=noise, n_anchors=64)
    wrong_fit = shape_observation_log_likelihood(target, wrong_field, domain, spacing_mm=2.0,
                                                 covariance_parameters=noise, n_anchors=64)

    assert true_fit["target_voxels"] == int(target.sum())
    assert true_fit["n_anchors"] == 64
    assert np.allclose(true_fit["residuals_mm"], 0.0)
    assert true_fit["log_likelihood"] > wrong_fit["log_likelihood"]


def test_parameter_grid_likelihood_is_aligned_to_joint_log_grid_nodes():
    initial = ellipsoid()
    domain = np.ones(initial.shape, bool)
    days = [0.0, 60.0, 120.0]
    masks = [initial, *fi.forward_masks(initial, domain, 0.04, 0.05, days, 2.0)]
    points = np.asarray([[np.log(0.04), np.log(0.05)], [np.log(0.04), np.log(0.15)]])
    values = evaluate_parameter_grid_log_likelihood(
        masks, days, points, h=2.0, domain=domain,
        covariance_parameters={"sigma_b_mm": 1.0, "tau_b_mm": 0.2, "ell_mm": 5.0}, n_anchors=64,
    )

    assert values.shape == (2,)
    assert np.isfinite(values).all()
    assert values[0] > values[1]


def test_batched_grid_likelihood_matches_serial_pde_runs():
    initial = ellipsoid()
    domain = np.ones(initial.shape, bool)
    days = [0.0, 30.0, 60.0]
    masks = [initial, *fi.forward_masks(initial, domain, 0.04, 0.05, days, 2.0)]
    points = np.asarray([
        [np.log(0.04), np.log(0.05)],
        [np.log(0.06), np.log(0.05)],
        [np.log(0.04), np.log(0.15)],
    ])
    noise = {"sigma_b_mm": 1.0, "tau_b_mm": 0.2, "ell_mm": 5.0}
    observed = [
        fgl._prepare_shape_observation(
            target, domain, spacing_mm=2.0, covariance_parameters=noise, n_anchors=32
        )
        for target in masks[1:]
    ]
    expected = []
    for log_rho, log_D in points:
        fields = fi.forward_concentrations(initial, domain, np.exp(log_rho), np.exp(log_D), days, 2.0)
        expected.append(sum(
            fgl._score_prepared_shape_observation(obs, field)["log_likelihood"]
            for obs, field in zip(observed, fields)
        ))

    actual = evaluate_parameter_grid_log_likelihood(
        masks, days, points, h=2.0, domain=domain,
        covariance_parameters=noise, n_anchors=32,
    )

    assert np.allclose(actual, expected, rtol=1e-10, atol=1e-8)


def test_parameter_grid_reuses_solver_across_rho_batches(monkeypatch):
    initial = ellipsoid()
    domain = np.ones(initial.shape, bool)
    days = [0.0, 10.0, 20.0]
    masks = [initial, *fi.forward_masks(initial, domain, 0.04, 0.05, days, 2.0)]
    log_D = np.log(0.05)
    points = np.asarray([
        [np.log(rho), log_D] for rho in (0.02, 0.03, 0.04, 0.05, 0.06)
    ])
    original_solver = fgl.fi.ria.TensorFK
    constructed = []

    class CountingTensorFK(original_solver):
        def __init__(self, *args, **kwargs):
            constructed.append(1)
            super().__init__(*args, **kwargs)

    monkeypatch.setattr(fgl.fi.ria, "TensorFK", CountingTensorFK)
    values = evaluate_parameter_grid_log_likelihood(
        masks, days, points, h=2.0, domain=domain,
        covariance_parameters={"sigma_b_mm": 1.0, "tau_b_mm": 0.2, "ell_mm": 5.0},
        n_anchors=32,
    )

    assert len(constructed) == 1
    assert values.shape == (5,)
    assert np.isfinite(values).all()


def test_parameter_grid_reuses_reference_covariance_factor(monkeypatch):
    initial = ellipsoid()
    domain = np.ones(initial.shape, bool)
    days = [0.0, 30.0, 60.0]
    masks = [initial, *fi.forward_masks(initial, domain, 0.04, 0.05, days, 2.0)]
    original = fgl.exponential_covariance_matrix
    calls = []

    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)

    monkeypatch.setattr(fgl, "exponential_covariance_matrix", counted)
    points = np.asarray([[np.log(0.04), np.log(0.05)], [np.log(0.04), np.log(0.15)]])
    evaluate_parameter_grid_log_likelihood(
        masks, days, points, h=2.0, domain=domain,
        covariance_parameters={"sigma_b_mm": 1.0, "tau_b_mm": 0.2, "ell_mm": 5.0}, n_anchors=32,
    )

    assert len(calls) == len(masks) - 1


def test_joint_posterior_wrapper_normalizes_refined_grid_and_draws(monkeypatch):
    center = np.log([0.01, 0.02])

    def likelihood(masks, days, points, **kwargs):
        return -0.5 * np.sum(((points - center) / 0.15) ** 2, axis=1)

    monkeypatch.setattr(fgl, "evaluate_parameter_grid_log_likelihood", likelihood)
    result = infer_joint_log_parameter_posterior(
        [np.ones((2, 2, 2), bool)] * 3, [0.0, 1.0, 2.0],
        prior={"mean": center, "covariance": np.eye(2) * 0.2}, h=2.0,
        covariance_parameters={"sigma_b_mm": 1.0, "tau_b_mm": 0.2, "ell_mm": 5.0},
        rho_bounds=(1e-3, 1e-1), D_bounds=(1e-3, 1e-1), n_coarse=5,
        n_refine=3, n_draws=40, seed=31, n_anchors=8,
    )

    assert np.isclose(result["posterior_probabilities"].sum(), 1.0)
    assert result["parameter_draws_log"].shape == (40, 2)
    assert np.linalg.norm(result["posterior_mean_log_parameters"] - center) < 0.5
    assert result["n_likelihood_evaluations"] == (
        result["coarse_grid_points"] + result["refined_grid_points"]
        - result["n_reused_coarse_likelihoods"]
    )
    assert "target excluded" in result["scope"]
