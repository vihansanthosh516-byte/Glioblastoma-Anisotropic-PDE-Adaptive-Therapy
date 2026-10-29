"""PDE-derived shape likelihood for joint Track B growth-parameter posteriors.

Each candidate PDE field is thresholded by matching its ranked voxels to the observed mask volume, then its
boundary displacement is scored at target-mask surface anchors under a supplied spatial Gaussian covariance.
The supplied segmentation is treated as a reference label, not as known ground truth or measured scanner noise.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from src import fk_mask_inverse as fi
from src.fk_joint_posterior import (log_parameter_prior_density,
                                    make_log_parameter_grid, normalize_grid_posterior,
                                    refine_log_parameter_grid, sample_joint_log_parameters)
from src.fk_spatial_uncertainty import deterministic_surface_anchors, exponential_covariance_matrix


def volume_match_concentration(concentration, domain, *, target_voxels: int) -> np.ndarray:
    """Select the highest-concentration voxels inside domain, with stable flat-index ties."""
    concentration = np.asarray(concentration, dtype=float)
    domain = np.asarray(domain, dtype=bool)
    if concentration.ndim != 3 or domain.shape != concentration.shape or not np.isfinite(concentration).all():
        raise ValueError("concentration and domain must be matching finite 3D arrays")
    if isinstance(target_voxels, (bool, np.bool_)) or int(target_voxels) != target_voxels:
        raise ValueError("target_voxels must be an integer")
    target_voxels = int(target_voxels)
    domain_indices = np.flatnonzero(domain)
    if target_voxels < 0 or target_voxels > len(domain_indices):
        raise ValueError("target_voxels must be between zero and domain size")
    result = np.zeros(concentration.shape, dtype=bool)
    if target_voxels == 0:
        return result
    values = concentration.ravel()[domain_indices]
    threshold_index = values.size - target_voxels
    threshold = np.partition(values, threshold_index)[threshold_index]
    above_threshold = values > threshold
    result.ravel()[domain_indices[above_threshold]] = True
    n_tied = target_voxels - int(above_threshold.sum())
    if n_tied:
        # domain_indices are sorted; preserve the prior lexsort's stable flat-index tie order.
        tied_indices = domain_indices[np.flatnonzero(values == threshold)]
        result.ravel()[tied_indices[:n_tied]] = True
    return result


def gaussian_residual_log_likelihood(residuals, covariance) -> float:
    """Evaluate the normalized zero-mean multivariate Gaussian log likelihood."""
    residuals = np.asarray(residuals, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if residuals.ndim != 1 or len(residuals) == 0 or not np.isfinite(residuals).all():
        raise ValueError("residuals must be a nonempty finite vector")
    if (covariance.shape != (len(residuals), len(residuals)) or not np.isfinite(covariance).all()
            or not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12)):
        raise ValueError("covariance must be finite, symmetric, and aligned to residuals")
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("covariance must be positive definite") from exc
    return _gaussian_log_likelihood_from_factor(residuals, factor)


def _gaussian_log_likelihood_from_factor(residuals, factor) -> float:
    residuals = np.asarray(residuals, dtype=float)
    factor = np.asarray(factor, dtype=float)
    if (residuals.ndim != 1 or factor.shape != (len(residuals), len(residuals))
            or not np.isfinite(residuals).all() or not np.isfinite(factor).all()
            or np.any(np.diag(factor) <= 0)):
        raise ValueError("residuals and Cholesky factor must be finite and aligned")
    whitened = np.linalg.solve(factor, residuals)
    log_determinant = 2.0 * np.log(np.diag(factor)).sum()
    return float(-0.5 * (len(residuals) * np.log(2.0 * np.pi) + log_determinant + whitened @ whitened))


def _signed_distance(mask, spacing):
    mask = np.asarray(mask, dtype=bool)
    return ndimage.distance_transform_edt(mask, sampling=spacing) - ndimage.distance_transform_edt(~mask, sampling=spacing)


def _prepare_shape_observation(reference_mask, domain, *, spacing_mm, covariance_parameters, n_anchors):
    reference = np.asarray(reference_mask, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if reference.ndim != 3 or domain.shape != reference.shape or not reference.any() or reference.all():
        raise ValueError("reference and domain must be valid matching 3D arrays")
    if np.any(reference & ~domain):
        raise ValueError("reference mask must be contained in domain")
    if not isinstance(covariance_parameters, dict):
        raise ValueError("covariance_parameters must be a mapping")
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    anchors = deterministic_surface_anchors(reference, spacing_mm=spacing, n_anchors=n_anchors)
    reference_distance = _signed_distance(reference, spacing)
    anchor_indices = (anchors / spacing).T
    observed = ndimage.map_coordinates(reference_distance, anchor_indices, order=1, mode="nearest")
    covariance = exponential_covariance_matrix(
        anchors,
        sigma_b_mm=covariance_parameters["sigma_b_mm"],
        tau_b_mm=covariance_parameters["tau_b_mm"],
        ell_mm=covariance_parameters["ell_mm"],
    )
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("boundary covariance must be positive definite") from exc
    return {
        "reference": reference,
        "domain": domain,
        "spacing": spacing,
        "anchor_indices": anchor_indices,
        "observed_distance_mm": observed,
        "factor": factor,
        "target_voxels": int(reference.sum()),
        "n_anchors": int(len(anchors)),
        "covariance_parameters": dict(covariance_parameters),
    }


def _score_prepared_shape_observation(prepared, concentration) -> dict:
    concentration = np.asarray(concentration, dtype=float)
    if concentration.shape != prepared["reference"].shape or not np.isfinite(concentration).all():
        raise ValueError("concentration must be finite and match reference shape")
    prediction = volume_match_concentration(
        concentration, prepared["domain"], target_voxels=prepared["target_voxels"]
    )
    prediction_distance = _signed_distance(prediction, prepared["spacing"])
    predicted = ndimage.map_coordinates(
        prediction_distance, prepared["anchor_indices"], order=1, mode="nearest"
    )
    residuals = np.asarray(predicted - prepared["observed_distance_mm"], dtype=float)
    return {
        "log_likelihood": _gaussian_log_likelihood_from_factor(residuals, prepared["factor"]),
        "residuals_mm": residuals,
        "predicted_mask": prediction,
        "target_voxels": prepared["target_voxels"],
        "n_anchors": prepared["n_anchors"],
        "covariance_parameters": prepared["covariance_parameters"],
    }


def shape_observation_log_likelihood(reference_mask, concentration, domain, *, spacing_mm,
                                     covariance_parameters: dict, n_anchors: int = 256) -> dict:
    """Score one predicted PDE field against one supplied reference mask.

    The signed-distance residual is predicted minus reference at deterministic inside-surface voxel centers.
    Covariance parameters are supplied by the outer-training process; they are never estimated from this target.
    """
    prepared = _prepare_shape_observation(
        reference_mask, domain, spacing_mm=spacing_mm,
        covariance_parameters=covariance_parameters, n_anchors=n_anchors,
    )
    return _score_prepared_shape_observation(prepared, concentration)


def evaluate_parameter_grid_log_likelihood(masks, days, log_parameter_points, *, h: float, domain=None,
                                           covariance_parameters: dict, n_anchors: int = 256) -> np.ndarray:
    """Evaluate a joint log(rho)-log(D) grid using training-series shape residuals only.

    The first mask seeds the PDE. Each later concentration is volume-matched to its corresponding supplied
    training mask before the spatial likelihood is evaluated. With only two masks, the single follow-up can leave
    a rho-D ridge; the joint posterior must retain that uncertainty and must not be presented as separate
    patient-specific parameter identification. Use only outer-training masks here.
    """
    masks = [np.asarray(mask, dtype=bool) for mask in masks]
    days = np.asarray(days, dtype=float)
    points = np.asarray(log_parameter_points, dtype=float)
    if len(masks) < 2 or days.shape != (len(masks),) or np.any(np.diff(days) <= 0) or not np.isfinite(days).all():
        raise ValueError("need at least two ordered masks and one finite day per mask")
    if any(mask.ndim != 3 or mask.shape != masks[0].shape or not mask.any() for mask in masks):
        raise ValueError("all training masks must be nonempty 3D arrays with one shape")
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("log_parameter_points must be finite with shape (n, 2)")
    if not np.isfinite(h) or h <= 0:
        raise ValueError("h must be positive and finite")
    domain = np.ones(masks[0].shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != masks[0].shape or any(np.any(mask & ~domain) for mask in masks):
        raise ValueError("domain must match and contain every training mask")
    observations = [
        _prepare_shape_observation(
            target, domain, spacing_mm=h, covariance_parameters=covariance_parameters,
            n_anchors=n_anchors,
        )
        for target in masks[1:]
    ]

    parameter_values = np.exp(points)
    if not np.isfinite(parameter_values).all() or np.any(parameter_values <= 0):
        raise ValueError("log parameter grid contains values outside the finite positive PDE range")
    values = np.empty(len(points), dtype=float)
    unique_D, D_index = np.unique(points[:, 1], return_inverse=True)
    batch_size = 4
    for diffusion_index, log_D in enumerate(unique_D):
        member_indices = np.flatnonzero(D_index == diffusion_index)
        diffusion = float(np.exp(log_D))
        tensor = np.zeros(masks[0].shape + (6,), dtype=np.float32)
        tensor[..., :3] = diffusion
        # The tensor operator depends on D and the domain, not on rho. Reuse its
        # precomputed face weights across all rho chunks for this diffusion value.
        solver = fi.ria.TensorFK(tensor, domain, h=h)
        for batch_start in range(0, len(member_indices), batch_size):
            indices = member_indices[batch_start:batch_start + batch_size]
            rhos = parameter_values[indices, 0]
            state = masks[0].astype(np.float32)
            previous_day = float(days[0])
            batch_values = np.zeros(len(indices), dtype=float)
            for observation, day in zip(observations, days[1:]):
                state, _ = solver.run(state, rhos, float(day - previous_day))
                for member_index, concentration in enumerate(state):
                    batch_values[member_index] += _score_prepared_shape_observation(
                        observation, concentration
                    )["log_likelihood"]
                previous_day = float(day)
            values[indices] = batch_values
    return values


def infer_joint_log_parameter_posterior(masks, days, *, prior: dict, h: float,
                                        covariance_parameters: dict, domain=None,
                                        rho_bounds=(1e-5, 1e-1), D_bounds=(1e-5, 10.0),
                                        n_coarse: int = 21, mass_fraction: float = 0.9,
                                        n_refine: int = 11, n_draws: int = 200,
                                        seed: int = 0, n_anchors: int = 256) -> dict:
    """Fit one training-history posterior with a fixed or outer-training joint prior.

    A data-derived prior and the spatial covariance must use outer-training patients only. A fixed, prespecified
    prior such as ``{"distribution": "log_uniform"}`` is also supported. ``masks`` must be the history available
    before the forecast target; no forecast target is accepted by this function.
    """
    coarse = make_log_parameter_grid(rho_bounds=rho_bounds, D_bounds=D_bounds, n_coarse=n_coarse)
    coarse_log_likelihood = evaluate_parameter_grid_log_likelihood(
        masks, days, coarse["points"], h=h, domain=domain,
        covariance_parameters=covariance_parameters, n_anchors=n_anchors,
    )
    coarse_log_prior = log_parameter_prior_density(
        coarse["points"], prior, rho_bounds=rho_bounds, D_bounds=D_bounds
    )
    coarse_log_posterior = coarse_log_likelihood + coarse_log_prior
    coarse_result = normalize_grid_posterior(
        coarse_log_likelihood, coarse_log_prior, coarse["quadrature_weights"]
    )
    refined = refine_log_parameter_grid(
        coarse, coarse_log_posterior, mass_fraction=mass_fraction, n_refine=n_refine
    )
    source_indices = refined["coarse_source_indices"]
    refined_log_likelihood = np.empty(len(refined["points"]), dtype=float)
    reused = source_indices >= 0
    refined_log_likelihood[reused] = coarse_log_likelihood[source_indices[reused]]
    new_points = ~reused
    if np.any(new_points):
        refined_log_likelihood[new_points] = evaluate_parameter_grid_log_likelihood(
            masks, days, refined["points"][new_points], h=h, domain=domain,
            covariance_parameters=covariance_parameters, n_anchors=n_anchors,
        )
    refined_log_prior = log_parameter_prior_density(
        refined["points"], prior, rho_bounds=rho_bounds, D_bounds=D_bounds
    )
    posterior = normalize_grid_posterior(
        refined_log_likelihood, refined_log_prior, refined["quadrature_weights"]
    )
    probabilities = posterior["probabilities"]
    mean = np.average(refined["points"], axis=0, weights=probabilities)
    centered = refined["points"] - mean
    covariance = (centered * probabilities[:, None]).T @ centered
    draws = sample_joint_log_parameters(refined["points"], probabilities, n_draws=n_draws, seed=seed)
    rho_log_bounds = np.log(np.asarray(rho_bounds, dtype=float))
    D_log_bounds = np.log(np.asarray(D_bounds, dtype=float))
    boundary_mass = {
        "rho_lower": float(probabilities[np.isclose(refined["points"][:, 0], rho_log_bounds[0])].sum()),
        "rho_upper": float(probabilities[np.isclose(refined["points"][:, 0], rho_log_bounds[1])].sum()),
        "D_lower": float(probabilities[np.isclose(refined["points"][:, 1], D_log_bounds[0])].sum()),
        "D_upper": float(probabilities[np.isclose(refined["points"][:, 1], D_log_bounds[1])].sum()),
    }
    return {
        "grid": refined,
        "posterior_probabilities": probabilities,
        "posterior_mean_log_parameters": mean,
        "posterior_covariance_log_parameters": covariance,
        "parameter_draws_log": draws,
        "coarse_log_evidence": coarse_result["log_evidence"],
        "refined_log_evidence": posterior["log_evidence"],
        "coarse_grid_points": int(len(coarse["points"])),
        "refined_grid_points": int(len(refined["points"])),
        "n_likelihood_evaluations": int(len(coarse["points"]) + np.count_nonzero(new_points)),
        "n_reused_coarse_likelihoods": int(np.count_nonzero(reused)),
        "posterior_mass_at_parameter_bounds": boundary_mass,
        "parameter_order": ("log_rho", "log_D"),
        "prior_distribution": str(prior.get("distribution", "gaussian")),
        "scope": "history only; target excluded; prior is fixed or development-only and covariance is outer-training",
    }
