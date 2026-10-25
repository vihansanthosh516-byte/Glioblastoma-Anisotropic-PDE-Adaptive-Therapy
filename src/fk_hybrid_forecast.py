"""Synthetic-only composition of Track B posterior, volume, and boundary uncertainty.

Callers must provide the joint prior, volume-population slope, shrinkage time constant, and boundary
covariance from outer-training patients only. This module accepts only pre-target history and has no data I/O.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from src import fk_mask_inverse as fi
from src.fk_grid_likelihood import infer_joint_log_parameter_posterior, volume_match_concentration
from src.fk_spatial_uncertainty import correlated_boundary_perturbation
from src.fk_volume_model import sample_positive_volumes, shrink_log_growth_slope


def estimate_log_volume_slope(volumes, days) -> float:
    """Estimate patient log-volume slope by ordinary least squares on history only."""
    volumes = np.asarray(volumes, dtype=float)
    days = np.asarray(days, dtype=float)
    if (volumes.ndim != 1 or days.shape != volumes.shape or len(volumes) < 2
            or not np.isfinite(volumes).all() or not np.isfinite(days).all()
            or np.any(volumes <= 0) or np.any(np.diff(days) <= 0)):
        raise ValueError("history volumes and days must be positive, finite, ordered, and contain at least two scans")
    centered_days = days - days.mean()
    denominator = float(centered_days @ centered_days)
    if denominator <= 0:
        raise ValueError("history scan days must span a positive interval")
    log_volumes = np.log(volumes)
    return float(centered_days @ (log_volumes - log_volumes.mean()) / denominator)


def _perturb_boundary_or_keep(mask, *, boundary_covariance, spacing_mm, seed, n_anchors):
    """Keep the valid pre-draw mask if boundary rasterization collapses a draw."""
    mask = np.asarray(mask, dtype=bool)
    perturbed = correlated_boundary_perturbation(
        mask,
        sigma_b_mm=boundary_covariance["sigma_b_mm"],
        tau_b_mm=boundary_covariance["tau_b_mm"],
        ell_mm=boundary_covariance["ell_mm"],
        spacing_mm=spacing_mm,
        seed=seed,
        n_anchors=n_anchors,
    )
    if not perturbed.any() or perturbed.all():
        return mask.copy(), True
    return perturbed, False


def signed_distance_volume_baseline_draws(last_mask, target_voxels, *, domain=None, spacing_mm=1.0,
                                          boundary_covariance: dict, seed: int = 0,
                                          n_anchors: int = 256,
                                          return_fallback_count: bool = False):
    """Draw matched-volume no-PDE forecasts by ranking voxels with last-mask signed distance."""
    last_mask = np.asarray(last_mask, dtype=bool)
    if last_mask.ndim != 3 or not last_mask.any() or last_mask.all():
        raise ValueError("last_mask must be a nonempty, non-full 3D array")
    domain = np.ones(last_mask.shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != last_mask.shape or np.any(last_mask & ~domain):
        raise ValueError("domain must match and contain last_mask")
    targets = np.asarray(target_voxels)
    if (targets.ndim != 1 or len(targets) < 1 or not np.issubdtype(targets.dtype, np.integer)
            or np.any(targets < 1) or np.any(targets >= int(domain.sum()))):
        raise ValueError("target_voxels must be positive integer counts below domain capacity")
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    signed_distance = (ndimage.distance_transform_edt(last_mask, sampling=spacing)
                       - ndimage.distance_transform_edt(~last_mask, sampling=spacing))
    seeds = np.random.SeedSequence(seed).spawn(len(targets))
    draws = []
    fallback_count = 0
    for target, child_seed in zip(targets, seeds):
        ranked = volume_match_concentration(signed_distance, domain, target_voxels=int(target))
        draw, collapsed = _perturb_boundary_or_keep(
            ranked,
            boundary_covariance=boundary_covariance,
            spacing_mm=spacing,
            seed=int(child_seed.generate_state(1)[0]),
            n_anchors=n_anchors,
        )
        fallback_count += int(collapsed)
        draw &= domain
        draws.append(draw)
    result = np.stack(draws, axis=0)
    return (result, fallback_count) if return_fallback_count else result


def hybrid_posterior_predictive_forecast(masks, days, target_day, *, prior: dict, h: float,
                                         boundary_covariance: dict,
                                         population_log_volume_slope_per_day: float,
                                         shrinkage_time_constant_days: float,
                                         volume_variance_fit: dict, domain=None,
                                         rho_bounds=(1e-5, 1e-1), D_bounds=(1e-5, 10.0),
                                         n_coarse: int = 21, mass_fraction: float = 0.9,
                                         n_refine: int = 11, n_draws: int = 200,
                                         n_anchors: int = 256, seed: int = 0) -> dict:
    """Draw joint PDE-shape, log-volume, and correlated-boundary forecasts from history only.

    The prior and nuisance estimates must come from outer-training patients. The target mask is deliberately
    not an argument. Two masks support the planned first-two-to-third forecast, but one spatial follow-up can
    leave a joint rho-D ridge, so integrate posterior uncertainty and do not interpret the two parameters
    separately. A forecast that reaches the computational-domain edge fails instead of being clipped.
    """
    masks = [np.asarray(mask, dtype=bool) for mask in masks]
    days = np.asarray(days, dtype=float)
    target_day = float(target_day)
    h = float(h)
    if len(masks) < 2 or days.shape != (len(masks),) or not np.isfinite(days).all():
        raise ValueError("need at least two finite history masks and one day per mask")
    if any(mask.ndim != 3 or mask.shape != masks[0].shape or not mask.any() for mask in masks):
        raise ValueError("history masks must be nonempty 3D arrays with one common shape")
    if np.any(np.diff(days) <= 0) or not np.isfinite(target_day) or target_day <= days[-1]:
        raise ValueError("history days must increase and target_day must follow history")
    if not np.isfinite(h) or h <= 0:
        raise ValueError("h must be positive and finite")
    if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 2:
        raise ValueError("n_draws must be an integer of at least two")
    n_draws = int(n_draws)
    domain = np.ones(masks[0].shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != masks[0].shape or any(np.any(mask & ~domain) for mask in masks):
        raise ValueError("domain must match and contain every history mask")

    posterior = infer_joint_log_parameter_posterior(
        masks, days, prior=prior, h=h, covariance_parameters=boundary_covariance, domain=domain,
        rho_bounds=rho_bounds, D_bounds=D_bounds, n_coarse=n_coarse,
        mass_fraction=mass_fraction, n_refine=n_refine, n_draws=n_draws,
        seed=seed, n_anchors=n_anchors,
    )
    history_volumes = np.asarray([mask.sum() for mask in masks], dtype=float)
    patient_slope = estimate_log_volume_slope(history_volumes, days)
    shrunk_slope = shrink_log_growth_slope(
        patient_slope, population_log_volume_slope_per_day,
        observation_gap_days=float(days[-1] - days[0]),
        time_constant_days=shrinkage_time_constant_days,
    )
    horizon = target_day - float(days[-1])
    volume_seed, boundary_seed = np.random.SeedSequence(seed).spawn(2)
    raw_volume_draws = sample_positive_volumes(
        history_volumes[-1], horizon_days=horizon, log_growth_slope=shrunk_slope,
        variance_fit=volume_variance_fit, n_draws=n_draws,
        seed=int(volume_seed.generate_state(1)[0]),
    )
    requested_voxels = np.maximum(1, np.rint(raw_volume_draws).astype(np.int64))
    domain_voxels = int(domain.sum())
    if np.any(requested_voxels >= domain_voxels):
        raise ValueError("forecast volume reaches computational-domain capacity; enlarge the domain")

    seed_sequences = boundary_seed.spawn(n_draws)
    parameter_draws = posterior["parameter_draws_log"]
    masks_out = []
    hybrid_boundary_fallback_count = 0
    pre_boundary_probability = np.zeros(masks[0].shape, dtype=np.float32)
    for draw_index, (log_rho, log_D) in enumerate(parameter_draws):
        rho, diffusion = float(np.exp(log_rho)), float(np.exp(log_D))
        concentration = fi.forward_concentrations(
            masks[-1], domain, rho, diffusion, [float(days[-1]), target_day], h
        )[0]
        sized_mask = volume_match_concentration(
            concentration, domain, target_voxels=int(requested_voxels[draw_index])
        )
        if sized_mask.all() or not sized_mask.any():
            raise ValueError("volume-matched forecast has no usable boundary")
        pre_boundary_probability += sized_mask.astype(np.float32) / n_draws
        boundary_draw, collapsed = _perturb_boundary_or_keep(
            sized_mask,
            boundary_covariance=boundary_covariance,
            spacing_mm=h,
            seed=int(seed_sequences[draw_index].generate_state(1)[0]),
            n_anchors=n_anchors,
        )
        hybrid_boundary_fallback_count += int(collapsed)
        boundary_draw &= domain
        masks_out.append(boundary_draw)

    mask_draws = np.stack(masks_out, axis=0)
    final_volumes = mask_draws.reshape(n_draws, -1).sum(axis=1)
    no_pde_ranker_draws, baseline_boundary_fallback_count = signed_distance_volume_baseline_draws(
        masks[-1], requested_voxels, domain=domain, spacing_mm=h,
        boundary_covariance=boundary_covariance,
        seed=int(boundary_seed.generate_state(1)[0]), n_anchors=n_anchors,
        return_fallback_count=True,
    )
    return {
        "mask_draws": mask_draws,
        "pre_boundary_probability": pre_boundary_probability,
        "no_pde_ranker_mask_draws": no_pde_ranker_draws,
        "hybrid_boundary_draw_fallback_count": int(hybrid_boundary_fallback_count),
        "no_pde_boundary_draw_fallback_count": int(baseline_boundary_fallback_count),
        "parameter_draws_log": parameter_draws,
        "requested_volume_draws_voxels": requested_voxels,
        "final_volume_draws_voxels": final_volumes,
        "patient_log_volume_slope_per_day": patient_slope,
        "shrunk_log_volume_slope_per_day": shrunk_slope,
        "horizon_days": float(horizon),
        "posterior": posterior,
        "boundary_covariance_assumed": dict(boundary_covariance),
        "scope": "synthetic utility; training-only prior and nuisance inputs required; no target mask used",
    }
