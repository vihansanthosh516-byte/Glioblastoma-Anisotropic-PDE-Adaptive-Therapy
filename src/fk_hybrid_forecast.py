"""Synthetic-only composition of Track B posterior, volume, and boundary uncertainty.

Callers must provide a fixed, prespecified joint prior or one derived from outer-training patients. Volume
population slope, shrinkage time constant, predictive variance, and boundary covariance must use outer-training
patients only. This module accepts only pre-target history and has no data I/O.
"""
from __future__ import annotations

import hashlib
import json
import numpy as np
from scipy import ndimage

from src import fk_mask_inverse as fi
from src.fk_grid_likelihood import infer_joint_log_parameter_posterior, volume_match_concentration
from src.fk_spatial_uncertainty import correlated_boundary_perturbation
from src.fk_volume_model import sample_positive_volumes, shrink_log_growth_slope


class ExpectedForecastFailure(ValueError):
    """A named numerical forecast failure that a predeclared fallback may handle."""


class ForecastVolumeCapacityError(ExpectedForecastFailure):
    """Requested history-only volume cannot fit inside the frozen support."""


class ForecastInvalidVolumeError(ExpectedForecastFailure):
    """History-only volume prediction is nonfinite, nonpositive, or malformed."""


class ForecastNoUsableBoundaryError(ExpectedForecastFailure):
    """A PDE draw cannot produce a valid boundary at its requested volume."""


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


def _volume_forecast_signature(masks, days, target_day, *, population_log_volume_slope_per_day,
                               shrinkage_time_constant_days, volume_variance_fit, n_draws, seed,
                               domain_voxels):
    return {
        "history_volumes_voxels": tuple(int(np.asarray(mask, dtype=bool).sum()) for mask in masks),
        "history_days": tuple(float(day) for day in days),
        "target_day": float(target_day),
        "population_log_volume_slope_per_day": float(population_log_volume_slope_per_day),
        "shrinkage_time_constant_days": float(shrinkage_time_constant_days),
        "volume_variance_fit": {str(k): float(v) for k, v in volume_variance_fit.items()},
        "n_draws": int(n_draws),
        "seed": int(seed),
        "domain_voxels": int(domain_voxels),
    }


def _volume_forecast_output_sha256(requested_volume_draws_voxels, geometric_target_voxels,
                                  patient_log_volume_slope_per_day, shrunk_log_volume_slope_per_day,
                                  horizon_days) -> str:
    payload = {
        "requested_volume_draws_voxels": [int(value) for value in requested_volume_draws_voxels],
        "geometric_target_voxels": int(geometric_target_voxels),
        "patient_log_volume_slope_per_day": float(patient_log_volume_slope_per_day),
        "shrunk_log_volume_slope_per_day": float(shrunk_log_volume_slope_per_day),
        "horizon_days": float(horizon_days),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def history_only_volume_forecast(masks, days, target_day, *,
                                 population_log_volume_slope_per_day: float,
                                 shrinkage_time_constant_days: float,
                                 volume_variance_fit: dict, n_draws: int, seed: int,
                                 domain=None) -> dict:
    """Fit and sample the history-only volume forecast without fitting or running the PDE.

    The seed derivation matches ``hybrid_posterior_predictive_forecast``. The returned object can be passed
    back to that function so the PDE candidate and another shape arm use the identical requested volumes.
    """
    masks = [np.asarray(mask, dtype=bool) for mask in masks]
    days = np.asarray(days, dtype=float)
    target_day = float(target_day)
    if len(masks) < 2 or days.shape != (len(masks),) or not np.isfinite(days).all():
        raise ValueError("need at least two finite history masks and one day per mask")
    if any(mask.ndim != 3 or mask.shape != masks[0].shape or not mask.any() for mask in masks):
        raise ValueError("history masks must be nonempty 3D arrays with one common shape")
    if np.any(np.diff(days) <= 0) or not np.isfinite(target_day) or target_day <= days[-1]:
        raise ValueError("history days must increase and target_day must follow history")
    if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 2:
        raise ValueError("n_draws must be an integer of at least two")
    if isinstance(seed, (bool, np.bool_)) or int(seed) != seed:
        raise ValueError("seed must be an integer")
    n_draws, seed = int(n_draws), int(seed)
    domain = np.ones(masks[0].shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != masks[0].shape or any(np.any(mask & ~domain) for mask in masks):
        raise ValueError("domain must match and contain every history mask")
    if domain.sum() < 2:
        raise ForecastVolumeCapacityError("valid support has insufficient capacity for a forecast mask")

    history_volumes = np.asarray([mask.sum() for mask in masks], dtype=float)
    patient_slope = estimate_log_volume_slope(history_volumes, days)
    shrunk_slope = shrink_log_growth_slope(
        patient_slope, population_log_volume_slope_per_day,
        observation_gap_days=float(days[-1] - days[0]),
        time_constant_days=shrinkage_time_constant_days,
    )
    horizon = target_day - float(days[-1])
    volume_seed, _ = np.random.SeedSequence(seed).spawn(2)
    raw_volume_draws = np.asarray(sample_positive_volumes(
        history_volumes[-1], horizon_days=horizon, log_growth_slope=shrunk_slope,
        variance_fit=volume_variance_fit, n_draws=n_draws,
        seed=int(volume_seed.generate_state(1)[0]),
    ), dtype=float)
    if raw_volume_draws.shape != (n_draws,) or not np.isfinite(raw_volume_draws).all() or np.any(raw_volume_draws <= 0):
        raise ForecastInvalidVolumeError("history-only volume forecast produced nonfinite or nonpositive draws")
    domain_voxels = int(domain.sum())
    if np.any(raw_volume_draws >= domain_voxels):
        raise ForecastVolumeCapacityError("forecast volume reaches computational-domain capacity; enlarge the domain")
    requested_voxels = np.maximum(1, np.rint(raw_volume_draws).astype(np.int64))
    if np.any(requested_voxels >= domain_voxels):
        raise ForecastVolumeCapacityError("forecast volume reaches computational-domain capacity; enlarge the domain")

    geometric_volume = history_volumes[-1] * np.exp(float(population_log_volume_slope_per_day) * horizon)
    if not np.isfinite(geometric_volume) or geometric_volume <= 0:
        raise ForecastInvalidVolumeError("geometric train-rate baseline produced an invalid target volume")
    geometric_target_voxels = int(np.rint(geometric_volume))
    if geometric_target_voxels < 1:
        raise ForecastInvalidVolumeError("geometric train-rate baseline produced an invalid target volume")
    if geometric_target_voxels >= domain_voxels:
        raise ForecastVolumeCapacityError("forecast volume reaches computational-domain capacity; enlarge the domain")
    signature = _volume_forecast_signature(
        masks, days, target_day,
        population_log_volume_slope_per_day=population_log_volume_slope_per_day,
        shrinkage_time_constant_days=shrinkage_time_constant_days,
        volume_variance_fit=volume_variance_fit, n_draws=n_draws, seed=seed,
        domain_voxels=domain_voxels,
    )
    signature["output_sha256"] = _volume_forecast_output_sha256(
        requested_voxels, geometric_target_voxels, patient_slope, shrunk_slope, horizon
    )
    return {
        "schema": "track_b_history_only_volume_forecast_v1",
        "signature": signature,
        "requested_volume_draws_voxels": requested_voxels,
        "geometric_target_voxels": geometric_target_voxels,
        "patient_log_volume_slope_per_day": patient_slope,
        "shrunk_log_volume_slope_per_day": float(shrunk_slope),
        "horizon_days": float(horizon),
    }


def _validated_precomputed_volume_forecast(volume_forecast, masks, days, target_day, *,
                                           population_log_volume_slope_per_day,
                                           shrinkage_time_constant_days, volume_variance_fit,
                                           n_draws, seed, domain):
    if not isinstance(volume_forecast, dict) or volume_forecast.get("schema") != "track_b_history_only_volume_forecast_v1":
        raise ValueError("precomputed_volume_forecast must come from history_only_volume_forecast")
    expected = _volume_forecast_signature(
        masks, days, target_day,
        population_log_volume_slope_per_day=population_log_volume_slope_per_day,
        shrinkage_time_constant_days=shrinkage_time_constant_days,
        volume_variance_fit=volume_variance_fit, n_draws=n_draws, seed=seed,
        domain_voxels=int(np.asarray(domain, dtype=bool).sum()),
    )
    signature = volume_forecast.get("signature")
    if not isinstance(signature, dict) or any(signature.get(key) != value for key, value in expected.items()):
        raise ValueError("precomputed_volume_forecast inputs do not match this candidate forecast")
    requested = np.asarray(volume_forecast.get("requested_volume_draws_voxels"))
    geometric = volume_forecast.get("geometric_target_voxels")
    capacity = int(np.asarray(domain, dtype=bool).sum())
    if (requested.shape != (int(n_draws),) or not np.issubdtype(requested.dtype, np.integer)
            or np.any(requested < 1) or np.any(requested >= capacity)
            or isinstance(geometric, bool) or not isinstance(geometric, (int, np.integer))
            or geometric < 1 or geometric >= capacity):
        raise ValueError("precomputed_volume_forecast contains invalid voxel counts")
    try:
        observed_output_sha256 = _volume_forecast_output_sha256(
            requested, geometric,
            volume_forecast["patient_log_volume_slope_per_day"],
            volume_forecast["shrunk_log_volume_slope_per_day"],
            volume_forecast["horizon_days"],
        )
    except (KeyError, TypeError, ValueError, OverflowError) as exc:
        raise ValueError("precomputed_volume_forecast contains invalid output fields") from exc
    if signature.get("output_sha256") != observed_output_sha256:
        raise ValueError("precomputed_volume_forecast output fingerprint does not match its draws")
    return volume_forecast


def _perturb_boundary_or_keep(mask, *, boundary_covariance, spacing_mm, seed, n_anchors,
                             domain=None, target_voxels=None):
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
        domain=domain,
        target_voxels=target_voxels,
    )
    if not perturbed.any() or perturbed.all():
        return mask.copy(), True
    return perturbed, False


def _batched_parameter_concentrations(initial_mask, domain, parameter_pairs, *, horizon_days, spacing_mm,
                                      batch_size=4):
    """Yield PDE fields in small batches sharing D and scalar-equivalent CFL step counts."""
    initial_mask = np.asarray(initial_mask, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    parameters = np.asarray(parameter_pairs, dtype=np.float32)
    horizon_days = float(horizon_days)
    spacing_mm = float(spacing_mm)
    if (initial_mask.ndim != 3 or domain.shape != initial_mask.shape
            or parameters.ndim != 2 or parameters.shape[1] != 2 or not len(parameters)
            or not np.isfinite(parameters).all() or np.any(parameters <= 0)
            or not np.isfinite(horizon_days) or horizon_days <= 0
            or not np.isfinite(spacing_mm) or spacing_mm <= 0):
        raise ValueError("batched PDE inputs must be finite positive parameter pairs on a matching 3D grid")
    if np.any(initial_mask & ~domain):
        raise ValueError("initial PDE mask must be contained in the computational domain")
    if isinstance(batch_size, (bool, np.bool_)) or int(batch_size) != batch_size or batch_size < 1:
        raise ValueError("batch_size must be a positive integer")

    unique_diffusion, diffusion_index = np.unique(parameters[:, 1], return_inverse=True)
    for d_index, diffusion in enumerate(unique_diffusion):
        parameter_indices = np.flatnonzero(diffusion_index == d_index)
        diffusion_tensor = np.zeros(initial_mask.shape + (6,), dtype=np.float32)
        diffusion_tensor[..., :3] = diffusion
        solver = fi.ria.TensorFK(diffusion_tensor, domain, h=spacing_mm)

        indices_by_step_count = {}
        for parameter_index in parameter_indices:
            rho = np.asarray([parameters[parameter_index, 0]], dtype=np.float32)
            step_count = solver.step_count(rho, horizon_days)
            indices_by_step_count.setdefault(step_count, []).append(int(parameter_index))

        for step_count, same_step_indices in sorted(indices_by_step_count.items()):
            for start in range(0, len(same_step_indices), int(batch_size)):
                batch_indices = np.asarray(same_step_indices[start:start + int(batch_size)], dtype=int)
                fields, observed_step_count = solver.run(
                    initial_mask, parameters[batch_indices, 0], horizon_days,
                )
                if observed_step_count != step_count:
                    raise RuntimeError("batched PDE integration changed the scalar-equivalent CFL step count")
                yield batch_indices, fields


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
            domain=domain,
            target_voxels=int(target),
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
                                         n_anchors: int = 256, seed: int = 0,
                                         calibration_only: bool = False,
                                         pde_only_concentration_cutoff: float | None = None,
                                         pde_only_cutoff_grid=None,
                                         precomputed_volume_forecast: dict | None = None) -> dict:
    """Draw joint PDE-shape, log-volume, and correlated-boundary forecasts from history only.

    A fixed prespecified prior or training-only prior is allowed; all fitted nuisance estimates must come from
    outer-training patients. The target mask is deliberately not an argument. Two masks support the planned
    first-two-to-third forecast, but one spatial follow-up can
    leave a joint rho-D ridge, so integrate posterior uncertainty and do not interpret the two parameters
    separately. A forecast that reaches the computational-domain edge fails instead of being clipped.

    ``pde_only_concentration_cutoff`` (PROPOSED LUMIERE PDE-only comparator; default None = unchanged output):
    each raw PDE concentration draw is thresholded at the cutoff inside the domain, before volume matching and
    boundary perturbation, from the same posterior draws and solver run as the candidate.
    ``pde_only_cutoff_grid`` (PROPOSED; default None = unchanged output): the same raw-field thresholding for
    every cutoff in the grid, in both full and calibration-only runs, so training patients' raw cutoff masks
    come from their existing solver runs. Output: inclusion probability and empty-draw count per cutoff.
    ``precomputed_volume_forecast`` (optional): a history-only volume forecast from
    ``history_only_volume_forecast``. Supplying it lets non-PDE shape arms use the exact same volume draws.
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
    if not isinstance(calibration_only, bool):
        raise ValueError("calibration_only must be a boolean")
    n_draws = int(n_draws)
    domain = np.ones(masks[0].shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != masks[0].shape or any(np.any(mask & ~domain) for mask in masks):
        raise ValueError("domain must match and contain every history mask")

    if precomputed_volume_forecast is None:
        volume_forecast = history_only_volume_forecast(
            masks, days, target_day,
            population_log_volume_slope_per_day=population_log_volume_slope_per_day,
            shrinkage_time_constant_days=shrinkage_time_constant_days,
            volume_variance_fit=volume_variance_fit, n_draws=n_draws, seed=seed, domain=domain,
        )
    else:
        volume_forecast = _validated_precomputed_volume_forecast(
            precomputed_volume_forecast, masks, days, target_day,
            population_log_volume_slope_per_day=population_log_volume_slope_per_day,
            shrinkage_time_constant_days=shrinkage_time_constant_days,
            volume_variance_fit=volume_variance_fit, n_draws=n_draws, seed=seed, domain=domain,
        )

    posterior = infer_joint_log_parameter_posterior(
        masks, days, prior=prior, h=h, covariance_parameters=boundary_covariance, domain=domain,
        rho_bounds=rho_bounds, D_bounds=D_bounds, n_coarse=n_coarse,
        mass_fraction=mass_fraction, n_refine=n_refine, n_draws=n_draws,
        seed=seed, n_anchors=n_anchors,
    )
    requested_voxels = np.asarray(volume_forecast["requested_volume_draws_voxels"], dtype=np.int64)
    geometric_target_voxels = int(volume_forecast["geometric_target_voxels"])
    patient_slope = float(volume_forecast["patient_log_volume_slope_per_day"])
    shrunk_slope = float(volume_forecast["shrunk_log_volume_slope_per_day"])
    horizon = float(volume_forecast["horizon_days"])
    domain_voxels = int(domain.sum())
    _, boundary_seed = np.random.SeedSequence(seed).spawn(2)

    seed_sequences = None if calibration_only else boundary_seed.spawn(n_draws)
    parameter_draws_log = np.asarray(posterior["parameter_draws_log"], dtype=float)
    if (parameter_draws_log.shape != (n_draws, 2) or not np.isfinite(parameter_draws_log).all()):
        raise RuntimeError("joint posterior returned malformed log-parameter draws")
    # TensorFK receives float32 rho and D values. Deduplicate on those effective solver values.
    parameter_draws = np.exp(parameter_draws_log).astype(np.float32)
    if not np.isfinite(parameter_draws).all() or np.any(parameter_draws <= 0):
        raise RuntimeError("joint posterior draws are outside the finite positive PDE parameter range")
    unique_parameter_draws, draw_to_parameter = np.unique(parameter_draws, axis=0, return_inverse=True)
    draw_indices_by_parameter = [
        np.flatnonzero(draw_to_parameter == parameter_index)
        for parameter_index in range(len(unique_parameter_draws))
    ]
    masks_out = None if calibration_only else [None] * n_draws
    hybrid_boundary_fallback_count = 0
    pre_boundary_inclusion_counts = np.zeros(masks[0].shape, dtype=np.uint32)
    want_pde_only = pde_only_concentration_cutoff is not None and not calibration_only
    if want_pde_only:
        cutoff = float(pde_only_concentration_cutoff)
        if not 0.0 < cutoff < 1.0:
            raise ValueError("pde_only_concentration_cutoff must lie strictly between 0 and 1")
        pde_only_inclusion_counts = np.zeros(masks[0].shape, dtype=np.uint32)
        pde_only_empty_draws = pde_only_full_draws = 0
    grid_cutoffs = None if pde_only_cutoff_grid is None else tuple(float(c) for c in pde_only_cutoff_grid)
    if grid_cutoffs is not None:
        if not grid_cutoffs or any(not 0.0 < c < 1.0 for c in grid_cutoffs) or len(set(grid_cutoffs)) != len(grid_cutoffs):
            raise ValueError("pde_only_cutoff_grid must hold unique cutoffs strictly between 0 and 1")
        grid_counts = {c: np.zeros(masks[0].shape, dtype=np.uint32) for c in grid_cutoffs}
        grid_empty = {c: 0 for c in grid_cutoffs}
    for parameter_indices, concentration_batch in _batched_parameter_concentrations(
        masks[-1], domain, unique_parameter_draws, horizon_days=horizon, spacing_mm=h,
    ):
        for lane, parameter_index in enumerate(parameter_indices):
            concentration = concentration_batch[lane]
            if want_pde_only:
                raw_mask = (concentration >= cutoff) & domain
                n_same = len(draw_indices_by_parameter[int(parameter_index)])
                pde_only_inclusion_counts += raw_mask.astype(np.uint32) * np.uint32(n_same)
                pde_only_empty_draws += n_same * int(not raw_mask.any())
                pde_only_full_draws += n_same * int(raw_mask.sum() == domain_voxels)
            if grid_cutoffs is not None:
                n_same = len(draw_indices_by_parameter[int(parameter_index)])
                for c in grid_cutoffs:
                    raw_c = (concentration >= c) & domain
                    grid_counts[c] += raw_c.astype(np.uint32) * np.uint32(n_same)
                    grid_empty[c] += n_same * int(not raw_c.any())
            for draw_index in draw_indices_by_parameter[int(parameter_index)]:
                sized_mask = volume_match_concentration(
                    concentration, domain, target_voxels=int(requested_voxels[draw_index])
                )
                if sized_mask.all() or not sized_mask.any():
                    raise ForecastNoUsableBoundaryError("volume-matched forecast has no usable boundary")
                pre_boundary_inclusion_counts += sized_mask
                if calibration_only:
                    continue
                boundary_draw, collapsed = _perturb_boundary_or_keep(
                    sized_mask,
                    boundary_covariance=boundary_covariance,
                    spacing_mm=h,
                    seed=int(seed_sequences[draw_index].generate_state(1)[0]),
                    n_anchors=n_anchors,
                    domain=domain,
                    target_voxels=int(requested_voxels[draw_index]),
                )
                hybrid_boundary_fallback_count += int(collapsed)
                boundary_draw &= domain
                masks_out[draw_index] = boundary_draw

    pre_boundary_probability = pre_boundary_inclusion_counts.astype(np.float32) / n_draws
    grid_out = {} if grid_cutoffs is None else {
        "pde_only_inclusion_probability_by_cutoff": {c: grid_counts[c].astype(np.float32) / n_draws for c in grid_cutoffs},
        "pde_only_empty_draws_by_cutoff": dict(grid_empty),
        "pde_only_horizon_days": float(horizon),
    }
    if calibration_only:
        return {
            **grid_out,
            "pre_boundary_probability": pre_boundary_probability,
            "requested_volume_draws_voxels": requested_voxels,
            "volume_forecast_signature": volume_forecast["signature"],
            "posterior": posterior,
            "calibration_only": True,
            "scope": "inner boundary calibration; volume-matched PDE marginals only; no final mask or comparator draws",
        }

    mask_draws = np.stack(masks_out, axis=0)
    final_volumes = mask_draws.reshape(n_draws, -1).sum(axis=1)
    no_pde_seed, geometric_seed = boundary_seed.spawn(2)
    no_pde_ranker_draws, baseline_boundary_fallback_count = signed_distance_volume_baseline_draws(
        masks[-1], requested_voxels, domain=domain, spacing_mm=h,
        boundary_covariance=boundary_covariance,
        seed=int(no_pde_seed.generate_state(1)[0]), n_anchors=n_anchors,
        return_fallback_count=True,
    )
    geometric_train_rate_draws, geometric_boundary_fallback_count = signed_distance_volume_baseline_draws(
        masks[-1], np.asarray([geometric_target_voxels], dtype=np.int64),
        domain=domain, spacing_mm=h,
        boundary_covariance={"sigma_b_mm": 0.0, "tau_b_mm": 0.0,
                             "ell_mm": float(boundary_covariance["ell_mm"])},
        seed=int(geometric_seed.generate_state(1)[0]), n_anchors=n_anchors,
        return_fallback_count=True,
    )
    no_pde_final_volumes = no_pde_ranker_draws.reshape(n_draws, -1).sum(axis=1)
    geometric_final_volume = int(geometric_train_rate_draws[0].sum())
    if (not np.array_equal(final_volumes, requested_voxels)
            or not np.array_equal(no_pde_final_volumes, requested_voxels)):
        raise RuntimeError("hybrid and matched-volume no-PDE draws must preserve every requested voxel count")
    if geometric_final_volume != geometric_target_voxels:
        raise RuntimeError("geometric train-rate draw did not preserve its requested voxel count")
    pde_only = {}
    if want_pde_only:
        probability = pde_only_inclusion_counts.astype(np.float32) / n_draws
        pde_only = {"pde_only_inclusion_probability": probability,
                    "pde_only_mask": (probability >= 0.5) & domain,  # same 0.5 inclusion rule as the candidate scorer
                    "pde_only_concentration_cutoff": cutoff,
                    "pde_only_empty_draw_count": int(pde_only_empty_draws),
                    "pde_only_full_domain_draw_count": int(pde_only_full_draws)}
    return {
        **pde_only,
        **grid_out,
        "mask_draws": mask_draws,
        "pre_boundary_probability": pre_boundary_probability,
        "no_pde_ranker_mask_draws": no_pde_ranker_draws,
        "no_pde_ranker_final_volume_draws_voxels": no_pde_final_volumes,
        "geometric_train_rate_mask": geometric_train_rate_draws[0],
        "geometric_train_rate_requested_volume_voxels": geometric_target_voxels,
        "geometric_train_rate_boundary_uncertainty": "none; deterministic signed-distance baseline",
        "hybrid_boundary_draw_fallback_count": int(hybrid_boundary_fallback_count),
        "no_pde_boundary_draw_fallback_count": int(baseline_boundary_fallback_count),
        "geometric_train_rate_boundary_draw_fallback_count": int(geometric_boundary_fallback_count),
        "parameter_draws_log": parameter_draws,
        "requested_volume_draws_voxels": requested_voxels,
        "volume_forecast_signature": volume_forecast["signature"],
        "final_volume_draws_voxels": final_volumes,
        "patient_log_volume_slope_per_day": patient_slope,
        "shrunk_log_volume_slope_per_day": shrunk_slope,
        "horizon_days": float(horizon),
        "posterior": posterior,
        "calibration_only": False,
        "boundary_covariance_assumed": dict(boundary_covariance),
        "scope": "synthetic utility; training-only prior and nuisance inputs required; no target mask used",
    }
