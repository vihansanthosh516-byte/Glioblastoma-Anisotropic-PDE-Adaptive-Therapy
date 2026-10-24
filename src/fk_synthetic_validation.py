"""Synthetic-only Track B forecast checks; no patient data are read by this module.

The boundary perturbation is a stated simulation assumption: a spatially correlated random displacement of
the signed-distance boundary. It is not a measured segmentation-error distribution. The forecast comparison
keeps the later target out of parameter fitting and compares a PDE fit with persistence and a volume-trend
shape-preserving baseline.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage

from src import fk_mask_inverse as fi
from src.forecast_metrics import all_metrics, dice


def _spacing3(spacing_mm) -> np.ndarray:
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    return spacing


def perturb_boundary(mask, *, sigma_mm: float, correlation_mm: float, spacing_mm=1.0, seed: int = 0) -> np.ndarray:
    """Perturb a 3D binary boundary with a smooth signed-distance noise field.

    ``sigma_mm`` sets the global standard deviation of the smoothed displacement field; ``correlation_mm``
    sets its approximate spatial correlation scale. These are simulation controls, not empirical estimates.
    """
    mask = np.asarray(mask, dtype=bool)
    spacing = _spacing3(spacing_mm)
    if mask.ndim != 3:
        raise ValueError("mask must be 3D")
    if not np.isfinite(sigma_mm) or sigma_mm < 0 or not np.isfinite(correlation_mm) or correlation_mm <= 0:
        raise ValueError("sigma_mm must be nonnegative and correlation_mm must be positive")
    if not mask.any() or mask.all():
        raise ValueError("mask must contain both inside and outside voxels")
    if sigma_mm == 0:
        return mask.copy()

    signed_distance = (ndimage.distance_transform_edt(mask, sampling=spacing)
                       - ndimage.distance_transform_edt(~mask, sampling=spacing))
    rng = np.random.default_rng(seed)
    white = rng.normal(size=mask.shape)
    field = ndimage.gaussian_filter(white, sigma=correlation_mm / spacing, mode="reflect")
    field_sd = float(field.std())
    if field_sd == 0:
        return mask.copy()
    return (signed_distance + field * (float(sigma_mm) / field_sd)) > 0


def project_mask_volume(mask, target_voxels: int, *, domain=None, spacing_mm=1.0) -> np.ndarray:
    """Make the smallest-boundary-change mask with exactly ``target_voxels`` voxels.

    Growth adds domain voxels nearest to the current mask. Shrinkage retains the deepest interior voxels.
    Stable flat-index tie-breaking makes the projection deterministic.
    """
    mask = np.asarray(mask, dtype=bool)
    spacing = _spacing3(spacing_mm)
    if mask.ndim != 3:
        raise ValueError("mask must be 3D")
    domain = np.ones(mask.shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != mask.shape or np.any(mask & ~domain):
        raise ValueError("domain must match mask shape and contain the mask")
    if isinstance(target_voxels, (bool, np.bool_)) or int(target_voxels) != target_voxels:
        raise ValueError("target_voxels must be an integer")
    target_voxels = int(target_voxels)
    if target_voxels < 0 or target_voxels > int(domain.sum()):
        raise ValueError("target_voxels must be between zero and the domain size")
    current = int(mask.sum())
    if target_voxels == current:
        return mask.copy()
    if target_voxels == 0:
        return np.zeros_like(mask)
    if current == 0:
        raise ValueError("cannot expand an empty mask by nearest-boundary projection")

    flat_indices = np.arange(mask.size)
    if target_voxels > current:
        candidates = domain & ~mask
        distance = ndimage.distance_transform_edt(~mask, sampling=spacing)
        indices = flat_indices[candidates.ravel()]
        order = np.lexsort((indices, distance.ravel()[indices]))
        selected = indices[order[:target_voxels - current]]
        result = mask.copy()
        result.ravel()[selected] = True
        return result

    depth = ndimage.distance_transform_edt(mask, sampling=spacing)
    indices = flat_indices[mask.ravel()]
    order = np.lexsort((indices, -depth.ravel()[indices]))
    selected = indices[order[:target_voxels]]
    result = np.zeros_like(mask)
    result.ravel()[selected] = True
    return result


def ensemble_boundary_disagreement(forecasts) -> float:
    """Mean normalized Bernoulli disagreement across a forecast-mask ensemble.

    The score is in [0, 1] and uses a one-voxel dilation of the forecast union as its support, so empty
    background does not dilute boundary uncertainty. It is a forecast-spread diagnostic, not an error rate.
    """
    masks = [np.asarray(mask, dtype=bool) for mask in forecasts]
    if len(masks) < 2:
        raise ValueError("need at least two forecast masks")
    if masks[0].ndim != 3 or any(mask.shape != masks[0].shape for mask in masks):
        raise ValueError("forecast masks must share one 3D shape")
    probability = np.mean(np.stack(masks, axis=0), axis=0)
    support = ndimage.binary_dilation(np.any(probability > 0, axis=0), iterations=1)
    if not support.any():
        return 0.0
    disagreement = 4.0 * probability * (1.0 - probability)
    return float(np.clip(np.mean(disagreement[support]), 0.0, 1.0))


def bootstrap_forecast_ensemble(training_masks, training_days, target_day: float, *, h: float,
                                sigma_mm: float, correlation_mm: float, n_draws: int = 20,
                                seed: int = 0, domain=None) -> dict:
    """Refit forecasts after independent assumed boundary perturbations of training masks.

    The target is represented only by its day; no held-out mask or score enters fitting or forecast spread.
    Each draw perturbs every training mask with independent seeded boundary noise, refits the isotropic PDE,
    and forecasts from that draw's first mask. The noise law is an explicit simulation assumption.
    """
    masks = [np.asarray(mask, dtype=bool) for mask in training_masks]
    days = np.asarray(training_days, dtype=float)
    if len(masks) < 3 or days.ndim != 1 or len(days) != len(masks):
        raise ValueError("need at least three training masks and one day per mask")
    if any(mask.ndim != 3 or mask.shape != masks[0].shape for mask in masks):
        raise ValueError("training masks must share one 3D shape")
    if not np.isfinite(days).all() or np.any(np.diff(days) <= 0):
        raise ValueError("training days must be finite and strictly increasing")
    if not np.isfinite(target_day) or target_day <= days[-1]:
        raise ValueError("target_day must follow the final training day")
    if not np.isfinite(h) or h <= 0:
        raise ValueError("h must be positive and finite")
    if not np.isfinite(sigma_mm) or sigma_mm < 0 or not np.isfinite(correlation_mm) or correlation_mm <= 0:
        raise ValueError("sigma_mm must be nonnegative and correlation_mm must be positive")
    if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 2:
        raise ValueError("n_draws must be an integer of at least two")
    n_draws = int(n_draws)
    domain = np.ones(masks[0].shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != masks[0].shape or any(np.any(mask & ~domain) for mask in masks):
        raise ValueError("domain must match and contain all training masks")

    seed_sequences = np.random.SeedSequence(seed).spawn(n_draws * len(masks))
    forecasts, fits = [], []
    for draw in range(n_draws):
        noisy_masks = []
        for scan, mask in enumerate(masks):
            noise_seed = int(seed_sequences[draw * len(masks) + scan].generate_state(1)[0])
            noisy = perturb_boundary(
                mask, sigma_mm=sigma_mm, correlation_mm=correlation_mm, spacing_mm=h, seed=noise_seed
            ) & domain
            if not noisy.any():
                raise ValueError("boundary perturbation emptied a training mask inside the domain")
            noisy_masks.append(noisy)
        fit = fi.fit_fk_masks(noisy_masks, days, h, domain=domain)
        prediction = fi.forward_masks(
            noisy_masks[0], domain, fit["rho"], fit["D"], [days[0], float(target_day)], h
        )[0]
        fits.append(fit)
        forecasts.append(prediction)
    return {"fits": fits, "forecasts": forecasts,
            "disagreement": ensemble_boundary_disagreement(forecasts),
            "assumed_sigma_mm": float(sigma_mm), "assumed_correlation_mm": float(correlation_mm),
            "seed": int(seed)}


def conformal_error_upper_bounds(calibration_scores, calibration_errors, prediction_scores, *,
                                 alpha: float = 0.1, score_floor: float = 1e-6) -> dict:
    """Calibrate upper bounds for [0, 1] forecast error from separate calibration cases.

    Uses split-conformal ratios error / (score + score_floor). Under exchangeable calibration and future
    cases, this gives marginal finite-sample coverage; it does not give conditional coverage or handle shift.
    If the calibration set is too small for the requested alpha, the returned bound is the uninformative 1.
    Scores must be fixed before observing their held-out errors. ``score_floor`` is fixed in normalized units.
    """
    scores = np.asarray(calibration_scores, dtype=float)
    errors = np.asarray(calibration_errors, dtype=float)
    future = np.asarray(prediction_scores, dtype=float)
    if scores.ndim != 1 or errors.ndim != 1 or len(scores) != len(errors) or len(scores) == 0:
        raise ValueError("calibration scores and errors must be nonempty 1D arrays of equal length")
    if future.ndim != 1:
        raise ValueError("prediction_scores must be a 1D array")
    if (not np.isfinite(scores).all() or not np.isfinite(errors).all() or not np.isfinite(future).all()
            or np.any((scores < 0) | (scores > 1)) or np.any((errors < 0) | (errors > 1))
            or np.any((future < 0) | (future > 1))):
        raise ValueError("scores and errors must be finite values in [0, 1]")
    if not np.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be strictly between zero and one")
    if not np.isfinite(score_floor) or score_floor <= 0:
        raise ValueError("score_floor must be positive and finite")

    n = len(errors)
    rank = int(np.ceil((n + 1) * (1.0 - alpha)))
    if rank > n:
        quantile = float("inf")
    else:
        ratios = errors / (scores + float(score_floor))
        quantile = float(np.partition(ratios, rank - 1)[rank - 1])
    bounds = np.minimum(1.0, quantile * (future + float(score_floor)))
    return {"quantile": quantile, "bounds": bounds, "alpha": float(alpha),
            "calibration_n": int(n), "score_floor": float(score_floor),
            "coverage_type": "marginal under exchangeability"}


def compare_next_scan(training_masks, training_days, target_day: float, target_mask, *, h: float, domain=None) -> dict:
    """Fit using training scans only, then compare PDE, persistence, and volume-trend forecasts.

    Requires at least two follow-up scans for a two-parameter fit. The target scan is used only for scoring.
    Returned masks are included to support transparent review of each forecast.
    """
    masks = [np.asarray(mask, dtype=bool) for mask in training_masks]
    days = np.asarray(training_days, dtype=float)
    target = np.asarray(target_mask, dtype=bool)
    if len(masks) < 3 or days.ndim != 1 or len(days) != len(masks):
        raise ValueError("need at least three training masks and one day per mask")
    if any(mask.shape != masks[0].shape for mask in masks) or target.shape != masks[0].shape or masks[0].ndim != 3:
        raise ValueError("training masks and target must share one 3D shape")
    if not np.isfinite(days).all() or np.any(np.diff(days) <= 0) or not np.isfinite(target_day) or target_day <= days[-1]:
        raise ValueError("training days must increase and target_day must follow the final training day")
    domain = np.ones(masks[0].shape, dtype=bool) if domain is None else np.asarray(domain, dtype=bool)
    if domain.shape != masks[0].shape or any(np.any(mask & ~domain) for mask in masks) or np.any(target & ~domain):
        raise ValueError("all masks must be contained in a matching domain")
    if any(not mask.any() for mask in masks[-2:]):
        raise ValueError("the final two training masks must be nonempty for the volume-trend baseline")

    fit = fi.fit_fk_masks(masks, days, h, domain=domain)
    pde = fi.forward_masks(masks[0], domain, fit["rho"], fit["D"], [days[0], float(target_day)], h)[0]
    persistence = masks[-1].copy()
    previous_volume, last_volume = int(masks[-2].sum()), int(masks[-1].sum())
    growth_per_day = np.log(last_volume / previous_volume) / float(days[-1] - days[-2])
    horizon = float(target_day - days[-1])
    trend_volume = int(np.rint(last_volume * np.exp(growth_per_day * horizon)))
    trend_volume = int(np.clip(trend_volume, 0, int(domain.sum())))
    volume_trend = project_mask_volume(masks[-1], trend_volume, domain=domain, spacing_mm=h)

    training_spatial_residual = float(np.mean(1.0 - np.asarray(fit["dice"], dtype=float)))
    volume_sse = float(fit["sse_log_volume"])
    if not np.isfinite(volume_sse) or volume_sse < 0:
        raise ValueError("fit sse_log_volume must be finite and nonnegative")
    training_volume_residual_score = volume_sse / (1.0 + volume_sse)
    predictions = {"pde": pde, "persistence": persistence, "volume_trend": volume_trend}
    metrics = {name: all_metrics(pred, target, spacing=h) for name, pred in predictions.items()}
    return {"fit": fit, "target_day": float(target_day), "target_volume_trend_voxels": trend_volume,
            "training_spatial_residual": training_spatial_residual,
            "training_volume_residual_score": training_volume_residual_score,
            "predictions": predictions, "metrics": metrics}
