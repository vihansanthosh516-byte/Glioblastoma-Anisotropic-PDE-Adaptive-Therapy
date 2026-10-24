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
from src.forecast_metrics import all_metrics


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

    predictions = {"pde": pde, "persistence": persistence, "volume_trend": volume_trend}
    metrics = {name: all_metrics(pred, target, spacing=h) for name, pred in predictions.items()}
    return {"fit": fit, "target_day": float(target_day), "target_volume_trend_voxels": trend_volume,
            "predictions": predictions, "metrics": metrics}
