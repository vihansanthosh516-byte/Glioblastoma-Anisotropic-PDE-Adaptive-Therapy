"""Patient-level split-conformal bands for binary tumor-boundary forecasts."""
from __future__ import annotations

import math

import numpy as np
from scipy import ndimage


def _spacing3(spacing_mm) -> np.ndarray:
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    return spacing


def boundary_surface(mask) -> np.ndarray:
    """Return the six-connected inside boundary of a nonempty 3D mask."""
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim != 3 or not mask.any() or mask.all():
        raise ValueError("mask must be a nonempty, non-full 3D array")
    eroded = ndimage.binary_erosion(
        mask, structure=ndimage.generate_binary_structure(3, 1), border_value=0
    )
    return mask & ~eroded


def symmetric_boundary_hausdorff_mm(predicted_mask, target_mask, *, spacing_mm=1.0) -> float:
    """Compute the symmetric Hausdorff distance between two voxelized surfaces in millimeters."""
    spacing = _spacing3(spacing_mm)
    predicted_surface = boundary_surface(predicted_mask)
    target_surface = boundary_surface(target_mask)
    distance_to_predicted = ndimage.distance_transform_edt(~predicted_surface, sampling=spacing)
    distance_to_target = ndimage.distance_transform_edt(~target_surface, sampling=spacing)
    target_to_predicted = float(distance_to_predicted[target_surface].max())
    predicted_to_target = float(distance_to_target[predicted_surface].max())
    return max(target_to_predicted, predicted_to_target)


def split_conformal_radius(scores, *, alpha: float) -> dict:
    """Return the finite-sample split-conformal order statistic for patient-level scores.

    For exchangeable calibration and future patients, the returned radius covers the future patient's
    entire target surface with marginal probability at least 1-alpha, provided the predictor and score
    definition were fixed before calibration. A non-finite result means the calibration set is too small.
    """
    values = np.asarray(scores, dtype=float)
    alpha = float(alpha)
    if (values.ndim != 1 or not len(values) or not np.isfinite(values).all()
            or np.any(values < 0)):
        raise ValueError("scores must be a nonempty vector of finite nonnegative values")
    if not np.isfinite(alpha) or not 0.0 < alpha < 1.0:
        raise ValueError("alpha must be finite and strictly between zero and one")
    rank = int(math.ceil((len(values) + 1) * (1.0 - alpha)))
    if rank > len(values):
        return {"radius_mm": float("inf"), "rank": rank, "n_calibration_patients": len(values),
                "finite": False}
    radius = float(np.partition(values, rank - 1)[rank - 1])
    return {"radius_mm": radius, "rank": rank, "n_calibration_patients": len(values), "finite": True}


def boundary_band(predicted_mask, *, radius_mm: float, spacing_mm=1.0) -> np.ndarray:
    """Return voxels no farther than radius_mm from the predicted boundary surface."""
    radius_mm = float(radius_mm)
    if np.isnan(radius_mm) or radius_mm < 0:
        raise ValueError("radius_mm must be nonnegative")
    spacing = _spacing3(spacing_mm)
    surface = boundary_surface(predicted_mask)
    if np.isinf(radius_mm):
        return np.ones(surface.shape, dtype=bool)
    distance = ndimage.distance_transform_edt(~surface, sampling=spacing)
    return distance <= radius_mm


def target_boundary_coverage(predicted_mask, target_mask, *, radius_mm: float, spacing_mm=1.0) -> dict:
    """Measure whether and what fraction of the target surface lies inside a predicted boundary band."""
    band = boundary_band(predicted_mask, radius_mm=radius_mm, spacing_mm=spacing_mm)
    target_surface = boundary_surface(target_mask)
    covered = band[target_surface]
    return {
        "whole_boundary_covered": bool(covered.all()),
        "target_surface_fraction_covered": float(covered.mean()),
        "band_voxel_fraction": float(band.mean()),
    }
