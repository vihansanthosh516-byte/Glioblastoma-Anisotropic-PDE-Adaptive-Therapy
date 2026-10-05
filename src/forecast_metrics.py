"""Spatial and volume forecast metrics (analysis_plan_v1.md; masterplan §22).

All functions take boolean 3D masks (same shape) and the voxel spacing in mm. Empty-mask rules are
explicit: overlap metrics return 0.0 when exactly one mask is empty and 1.0 when both are empty;
distance metrics return nan when either mask is empty.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage


def _check(a, b):
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    assert a.shape == b.shape, f"shape mismatch {a.shape} vs {b.shape}"
    return a, b


def dice(a, b) -> float:
    a, b = _check(a, b)
    s = a.sum() + b.sum()
    return 1.0 if s == 0 else float(2.0 * (a & b).sum() / s)


def jaccard(a, b) -> float:
    a, b = _check(a, b)
    u = (a | b).sum()
    return 1.0 if u == 0 else float((a & b).sum() / u)


def _surface(m):
    return m & ~ndimage.binary_erosion(m, border_value=0)


def _surface_distances(a, b, spacing):
    """Distances (mm) from each surface voxel of a to the surface of b."""
    sb = _surface(b)
    dist_to_b = ndimage.distance_transform_edt(~sb, sampling=spacing)
    return dist_to_b[_surface(a)]


def hd95(a, b, spacing=1.0) -> float:
    a, b = _check(a, b)
    if not a.any() or not b.any():
        return float("nan")
    d = np.concatenate([_surface_distances(a, b, spacing), _surface_distances(b, a, spacing)])
    return float(np.percentile(d, 95))


def asd(a, b, spacing=1.0) -> float:
    """Average symmetric surface distance (mm)."""
    a, b = _check(a, b)
    if not a.any() or not b.any():
        return float("nan")
    d1, d2 = _surface_distances(a, b, spacing), _surface_distances(b, a, spacing)
    return float((d1.sum() + d2.sum()) / (len(d1) + len(d2)))


def volume_mm3(m, spacing=1.0) -> float:
    vox = float(np.prod(spacing)) if np.ndim(spacing) else float(spacing) ** 3
    return float(np.asarray(m, bool).sum() * vox)


def symmetric_pct_error(v_pred: float, v_true: float) -> float:
    """100 * |pred - true| / ((pred + true) / 2); 0 when both are 0."""
    den = (v_pred + v_true) / 2.0
    return 0.0 if den == 0 else float(100.0 * abs(v_pred - v_true) / den)


def log_volume_ratio(v_pred: float, v_true: float, eps: float = 1.0) -> float:
    """ln((pred + eps) / (true + eps)); eps keeps empty masks finite."""
    return float(np.log((v_pred + eps) / (v_true + eps)))


def centroid_displacement(a, b, spacing=1.0) -> float:
    a, b = _check(a, b)
    if not a.any() or not b.any():
        return float("nan")
    sp = np.broadcast_to(np.asarray(spacing, float), (3,))
    ca = np.argwhere(a).mean(0) * sp
    cb = np.argwhere(b).mean(0) * sp
    return float(np.linalg.norm(ca - cb))


def all_metrics(pred, truth, spacing=1.0) -> dict:
    return {"dice": dice(pred, truth), "jaccard": jaccard(pred, truth), "hd95_mm": hd95(pred, truth, spacing),
            "asd_mm": asd(pred, truth, spacing),
            "spe_pct": symmetric_pct_error(volume_mm3(pred, spacing), volume_mm3(truth, spacing)),
            "log_vol_ratio": log_volume_ratio(volume_mm3(pred, spacing), volume_mm3(truth, spacing)),
            "centroid_mm": centroid_displacement(pred, truth, spacing)}
