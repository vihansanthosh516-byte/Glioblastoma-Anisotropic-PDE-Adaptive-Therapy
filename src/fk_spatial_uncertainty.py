"""Spatially correlated boundary draws and proper scores for Track B forecasts.

All distances use physical millimeters. The displacement law is a simulation or fitted-model assumption;
these utilities do not turn supplied segmentations into scanner-noise measurements.
"""
from __future__ import annotations

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree


STRATUM_NAMES = ("0-2 mm", "2-5 mm", "5-10 mm")


def _spacing3(spacing_mm) -> np.ndarray:
    spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    return spacing


def _nonnegative_finite(value, name: str) -> float:
    value = float(value)
    if not np.isfinite(value) or value < 0:
        raise ValueError(f"{name} must be finite and nonnegative")
    return value


def deterministic_surface_anchors(mask, *, spacing_mm=1.0, n_anchors: int = 256) -> np.ndarray:
    """Choose reproducible, spatially spread inside-surface voxel centers in millimeters."""
    mask = np.asarray(mask, dtype=bool)
    _spacing3(spacing_mm)
    if mask.ndim != 3 or not mask.any() or mask.all():
        raise ValueError("mask must be a nonempty, non-full 3D binary array")
    if isinstance(n_anchors, (bool, np.bool_)) or int(n_anchors) != n_anchors or n_anchors < 1:
        raise ValueError("n_anchors must be a positive integer")
    surface = mask & ~ndimage.binary_erosion(
        mask, structure=ndimage.generate_binary_structure(3, 1), border_value=0
    )
    coords = np.argwhere(surface).astype(float)
    if not len(coords):
        raise ValueError("mask has no surface voxels")
    spacing = _spacing3(spacing_mm)
    coords *= spacing
    count = min(int(n_anchors), len(coords))
    if count == len(coords):
        return coords

    center = coords.mean(axis=0)
    first = int(np.argmax(np.sum((coords - center) ** 2, axis=1)))
    chosen = np.empty(count, dtype=int)
    chosen[0] = first
    min_sq_dist = np.sum((coords - coords[first]) ** 2, axis=1)
    for index in range(1, count):
        next_index = int(np.argmax(min_sq_dist))
        chosen[index] = next_index
        sq_dist = np.sum((coords - coords[next_index]) ** 2, axis=1)
        min_sq_dist = np.minimum(min_sq_dist, sq_dist)
    return coords[chosen]


def exponential_covariance_matrix(anchor_points_mm, *, sigma_b_mm: float, tau_b_mm: float,
                                 ell_mm: float) -> np.ndarray:
    """Build Sigma_kl = sigma_b^2 exp(-distance/ell) + tau_b^2 I."""
    points = np.asarray(anchor_points_mm, dtype=float)
    sigma_b_mm = _nonnegative_finite(sigma_b_mm, "sigma_b_mm")
    tau_b_mm = _nonnegative_finite(tau_b_mm, "tau_b_mm")
    ell_mm = float(ell_mm)
    if points.ndim != 2 or points.shape[1] != 3 or len(points) < 1 or not np.isfinite(points).all():
        raise ValueError("anchor_points_mm must be a nonempty finite array with shape (n, 3)")
    if not np.isfinite(ell_mm) or ell_mm <= 0:
        raise ValueError("ell_mm must be positive and finite")
    distances = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    covariance = sigma_b_mm**2 * np.exp(-distances / ell_mm)
    covariance.flat[::len(points) + 1] += tau_b_mm**2
    return covariance


def sample_anchor_displacements(anchor_points_mm, *, sigma_b_mm: float, tau_b_mm: float,
                                ell_mm: float, seed: int = 0) -> np.ndarray:
    """Draw one zero-mean Gaussian displacement vector at the supplied surface anchors."""
    covariance = exponential_covariance_matrix(
        anchor_points_mm, sigma_b_mm=sigma_b_mm, tau_b_mm=tau_b_mm, ell_mm=ell_mm
    )
    if not np.any(covariance):
        return np.zeros(len(covariance), dtype=float)
    try:
        factor = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("boundary covariance is not positive definite") from exc
    return factor @ np.random.default_rng(seed).normal(size=len(covariance))


def correlated_boundary_perturbation(mask, *, sigma_b_mm: float, tau_b_mm: float, ell_mm: float,
                                     spacing_mm=1.0, seed: int = 0, n_anchors: int = 256,
                                     anchor_points_mm=None, band_mm: float = 10.0) -> np.ndarray:
    """Move a signed-distance contour using a correlated Gaussian field on surface anchors.

    Anchor displacements are interpolated to voxels in the boundary band from their three nearest anchors,
    with weights proportional to (distance + 0.1 mm)^-2. Updating the signed-distance level set approximates
    moving the contour along its local normals, then rasterizes it back to a binary mask.
    """
    mask = np.asarray(mask, dtype=bool)
    spacing = _spacing3(spacing_mm)
    if mask.ndim != 3 or not mask.any() or mask.all():
        raise ValueError("mask must be a nonempty, non-full 3D binary array")
    sigma_b_mm = _nonnegative_finite(sigma_b_mm, "sigma_b_mm")
    tau_b_mm = _nonnegative_finite(tau_b_mm, "tau_b_mm")
    ell_mm = float(ell_mm)
    band_mm = float(band_mm)
    if not np.isfinite(ell_mm) or ell_mm <= 0 or not np.isfinite(band_mm) or band_mm <= 0:
        raise ValueError("ell_mm and band_mm must be positive and finite")
    if anchor_points_mm is None:
        anchors = deterministic_surface_anchors(mask, spacing_mm=spacing, n_anchors=n_anchors)
    else:
        anchors = np.asarray(anchor_points_mm, dtype=float)
        if anchors.ndim != 2 or anchors.shape[1] != 3 or len(anchors) < 1 or not np.isfinite(anchors).all():
            raise ValueError("anchor_points_mm must be a nonempty finite array with shape (n, 3)")
    if sigma_b_mm == 0 and tau_b_mm == 0:
        return mask.copy()

    displacement = sample_anchor_displacements(
        anchors, sigma_b_mm=sigma_b_mm, tau_b_mm=tau_b_mm, ell_mm=ell_mm, seed=seed
    )
    signed_distance = (ndimage.distance_transform_edt(mask, sampling=spacing)
                       - ndimage.distance_transform_edt(~mask, sampling=spacing))
    in_band = np.abs(signed_distance) <= band_mm
    voxel_indices = np.flatnonzero(in_band)
    if not len(voxel_indices):
        return mask.copy()
    voxel_coords_mm = np.column_stack(np.unravel_index(voxel_indices, mask.shape)).astype(float) * spacing
    k = min(3, len(anchors))
    distances, nearest = cKDTree(anchors).query(voxel_coords_mm, k=k)
    distances = np.asarray(distances, dtype=float).reshape(len(voxel_indices), k)
    nearest = np.asarray(nearest, dtype=int).reshape(len(voxel_indices), k)
    weights = (distances + 0.1) ** -2
    weights /= weights.sum(axis=1, keepdims=True)
    interpolated = np.sum(weights * displacement[nearest], axis=1)
    deformed = signed_distance.copy()
    deformed.flat[voxel_indices] += interpolated
    return deformed > 0


def fixed_boundary_score_pairs(history_mask, *, spacing_mm=1.0, n_voxels: int = 256,
                               band_mm: float = 10.0, max_pair_distance_mm: float = 10.0) -> dict:
    """Select a target-independent score design from the last observed history mask.

    Reuse this design for every forecast arm and held-out target in the corresponding case. Do not build it
    from the target mask: choosing scored locations after observing the outcome changes the scoring rule.
    """
    mask = np.asarray(history_mask, dtype=bool)
    spacing = _spacing3(spacing_mm)
    band_mm = float(band_mm)
    max_pair_distance_mm = float(max_pair_distance_mm)
    if mask.ndim != 3 or not mask.any() or mask.all():
        raise ValueError("history_mask must be a nonempty, non-full 3D binary array")
    if isinstance(n_voxels, (bool, np.bool_)) or int(n_voxels) != n_voxels or n_voxels < 2:
        raise ValueError("n_voxels must be an integer of at least two")
    if not np.isfinite(band_mm) or band_mm <= 0 or not np.isfinite(max_pair_distance_mm) or max_pair_distance_mm <= 0:
        raise ValueError("band and maximum pair distance must be positive and finite")
    signed_distance = (ndimage.distance_transform_edt(mask, sampling=spacing)
                       - ndimage.distance_transform_edt(~mask, sampling=spacing))
    candidates = np.flatnonzero(np.abs(signed_distance) <= band_mm)
    if len(candidates) < 2:
        raise ValueError("boundary band contains fewer than two voxels")
    if len(candidates) > int(n_voxels):
        positions = np.linspace(0, len(candidates) - 1, int(n_voxels), dtype=int)
        voxel_indices = candidates[positions]
    else:
        voxel_indices = candidates
    voxel_coords_mm = np.column_stack(np.unravel_index(voxel_indices, mask.shape)).astype(float) * spacing
    pair_i, pair_j = np.triu_indices(len(voxel_indices), k=1)
    distances = np.linalg.norm(voxel_coords_mm[pair_i] - voxel_coords_mm[pair_j], axis=1)
    keep = (distances > 0) & (distances <= max_pair_distance_mm)
    pair_i, pair_j, distances = pair_i[keep], pair_j[keep], distances[keep]
    strata = np.full(len(distances), -1, dtype=int)
    strata[distances <= 2.0] = 0
    strata[(distances > 2.0) & (distances <= 5.0)] = 1
    strata[(distances > 5.0) & (distances <= 10.0)] = 2
    return {
        "voxel_indices": np.asarray(voxel_indices, dtype=int),
        "pairs": np.column_stack((pair_i, pair_j)).astype(int),
        "distance_mm": distances,
        "strata": strata,
        "stratum_names": STRATUM_NAMES,
    }


def spatial_variogram_score(draw_masks, target_mask, pair_design: dict) -> dict:
    """Score mask-draw dependence using the preselected, equally weighted distance strata."""
    draws = np.asarray(draw_masks, dtype=bool)
    target = np.asarray(target_mask, dtype=bool)
    if draws.ndim != 4 or draws.shape[0] < 1 or target.shape != draws.shape[1:]:
        raise ValueError("draw_masks must be a nonempty (draw, x, y, z) array matching target_mask")
    voxel_indices = np.asarray(pair_design.get("voxel_indices"), dtype=int)
    pairs = np.asarray(pair_design.get("pairs"), dtype=int)
    strata = np.asarray(pair_design.get("strata"), dtype=int)
    names = tuple(pair_design.get("stratum_names", STRATUM_NAMES))
    if voxel_indices.ndim != 1 or np.any(voxel_indices < 0) or np.any(voxel_indices >= target.size):
        raise ValueError("pair-design voxel indices are invalid")
    if pairs.ndim != 2 or pairs.shape[1] != 2 or strata.ndim != 1 or len(strata) != len(pairs):
        raise ValueError("pair design must contain aligned (n, 2) pairs and stratum labels")
    if len(pairs) and (np.any(pairs < 0) or np.any(pairs >= len(voxel_indices))):
        raise ValueError("pair indices are outside the selected voxel list")
    if any(label < -1 or label >= len(names) for label in strata):
        raise ValueError("pair design contains an unknown stratum label")

    draw_flat = draws.reshape(len(draws), -1)[:, voxel_indices]
    target_flat = target.ravel()[voxel_indices]
    if len(pairs):
        expected_difference = np.mean(draw_flat[:, pairs[:, 0]] != draw_flat[:, pairs[:, 1]], axis=0)
        target_difference = (target_flat[pairs[:, 0]] != target_flat[pairs[:, 1]]).astype(float)
        losses = (expected_difference - target_difference) ** 2
    else:
        losses = np.empty(0, dtype=float)

    per_stratum = {}
    populated_scores = []
    empty = []
    for index, name in enumerate(names):
        selected = losses[strata == index]
        if len(selected):
            score = float(np.mean(selected))
            populated_scores.append(score)
            per_stratum[name] = {"n_pairs": int(len(selected)), "score": score}
        else:
            empty.append(name)
            per_stratum[name] = {"n_pairs": 0, "score": None}
    return {
        "score": float(np.mean(populated_scores)) if populated_scores else None,
        "strata": per_stratum,
        "empty_strata": empty,
        "aggregation": "equal weight over nonempty distance strata",
    }


def marginal_brier_score(draw_masks, target_mask, voxel_indices) -> float:
    """Compute the marginal Brier score on a fixed, target-independent voxel set."""
    draws = np.asarray(draw_masks, dtype=bool)
    target = np.asarray(target_mask, dtype=bool)
    indices = np.asarray(voxel_indices, dtype=int)
    if draws.ndim != 4 or draws.shape[0] < 1 or target.shape != draws.shape[1:]:
        raise ValueError("draw_masks must be a nonempty (draw, x, y, z) array matching target_mask")
    if indices.ndim != 1 or len(indices) == 0 or np.any(indices < 0) or np.any(indices >= target.size):
        raise ValueError("voxel_indices must be a nonempty valid flat-index array")
    probability = draws.reshape(len(draws), -1)[:, indices].mean(axis=0)
    observed = target.ravel()[indices].astype(float)
    return float(np.mean((probability - observed) ** 2))


def energy_score_mask(draw_masks, target_mask, voxel_indices) -> float:
    """Compute the empirical energy score on fixed selected voxels."""
    draws = np.asarray(draw_masks, dtype=bool)
    target = np.asarray(target_mask, dtype=bool)
    indices = np.asarray(voxel_indices, dtype=int)
    if draws.ndim != 4 or draws.shape[0] < 1 or target.shape != draws.shape[1:]:
        raise ValueError("draw_masks must be a nonempty (draw, x, y, z) array matching target_mask")
    if indices.ndim != 1 or len(indices) == 0 or np.any(indices < 0) or np.any(indices >= target.size):
        raise ValueError("voxel_indices must be a nonempty valid flat-index array")
    samples = draws.reshape(len(draws), -1)[:, indices].astype(float)
    observed = target.ravel()[indices].astype(float)
    first_term = np.linalg.norm(samples - observed[None, :], axis=1).mean()
    pairwise = np.linalg.norm(samples[:, None, :] - samples[None, :, :], axis=2).mean()
    return float(first_term - 0.5 * pairwise)


def volume_crps(draw_volumes, target_volume: float) -> float:
    """Compute the empirical CRPS for predictive volumes in a common physical unit."""
    draws = np.asarray(draw_volumes, dtype=float)
    target_volume = float(target_volume)
    if draws.ndim != 1 or len(draws) == 0 or not np.isfinite(draws).all():
        raise ValueError("draw_volumes must be a nonempty finite 1D array")
    if not np.isfinite(target_volume):
        raise ValueError("target_volume must be finite")
    return float(np.mean(np.abs(draws - target_volume))
                 - 0.5 * np.mean(np.abs(draws[:, None] - draws[None, :])))
