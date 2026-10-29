"""Target-only scoring callback for the first-two-to-third Track B forecasts.

The forecast code must not receive the target. This module receives the forecast output and target only after
the forecast is complete, thresholds hybrid inclusion probabilities at 0.5, volume-matches the no-PDE point
mask to that hybrid mask, and returns one patient score row.
Use :func:`make_first_third_score_callback` with the grid's voxel spacing for physical volume and surface metrics.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from scipy.ndimage import binary_erosion, distance_transform_edt, generate_binary_structure

from src.fk_spatial_uncertainty import (
    energy_score_mask,
    fixed_boundary_score_pairs,
    marginal_brier_score,
    spatial_variogram_score,
)

BOUNDARY_SCORE_VOXELS = 128
BOUNDARY_SCORE_BAND_MM = 10.0
BOUNDARY_SCORE_MAX_PAIR_DISTANCE_MM = 10.0


def _binary_mask(value, *, name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim != 3:
        raise ValueError(f"{name} must be a 3D binary mask")
    if array.dtype == np.bool_:
        return array.astype(bool, copy=False)
    if (not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all()
            or not np.isin(array, (0, 1)).all()):
        raise ValueError(f"{name} must contain only finite binary values")
    return array.astype(bool, copy=False)


def _ensemble_median_mask(draws, *, name: str, expected_shape: tuple[int, int, int]) -> np.ndarray:
    array = np.asarray(draws)
    if array.ndim != 4 or tuple(array.shape[1:]) != expected_shape or array.shape[0] < 1:
        raise ValueError(f"{name} must have shape (n_draws, *target_shape) with at least one draw")
    if array.dtype != np.bool_:
        if (not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all()
                or not np.isin(array, (0, 1)).all()):
            raise ValueError(f"{name} must contain only finite binary masks")
    inclusion_probability = array.mean(axis=0, dtype=np.float32)
    if not np.isfinite(inclusion_probability).all():
        raise ValueError(f"{name} inclusion probabilities must be finite")
    return inclusion_probability >= 0.5


def _ensemble_probabilities_and_volumes(draws, *, name: str, expected_shape, support) -> tuple[np.ndarray, np.ndarray]:
    array = np.asarray(draws)
    if array.ndim != 4 or tuple(array.shape[1:]) != expected_shape or array.shape[0] < 1:
        raise ValueError(f"{name} must have shape (n_draws, *target_shape) with at least one draw")
    if array.dtype != np.bool_:
        if (not np.issubdtype(array.dtype, np.number) or not np.isfinite(array).all()
                or not np.isin(array, (0, 1)).all()):
            raise ValueError(f"{name} must contain only finite binary masks")
    probabilities = array.mean(axis=0, dtype=np.float32)
    if not np.isfinite(probabilities).all():
        raise ValueError(f"{name} inclusion probabilities must be finite")
    volumes = np.fromiter(
        (np.count_nonzero(draw & support) for draw in array), dtype=np.int64, count=len(array),
    )
    return probabilities, volumes


def _ranked_probability_mask(probabilities, *, support, target_voxels: int) -> np.ndarray:
    """Select top probabilities inside support, breaking ties by ascending flat index."""
    probabilities = np.asarray(probabilities, dtype=float)
    support = np.asarray(support, dtype=bool)
    if (probabilities.ndim != 3 or support.shape != probabilities.shape
            or not np.isfinite(probabilities).all()):
        raise ValueError("probabilities and support must be matching finite 3D arrays")
    if isinstance(target_voxels, (bool, np.bool_)) or int(target_voxels) != target_voxels:
        raise ValueError("target_voxels must be an integer")
    target_voxels = int(target_voxels)
    indices = np.flatnonzero(support)
    if target_voxels < 0 or target_voxels > len(indices):
        raise ValueError("target_voxels must be between zero and support size")
    selected = np.zeros(probabilities.shape, dtype=bool)
    if target_voxels == 0:
        return selected
    values = probabilities.ravel()[indices]
    threshold_index = values.size - target_voxels
    threshold = np.partition(values, threshold_index)[threshold_index]
    above = values > threshold
    selected.ravel()[indices[above]] = True
    tied_count = target_voxels - int(above.sum())
    if tied_count:
        tied_indices = indices[np.flatnonzero(values == threshold)]
        selected.ravel()[tied_indices[:tied_count]] = True
    return selected


def dice_with_empty_rule(prediction, target) -> float:
    """Return Dice; define both-empty as 1 and one-empty as 0, without dropping the case."""
    prediction = _binary_mask(prediction, name="prediction")
    target = _binary_mask(target, name="target")
    if prediction.shape != target.shape:
        raise ValueError("prediction and target masks must have the same shape")
    prediction_count = int(prediction.sum())
    target_count = int(target.sum())
    if prediction_count + target_count == 0:
        return 1.0
    return float(2 * int(np.count_nonzero(prediction & target)) / (prediction_count + target_count))


def _surface_distance_metrics(prediction, target, *, spacing_mm, support) -> dict:
    if not prediction.any() and not target.any():
        return {"mean_surface_distance_mm": 0.0, "hausdorff95_surface_distance_mm": 0.0,
                "surface_distance_status": "both_empty_zero_rule"}
    if not prediction.any() or not target.any():
        return {"mean_surface_distance_mm": None, "hausdorff95_surface_distance_mm": None,
                "surface_distance_status": "undefined_one_empty"}
    if not support.all():
        return {"mean_surface_distance_mm": None, "hausdorff95_surface_distance_mm": None,
                "surface_distance_status": "undefined_partial_support"}

    structure = generate_binary_structure(3, 1)
    prediction_surface = prediction & ~binary_erosion(prediction, structure=structure, border_value=0)
    target_surface = target & ~binary_erosion(target, structure=structure, border_value=0)
    distance_to_prediction = distance_transform_edt(~prediction_surface, sampling=spacing_mm)
    distance_to_target = distance_transform_edt(~target_surface, sampling=spacing_mm)
    distances = np.concatenate((
        distance_to_target[prediction_surface],
        distance_to_prediction[target_surface],
    ))
    return {
        "mean_surface_distance_mm": float(distances.mean()),
        "hausdorff95_surface_distance_mm": float(np.quantile(distances, 0.95)),
        "surface_distance_status": "defined; 6-connected boundary voxels; symmetric directed distances",
    }


def _volume_error_metrics(prediction_voxels: int, target_voxels: int, *, voxel_volume_mm3: float) -> dict:
    error_mm3 = float((prediction_voxels - target_voxels) * voxel_volume_mm3)
    return {
        "volume_error_mm3": error_mm3,
        "absolute_volume_error_mm3": abs(error_mm3),
        "relative_volume_error": float(prediction_voxels / target_voxels - 1.0) if target_voxels else None,
    }


def _volume_distribution_metrics(sample_volumes_voxels, target_voxels: int, *, voxel_volume_mm3: float) -> dict:
    samples = np.asarray(sample_volumes_voxels, dtype=float)
    if (samples.ndim != 1 or not len(samples) or not np.isfinite(samples).all()
            or np.any(samples < 0)):
        raise ValueError("volume forecast samples must be finite nonnegative voxel counts")
    samples *= voxel_volume_mm3
    target_volume = float(target_voxels * voxel_volume_mm3)
    pairwise_distance = np.abs(samples[:, None] - samples[None, :]).mean()
    crps = float(np.abs(samples - target_volume).mean() - 0.5 * pairwise_distance)
    lower, upper = np.quantile(samples, [0.05, 0.95])
    return {
        "volume_crps_mm3": crps,
        "volume_interval_90_lower_mm3": float(lower),
        "volume_interval_90_upper_mm3": float(upper),
        "volume_interval_90_covered": bool(lower <= target_volume <= upper),
        "volume_interval_90_sample_count": int(len(samples)),
    }


def _history_only_boundary_design(history, *, support, spacing_mm):
    if history is None:
        return None, "unavailable_no_history_mask"
    if not support.all():
        return None, "undefined_partial_support"
    if not history.any() or history.all():
        return None, "unavailable_empty_or_full_history_mask"
    try:
        design = fixed_boundary_score_pairs(
            history, spacing_mm=spacing_mm, n_voxels=BOUNDARY_SCORE_VOXELS,
            band_mm=BOUNDARY_SCORE_BAND_MM,
            max_pair_distance_mm=BOUNDARY_SCORE_MAX_PAIR_DISTANCE_MM,
        )
    except ValueError as exc:
        return None, f"unavailable_history_only_boundary_design: {exc}"
    return design, "selected from history only; fixed before target scoring"


def _get_target_mask(target) -> np.ndarray:
    if isinstance(target, Mapping):
        if "mask" not in target:
            raise ValueError("target mapping must contain a mask")
        target = target["mask"]
    elif hasattr(target, "mask"):
        target = target.mask
    return _binary_mask(target, name="target")


def _get_prediction_masks(prediction, target_shape, *, support, spacing_mm):
    """Resolve a successful ensemble or explicit persistence-fallback masks."""
    if not isinstance(prediction, Mapping):
        fallback = _binary_mask(prediction, name="explicit fallback mask")
        if fallback.shape != target_shape:
            raise ValueError("fallback mask and target must have the same shape")
        design, design_status = _history_only_boundary_design(
            None, support=support, spacing_mm=spacing_mm,
        )
        return (fallback, fallback, fallback, fallback, "persistence_fallback", {}, {}, {}, design,
                design_status, "identical persistence fallback masks; paired comparator difference is zero")

    history_raw = prediction.get("history_mask", prediction.get("persistence_mask"))
    history = None if history_raw is None else _binary_mask(history_raw, name="history/persistence mask")
    if history is not None and history.shape != target_shape:
        raise ValueError("history/persistence mask and target must have the same shape")

    fallback_masks = prediction.get("fallback_masks")
    probability_fields = {}
    volume_samples = {}
    spatial_draws = {}
    point_volume_rule = "identical persistence fallback masks; paired comparator difference is zero"
    if fallback_masks is not None:
        if not isinstance(fallback_masks, Mapping):
            raise ValueError("fallback_masks must map hybrid and no_pde to explicit masks")
        if "hybrid" not in fallback_masks or "no_pde" not in fallback_masks:
            raise ValueError("fallback_masks must contain hybrid and no_pde masks")
        hybrid = _binary_mask(fallback_masks["hybrid"], name="hybrid fallback mask")
        no_pde = _binary_mask(fallback_masks["no_pde"], name="no-PDE fallback mask")
        persistence_raw = fallback_masks.get("persistence", history)
        if persistence_raw is None:
            raise ValueError("persistence mask is required for fallback scoring")
        persistence = _binary_mask(persistence_raw, name="persistence fallback mask")
        geometric = _binary_mask(
            fallback_masks.get("geometric_train_rate", persistence),
            name="geometric train-rate fallback mask",
        )
        geometric_status = "persistence_fallback"
        if int(hybrid.sum()) != int(no_pde.sum()):
            raise ValueError("fallback masks must have equal volumes")
        point_volume_rule = "explicit fallback masks supplied; equal voxel counts required"
    elif "fallback_mask" in prediction:
        hybrid = _binary_mask(prediction["fallback_mask"], name="explicit fallback mask")
        no_pde = hybrid
        if history is None:
            raise ValueError("history_mask is required to score persistence separately")
        persistence = history
        geometric = persistence
        geometric_status = "persistence_fallback"
    else:
        forecast = prediction.get("forecast", prediction)
        if not isinstance(forecast, Mapping):
            raise ValueError("prediction must contain a hybrid forecast or explicit fallback masks")
        if "mask_draws" not in forecast or "no_pde_ranker_mask_draws" not in forecast:
            raise ValueError("hybrid forecast must contain mask_draws and no_pde_ranker_mask_draws")
        if history is None:
            raise ValueError("history_mask is required to score persistence")
        probability_fields["hybrid"], volume_samples["hybrid"] = _ensemble_probabilities_and_volumes(
            forecast["mask_draws"], name="mask_draws", expected_shape=target_shape, support=support,
        )
        spatial_draws["hybrid"] = forecast["mask_draws"]
        probability_fields["no_pde"], volume_samples["no_pde"] = _ensemble_probabilities_and_volumes(
            forecast["no_pde_ranker_mask_draws"], name="no_pde_ranker_mask_draws",
            expected_shape=target_shape, support=support,
        )
        spatial_draws["no_pde"] = forecast["no_pde_ranker_mask_draws"]
        hybrid = probability_fields["hybrid"] >= 0.5
        hybrid &= support
        no_pde = _ranked_probability_mask(
            probability_fields["no_pde"], support=support,
            target_voxels=int(hybrid.sum()),
        )
        point_volume_rule = (
            "rank no-PDE inclusion probabilities to hybrid point-mask voxel count; stable flat-index ties"
        )
        persistence = history
        if "geometric_train_rate_mask" in forecast:
            geometric = _binary_mask(
                forecast["geometric_train_rate_mask"], name="geometric_train_rate_mask",
            )
            geometric_status = "forecast"
        elif "geometric_train_rate_mask_draws" in forecast:
            geometric = _ensemble_median_mask(
                forecast["geometric_train_rate_mask_draws"],
                name="geometric_train_rate_mask_draws", expected_shape=target_shape,
            )
            geometric_status = "forecast"
        else:
            geometric = persistence
            geometric_status = "unavailable_persistence_fallback"

    for name, mask in (("hybrid", hybrid), ("no-PDE", no_pde),
                       ("persistence", persistence), ("geometric train-rate", geometric)):
        if mask.shape != target_shape:
            raise ValueError(f"{name} prediction and target must have the same shape")
    design, design_status = _history_only_boundary_design(
        history, support=support, spacing_mm=spacing_mm,
    )
    return (hybrid, no_pde, persistence, geometric, geometric_status, probability_fields,
            volume_samples, spatial_draws, design, design_status, point_volume_rule)


def score_first_third_target(prediction, target, *, voxel_volume_mm3: float = 1.0,
                             voxel_spacing_mm=None, valid_support=None) -> dict:
    """Score one target after forecasting; return finite, patient-level paired metrics.

    ``prediction`` may be the hybrid result directly, a wrapper containing ``forecast`` and ``history_mask``,
    a mapping with explicit ``fallback_masks`` (hybrid/no_pde and optionally persistence), a single
    ``fallback_mask`` plus ``history_mask``, or a hard fallback mask used for all three arms. For an ordinary
    PDE failure, pass the persistence mask as the same hybrid and no-PDE fallback; its paired delta is then zero.

    Dice treats an empty target and empty prediction as 1, and exactly one empty mask as 0. Every target remains
    in the score output, including empty targets. Volume errors and mean/95th-percentile surface distances are
    secondary metrics; surface distances are undefined for one empty mask or partial support.
    """
    try:
        voxel_volume_mm3 = float(voxel_volume_mm3)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("voxel_volume_mm3 must be a positive finite scalar") from exc
    if not np.isfinite(voxel_volume_mm3) or voxel_volume_mm3 <= 0:
        raise ValueError("voxel_volume_mm3 must be a positive finite scalar")
    if voxel_spacing_mm is None:
        spacing_mm = np.full(3, voxel_volume_mm3 ** (1.0 / 3.0), dtype=float)
    else:
        try:
            spacing_mm = np.broadcast_to(np.asarray(voxel_spacing_mm, dtype=float), (3,)).copy()
        except (TypeError, ValueError) as exc:
            raise ValueError("voxel_spacing_mm must be one positive value or three positive finite values") from exc
        if not np.isfinite(spacing_mm).all() or np.any(spacing_mm <= 0):
            raise ValueError("voxel_spacing_mm must be one positive value or three positive finite values")
        if not np.isclose(float(np.prod(spacing_mm)), voxel_volume_mm3, rtol=1e-9, atol=1e-12):
            raise ValueError("voxel_volume_mm3 must match the product of voxel_spacing_mm")
    target_mask = _get_target_mask(target)
    if valid_support is None:
        support_mask = np.ones(target_mask.shape, dtype=bool)
    else:
        support_mask = _binary_mask(valid_support, name="valid target support")
        if support_mask.shape != target_mask.shape or not support_mask.any():
            raise ValueError("valid target support must be nonempty and match the target grid")
    (hybrid_mask, no_pde_mask, persistence_mask, geometric_mask, geometric_status,
     probability_fields, volume_samples, spatial_draws, boundary_design, boundary_design_status,
     point_volume_rule) = (
        _get_prediction_masks(
            prediction, target_mask.shape, support=support_mask, spacing_mm=spacing_mm,
        )
    )
    if np.any(target_mask & ~support_mask):
        raise ValueError("target foreground falls outside validated target support")
    for name, mask in (("hybrid", hybrid_mask), ("no-PDE", no_pde_mask),
                       ("persistence", persistence_mask), ("geometric train-rate", geometric_mask)):
        if np.any(mask & ~support_mask):
            raise ValueError(f"{name} prediction falls outside validated target support")
    target_mask = target_mask & support_mask
    hybrid_mask = hybrid_mask & support_mask
    no_pde_mask = no_pde_mask & support_mask
    persistence_mask = persistence_mask & support_mask
    geometric_mask = geometric_mask & support_mask

    target_voxels = int(target_mask.sum())
    hybrid_voxels = int(hybrid_mask.sum())
    no_pde_voxels = int(no_pde_mask.sum())
    persistence_voxels = int(persistence_mask.sum())
    geometric_voxels = int(geometric_mask.sum())
    dice_hybrid = dice_with_empty_rule(hybrid_mask, target_mask)
    dice_no_pde = dice_with_empty_rule(no_pde_mask, target_mask)
    dice_persistence = dice_with_empty_rule(persistence_mask, target_mask)
    paired_delta = float(dice_hybrid - dice_no_pde)
    result = {
        "dice_hybrid": float(dice_hybrid),
        "dice_no_pde_matched_volume": float(dice_no_pde),
        "dice_persistence": float(dice_persistence),
        "dice_geometric_train_rate": float(dice_with_empty_rule(geometric_mask, target_mask)),
        "paired_delta_dice_hybrid_minus_no_pde": paired_delta,
        "paired_delta_dice_hybrid_minus_geometric_train_rate": float(
            dice_hybrid - dice_with_empty_rule(geometric_mask, target_mask)
        ),
        "target_voxels": target_voxels,
        "hybrid_prediction_voxels": hybrid_voxels,
        "no_pde_prediction_voxels": no_pde_voxels,
        "persistence_prediction_voxels": persistence_voxels,
        "geometric_train_rate_prediction_voxels": geometric_voxels,
        "target_volume_mm3": float(target_voxels * voxel_volume_mm3),
        "hybrid_prediction_volume_mm3": float(hybrid_voxels * voxel_volume_mm3),
        "no_pde_prediction_volume_mm3": float(no_pde_voxels * voxel_volume_mm3),
        "persistence_prediction_volume_mm3": float(persistence_voxels * voxel_volume_mm3),
        "geometric_train_rate_prediction_volume_mm3": float(geometric_voxels * voxel_volume_mm3),
        "target_empty": bool(target_voxels == 0),
        "valid_support_voxels": int(support_mask.sum()),
        "hybrid_prediction_empty": bool(hybrid_voxels == 0),
        "no_pde_prediction_empty": bool(no_pde_voxels == 0),
        "persistence_prediction_empty": bool(persistence_voxels == 0),
        "geometric_train_rate_prediction_empty": bool(geometric_voxels == 0),
        "geometric_train_rate_status": geometric_status,
        "ensemble_inclusion_threshold": 0.5,
        "no_pde_matched_volume_point_forecast_rule": point_volume_rule,
        "empty_mask_dice_rule": "both_empty=1; one_empty=0",
        "boundary_spatial_score_design_status": boundary_design_status,
        "boundary_spatial_score_design_voxels": (
            int(len(boundary_design["voxel_indices"])) if boundary_design is not None else 0
        ),
        "boundary_spatial_score_design_pairs": (
            int(len(boundary_design["pairs"])) if boundary_design is not None else 0
        ),
    }
    for name, predicted_mask in (
        ("hybrid", hybrid_mask), ("no_pde", no_pde_mask), ("persistence", persistence_mask),
        ("geometric_train_rate", geometric_mask),
    ):
        result.update({
            f"{key}_{name}": value
            for key, value in _volume_error_metrics(
                int(predicted_mask.sum()), target_voxels, voxel_volume_mm3=voxel_volume_mm3,
            ).items()
        })
        result.update({
            f"{key}_{name}": value
            for key, value in _surface_distance_metrics(
                predicted_mask, target_mask, spacing_mm=spacing_mm, support=support_mask,
            ).items()
        })
        probabilities = probability_fields.get(name, predicted_mask.astype(np.float32))
        result[f"marginal_voxel_brier_{name}"] = float(
            np.mean((probabilities[support_mask] - target_mask[support_mask]) ** 2)
        )
        samples = volume_samples.get(name, np.asarray([predicted_mask.sum()], dtype=np.int64))
        result.update({
            f"{key}_{name}": value
            for key, value in _volume_distribution_metrics(
                samples, target_voxels, voxel_volume_mm3=voxel_volume_mm3,
            ).items()
        })
        if boundary_design is None:
            result[f"boundary_energy_score_{name}"] = None
            result[f"boundary_variogram_score_{name}"] = None
            result[f"boundary_variogram_strata_{name}"] = None
            result[f"boundary_marginal_brier_{name}"] = None
            result[f"boundary_spatial_score_status_{name}"] = boundary_design_status
        else:
            draws = spatial_draws.get(name)
            if draws is None:
                draws = predicted_mask[None, ...]
            variogram = spatial_variogram_score(draws, target_mask, boundary_design)
            result[f"boundary_energy_score_{name}"] = float(
                energy_score_mask(draws, target_mask, boundary_design["voxel_indices"])
            )
            result[f"boundary_variogram_score_{name}"] = variogram["score"]
            result[f"boundary_variogram_strata_{name}"] = variogram["strata"]
            result[f"boundary_marginal_brier_{name}"] = float(
                marginal_brier_score(draws, target_mask, boundary_design["voxel_indices"])
            )
            result[f"boundary_spatial_score_status_{name}"] = "scored on locked history-only boundary design"
    for key, value in result.items():
        if isinstance(value, (float, np.floating)) and not np.isfinite(value):
            raise RuntimeError(f"non-finite score generated for {key}")
    return result


def make_first_third_score_callback(spacing_mm=2.0, *, valid_support_by_patient=None):
    """Create a ``score(prediction, target)`` callback with physical volumes in mm^3."""
    try:
        spacing = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    except (TypeError, ValueError) as exc:
        raise ValueError("spacing_mm must be one positive value or three positive finite values") from exc
    if not np.isfinite(spacing).all() or np.any(spacing <= 0):
        raise ValueError("spacing_mm must be one positive value or three positive finite values")
    voxel_volume_mm3 = float(np.prod(spacing))
    support_mapping = None
    if valid_support_by_patient is not None:
        if not isinstance(valid_support_by_patient, Mapping):
            raise ValueError("valid_support_by_patient must map patient IDs to support masks")
        support_mapping = {}
        for patient_id, support in valid_support_by_patient.items():
            key = str(patient_id).strip()
            if not key or key in support_mapping:
                raise ValueError("valid-support patient IDs must be nonempty and unique")
            support_mapping[key] = _binary_mask(support, name=f"valid support for {key}")

    def score(prediction, target):
        valid_support = None
        if support_mapping is not None:
            patient_id = (target.get("patient_id") if isinstance(target, Mapping)
                          else getattr(target, "patient_id", None))
            patient_id = str(patient_id or "").strip()
            if patient_id not in support_mapping:
                raise ValueError("no validated target support is available for this patient")
            valid_support = support_mapping[patient_id]
        return score_first_third_target(
            prediction, target, voxel_volume_mm3=voxel_volume_mm3,
            voxel_spacing_mm=tuple(spacing),
            valid_support=valid_support,
        )

    return score
