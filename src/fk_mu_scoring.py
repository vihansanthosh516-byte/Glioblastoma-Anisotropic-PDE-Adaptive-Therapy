"""Target-only scoring callback for the first-two-to-third Track B forecasts.

The forecast code must not receive the target. This module receives the forecast output and target only after
the forecast is complete, thresholds ensemble inclusion probabilities at 0.5, and returns one paired score row.
Use :func:`make_first_third_score_callback` with the grid's voxel spacing to report volumes in mm^3.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np


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
    inclusion_probability = np.asarray(array, dtype=np.float64).mean(axis=0)
    if not np.isfinite(inclusion_probability).all():
        raise ValueError(f"{name} inclusion probabilities must be finite")
    return inclusion_probability >= 0.5


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


def _get_target_mask(target) -> np.ndarray:
    if isinstance(target, Mapping):
        if "mask" not in target:
            raise ValueError("target mapping must contain a mask")
        target = target["mask"]
    elif hasattr(target, "mask"):
        target = target.mask
    return _binary_mask(target, name="target")


def _get_prediction_masks(prediction, target_shape) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Resolve a successful ensemble or explicit persistence-fallback masks."""
    if not isinstance(prediction, Mapping):
        fallback = _binary_mask(prediction, name="explicit fallback mask")
        if fallback.shape != target_shape:
            raise ValueError("fallback mask and target must have the same shape")
        return fallback, fallback, fallback

    history_raw = prediction.get("history_mask", prediction.get("persistence_mask"))
    history = None if history_raw is None else _binary_mask(history_raw, name="history/persistence mask")
    if history is not None and history.shape != target_shape:
        raise ValueError("history/persistence mask and target must have the same shape")

    fallback_masks = prediction.get("fallback_masks")
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
    elif "fallback_mask" in prediction:
        hybrid = _binary_mask(prediction["fallback_mask"], name="explicit fallback mask")
        no_pde = hybrid
        if history is None:
            raise ValueError("history_mask is required to score persistence separately")
        persistence = history
    else:
        forecast = prediction.get("forecast", prediction)
        if not isinstance(forecast, Mapping):
            raise ValueError("prediction must contain a hybrid forecast or explicit fallback masks")
        if "mask_draws" not in forecast or "no_pde_ranker_mask_draws" not in forecast:
            raise ValueError("hybrid forecast must contain mask_draws and no_pde_ranker_mask_draws")
        if history is None:
            raise ValueError("history_mask is required to score persistence")
        hybrid = _ensemble_median_mask(
            forecast["mask_draws"], name="mask_draws", expected_shape=target_shape,
        )
        no_pde = _ensemble_median_mask(
            forecast["no_pde_ranker_mask_draws"],
            name="no_pde_ranker_mask_draws", expected_shape=target_shape,
        )
        persistence = history

    for name, mask in (("hybrid", hybrid), ("no-PDE", no_pde), ("persistence", persistence)):
        if mask.shape != target_shape:
            raise ValueError(f"{name} prediction and target must have the same shape")
    return hybrid, no_pde, persistence


def score_first_third_target(prediction, target, *, voxel_volume_mm3: float = 1.0,
                             valid_support=None) -> dict:
    """Score one target after forecasting; return finite, patient-level paired metrics.

    ``prediction`` may be the hybrid result directly, a wrapper containing ``forecast`` and ``history_mask``,
    a mapping with explicit ``fallback_masks`` (hybrid/no_pde and optionally persistence), a single
    ``fallback_mask`` plus ``history_mask``, or a hard fallback mask used for all three arms. For an ordinary
    PDE failure, pass the persistence mask as the same hybrid and no-PDE fallback; its paired delta is then zero.

    Dice treats an empty target and empty prediction as 1, and exactly one empty mask as 0. Every target remains
    in the score output, including empty targets. Volumes use voxel counts multiplied by voxel_volume_mm3.
    """
    try:
        voxel_volume_mm3 = float(voxel_volume_mm3)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("voxel_volume_mm3 must be a positive finite scalar") from exc
    if not np.isfinite(voxel_volume_mm3) or voxel_volume_mm3 <= 0:
        raise ValueError("voxel_volume_mm3 must be a positive finite scalar")
    target_mask = _get_target_mask(target)
    hybrid_mask, no_pde_mask, persistence_mask = _get_prediction_masks(prediction, target_mask.shape)
    if valid_support is None:
        support_mask = np.ones(target_mask.shape, dtype=bool)
    else:
        support_mask = _binary_mask(valid_support, name="valid target support")
        if support_mask.shape != target_mask.shape or not support_mask.any():
            raise ValueError("valid target support must be nonempty and match the target grid")
    if np.any(target_mask & ~support_mask):
        raise ValueError("target foreground falls outside validated target support")
    for name, mask in (("hybrid", hybrid_mask), ("no-PDE", no_pde_mask), ("persistence", persistence_mask)):
        if np.any(mask & ~support_mask):
            raise ValueError(f"{name} prediction falls outside validated target support")
    target_mask = target_mask & support_mask
    hybrid_mask = hybrid_mask & support_mask
    no_pde_mask = no_pde_mask & support_mask
    persistence_mask = persistence_mask & support_mask

    target_voxels = int(target_mask.sum())
    hybrid_voxels = int(hybrid_mask.sum())
    no_pde_voxels = int(no_pde_mask.sum())
    persistence_voxels = int(persistence_mask.sum())
    dice_hybrid = dice_with_empty_rule(hybrid_mask, target_mask)
    dice_no_pde = dice_with_empty_rule(no_pde_mask, target_mask)
    dice_persistence = dice_with_empty_rule(persistence_mask, target_mask)
    paired_delta = float(dice_hybrid - dice_no_pde)
    result = {
        "dice_hybrid": float(dice_hybrid),
        "dice_no_pde_matched_volume": float(dice_no_pde),
        "dice_persistence": float(dice_persistence),
        "paired_delta_dice_hybrid_minus_no_pde": paired_delta,
        "target_voxels": target_voxels,
        "hybrid_prediction_voxels": hybrid_voxels,
        "no_pde_prediction_voxels": no_pde_voxels,
        "persistence_prediction_voxels": persistence_voxels,
        "target_volume_mm3": float(target_voxels * voxel_volume_mm3),
        "hybrid_prediction_volume_mm3": float(hybrid_voxels * voxel_volume_mm3),
        "no_pde_prediction_volume_mm3": float(no_pde_voxels * voxel_volume_mm3),
        "persistence_prediction_volume_mm3": float(persistence_voxels * voxel_volume_mm3),
        "target_empty": bool(target_voxels == 0),
        "valid_support_voxels": int(support_mask.sum()),
        "hybrid_prediction_empty": bool(hybrid_voxels == 0),
        "no_pde_prediction_empty": bool(no_pde_voxels == 0),
        "persistence_prediction_empty": bool(persistence_voxels == 0),
        "ensemble_inclusion_threshold": 0.5,
        "empty_mask_dice_rule": "both_empty=1; one_empty=0",
    }
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
            valid_support=valid_support,
        )

    return score
