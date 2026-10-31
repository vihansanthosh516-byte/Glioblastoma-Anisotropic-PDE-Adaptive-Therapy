"""Fixed-growth mask comparator for the Track B day-90 forecast.

The forecast dilates the last history tumor-core mask by one physical radius and clips it to an approved
valid-support mask. The radius is chosen per outer fold from outer-training patients only. The comparator has
one free choice, so it uses no inner folds: the radius with the best training patient-macro Dice wins, and ties
go to the smaller radius. Radius zero is copy-last clipped to support. This module does not decide valid support;
the caller must pass a source-backed support mask, and scoring stays blocked until that gate passes.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import distance_transform_edt

from src.fk_mu_scoring import dice_with_empty_rule

RADIUS_GRID_MM = tuple(float(r) for r in range(0, 22, 2))
TIE_TOLERANCE = 1e-9


@dataclass(frozen=True)
class FixedGrowthCase:
    patient_id: str
    last_mask: np.ndarray
    support: np.ndarray
    spacing_mm: tuple[float, float, float]
    target_mask: np.ndarray | None = None


def _mask(value, *, name: str) -> np.ndarray:
    array = np.asarray(value)
    if array.ndim != 3:
        raise ValueError(f"{name} must be a 3D binary mask")
    if array.dtype != np.bool_ and not np.isin(array, (0, 1)).all():
        raise ValueError(f"{name} must contain only binary values")
    return array.astype(bool, copy=False)


def _spacing(value) -> tuple[float, float, float]:
    spacing = tuple(float(s) for s in value)
    if len(spacing) != 3 or not all(np.isfinite(s) and s > 0 for s in spacing):
        raise ValueError("spacing_mm must be three positive finite values")
    return spacing


def fixed_growth_forecast(last_mask, support, *, spacing_mm, radius_mm: float) -> np.ndarray:
    """Return {v in support : distance(v, last core) <= radius_mm}, distance in mm on the voxel grid."""
    last = _mask(last_mask, name="last_mask")
    valid = _mask(support, name="support")
    if last.shape != valid.shape:
        raise ValueError("last_mask and support must have the same shape")
    radius = float(radius_mm)
    if not np.isfinite(radius) or radius < 0:
        raise ValueError("radius_mm must be finite and nonnegative")
    if not last.any():
        raise ValueError("empty last history core is ineligible for the fixed-growth comparator")
    if radius == 0:
        return last & valid
    distance = distance_transform_edt(~last, sampling=_spacing(spacing_mm))
    return (distance <= radius + 1e-9) & valid


def select_radius(training_cases: Sequence[FixedGrowthCase], *, radius_grid_mm=RADIUS_GRID_MM) -> dict:
    """Pick the radius with the best patient-macro Dice on training cases; ties go to the smaller radius."""
    grid = sorted(float(r) for r in radius_grid_mm)
    if not grid or grid[0] != 0.0 or len(set(grid)) != len(grid):
        raise ValueError("radius grid must be unique, nonnegative and include 0")
    if not training_cases:
        raise ValueError("radius selection needs at least one training case")
    ids = [case.patient_id for case in training_cases]
    if len(set(ids)) != len(ids):
        raise ValueError("each training patient must appear once")
    mean_dice = []
    for radius in grid:
        scores = []
        for case in training_cases:
            if case.target_mask is None:
                raise ValueError(f"training case {case.patient_id} has no target mask")
            forecast = fixed_growth_forecast(case.last_mask, case.support, spacing_mm=case.spacing_mm,
                                             radius_mm=radius)
            # Target core outside support stays in the target, so it counts as a miss.
            scores.append(dice_with_empty_rule(forecast, _mask(case.target_mask, name="target_mask")))
        mean_dice.append(float(np.mean(scores)))
    best = max(mean_dice)
    chosen = next(r for r, d in zip(grid, mean_dice) if d >= best - TIE_TOLERANCE)
    return {"radius_mm": chosen, "grid_mm": grid, "training_mean_dice": mean_dice,
            "n_training": len(training_cases), "boundary_hit": chosen == grid[-1] and len(grid) > 1}


def run_outer_fixed_growth(cases: Sequence[FixedGrowthCase], fold_by_patient: Mapping[str, int], *,
                           radius_grid_mm=RADIUS_GRID_MM) -> dict:
    """Select a radius on each outer fold's training patients, then forecast that fold's test patients.

    Test-patient targets are never passed to selection. Forecast masks are built before any test target is read.
    """
    ids = [case.patient_id for case in cases]
    if len(set(ids)) != len(ids):
        raise ValueError("each patient must appear once")
    if set(fold_by_patient) != set(ids):
        raise ValueError("fold_by_patient must cover exactly the case patients")
    rows, selections = [], {}
    for fold in sorted(set(fold_by_patient.values())):
        train = [c for c in cases if fold_by_patient[c.patient_id] != fold]
        test = [c for c in cases if fold_by_patient[c.patient_id] == fold]
        selection = select_radius(train, radius_grid_mm=radius_grid_mm)
        selections[fold] = selection
        forecasts = {c.patient_id: fixed_growth_forecast(c.last_mask, c.support, spacing_mm=c.spacing_mm,
                                                         radius_mm=selection["radius_mm"]) for c in test}
        for case in test:
            row = {"patient_id": case.patient_id, "fold": fold, "radius_mm": selection["radius_mm"],
                   "equals_copy_last": selection["radius_mm"] == 0.0, "forecast": forecasts[case.patient_id]}
            if case.target_mask is not None:
                target = _mask(case.target_mask, name="target_mask")
                valid = _mask(case.support, name="support")
                row["dice"] = dice_with_empty_rule(forecasts[case.patient_id], target)
                row["target_voxels_outside_support"] = int(np.count_nonzero(target & ~valid))
            rows.append(row)
    return {"rows": rows, "selections": selections}


def add_fixed_growth_arm(fit_nuisance, forecast, fallback_forecast, score, *, domain_by_patient: Mapping,
                         spacing_mm, valid_support_by_patient: Mapping | None = None,
                         radius_grid_mm=RADIUS_GRID_MM):
    """Wrap patient-CV callbacks so every forecast also carries the fixed-growth comparator.

    ``fit_nuisance`` receives outer-training records only (``run_outer_patient_cv`` contract). The radius is chosen
    there from each training patient's last history mask (visit 2) and target (visit 3), clipped to that patient's
    origin-only domain. Training patients with an empty last history core are ineligible for selection. Test
    forecasts use only the test history and the origin-only domain. Empty test history gets copy-last (empty).
    The score wrapper adds plain Dice for fixed growth after the base score is complete. If target support is
    given, it applies the base scorer's rule: a forecast outside validated support is an error, not a clip.
    """
    spacing = _spacing(spacing_mm)

    def fit(training_records):
        training_records = tuple(training_records)
        cases = []
        for record in training_records:
            last = _mask(record.masks[-2], name="training last history mask")
            if not last.any():
                continue
            cases.append(FixedGrowthCase(str(record.patient_id), last, domain_by_patient[record.patient_id],
                                         spacing, _mask(record.masks[-1], name="training target mask")))
        selection = select_radius(cases, radius_grid_mm=radius_grid_mm)
        selection["n_training_empty_history_excluded"] = len(training_records) - len(cases)
        return {"base": fit_nuisance(training_records), "fixed_growth": selection}

    def _with_fixed_growth(prediction, history, wrapped):
        last = _mask(history.masks[-1], name="history last mask")
        domain = _mask(domain_by_patient[history.patient_id], name="origin domain")
        if last.any():
            mask = fixed_growth_forecast(last, domain, spacing_mm=spacing,
                                         radius_mm=wrapped["fixed_growth"]["radius_mm"])
            status = "forecast" if wrapped["fixed_growth"]["radius_mm"] > 0 else "equals_copy_last"
        else:
            mask, status = last.copy(), "empty_history_copy_last"
        if not isinstance(prediction, Mapping):
            prediction = {"fallback_masks": {"hybrid": prediction, "no_pde": prediction}, "history_mask": last}
        return {**prediction, "fixed_growth_mask": mask, "fixed_growth_status": status,
                "fixed_growth_radius_mm": wrapped["fixed_growth"]["radius_mm"]}

    def wrapped_forecast(history, wrapped):
        return _with_fixed_growth(forecast(history, wrapped["base"]), history, wrapped)

    def wrapped_fallback(history, wrapped):
        return _with_fixed_growth(fallback_forecast(history, wrapped["base"]), history, wrapped)

    def wrapped_score(prediction, target):
        base = {k: v for k, v in prediction.items()
                if k not in ("fixed_growth_mask", "fixed_growth_status", "fixed_growth_radius_mm")}
        scores = dict(score(base, target))
        target_mask = _mask(target.mask if hasattr(target, "mask") else target["mask"], name="target")
        forecast_mask = _mask(prediction["fixed_growth_mask"], name="fixed-growth mask")
        if valid_support_by_patient is not None:
            support = _mask(valid_support_by_patient[target.patient_id], name="valid target support")
            if np.any(forecast_mask & ~support):
                raise ValueError("fixed-growth prediction falls outside validated target support")
            target_mask = target_mask & support
        dice = dice_with_empty_rule(forecast_mask, target_mask)
        scores.update({"dice_fixed_growth": float(dice),
                       "fixed_growth_status": prediction["fixed_growth_status"],
                       "fixed_growth_radius_mm": float(prediction["fixed_growth_radius_mm"])})
        if "dice_hybrid" in scores:
            scores["paired_delta_dice_hybrid_minus_fixed_growth"] = float(scores["dice_hybrid"] - dice)
        return scores

    return fit, wrapped_forecast, wrapped_fallback, wrapped_score
