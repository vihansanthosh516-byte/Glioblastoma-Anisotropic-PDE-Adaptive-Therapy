"""Outer patient-fold, rolling-origin evaluation harness for Track B.

The nuisance-fitting callback sees outer-training patients only. The forecast callback receives history only;
the target mask is created separately and passed only to the scoring callback. This harness does not implement
the nuisance model or prove that a callback internally avoids leakage.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

import numpy as np


class ForecastFailure(RuntimeError):
    """Raise only for a forecast that should use the prespecified history-only fallback."""


@dataclass(frozen=True)
class PatientSeries:
    patient_id: str
    days: object
    masks: Sequence[object]


@dataclass(frozen=True)
class PatientHistory:
    patient_id: str
    target_index: int
    forecast_day: float
    days: np.ndarray
    masks: tuple[np.ndarray, ...]


@dataclass(frozen=True)
class PatientTarget:
    patient_id: str
    target_index: int
    day: float
    mask: np.ndarray


def _readonly_array(value, *, dtype, name: str) -> np.ndarray:
    result = np.array(value, dtype=dtype, copy=True)
    if result.ndim == 0:
        raise ValueError(f"{name} must be an array")
    result.setflags(write=False)
    return result


def _normalize_series(records, *, minimum_history: int,
                      allow_empty_history_masks: bool = False) -> tuple[PatientSeries, ...]:
    normalized = []
    seen = set()
    for record in records:
        if not isinstance(record, PatientSeries):
            raise ValueError("records must contain PatientSeries values")
        patient_id = str(record.patient_id).strip()
        if not patient_id:
            raise ValueError("patient IDs must be nonempty")
        if patient_id in seen:
            raise ValueError("duplicate patient in series records")
        seen.add(patient_id)
        days = np.asarray(record.days, dtype=float)
        masks = tuple(record.masks)
        if (days.ndim != 1 or len(days) != len(masks) or len(days) < minimum_history + 1
                or not np.isfinite(days).all() or np.any(np.diff(days) <= 0)):
            raise ValueError("each patient needs finite strictly increasing days and enough scans for a target")
        checked_masks = []
        shape = None
        for index, raw_mask in enumerate(masks):
            raw = np.asarray(raw_mask)
            if raw.ndim != 3 or not (raw.dtype == np.bool_ or np.isin(raw, (0, 1)).all()):
                raise ValueError("masks must be binary 3D arrays")
            mask = _readonly_array(raw, dtype=bool, name="mask")
            if shape is None:
                shape = mask.shape
            elif mask.shape != shape:
                raise ValueError("all masks for a patient must share one shape")
            if index < len(masks) - 1:
                if mask.all() or (not mask.any() and not allow_empty_history_masks):
                    raise ValueError("history masks must be nonempty unless the explicit fallback path is enabled, and must not fill the computational domain")
            checked_masks.append(mask)
        normalized.append(PatientSeries(patient_id, _readonly_array(days, dtype=float, name="days"),
                                        tuple(checked_masks)))
    if not normalized:
        raise ValueError("records must contain at least one patient")
    return tuple(normalized)


def _normalize_folds(fold_by_patient: Mapping, patient_ids: set[str], n_folds: int) -> dict[str, int]:
    if not isinstance(fold_by_patient, Mapping):
        raise ValueError("fold_by_patient must map patient IDs to integer folds")
    normalized = {}
    for raw_patient_id, raw_fold in fold_by_patient.items():
        patient_id = str(raw_patient_id).strip()
        if not patient_id or patient_id in normalized:
            raise ValueError("fold map has blank or duplicate normalized patient IDs")
        if isinstance(raw_fold, (bool, np.bool_)):
            raise ValueError("folds must be integers, not booleans")
        try:
            fold = int(raw_fold)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("folds must be integers") from exc
        if fold != raw_fold or fold < 0 or fold >= n_folds:
            raise ValueError("fold value is outside the declared outer folds")
        normalized[patient_id] = fold
    if set(normalized) != patient_ids:
        raise ValueError("fold map must cover exactly the patient series records")
    if set(normalized.values()) != set(range(n_folds)):
        raise ValueError("every outer fold must contain at least one patient")
    return normalized


def run_outer_patient_cv(records: Sequence[PatientSeries], fold_by_patient: Mapping,
                         fit_nuisance: Callable, forecast: Callable, score: Callable, *,
                         fallback_forecast: Callable | None = None, n_folds: int = 5,
                         minimum_history: int = 2,
                         forecast_horizon_days: float | None = None,
                         allow_empty_history_masks: bool = False) -> dict:
    """Run rolling-origin forecasts with patient-level outer folds and target-isolated callbacks.

    ``fit_nuisance`` is called once per fold with training patients only. ``forecast`` receives one patient's
    history, the scheduled forecast day, and the fitted nuisance object, never its target mask. ``score`` receives
    the forecast and that target.
    ``minimum_history`` counts observed masks before each target; two masks allow the planned first-two-to-third
    forecast. If ``forecast`` raises :class:`ForecastFailure`, the required fallback uses same history only.
    ``forecast_horizon_days`` optionally schedules the forecast at that fixed offset after the latest history day;
    the actual target day remains available only to the scoring callback. ``None`` preserves the historical behavior
    of scheduling at the observed target day. For a one-target-per-patient protocol, callers must pass each
    patient's first two history visits plus only the immediate third target (truncate the series to three visits).
    Inner selection and cross-fitting used by ``fit_nuisance`` remain the caller's responsibility.
    """
    if isinstance(n_folds, (bool, np.bool_)) or int(n_folds) != n_folds or n_folds < 2:
        raise ValueError("n_folds must be an integer of at least two")
    if (isinstance(minimum_history, (bool, np.bool_)) or int(minimum_history) != minimum_history
            or minimum_history < 2):
        raise ValueError("minimum_history must be an integer of at least two scans")
    if forecast_horizon_days is not None:
        if (isinstance(forecast_horizon_days, (bool, np.bool_))
                or not np.isscalar(forecast_horizon_days)
                or isinstance(forecast_horizon_days, (str, bytes))):
            raise ValueError("forecast_horizon_days must be a finite positive numeric scalar")
        try:
            forecast_horizon_days = float(forecast_horizon_days)
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("forecast_horizon_days must be a finite positive numeric scalar") from exc
        if not np.isfinite(forecast_horizon_days) or forecast_horizon_days <= 0:
            raise ValueError("forecast_horizon_days must be a finite positive numeric scalar")
    n_folds, minimum_history = int(n_folds), int(minimum_history)
    if not isinstance(allow_empty_history_masks, (bool, np.bool_)):
        raise ValueError("allow_empty_history_masks must be a boolean")
    series = _normalize_series(
        records, minimum_history=minimum_history,
        allow_empty_history_masks=bool(allow_empty_history_masks),
    )
    folds = _normalize_folds(fold_by_patient, {record.patient_id for record in series}, n_folds)

    rows = []
    fold_audit = []
    for fold in range(n_folds):
        test_records = tuple(sorted((record for record in series if folds[record.patient_id] == fold),
                                    key=lambda record: record.patient_id))
        train_records = tuple(record for record in series if folds[record.patient_id] != fold)
        if not train_records or not test_records:
            raise ValueError("outer fold has an empty training or test patient set")
        train_ids = {record.patient_id for record in train_records}
        test_ids = {record.patient_id for record in test_records}
        if train_ids & test_ids:
            raise RuntimeError("patient crossed the outer train/test boundary")
        nuisance = fit_nuisance(train_records)
        fold_rows = []
        fallback_count = 0
        for record in test_records:
            for target_index in range(minimum_history, len(record.days)):
                history_days = _readonly_array(record.days[:target_index], dtype=float, name="history days")
                scheduled_forecast_day = (
                    float(record.days[target_index]) if forecast_horizon_days is None
                    else float(history_days[-1] + forecast_horizon_days)
                )
                history = PatientHistory(
                    patient_id=record.patient_id,
                    target_index=target_index,
                    forecast_day=scheduled_forecast_day,
                    days=history_days,
                    masks=tuple(record.masks[:target_index]),
                )
                target = PatientTarget(
                    patient_id=record.patient_id,
                    target_index=target_index,
                    day=float(record.days[target_index]),
                    mask=record.masks[target_index],
                )
                try:
                    prediction = forecast(history, nuisance)
                    forecast_status = "ok"
                except ForecastFailure:
                    if fallback_forecast is None:
                        raise RuntimeError("forecast failed but no history-only fallback was supplied")
                    prediction = fallback_forecast(history, nuisance)
                    forecast_status = "fallback"
                    fallback_count += 1
                if prediction is None:
                    raise ValueError("forecast callbacks must return a prediction")
                scores = score(prediction, target)
                if not isinstance(scores, Mapping):
                    raise ValueError("score callback must return a mapping of named scores")
                row = {
                    "patient_id": record.patient_id,
                    "fold": fold,
                    "target_index": target_index,
                    "target_day": target.day,
                    "scheduled_forecast_day": scheduled_forecast_day,
                    "history_scan_count": target_index,
                    "forecast_status": forecast_status,
                    "scores": dict(scores),
                }
                rows.append(row)
                fold_rows.append(row)
        fold_audit.append({
            "fold": fold,
            "n_training_patients": len(train_records),
            "n_test_patients": len(test_records),
            "n_forecasts": len(fold_rows),
            "fallback_count": fallback_count,
        })
    return {
        "rows": rows,
        "fold_audit": fold_audit,
        "scope": "patient-held-out rolling-origin evaluation; target masks reach scoring only",
    }


def run_first_third_90day_patient_cv(
    records: Sequence[PatientSeries], fold_by_patient: Mapping,
    fit_nuisance: Callable, forecast: Callable, score: Callable, *,
    fallback_forecast: Callable | None = None,
    n_folds: int = 5,
) -> dict:
    """Run exactly one day-90 forecast from visits 1–2 to visit 3 per patient.

    The observed third-visit day is available only to ``score``. Empty history masks are admitted so the
    explicitly supplied history-only fallback can handle those patients; the forecast callback must raise
    ``ForecastFailure`` when its PDE history is unusable. This wrapper prevents the protocol call site from
    silently defaulting to the realized target day or scoring later visits.
    """
    if not isinstance(records, Sequence) or isinstance(records, (str, bytes)) or not records:
        raise ValueError("records must be a nonempty sequence of patient series")
    for record in records:
        if not isinstance(record, PatientSeries):
            raise ValueError("records must contain PatientSeries values")
        try:
            n_days, n_masks = len(record.days), len(record.masks)
        except TypeError as exc:
            raise ValueError("first-third CV requires exactly three visits per patient") from exc
        if n_days != 3 or n_masks != 3:
            raise ValueError("first-third CV requires exactly three visits per patient")
    if fallback_forecast is None:
        raise ValueError("the first-third protocol requires an explicit history-only fallback")
    result = run_outer_patient_cv(
        records, fold_by_patient, fit_nuisance, forecast, score,
        fallback_forecast=fallback_forecast,
        n_folds=n_folds,
        minimum_history=2,
        forecast_horizon_days=90.0,
        allow_empty_history_masks=True,
    )
    if len(result["rows"]) != len(records) or any(row["target_index"] != 2 for row in result["rows"]):
        raise RuntimeError("first-third wrapper did not produce exactly one third-visit target per patient")
    result["scope"] = "patient-held-out day-90 forecast from first two visits to third; target masks reach scoring only"
    return result
