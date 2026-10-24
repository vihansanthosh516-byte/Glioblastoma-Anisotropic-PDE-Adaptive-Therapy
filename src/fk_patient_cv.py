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


def _normalize_series(records, *, minimum_history: int) -> tuple[PatientSeries, ...]:
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
            if index < len(masks) - 1 and (not mask.any() or mask.all()):
                raise ValueError("history masks must be nonempty and not fill the computational domain")
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
                         minimum_history: int = 3) -> dict:
    """Run rolling-origin forecasts with patient-level outer folds and target-isolated callbacks.

    ``fit_nuisance`` is called once per fold with training patients only. ``forecast`` receives one patient's
    history, the scheduled forecast day, and the fitted nuisance object, never its target mask. ``score`` receives
    the forecast and that target.
    If ``forecast`` raises :class:`ForecastFailure`, the required fallback is run on the same history only.
    Inner selection and cross-fitting used by ``fit_nuisance`` remain the caller's responsibility.
    """
    if isinstance(n_folds, (bool, np.bool_)) or int(n_folds) != n_folds or n_folds < 2:
        raise ValueError("n_folds must be an integer of at least two")
    if (isinstance(minimum_history, (bool, np.bool_)) or int(minimum_history) != minimum_history
            or minimum_history < 3):
        raise ValueError("minimum_history must be an integer of at least three scans")
    n_folds, minimum_history = int(n_folds), int(minimum_history)
    series = _normalize_series(records, minimum_history=minimum_history)
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
                history = PatientHistory(
                    patient_id=record.patient_id,
                    target_index=target_index,
                    forecast_day=float(record.days[target_index]),
                    days=_readonly_array(record.days[:target_index], dtype=float, name="history days"),
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
