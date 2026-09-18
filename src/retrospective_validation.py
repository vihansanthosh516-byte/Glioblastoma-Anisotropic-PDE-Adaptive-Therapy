"""Metrics for retrospective longitudinal tumor-volume validation."""
from __future__ import annotations

from typing import Iterable, Sequence
import numpy as np


def volume_metrics(predicted: Sequence[float], observed: Sequence[float]) -> dict[str, float]:
    pred = np.asarray(predicted, dtype=float)
    obs = np.asarray(observed, dtype=float)
    if pred.shape != obs.shape or pred.size == 0:
        raise ValueError("predicted and observed must be non-empty and equal length")
    err = pred - obs
    result = {"mae_mm3": float(np.mean(np.abs(err))), "rmse_mm3": float(np.sqrt(np.mean(err ** 2)))}
    if pred.size < 2 or np.std(pred) == 0 or np.std(obs) == 0:
        result["pearson_r"] = float("nan")
    else:
        result["pearson_r"] = float(np.corrcoef(pred, obs)[0, 1])
    return result


def dice_score(predicted_mask: np.ndarray, observed_mask: np.ndarray) -> float:
    pred = np.asarray(predicted_mask, dtype=bool)
    obs = np.asarray(observed_mask, dtype=bool)
    if pred.shape != obs.shape:
        raise ValueError("masks must have identical shapes")
    denom = pred.sum() + obs.sum()
    return float(2 * np.logical_and(pred, obs).sum() / denom) if denom else 1.0


def validate_records(records: Iterable[dict]) -> list[dict]:
    """Compute metrics for records containing patient_id, predicted, observed."""
    return [{"patient_id": r["patient_id"], **volume_metrics(r["predicted"], r["observed"])} for r in records]
