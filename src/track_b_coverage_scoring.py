"""Coverage-aware Track B scoring (amendment paper/trackB_coverage_scoring_amendment_2026-11-10.md).

Support is built only from scan geometry, never from labels. Locations outside the
common support of all required visits are unobserved: excluded, not negative. Every
method in a window is scored on the same support. Missing or invalid coverage fails
closed. Historical scorers are unchanged.

Boundary enforcement: scoring functions accept only an `ObservableSupport` returned by
`common_support()`, which accepts only `ScanCoverage` objects built by
`coverage_from_geometry()`. Raw Boolean arrays are rejected. Masks are copied and made
read-only at creation.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from src.fk_mu_scoring import dice_with_empty_rule
from src.track_b_transition_estimand import geometric_coverage, summarize_window

_GEOMETRY_TOKEN = object()
_SUPPORT_TOKEN = object()


def _frozen_copy(mask) -> np.ndarray:
    arr = np.array(mask, dtype=bool, copy=True)
    arr.setflags(write=False)
    return arr


@dataclass(frozen=True)
class ScanCoverage:
    """Boolean coverage of one visit on the scoring grid, derived from scan geometry."""
    mask: np.ndarray
    visit: str
    _token: object = None

    def __post_init__(self):
        if self._token is not _GEOMETRY_TOKEN:
            raise TypeError("ScanCoverage must be built with coverage_from_geometry")
        object.__setattr__(self, "mask", _frozen_copy(self.mask))


@dataclass(frozen=True)
class ObservableSupport:
    """Common support of all required visits, with its geometry provenance."""
    mask: np.ndarray
    visits: tuple
    provenance: str
    _token: object = None

    def __post_init__(self):
        if self._token is not _SUPPORT_TOKEN:
            raise TypeError("ObservableSupport must be built with common_support")
        object.__setattr__(self, "mask", _frozen_copy(self.mask))


def coverage_from_geometry(visit: str, moving_shape, output_shape, output_to_input) -> ScanCoverage:
    """Coverage from grid shapes and the exact resampling transform; no image or label input."""
    mask = geometric_coverage(moving_shape, output_shape, output_to_input)
    return ScanCoverage(mask=mask, visit=str(visit), _token=_GEOMETRY_TOKEN)


def common_support(coverages, required_visits) -> ObservableSupport:
    """AND of coverage over every required visit. Fails closed."""
    required = tuple(required_visits)
    if not required:
        raise ValueError("at least one required visit")
    if len(set(required)) != len(required):
        raise ValueError("required visits must be distinct")
    by_visit = {}
    for c in coverages:
        if not isinstance(c, ScanCoverage):
            raise TypeError("coverage must be ScanCoverage objects")
        if c.visit in by_visit:
            raise ValueError(f"duplicate coverage entry for visit {c.visit!r}")
        by_visit[c.visit] = c
    missing = [v for v in required if v not in by_visit]
    if missing:
        raise ValueError(f"no scan coverage for required visits: {missing}")
    masks = [by_visit[v].mask for v in required]
    shape = masks[0].shape
    for m in masks:
        if m.dtype != np.bool_ or m.shape != shape or m.ndim != 3:
            raise ValueError("coverage masks must be same-shape 3D boolean arrays")
    support = np.logical_and.reduce(masks)
    if not support.any():
        raise ValueError("common observable support is empty")
    return ObservableSupport(mask=support, visits=required, provenance="scan_geometry", _token=_SUPPORT_TOKEN)


def _require_support(support, shape) -> np.ndarray:
    if not isinstance(support, ObservableSupport) or support.provenance != "scan_geometry":
        raise TypeError("support must be an ObservableSupport from common_support; raw masks are rejected")
    if support.mask.shape != tuple(shape):
        raise ValueError("support shape does not match the scored arrays")
    return support.mask


def _support_report(support):
    mask = support.mask
    return {"support_voxels": int(mask.sum()), "support_fraction": float(mask.mean()),
            "support_visits": list(support.visits)}


def transition_window_on_support(state, bands, target, support, *, min_voxels=10):
    """Transition contrast using only risk voxels inside support. None if no band is eligible."""
    state, bands, target = map(np.asarray, (state, bands, target))
    mask = _require_support(support, state.shape)
    restricted = np.where(mask, state, 0).astype(state.dtype)
    target_on_support = np.where(mask, target, 0)
    result = summarize_window(restricted, bands, target_on_support, min_voxels=min_voxels)
    if result is None:
        return None
    return {**result, **_support_report(support)}


def forecast_dice_on_support(predictions: dict, target, support) -> dict:
    """Dice for every method on one common support; outside support ignored for all."""
    target = np.asarray(target, bool)
    mask = _require_support(support, target.shape)
    if not predictions:
        raise ValueError("no methods to score")
    scores = {}
    for name, pred in predictions.items():
        pred = np.asarray(pred, bool)
        if pred.shape != target.shape:
            raise ValueError(f"prediction {name} has the wrong shape")
        scores[name] = dice_with_empty_rule(pred & mask, target & mask)
    return {"dice": scores, **_support_report(support)}
