"""Synthetic-only two-history, time-aware NON-PDE shape comparator for Track B development.

Rule: signed-distance (SD) maps of the earlier history mask h1 and the origin mask h2 are extrapolated,
    SD_f = SD_h2 + f * clip(SD_h2 - SD_h1, +-MAX_SD_DELTA_MM),   f = min(w * dt_forecast / dt_history, MAX_FACTOR),
    w = dt_history / (dt_history + SHRINK_TAU_DAYS)   (shrinks the extrapolation for short history gaps).
For each volume draw, the top-k voxels of SD_f inside the support S match that draw's volume exactly. Draws and
boundary seeds can be shared with the PDE-hybrid candidate. No patient I/O; not wired into any runner.

PROPOSED (not adopted): every constant is a development default, untuned, and the strength of this control is
unverified. Before adoption a development sensitivity grid over tau, cap and clip, compared with copy-last and the
last-mask signed-distance arm, is required. Note f = min(T / (g + 7), MAX_FACTOR) (g = history gap days, T = forecast
days): g cancels in w * T / g. Ties among equal SD_f values break by flat index (an arbitrary, grid-order
dependent rule). The control is not in the max-T comparator family; an amendment must define that family. If the
control fails for a patient (a ForecastFailure-style event), its fallback is not defined here and must be defined
before wiring. The boundary covariance used with apply_boundary=True was calibrated on the PDE arm.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
from scipy import ndimage

from src.fk_grid_likelihood import volume_match_concentration
from src.fk_hybrid_forecast import _perturb_boundary_or_keep

# Proposed (not adopted) development defaults. The explicit settings object lets a future
# outer-training-only selector evaluate the proposed grid without changing module globals.
MAX_FACTOR = 4.0
MIN_HISTORY_GAP_DAYS = 1.0
SHRINK_TAU_DAYS = 7.0
MAX_SD_DELTA_MM = 30.0


@dataclass(frozen=True)
class TwoHistoryShapeSettings:
    """Validated settings for the synthetic-only two-history shape prototype."""

    max_factor: float = MAX_FACTOR
    min_history_gap_days: float = MIN_HISTORY_GAP_DAYS
    shrink_tau_days: float = SHRINK_TAU_DAYS
    max_sd_delta_mm: float = MAX_SD_DELTA_MM

    def __post_init__(self):
        values = np.asarray(
            [self.max_factor, self.min_history_gap_days, self.shrink_tau_days, self.max_sd_delta_mm],
            dtype=float,
        )
        if not np.isfinite(values).all():
            raise ValueError("two-history settings must be finite")
        if self.max_factor <= 0:
            raise ValueError("max_factor must be positive; use f=0 for the last-mask rule")
        if self.min_history_gap_days < 0 or self.shrink_tau_days < 0:
            raise ValueError("history-gap threshold and shrinkage time must be nonnegative")
        if self.max_sd_delta_mm <= 0:
            raise ValueError("max_sd_delta_mm must be positive")

    def as_dict(self) -> dict[str, float]:
        """Return a stable, JSON-friendly representation for result records."""
        return {
            "max_factor": float(self.max_factor),
            "min_history_gap_days": float(self.min_history_gap_days),
            "shrink_tau_days": float(self.shrink_tau_days),
            "max_sd_delta_mm": float(self.max_sd_delta_mm),
        }


DEFAULT_SETTINGS = TwoHistoryShapeSettings()


def _spacing3(spacing_mm) -> np.ndarray:
    s = np.broadcast_to(np.asarray(spacing_mm, dtype=float), (3,)).copy()
    if not np.isfinite(s).all() or np.any(s <= 0):
        raise ValueError("spacing_mm must contain three positive finite values")
    return s


def signed_distance(mask, spacing_mm) -> np.ndarray:
    mask = np.asarray(mask, dtype=bool)
    s = _spacing3(spacing_mm)
    return ndimage.distance_transform_edt(mask, sampling=s) - ndimage.distance_transform_edt(~mask, sampling=s)


def extrapolation_factor(
    gap_days: float,
    forecast_days: float,
    *,
    settings: TwoHistoryShapeSettings = DEFAULT_SETTINGS,
) -> float:
    """Proposed factor f = min(T / (g + tau), cap); zero for tiny history gaps."""
    gap_days, forecast_days = float(gap_days), float(forecast_days)
    if not np.isfinite(gap_days) or not np.isfinite(forecast_days) or forecast_days <= 0:
        raise ValueError("forecast interval must be positive and finite; history gap must be finite")
    if gap_days < 0:
        raise ValueError("history gap must be nonnegative")
    if gap_days == 0 or gap_days < settings.min_history_gap_days:
        return 0.0
    weight = gap_days / (gap_days + settings.shrink_tau_days)
    return float(min(weight * forecast_days / gap_days, settings.max_factor))


def forecast_signed_distance(
    h1,
    h2,
    gap_days,
    forecast_days,
    spacing_mm=1.0,
    *,
    settings: TwoHistoryShapeSettings = DEFAULT_SETTINGS,
) -> tuple[np.ndarray, dict]:
    """Return SD_f and details. Empty/full h1 falls back to the last-mask rule (f=0); empty/full h2 raises."""
    h1, h2 = np.asarray(h1, dtype=bool), np.asarray(h2, dtype=bool)
    if h2.ndim != 3 or h1.shape != h2.shape:
        raise ValueError("h1 and h2 must be 3D arrays on one grid")
    if not h2.any() or h2.all():
        raise ValueError("origin mask h2 must be nonempty and non-full")
    sd2 = signed_distance(h2, spacing_mm)
    h1_unusable = (not h1.any()) or bool(h1.all())
    requested_factor = extrapolation_factor(gap_days, forecast_days, settings=settings)
    f = 0.0 if h1_unusable else requested_factor
    setting_record = settings.as_dict()
    if f == 0.0:
        return sd2, {"factor": 0.0, "h1_unusable_fallback": h1_unusable,
                     "tiny_gap_fallback": (not h1_unusable), "settings": setting_record}
    delta = np.clip(
        sd2 - signed_distance(h1, spacing_mm),
        -settings.max_sd_delta_mm,
        settings.max_sd_delta_mm,
    )
    return sd2 + f * delta, {
        "factor": f,
        "h1_unusable_fallback": False,
        "tiny_gap_fallback": False,
        "settings": setting_record,
    }


def candidate_boundary_seeds(seed: int, n_draws: int) -> list[int]:
    """Per-draw boundary seeds exactly as hybrid_posterior_predictive_forecast derives them."""
    _, boundary_seed = np.random.SeedSequence(seed).spawn(2)
    return [int(s.generate_state(1)[0]) for s in boundary_seed.spawn(int(n_draws))]


def two_history_shape_draws(h1, h2, gap_days, forecast_days, target_voxels, *, domain,
                            spacing_mm=1.0, boundary_covariance: dict | None = None,
                            boundary_seeds=None, seed: int = 0, n_anchors: int = 256,
                            apply_boundary: bool = True,
                            settings: TwoHistoryShapeSettings = DEFAULT_SETTINGS) -> dict:
    """Matched-volume draws from the two-history rule inside ``domain`` (the shared support S).

    target_voxels: the candidate's requested_volume_draws_voxels (same draws). boundary_seeds: pass
    candidate_boundary_seeds(candidate_seed, n) for common random numbers; default derives them from ``seed``.
    apply_boundary=False gives the deterministic pre-boundary matched-volume draws.
    """
    domain = np.asarray(domain, dtype=bool)
    h2b = np.asarray(h2, dtype=bool)
    if domain.shape != h2b.shape or np.any(h2b & ~domain):
        raise ValueError("domain must match and contain h2")
    targets = np.asarray(target_voxels)
    if (targets.ndim != 1 or len(targets) < 1 or not np.issubdtype(targets.dtype, np.integer)
            or np.any(targets < 1) or np.any(targets >= int(domain.sum()))):
        raise ValueError("target_voxels must be positive integer counts below domain capacity")
    spacing = _spacing3(spacing_mm)
    sd_f, details = forecast_signed_distance(
        h1, h2, gap_days, forecast_days, spacing, settings=settings
    )
    if apply_boundary:
        if boundary_covariance is None:
            raise ValueError("boundary_covariance is required when apply_boundary is True")
        seeds = (candidate_boundary_seeds(seed, len(targets)) if boundary_seeds is None
                 else [int(s) for s in boundary_seeds])
        if len(seeds) != len(targets):
            raise ValueError("need one boundary seed per volume draw")
    pre, post, fallbacks = [], [], 0
    for i, k in enumerate(targets):
        ranked = volume_match_concentration(sd_f, domain, target_voxels=int(k))
        pre.append(ranked)
        if not apply_boundary:
            continue
        draw, collapsed = _perturb_boundary_or_keep(
            ranked, boundary_covariance=boundary_covariance, spacing_mm=spacing, seed=seeds[i],
            n_anchors=n_anchors, domain=domain, target_voxels=int(k))
        fallbacks += int(collapsed)
        post.append(draw & domain)
    pre = np.stack(pre, axis=0)
    out = pre if not apply_boundary else np.stack(post, axis=0)
    if not np.array_equal(out.reshape(len(targets), -1).sum(axis=1), targets):
        raise RuntimeError("two-history draws must preserve every requested voxel count")
    return {"mask_draws": out, "pre_boundary_mask_draws": pre, "boundary_fallback_count": fallbacks,
            "boundary_applied": bool(apply_boundary), **details}


def consensus_mask(draws, domain, threshold: float = 0.5):
    """Inclusion probability >= threshold inside domain (the candidate scorer's combination rule)."""
    p = np.asarray(draws).mean(axis=0)
    return (p >= threshold) & np.asarray(domain, dtype=bool), p
