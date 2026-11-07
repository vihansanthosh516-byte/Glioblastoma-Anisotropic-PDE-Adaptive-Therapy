"""World-level summaries for paired simulated policy differences.

This development helper does not define an adopted scoring protocol.
"""
from __future__ import annotations

import numpy as np


def paired_world_summary(candidate, comparator, world_ids, *, resamples, seed):
    """Average patients within worlds, then bootstrap independent worlds equally.

    Inputs are matching two-dimensional world-by-patient arrays. Patient positions
    must represent paired simulated patients. World identifiers must be unique.
    """
    candidate = np.asarray(candidate, dtype=float)
    comparator = np.asarray(comparator, dtype=float)
    if candidate.ndim != 2 or candidate.shape != comparator.shape:
        raise ValueError("paired outcomes must be matching world-by-patient arrays")
    n_worlds, n_patients = candidate.shape
    if n_worlds < 2 or n_patients < 1:
        raise ValueError("at least two worlds and one patient per world are required")
    if not np.isfinite(candidate).all() or not np.isfinite(comparator).all():
        raise ValueError("paired outcomes must be finite")
    ids = tuple(world_ids)
    if len(ids) != n_worlds or any(not isinstance(x, str) or not x.strip() for x in ids) or len(set(ids)) != n_worlds:
        raise ValueError("world identifiers must be unique nonempty strings")
    if isinstance(resamples, (bool, np.bool_)) or not isinstance(resamples, (int, np.integer)) or resamples < 2:
        raise ValueError("resamples must be an integer of at least 2")
    differences = (candidate - comparator).mean(axis=1)
    rng = np.random.default_rng(seed)
    boot = differences[rng.integers(0, n_worlds, size=(resamples, n_worlds))].mean(axis=1)
    return {
        "n_worlds": n_worlds,
        "n_patients_per_world": n_patients,
        "world_mean_differences": differences,
        "equal_world_mean_difference": float(differences.mean()),
        "bootstrap_percentile_95": np.quantile(boot, [0.025, 0.975]),
        "bootstrap_unit": "world",
        "resamples": int(resamples),
        "seed": seed,
    }
