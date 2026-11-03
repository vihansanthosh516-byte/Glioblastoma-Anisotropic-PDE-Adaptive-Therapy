"""Paired patient bootstrap with single-step max-T intervals (draft amendment A-2026-11-02-T §5e). Not adopted.

Choices frozen in the draft before coding: B = 10,000, seed 20261110, the same resampled patients for every
contrast, studentized by the bootstrap SD (ddof 1), c = 95th percentile of the max statistic, interval
mean +/- c * se. Pass only if every adjusted lower bound is strictly above zero; missing arms -> incomplete.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

N_BOOTSTRAP = 10_000
BOOTSTRAP_SEED = 20261110
LEVEL = 0.95


def max_t_pass_test(dice: Mapping[str, Mapping[str, float | None]], comparators: Sequence[str], *,
                    candidate: str = "candidate", n_bootstrap: int = N_BOOTSTRAP, seed: int = BOOTSTRAP_SEED) -> dict:
    ids = sorted(dice[candidate])
    missing = {k: [p for p in ids if dice.get(k, {}).get(p) is None] for k in comparators}
    missing[candidate] = [p for p in ids if dice[candidate][p] is None]
    if any(missing.values()) or any(set(dice.get(k, {})) != set(ids) for k in comparators):
        return {"status": "incomplete", "reason": "missing arm for at least one scored patient",
                "missing": {k: v for k, v in missing.items() if v}}
    if len(ids) < 2:
        return {"status": "incomplete", "reason": "fewer than two patients"}
    d = np.array([[dice[candidate][p] - dice[k][p] for k in comparators] for p in ids], float)
    means = d.mean(axis=0)
    rng = np.random.default_rng(int(seed))
    idx = rng.integers(0, len(ids), size=(int(n_bootstrap), len(ids)))
    boot = d[idx].mean(axis=1)                       # (B, K), same patients for every contrast
    se = boot.std(axis=0, ddof=1)
    if not np.isfinite(se).all() or np.any(se <= 0):
        return {"status": "incomplete", "reason": "zero or non-finite bootstrap SE", "mean_difference": means.tolist(),
                "se": se.tolist()}
    t_max = np.max(np.abs(boot - means) / se, axis=1)
    c = float(np.quantile(t_max, LEVEL))
    lower, upper = means - c * se, means + c * se
    rows = {k: {"mean_difference": float(m), "se": float(s), "adjusted_interval": [float(lo), float(hi)]}
            for k, m, s, lo, hi in zip(comparators, means, se, lower, upper)}
    return {"status": "pass" if bool(np.all(lower > 0)) else "not passed", "n_patients": len(ids),
            "critical_value": c, "contrasts": rows, "n_bootstrap": int(n_bootstrap), "seed": int(seed),
            "bootstrap_unit": "patient", "adjustment": "single-step max-T, studentized by bootstrap SD"}
