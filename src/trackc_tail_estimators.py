"""Tail estimators for Track C outcome arrays (SIMULATION statistics only; runs no simulation).

Input arrays X have shape (R, n): R noise repeats of the same n scenarios. Lower outcome = worse
(time to progression). The tail score is the mean of the worst alpha fraction (CVaR_alpha).

Why this module exists (record of 2026-11-03): in the 167 per-scenario file, the worst-decile score
of the mean over repeats ("plug-in") rises with the number of repeats for rules that read the noisy
scans (wide containment 171.1, 174.5, 176.0 at 1, 2, 4 repeats) and stays flat for a rule that
ignores the scans (front-loaded, 182.2). Selecting the worst scenarios on noisy means picks scenarios
that look bad partly by luck. The estimators below separate three things:

    plugin       CVaR of the per-scenario mean over R repeats. Biased low at small R for closed-loop rules.
    split        choose the worst scenarios on one half of the repeats, score them on the other half.
                 Biased high (the chosen set is not exactly the true worst set). Needs R >= 2.
    joint        CVaR over all R x n single outcomes (what one patient would experience, noise included).
                 Not selection-biased; it measures a different thing (individual realised tail).
    paired       CVaR of (policy - baseline) differences, paired by scenario (and by noise key in the joint form).

The inverse-R fit  plugin(R) = c - a / R  extrapolates the plug-in curve to R -> infinity. It is a model
of the bias, not a theorem; check the residuals before using c.
"""
from __future__ import annotations

import numpy as np

import trackc_sim as tc


def plugin(X, alpha=0.10):
    """CVaR_alpha of the per-scenario mean over repeats."""
    return float(tc.cvar_lower(np.asarray(X, float).mean(0), alpha))


def joint(X, alpha=0.10):
    """CVaR_alpha over every single (repeat, scenario) outcome."""
    return float(tc.cvar_lower(np.ravel(np.asarray(X, float)), alpha))


def _weights(n, alpha):
    """Weights of the sorted order that reproduce tc.cvar_lower: 1 for the lowest floor(alpha n),
    a fraction for the next one, 0 after. Sum = alpha n."""
    m = alpha * n
    k = int(np.floor(m))
    w = np.zeros(n)
    w[:k] = 1.0
    if k < n:
        w[k] = m - k
    return w


def split_sample(X, alpha=0.10, rng=None, n_splits=200):
    """Select the worst scenarios on one random half of the repeats and average their outcomes on the
    other half. Mean over n_splits random splits. Returns NaN if R < 2."""
    X = np.asarray(X, float)
    R, n = X.shape
    if R < 2:
        return float("nan")
    rng = np.random.default_rng(0) if rng is None else rng
    w = _weights(n, alpha)
    out = np.empty(n_splits)
    for s in range(n_splits):
        perm = rng.permutation(R)
        a, b = perm[: R // 2], perm[R // 2:]
        order = np.argsort(X[a].mean(0), kind="stable")
        out[s] = float((w * X[b].mean(0)[order]).sum() / (alpha * n))
    return float(out.mean())


def paired_regret(X_policy, X_base, alpha=0.10):
    """Paired differences (policy - baseline), same scenario and same noise key.
    scenario_cvar : CVaR_alpha over per-scenario mean differences (how bad is the worst decile of scenarios)
    joint_cvar    : CVaR_alpha over every single paired difference
    share_worse   : fraction of scenarios where the policy mean is below the baseline mean
    mean_diff     : mean difference"""
    A, B = np.asarray(X_policy, float), np.asarray(X_base, float)
    if A.shape != B.shape:
        raise ValueError(f"shape mismatch {A.shape} vs {B.shape}")
    d = A - B
    per_scen = d.mean(0)
    return {"scenario_cvar": float(tc.cvar_lower(per_scen, alpha)),
            "joint_cvar": float(tc.cvar_lower(d.ravel(), alpha)),
            "share_worse": float((per_scen < 0).mean()), "mean_diff": float(per_scen.mean())}


def curve(X, Rs, estimator, rng=None, n_subsets=100, **kw):
    """Estimator value as a function of the number of repeats. For each R in Rs (R <= R_max) draw
    n_subsets random subsets of R repeats (one subset if R = R_max) and return {R: (mean, sd)}."""
    X = np.asarray(X, float)
    r_max = X.shape[0]
    rng = np.random.default_rng(0) if rng is None else rng
    out = {}
    for R in Rs:
        if R > r_max:
            raise ValueError(f"R={R} exceeds the {r_max} repeats available")
        reps = 1 if R == r_max else n_subsets
        vals = np.array([estimator(X[rng.choice(r_max, R, replace=False)], **kw) for _ in range(reps)])
        vals = vals[np.isfinite(vals)]
        out[int(R)] = (float(vals.mean()), float(vals.std(ddof=1)) if vals.size > 1 else 0.0) if vals.size \
            else (float("nan"), float("nan"))
    return out


def fit_inverse_R(Rs, values):
    """Least squares fit  value = c - a / R.  Returns (c, a, rms_residual)."""
    Rs = np.asarray(Rs, float)
    y = np.asarray(values, float)
    A = np.column_stack([np.ones_like(Rs), -1.0 / Rs])
    (c, a), *_ = np.linalg.lstsq(A, y, rcond=None)
    resid = y - A @ np.array([c, a])
    return float(c), float(a), float(np.sqrt((resid ** 2).mean()))


def scenario_bootstrap(X_a, X_b, stat, n_boot=1000, seed=0):
    """Point estimate and 95% interval of stat(X_a) - stat(X_b) when scenarios are resampled
    (columns kept together across repeats and across the two arms)."""
    A, B = np.asarray(X_a, float), np.asarray(X_b, float)
    n = A.shape[1]
    rng = np.random.default_rng(seed)
    d = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        d[i] = stat(A[:, idx]) - stat(B[:, idx])
    return float(stat(A) - stat(B)), np.percentile(d, [2.5, 97.5]).tolist()
