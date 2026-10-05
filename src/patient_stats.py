"""Patient-level statistics (analysis_plan_v1.md A1.1, A1.6, A1.16; masterplan §23, §73-77).

The independent unit is the patient. A patient with several forecasts contributes one value
(`patient_means`). Every CI resamples patients, never scans or voxels.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

DEFAULT_N_BOOT = 10_000
DEFAULT_SEED = 20261004


def patient_means(df: pd.DataFrame, value: str, patient: str = "patient_id") -> pd.Series:
    """One value per patient: the equal-weight mean of that patient's forecast-level values (nan rows dropped)."""
    return df.dropna(subset=[value]).groupby(patient)[value].mean()


def bootstrap_ci(x, stat=np.mean, n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED, alpha: float = 0.05):
    """Percentile bootstrap CI of stat(x), resampling the entries of x (each entry = one patient)."""
    x = np.asarray(x, float)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    vals = np.array([stat(x[i]) for i in idx]) if stat not in (np.mean,) else x[idx].mean(1)
    return float(np.percentile(vals, 100 * alpha / 2)), float(np.percentile(vals, 100 * (1 - alpha / 2)))


def cluster_bootstrap(df: pd.DataFrame, fn, patient: str = "patient_id", n_boot: int = DEFAULT_N_BOOT,
                      seed: int = DEFAULT_SEED, alpha: float = 0.05):
    """CI of fn(dataframe) when patients are resampled with ALL of their rows (for forecast-level analyses)."""
    groups = {p: g for p, g in df.groupby(patient)}
    ids = np.array(list(groups))
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(ids, size=len(ids), replace=True)
        vals.append(fn(pd.concat([groups[p] for p in pick], ignore_index=True)))
    vals = np.asarray(vals, float)
    return float(np.nanpercentile(vals, 100 * alpha / 2)), float(np.nanpercentile(vals, 100 * (1 - alpha / 2)))


def paired_sign_flip_p(delta, n_perm: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED) -> float:
    """Two-sided paired permutation test of mean(delta) = 0 by random sign flips (+1 smoothing)."""
    d = np.asarray(delta, float)
    d = d[d != 0] if (d != 0).any() else d
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(n_perm, len(d)))
    null = np.abs((signs * d).mean(1))
    obs = abs(d.mean())
    return float((1 + (null >= obs - 1e-15).sum()) / (n_perm + 1))


def summarize_delta(delta, n_boot: int = DEFAULT_N_BOOT, seed: int = DEFAULT_SEED) -> dict:
    """Primary-endpoint summary of per-patient paired differences delta_i = metric_model - metric_baseline.

    'broadly_positive' (A1.6) = CI lower bound of the mean > 0 AND median > 0.
    """
    d = np.asarray(delta, float)
    d = d[~np.isnan(d)]
    n = len(d)
    lo, hi = bootstrap_ci(d, np.mean, n_boot, seed)
    mlo, mhi = bootstrap_ci(d, np.median, n_boot, seed)
    k = max(1, int(round(0.10 * n)))
    srt = np.sort(d)
    try:
        wil = float(stats.wilcoxon(d[d != 0]).pvalue) if (d != 0).any() else 1.0
    except ValueError:
        wil = float("nan")
    sd = d.std(ddof=1) if n > 1 else float("nan")
    return {"n": n, "mean": float(d.mean()), "mean_ci95": [lo, hi], "median": float(np.median(d)),
            "median_ci95": [mlo, mhi], "iqr": [float(np.percentile(d, 25)), float(np.percentile(d, 75))],
            "pct_improved": float((d > 0).mean() * 100), "pct_worse": float((d < 0).mean() * 100),
            "pct_tied": float((d == 0).mean() * 100),
            "worst10pct_mean": float(srt[:k].mean()), "best10pct_mean": float(srt[-k:].mean()),
            "effect_size_dz": float(d.mean() / sd) if sd and sd > 0 else float("nan"),
            "permutation_p": paired_sign_flip_p(d, n_boot, seed), "wilcoxon_p": wil,
            "broadly_positive": bool(lo > 0 and np.median(d) > 0),
            "mixed": bool((lo > 0) != (np.median(d) > 0))}


def holm(pvals: dict) -> dict:
    """Holm step-down adjusted p-values for a dict name -> p."""
    names = sorted(pvals, key=lambda k: pvals[k])
    m = len(names)
    adj, running = {}, 0.0
    for i, k in enumerate(names):
        running = max(running, min(1.0, (m - i) * pvals[k]))
        adj[k] = running
    return adj


def benjamini_hochberg(pvals: dict) -> dict:
    """BH adjusted p-values (q-values) for a dict name -> p."""
    names = sorted(pvals, key=lambda k: pvals[k])
    m = len(names)
    q, prev = {}, 1.0
    for i in range(m - 1, -1, -1):
        k = names[i]
        prev = min(prev, pvals[k] * m / (i + 1))
        q[k] = prev
    return q


def precision_and_power(delta, target_half_width: float | None = None, effects=None, n_grid=(25, 50, 100, 200, 400),
                        n_sim: int = 2000, seed: int = DEFAULT_SEED) -> dict:
    """Precision / power from the observed spread of per-patient differences (masterplan §77).

    * half-width of a 95% CI for the mean at n patients: 1.96 sd / sqrt(n)
    * n needed for target_half_width
    * power to detect each effect (default: observed mean) with a sign-flip test at alpha 0.05, drawing patients
      with replacement from the observed delta shifted to the effect.
    """
    d = np.asarray(delta, float)
    d = d[~np.isnan(d)]
    sd = float(d.std(ddof=1))
    out = {"sd": sd, "half_width_at_n": {int(n): float(1.96 * sd / np.sqrt(n)) for n in n_grid}}
    if target_half_width:
        out["n_for_target_half_width"] = int(np.ceil((1.96 * sd / target_half_width) ** 2))
    rng = np.random.default_rng(seed)
    effects = [float(d.mean())] if effects is None else list(effects)
    out["power"] = {}
    for eff in effects:
        base = d - d.mean() + eff
        row = {}
        for n in n_grid:
            hit = 0
            for _ in range(n_sim):
                x = rng.choice(base, size=n, replace=True)
                se = x.std(ddof=1) / np.sqrt(n)
                hit += abs(x.mean() / se) > 1.96 if se > 0 else 0
            row[int(n)] = hit / n_sim
        out["power"][f"{eff:.4f}"] = row
    return out
