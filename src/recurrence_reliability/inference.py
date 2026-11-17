"""Patient-level change, cell-bootstrap indeterminate zone, discordance endpoint, intervals."""
from __future__ import annotations

import numpy as np
from scipy.stats import beta

Z = 1.959963984540054   # two-sided 95% normal quantile, fixed


def boot_se(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> float:
    """SD of the nucleus-bootstrap distribution of the specimen mean."""
    v = np.asarray(values, dtype=float)
    idx = rng.integers(0, len(v), size=(n_boot, len(v)))
    return float(v[idx].mean(axis=1).std(ddof=1))


def direction(delta: float, se: float, z: float = Z) -> int:
    """+1 increase, -1 decrease, 0 indeterminate (|delta| <= z * se). Exact ties fall in 0."""
    if not np.isfinite(delta) or not np.isfinite(se):
        raise ValueError("non-finite delta or se")
    if delta > z * se:
        return 1
    if delta < -z * se:
        return -1
    return 0


def clopper_pearson(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    a = (1 - level) / 2
    lo = 0.0 if k == 0 else float(beta.ppf(a, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - a, k + 1, n - k))
    return lo, hi


def one_sided_upper_zero(n: int, level: float = 0.95) -> float:
    return 1 - (1 - level) ** (1 / n)


def agreement_table(dir_a: list[int], dir_b: list[int]) -> dict:
    if len(dir_a) != len(dir_b):
        raise ValueError("direction lists differ in length")
    n = len(dir_a)
    cells = {f"{a:+d}|{b:+d}".replace("+0", "0"): 0 for a in (1, 0, -1) for b in (1, 0, -1)}
    for a, b in zip(dir_a, dir_b):
        cells[f"{a:+d}|{b:+d}".replace("+0", "0")] += 1
    opp = sum(1 for a, b in zip(dir_a, dir_b) if a * b == -1)
    one_ind = sum(1 for a, b in zip(dir_a, dir_b) if (a == 0) != (b == 0))
    both_ind = sum(1 for a, b in zip(dir_a, dir_b) if a == 0 and b == 0)
    same = sum(1 for a, b in zip(dir_a, dir_b) if a == b and a != 0)
    lo, hi = clopper_pearson(opp, n)
    return {"n": n, "table_R_by_A": cells, "opposite_direction": opp, "opposite_fraction": opp / n,
            "opposite_ci95_clopper_pearson": [lo, hi], "one_indeterminate": one_ind, "both_indeterminate": both_ind,
            "same_determinate_direction": same,
            "any_disagreement": opp + one_ind, "any_disagreement_ci95": list(clopper_pearson(opp + one_ind, n))}


def patient_bootstrap_mean(x: np.ndarray, rng: np.random.Generator, n_boot: int) -> dict:
    x = np.asarray(x, dtype=float)
    idx = rng.integers(0, len(x), size=(n_boot, len(x)))
    m = x[idx].mean(axis=1)
    return {"mean": float(x.mean()), "median": float(np.median(x)),
            "ci95_percentile": [float(np.quantile(m, 0.025)), float(np.quantile(m, 0.975))], "n": len(x)}
