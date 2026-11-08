"""Outcome-blind coverage screen for a proposed two-family t-Bonferroni method."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from collections.abc import Mapping

import numpy as np
from scipy.stats import t


FAMILIES = ("single-response", "mixture-response")
DISTRIBUTION_PAIRS = (
    ("normal_normal", "normal", "normal"),
    ("gamma_gamma", "centered_gamma2", "centered_gamma2"),
    ("t3_t3", "student_t3", "student_t3"),
    ("uniform_uniform", "uniform", "uniform"),
    ("normal_gamma", "normal", "centered_gamma2"),
    ("gamma_t3", "centered_gamma2", "student_t3"),
)
WORLD_COUNTS = (5, 10, 14, 20)
MONTE_CARLO_REPETITIONS = 20_000
MASTER_SEED = 2026110809
MIN_WILSON_LOWER = 0.93
CONFIDENCE_LEVEL = 0.95
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = PROJECT_ROOT / "paper" / "trackc_t_bonferroni_coverage_protocol_2026-11-08.md"
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackc_t_bonferroni_coverage_audit_20261108" / "coverage.json"


def draw_standardized(name: str, size: int, rng: np.random.Generator) -> np.ndarray:
    """Draw from one fixed mean-zero, variance-one audit distribution."""
    if isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size < 1:
        raise ValueError("size must be a positive integer")
    if not isinstance(rng, np.random.Generator):
        raise TypeError("rng must be a numpy Generator")
    if name == "normal":
        values = rng.normal(0.0, 1.0, size=size)
    elif name == "centered_gamma2":
        values = (rng.gamma(shape=2.0, scale=1.0, size=size) - 2.0) / math.sqrt(2.0)
    elif name == "student_t3":
        values = rng.standard_t(df=3, size=size) / math.sqrt(3.0)
    elif name == "uniform":
        values = rng.uniform(-math.sqrt(3.0), math.sqrt(3.0), size=size)
    else:
        raise ValueError(f"unknown audit distribution: {name}")
    return np.asarray(values, dtype=float)


def bonferroni_t_intervals(
    samples_by_family: Mapping[str, object], *, confidence_level: float = CONFIDENCE_LEVEL
) -> dict:
    """Return two-sided simultaneous Student-t intervals for two family means.

    Exact marginal coverage requires iid normally distributed world differences
    within each family. Bonferroni does not require dependence between families.
    """
    if not isinstance(samples_by_family, Mapping) or set(samples_by_family) != set(FAMILIES):
        raise ValueError("samples must contain exactly the two required response families")
    if isinstance(confidence_level, (bool, np.bool_)):
        raise ValueError("confidence_level must lie strictly between zero and one")
    try:
        confidence_level = float(confidence_level)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError("confidence_level must be numeric") from exc
    if not np.isfinite(confidence_level) or not 0.0 < confidence_level < 1.0:
        raise ValueError("confidence_level must lie strictly between zero and one")

    arrays = {}
    for family in FAMILIES:
        values = np.asarray(samples_by_family[family], dtype=float)
        if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
            raise ValueError(f"{family} needs at least two finite world differences")
        if float(values.std(ddof=1)) <= 0.0:
            raise ValueError(f"{family} needs positive observed variance")
        arrays[family] = values
    sizes = {values.size for values in arrays.values()}
    if len(sizes) != 1:
        raise ValueError("both response families must have the same world count")

    n_worlds = sizes.pop()
    family_count = len(FAMILIES)
    alpha = 1.0 - confidence_level
    critical = float(t.ppf(1.0 - alpha / (2.0 * family_count), df=n_worlds - 1))
    intervals = {}
    for family, values in arrays.items():
        mean = float(values.mean())
        standard_error = float(values.std(ddof=1) / math.sqrt(n_worlds))
        half_width = critical * standard_error
        intervals[family] = {
            "n_worlds": int(n_worlds),
            "mean": mean,
            "standard_error": standard_error,
            "critical_value": critical,
            "lower": mean - half_width,
            "upper": mean + half_width,
            "half_width": half_width,
        }
    return {
        "method": "two-sided Student-t intervals with Bonferroni adjustment",
        "confidence_level": confidence_level,
        "familywise_alpha": alpha,
        "family_count": family_count,
        "critical_value": critical,
        "intervals": intervals,
    }


def wilson_interval(successes: int, total: int, *, z: float = 1.959963984540054) -> tuple[float, float]:
    """Return a two-sided Wilson score interval for a binomial proportion."""
    if isinstance(successes, bool) or not isinstance(successes, (int, np.integer)):
        raise ValueError("successes must be an integer")
    if isinstance(total, bool) or not isinstance(total, (int, np.integer)) or total < 1:
        raise ValueError("total must be a positive integer")
    if not 0 <= successes <= total:
        raise ValueError("successes must lie between zero and total")
    if not np.isfinite(z) or z <= 0.0:
        raise ValueError("z must be finite and positive")
    if successes == 0:
        z2 = z * z
        return 0.0, z2 / (total + z2)
    if successes == total:
        z2 = z * z
        return total / (total + z2), 1.0
    p = successes / total
    z2 = z * z
    denominator = 1.0 + z2 / total
    center = (p + z2 / (2.0 * total)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / total + z2 / (4.0 * total * total)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def _draw_matrix(name: str, repetitions: int, n_worlds: int, rng: np.random.Generator) -> np.ndarray:
    """Vectorized version of draw_standardized for the fixed audit grid."""
    if name == "normal":
        return rng.normal(0.0, 1.0, size=(repetitions, n_worlds))
    if name == "centered_gamma2":
        return (rng.gamma(shape=2.0, scale=1.0, size=(repetitions, n_worlds)) - 2.0) / math.sqrt(2.0)
    if name == "student_t3":
        return rng.standard_t(df=3, size=(repetitions, n_worlds)) / math.sqrt(3.0)
    if name == "uniform":
        return rng.uniform(-math.sqrt(3.0), math.sqrt(3.0), size=(repetitions, n_worlds))
    raise ValueError(f"unknown audit distribution: {name}")


def run_coverage_screen() -> dict:
    """Run fixed, outcome-blind coverage screen without project data or policies."""
    cells = []
    failures_total = 0
    alpha = 1.0 - CONFIDENCE_LEVEL
    critical_cache = {
        n_worlds: float(t.ppf(1.0 - alpha / (2.0 * len(FAMILIES)), df=n_worlds - 1))
        for n_worlds in WORLD_COUNTS
    }
    for pair_index, (pair_name, first_name, second_name) in enumerate(DISTRIBUTION_PAIRS):
        for n_worlds in WORLD_COUNTS:
            draws = {}
            for family_index, (family, distribution) in enumerate(zip(FAMILIES, (first_name, second_name))):
                rng = np.random.Generator(
                    np.random.PCG64DXSM(
                        np.random.SeedSequence([MASTER_SEED, pair_index, n_worlds, family_index])
                    )
                )
                draws[family] = _draw_matrix(distribution, MONTE_CARLO_REPETITIONS, n_worlds, rng)

            critical = critical_cache[n_worlds]
            covers = []
            positives = []
            valid_all = np.ones(MONTE_CARLO_REPETITIONS, dtype=bool)
            for family in FAMILIES:
                values = draws[family]
                means = values.mean(axis=1)
                standard_deviations = values.std(axis=1, ddof=1)
                valid = np.isfinite(means) & np.isfinite(standard_deviations) & (standard_deviations > 0.0)
                valid_all &= valid
                half_widths = critical * standard_deviations / math.sqrt(n_worlds)
                lower = means - half_widths
                upper = means + half_widths
                covers.append(valid & (lower <= 0.0) & (upper >= 0.0))
                positives.append(valid & (lower > 0.0))

            joint_cover = covers[0] & covers[1]
            joint_positive = positives[0] & positives[1]
            failures = int((~valid_all).sum())
            failures_total += failures
            coverage_count = int(joint_cover.sum())
            positive_count = int(joint_positive.sum())
            coverage_lower, coverage_upper = wilson_interval(coverage_count, MONTE_CARLO_REPETITIONS)
            positive_lower, positive_upper = wilson_interval(positive_count, MONTE_CARLO_REPETITIONS)
            cells.append({
                "pair": pair_name,
                "family_distributions": {
                    FAMILIES[0]: first_name,
                    FAMILIES[1]: second_name,
                },
                "world_count_per_family": n_worlds,
                "repetitions": MONTE_CARLO_REPETITIONS,
                "joint_coverage_count": coverage_count,
                "joint_coverage": coverage_count / MONTE_CARLO_REPETITIONS,
                "joint_coverage_wilson_95": [coverage_lower, coverage_upper],
                "both_lower_bounds_above_zero_count": positive_count,
                "both_lower_bounds_above_zero_rate": positive_count / MONTE_CARLO_REPETITIONS,
                "false_positive_wilson_95": [positive_lower, positive_upper],
                "inference_failures": failures,
                "screen_cell_pass": coverage_lower >= MIN_WILSON_LOWER and failures == 0,
            })

    all_pass = all(cell["screen_cell_pass"] for cell in cells) and failures_total == 0
    return {
        "record_type": "outcome_blind_trackc_t_bonferroni_coverage_screen",
        "status": "candidate_screen_passed" if all_pass else "candidate_screen_failed",
        "method": "two-sided Student-t intervals with Bonferroni adjustment across two response families",
        "estimand": "joint coverage of two zero-mean response-family means from independent world-level differences",
        "confidence_level": CONFIDENCE_LEVEL,
        "world_counts_per_family": list(WORLD_COUNTS),
        "monte_carlo_repetitions_per_cell": MONTE_CARLO_REPETITIONS,
        "master_seed": MASTER_SEED,
        "rng": "PCG64DXSM with SeedSequence([master_seed, pair_index, world_count, family_index])",
        "distribution_pairs": [
            {"name": name, FAMILIES[0]: first, FAMILIES[1]: second}
            for name, first, second in DISTRIBUTION_PAIRS
        ],
        "minimum_wilson_lower_bound": MIN_WILSON_LOWER,
        "inference_failures_total": failures_total,
        "cells_passed": sum(bool(cell["screen_cell_pass"]) for cell in cells),
        "cells_total": len(cells),
        "policy_outcomes_used": False,
        "patient_data_used": False,
        "protocol_status": "proposed_screen_only_not_frozen",
        "normality_assumption": "Exact marginal Student-t coverage requires independent normally distributed world differences within each family; non-normal coverage is empirical only.",
        "limitations": "Passing this generic distribution screen would not validate the actual Track C world-difference distribution, freeze inference, establish benefit, or authorize policy evaluation.",
        "cells": cells,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    if OUTPUT_PATH.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT_PATH}")
    if not PROTOCOL_PATH.is_file():
        raise FileNotFoundError(f"missing protocol: {PROTOCOL_PATH}")
    result = run_coverage_screen()
    result["provenance"] = {
        "code_sha256": _sha256(Path(__file__).resolve()),
        "protocol_sha256": _sha256(PROTOCOL_PATH),
        "python": __import__("sys").version,
        "numpy": np.__version__,
    }
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"status={result['status']} cells={result['cells_passed']}/{result['cells_total']} output={OUTPUT_PATH}")


if __name__ == "__main__":
    main()
