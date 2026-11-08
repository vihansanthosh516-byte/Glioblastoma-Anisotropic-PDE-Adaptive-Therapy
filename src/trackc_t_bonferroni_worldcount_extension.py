"""Outcome-blind larger-world-count screen for the fixed t-Bonferroni candidate."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import scipy
from scipy.stats import t

from trackc_t_bonferroni_coverage_audit import (
    CONFIDENCE_LEVEL,
    DISTRIBUTION_PAIRS as BASE_DISTRIBUTION_PAIRS,
    FAMILIES,
    MIN_WILSON_LOWER,
    _draw_matrix as _draw_base_matrix,
    wilson_interval,
)


DISTRIBUTION_PAIRS = BASE_DISTRIBUTION_PAIRS + (
    ("beta_2_5_pair", "beta_2_5", "beta_2_5"),
    ("beta_mirror_pair", "beta_2_5", "beta_5_2"),
    ("bernoulli_0_1_pair", "bernoulli_0_1", "bernoulli_0_1"),
    ("gamma_bernoulli_pair", "centered_gamma2", "bernoulli_0_1"),
    ("t3_bernoulli_pair", "student_t3", "bernoulli_0_1"),
)
WORLD_COUNTS = (30, 50, 100, 200)
MONTE_CARLO_REPETITIONS = 20_000
MASTER_SEED = 2026110811
PROJECT_ROOT = Path(__file__).resolve().parents[1]
PROTOCOL_PATH = PROJECT_ROOT / "paper" / "trackc_t_bonferroni_worldcount_extension_protocol_2026-11-08.md"
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackc_t_bonferroni_worldcount_extension_20261108" / "coverage.json"


def draw_standardized(name: str, size: int, rng: np.random.Generator) -> np.ndarray:
    """Draw a fixed mean-zero, variance-one distribution used in the extension."""
    if isinstance(size, bool) or not isinstance(size, (int, np.integer)) or size < 1:
        raise ValueError("size must be a positive integer")
    if name in ("normal", "centered_gamma2", "student_t3", "uniform"):
        return _draw_base_matrix(name, 1, int(size), rng)[0]
    if name in ("beta_2_5", "beta_5_2"):
        a, b = (2.0, 5.0) if name == "beta_2_5" else (5.0, 2.0)
        mean = a / (a + b)
        variance = a * b / ((a + b) ** 2 * (a + b + 1.0))
        return (rng.beta(a, b, size=size) - mean) / math.sqrt(variance)
    if name == "bernoulli_0_1":
        p = 0.1
        return (rng.binomial(1, p, size=size).astype(float) - p) / math.sqrt(p * (1.0 - p))
    raise ValueError(f"unknown extension distribution: {name}")


def _draw_audit_matrix(name: str, repetitions: int, n_worlds: int, rng: np.random.Generator) -> np.ndarray:
    """Vectorized draws for one family and one Monte Carlo cell."""
    if name in ("normal", "centered_gamma2", "student_t3", "uniform"):
        return _draw_base_matrix(name, repetitions, n_worlds, rng)
    if name in ("beta_2_5", "beta_5_2"):
        a, b = (2.0, 5.0) if name == "beta_2_5" else (5.0, 2.0)
        mean = a / (a + b)
        variance = a * b / ((a + b) ** 2 * (a + b + 1.0))
        return (rng.beta(a, b, size=(repetitions, n_worlds)) - mean) / math.sqrt(variance)
    if name == "bernoulli_0_1":
        p = 0.1
        draws = rng.binomial(1, p, size=(repetitions, n_worlds)).astype(float)
        return (draws - p) / math.sqrt(p * (1.0 - p))
    raise ValueError(f"unknown extension distribution: {name}")


def run_coverage_screen() -> dict:
    """Run the fixed large-world-count screen; do not load project outcomes."""
    alpha = 1.0 - CONFIDENCE_LEVEL
    critical_by_n = {
        n_worlds: float(t.ppf(1.0 - alpha / (2.0 * len(FAMILIES)), df=n_worlds - 1))
        for n_worlds in WORLD_COUNTS
    }
    cells = []
    total_failures = 0
    for pair_index, (pair_name, first_distribution, second_distribution) in enumerate(DISTRIBUTION_PAIRS):
        for n_worlds in WORLD_COUNTS:
            family_draws = {}
            for family_index, (family, distribution) in enumerate(
                zip(FAMILIES, (first_distribution, second_distribution))
            ):
                rng = np.random.Generator(
                    np.random.PCG64DXSM(
                        np.random.SeedSequence([MASTER_SEED, pair_index, n_worlds, family_index])
                    )
                )
                family_draws[family] = _draw_audit_matrix(
                    distribution, MONTE_CARLO_REPETITIONS, n_worlds, rng
                )

            critical = critical_by_n[n_worlds]
            cover_by_family = []
            lower_positive_by_family = []
            valid_all = np.ones(MONTE_CARLO_REPETITIONS, dtype=bool)
            for family in FAMILIES:
                values = family_draws[family]
                means = values.mean(axis=1)
                standard_deviations = values.std(axis=1, ddof=1)
                valid = np.isfinite(means) & np.isfinite(standard_deviations) & (standard_deviations > 0.0)
                valid_all &= valid
                half_widths = critical * standard_deviations / math.sqrt(n_worlds)
                lower = means - half_widths
                upper = means + half_widths
                cover_by_family.append(valid & (lower <= 0.0) & (upper >= 0.0))
                lower_positive_by_family.append(valid & (lower > 0.0))

            joint_cover = cover_by_family[0] & cover_by_family[1]
            joint_positive = lower_positive_by_family[0] & lower_positive_by_family[1]
            failures = int((~valid_all).sum())
            total_failures += failures
            cover_count = int(joint_cover.sum())
            positive_count = int(joint_positive.sum())
            cover_wilson = wilson_interval(cover_count, MONTE_CARLO_REPETITIONS)
            positive_wilson = wilson_interval(positive_count, MONTE_CARLO_REPETITIONS)
            cells.append({
                "pair": pair_name,
                "family_distributions": {
                    FAMILIES[0]: first_distribution,
                    FAMILIES[1]: second_distribution,
                },
                "world_count_per_family": n_worlds,
                "repetitions": MONTE_CARLO_REPETITIONS,
                "joint_coverage_count": cover_count,
                "joint_coverage": cover_count / MONTE_CARLO_REPETITIONS,
                "joint_coverage_wilson_95": list(cover_wilson),
                "both_lower_bounds_above_zero_count": positive_count,
                "both_lower_bounds_above_zero_rate": positive_count / MONTE_CARLO_REPETITIONS,
                "false_positive_wilson_95": list(positive_wilson),
                "inference_failures": failures,
                "screen_cell_pass": cover_wilson[0] >= MIN_WILSON_LOWER and failures == 0,
            })

    passed_by_n = {}
    for n_worlds in WORLD_COUNTS:
        group = [cell for cell in cells if cell["world_count_per_family"] == n_worlds]
        passed_by_n[str(n_worlds)] = sum(bool(cell["screen_cell_pass"]) for cell in group)
    cells_passed = sum(bool(cell["screen_cell_pass"]) for cell in cells)
    minimum_coverage = min(cell["joint_coverage"] for cell in cells)
    maximum_coverage = max(cell["joint_coverage"] for cell in cells)
    return {
        "record_type": "outcome_blind_trackc_t_bonferroni_worldcount_extension",
        "status": "candidate_screen_passed" if cells_passed == len(cells) and total_failures == 0 else "candidate_screen_failed",
        "method": "two-sided Student-t intervals with Bonferroni adjustment across two response families",
        "confidence_level": CONFIDENCE_LEVEL,
        "world_counts_per_family": list(WORLD_COUNTS),
        "monte_carlo_repetitions_per_cell": MONTE_CARLO_REPETITIONS,
        "master_seed": MASTER_SEED,
        "rng": "PCG64DXSM with SeedSequence([master_seed, pair_index, world_count, family_index])",
        "critical_value_by_world_count": {str(n): critical_by_n[n] for n in WORLD_COUNTS},
        "distribution_pairs": [
            {"name": name, FAMILIES[0]: first, FAMILIES[1]: second}
            for name, first, second in DISTRIBUTION_PAIRS
        ],
        "cells_passed": cells_passed,
        "cells_total": len(cells),
        "cells_passed_by_world_count": passed_by_n,
        "minimum_joint_coverage": minimum_coverage,
        "maximum_joint_coverage": maximum_coverage,
        "inference_failures_total": total_failures,
        "minimum_wilson_lower_bound": min(cell["joint_coverage_wilson_95"][0] for cell in cells),
        "maximum_false_positive_rate": max(cell["both_lower_bounds_above_zero_rate"] for cell in cells),
        "patient_data_used": False,
        "policy_outcomes_used": False,
        "protocol_status": "proposed_outcome_blind_feasibility_screen_only",
        "normality_assumption": "Exact marginal Student-t coverage requires independent normally distributed world differences within each family; all other coverage is empirical only.",
        "interpretation_limit": "Even a generic screen pass would not validate the actual Track C world-difference distribution or authorize policy evaluation.",
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
        "helper_sha256": _sha256(PROJECT_ROOT / "src" / "trackc_t_bonferroni_coverage_audit.py"),
        "protocol_sha256": _sha256(PROTOCOL_PATH),
        "python": __import__("sys").version,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"status={result['status']} cells={result['cells_passed']}/{result['cells_total']} output={OUTPUT_PATH}")


if __name__ == "__main__":
    main()
