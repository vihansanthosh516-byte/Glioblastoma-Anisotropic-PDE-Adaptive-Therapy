"""Outcome-free coverage screen for fixed-horizon bounded-mean intervals.

Protocol: paper/trackc_ge_betting_coverage_audit_protocol_2026-11-09.md.
Synthetic bounded draws only. No Track C policy outcomes, patient data, or
evaluation seeds are used. Passing this screen adopts no method.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from multiprocessing import Pool
from pathlib import Path
from statistics import NormalDist
from time import perf_counter

import numpy as np
from scipy import stats

from trackc_ge_betting import (
    _accepts_mean,
    _variance_path,
    ge_betting_interval_unit,
)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output" / "trackc_ge_betting_coverage_audit_20261109"

AUDIT_NAMESPACE = "trackc-ge-betting-coverage-audit-v1"
MASTER_SEED = 20261109
SUPPORT = (-790.0, 790.0)
SUPPORT_WIDTH = SUPPORT[1] - SUPPORT[0]
DELTA_PER_FAMILY = 0.025
WORLD_COUNTS = (7, 14, 30, 100, 207)
MONTE_CARLO_REPETITIONS = 2000
WIDTH_REPETITIONS = 200
JOINT_COVERAGE_LOWER_BOUND = 0.93
PRECISION_TARGET_DAYS = 14.0
CONSISTENCY_TOLERANCE_UNIT = 1e-8
DISTRIBUTION_PAIRS = (
    ("rare_pos_rare_pos", "rare_pos", "rare_pos"),
    ("rare_pos_rare_neg", "rare_pos", "rare_neg"),
    ("two_point_two_point", "two_point", "two_point"),
    ("uniform_uniform", "uniform", "uniform"),
    ("narrow_normal_narrow_normal", "narrow_normal", "narrow_normal"),
    ("skew_beta_narrow_normal", "skew_beta", "narrow_normal"),
)
METHODS = ("ge_betting", "hoeffding", "empirical_bernstein", "student_t")
NARROW_MEAN, NARROW_SD = 10.0, 30.0
BETA_A, BETA_B = 0.5, 5.0


def rare_event_probability(n: int) -> float:
    return 1.0 - DELTA_PER_FAMILY ** (1.0 / n)


def draw_days(distribution: str, rng: np.random.Generator, shape: tuple[int, int]) -> tuple[np.ndarray, float]:
    """Return draws on the day scale and the exact distribution mean."""
    n = shape[1]
    if distribution == "rare_pos":
        p = rare_event_probability(n)
        return np.where(rng.random(shape) < p, SUPPORT[1], 0.0), SUPPORT[1] * p
    if distribution == "rare_neg":
        p = rare_event_probability(n)
        return np.where(rng.random(shape) < p, SUPPORT[0], 0.0), SUPPORT[0] * p
    if distribution == "two_point":
        return np.where(rng.random(shape) < 0.5, SUPPORT[1], SUPPORT[0]), 0.0
    if distribution == "uniform":
        return rng.uniform(SUPPORT[0], SUPPORT[1], size=shape), 0.0
    if distribution == "narrow_normal":
        draws = rng.normal(NARROW_MEAN, NARROW_SD, size=shape)
        if (draws < SUPPORT[0]).any() or (draws > SUPPORT[1]).any():
            raise ArithmeticError("narrow_normal draw left the support")
        return draws, NARROW_MEAN
    if distribution == "skew_beta":
        draws = SUPPORT_WIDTH * rng.beta(BETA_A, BETA_B, size=shape) + SUPPORT[0]
        return draws, SUPPORT_WIDTH * BETA_A / (BETA_A + BETA_B) + SUPPORT[0]
    raise ValueError(f"unknown distribution: {distribution}")


def wilson_interval(successes: int, trials: int, confidence_level: float = 0.95) -> tuple[float, float]:
    if trials < 1 or successes < 0 or successes > trials:
        raise ValueError("invalid binomial counts")
    z = NormalDist().inv_cdf(0.5 + confidence_level / 2.0)
    proportion = successes / trials
    denominator = 1.0 + z * z / trials
    center = (proportion + z * z / (2.0 * trials)) / denominator
    half = z * math.sqrt(proportion * (1.0 - proportion) / trials + z * z / (4.0 * trials * trials)) / denominator
    return center - half, center + half


def hoeffding_half_width_days(n: int, delta: float = DELTA_PER_FAMILY) -> float:
    return SUPPORT_WIDTH * math.sqrt(math.log(2.0 / delta) / (2.0 * n))


def empirical_bernstein_half_width_days(days: np.ndarray, delta: float = DELTA_PER_FAMILY) -> np.ndarray:
    """Maurer-Pontil two-sided half-width per row, on the day scale."""
    n = days.shape[-1]
    unit_variance = np.var(days / SUPPORT_WIDTH, axis=-1, ddof=1)
    log_term = math.log(4.0 / delta)
    return SUPPORT_WIDTH * (np.sqrt(2.0 * unit_variance * log_term / n) + 7.0 * log_term / (3.0 * (n - 1)))


def student_t_half_width_days(days: np.ndarray, delta: float = DELTA_PER_FAMILY) -> np.ndarray:
    n = days.shape[-1]
    critical = float(stats.t.ppf(1.0 - delta / 2.0, df=n - 1))
    return critical * np.std(days, axis=-1, ddof=1) / math.sqrt(n)


def _clipped_bounds(center: np.ndarray, half: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.maximum(center - half, SUPPORT[0]), np.minimum(center + half, SUPPORT[1])


def _family_arrays(days: np.ndarray, true_mean: float) -> dict:
    """Per-repetition coverage flags and half-widths for one family."""
    repetitions, n = days.shape
    center = days.mean(axis=1)
    out: dict = {}

    hoeffding = np.full(repetitions, hoeffding_half_width_days(n))
    for name, half in (
        ("hoeffding", hoeffding),
        ("empirical_bernstein", empirical_bernstein_half_width_days(days)),
        ("student_t", student_t_half_width_days(days)),
    ):
        lower, upper = _clipped_bounds(center, half)
        out[name] = {
            "covered": (lower <= true_mean) & (true_mean <= upper),
            "half_width": (upper - lower) / 2.0,
            "failed": np.zeros(repetitions, dtype=bool),
        }

    unit = (days - SUPPORT[0]) / SUPPORT_WIDTH
    unit_mean = (true_mean - SUPPORT[0]) / SUPPORT_WIDTH
    covered = np.zeros(repetitions, dtype=bool)
    failed = np.zeros(repetitions, dtype=bool)
    half_width = np.full(WIDTH_REPETITIONS, np.nan)
    inconsistent = 0
    for index in range(repetitions):
        row = unit[index]
        covered[index] = _accepts_mean(row, unit_mean, DELTA_PER_FAMILY, _variance_path(row))
        if index < WIDTH_REPETITIONS:
            interval = ge_betting_interval_unit(row, delta=DELTA_PER_FAMILY)
            if interval is None:
                failed[index] = True
                continue
            half_width[index] = SUPPORT_WIDTH * (interval[1] - interval[0]) / 2.0
            inside = (
                interval[0] - CONSISTENCY_TOLERANCE_UNIT <= unit_mean <= interval[1] + CONSISTENCY_TOLERANCE_UNIT
            )
            strictly_inside = (
                interval[0] + CONSISTENCY_TOLERANCE_UNIT <= unit_mean <= interval[1] - CONSISTENCY_TOLERANCE_UNIT
            )
            if (covered[index] and not inside) or (strictly_inside and not covered[index]):
                inconsistent += 1
    out["ge_betting"] = {
        "covered": covered,
        "half_width": half_width,
        "failed": failed,
        "inconsistent": inconsistent,
    }
    return out


def run_cell(task: tuple[int, int, str, str, str, int]) -> dict:
    cell_index, n, pair_name, first, second, repetitions = task
    families = []
    means = []
    for family_index, distribution in enumerate((first, second)):
        rng = np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence([MASTER_SEED, cell_index, family_index])))
        days, true_mean = draw_days(distribution, rng, (repetitions, n))
        means.append(true_mean)
        families.append(_family_arrays(days, true_mean))

    methods = {}
    for method in METHODS:
        a, b = families[0][method], families[1][method]
        joint = int(np.count_nonzero(a["covered"] & b["covered"]))
        failures = int(np.count_nonzero(a["failed"]) + np.count_nonzero(b["failed"]))
        lower, upper = wilson_interval(joint, repetitions)
        width_count = min(len(a["half_width"]), len(b["half_width"]), repetitions)
        larger = np.maximum(a["half_width"][:width_count], b["half_width"][:width_count])
        larger = larger[np.isfinite(larger)]
        inconsistent = int(a.get("inconsistent", 0) + b.get("inconsistent", 0))
        methods[method] = {
            "joint_coverage": joint / repetitions,
            "wilson_lower": lower,
            "wilson_upper": upper,
            "marginal_coverage": [float(a["covered"].mean()), float(b["covered"].mean())],
            "inference_failures": failures,
            "inversion_inconsistencies": inconsistent,
            "half_width_repetitions": int(larger.size),
            "median_larger_half_width_days": float(np.median(larger)) if larger.size else None,
            "p90_larger_half_width_days": float(np.quantile(larger, 0.9)) if larger.size else None,
            "fraction_both_half_widths_within_target": float(np.mean(larger <= PRECISION_TARGET_DAYS)) if larger.size else None,
            "passed": bool(lower >= JOINT_COVERAGE_LOWER_BOUND and failures == 0 and inconsistent == 0),
        }
    return {
        "cell_index": cell_index,
        "worlds_per_family": n,
        "distribution_pair": pair_name,
        "true_means_days": means,
        "methods": methods,
    }


def build_tasks(repetitions: int = MONTE_CARLO_REPETITIONS, world_counts=WORLD_COUNTS) -> list:
    tasks = []
    for n in world_counts:
        for pair_name, first, second in DISTRIBUTION_PAIRS:
            tasks.append((len(tasks), n, pair_name, first, second, repetitions))
    return tasks


def summarize(cells: list) -> dict:
    summary = {}
    for method in METHODS:
        passed = [cell["methods"][method]["passed"] for cell in cells]
        coverage = [cell["methods"][method]["joint_coverage"] for cell in cells]
        summary[method] = {
            "cells": len(cells),
            "cells_passed": int(sum(passed)),
            "screen_passed": bool(all(passed)),
            "minimum_joint_coverage": min(coverage),
            "maximum_joint_coverage": max(coverage),
        }
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()

    started = perf_counter()
    tasks = build_tasks()
    with Pool(processes=args.processes) as pool:
        cells = pool.map(run_cell, tasks, chunksize=1)
    cells.sort(key=lambda cell: cell["cell_index"])

    sources = {
        name: hashlib.sha256((Path(__file__).resolve().parent / name).read_bytes()).hexdigest()
        for name in ("trackc_ge_betting.py", "trackc_ge_betting_coverage_audit.py")
    }
    result = {
        "audit_namespace": AUDIT_NAMESPACE,
        "protocol": "paper/trackc_ge_betting_coverage_audit_protocol_2026-11-09.md",
        "policy_outcomes_used": False,
        "patient_data_used": False,
        "evaluation_seeds_used": False,
        "settings": {
            "master_seed": MASTER_SEED,
            "support_days": list(SUPPORT),
            "delta_per_family": DELTA_PER_FAMILY,
            "world_counts_per_family": list(WORLD_COUNTS),
            "monte_carlo_repetitions_per_cell": MONTE_CARLO_REPETITIONS,
            "ge_betting_width_repetitions_per_cell": WIDTH_REPETITIONS,
            "joint_coverage_wilson_lower_bound_minimum": JOINT_COVERAGE_LOWER_BOUND,
            "precision_target_days": PRECISION_TARGET_DAYS,
            "bit_generator": "PCG64DXSM",
        },
        "source_sha256": sources,
        "software": {"python": platform.python_version(), "numpy": np.__version__},
        "runtime_seconds": None,
        "summary": summarize(cells),
        "cells": cells,
    }
    result["runtime_seconds"] = perf_counter() - started
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUTPUT_DIR / "coverage.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
