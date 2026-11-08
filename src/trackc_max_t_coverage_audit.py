"""Outcome-blind coverage screen for the proposed Track C max-T intervals."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
from pathlib import Path
from statistics import NormalDist
from time import perf_counter
from typing import Callable

import numpy as np

from trackc_world_summary import simultaneous_max_t_world_intervals


AUDIT_NAMESPACE = "trackc-max-t-coverage-audit-v1"
MASTER_SEED = 20261108
WORLD_COUNTS = (5, 10, 14, 20)
MONTE_CARLO_REPETITIONS = 2000
BOOTSTRAP_RESAMPLES = 2000
JOINT_COVERAGE_LOWER_BOUND = 0.93
CONFIDENCE_LEVEL = 0.95
DISTRIBUTION_PAIRS = (
    ("normal_normal", "normal", "normal"),
    ("gamma2_gamma2", "gamma2", "gamma2"),
    ("student_t3_student_t3", "student_t3", "student_t3"),
    ("uniform_uniform", "uniform", "uniform"),
    ("normal_gamma2", "normal", "gamma2"),
    ("gamma2_student_t3", "gamma2", "student_t3"),
)
FAMILIES = ("single-response", "mixture-response")


def _draw_standardized(distribution: str, rng: np.random.Generator, size: int) -> np.ndarray:
    if distribution == "normal":
        return rng.standard_normal(size=size)
    if distribution == "gamma2":
        return (rng.gamma(shape=2.0, scale=1.0, size=size) - 2.0) / math.sqrt(2.0)
    if distribution == "student_t3":
        return rng.standard_t(df=3, size=size) / math.sqrt(3.0)
    if distribution == "uniform":
        return rng.uniform(-math.sqrt(3.0), math.sqrt(3.0), size=size)
    raise ValueError(f"unknown distribution: {distribution}")


def _wilson_interval(successes: int, trials: int, confidence_level: float) -> tuple[float, float]:
    if trials < 1 or successes < 0 or successes > trials:
        raise ValueError("invalid binomial counts")
    z = NormalDist().inv_cdf(0.5 + confidence_level / 2.0)
    proportion = successes / trials
    denominator = 1.0 + (z * z / trials)
    center = (proportion + z * z / (2.0 * trials)) / denominator
    half_width = z * math.sqrt(
        proportion * (1.0 - proportion) / trials + (z * z / (4.0 * trials * trials))
    ) / denominator
    return center - half_width, center + half_width


def _bootstrap_seed(pair_name: str, n_worlds: int, repetition: int) -> int:
    key = f"{AUDIT_NAMESPACE}|bootstrap|{MASTER_SEED}|{pair_name}|{n_worlds}|{repetition}"
    return int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:8], "little")


def run_coverage_audit(
    *,
    repetitions: int = MONTE_CARLO_REPETITIONS,
    bootstrap_resamples: int = BOOTSTRAP_RESAMPLES,
    world_counts: tuple[int, ...] = WORLD_COUNTS,
    progress_callback: Callable[[dict], None] | None = None,
) -> dict:
    """Run fixed null-coverage scenarios; parameters allow fast mechanics tests."""
    for label, value in (("repetitions", repetitions), ("bootstrap_resamples", bootstrap_resamples)):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
            raise ValueError(f"{label} must be a positive integer")
    counts = tuple(world_counts)
    if not counts or any(
        isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 2
        for value in counts
    ):
        raise ValueError("world_counts must contain integers of at least 2")

    cells = []
    for pair_index, (pair_name, first_distribution, second_distribution) in enumerate(DISTRIBUTION_PAIRS):
        for n_worlds in counts:
            joint_covered = 0
            positive_rule_false_positives = 0
            inference_failures = 0
            for repetition in range(repetitions):
                data_rng = np.random.Generator(
                    np.random.PCG64DXSM(
                        np.random.SeedSequence([MASTER_SEED, pair_index, int(n_worlds), repetition])
                    )
                )
                differences = {
                    FAMILIES[0]: _draw_standardized(first_distribution, data_rng, int(n_worlds)),
                    FAMILIES[1]: _draw_standardized(second_distribution, data_rng, int(n_worlds)),
                }
                try:
                    result = simultaneous_max_t_world_intervals(
                        differences,
                        resamples=bootstrap_resamples,
                        seed=_bootstrap_seed(pair_name, int(n_worlds), repetition),
                        confidence_level=CONFIDENCE_LEVEL,
                        precision_target_half_width=1.0e9,
                    )
                except (ValueError, FloatingPointError):
                    inference_failures += 1
                    continue
                intervals = result["intervals_by_family"]
                joint_covered += int(
                    all(row["adjusted_lower"] <= 0.0 <= row["adjusted_upper"] for row in intervals.values())
                )
                positive_rule_false_positives += int(
                    all(row["adjusted_lower"] > 0.0 for row in intervals.values())
                )

            lower, upper = _wilson_interval(joint_covered, repetitions, CONFIDENCE_LEVEL)
            false_positive_lower, false_positive_upper = _wilson_interval(
                positive_rule_false_positives, repetitions, CONFIDENCE_LEVEL
            )
            observed_coverage = joint_covered / repetitions
            cell = {
                "distribution_pair": pair_name,
                "single_response_distribution": first_distribution,
                "mixture_response_distribution": second_distribution,
                "worlds_per_family": int(n_worlds),
                "monte_carlo_repetitions": int(repetitions),
                "bootstrap_resamples": int(bootstrap_resamples),
                "joint_coverage_successes": joint_covered,
                "joint_coverage_estimate": observed_coverage,
                "joint_coverage_wilson_95": [lower, upper],
                "positive_rule_false_positives_when_both_means_are_zero": positive_rule_false_positives,
                "positive_rule_false_positive_rate": positive_rule_false_positives / repetitions,
                "positive_rule_false_positive_rate_wilson_95": [false_positive_lower, false_positive_upper],
                "inference_failures": inference_failures,
                "screen_pass": lower >= JOINT_COVERAGE_LOWER_BOUND and inference_failures == 0,
            }
            cells.append(cell)
            if progress_callback is not None:
                progress_callback(cell)

    any_inference_failure = any(cell["inference_failures"] for cell in cells)
    passes = all(cell["screen_pass"] for cell in cells)
    if any_inference_failure:
        screen_result = "incomplete_inference_failure"
    elif passes:
        screen_result = "passes_predeclared_coverage_screen"
    else:
        screen_result = "fails_predeclared_coverage_screen"
    return {
        "record_type": "outcome-blind synthetic interval coverage screen",
        "audit_namespace": AUDIT_NAMESPACE,
        "protocol_status": "predeclared_screen_not_protocol_freeze",
        "patient_data_used": False,
        "trackc_policy_outcomes_used": False,
        "method": "world-resampled observed-SE max-T helper",
        "target": "joint coverage of both two-sided 95% intervals when both true family means are zero",
        "master_seed": MASTER_SEED,
        "rng": "numpy.random.Generator(PCG64DXSM); SeedSequence keyed by distribution-pair index, world count, and repetition",
        "bootstrap_seed_rule": "first 8 SHA-256 bytes, little-endian, from fixed audit namespace and scenario identifiers",
        "distribution_pairs": [
            {
                "name": name,
                "single_response": first,
                "mixture_response": second,
            }
            for name, first, second in DISTRIBUTION_PAIRS
        ],
        "world_counts": [int(value) for value in counts],
        "monte_carlo_repetitions": int(repetitions),
        "bootstrap_resamples": int(bootstrap_resamples),
        "confidence_level": CONFIDENCE_LEVEL,
        "screen_rule": f"every cell's Wilson lower bound must be >= {JOINT_COVERAGE_LOWER_BOUND:.2f} and have zero inference failures",
        "screen_result": screen_result,
        "all_cells_pass": passes,
        "cells": cells,
        "interpretation_limit": "Passing screens only these assumed distributions and world counts; it does not prove future family independence, validate the simulator, or imply patient or treatment benefit.",
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_record(progress_callback: Callable[[dict], None] | None = None) -> dict:
    root = Path(__file__).resolve().parents[1]
    started = perf_counter()
    record = run_coverage_audit(progress_callback=progress_callback)
    record["elapsed_seconds"] = perf_counter() - started
    record["fingerprints_sha256"] = {
        "coverage_audit_code": _sha256(Path(__file__).resolve()),
        "interval_helper_code": _sha256(root / "src" / "trackc_world_summary.py"),
        "protocol": _sha256(root / "paper" / "trackc_max_t_coverage_audit_protocol_2026-11-08.md"),
        "integrated_config": _sha256(root / "configs" / "trackABC_integrated_v2.json"),
    }
    return record


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    default_output = root / "output" / "trackc_max_t_coverage_audit_20261108" / "coverage.json"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=default_output)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    def report_cell(cell: dict) -> None:
        print(
            json.dumps(
                {
                    "distribution_pair": cell["distribution_pair"],
                    "worlds_per_family": cell["worlds_per_family"],
                    "joint_coverage": cell["joint_coverage_estimate"],
                    "wilson_lower": cell["joint_coverage_wilson_95"][0],
                    "screen_pass": cell["screen_pass"],
                }
            ),
            flush=True,
        )

    record = build_record(progress_callback=report_cell)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(record, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"screen_result": record["screen_result"], "cells": len(record["cells"]), "output": str(output)}))


if __name__ == "__main__":
    main()
