"""Run the amended ACRIN area forecast from the two published S0 tables."""

import argparse
import csv
import json
import sys
from collections import defaultdict

from itertools import groupby
from math import isfinite, log
from pathlib import Path
from random import Random
from statistics import mean, median, quantiles

# Isolated Python omits the script directory from its import path.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from track_b_acrin6677_forecast import predict_b1, predict_b2, predict_b3, predict_m, score

ROOT = Path(__file__).resolve().parents[1]
METHODS = ("B1", "B2", "B3", "M")
MIN_ELIGIBLE = 25


def number(value):
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if isfinite(numeric) else None


def build_windows(rows):
    visits = defaultdict(lambda: defaultdict(list))
    for row in rows:
        day = number(row["S0e2d"])
        patient = row["cn"].strip()
        if not patient or day is None or row["s0e5"].strip() == "1":
            continue
        value = number(row["s0e23"])
        if value is None and row["s0e13"].strip() == "1":
            value = 0.0
        if value is not None:
            visits[patient][day].append(value)
    return {
        patient: tuple((day, mean(days[day])) for day in sorted(days)[:3])
        for patient, days in sorted(visits.items())
        if len(days) >= 3
    }


def patient_errors(window):
    (t1, a1), (t2, a2), (t3, a3) = window
    y1, y2, y3 = (log(area + 0.1) for area in (a1, a2, a3))
    args = t1, y1, t2, y2, t3
    predictions = predict_b1(y1, y2), predict_b2(*args), predict_b3(*args), predict_m(*args)
    return dict(zip(METHODS, (score(pred, y3) for pred in predictions)))


def exact_wilcoxon(differences):
    """Drop zeros, assign average ranks to ties, and count sign assignments exactly."""
    nonzero = sorted((abs(d), d > 0) for d in differences if d != 0)
    counts = {0: 1}
    observed = 0
    first_rank = 1
    for _, group in groupby(nonzero, key=lambda pair: pair[0]):
        signs = [positive for _, positive in group]
        # Double average ranks so the permutation sums are integers.
        rank = 2 * first_rank + len(signs) - 1
        first_rank += len(signs)
        observed += rank * sum(signs)
        for _ in signs:
            for subtotal, count in list(counts.items()):
                counts[subtotal + rank] = counts.get(subtotal + rank, 0) + count
    total = max(counts)
    cutoff = min(observed, total - observed)
    extreme = sum(count for subtotal, count in counts.items()
                  if min(subtotal, total - subtotal) <= cutoff)
    return extreme / (1 << len(nonzero)), 2 * observed > total


def summarize(errors):
    n = len(errors)
    medians = {method: median(row[method] for row in errors) for method in METHODS}
    # Protocol lists B1, B2, B3 in this order. Use that order for equal medians.
    baseline = min(METHODS[:-1], key=medians.__getitem__)
    differences = [row["M"] - row[baseline] for row in errors]
    paired_median = median(differences)
    p, baseline_direction = exact_wilcoxon(differences)
    verdict = "inconclusive"
    if paired_median > 0 or (baseline_direction and p < 0.05):
        verdict = "failure"
    elif paired_median < 0 and p < 0.05:
        verdict = "success"
    rng = Random(20261112)
    boot = {method: [] for method in METHODS}
    for _ in range(10_000):
        sample = rng.choices(errors, k=n)
        for method in METHODS:
            boot[method].append(median(row[method] for row in sample))
    intervals = {}
    for method, estimates in boot.items():
        percentiles = quantiles(estimates, n=200, method="inclusive")
        intervals[method] = {"median_error": medians[method],
                             "bootstrap_95_ci": [percentiles[4], percentiles[194]]}
    return {"n": n, "methods": intervals, "best_baseline": baseline,
            "wilcoxon_p": p, "median_paired_difference": paired_median,
            "win_count": sum(d < 0 for d in differences), "verdict": verdict}


def read_rows():
    paths = sorted((ROOT / "data" / "acrin6677" / "x").rglob("S0.csv"))
    for split in ("Training-75", "Test-25"):
        matches = [path for path in paths if path.parent.name.endswith(split)]
        if len(matches) != 1:
            raise RuntimeError(f"Expected one S0 table for {split}")
        with matches[0].open(encoding="utf-8-sig", newline="") as stream:
            yield from csv.DictReader(stream)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    windows = build_windows(read_rows())
    if args.dry_run:
        print(len(windows))
        return
    if len(windows) < MIN_ELIGIBLE:
        raise SystemExit("Eligibility minimum not met.")
    today = "2026-11-11"
    output = ROOT / "output" / f"track_b_acrin6677_{today}"
    errors = [patient_errors(window) for window in windows.values()]
    result = summarize(errors)
    output.mkdir(parents=True, exist_ok=False)
    with (output / "per_patient.csv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("cn", *METHODS))
        writer.writerows((patient, *(row[method] for method in METHODS))
                         for patient, row in zip(windows, errors))
    (output / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

