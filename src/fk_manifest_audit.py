"""Fail-closed, aggregate-only audit for the Track B MU cohort manifests.

This module reads tabular manifests only. It does not open MRI or segmentation files and never includes
patient identifiers in its report. A blocked report means the cohort manifest must be reconciled before scoring.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_DIR = PROJECT_ROOT / "data" / "manifests"
_SCAN_FIELDS = {"patient_id", "dataset", "scan_number", "day", "mask_exists"}
_PATIENT_FIELDS = {"patient_id", "dataset", "eligible", "n_forecast_pairs", "n_primary_pairs"}
_PAIR_FIELDS = {"patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days",
                "input_core_nonempty", "in_primary"}
_SPLIT_FIELDS = {"patient_id", "fold"}


def _require_fields(rows, required: set[str], label: str) -> None:
    if not rows:
        raise ValueError(f"{label} rows must be nonempty")
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(sorted(missing))}")


def _identifier(row) -> str | None:
    value = row.get("patient_id")
    if value is None:
        return None
    value = str(value).strip()
    return value or None


def _number(value, *, integer: bool = False) -> float | int | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(result) or (integer and not result.is_integer()):
        return None
    return int(result) if integer else result


def _boolean(value) -> bool | None:
    if isinstance(value, bool):
        return value
    if value is None:
        return None
    normalized = str(value).strip().lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n"}:
        return False
    return None


def _mu_rows(rows, label: str):
    selected = [dict(row) for row in rows if str(row.get("dataset", "")).strip().upper() == "MU"]
    if not selected:
        raise ValueError(f"{label} has no MU rows")
    return selected


def audit_mu_manifests(scan_rows, patient_rows, pair_rows, split_rows, *,
                       source_reported_scan_count: int | None = None, n_folds: int = 5) -> dict:
    """Check aggregate consistency across scan, patient, pair, and fold manifests.

    The report contains counts and issue codes only. No IDs or row-level values leave this function.
    """
    scans = _mu_rows(list(scan_rows), "scan manifest")
    patients = _mu_rows(list(patient_rows), "patient manifest")
    pairs = [dict(row) for row in pair_rows]
    splits = [dict(row) for row in split_rows]
    _require_fields(scans, _SCAN_FIELDS, "MU scan manifest")
    _require_fields(patients, _PATIENT_FIELDS, "MU patient manifest")
    _require_fields(pairs, _PAIR_FIELDS, "MU forecast-pair manifest")
    _require_fields(splits, _SPLIT_FIELDS, "MU split manifest")
    if isinstance(n_folds, bool) or int(n_folds) != n_folds or n_folds < 2:
        raise ValueError("n_folds must be an integer of at least two")
    n_folds = int(n_folds)
    if source_reported_scan_count is not None and (
            isinstance(source_reported_scan_count, bool) or int(source_reported_scan_count) != source_reported_scan_count
            or source_reported_scan_count < 1):
        raise ValueError("source_reported_scan_count must be a positive integer")

    issues: Counter[str] = Counter()
    scan_map = {}
    scan_patient_ids = set()
    scan_days_missing = 0
    duplicate_scan_keys = 0
    for row in scans:
        patient_id = _identifier(row)
        scan_number = _number(row.get("scan_number"), integer=True)
        if patient_id is None or scan_number is None:
            issues["invalid_scan_rows"] += 1
            continue
        key = (patient_id, scan_number)
        if key in scan_map:
            duplicate_scan_keys += 1
        else:
            scan_map[key] = row
        scan_patient_ids.add(patient_id)
        if _number(row.get("day")) is None:
            scan_days_missing += 1
        if _boolean(row.get("mask_exists")) is None:
            issues["invalid_scan_mask_flags"] += 1
    if duplicate_scan_keys:
        issues["duplicate_scan_keys"] = duplicate_scan_keys
    if scan_days_missing:
        issues["scan_records_missing_day"] = scan_days_missing

    patient_map = {}
    duplicate_patient_rows = 0
    for row in patients:
        patient_id = _identifier(row)
        if patient_id is None:
            issues["invalid_patient_rows"] += 1
            continue
        if patient_id in patient_map:
            duplicate_patient_rows += 1
        else:
            patient_map[patient_id] = row
        if _boolean(row.get("eligible")) is None:
            issues["invalid_patient_eligibility_flags"] += 1
    if duplicate_patient_rows:
        issues["duplicate_patient_rows"] = duplicate_patient_rows

    split_map = {}
    duplicate_split_patients = 0
    invalid_folds = 0
    for row in splits:
        patient_id = _identifier(row)
        fold = _number(row.get("fold"), integer=True)
        if patient_id is None:
            issues["invalid_split_rows"] += 1
            continue
        if patient_id in split_map:
            duplicate_split_patients += 1
        else:
            split_map[patient_id] = fold
        if fold is None or fold < 0 or fold >= n_folds:
            invalid_folds += 1
    if duplicate_split_patients:
        issues["duplicate_split_patients"] = duplicate_split_patients
    if invalid_folds:
        issues["invalid_fold_values"] = invalid_folds

    observed_pairs = Counter()
    observed_primary = Counter()
    pair_keys = set()
    duplicate_pairs = 0
    missing_references = 0
    interval_mismatches = 0
    pair_day_mismatches = 0
    invalid_pair_rows = 0
    scan_order_errors = 0
    nonpositive_intervals = 0
    primary_pair_rows = 0
    primary_patients = set()
    primary_on_ineligible_rows = 0
    primary_on_ineligible_patients = set()
    primary_criteria_failures = 0
    pair_patient_ids = set()

    for row in pairs:
        patient_id = _identifier(row)
        tp_in = _number(row.get("tp_in"), integer=True)
        tp_out = _number(row.get("tp_out"), integer=True)
        day_in, day_out, delta = (_number(row.get(key)) for key in ("day_in", "day_out", "dt_days"))
        primary = _boolean(row.get("in_primary"))
        input_nonempty = _boolean(row.get("input_core_nonempty"))
        if (patient_id is None or tp_in is None or tp_out is None or day_in is None
                or day_out is None or delta is None or primary is None or input_nonempty is None):
            invalid_pair_rows += 1
            continue
        pair_patient_ids.add(patient_id)
        observed_pairs[patient_id] += 1
        if tp_out <= tp_in:
            scan_order_errors += 1
        if day_out <= day_in or delta <= 0:
            nonpositive_intervals += 1
        key = (patient_id, tp_in, tp_out)
        if key in pair_keys:
            duplicate_pairs += 1
        pair_keys.add(key)

        input_scan = scan_map.get((patient_id, tp_in))
        output_scan = scan_map.get((patient_id, tp_out))
        if input_scan is None or output_scan is None:
            missing_references += 1
        else:
            for supplied_day, scan_row in ((day_in, input_scan), (day_out, output_scan)):
                manifest_day = _number(scan_row.get("day"))
                if manifest_day is None or abs(supplied_day - manifest_day) > 1e-6:
                    pair_day_mismatches += 1
                    break
        if abs((day_out - day_in) - delta) > 1e-6:
            interval_mismatches += 1

        if primary:
            primary_pair_rows += 1
            primary_patients.add(patient_id)
            observed_primary[patient_id] += 1
            patient = patient_map.get(patient_id)
            if patient is not None and _boolean(patient.get("eligible")) is False:
                primary_on_ineligible_rows += 1
                primary_on_ineligible_patients.add(patient_id)
            input_scan_ok = input_scan is not None and _boolean(input_scan.get("mask_exists")) is True
            output_scan_ok = output_scan is not None and _boolean(output_scan.get("mask_exists")) is True
            if (not input_nonempty or not (14.0 <= delta <= 365.0)
                    or not input_scan_ok or not output_scan_ok):
                primary_criteria_failures += 1

    if invalid_pair_rows:
        issues["invalid_forecast_pair_rows"] = invalid_pair_rows
    if scan_order_errors:
        issues["forecast_pair_scan_order_errors"] = scan_order_errors
    if nonpositive_intervals:
        issues["nonpositive_pair_intervals"] = nonpositive_intervals
    if duplicate_pairs:
        issues["duplicate_forecast_pairs"] = duplicate_pairs
    if missing_references:
        issues["missing_pair_scan_references"] = missing_references
    if interval_mismatches:
        issues["pair_interval_mismatches"] = interval_mismatches
    if pair_day_mismatches:
        issues["pair_scan_day_mismatches"] = pair_day_mismatches
    if primary_criteria_failures:
        issues["primary_pairs_fail_declared_pair_criteria"] = primary_criteria_failures
    if primary_on_ineligible_rows:
        issues["primary_pairs_for_ineligible_patients"] = primary_on_ineligible_rows

    manifest_forecast_sum = 0
    manifest_primary_sum = 0
    forecast_count_mismatches = 0
    primary_count_mismatches = 0
    eligible_without_primary = set()
    for patient_id, row in patient_map.items():
        forecast_count = _number(row.get("n_forecast_pairs"), integer=True)
        primary_count = _number(row.get("n_primary_pairs"), integer=True)
        if forecast_count is None or primary_count is None:
            issues["invalid_patient_pair_counts"] += 1
            continue
        manifest_forecast_sum += forecast_count
        manifest_primary_sum += primary_count
        if forecast_count != observed_pairs[patient_id]:
            forecast_count_mismatches += 1
        if primary_count != observed_primary[patient_id]:
            primary_count_mismatches += 1
        if _boolean(row.get("eligible")) is True and observed_primary[patient_id] == 0:
            eligible_without_primary.add(patient_id)
    if forecast_count_mismatches:
        issues["patient_forecast_pair_count_mismatches"] = forecast_count_mismatches
    if primary_count_mismatches:
        issues["patient_primary_pair_count_mismatches"] = primary_count_mismatches
    if eligible_without_primary:
        issues["eligible_patients_without_primary_pairs"] = len(eligible_without_primary)
    missing_patient_rows = pair_patient_ids - set(patient_map)
    if missing_patient_rows:
        issues["pair_patients_missing_patient_manifest"] = len(missing_patient_rows)
    scans_missing_patient_rows = scan_patient_ids - set(patient_map)
    if scans_missing_patient_rows:
        issues["scan_patients_missing_patient_manifest"] = len(scans_missing_patient_rows)
    patients_without_scans = set(patient_map) - scan_patient_ids
    if patients_without_scans:
        issues["patient_manifest_rows_without_scans"] = len(patients_without_scans)
    if primary_patients - set(patient_map):
        issues["primary_pair_patients_missing_patient_manifest"] = len(primary_patients - set(patient_map))

    missing_split_rows = set(patient_map) - set(split_map)
    extra_split_rows = set(split_map) - set(patient_map)
    if missing_split_rows:
        issues["patients_missing_split"] = len(missing_split_rows)
    if extra_split_rows:
        issues["split_rows_without_patient_manifest"] = len(extra_split_rows)

    if source_reported_scan_count is not None and int(source_reported_scan_count) != len(scans):
        issues["local_source_scan_count_mismatch"] = 1

    ordered_issues = [{"code": code, "count": int(count)} for code, count in sorted(issues.items()) if count]
    summary = {
        "mu_scan_rows": len(scans),
        "mu_unique_scan_keys": len(scan_map),
        "mu_unique_scan_patients": len(scan_patient_ids),
        "scan_records_missing_day": scan_days_missing,
        "mu_patient_rows": len(patients),
        "mu_eligible_patients": sum(_boolean(row.get("eligible")) is True for row in patients),
        "forecast_pair_rows": len(pairs),
        "primary_pair_rows": primary_pair_rows,
        "patients_with_primary_pairs": len(primary_patients),
        "primary_pairs_on_ineligible_patients": primary_on_ineligible_rows,
        "ineligible_patients_with_primary_pairs": len(primary_on_ineligible_patients),
        "manifest_forecast_pair_sum": manifest_forecast_sum,
        "manifest_primary_pair_sum": manifest_primary_sum,
        "split_rows": len(splits),
        "n_folds": n_folds,
        "source_reported_scan_count": source_reported_scan_count,
    }
    return {"status": "blocked" if ordered_issues else "ready", "summary": summary, "issues": ordered_issues}


def _read_csv(path: Path) -> list[dict]:
    with Path(path).open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def audit_manifest_files(manifest_dir: Path = MANIFEST_DIR, *, source_reported_scan_count: int | None = None) -> dict:
    """Read the four local manifest files and return an aggregate-only audit report."""
    manifest_dir = Path(manifest_dir)
    return audit_mu_manifests(
        _read_csv(manifest_dir / "scan_manifest.csv"),
        _read_csv(manifest_dir / "patient_manifest.csv"),
        _read_csv(manifest_dir / "forecast_pairs_mu.csv"),
        _read_csv(manifest_dir / "split_mu.csv"),
        source_reported_scan_count=source_reported_scan_count,
    )


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--source-reported-scan-count", type=int)
    args = parser.parse_args(argv)
    report = audit_manifest_files(args.manifest_dir,
                                  source_reported_scan_count=args.source_reported_scan_count)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["status"] == "ready" else 2


if __name__ == "__main__":
    raise SystemExit(main())
