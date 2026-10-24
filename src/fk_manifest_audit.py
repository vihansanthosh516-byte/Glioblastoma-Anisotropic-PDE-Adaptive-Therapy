"""Fail-closed, aggregate-only audit for the Track B MU cohort manifests.

This module reads tabular manifests only. It does not open MRI or segmentation files and never includes
patient identifiers in its report. A blocked report means the cohort manifest must be reconciled before scoring.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent.parent
MANIFEST_DIR = PROJECT_ROOT / "data" / "manifests"
_SCAN_FIELDS = {"patient_id", "dataset", "scan_number", "day", "mask_exists"}
_PATIENT_FIELDS = {"patient_id", "dataset", "eligible", "n_forecast_pairs", "n_primary_pairs"}
_PAIR_FIELDS = {"patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days",
                "input_core_nonempty", "in_primary"}
_SPLIT_FIELDS = {"patient_id", "fold"}
_MRI_FILE_KINDS = {"t1n", "t1c", "t2w", "t2f"}
_REQUIRED_TIMEPOINT_FILE_KINDS = _MRI_FILE_KINDS | {"segmentation"}


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
    patient_usable_scan_counts = Counter()
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
            if (_number(row.get("day")) is not None
                    and _boolean(row.get("mask_exists")) is True):
                patient_usable_scan_counts[patient_id] += 1
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
        "mu_patients_with_4plus_day_and_mask_flag_scans": sum(
            count >= 4 for count in patient_usable_scan_counts.values()
        ),
        "eligible_mu_patients_with_4plus_manifest_flagged_scans": sum(
            _boolean(row.get("eligible")) is True and patient_usable_scan_counts[patient_id] >= 4
            for patient_id, row in patient_map.items()
        ),
        "eligible_mu_patients_with_5plus_manifest_flagged_scans": sum(
            _boolean(row.get("eligible")) is True and patient_usable_scan_counts[patient_id] >= 5
            for patient_id, row in patient_map.items()
        ),
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


def _local_mri_file_kind(filename: str) -> str:
    name = str(filename).lower()
    if re.search(r"seg|mask|label", name):
        return "segmentation"
    if re.search(r"t1n|native.?t1", name):
        return "t1n"
    if re.search(r"t1c|t1.?ce|contrast", name):
        return "t1c"
    if re.search(r"t2f|flair", name):
        return "t2f"
    if re.search(r"t2w|t2", name):
        return "t2w"
    return "unknown"


def audit_mu_local_image_tree(scan_rows, image_dir: Path, *,
                              source_reported_timepoint_count: int | None = None,
                              source_reported_image_count: int | None = None) -> dict:
    """Check MU timepoint paths and NIfTI filenames without opening image contents.

    A complete manifest-listed timepoint must have four named MRI sequences plus one segmentation. An
    unmanifested local timepoint is accepted only when it has exactly the four MRI sequence files and no
    segmentation. The report contains aggregate counts and issue codes only; it cannot validate NIfTI bytes,
    image geometry, segmentation values, or clinical suitability.
    """
    scans = _mu_rows(list(scan_rows), "scan manifest")
    _require_fields(scans, _SCAN_FIELDS, "MU scan manifest")
    image_dir = Path(image_dir)
    for name, value in (("source_reported_timepoint_count", source_reported_timepoint_count),
                        ("source_reported_image_count", source_reported_image_count)):
        if value is not None and (isinstance(value, bool) or int(value) != value or value < 0):
            raise ValueError(f"{name} must be a nonnegative integer")

    issues: Counter[str] = Counter()
    if not image_dir.is_dir():
        return {
            "status": "blocked",
            "summary": {"image_directory_present": False},
            "issues": [{"code": "local_image_directory_missing", "count": 1}],
            "interpretation_limit": "no image directory metadata was audited",
        }

    manifest_by_key = {}
    manifest_patient_ids = set()
    for row in scans:
        patient_id = _identifier(row)
        scan_number = _number(row.get("scan_number"), integer=True)
        if patient_id is None or scan_number is None or scan_number < 1:
            issues["invalid_manifest_patient_or_scan_number"] += 1
            continue
        key = (patient_id, int(scan_number))
        if key in manifest_by_key:
            issues["duplicate_manifest_patient_scan_keys"] += 1
            continue
        manifest_by_key[key] = row
        manifest_patient_ids.add(patient_id)

    patient_dirs = [path for path in image_dir.iterdir() if path.is_dir()]
    patient_dirs_by_id = {path.name: path for path in patient_dirs}
    if len(patient_dirs_by_id) != len(patient_dirs):
        issues["duplicate_local_patient_directory_names"] += len(patient_dirs) - len(patient_dirs_by_id)
    missing_patient_dirs = manifest_patient_ids - set(patient_dirs_by_id)
    extra_patient_dirs = set(patient_dirs_by_id) - manifest_patient_ids
    if missing_patient_dirs:
        issues["manifest_patients_without_local_directories"] += len(missing_patient_dirs)
    if extra_patient_dirs:
        issues["local_patient_directories_without_manifest_rows"] += len(extra_patient_dirs)

    local_timepoints = {}
    file_counts_by_timepoint: Counter[int] = Counter()
    modality_file_counts: Counter[str] = Counter()
    total_nifti_files = 0
    for patient_id, patient_dir in patient_dirs_by_id.items():
        for timepoint_dir in patient_dir.rglob("Timepoint_*"):
            match = re.fullmatch(r"Timepoint_(\d+)", timepoint_dir.name)
            if not timepoint_dir.is_dir() or match is None:
                continue
            timepoint_number = int(match.group(1))
            key = (patient_id, timepoint_number)
            if key in local_timepoints:
                issues["duplicate_local_timepoint_directories"] += 1
                continue
            nifti_files = [path for path in timepoint_dir.iterdir()
                           if path.is_file() and path.name.lower().endswith(".nii.gz")]
            kinds = [_local_mri_file_kind(path.name) for path in nifti_files]
            kind_counts = Counter(kinds)
            local_timepoints[key] = {"nifti_files": len(nifti_files), "kind_counts": kind_counts}
            file_counts_by_timepoint[len(nifti_files)] += 1
            modality_file_counts.update(kinds)
            total_nifti_files += len(nifti_files)

    manifest_rows_matching = 0
    manifest_rows_complete = 0
    for key, row in manifest_by_key.items():
        item = local_timepoints.get(key)
        if item is None:
            issues["manifest_timepoint_directory_missing"] += 1
            continue
        manifest_rows_matching += 1
        kind_counts = item["kind_counts"]
        if any(kind_counts[kind] != 1 for kind in _MRI_FILE_KINDS):
            issues["manifest_timepoint_missing_required_mri_sequence"] += 1
        if _boolean(row.get("mask_exists")) is True and kind_counts["segmentation"] != 1:
            issues["manifest_flagged_segmentation_missing"] += 1
        if (item["nifti_files"] == 5
                and all(kind_counts[kind] == 1 for kind in _REQUIRED_TIMEPOINT_FILE_KINDS)):
            manifest_rows_complete += 1
        if kind_counts["unknown"] or any(value > 1 for value in kind_counts.values()):
            issues["manifest_timepoint_has_unrecognized_or_duplicate_files"] += 1

    extra_keys = set(local_timepoints) - set(manifest_by_key)
    unmanifested_incomplete = 0
    for key in extra_keys:
        item = local_timepoints[key]
        kind_counts = item["kind_counts"]
        if (item["nifti_files"] == 4 and kind_counts["segmentation"] == 0
                and kind_counts["unknown"] == 0
                and all(kind_counts[kind] == 1 for kind in _MRI_FILE_KINDS)):
            unmanifested_incomplete += 1
        else:
            issues["unmanifested_local_timepoint_not_expected_incomplete"] += 1

    if (source_reported_timepoint_count is not None
            and source_reported_timepoint_count != len(local_timepoints)):
        issues["local_timepoint_directory_count_mismatch"] += 1
    if source_reported_image_count is not None and source_reported_image_count != total_nifti_files:
        issues["local_nifti_file_count_mismatch"] += 1
    if (source_reported_timepoint_count is not None
            and source_reported_timepoint_count != len(manifest_by_key) + unmanifested_incomplete):
        issues["source_timepoint_count_not_explained_by_manifest_and_incomplete_dirs"] += 1

    ordered_issues = [{"code": code, "count": int(count)} for code, count in sorted(issues.items()) if count]
    return {
        "status": "blocked" if ordered_issues else "consistent_by_file_metadata",
        "summary": {
            "image_directory_present": True,
            "local_patient_directories": len(patient_dirs),
            "manifest_mu_patients": len(manifest_patient_ids),
            "manifest_mu_rows": len(manifest_by_key),
            "local_timepoint_directories": len(local_timepoints),
            "manifest_rows_matching_local_timepoints": manifest_rows_matching,
            "manifest_rows_with_all_four_mri_plus_seg": manifest_rows_complete,
            "unmanifested_local_timepoints": len(extra_keys),
            "unmanifested_timepoints_with_four_mri_no_seg": unmanifested_incomplete,
            "local_nifti_file_count": total_nifti_files,
            "nifti_file_counts_per_timepoint": dict(sorted(file_counts_by_timepoint.items())),
            "modality_file_counts": dict(sorted(modality_file_counts.items())),
            "source_reported_timepoint_count": source_reported_timepoint_count,
            "source_reported_image_count": source_reported_image_count,
        },
        "issues": ordered_issues,
        "interpretation_limit": (
            "file and directory names only; no NIfTI bytes, image contents, geometry, segmentation values, "
            "or clinical suitability were checked"
        ),
    }


def audit_manifest_files(manifest_dir: Path = MANIFEST_DIR, *,
                         source_reported_scan_count: int | None = None,
                         source_reported_image_count: int | None = None,
                         local_image_dir: Path | None = None) -> dict:
    """Read local manifests; optionally reconcile them to local TCIA file metadata."""
    manifest_dir = Path(manifest_dir)
    scan_rows = _read_csv(manifest_dir / "scan_manifest.csv")
    report = audit_mu_manifests(
        scan_rows,
        _read_csv(manifest_dir / "patient_manifest.csv"),
        _read_csv(manifest_dir / "forecast_pairs_mu.csv"),
        _read_csv(manifest_dir / "split_mu.csv"),
        source_reported_scan_count=source_reported_scan_count,
    )
    if local_image_dir is not None:
        release_audit = audit_mu_local_image_tree(
            scan_rows, local_image_dir,
            source_reported_timepoint_count=source_reported_scan_count,
            source_reported_image_count=source_reported_image_count,
        )
        report["local_release_audit"] = release_audit
        if release_audit["status"] == "consistent_by_file_metadata":
            report["issues"] = [issue for issue in report["issues"]
                                if issue["code"] != "local_source_scan_count_mismatch"]
            report["summary"]["source_unsegmented_timepoint_count"] = (
                release_audit["summary"]["unmanifested_timepoints_with_four_mri_no_seg"]
            )
            report["status"] = "blocked" if report["issues"] else "ready"
        else:
            report["issues"].append({"code": "local_release_file_metadata_mismatch", "count": 1})
            report["status"] = "blocked"
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path, default=MANIFEST_DIR)
    parser.add_argument("--source-reported-scan-count", type=int,
                        help="source total timepoint count (legacy option name)")
    parser.add_argument("--source-reported-image-count", type=int)
    parser.add_argument("--local-image-dir", type=Path)
    args = parser.parse_args(argv)
    report = audit_manifest_files(args.manifest_dir,
                                  source_reported_scan_count=args.source_reported_scan_count,
                                  source_reported_image_count=args.source_reported_image_count,
                                  local_image_dir=args.local_image_dir)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["status"] == "ready" else 2


if __name__ == "__main__":
    raise SystemExit(main())
