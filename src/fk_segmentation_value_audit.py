"""Aggregate-only audit of MU-Glioma-Post segmentation labels and manifest core counts.

Unlike the header audit, this module reads one segmentation array at a time. It never reads MRI intensity
arrays and never returns patient identifiers, filenames, or row-level values. The report checks file contents;
it does not establish segmentation accuracy or clinical validity.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from pathlib import Path
import re
import sys

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.fk_manifest_audit import _boolean, _identifier, _local_mri_file_kind, _mu_rows, _number, _read_csv


DEFAULT_MANIFEST_DIR = PROJECT_ROOT / "data" / "manifests"
DEFAULT_IMAGE_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "track_b_segmentation_value_audit_20261024.json"
_EXPECTED_LABELS = frozenset({0, 1, 2, 3, 4})
_CORE_LABELS = frozenset({1, 3})
_INTEGER_TOLERANCE = 1e-6
_PRIMARY_MIN_DAYS = 14.0
_PRIMARY_MAX_DAYS = 365.0


def _local_timepoints(image_dir: Path, issues: Counter[str]) -> dict[tuple[str, int], Path]:
    result: dict[tuple[str, int], Path] = {}
    if not image_dir.is_dir():
        issues["local_image_directory_missing"] += 1
        return result
    for patient_dir in (path for path in image_dir.iterdir() if path.is_dir()):
        for timepoint_dir in patient_dir.rglob("Timepoint_*"):
            match = re.fullmatch(r"Timepoint_(\d+)", timepoint_dir.name)
            if not timepoint_dir.is_dir() or match is None:
                continue
            key = (patient_dir.name, int(match.group(1)))
            if key in result:
                issues["duplicate_local_timepoint_keys"] += 1
            else:
                result[key] = timepoint_dir
    return result


def _read_mask(path: Path, nib):
    """Read one label volume and return only aggregate diagnostics and its core count."""
    image = nib.load(str(path), mmap=False)
    if len(image.shape) != 3:
        return {"issue": "segmentation_not_3d", "core_voxels": None, "label_counts": {}}
    data = np.asanyarray(image.dataobj)
    if not np.isfinite(data).all():
        return {"issue": "segmentation_has_nonfinite_values", "core_voxels": None, "label_counts": {}}

    if np.issubdtype(data.dtype, np.integer):
        labels = data
        fractional_voxels = 0
    else:
        rounded = np.rint(data)
        fractional_voxels = int(np.count_nonzero(np.abs(data - rounded) > _INTEGER_TOLERANCE))
        if fractional_voxels:
            return {
                "issue": "segmentation_has_noninteger_values",
                "core_voxels": None,
                "label_counts": {},
                "fractional_voxels": fractional_voxels,
            }
        labels = rounded

    minimum, maximum = int(labels.min()), int(labels.max())
    if 0 <= minimum <= maximum <= 10_000:
        # Label volumes normally contain only a few small nonnegative integers. A histogram is linear-time;
        # sorting every voxel with unique() would make a whole-cohort audit needlessly slow.
        histogram = np.bincount(labels.reshape(-1).astype(np.int64, copy=False), minlength=maximum + 1)
        label_counts = {label: int(count) for label, count in enumerate(histogram) if count}
    else:
        unique, counts = np.unique(labels, return_counts=True)
        label_counts = {int(label): int(count) for label, count in zip(unique, counts)}
    unexpected = sorted(set(label_counts) - _EXPECTED_LABELS)
    core_voxels = sum(label_counts.get(label, 0) for label in _CORE_LABELS)
    del data, labels
    if unexpected:
        return {
            "issue": "segmentation_has_unexpected_labels",
            "core_voxels": core_voxels,
            "label_counts": label_counts,
            "unexpected_label_values": unexpected[:16],
            "unexpected_label_value_count": len(unexpected),
        }
    return {"issue": None, "core_voxels": core_voxels, "label_counts": label_counts}


def audit_mu_segmentation_values(manifest_dir: Path = DEFAULT_MANIFEST_DIR,
                                 image_dir: Path = DEFAULT_IMAGE_DIR) -> dict:
    """Audit local segmentation labels against the declared 1+3 core count in the scan manifest.

    The function reads segmentation arrays only, one at a time. All returned fields are aggregate counts;
    identifiers and paths are used internally for matching but are not emitted.
    """
    try:
        import nibabel as nib
    except ImportError as exc:
        raise RuntimeError("nibabel is required for the segmentation value audit") from exc

    manifest_dir = Path(manifest_dir)
    image_dir = Path(image_dir)
    scan_rows = _mu_rows(_read_csv(manifest_dir / "scan_manifest.csv"), "scan manifest")
    patient_rows = _mu_rows(_read_csv(manifest_dir / "patient_manifest.csv"), "patient manifest")
    pair_rows = _read_csv(manifest_dir / "forecast_pairs_mu.csv")
    issues: Counter[str] = Counter()
    timepoint_paths = _local_timepoints(image_dir, issues)

    scans_by_key: dict[tuple[str, int], dict] = {}
    scans_by_patient: dict[str, list[dict]] = defaultdict(list)
    for row in scan_rows:
        patient_id = _identifier(row)
        scan_number = _number(row.get("scan_number"), integer=True)
        if patient_id is None or scan_number is None or scan_number < 1:
            issues["invalid_manifest_scan_rows"] += 1
            continue
        key = (patient_id, int(scan_number))
        if key in scans_by_key:
            issues["duplicate_manifest_scan_keys"] += 1
            continue
        item = {
            "patient_id": patient_id,
            "scan_number": int(scan_number),
            "day": _number(row.get("day")),
            "mask_flag": _boolean(row.get("mask_exists")),
            "stored_core_voxels": _number(row.get("core_voxels_native"), integer=True),
            "local_core_voxels": None,
            "valid_labels": False,
        }
        scans_by_key[key] = item
        scans_by_patient[patient_id].append(item)
        if item["mask_flag"] is None:
            issues["invalid_manifest_mask_flags"] += 1

    patient_eligibility: dict[str, bool] = {}
    for row in patient_rows:
        patient_id = _identifier(row)
        if patient_id is None:
            issues["invalid_manifest_patient_rows"] += 1
            continue
        eligible = _boolean(row.get("eligible"))
        if patient_id in patient_eligibility:
            issues["duplicate_manifest_patient_rows"] += 1
        else:
            patient_eligibility[patient_id] = eligible is True
        if eligible is None:
            issues["invalid_manifest_patient_eligibility"] += 1

    pair_primary: dict[tuple[str, int, int], bool] = {}
    for row in pair_rows:
        patient_id = _identifier(row)
        tp_in = _number(row.get("tp_in"), integer=True)
        tp_out = _number(row.get("tp_out"), integer=True)
        primary = _boolean(row.get("in_primary"))
        if patient_id is None or tp_in is None or tp_out is None or primary is None:
            issues["invalid_manifest_pair_rows"] += 1
            continue
        key = (patient_id, int(tp_in), int(tp_out))
        if key in pair_primary:
            issues["duplicate_manifest_pair_keys"] += 1
        else:
            pair_primary[key] = primary

    counters: Counter[str] = Counter()
    total_label_counts: Counter[int] = Counter()
    unexpected_values: Counter[int] = Counter()
    fractional_voxels_total = 0
    for key, item in scans_by_key.items():
        if item["mask_flag"] is not True:
            continue
        counters["manifest_scans_flagged_with_segmentation"] += 1
        timepoint_dir = timepoint_paths.get(key)
        if timepoint_dir is None:
            issues["flagged_scan_timepoint_directory_missing"] += 1
            continue
        segmentation_paths = [
            path for path in timepoint_dir.iterdir()
            if path.is_file() and path.name.lower().endswith(".nii.gz")
            and _local_mri_file_kind(path.name) == "segmentation"
        ]
        if len(segmentation_paths) != 1:
            issues["flagged_scan_missing_or_duplicate_segmentation_files"] += 1
            continue
        counters["segmentation_arrays_attempted"] += 1
        try:
            result = _read_mask(segmentation_paths[0], nib)
        except Exception:
            issues["unreadable_segmentation_arrays"] += 1
            continue
        if result.get("issue"):
            issues[result["issue"]] += 1
        fractional_voxels_total += int(result.get("fractional_voxels", 0))
        total_label_counts.update(result.get("label_counts", {}))
        for label in result.get("unexpected_label_values", []):
            unexpected_values[label] += result["label_counts"].get(label, 0)
        if result["core_voxels"] is None:
            continue
        item["local_core_voxels"] = int(result["core_voxels"])
        item["valid_labels"] = result.get("issue") is None
        if item["stored_core_voxels"] is None:
            issues["manifest_core_voxel_counts_missing"] += 1
        else:
            counters["manifest_core_voxel_counts_compared"] += 1
            if item["stored_core_voxels"] != item["local_core_voxels"]:
                issues["manifest_core_voxel_count_mismatches"] += 1
            else:
                counters["manifest_core_voxel_counts_matched"] += 1
        if item["local_core_voxels"] == 0:
            counters["empty_source_core_masks"] += 1
        if item["valid_labels"]:
            counters["segmentation_arrays_with_integer_expected_label_map"] += 1

    # Count first-two-to-third candidates without skipping a broken or empty visit.
    eligible_patients = {pid for pid, eligible in patient_eligibility.items() if eligible}
    candidate_counts = Counter()
    primary_flag_mismatches = 0
    for patient_id in eligible_patients:
        visits = sorted(scans_by_patient.get(patient_id, []), key=lambda row: row["scan_number"])
        if len(visits) < 3:
            continue
        candidate_counts["eligible_patients_with_at_least_three_manifest_visits"] += 1
        first_three = visits[:3]
        if any(row["mask_flag"] is not True or row["local_core_voxels"] is None
               or not row["valid_labels"] or row["local_core_voxels"] <= 0
               for row in first_three):
            continue
        candidate_counts["eligible_patients_with_three_valid_nonempty_first_visits"] += 1
        days = [row["day"] for row in first_three]
        if any(day is None for day in days) or any(b <= a for a, b in zip(days, days[1:])):
            candidate_counts["three_valid_masks_with_missing_or_nonincreasing_days"] += 1
            continue
        gap = days[2] - days[1]
        pair_key = (patient_id, first_three[1]["scan_number"], first_three[2]["scan_number"])
        expected_primary = (
            first_three[1]["local_core_voxels"] > 0
            and _PRIMARY_MIN_DAYS <= gap <= _PRIMARY_MAX_DAYS
            and first_three[2]["local_core_voxels"] is not None
        )
        supplied_primary = pair_primary.get(pair_key)
        if supplied_primary is not None and supplied_primary != expected_primary:
            primary_flag_mismatches += 1
        if not (_PRIMARY_MIN_DAYS <= gap <= _PRIMARY_MAX_DAYS):
            candidate_counts["three_valid_masks_with_target_gap_outside_primary_window"] += 1
            continue
        candidate_counts["eligible_patients_with_first_two_to_third_in_primary_window"] += 1
        if supplied_primary is True:
            candidate_counts["eligible_patients_with_first_two_to_third_manifest_primary_pair"] += 1

    if primary_flag_mismatches:
        issues["first_two_to_third_primary_flag_mismatches"] += primary_flag_mismatches

    blocking_prefixes = (
        "local_image_directory_missing", "invalid_manifest", "duplicate_manifest",
        "flagged_scan", "unreadable_segmentation", "segmentation_", "manifest_core_voxel",
        "first_two_to_third_primary_flag_mismatches",
    )
    blocking_issues = [
        {"code": code, "count": int(count)}
        for code, count in sorted(issues.items())
        if count and code.startswith(blocking_prefixes)
    ]
    all_issues = [{"code": code, "count": int(count)}
                  for code, count in sorted(issues.items()) if count]
    return {
        "status": "blocked" if blocking_issues else "consistent_with_local_manifest_counts",
        "summary": {
            **{key: int(value) for key, value in sorted(counters.items())},
            "local_timepoint_directories": len(timepoint_paths),
            "mu_manifest_scan_rows": len(scans_by_key),
            "mu_manifest_patient_rows": len(patient_eligibility),
            "mu_manifest_eligible_patients": len(eligible_patients),
            **{key: int(value) for key, value in sorted(candidate_counts.items())},
            "total_voxel_counts_by_label": {
                str(label): int(total_label_counts.get(label, 0))
                for label in sorted(_EXPECTED_LABELS)
            },
            "unexpected_label_value_count": len(unexpected_values),
            "unexpected_label_voxel_count": int(sum(unexpected_values.values())),
            "unexpected_label_values_first_16": [int(v) for v in sorted(unexpected_values)[:16]],
            "fractional_voxels_total": int(fractional_voxels_total),
            "primary_interval_days": [_PRIMARY_MIN_DAYS, _PRIMARY_MAX_DAYS],
        },
        "issues": all_issues,
        "forecast_cohort_status": "provisional_not_ready_for_patient_scoring",
        "interpretation_limit": (
            "Segmentation values and labels 1+3 counts were checked against local manifest values only. "
            "This does not validate segmentation accuracy, biological boundary truth, source-release identity, "
            "clinical suitability, patient-fold provenance, or forecast performance. The strict cohort gate "
            "remains separate and unresolved."
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest-dir", type=Path, default=DEFAULT_MANIFEST_DIR)
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGE_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    report = audit_mu_segmentation_values(args.manifest_dir, args.image_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["status"] == "consistent_with_local_manifest_counts" else 2


if __name__ == "__main__":
    raise SystemExit(main())
