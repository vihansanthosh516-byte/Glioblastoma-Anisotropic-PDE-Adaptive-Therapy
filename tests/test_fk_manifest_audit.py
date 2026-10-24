from __future__ import annotations

import json
import csv
from pathlib import Path

import pytest

from src.fk_manifest_audit import audit_manifest_files, audit_mu_local_image_tree, audit_mu_manifests


def valid_rows():
    scans = [
        {"patient_id": "MU-1", "dataset": "MU", "scan_number": "1", "day": "0", "mask_exists": "True"},
        {"patient_id": "MU-1", "dataset": "MU", "scan_number": "2", "day": "30", "mask_exists": "True"},
    ]
    patients = [{
        "patient_id": "MU-1", "dataset": "MU", "eligible": "True",
        "n_scans": "2", "n_forecast_pairs": "1", "n_primary_pairs": "1",
    }]
    pairs = [{
        "patient_id": "MU-1", "tp_in": "1", "tp_out": "2", "day_in": "0",
        "day_out": "30", "dt_days": "30", "input_core_nonempty": "True", "in_primary": "True",
    }]
    splits = [{"patient_id": "MU-1", "fold": "0"}]
    return scans, patients, pairs, splits


def test_consistent_manifests_are_ready_and_report_only_aggregate_counts():
    report = audit_mu_manifests(*valid_rows())

    assert report["status"] == "ready"
    assert report["issues"] == []
    assert report["summary"]["mu_scan_rows"] == 2
    assert report["summary"]["mu_patient_rows"] == 1
    assert report["summary"]["primary_pair_rows"] == 1
    assert "MU-1" not in json.dumps(report)


def test_primary_pair_on_ineligible_patient_blocks_scoring():
    scans, patients, pairs, splits = valid_rows()
    patients[0]["eligible"] = "False"
    patients[0]["reason_excluded"] = "missing scan day"

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert report["status"] == "blocked"
    assert codes["primary_pairs_for_ineligible_patients"] == 1


def test_missing_pair_scan_link_and_bad_interval_both_block():
    scans, patients, pairs, splits = valid_rows()
    pairs[0]["tp_out"] = "3"
    pairs[0]["day_out"] = "31"
    pairs[0]["dt_days"] = "30"

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert report["status"] == "blocked"
    assert codes["missing_pair_scan_references"] == 1
    assert codes["pair_interval_mismatches"] == 1


def test_duplicate_scans_and_split_rows_block():
    scans, patients, pairs, splits = valid_rows()
    scans.append(dict(scans[0]))
    splits.append(dict(splits[0]))

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert report["status"] == "blocked"
    assert codes["duplicate_scan_keys"] == 1
    assert codes["duplicate_split_patients"] == 1


def test_source_scan_count_difference_is_reported_as_a_gate():
    report = audit_mu_manifests(*valid_rows(), source_reported_scan_count=4)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert report["status"] == "blocked"
    assert report["summary"]["source_reported_scan_count"] == 4
    assert report["summary"]["mu_scan_rows"] == 2
    assert codes["local_source_scan_count_mismatch"] == 1


def test_primary_pair_with_empty_input_is_blocked():
    scans, patients, pairs, splits = valid_rows()
    pairs[0]["input_core_nonempty"] = "False"

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["primary_pairs_fail_declared_pair_criteria"] == 1


def test_patient_pair_counts_must_match_pair_manifest():
    scans, patients, pairs, splits = valid_rows()
    patients[0]["n_primary_pairs"] = "0"

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["patient_primary_pair_count_mismatches"] == 1


def test_audit_reports_aggregate_count_of_eligible_four_visit_patients():
    scans, patients, pairs, splits = valid_rows()
    scans.extend([
        {"patient_id": "MU-1", "dataset": "MU", "scan_number": "3", "day": "60", "mask_exists": "True"},
        {"patient_id": "MU-1", "dataset": "MU", "scan_number": "4", "day": "90", "mask_exists": "True"},
    ])
    patients[0].update({"n_scans": "4", "n_forecast_pairs": "3", "n_primary_pairs": "3"})
    pairs.extend([
        {"patient_id": "MU-1", "tp_in": "2", "tp_out": "3", "day_in": "30",
         "day_out": "60", "dt_days": "30", "input_core_nonempty": "True", "in_primary": "True"},
        {"patient_id": "MU-1", "tp_in": "3", "tp_out": "4", "day_in": "60",
         "day_out": "90", "dt_days": "30", "input_core_nonempty": "True", "in_primary": "True"},
    ])

    report = audit_mu_manifests(scans, patients, pairs, splits)

    assert report["summary"]["eligible_mu_patients_with_4plus_manifest_flagged_scans"] == 1
    assert "MU-1" not in json.dumps(report)


def test_missing_days_and_unmatched_split_rows_are_reported():
    scans, patients, pairs, splits = valid_rows()
    scans[1]["day"] = ""
    splits[0]["fold"] = "8"
    splits.append({"patient_id": "MU-extra", "fold": "0"})

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["scan_records_missing_day"] == 1
    assert codes["pair_scan_day_mismatches"] == 1
    assert codes["invalid_fold_values"] == 1
    assert codes["split_rows_without_patient_manifest"] == 1


def test_eligible_patient_without_primary_pair_is_reported():
    scans, patients, pairs, splits = valid_rows()
    pairs[0]["in_primary"] = "False"
    patients[0]["n_primary_pairs"] = "0"

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["eligible_patients_without_primary_pairs"] == 1


def test_duplicate_patients_and_invalid_eligibility_flags_are_blocked():
    scans, patients, pairs, splits = valid_rows()
    patients[0]["eligible"] = "maybe"
    patients.append(dict(patients[0]))

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["duplicate_patient_rows"] == 1
    assert codes["invalid_patient_eligibility_flags"] == 2


def test_pair_for_unlisted_patient_is_not_silent():
    scans, patients, pairs, splits = valid_rows()
    orphan = dict(pairs[0])
    orphan["patient_id"] = "MU-orphan"
    pairs.append(orphan)

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["pair_patients_missing_patient_manifest"] == 1
    assert codes["missing_pair_scan_references"] == 1


def _make_timepoint(image_root, patient_id, number, *, include_segmentation=True):
    folder = image_root / patient_id / f"Timepoint_{number}"
    folder.mkdir(parents=True)
    names = ["t1n.nii.gz", "t1c.nii.gz", "t2w.nii.gz", "t2f.nii.gz"]
    if include_segmentation:
        names.append("segmentation.nii.gz")
    for name in names:
        (folder / name).touch()


def _write_manifest_files(manifest_dir, rows):
    manifest_dir.mkdir()
    scans, patients, pairs, splits = rows
    for name, values in (
        ("scan_manifest.csv", scans), ("patient_manifest.csv", patients),
        ("forecast_pairs_mu.csv", pairs), ("split_mu.csv", splits),
    ):
        with (manifest_dir / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(values[0]))
            writer.writeheader()
            writer.writerows(values)


def test_local_release_metadata_resolves_two_unsegmented_timepoints(tmp_path):
    rows = valid_rows()
    image_root = tmp_path / "images"
    _make_timepoint(image_root, "MU-1", 1)
    _make_timepoint(image_root, "MU-1", 2)
    _make_timepoint(image_root, "MU-1", 3, include_segmentation=False)
    manifest_dir = tmp_path / "manifests"
    _write_manifest_files(manifest_dir, rows)

    report = audit_manifest_files(
        manifest_dir, source_reported_scan_count=3, source_reported_image_count=14,
        local_image_dir=image_root,
    )

    assert report["status"] == "ready"
    assert report["issues"] == []
    assert report["summary"]["mu_scan_rows"] == 2
    assert report["summary"]["source_unsegmented_timepoint_count"] == 1
    assert report["local_release_audit"]["status"] == "consistent_by_file_metadata"
    assert report["local_release_audit"]["summary"]["manifest_rows_with_all_four_mri_plus_seg"] == 2
    assert "MU-1" not in json.dumps(report)


def test_local_release_metadata_audit_blocks_manifest_mask_without_segmentation(tmp_path):
    scans, _, _, _ = valid_rows()
    image_root = tmp_path / "images"
    _make_timepoint(image_root, "MU-1", 1, include_segmentation=False)

    report = audit_mu_local_image_tree(
        scans, image_root, source_reported_timepoint_count=1, source_reported_image_count=4,
    )
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert report["status"] == "blocked"
    assert codes["manifest_flagged_segmentation_missing"] == 1
    assert "MU-1" not in json.dumps(report)


def test_missing_required_columns_and_invalid_fold_count_raise():
    scans, patients, pairs, splits = valid_rows()
    with pytest.raises(ValueError, match="missing required columns"):
        audit_mu_manifests([{"patient_id": "MU-1", "dataset": "MU"}], patients, pairs, splits)
    with pytest.raises(ValueError, match="n_folds"):
        audit_mu_manifests(scans, patients, pairs, splits, n_folds=1)
    with pytest.raises(ValueError, match="source_reported_scan_count"):
        audit_mu_manifests(scans, patients, pairs, splits, source_reported_scan_count=0)


def test_manifest_file_reader_loads_all_four_tables(tmp_path: Path):
    scans, patients, pairs, splits = valid_rows()
    for name, rows in (
        ("scan_manifest.csv", scans),
        ("patient_manifest.csv", patients),
        ("forecast_pairs_mu.csv", pairs),
        ("split_mu.csv", splits),
    ):
        with (tmp_path / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    from src.fk_manifest_audit import audit_manifest_files, main

    report = audit_manifest_files(tmp_path)
    assert report["status"] == "ready"
    assert report["summary"]["split_rows"] == 1
    assert main(["--manifest-dir", str(tmp_path), "--source-reported-scan-count", "2"]) == 0


def test_scan_and_patient_manifests_must_cover_each_other():
    scans, patients, pairs, splits = valid_rows()
    scan_only = dict(scans[0])
    scan_only["patient_id"] = "MU-scan-only"
    scans.append(scan_only)
    patient_only = dict(patients[0])
    patient_only.update({"patient_id": "MU-patient-only", "eligible": "False",
                         "n_scans": "0", "n_forecast_pairs": "0", "n_primary_pairs": "0"})
    patients.append(patient_only)
    splits.extend([{"patient_id": "MU-scan-only", "fold": "1"},
                   {"patient_id": "MU-patient-only", "fold": "1"}])

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["scan_patients_missing_patient_manifest"] == 1
    assert codes["patient_manifest_rows_without_scans"] == 1


def test_nonincreasing_pair_times_and_scan_order_are_blocked():
    scans, patients, pairs, splits = valid_rows()
    patients[0]["eligible"] = "False"
    patients[0]["n_primary_pairs"] = "0"
    pairs[0].update({"tp_out": "1", "day_out": "-5", "dt_days": "-5", "in_primary": "False"})

    report = audit_mu_manifests(scans, patients, pairs, splits)
    codes = {item["code"]: item["count"] for item in report["issues"]}

    assert codes["forecast_pair_scan_order_errors"] == 1
    assert codes["nonpositive_pair_intervals"] == 1


def test_invalid_pair_booleans_are_reported_without_exposing_row_data():
    scans, patients, pairs, splits = valid_rows()
    pairs[0]["input_core_nonempty"] = "unknown"

    report = audit_mu_manifests(scans, patients, pairs, splits)

    assert report["status"] == "blocked"
    assert {"code": "invalid_forecast_pair_rows", "count": 1} in report["issues"]
    assert "MU-1" not in json.dumps(report)
