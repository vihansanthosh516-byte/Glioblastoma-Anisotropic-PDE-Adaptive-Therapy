from __future__ import annotations

import csv
import io
import json
import zipfile

import openpyxl
import pytest

from src.fk_duke_metadata_census import EXPECTED_MODALITIES, build_metadata_census


def _workbook(event_rows):
    stream = io.BytesIO()
    workbook = openpyxl.Workbook()
    sheet = workbook.active
    sheet.append([
        "Patient", "Event", "Notes", "Relative Start date (Days from Scan 1)",
        "Relative End date (Days from Scan 1)",
    ])
    for row in event_rows:
        sheet.append(row)
    workbook.save(stream)
    workbook.close()
    return stream.getvalue()


def _archive(path, visits, events, *, recurrence_progression_day=None):
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_STORED) as archive:
        csv_stream = io.StringIO(newline="")
        writer = csv.DictWriter(
            csv_stream,
            fieldnames=["PatientNumber", "ScanNumber", "ScanDate", "IDH", "Diagnosis"],
        )
        writer.writeheader()
        for scan, day in visits:
            writer.writerow({
                "PatientNumber": "Patient1", "ScanNumber": scan, "ScanDate": day,
                "IDH": "Negative", "Diagnosis": "Glioblastoma",
            })
        archive.writestr("DATASET_FINAL/patient_demographics_and_BTRADS.csv", csv_stream.getvalue())
        archive.writestr("DATASET_FINAL/patient_treatment_info.xlsx", _workbook(events))
        archive.writestr("DATASET_FINAL/Patient1/scan1/Patient1_scan1_landmarks_near.csv", "x,y,z\n1,2,3\n")

        sessions = [(1, "preoperative"), *[(scan, "regular") for scan, _ in visits if scan not in (1, 4)]]
        if any(scan == 4 for scan, _ in visits):
            sessions.append((4, "recurrence"))
        day_by_scan = dict(visits)
        for scan, suffix in sessions:
            folder = f"scan{scan}" if suffix == "regular" else f"scan{scan}_{suffix}"
            for modality in EXPECTED_MODALITIES:
                archive.writestr(
                    f"DATASET_FINAL/Patient1/{folder}/Patient1_scan{scan}_{modality}.nii.gz",
                    b"synthetic-member-body-must-not-be-opened",
                )
        assert set(day_by_scan) == {scan for scan, _ in visits}


def test_census_joins_metadata_without_opening_image_or_landmark_members(tmp_path, monkeypatch):
    path = tmp_path / "synthetic.zip"
    _archive(
        path,
        [(1, 0), (2, 10), (3, 90), (4, 180)],
        [
            ["Patient1", "Surgery", "ignored note", 1, 1],
            [None, None, None, None, None],
            ["Patient1", "Radiation + TMZ", "ignored note", 5, 80],
            ["Patient1", "Progression", "ignored note", 180, 180],
        ],
    )
    opened = []
    original_open = zipfile.ZipFile.open

    def track_open(self, name, *args, **kwargs):
        member_name = name.filename if isinstance(name, zipfile.ZipInfo) else str(name)
        opened.append(member_name)
        if member_name.lower().endswith((".nii.gz", "landmarks_near.csv")):
            raise AssertionError("census must not open NIfTI or landmark members")
        return original_open(self, name, *args, **kwargs)

    monkeypatch.setattr(zipfile.ZipFile, "open", track_open)
    result = build_metadata_census(path)

    assert result["metadata"]["scan_row_count"] == 4
    assert result["metadata"]["unique_patient_count"] == 1
    assert result["filename_join"]["metadata_pairs_missing_from_filenames"] == 0
    assert result["filename_join"]["filename_pairs_missing_from_metadata"] == 0
    assert result["filename_join"]["session_groups_with_wrong_modality_set"] == 0
    assert result["treatment_metadata"]["blank_rows_skipped"] == 1
    assert result["treatment_metadata"]["event_row_count"] == 3
    assert result["timeline_screen"]["patients_with_three_followup_scans"] == 1
    assert result["timeline_screen"]["first_triplet_v2_v3_gap_counts"]["76_to_105_days"] == 1
    assert result["timeline_screen"]["recurrence_vs_progression_date_counts"] == {
        "same_day_as_progression_event": 1
    }
    assert result["read_scope"]["patient_identifiers_in_output"] is False
    assert not any(name.lower().endswith((".nii.gz", "landmarks_near.csv")) for name in opened)
    assert "Patient1" not in json.dumps(result)


def test_census_stops_triplet_at_next_surgery(tmp_path):
    path = tmp_path / "repeat_surgery.zip"
    _archive(
        path,
        [(1, 0), (2, 10), (3, 80), (4, 150)],
        [
            ["Patient1", "Surgery", "", 1, 1],
            ["Patient1", "Surgery", "", 100, 100],
        ],
    )
    result = build_metadata_census(path)
    assert result["timeline_screen"]["patients_with_three_followup_scans"] == 0
    assert result["timeline_screen"]["first_triplet_v2_v3_gap_counts"]["76_to_105_days"] == 0


def test_census_requires_unique_metadata_members(tmp_path):
    path = tmp_path / "missing_metadata.zip"
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("a.csv", "x\n1\n")
    with pytest.raises(ValueError, match="exactly one archive member"):
        build_metadata_census(path)


def test_census_refuses_populated_treatment_row_without_patient(tmp_path):
    path = tmp_path / "missing_patient.zip"
    _archive(
        path,
        [(1, 0), (2, 10), (3, 90)],
        [[None, "Surgery", "", 1, 1]],
    )
    with pytest.raises(ValueError, match="invalid Patient"):
        build_metadata_census(path)
