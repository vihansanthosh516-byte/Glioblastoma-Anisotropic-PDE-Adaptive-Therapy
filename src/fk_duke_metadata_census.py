"""Aggregate-only metadata census for the restricted Duke GBM archive.

Reads the demographics CSV and treatment workbook plus ZIP central-directory
filenames. It does not open NIfTI or landmark CSV members and returns no patient
identifiers. The triplet screen is for feasibility only; it does not freeze a
cohort or authorize outcome scoring.
"""
from __future__ import annotations

import csv
import io
import math
import re
import zipfile
from collections import Counter, defaultdict
from pathlib import Path

import openpyxl


DEMOGRAPHICS_BASENAME = "patient_demographics_and_BTRADS.csv"
TREATMENT_BASENAME = "patient_treatment_info.xlsx"
EXPECTED_MODALITIES = frozenset({"adc", "dwi", "flair", "t1", "t1c", "t2", "tumor-seg"})
_NIFTI_PATH = re.compile(
    r"(?:^|/)Patient(\d+)/scan(\d+)(?:_(preoperative|recurrence))?/"
    r"Patient(\d+)_scan(\d+)_(adc|dwi|flair|t1|t1c|t2|tumor-seg)\.nii\.gz$",
    re.IGNORECASE,
)


def _unique_basename(infos, basename: str):
    matches = [info for info in infos if info.filename.rsplit("/", 1)[-1].casefold() == basename.casefold()]
    if len(matches) != 1 or matches[0].is_dir():
        raise ValueError(f"expected exactly one archive member named {basename}")
    return matches[0]


def _numeric_id(value, field_name: str) -> str:
    match = re.search(r"\d+", str(value or ""))
    if match is None:
        raise ValueError(f"invalid {field_name} in metadata")
    return str(int(match.group(0)))


def _number_or_none(value):
    if value is None or str(value).strip() == "":
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def _read_demographics(archive, info):
    with archive.open(info, "r") as binary:
        with io.TextIOWrapper(binary, encoding="utf-8-sig", errors="replace", newline="") as text:
            reader = csv.DictReader(text)
            required = {"PatientNumber", "ScanNumber", "ScanDate", "IDH", "Diagnosis"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValueError("demographics CSV is missing required metadata columns")
            rows = {}
            for row in reader:
                patient = _numeric_id(row.get("PatientNumber"), "PatientNumber")
                scan = _numeric_id(row.get("ScanNumber"), "ScanNumber")
                day = _number_or_none(row.get("ScanDate"))
                key = (patient, scan)
                if key in rows:
                    raise ValueError("duplicate patient/scan row in demographics metadata")
                rows[key] = {
                    "patient": patient,
                    "scan": int(scan),
                    "day": day,
                    "idh": (row.get("IDH") or "").strip() or "[missing]",
                    "diagnosis": (row.get("Diagnosis") or "").strip() or "[missing]",
                }
    return rows


def _read_treatment_events(archive, info):
    events = defaultdict(list)
    event_patients = defaultdict(set)
    missing_start_dates = 0
    row_count = 0
    blank_row_count = 0
    with archive.open(info, "r") as binary:
        workbook = openpyxl.load_workbook(binary, read_only=True, data_only=True)
        try:
            sheet = workbook[workbook.sheetnames[0]]
            iterator = sheet.iter_rows(values_only=True)
            header = next(iterator, None)
            if header is None:
                raise ValueError("treatment workbook has no header")
            columns = {str(value).strip(): index for index, value in enumerate(header) if value is not None}
            required = {
                "Patient", "Event", "Relative Start date (Days from Scan 1)",
                "Relative End date (Days from Scan 1)",
            }
            if not required.issubset(columns):
                raise ValueError("treatment workbook is missing required metadata columns")
            start_index = columns["Relative Start date (Days from Scan 1)"]
            end_index = columns["Relative End date (Days from Scan 1)"]
            patient_index = columns["Patient"]
            event_index = columns["Event"]
            for row in iterator:
                if not any(value is not None and str(value).strip() for value in row):
                    blank_row_count += 1
                    continue
                row_count += 1
                patient = _numeric_id(row[patient_index] if patient_index < len(row) else None, "Patient")
                raw_event = row[event_index] if event_index < len(row) else None
                event = str(raw_event).strip() if raw_event is not None else "[missing]"
                start = _number_or_none(row[start_index] if start_index < len(row) else None)
                end = _number_or_none(row[end_index] if end_index < len(row) else None)
                if start is None:
                    missing_start_dates += 1
                else:
                    events[patient].append((event.casefold(), start, end))
                if event != "[missing]":
                    event_patients[event.casefold()].add(patient)
        finally:
            workbook.close()
    return events, event_patients, row_count, missing_start_dates, blank_row_count


def _read_filename_inventory(infos):
    groups = defaultdict(list)
    suffixes = defaultdict(set)
    unparsed_nifti = 0
    id_mismatch_nifti = 0
    nifti_count = 0
    csv_count = 0
    xlsx_count = 0
    for info in infos:
        if info.is_dir():
            continue
        lower = info.filename.casefold()
        if lower.endswith(".csv"):
            csv_count += 1
        if lower.endswith(".xlsx"):
            xlsx_count += 1
        if not lower.endswith(".nii.gz"):
            continue
        nifti_count += 1
        match = _NIFTI_PATH.search(info.filename)
        if match is None:
            unparsed_nifti += 1
            continue
        folder_patient, folder_scan, suffix, file_patient, file_scan, modality = match.groups()
        if folder_patient != file_patient or folder_scan != file_scan:
            id_mismatch_nifti += 1
            continue
        patient = str(int(folder_patient))
        scan = str(int(folder_scan))
        key = (patient, scan, (suffix or "regular").casefold())
        groups[key].append(modality.casefold())
        suffixes[(patient, scan)].add((suffix or "regular").casefold())

    bad_modality_groups = sum(Counter(modalities) != Counter(EXPECTED_MODALITIES) for modalities in groups.values())
    metadata_pairs = {(patient, scan) for patient, scan, _ in groups}
    return {
        "nifti_member_count": nifti_count,
        "csv_member_count": csv_count,
        "xlsx_member_count": xlsx_count,
        "session_folder_group_count": len(groups),
        "session_groups_with_wrong_modality_set": bad_modality_groups,
        "unparsed_nifti_path_count": unparsed_nifti,
        "patient_scan_filename_mismatch_count": id_mismatch_nifti,
        "filename_scan_pair_count": len(metadata_pairs),
        "suffix_session_counts": dict(sorted(Counter(suffix for _, _, suffix in groups).items())),
        "suffixes_by_patient_scan": suffixes,
        "groups": groups,
    }


def _summarize(
    rows, events, event_patients, event_row_count, missing_start_dates, blank_event_rows, filenames
):
    metadata_pairs = set(rows)
    filename_pairs = {(patient, scan) for patient, scan, _ in filenames["groups"]}
    patients = defaultdict(list)
    for row in rows.values():
        patients[row["patient"]].append(row)

    idh_by_patient = Counter()
    diagnosis_by_patient = Counter()
    for visits in patients.values():
        idh_by_patient[";".join(sorted({row["idh"] for row in visits}))] += 1
        diagnosis_by_patient[";".join(sorted({row["diagnosis"] for row in visits}))] += 1

    surgery_by_patient = defaultdict(list)
    progression_by_patient = defaultdict(list)
    for patient, items in events.items():
        for event, start, _end in items:
            if event == "surgery":
                surgery_by_patient[patient].append(start)
            elif event == "progression":
                progression_by_patient[patient].append(start)

    preop_sessions = []
    recurrence_sessions = []
    for (patient, scan), tags in filenames["suffixes_by_patient_scan"].items():
        day = rows.get((patient, scan), {}).get("day")
        if day is None:
            continue
        if "preoperative" in tags:
            preop_sessions.append((patient, scan, day))
        if "recurrence" in tags:
            recurrence_sessions.append((patient, scan, day))

    preop_surgery_relation = Counter()
    between_surgery_count = 0
    candidates = []
    candidate_gap_counts = Counter({"76_to_105_days": 0, "76_to_104_days": 0})
    candidate_idh = Counter()
    ambiguous_preop = 0
    no_surgery_after_preop = 0
    for patient, preop_scan, preop_day in preop_sessions:
        patient_preop = [item for item in preop_sessions if item[0] == patient]
        if len(patient_preop) != 1:
            ambiguous_preop += 1
            continue
        surgeries = sorted(set(surgery_by_patient.get(patient, [])))
        future_surgeries = [day for day in surgeries if day >= preop_day]
        if not future_surgeries:
            no_surgery_after_preop += 1
            continue
        index_surgery = future_surgeries[0]
        if index_surgery == preop_day:
            preop_surgery_relation["same_day"] += 1
        else:
            preop_surgery_relation["before_later_surgery"] += 1
        if any(day < preop_day for day in surgeries) and any(day > preop_day for day in surgeries):
            between_surgery_count += 1
        next_surgery = min((day for day in surgeries if day > index_surgery), default=math.inf)
        followup = sorted(
            (row for row in patients[patient] if row["day"] is not None
             and row["day"] > index_surgery and row["day"] < next_surgery),
            key=lambda row: (row["day"], row["scan"]),
        )
        if len(followup) >= 3:
            first_three = followup[:3]
            gap = first_three[2]["day"] - first_three[1]["day"]
            candidates.append((patient, gap))
            candidate_gap_counts["76_to_105_days"] += 76 <= gap <= 105
            candidate_gap_counts["76_to_104_days"] += 76 <= gap <= 104
            candidate_idh[";".join(sorted({row["idh"] for row in patients[patient]}))] += 1

    recurrence_vs_progression = Counter()
    for patient, _scan, day in recurrence_sessions:
        dates = sorted(set(progression_by_patient.get(patient, [])))
        if not dates:
            relation = "no_dated_progression_event"
        elif day in dates:
            relation = "same_day_as_progression_event"
        elif any(value < day for value in dates) and any(value > day for value in dates):
            relation = "between_progression_events"
        elif all(day < value for value in dates):
            relation = "before_all_progression_events"
        else:
            relation = "after_all_progression_events"
        recurrence_vs_progression[relation] += 1

    selected_event_names = ("surgery", "radiation + tmz", "progression")
    selected_event_counts = {
        name: len(event_patients.get(name, set())) for name in selected_event_names
    }
    return {
        "metadata": {
            "scan_row_count": len(rows),
            "unique_patient_count": len(patients),
            "rows_with_missing_or_invalid_scan_date": sum(row["day"] is None for row in rows.values()),
            "duplicate_scan_rows": 0,
            "idh_status_patient_counts": dict(sorted(idh_by_patient.items())),
            "diagnosis_status_patient_counts": dict(sorted(diagnosis_by_patient.items())),
        },
        "filename_join": {
            "metadata_patient_scan_pair_count": len(metadata_pairs),
            "filename_patient_scan_pair_count": len(filename_pairs),
            "metadata_pairs_missing_from_filenames": len(metadata_pairs - filename_pairs),
            "filename_pairs_missing_from_metadata": len(filename_pairs - metadata_pairs),
            "nifti_member_count": filenames["nifti_member_count"],
            "csv_member_count": filenames["csv_member_count"],
            "xlsx_member_count": filenames["xlsx_member_count"],
            "session_folder_group_count": filenames["session_folder_group_count"],
            "session_groups_with_wrong_modality_set": filenames["session_groups_with_wrong_modality_set"],
            "unparsed_nifti_path_count": filenames["unparsed_nifti_path_count"],
            "patient_scan_filename_mismatch_count": filenames["patient_scan_filename_mismatch_count"],
            "suffix_session_counts": filenames["suffix_session_counts"],
        },
        "treatment_metadata": {
            "event_row_count": event_row_count,
            "blank_rows_skipped": blank_event_rows,
            "rows_with_missing_or_invalid_start_date": missing_start_dates,
            "selected_event_patient_counts": selected_event_counts,
        },
        "timeline_screen": {
            "status": "planning-only; not a frozen cohort",
            "preoperative_sessions": len(preop_sessions),
            "preoperative_vs_surgery_date_counts": dict(sorted(preop_surgery_relation.items())),
            "preoperative_between_surgery_dates": between_surgery_count,
            "preoperative_without_dated_future_surgery": no_surgery_after_preop,
            "ambiguous_preoperative_session_count": ambiguous_preop,
            "recurrence_sessions": len(recurrence_sessions),
            "recurrence_vs_progression_date_counts": dict(sorted(recurrence_vs_progression.items())),
            "candidate_rule": (
                "For each unique preoperative-tagged scan, choose the first dated surgery on or after it; "
                "take the first three later scan sessions before the next surgery."
            ),
            "patients_with_three_followup_scans": len(candidates),
            "first_triplet_v2_v3_gap_counts": dict(sorted(candidate_gap_counts.items())),
            "candidate_idh_status_counts": dict(sorted(candidate_idh.items())),
        },
        "read_scope": {
            "metadata_members_opened": [DEMOGRAPHICS_BASENAME, TREATMENT_BASENAME],
            "nifti_member_bodies_opened": False,
            "landmark_csv_member_bodies_opened": False,
            "patient_identifiers_in_output": False,
            "voxel_arrays_read": False,
        },
    }


def build_metadata_census(archive_path) -> dict:
    """Return aggregate metadata and filename checks; never open image or landmark members."""
    archive_path = Path(archive_path)
    if not archive_path.is_file():
        raise FileNotFoundError(f"archive does not exist: {archive_path}")
    try:
        with zipfile.ZipFile(archive_path, "r") as archive:
            infos = archive.infolist()
            demographics_info = _unique_basename(infos, DEMOGRAPHICS_BASENAME)
            treatment_info = _unique_basename(infos, TREATMENT_BASENAME)
            rows = _read_demographics(archive, demographics_info)
            events, event_patients, event_rows, missing_start_dates, blank_event_rows = (
                _read_treatment_events(archive, treatment_info)
            )
            filenames = _read_filename_inventory(infos)
    except (OSError, zipfile.BadZipFile, zipfile.LargeZipFile) as exc:
        raise ValueError(f"could not inspect Duke metadata archive: {archive_path.name}") from exc
    return _summarize(
        rows, events, event_patients, event_rows, missing_start_dates, blank_event_rows, filenames
    )
