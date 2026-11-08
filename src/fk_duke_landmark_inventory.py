"""Aggregate-only audit of Duke preoperative/recurrence vessel landmarks.

This reads landmark CSV coordinates and ZIP member names only. It never opens
MRI or segmentation members, and it returns no patient identifiers or point
coordinates. The audit checks file pairing and coordinate integrity; it does
not assess registration quality or select a study cohort.
"""
from __future__ import annotations

import csv
import io
import math
import re
import statistics
import zipfile
from collections import Counter, defaultdict
from pathlib import Path


_LANDMARK_PATH = re.compile(
    r"(?:^|/)Patient(?P<patient>\d+)/Patient(?P<file_patient>\d+)_scan(?P<scan>\d+)_"
    r"landmarks_(?P<zone>near|far)\.csv$",
    re.IGNORECASE,
)
_T1C_PATH = re.compile(
    r"(?:^|/)Patient(?P<patient>\d+)/scan(?P<scan>\d+)(?:_(?P<visit>preoperative|recurrence))?/"
    r"Patient(?P<file_patient>\d+)_scan(?P<file_scan>\d+)_t1c\.nii\.gz$",
    re.IGNORECASE,
)
_VISITS = ("preoperative", "recurrence")
_ZONES = ("near", "far")


def _summary(values):
    values = sorted(values)
    if not values:
        return {"files": 0, "min_points": None, "median_points": None, "max_points": None}
    return {
        "files": len(values),
        "min_points": int(values[0]),
        "median_points": float(statistics.median(values)),
        "max_points": int(values[-1]),
    }


def _read_landmark_member(archive, info, volume_shape):
    with archive.open(info, "r") as binary:
        with io.TextIOWrapper(binary, encoding="utf-8-sig", errors="replace", newline="") as text:
            reader = csv.DictReader(text)
            if reader.fieldnames is None:
                raise ValueError("landmark CSV has no header")
            columns = {name.strip().casefold(): name for name in reader.fieldnames if name}
            if not {"x", "y", "z"}.issubset(columns):
                raise ValueError("landmark CSV must contain x, y, and z columns")
            points = []
            invalid = 0
            fractional = 0
            out_of_bounds = 0
            for row in reader:
                raw = [row.get(columns[axis]) for axis in ("x", "y", "z")]
                try:
                    point = tuple(float(value) for value in raw)
                except (TypeError, ValueError):
                    invalid += 1
                    continue
                if not all(math.isfinite(value) for value in point):
                    invalid += 1
                    continue
                if any(value % 1 != 0 for value in point):
                    fractional += 1
                if any(value < 0 or value >= volume_shape[axis] for axis, value in enumerate(point)):
                    out_of_bounds += 1
                    continue
                points.append(point)
    duplicate_count = len(points) - len(set(points))
    return {
        "point_count": len(points),
        "invalid_coordinate_rows": invalid,
        "fractional_coordinate_rows": fractional,
        "out_of_bounds_rows": out_of_bounds,
        "duplicate_points": duplicate_count,
    }


def build_landmark_inventory(archive_path, *, volume_shape=(182, 218, 182)) -> dict:
    """Return aggregate checks for paired landmark CSVs; no images or masks are opened."""
    archive_path = Path(archive_path)
    if not archive_path.is_file():
        raise FileNotFoundError(f"archive does not exist: {archive_path}")
    if (len(volume_shape) != 3 or any(isinstance(value, bool) or not isinstance(value, int) or value < 1
                                      for value in volume_shape)):
        raise ValueError("volume_shape must contain three positive integers")

    try:
        with zipfile.ZipFile(archive_path, "r") as archive:
            infos = archive.infolist()
            scan_visit_tags = defaultdict(set)
            scan_filename_id_mismatches = 0
            for info in infos:
                if info.is_dir() or not info.filename.casefold().endswith("_t1c.nii.gz"):
                    continue
                match = _T1C_PATH.search(info.filename)
                if match is None:
                    continue
                patient, scan, visit, file_patient, file_scan = match.groups()
                if patient != file_patient or scan != file_scan:
                    scan_filename_id_mismatches += 1
                    continue
                key = (str(int(patient)), str(int(scan)))
                scan_visit_tags[key].add((visit or "regular").casefold())

            groups = defaultdict(dict)
            unparsed_landmark_files = 0
            duplicate_group_files = 0
            landmark_scan_without_t1c = 0
            ambiguous_visit_tag = 0
            landmark_bodies_opened = 0
            files_seen = 0
            issue_totals = Counter()
            issue_totals["scan_filename_id_mismatches"] = scan_filename_id_mismatches
            point_counts = defaultdict(list)

            for info in archive.infolist():
                if info.is_dir() or not info.filename.casefold().endswith(".csv"):
                    continue
                if "landmarks_" not in info.filename.casefold():
                    continue
                files_seen += 1
                match = _LANDMARK_PATH.search(info.filename)
                if match is None:
                    unparsed_landmark_files += 1
                    continue
                patient, file_patient, scan, zone = match.groups()
                if patient != file_patient:
                    issue_totals["filename_id_mismatches"] += 1
                    continue
                patient, scan = str(int(patient)), str(int(scan))
                visits = scan_visit_tags.get((patient, scan), set())
                if not visits:
                    landmark_scan_without_t1c += 1
                    continue
                special_visits = visits & {"preoperative", "recurrence"}
                if len(visits) != 1 or len(special_visits) != 1:
                    ambiguous_visit_tag += 1
                    continue
                visit = next(iter(special_visits))
                key = (patient, visit, zone.casefold())
                if key in groups:
                    duplicate_group_files += 1
                    continue
                record = _read_landmark_member(archive, info, volume_shape)
                landmark_bodies_opened += 1
                groups[key] = record
                point_counts[(visit.casefold(), zone.casefold())].append(record["point_count"])
                for field in ("invalid_coordinate_rows", "fractional_coordinate_rows", "out_of_bounds_rows",
                              "duplicate_points"):
                    issue_totals[field] += record[field]

            patients = {patient for patient, _visit, _zone in groups}
            complete = 0
            missing_groups = 0
            near_pair_count_mismatches = 0
            far_pair_count_mismatches = 0
            for patient in patients:
                required = {(visit, zone) for visit in _VISITS for zone in _ZONES}
                found = {(visit, zone) for p, visit, zone in groups if p == patient}
                missing_groups += len(required - found)
                complete += int(required == found)
                for zone in _ZONES:
                    before = groups.get((patient, "preoperative", zone))
                    after = groups.get((patient, "recurrence", zone))
                    if before is None or after is None:
                        continue
                    if before["point_count"] != after["point_count"]:
                        if zone == "near":
                            near_pair_count_mismatches += 1
                        else:
                            far_pair_count_mismatches += 1

            return {
                "archive_landmark_csv_count": files_seen,
                "unparsed_landmark_path_count": unparsed_landmark_files,
                "filename_id_mismatch_count": issue_totals["filename_id_mismatches"],
                "t1c_scan_filename_id_mismatch_count": issue_totals["scan_filename_id_mismatches"],
                "landmark_scans_without_t1c": landmark_scan_without_t1c,
                "ambiguous_landmark_visit_tags": ambiguous_visit_tag,
                "duplicate_patient_visit_zone_files": duplicate_group_files,
                "patient_count_with_any_parsed_landmarks": len(patients),
                "patients_with_all_four_visit_zone_files": complete,
                "missing_patient_visit_zone_groups": missing_groups,
                "near_pair_point_count_mismatches": near_pair_count_mismatches,
                "far_pair_point_count_mismatches": far_pair_count_mismatches,
                "point_count_summaries": {
                    f"{visit}_{zone}": _summary(point_counts[(visit, zone)])
                    for visit in _VISITS for zone in _ZONES
                },
                "coordinate_quality_totals": {
                    "invalid_coordinate_rows": issue_totals["invalid_coordinate_rows"],
                    "fractional_coordinate_rows": issue_totals["fractional_coordinate_rows"],
                    "out_of_bounds_rows": issue_totals["out_of_bounds_rows"],
                    "duplicate_points_within_file": issue_totals["duplicate_points"],
                },
                "read_scope": {
                    "landmark_coordinate_csv_bodies_opened": landmark_bodies_opened > 0,
                    "nifti_member_bodies_opened": False,
                    "segmentation_member_bodies_opened": False,
                    "patient_identifiers_or_coordinates_in_output": False,
                    "registration_quality_scored": False,
                    "cohort_selected_or_frozen": False,
                },
            }
    except (OSError, zipfile.BadZipFile, zipfile.LargeZipFile) as exc:
        raise ValueError(f"could not inspect Duke landmark archive: {archive_path.name}") from exc
