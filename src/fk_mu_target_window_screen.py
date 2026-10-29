"""Outcome-blind metadata and header screen for the MU Track B development cohort.

This screen selects by diagnosis, visit dates, and file availability only. It reads NIfTI headers for selected
files but never accesses target-mask voxel arrays or outcomes. Header alignment is not registration validation.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SCAN_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "scan_manifest.csv"
DEFAULT_PATIENT_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "patient_manifest.csv"
DEFAULT_DATA_ROOT = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "fk_mu_target_window_screen_20261028"
GRID_TOLERANCE_MM = 1e-4
TARGET_GAP_MIN_DAYS = 76.0
TARGET_GAP_MAX_DAYS = 104.0


def _read_csv(path: Path, fields: set[str]) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"manifest is empty: {path}") from None
        indices = {name: index for index, name in enumerate(header)}
        if not fields.issubset(indices):
            raise ValueError(f"manifest is missing required columns: {sorted(fields - set(indices))}")
        rows = []
        for raw in reader:
            if len(raw) < len(header):
                raise ValueError(f"manifest has a malformed row: {path}")
            rows.append({name: raw[indices[name]] for name in fields})
    if not rows:
        raise ValueError(f"manifest is empty: {path}")
    return rows


def _truth(value: str) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes"}


def _projected_rows_sha256(rows: list[dict[str, str]]) -> str:
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _header(path: Path, data_root: Path) -> dict:
    try:
        path.resolve().relative_to(data_root.resolve())
    except ValueError as exc:
        raise ValueError("selected file path escapes the MU data root") from exc
    if not path.is_file():
        raise FileNotFoundError(path)
    import nibabel as nib

    image = nib.load(str(path), mmap=False)
    shape = tuple(int(value) for value in image.shape)
    affine = np.asarray(image.affine, dtype=float)
    spacing = tuple(float(value) for value in image.header.get_zooms()[:3])
    xyz_unit = image.header.get_xyzt_units()[0]
    if len(shape) != 3 or affine.shape != (4, 4) or not np.isfinite(affine).all():
        raise ValueError("selected NIfTI has an invalid 3D header or affine")
    determinant = float(np.linalg.det(affine[:3, :3]))
    if not np.isfinite(determinant) or abs(determinant) < 1e-12:
        raise ValueError("selected NIfTI has a singular spatial affine")
    if len(spacing) != 3 or not np.isfinite(spacing).all() or any(value <= 0 for value in spacing):
        raise ValueError("selected NIfTI has invalid voxel spacing")
    if xyz_unit != "mm":
        raise ValueError(f"selected NIfTI spatial unit is not explicitly millimeters: {xyz_unit}")
    return {"shape": shape, "spacing_mm": spacing, "affine": affine, "xyz_unit": xyz_unit}


def _read_selected_history_core_counts(path: Path, selected_pairs: set[tuple[str, int]]) -> dict:
    """Read manifest core counts only for selected first/second history visits."""
    if not selected_pairs:
        return {}
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"manifest is empty: {path}") from None
        required = ("patient_id", "dataset", "scan_number", "core_voxels_native")
        indices = {name: index for index, name in enumerate(header)}
        if not set(required).issubset(indices):
            raise ValueError("scan manifest lacks history-core metadata columns")
        result = {}
        for raw in reader:
            if len(raw) < len(header):
                raise ValueError(f"manifest has a malformed row: {path}")
            if raw[indices["dataset"]] != "MU":
                continue
            try:
                key = (raw[indices["patient_id"]], int(raw[indices["scan_number"]]))
            except (TypeError, ValueError):
                continue
            if key in selected_pairs:
                result[key] = raw[indices["core_voxels_native"]]
    return result


def _selected_paths(patient_id: str, scan_numbers: list[int], data_root: Path) -> tuple[list[Path], Path]:
    patient_dir = data_root / patient_id
    mask_paths = [
        patient_dir / f"Timepoint_{number}" / f"{patient_id}_Timepoint_{number}_tumorMask.nii.gz"
        for number in scan_numbers
    ]
    origin_number = scan_numbers[1]
    origin_brain = (
        patient_dir / f"Timepoint_{origin_number}" /
        f"{patient_id}_Timepoint_{origin_number}_brain_t1n.nii.gz"
    )
    return mask_paths, origin_brain


def build_mu_window_screen(
    scan_manifest_path: Path = DEFAULT_SCAN_MANIFEST,
    patient_manifest_path: Path = DEFAULT_PATIENT_MANIFEST,
    data_root: Path = DEFAULT_DATA_ROOT,
) -> dict:
    """Build an aggregate-only first-two-to-third eligibility and header report."""
    scan_manifest_path = Path(scan_manifest_path).resolve()
    patient_manifest_path = Path(patient_manifest_path).resolve()
    data_root = Path(data_root).resolve()
    required_scan = {"patient_id", "dataset", "scan_number", "day", "mask_exists"}
    required_patient = {"patient_id", "dataset", "is_gbm"}
    scan_rows = _read_csv(scan_manifest_path, required_scan)
    patient_rows = _read_csv(patient_manifest_path, required_patient)

    gbm_rows = [row for row in patient_rows if row["dataset"] == "MU" and _truth(row["is_gbm"])]
    gbm_by_id = {row["patient_id"]: row for row in gbm_rows}
    if len(gbm_by_id) != len(gbm_rows):
        raise ValueError("MU patient manifest contains duplicate GBM patient rows")

    visits_by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in scan_rows:
        if row["dataset"] == "MU" and _truth(row["mask_exists"]) and row["patient_id"] in gbm_by_id:
            visits_by_patient[row["patient_id"]].append(row)

    counts: Counter[str] = Counter()
    header_shapes: Counter[str] = Counter()
    header_spacings: Counter[str] = Counter()
    header_files_read = 0
    header_pass_candidates = 0
    eligible_candidates = []
    history_pairs: set[tuple[str, int]] = set()

    for patient_id, rows in visits_by_patient.items():
        try:
            ordered = sorted(rows, key=lambda row: int(row["scan_number"]))
        except (TypeError, ValueError) as exc:
            counts["invalid_scan_number_rows"] += 1
            continue
        if len(ordered) < 3:
            continue
        counts["gbm_patients_with_three_mask_available_visits"] += 1
        first_three = ordered[:3]
        try:
            scan_numbers = [int(row["scan_number"]) for row in first_three]
            days = [float(row["day"]) for row in first_three]
        except (TypeError, ValueError):
            counts["first_three_dates_or_scan_numbers_invalid"] += 1
            continue
        if (len(set(scan_numbers)) != 3 or not np.isfinite(days).all()
                or not (days[0] < days[1] < days[2])):
            counts["first_three_dates_or_scan_numbers_invalid"] += 1
            continue
        counts["gbm_with_valid_first_three_dates"] += 1
        target_gap = days[2] - days[1]
        if not TARGET_GAP_MIN_DAYS <= target_gap <= TARGET_GAP_MAX_DAYS:
            continue
        counts["gbm_in_76_104_day_target_window"] += 1

        mask_paths, origin_brain_path = _selected_paths(patient_id, scan_numbers, data_root)
        if not all(path.is_file() for path in mask_paths) or not origin_brain_path.is_file():
            counts["target_window_candidates_missing_selected_files"] += 1
            continue
        counts["gbm_target_window_mask_and_origin_brain_files_present"] += 1

        headers = []
        try:
            for path in mask_paths:
                headers.append(_header(path, data_root))
                header_files_read += 1
            origin_brain_header = _header(origin_brain_path, data_root)
            header_files_read += 1
        except Exception as exc:
            counts[f"target_window_header_error_{type(exc).__name__}"] += 1
            continue
        reference = headers[0]
        aligned = all(
            header["shape"] == reference["shape"]
            and np.allclose(header["affine"], reference["affine"], rtol=0.0, atol=GRID_TOLERANCE_MM)
            and np.allclose(header["spacing_mm"], reference["spacing_mm"], rtol=0.0, atol=GRID_TOLERANCE_MM)
            for header in headers[1:] + [origin_brain_header]
        )
        if not aligned:
            counts["target_window_header_grid_mismatch"] += 1
            continue
        header_pass_candidates += 1
        shape_key = "x".join(str(value) for value in reference["shape"])
        spacing_key = "x".join(f"{value:.6g}" for value in reference["spacing_mm"])
        header_shapes[shape_key] += 1
        header_spacings[spacing_key] += 1
        history_pairs.update((patient_id, number) for number in scan_numbers[:2])
        eligible_candidates.append({
            "patient_id": patient_id,
            "history_scan_numbers": scan_numbers[:2],
            "days": days,
            "gap_days": target_gap,
        })

    history_core_counts = _read_selected_history_core_counts(scan_manifest_path, history_pairs)
    empty_or_unknown_history_patients = 0
    for candidate in eligible_candidates:
        history_values = [
            history_core_counts.get((candidate["patient_id"], scan_number), "")
            for scan_number in candidate["history_scan_numbers"]
        ]
        empty_or_unknown = False
        for raw_value in history_values:
            try:
                count = float(raw_value)
                empty_or_unknown |= not np.isfinite(count) or count <= 0
            except (TypeError, ValueError):
                empty_or_unknown = True
        empty_or_unknown_history_patients += int(empty_or_unknown)
    counts["header_pass_candidates_with_empty_or_unknown_history_core_in_manifest"] = (
        empty_or_unknown_history_patients
    )
    history_core_payload = json.dumps(
        sorted((patient_id, scan_number, history_core_counts.get((patient_id, scan_number), ""))
               for patient_id, scan_number in history_pairs),
        separators=(",", ":"), ensure_ascii=False,
    )

    return {
        "schema_version": 1,
        "study": "MU exploratory first-two-to-third 90-day-window development screen",
        "screen_rule": {
            "diagnosis": "MU rows with is_gbm=true",
            "visits": "first three mask-available visits ordered by scan number",
            "date_rule": "finite strictly increasing day values; no duplicate scan numbers",
            "target_window_days": [TARGET_GAP_MIN_DAYS, TARGET_GAP_MAX_DAYS],
            "selection_uses_target_mask_values": False,
        },
        "counts": dict(sorted(counts.items())),
        "header_check": {
            "candidate_count_with_mask_and_origin_brain_files_present": int(
                counts["gbm_target_window_mask_and_origin_brain_files_present"]
            ),
            "header_pass_candidates": int(header_pass_candidates),
            "nifti_headers_read": int(header_files_read),
            "shape_counts": dict(sorted(header_shapes.items())),
            "voxel_spacing_mm_counts": dict(sorted(header_spacings.items())),
            "affine_tolerance_mm": GRID_TOLERANCE_MM,
            "spatial_units_required": "mm",
            "voxel_arrays_read": False,
        },
        "selected_patient_count_after_header_check": int(len(eligible_candidates)),
        "access": {
            "manifest_columns_read": {
                "scan_manifest_for_selection": [
                    "patient_id", "dataset", "scan_number", "day", "mask_exists"
                ],
                "patient_manifest_for_selection": ["patient_id", "dataset", "is_gbm"],
                "history_only_summary_after_selection": [
                    "patient_id", "dataset", "scan_number", "core_voxels_native"
                ],
            },
            "nifti_headers_read": True,
            "target_mask_voxel_values_read": False,
            "MRI_intensities_read": False,
        },
        "source_hashes": {
            "scan_selection_columns_sha256": _projected_rows_sha256(scan_rows),
            "patient_selection_columns_sha256": _projected_rows_sha256(patient_rows),
            "selected_history_core_counts_sha256": hashlib.sha256(
                history_core_payload.encode("utf-8")
            ).hexdigest(),
            "screen_script_sha256": hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        },
        "interpretation_limit": (
            "MU is previously explored development data. The file-availability count only covers selected masks "
            "and the origin brain mask; it does not establish target-MRI availability or registration quality. "
            "Header agreement does not prove anatomical registration, valid annotation support, or external "
            "performance. This screen does not score any patient."
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-manifest", type=Path, default=DEFAULT_SCAN_MANIFEST)
    parser.add_argument("--patient-manifest", type=Path, default=DEFAULT_PATIENT_MANIFEST)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    output_dir = args.out_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing screen output: {output_dir}")
    report = build_mu_window_screen(args.scan_manifest, args.patient_manifest, args.data_root)
    output_dir.mkdir(parents=True)
    report_path = output_dir / "metadata_screen.json"
    report_path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "report": str(report_path), **report["counts"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
