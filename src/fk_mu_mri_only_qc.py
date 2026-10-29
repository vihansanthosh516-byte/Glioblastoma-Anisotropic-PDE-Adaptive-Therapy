"""MRI-only geometry and intensity screen for the MU Track B first-to-third cohort.

This screen reads skull-stripped T1 non-contrast MRI only. It never opens tumor-mask arrays. The reported image
support and intensity correlations are descriptive checks for gross longitudinal mismatch, not registration
certification or a basis for selecting target masks.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.fk_mu_target_window_screen import (
    DEFAULT_DATA_ROOT,
    DEFAULT_PATIENT_MANIFEST,
    DEFAULT_SCAN_MANIFEST,
    TARGET_GAP_MAX_DAYS,
    TARGET_GAP_MIN_DAYS,
    _read_csv,
    _truth,
    build_mu_window_screen,
)


DEFAULT_METADATA_REPORT = PROJECT_ROOT / "output" / "fk_mu_target_window_screen_20261028" / "metadata_screen.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "fk_mu_mri_only_qc_20261028"
SEQUENCES = ("t1n", "t1c", "t2w", "t2f")
GRID_TOLERANCE_MM = 1e-4


def _candidate_visits(scan_manifest: Path, patient_manifest: Path, data_root: Path) -> list[dict]:
    scans = _read_csv(scan_manifest, {"patient_id", "dataset", "scan_number", "day", "mask_exists"})
    patients = _read_csv(patient_manifest, {"patient_id", "dataset", "is_gbm"})
    gbm_ids = {row["patient_id"] for row in patients
               if row["dataset"] == "MU" and _truth(row["is_gbm"])}
    by_patient: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in scans:
        if row["dataset"] == "MU" and _truth(row["mask_exists"]) and row["patient_id"] in gbm_ids:
            by_patient[row["patient_id"]].append(row)

    selected = []
    for patient_id, rows in by_patient.items():
        try:
            first_three = sorted(rows, key=lambda row: int(row["scan_number"]))[:3]
            if len(first_three) != 3:
                continue
            scan_numbers = [int(row["scan_number"]) for row in first_three]
            days = [float(row["day"]) for row in first_three]
        except (TypeError, ValueError):
            continue
        if (len(set(scan_numbers)) != 3 or not np.isfinite(days).all()
                or not (days[0] < days[1] < days[2])
                or not TARGET_GAP_MIN_DAYS <= days[2] - days[1] <= TARGET_GAP_MAX_DAYS):
            continue
        patient_dir = data_root / patient_id
        mask_paths = [patient_dir / f"Timepoint_{number}" / f"{patient_id}_Timepoint_{number}_tumorMask.nii.gz"
                      for number in scan_numbers]
        mri_paths = {
            sequence: [patient_dir / f"Timepoint_{number}" /
                       f"{patient_id}_Timepoint_{number}_brain_{sequence}.nii.gz"
                       for number in scan_numbers]
            for sequence in SEQUENCES
        }
        if not all(path.is_file() for path in mask_paths) or not mri_paths["t1n"][1].is_file():
            continue
        selected.append({
            "patient_id": patient_id,
            "scan_numbers": scan_numbers,
            "days": days,
            "mri_paths": mri_paths,
        })
    return sorted(selected, key=lambda row: row["patient_id"])


def _header(path: Path) -> dict:
    import nibabel as nib

    image = nib.load(str(path), mmap=False)
    affine = np.asarray(image.affine, dtype=float)
    unit = image.header.get_xyzt_units()[0]
    if (len(image.shape) != 3 or not np.isfinite(affine).all() or unit != "mm"
            or not np.isfinite(image.header.get_zooms()[:3]).all()
            or np.any(np.asarray(image.header.get_zooms()[:3], dtype=float) <= 0)):
        raise ValueError(f"invalid 3D MRI geometry or non-mm unit: {path.name}")
    determinant = float(np.linalg.det(affine[:3, :3]))
    if not np.isfinite(determinant) or abs(determinant) < 1e-12:
        raise ValueError(f"singular MRI affine: {path.name}")
    return {
        "shape": tuple(int(value) for value in image.shape),
        "affine": affine,
        "spacing": tuple(float(value) for value in image.header.get_zooms()[:3]),
    }


def _pair_metrics(first: np.ndarray, second: np.ndarray) -> dict:
    support_first = np.isfinite(first) & (first > 0)
    support_second = np.isfinite(second) & (second > 0)
    intersection = support_first & support_second
    denominator = int(support_first.sum() + support_second.sum())
    support_dice = 2.0 * int(intersection.sum()) / denominator if denominator else float("nan")
    if int(intersection.sum()) < 100:
        correlation = float("nan")
    else:
        x = first[intersection].astype(np.float64, copy=False)
        y = second[intersection].astype(np.float64, copy=False)
        x -= x.mean()
        y -= y.mean()
        denominator_corr = float(np.sqrt(x @ x) * np.sqrt(y @ y))
        correlation = float((x @ y) / denominator_corr) if denominator_corr > 0 else float("nan")
    return {
        "nonzero_support_dice": support_dice,
        "t1n_intensity_correlation_on_common_nonzero_support": correlation,
        "common_nonzero_voxel_count": int(intersection.sum()),
    }


def _summarize(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    finite = array[np.isfinite(array)]
    if not len(finite):
        return {"n_finite": 0, "minimum": None, "p05": None, "median": None, "p95": None, "maximum": None}
    return {
        "n_finite": int(len(finite)),
        "minimum": float(finite.min()),
        "p05": float(np.quantile(finite, 0.05)),
        "median": float(np.median(finite)),
        "p95": float(np.quantile(finite, 0.95)),
        "maximum": float(finite.max()),
    }


def build_mri_only_qc(
    scan_manifest: Path = DEFAULT_SCAN_MANIFEST,
    patient_manifest: Path = DEFAULT_PATIENT_MANIFEST,
    data_root: Path = DEFAULT_DATA_ROOT,
    metadata_report: Path = DEFAULT_METADATA_REPORT,
) -> dict:
    scan_manifest = Path(scan_manifest).resolve()
    patient_manifest = Path(patient_manifest).resolve()
    data_root = Path(data_root).resolve()
    metadata_report = Path(metadata_report).resolve()
    prior = json.loads(metadata_report.read_text(encoding="utf-8"))
    if prior.get("selected_patient_count_after_header_check") != 18:
        raise ValueError("the frozen MU metadata screen does not report the expected 18 candidates")
    current_metadata_screen = build_mu_window_screen(scan_manifest, patient_manifest, data_root)
    if current_metadata_screen["source_hashes"] != prior.get("source_hashes"):
        raise ValueError("the frozen metadata report no longer matches the current manifests and screen code")
    candidates = _candidate_visits(scan_manifest, patient_manifest, data_root)
    if len(candidates) != prior["selected_patient_count_after_header_check"]:
        raise ValueError("MRI-only candidates do not reproduce the frozen metadata-screen count")

    pair_dice = []
    pair_correlation = []
    pair_count = 0
    mri_volumes_read = 0
    header_checks = 0
    sequence_files = 0
    all_t1n_headers_aligned = True
    for candidate in candidates:
        t1n_paths = candidate["mri_paths"]["t1n"]
        for paths in candidate["mri_paths"].values():
            sequence_files += sum(path.is_file() for path in paths)
        if not all(path.is_file() for path in t1n_paths):
            all_t1n_headers_aligned = False
            continue
        try:
            headers = [_header(path) for path in t1n_paths]
        except (OSError, ValueError):
            all_t1n_headers_aligned = False
            continue
        header_checks += len(headers)
        reference = headers[0]
        if not all(
            item["shape"] == reference["shape"]
            and np.allclose(item["affine"], reference["affine"], rtol=0.0, atol=GRID_TOLERANCE_MM)
            and np.allclose(item["spacing"], reference["spacing"], rtol=0.0, atol=GRID_TOLERANCE_MM)
            for item in headers[1:]
        ):
            all_t1n_headers_aligned = False
            continue
        import nibabel as nib

        image_arrays = []
        for path in t1n_paths:
            image = nib.load(str(path), mmap=False)
            image_arrays.append(np.asarray(image.dataobj, dtype=np.float32))
            mri_volumes_read += 1
        for earlier, later in ((0, 1), (1, 2)):
            metrics = _pair_metrics(image_arrays[earlier], image_arrays[later])
            pair_dice.append(metrics["nonzero_support_dice"])
            pair_correlation.append(metrics["t1n_intensity_correlation_on_common_nonzero_support"])
            pair_count += 1
        del image_arrays

    if not all_t1n_headers_aligned:
        geometry_status = "not_all_selected_t1n_headers_align"
    else:
        geometry_status = "all_selected_t1n_headers_align"
    return {
        "schema_version": 1,
        "study": "MU MRI-only longitudinal geometry/support screen; exploratory and descriptive",
        "selection": {
            "rule": "exactly the frozen first-three-mask-available GBM visits with a 76–104-day second-to-third gap",
            "selected_patient_count": len(candidates),
            "target_mask_values_used_for_selection": False,
        },
        "file_availability": {
            "all_four_skull_stripped_sequences_t1n_t1c_t2w_t2f_at_each_selected_visit": bool(
                sequence_files == len(candidates) * 3 * len(SEQUENCES)
            ),
            "sequence_files_present": int(sequence_files),
            "sequence_files_expected": int(len(candidates) * 3 * len(SEQUENCES)),
        },
        "t1n_geometry": {
            "header_comparisons": int(header_checks),
            "all_first_three_visit_headers_align_within_mm_tolerance": bool(all_t1n_headers_aligned),
            "affine_tolerance_mm": GRID_TOLERANCE_MM,
            "status": geometry_status,
        },
        "pairwise_mri_only_metrics": {
            "n_adjacent_visit_pairs": int(pair_count),
            "nonzero_support_dice": _summarize(pair_dice),
            "t1n_intensity_correlation_on_common_nonzero_support": _summarize(pair_correlation),
            "support_definition": "finite, strictly positive voxels in the skull-stripped brain_t1n image; a support proxy, not a validated acquisition mask",
        },
        "data_access": {
            "t1n_mri_volumes_read": int(mri_volumes_read),
            "tumor_mask_files_opened": 0,
            "tumor_mask_voxel_arrays_read": 0,
            "target_mri_used_to_fit_or_accept_a_registration": False,
        },
        "source_hashes": {
            "metadata_screen_sha256": hashlib.sha256(metadata_report.read_bytes()).hexdigest(),
            "mri_qc_script_sha256": hashlib.sha256(Path(__file__).resolve().read_bytes()).hexdigest(),
        },
        "gate_status": "descriptive_screen_only; no registration or valid-support gate cleared",
        "interpretation_limit": (
            "The MU paper reports atlas registration and manual image review but does not provide per-case transforms "
            "or registration metrics. Header agreement, nonzero-image-support overlap, and same-sequence correlation "
            "can flag gross mismatches but do not certify anatomical or tumor-boundary alignment, nor establish that "
            "zero mask labels mean valid negative coverage. Do not score patient masks from this screen alone."
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan-manifest", type=Path, default=DEFAULT_SCAN_MANIFEST)
    parser.add_argument("--patient-manifest", type=Path, default=DEFAULT_PATIENT_MANIFEST)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--metadata-report", type=Path, default=DEFAULT_METADATA_REPORT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    output_dir = args.out_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite MRI-only QC output: {output_dir}")
    report = build_mri_only_qc(args.scan_manifest, args.patient_manifest, args.data_root, args.metadata_report)
    output_dir.mkdir(parents=True)
    path = output_dir / "mri_only_qc.json"
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "complete",
        "report": str(path),
        "selected_patient_count": report["selection"]["selected_patient_count"],
        "t1n_mri_volumes_read": report["data_access"]["t1n_mri_volumes_read"],
        "tumor_mask_voxel_arrays_read": report["data_access"]["tumor_mask_voxel_arrays_read"],
        "gate_status": report["gate_status"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
