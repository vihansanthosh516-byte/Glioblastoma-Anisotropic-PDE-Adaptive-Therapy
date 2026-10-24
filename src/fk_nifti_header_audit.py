"""Aggregate-only NIfTI header audit for Track B; never reads image voxel arrays."""
from __future__ import annotations

import argparse
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np

from src.fk_manifest_audit import _local_mri_file_kind


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_IMAGE_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
_MRI_KINDS = ("t1n", "t1c", "t2w", "t2f")
_GRID_TOLERANCE_MM = 1e-4


def audit_nifti_headers(image_dir: Path = DEFAULT_IMAGE_DIR) -> dict:
    """Check shape, voxel spacing, and affine agreement using NIfTI headers only.

    This function calls ``nib.load`` and reads ``shape``, ``affine``, and header zooms. It never accesses
    ``dataobj``, ``get_fdata``, or voxel values. The aggregate report contains no patient or file names.
    """
    try:
        import nibabel as nib
    except ImportError as exc:
        raise RuntimeError("nibabel is required for the NIfTI header audit") from exc

    image_dir = Path(image_dir)
    if not image_dir.is_dir():
        return {
            "status": "blocked",
            "summary": {"image_directory_present": False},
            "issues": [{"code": "local_image_directory_missing", "count": 1}],
            "interpretation_limit": "NIfTI voxel arrays and pixel values were not read",
        }

    issues: Counter[str] = Counter()
    header_counts: Counter[str] = Counter()
    t1c_shapes: Counter[str] = Counter()
    t1c_spacings: Counter[str] = Counter()
    patient_dirs = [path for path in image_dir.iterdir() if path.is_dir()]
    timepoint_dirs = []
    for patient_dir in patient_dirs:
        timepoint_dirs.extend(
            path for path in patient_dir.rglob("Timepoint_*")
            if path.is_dir() and re.fullmatch(r"Timepoint_\d+", path.name)
        )

    total_files = 0
    matched_timepoints = 0
    for timepoint_dir in timepoint_dirs:
        files = [path for path in timepoint_dir.iterdir()
                 if path.is_file() and path.name.lower().endswith(".nii.gz")]
        total_files += len(files)
        paths_by_kind: dict[str, list[Path]] = {}
        for path in files:
            kind = _local_mri_file_kind(path.name)
            header_counts[kind] += 1
            paths_by_kind.setdefault(kind, []).append(path)
        if any(len(paths_by_kind.get(kind, ())) != 1 for kind in _MRI_KINDS):
            issues["timepoints_missing_or_duplicating_mri_sequences"] += 1
        if len(paths_by_kind.get("segmentation", ())) > 1:
            issues["timepoints_with_duplicate_segmentations"] += 1
        if paths_by_kind.get("unknown"):
            issues["unclassified_nifti_files"] += len(paths_by_kind["unknown"])

        loaded = {}
        for kind, paths in paths_by_kind.items():
            for index, path in enumerate(paths):
                try:
                    image = nib.load(str(path), mmap=False)
                    shape = tuple(int(value) for value in image.shape)
                    affine = np.asarray(image.affine, dtype=float)
                    spacing = tuple(float(value) for value in image.header.get_zooms()[:3])
                except Exception:
                    issues["unreadable_nifti_headers"] += 1
                    continue
                if len(shape) != 3:
                    issues["non_3d_nifti_headers"] += 1
                if (len(spacing) != 3 or not np.isfinite(spacing).all()
                        or any(value <= 0 for value in spacing)):
                    issues["invalid_voxel_spacing_headers"] += 1
                if affine.shape != (4, 4) or not np.isfinite(affine).all():
                    issues["invalid_affine_headers"] += 1
                if kind == "t1c":
                    t1c_shapes["x".join(map(str, shape))] += 1
                    t1c_spacings["x".join(f"{value:.6g}" for value in spacing)] += 1
                loaded.setdefault(kind, []).append((shape, spacing, affine))

        reference_kind = "t1c" if loaded.get("t1c") else next(
            (kind for kind in _MRI_KINDS if loaded.get(kind)), None
        )
        if reference_kind is None:
            issues["timepoints_without_readable_mri_header"] += 1
            continue
        reference = loaded[reference_kind][0]
        timepoint_has_grid_mismatch = False
        for kind, entries in loaded.items():
            for shape, spacing, affine in entries:
                if shape != reference[0]:
                    issues["within_timepoint_shape_mismatches"] += 1
                    timepoint_has_grid_mismatch = True
                if (len(spacing) == 3 and len(reference[1]) == 3
                        and not np.allclose(spacing, reference[1], rtol=0.0, atol=_GRID_TOLERANCE_MM)):
                    issues["within_timepoint_spacing_mismatches"] += 1
                    timepoint_has_grid_mismatch = True
                if (affine.shape == (4, 4) and reference[2].shape == (4, 4)
                        and not np.allclose(affine, reference[2], rtol=0.0,
                                            atol=_GRID_TOLERANCE_MM)):
                    issues["within_timepoint_affine_mismatches"] += 1
                    timepoint_has_grid_mismatch = True
        if not timepoint_has_grid_mismatch:
            matched_timepoints += 1

    ordered_issues = [{"code": code, "count": int(count)}
                      for code, count in sorted(issues.items()) if count]
    return {
        "status": "blocked" if ordered_issues else "consistent_by_header_metadata",
        "summary": {
            "local_patient_directories": len(patient_dirs),
            "local_timepoint_directories": len(timepoint_dirs),
            "nifti_files_seen": total_files,
            "nifti_headers_read": sum(header_counts.values()) - issues["unreadable_nifti_headers"],
            "timepoints_without_shape_spacing_or_affine_mismatch": matched_timepoints,
            "header_counts_by_file_kind": dict(sorted(header_counts.items())),
            "t1c_shape_counts": dict(sorted(t1c_shapes.items())),
            "t1c_voxel_spacing_counts_mm": dict(sorted(t1c_spacings.items())),
            "grid_tolerance_mm": _GRID_TOLERANCE_MM,
        },
        "issues": ordered_issues,
        "interpretation_limit": (
            "NIfTI header metadata only; voxel arrays, image intensities, segmentation labels, registration quality, "
            "and clinical validity were not checked"
        ),
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--image-dir", type=Path, default=DEFAULT_IMAGE_DIR)
    args = parser.parse_args(argv)
    report = audit_nifti_headers(args.image_dir)
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0 if report["status"] == "consistent_by_header_metadata" else 2


if __name__ == "__main__":
    raise SystemExit(main())
