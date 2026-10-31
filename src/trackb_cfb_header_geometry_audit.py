"""Check CFB-GBM NIfTI geometry from headers only; never read voxel arrays.

Run with --package-root pointing to the extracted CFB-GBM version 4 folder.
The selected patient IDs come from the existing paginated MRI review manifest.
This check verifies file presence, basic NIfTI-1 header validity, millimeter
units, and same-grid geometry. It does not establish registration quality,
mask provenance, valid-negative semantics, or scoring authorization.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import struct
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_MANIFEST = PROJECT_ROOT / "output" / "cfb_mri_paginated_review_full_20261031" / "manifest.json"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "trackb_cfb_header_geometry_audit_20261031.json"
AFFINE_ATOL = 1e-3


def _read_header(path: Path) -> dict[str, Any]:
    """Read and decode only the first 348 uncompressed NIfTI-1 header bytes."""
    try:
        with gzip.open(path, "rb") as stream:
            raw = stream.read(348)
    except (OSError, EOFError) as exc:
        raise ValueError(f"gzip_header_read_error:{type(exc).__name__}") from exc
    if len(raw) != 348:
        raise ValueError(f"short_nifti1_header:{len(raw)}")

    little_size = struct.unpack_from("<i", raw, 0)[0]
    big_size = struct.unpack_from(">i", raw, 0)[0]
    if little_size == 348:
        endian = "<"
    elif big_size == 348:
        endian = ">"
    else:
        raise ValueError("unsupported_nifti_header_size")

    magic = raw[344:348]
    if magic not in (b"n+1\x00", b"ni1\x00"):
        raise ValueError("unexpected_nifti1_magic")

    dims = struct.unpack_from(endian + "8h", raw, 40)
    if dims[0] < 3 or any(int(value) <= 0 for value in dims[1:4]):
        raise ValueError("not_positive_3d_image")
    if any(int(value) != 1 for value in dims[4:dims[0] + 1]):
        raise ValueError("non_singleton_extra_dimension")
    shape = [int(value) for value in dims[1:4]]

    unit_code = raw[123] & 0x07
    units = {1: "m", 2: "mm", 3: "um"}.get(unit_code, "unknown")
    if units != "mm":
        raise ValueError(f"spatial_units_not_mm:{units}")

    qform_code = struct.unpack_from(endian + "h", raw, 252)[0]
    sform_code = struct.unpack_from(endian + "h", raw, 254)[0]
    if sform_code > 0:
        rows = [
            struct.unpack_from(endian + "4f", raw, offset)
            for offset in (280, 296, 312)
        ]
        affine = [list(map(float, row)) for row in rows] + [[0.0, 0.0, 0.0, 1.0]]
    elif qform_code > 0:
        affine = _qform_affine(raw, endian)
    else:
        raise ValueError("no_sform_or_qform")

    if any(not math.isfinite(value) for row in affine for value in row):
        raise ValueError("nonfinite_affine")
    return {
        "shape": shape,
        "affine": affine,
        "units": units,
        "qform_code": int(qform_code),
        "sform_code": int(sform_code),
    }


def _qform_affine(raw: bytes, endian: str) -> list[list[float]]:
    pixdim = struct.unpack_from(endian + "8f", raw, 76)
    b, c, d = struct.unpack_from(endian + "3f", raw, 256)
    x, y, z = struct.unpack_from(endian + "3f", raw, 268)
    a2 = 1.0 - (b * b + c * c + d * d)
    if a2 < 0.0:
        norm = math.sqrt(b * b + c * c + d * d)
        if norm == 0.0:
            raise ValueError("invalid_qform_quaternion")
        b, c, d = b / norm, c / norm, d / norm
        a = 0.0
    else:
        a = math.sqrt(a2)

    rotation = [
        [a*a + b*b - c*c - d*d, 2*b*c - 2*a*d, 2*b*d + 2*a*c],
        [2*b*c + 2*a*d, a*a + c*c - b*b - d*d, 2*c*d - 2*a*b],
        [2*b*d - 2*a*c, 2*c*d + 2*a*b, a*a + d*d - c*c - b*b],
    ]
    scales = [abs(float(pixdim[1])), abs(float(pixdim[2])), abs(float(pixdim[3]))]
    if any(not math.isfinite(value) or value <= 0.0 for value in scales):
        raise ValueError("invalid_qform_spacing")
    if pixdim[0] < 0:
        scales[2] *= -1.0
    affine = [
        [rotation[row][col] * scales[col] for col in range(3)] + [float(offset)]
        for row, offset in enumerate((x, y, z))
    ]
    return affine + [[0.0, 0.0, 0.0, 1.0]]


def _same_grid(left: dict[str, Any], right: dict[str, Any]) -> bool:
    if left["shape"] != right["shape"]:
        return False
    return all(
        abs(float(left["affine"][row][col]) - float(right["affine"][row][col])) <= AFFINE_ATOL
        for row in range(4)
        for col in range(4)
    )


def run_audit(package_root: Path, manifest_path: Path, output_path: Path) -> dict[str, Any]:
    package_root = package_root.expanduser().resolve()
    manifest_path = manifest_path.expanduser().resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    patient_ids = [int(value) for value in manifest["selected_patient_ids"]]
    results: list[dict[str, Any]] = []
    files_checked = 0

    for patient_id in patient_ids:
        folder = f"{patient_id:03d}"
        visits: dict[str, dict[str, Any]] = {}
        case_errors: list[str] = []
        for visit in ("t0", "t1", "t2"):
            base = package_root / folder / visit
            paths = {
                "t1gd": base / f"{patient_id}_{visit}_t1gd.nii.gz",
                "gtv": base / f"{patient_id}_{visit}_gtv.nii.gz",
                "brain_mask": base / f"{patient_id}_{visit}_brain_mask.nii.gz",
            }
            for role, path in paths.items():
                if not path.is_file():
                    case_errors.append(f"{visit}_{role}_missing")
                    continue
                files_checked += 1
                try:
                    visits.setdefault(visit, {})[role] = _read_header(path)
                except ValueError as exc:
                    case_errors.append(f"{visit}_{role}_{exc}")

        reference = visits.get("t0", {}).get("t1gd")
        if reference is None:
            case_errors.append("t0_reference_header_missing")
        else:
            for visit, roles in visits.items():
                for role, header in roles.items():
                    if not _same_grid(reference, header):
                        case_errors.append(f"{visit}_{role}_not_on_t0_t1gd_grid")

        results.append({
            "patient_id": patient_id,
            "header_file_count": sum(len(roles) for roles in visits.values()),
            "same_grid_as_t0_t1gd": not any("not_on_t0_t1gd_grid" in item for item in case_errors),
            "errors": case_errors,
        })

    missing_cases = [row["patient_id"] for row in results if row["errors"]]
    try:
        manifest_display = manifest_path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        manifest_display = manifest_path.name
    summary = {
        "status": "complete" if len(results) == len(patient_ids) else "incomplete",
        "scope": "NIfTI-1 headers only; reads first 348 uncompressed bytes per file; no voxel arrays",
        "package_root_label": f"{package_root.parent.name}/{package_root.name}",
        "manifest_path": manifest_display,
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "candidate_count": len(patient_ids),
        "expected_headers": len(patient_ids) * 3 * 3,
        "headers_read": files_checked,
        "same_grid_cases": sum(bool(row["same_grid_as_t0_t1gd"]) and not row["errors"] for row in results),
        "cases_with_errors": len(missing_cases),
        "error_patient_ids": missing_cases,
        "gtv_voxel_values_read": False,
        "brain_mask_voxel_values_read": False,
        "mri_voxel_values_read": False,
        "registration_quality_pass": False,
        "valid_negative_support_pass": False,
        "scoring_authorized": False,
        "cases": results,
    }
    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--package-root", required=True, type=Path)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run_audit(args.package_root, args.manifest, args.output), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
