"""Outcome-blind CFB file, header, provenance, and origin-MRI preflight.

This tool does not authorize scoring. It reads NIfTI headers for all listed
files and reads voxel values only from each patient's t1 T1-Gd image to check
the frozen robust-normalization inputs. It never reads any GTV mask voxels or
t0/t2 MRI voxels. The t0, t1, and t2 GTV files in the crosswalk must already
be transformed into the t1 reference grid. Native t0/t2 MRI files need valid
headers but may retain their native affines.

Crosswalk CSV columns. CFB v4 paths must use this layout:
`<crosswalk directory>/<three-digit folder>/<t0|t1|t2>/`,
with exact filename `<numeric patient_id>_t0|1|2_t1gd|gtv.nii.gz`.
patient_id, identity_mapping_source,
t0_t1gd_path, t1_t1gd_path, t2_t1gd_path,
t0_gtv_path, t1_gtv_path, t2_gtv_path,
provenance_t0, provenance_source_t0,
provenance_t1, provenance_source_t1,
provenance_t2, provenance_source_t2,
review_status_t0, review_source_t0,
review_status_t1, review_source_t1,
review_status_t2, review_source_t2,
t0_to_t1_transform_documented, t0_to_t1_transform_source,
t0_to_t1_qc_pass, t0_to_t1_qc_source,
t2_to_t1_transform_documented, t2_to_t1_transform_source,
t2_to_t1_qc_pass, t2_to_t1_qc_source.

Paths may be absolute or relative to the crosswalk file. A successful report
means only that these preflight checks passed. It does not resolve source
release discrepancies, prove transform validity, prove clinical provenance,
or permit target-mask access.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Callable, Mapping

import numpy as np


IMAGE_COLUMNS = {
    "t0_t1gd_path": "t0_t1gd",
    "t1_t1gd_path": "t1_t1gd",
    "t2_t1gd_path": "t2_t1gd",
    "t0_gtv_path": "t0_gtv",
    "t1_gtv_path": "t1_gtv",
    "t2_gtv_path": "t2_gtv",
}
ALIGNMENT_REQUIRED = {"t1_t1gd", "t0_gtv", "t1_gtv", "t2_gtv"}
REQUIRED_COLUMNS = (
    "patient_id",
    "identity_mapping_source",
    *IMAGE_COLUMNS,
    "provenance_t0",
    "provenance_source_t0",
    "provenance_t1",
    "provenance_source_t1",
    "provenance_t2",
    "provenance_source_t2",
    "review_status_t0",
    "review_source_t0",
    "review_status_t1",
    "review_source_t1",
    "review_status_t2",
    "review_source_t2",
    "t0_to_t1_transform_documented",
    "t0_to_t1_transform_source",
    "t0_to_t1_qc_pass",
    "t0_to_t1_qc_source",
    "t2_to_t1_transform_documented",
    "t2_to_t1_transform_source",
    "t2_to_t1_qc_pass",
    "t2_to_t1_qc_source",
)
UNKNOWN_PROVENANCE = {"", "unknown", "unverified", "na", "n/a", "none"}
UNKNOWN_REVIEW = {"", "unknown", "unverified", "pending", "not_recorded", "na", "n/a", "none"}
HUMAN_PROVENANCE_ALIASES = {
    "humansegmentation", "human_segmentation", "human", "manual",
    "manual_segmentation", "clinical_manual", "clinical_manual_annotation",
}
GENERATED_PROVENANCE_ALIASES = {
    "gensegmentation", "gen_segmentation", "generated", "generated_segmentation",
    "nnunet_generated", "model_generated",
}
HUMAN_REVIEW_ALIASES = {
    "manual", "human_manual", "manual_annotation", "not_applicable",
    "not_applicable_manual", "not_applicable_human_manual",
}
GENERATED_REVIEW_ALIASES = {
    "reviewed", "reviewed_corrected", "reviewed_and_corrected",
    "reviewed_no_change", "not_reviewed", "unreviewed",
}


def _load_nifti(path: Path) -> Any:
    try:
        import nibabel as nib
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise RuntimeError("nibabel is required to inspect NIfTI headers") from exc
    return nib.load(str(path))


def _normalize_label(value: Any) -> str:
    return "_".join(str(value or "").strip().casefold().replace("-", " ").split())


def _parse_transform_documented(value: Any) -> bool | None:
    normalized = str(value or "").strip().lower()
    if normalized in {"1", "true", "yes", "documented"}:
        return True
    if normalized in {"0", "false", "no", "undocumented"}:
        return False
    return None


def _parse_qc_pass(value: Any) -> bool | None:
    normalized = str(value or "").strip().lower()
    if normalized in {"1", "true", "yes", "pass"}:
        return True
    if normalized in {"0", "false", "no", "fail"}:
        return False
    return None


def _path_from_value(value: Any, crosswalk_dir: Path) -> Path | None:
    raw = str(value or "").strip()
    if not raw:
        return None
    path = Path(raw).expanduser()
    return path.resolve() if path.is_absolute() else (crosswalk_dir / path).resolve()


def _path_identity_errors(path: Path, patient_id: str, role: str) -> list[str]:
    """Reject crosswalk swaps before loading any NIfTI header or voxel data."""
    visit, modality = role.split("_", maxsplit=1)
    visit_index = int(visit[1:])
    errors: list[str] = []
    if path.name != f"{patient_id}_{visit}_{modality}.nii.gz":
        errors.append(f"{role}_filename_identity_mismatch")
    if patient_id.isdigit():
        expected_patient_folder = patient_id.zfill(3)
        expected_visit_folder = visit
    else:
        # Non-numeric IDs are reserved for the isolated synthetic test fixtures.
        expected_patient_folder = patient_id
        expected_visit_folder = f"visit_{visit_index:02d}"
    if len(path.parents) < 2 or path.parent.parent.name != expected_patient_folder:
        errors.append(f"{role}_patient_directory_mismatch")
    if path.parent.name != expected_visit_folder:
        errors.append(f"{role}_visit_directory_mismatch")
    return errors


def _valid_header(image: Any) -> tuple[bool, dict[str, Any], str | None]:
    try:
        shape = tuple(int(axis) for axis in image.shape)
        affine = np.asarray(image.affine, dtype=np.float64)
        units = image.header.get_xyzt_units()[0]
    except Exception as exc:
        return False, {}, f"header_read_error:{type(exc).__name__}"

    if len(shape) != 3 or any(axis <= 0 for axis in shape):
        return False, {"shape": shape, "units": units}, "not_positive_3d"
    if affine.shape != (4, 4) or not np.all(np.isfinite(affine)):
        return False, {"shape": shape, "units": units}, "invalid_affine"
    determinant = float(np.linalg.det(affine[:3, :3]))
    spacing = np.linalg.norm(affine[:3, :3], axis=0)
    if not np.isfinite(determinant) or abs(determinant) <= 1e-10:
        return False, {"shape": shape, "units": units}, "singular_affine"
    if not np.all(np.isfinite(spacing)) or np.any(spacing <= 0):
        return False, {"shape": shape, "units": units}, "invalid_spacing"
    if units != "mm":
        return False, {"shape": shape, "units": units}, "spatial_units_not_mm"
    directions = affine[:3, :3] / spacing[np.newaxis, :]
    if not np.allclose(directions.T @ directions, np.eye(3), rtol=0.0, atol=1e-5):
        return False, {
            "shape": shape,
            "units": units,
            "spacing_mm": [float(value) for value in spacing],
        }, "non_orthogonal_affine"

    return True, {
        "shape": list(shape),
        "units": units,
        "spacing_mm": [float(value) for value in spacing],
        "affine": affine.tolist(),
    }, None


def _origin_t1_qc(image: Any) -> tuple[bool, dict[str, Any], str | None]:
    """Check t1 T1-Gd scale without using any mask or other visit's voxels."""
    try:
        values = np.asarray(image.dataobj, dtype=np.float32)
    except Exception as exc:
        return False, {}, f"origin_voxel_read_error:{type(exc).__name__}"
    if values.ndim != 3:
        return False, {"shape": list(values.shape)}, "origin_not_3d"
    if not np.all(np.isfinite(values)):
        return False, {"shape": list(values.shape)}, "origin_nonfinite_values"

    support = np.abs(values) > 1e-6
    count = int(np.count_nonzero(support))
    if count < 100:
        return False, {"nonzero_support_voxels": count}, "origin_too_few_supported_voxels"
    sample = values[support].astype(np.float64, copy=False)
    median = float(np.median(sample))
    scale = float(np.median(np.abs(sample - median))) * 1.4826
    scale_method = "median_absolute_deviation"
    if not np.isfinite(scale) or scale <= 1e-6:
        low, high = np.percentile(sample, [10.0, 90.0])
        scale = float((high - low) / 2.563)
        scale_method = "percentile_10_90_fallback"
    if not np.isfinite(scale) or scale <= 1e-6:
        return False, {
            "nonzero_support_voxels": count,
            "median": median,
            "scale": scale if math.isfinite(scale) else None,
            "scale_method": scale_method,
        }, "origin_no_usable_intensity_variation"
    return True, {
        "nonzero_support_voxels": count,
        "median": median,
        "robust_scale": scale,
        "scale_method": scale_method,
    }, None


def _field_nonempty(row: Mapping[str, Any], key: str) -> bool:
    return bool(str(row.get(key, "") or "").strip())


def _documented_check(row: Mapping[str, Any], prefix: str) -> tuple[bool, list[str]]:
    errors: list[str] = []
    documented = _parse_transform_documented(row.get(f"{prefix}_transform_documented"))
    qc_pass = _parse_qc_pass(row.get(f"{prefix}_qc_pass"))
    if documented is not True:
        errors.append(f"{prefix}_transform_not_documented")
    if not _field_nonempty(row, f"{prefix}_transform_source"):
        errors.append(f"{prefix}_transform_source_missing")
    if qc_pass is not True:
        errors.append(f"{prefix}_qc_not_passed")
    if not _field_nonempty(row, f"{prefix}_qc_source"):
        errors.append(f"{prefix}_qc_source_missing")
    return not errors, errors


def _provenance_check(row: Mapping[str, Any]) -> tuple[bool, list[str], list[str]]:
    values = [str(row.get(f"provenance_{visit}", "") or "").strip() for visit in ("t0", "t1", "t2")]
    normalized = [_normalize_label(value) for value in values]
    canonical: list[str | None] = []
    errors: list[str] = []
    for visit, value in zip(("t0", "t1", "t2"), normalized):
        if value in UNKNOWN_PROVENANCE:
            errors.append(f"provenance_{visit}_unknown")
            canonical.append(None)
        elif value in HUMAN_PROVENANCE_ALIASES:
            canonical.append("humanSegmentation")
        elif value in GENERATED_PROVENANCE_ALIASES:
            canonical.append("genSegmentation")
        else:
            canonical.append(None)
            errors.append(f"provenance_{visit}_invalid")
        if not _field_nonempty(row, f"provenance_source_{visit}"):
            errors.append(f"provenance_source_{visit}_missing")
    known = [value for value in canonical if value is not None]
    if len(set(known)) > 1:
        errors.append("provenance_classes_not_matched")
    return not errors, [value or "unknown_or_invalid" for value in canonical], errors


def _review_check(row: Mapping[str, Any]) -> tuple[bool, list[str], list[str]]:
    """Require explicit, source-backed review state for each visit mask."""
    values = [str(row.get(f"review_status_{visit}", "") or "").strip() for visit in ("t0", "t1", "t2")]
    normalized = [_normalize_label(value) for value in values]
    provenance = [_normalize_label(row.get(f"provenance_{visit}", "")) for visit in ("t0", "t1", "t2")]
    errors: list[str] = []
    canonical: list[str | None] = []
    for visit, value, provenance_value in zip(("t0", "t1", "t2"), normalized, provenance):
        if value in UNKNOWN_REVIEW:
            errors.append(f"mask_review_status_{visit}_unknown")
            canonical.append(None)
        elif provenance_value in HUMAN_PROVENANCE_ALIASES and value in HUMAN_REVIEW_ALIASES:
            canonical.append("manual")
        elif provenance_value in GENERATED_PROVENANCE_ALIASES and value in GENERATED_REVIEW_ALIASES:
            if value in {"reviewed_corrected", "reviewed_and_corrected", "reviewed_no_change"}:
                canonical.append("reviewed")
            else:
                canonical.append("not_reviewed" if value == "unreviewed" else value)
        else:
            canonical.append(None)
            errors.append(f"mask_review_status_{visit}_invalid_for_provenance")
        if not _field_nonempty(row, f"review_source_{visit}"):
            errors.append(f"mask_review_source_{visit}_missing")
    known = [value for value in canonical if value is not None]
    if len(set(known)) > 1:
        errors.append("mask_review_statuses_not_matched")
    return not errors, [value or "unknown_or_invalid" for value in canonical], errors


def preflight_rows(
    rows: list[dict[str, str]],
    crosswalk_dir: Path,
    package_root: Path | None = None,
    image_loader: Callable[[Path], Any] | None = None,
) -> list[dict[str, Any]]:
    """Evaluate explicit rows; no visit is selected or ranked from outcomes."""
    trusted_root = (package_root or crosswalk_dir).resolve()
    if not trusted_root.is_dir():
        raise ValueError(f"package root is not a directory: {trusted_root}")
    loader = _load_nifti if image_loader is None else image_loader
    patient_counts: dict[str, int] = {}
    for row in rows:
        patient = str(row.get("patient_id", "") or "").strip()
        patient_counts[patient] = patient_counts.get(patient, 0) + 1

    results: list[dict[str, Any]] = []
    for row_number, row in enumerate(rows, start=2):
        patient = str(row.get("patient_id", "") or "").strip()
        result: dict[str, Any] = {
            "source_row": row_number,
            "patient_id": patient,
            "candidate_id": hashlib.sha256(patient.encode("utf-8")).hexdigest()[:16] if patient else None,
            "file_headers": {},
            "errors": [],
            "target_mask_voxels_read": False,
            "t0_t2_mri_voxels_read": False,
            "scoring_authorized": False,
            "metadata_evidence": {
                key: str(row.get(key, "") or "").strip()
                for key in REQUIRED_COLUMNS
                if key.endswith("_source") or key == "identity_mapping_source"
            },
        }
        errors: list[str] = result["errors"]
        if not patient:
            errors.append("patient_id_missing")
        elif not (patient.isdigit() or re.fullmatch(r"deidentified-case-\d{3}", patient)):
            errors.append("patient_id_not_supported")
        if patient_counts.get(patient, 0) != 1:
            errors.append("patient_id_not_unique_in_crosswalk")
        if not _field_nonempty(row, "identity_mapping_source"):
            errors.append("identity_mapping_source_missing")

        provenance_pass, provenance_values, provenance_errors = _provenance_check(row)
        result["provenance"] = provenance_values
        result["provenance_pass"] = provenance_pass
        errors.extend(provenance_errors)
        review_pass, review_values, review_errors = _review_check(row)
        result["mask_review_status"] = review_values
        result["mask_review_pass"] = review_pass
        errors.extend(review_errors)
        for prefix in ("t0_to_t1", "t2_to_t1"):
            passed, doc_errors = _documented_check(row, prefix)
            result[f"{prefix}_documented_qc_pass"] = passed
            errors.extend(doc_errors)

        loaded: dict[str, Any] = {}
        headers: dict[str, dict[str, Any]] = result["file_headers"]
        path_keys: dict[str, Path] = {}
        path_errors: list[str] = []
        for column, role in IMAGE_COLUMNS.items():
            path = _path_from_value(row.get(column), crosswalk_dir)
            if path is None:
                path_errors.append(f"{role}_path_missing")
                continue
            if not path.name.lower().endswith(".nii.gz"):
                path_errors.append(f"{role}_not_nifti_gz")
                continue
            try:
                path.relative_to(trusted_root)
            except ValueError:
                path_errors.append(f"{role}_outside_package_root")
                continue
            path_errors.extend(_path_identity_errors(path, patient, role))
            if not path.is_file():
                path_errors.append(f"{role}_file_not_found")
                continue
            path_keys[role] = path

        # Validate every path before loading any file. A swapped role must not
        # trigger even a header read from another visit or modality.
        errors.extend(path_errors)
        if not path_errors and len(path_keys) == len(IMAGE_COLUMNS):
            for role, path in path_keys.items():
                try:
                    image = loader(path)
                except Exception as exc:
                    errors.append(f"{role}_header_load_error:{type(exc).__name__}")
                    continue
                valid, header, header_error = _valid_header(image)
                headers[role] = {**header, "header_pass": valid}
                loaded[role] = image
                if header_error:
                    errors.append(f"{role}_{header_error}")

        reference = headers.get("t1_t1gd")
        if reference and reference.get("header_pass"):
            ref_shape = reference["shape"]
            ref_affine = np.asarray(reference["affine"], dtype=np.float64)
            for role, header in headers.items():
                header["alignment_required"] = role in ALIGNMENT_REQUIRED
                if role not in ALIGNMENT_REQUIRED:
                    header["aligned_to_t1_reference"] = None
                    continue
                if not header.get("header_pass"):
                    continue
                aligned = header["shape"] == ref_shape and np.allclose(
                    np.asarray(header["affine"], dtype=np.float64), ref_affine,
                    rtol=0.0, atol=1e-3,
                )
                header["aligned_to_t1_reference"] = bool(aligned)
                if not aligned:
                    errors.append(f"{role}_not_aligned_to_t1_reference")

        if "t1_t1gd" in loaded and headers.get("t1_t1gd", {}).get("header_pass"):
            qc_pass, diagnostics, qc_error = _origin_t1_qc(loaded["t1_t1gd"])
            result["origin_t1_t1gd_qc"] = {"pass": qc_pass, **diagnostics}
            if qc_error:
                errors.append(qc_error)
        else:
            result["origin_t1_t1gd_qc"] = {"pass": False, "error": "origin_header_unavailable"}
            errors.append("origin_header_unavailable")

        result["headers_pass"] = (
            len(headers) == len(IMAGE_COLUMNS)
            and all(header.get("header_pass") for header in headers.values())
            and all(
                header.get("aligned_to_t1_reference") is True
                for role, header in headers.items()
                if role in ALIGNMENT_REQUIRED
            )
        )
        result["metadata_gates_pass"] = (
            _field_nonempty(row, "identity_mapping_source")
            and provenance_pass
            and review_pass
            and result["t0_to_t1_documented_qc_pass"]
            and result["t2_to_t1_documented_qc_pass"]
        )
        result["origin_qc_pass"] = bool(result["origin_t1_t1gd_qc"].get("pass"))
        result["preflight_pass"] = (
            not errors and result["headers_pass"] and result["metadata_gates_pass"] and result["origin_qc_pass"]
        )
        results.append(result)
    return results


def read_crosswalk(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if len(fieldnames) != len(set(fieldnames)):
            raise ValueError("crosswalk contains duplicate column names")
        columns = set(fieldnames)
        missing = sorted(set(REQUIRED_COLUMNS) - columns)
        if missing:
            raise ValueError(f"crosswalk missing required columns: {', '.join(missing)}")
        output: list[dict[str, str]] = []
        for row in reader:
            if None in row:
                raise ValueError("crosswalk contains a row with more values than columns")
            if not any(str(value or "").strip() for value in row.values()):
                raise ValueError("crosswalk contains a blank row")
            output.append({str(key): str(value or "") for key, value in row.items()})
        return output


def run_preflight(
    crosswalk_path: Path,
    output_dir: Path,
    package_root: Path,
    protocol_path: Path,
) -> dict[str, Any]:
    rows = read_crosswalk(crosswalk_path)
    if not rows:
        raise ValueError("crosswalk has no candidate rows")
    if output_dir.exists():
        raise FileExistsError(f"output directory already exists: {output_dir}")
    trusted_root = package_root.resolve()
    if not trusted_root.is_dir():
        raise ValueError(f"package root is not a directory: {trusted_root}")
    protocol_resolved = protocol_path.resolve()
    if not protocol_resolved.is_file():
        raise ValueError(f"protocol file does not exist: {protocol_resolved}")
    results = preflight_rows(rows, crosswalk_path.resolve().parent, package_root=trusted_root)
    crosswalk_digest = hashlib.sha256(crosswalk_path.read_bytes()).hexdigest()
    script_path = Path(__file__).resolve()
    summary = {
        "source_package_root": str(trusted_root),
        "source_crosswalk": str(crosswalk_path.resolve()),
        "source_crosswalk_sha256": crosswalk_digest,
        "preflight_script": str(script_path),
        "preflight_script_sha256": hashlib.sha256(script_path.read_bytes()).hexdigest(),
        "analysis_protocol": str(protocol_resolved),
        "analysis_protocol_sha256": hashlib.sha256(protocol_resolved.read_bytes()).hexdigest(),
        "candidate_rows": len(results),
        "preflight_pass_rows": sum(bool(row["preflight_pass"]) for row in results),
        "preflight_fail_rows": sum(not bool(row["preflight_pass"]) for row in results),
        "target_mask_voxels_read": False,
        "t0_t2_mri_voxels_read": False,
        "scoring_authorized": False,
        "limitations": [
            "Header alignment does not prove registration provenance or anatomical registration quality.",
            "All three GTV mask files must already be transformed to the t1 reference grid; this tool does not transform masks.",
            "Metadata source fields are recorded but are not independently verified by this tool.",
            "This preflight does not resolve TCIA release crosswalk or collection-count discrepancies.",
            "This preflight does not inspect target masks or compute outcomes.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    (output_dir / "rows.json").write_text(json.dumps(results, indent=2), encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", required=True, type=Path, help="explicit outcome-blind file crosswalk CSV")
    parser.add_argument("--out-dir", required=True, type=Path, help="new directory for preflight report files")
    parser.add_argument("--package-root", required=True, type=Path, help="trusted root of the downloaded CFB package")
    parser.add_argument("--protocol", required=True, type=Path, help="frozen protocol file to hash into the report")
    args = parser.parse_args()
    summary = run_preflight(
        args.crosswalk,
        args.out_dir,
        package_root=args.package_root,
        protocol_path=args.protocol,
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
