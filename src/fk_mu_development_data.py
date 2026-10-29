"""Load the frozen MU first-two-to-third PDE development cohort behind explicit QC gates.

Metadata selection is recomputed before any tumor-mask voxel array is accessed. The supplied evidence JSON
must separately establish source-backed longitudinal registration QC and valid target-label support. This
module does not infer either fact from matching headers or MRI intensity summaries.

Required evidence JSON fields (all values are checked before mask voxel access)::

    {
      "schema_version": 1,
      "manifest_sha256": {"scan_manifest": "...", "patient_manifest": "...",
                          "metadata_screen": "...", "mri_only_qc": "..."},
      "selected_patient_ids_sha256_sorted_nul_joined": "...",
      "registration": {"source_backed": true, "all_selected_patients_covered": true,
                       "qc_status": "pass", "method": "...", "reference": "...",
                       "evidence_reference": "...", "evidence_sha256": "..."},
      "target_support": {"source_backed": true, "all_selected_patients_covered": true,
                         "zero_is_valid_negative": true, "definition": "...",
                         "evidence_reference": "...", "evidence_sha256": "..."},
      "origin_support": {"source": "origin-only", "not_target_derived": true}
    }

Origin support must also be passed as ``origin_support_by_patient``. Each value is an ``OriginSupport`` with
``mask`` (boolean 3D array), ``affine`` (4x4 affine on the 2-mm output grid), and ``source`` (a nonempty
origin-only source description). No target mask is used to define or trim this support. Callers must also pass an
explicit ``allowed_evidence_roots`` allowlist; each evidence reference must resolve to an existing local file
under one of those roots and match its declared SHA-256.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
import csv
import hashlib
import json
from pathlib import Path
import re
from types import MappingProxyType
from typing import Iterable, Mapping
import weakref

import numpy as np

from src.fk_patient_cv import PatientSeries


PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_SCAN_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "scan_manifest.csv"
DEFAULT_PATIENT_MANIFEST = PROJECT_ROOT / "data" / "manifests" / "patient_manifest.csv"
DEFAULT_DATA_ROOT = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
DEFAULT_METADATA_REPORT = PROJECT_ROOT / "output" / "fk_mu_target_window_screen_20261028" / "metadata_screen.json"
DEFAULT_MRI_QC_REPORT = PROJECT_ROOT / "output" / "fk_mu_mri_only_qc_20261028_v2" / "mri_only_qc.json"

TARGET_GAP_MIN_DAYS = 76.0
TARGET_GAP_MAX_DAYS = 104.0
GRID_TOLERANCE_MM = 1e-4
EXPECTED_PATIENTS = 18

# Snapshot hashes bind this loader to the currently reviewed metadata and MRI-only reports.
EXPECTED_SHA256 = {
    "scan_manifest": "ed2307fc23f7e03fb833e4031124c2b83642604d192ccbb1528fe3050966dcf5",
    "patient_manifest": "24256d0fb5a352349725ffd73d7fb100d5ecad017cedbad32878b7e35311598a",
    "metadata_screen": "12bdd4b79d5bcbdc59cd0886eed4c0487a8cb618fbc54c60d7d66013c67547bf",
    "mri_only_qc": "b0825ad95843e8a6f66c36cf043754280e67197d034d0e48afeba45f06dd76c8",
}
EXPECTED_PROJECTED_MANIFEST_SHA256 = {
    "scan_selection_columns": "f5f2c1b10d5c30f97f5506230bf6ef6bb0a9dc2236fdba9b71ad17d9d835fb73",
    "patient_selection_columns": "be804e92deb8a24bd2e2269cfac1dced1e5c7b47233a8b0e3240011e341afbbc",
}
EXPECTED_MRI_QC = {
    "selected_patient_count": 18,
    "sequence_files_present": 216,
    "sequence_files_expected": 216,
    "t1n_headers_aligned": True,
    "t1n_volumes_read": 54,
    "tumor_mask_files_opened": 0,
    "tumor_mask_voxel_arrays_read": 0,
}


@dataclass(frozen=True)
class OriginSupport:
    """Caller-supplied support derived only from the forecast-origin scan."""

    mask: np.ndarray
    affine: np.ndarray
    source: str


@dataclass(frozen=True, eq=False)
class MUDevelopmentCohort:
    """Loaded histories, patient folds, origin domains, valid scoring support, and access audit."""

    records: tuple[PatientSeries, ...]
    fold_by_patient: Mapping[str, int]
    origin_domain_by_patient: Mapping[str, np.ndarray]
    target_valid_support_by_patient: Mapping[str, np.ndarray]
    audit: dict


_LOADER_ISSUED_COHORTS = weakref.WeakSet()


def is_loader_issued_cohort(value) -> bool:
    """Return whether this exact in-memory cohort object came from the gated loader."""
    return isinstance(value, MUDevelopmentCohort) and value in _LOADER_ISSUED_COHORTS


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _projected_rows_sha256(rows: list[dict[str, str]]) -> str:
    payload = json.dumps(rows, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _read_csv(path: Path, fields: set[str]) -> list[dict[str, str]]:
    with Path(path).open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration:
            raise ValueError(f"manifest is empty: {path}") from None
        indices = {name: index for index, name in enumerate(header)}
        missing = fields - set(indices)
        if missing:
            raise ValueError(f"manifest is missing required columns: {sorted(missing)}")
        rows = []
        for row_number, raw in enumerate(reader, start=2):
            if len(raw) < len(header):
                raise ValueError(f"malformed manifest row {row_number}: {path}")
            rows.append({name: raw[indices[name]] for name in fields})
    if not rows:
        raise ValueError(f"manifest has no data rows: {path}")
    return rows


def _truth(value: str) -> bool:
    return str(value or "").strip().lower() in {"true", "1", "yes"}


def _selected_id_hash(patient_ids) -> str:
    return hashlib.sha256("\0".join(sorted(map(str, patient_ids))).encode("utf-8")).hexdigest()


def _validate_snapshot(scan_path: Path, patient_path: Path, metadata_path: Path, mri_path: Path) -> tuple[dict, dict]:
    paths = {
        "scan_manifest": scan_path,
        "patient_manifest": patient_path,
        "metadata_screen": metadata_path,
        "mri_only_qc": mri_path,
    }
    actual = {name: _sha256_file(path) for name, path in paths.items()}
    if actual != EXPECTED_SHA256:
        changed = [name for name in actual if actual[name] != EXPECTED_SHA256[name]]
        raise ValueError(f"reviewed metadata/report snapshot mismatch: {changed}")

    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    mri_qc = json.loads(mri_path.read_text(encoding="utf-8"))
    counts = metadata.get("counts", {})
    if (metadata.get("selected_patient_count_after_header_check") != EXPECTED_PATIENTS
            or counts.get("gbm_in_76_104_day_target_window") != EXPECTED_PATIENTS
            or counts.get("gbm_target_window_mask_and_origin_brain_files_present") != EXPECTED_PATIENTS
            or metadata.get("header_check", {}).get("header_pass_candidates") != EXPECTED_PATIENTS
            or metadata.get("header_check", {}).get("nifti_headers_read") != EXPECTED_PATIENTS * 4
            or metadata.get("header_check", {}).get("voxel_arrays_read") is not False):
        raise ValueError("frozen metadata screen counts/access record do not match reviewed expectations")
    source_hashes = metadata.get("source_hashes", {})
    if (source_hashes.get("scan_selection_columns_sha256")
            != EXPECTED_PROJECTED_MANIFEST_SHA256["scan_selection_columns"]
            or source_hashes.get("patient_selection_columns_sha256")
            != EXPECTED_PROJECTED_MANIFEST_SHA256["patient_selection_columns"]):
        raise ValueError("frozen metadata screen selection-column hashes do not match")

    mri_selection = mri_qc.get("selection", {})
    mri_files = mri_qc.get("file_availability", {})
    mri_geometry = mri_qc.get("t1n_geometry", {})
    mri_access = mri_qc.get("data_access", {})
    mri_source_hashes = mri_qc.get("source_hashes", {})
    checks = {
        "selected_patient_count": mri_selection.get("selected_patient_count"),
        "sequence_files_present": mri_files.get("sequence_files_present"),
        "sequence_files_expected": mri_files.get("sequence_files_expected"),
        "t1n_headers_aligned": mri_geometry.get("all_first_three_visit_headers_align_within_mm_tolerance"),
        "t1n_volumes_read": mri_access.get("t1n_mri_volumes_read"),
        "tumor_mask_files_opened": mri_access.get("tumor_mask_files_opened"),
        "tumor_mask_voxel_arrays_read": mri_access.get("tumor_mask_voxel_arrays_read"),
    }
    if checks != EXPECTED_MRI_QC:
        raise ValueError("MRI-only QC report counts/access record do not match reviewed expectations")
    if mri_source_hashes.get("metadata_screen_sha256") != actual["metadata_screen"]:
        raise ValueError("MRI-only QC report is not tied to the reviewed metadata-screen report")
    if mri_selection.get("target_mask_values_used_for_selection") is not False:
        raise ValueError("MRI-only QC report does not attest outcome-blind cohort selection")
    if mri_qc.get("gate_status") != "descriptive_screen_only; no registration or valid-support gate cleared":
        raise ValueError("MRI-only report status changed; re-review its interpretation before loading masks")
    return metadata, mri_qc


def _evidence_hash(value, field: str) -> str:
    result = str(value or "").strip().lower()
    if not re.fullmatch(r"[0-9a-f]{64}", result):
        raise ValueError(f"registration/support evidence requires a SHA-256 field: {field}")
    return result


def _resolve_local_evidence(reference: str, *, allowed_roots: Iterable[Path], field: str) -> Path:
    reference = str(reference or "").strip()
    if not reference:
        raise PermissionError(f"evidence gate requires a local evidence file: {field}")
    raw_path = Path(reference)
    candidate = raw_path if raw_path.is_absolute() else PROJECT_ROOT / raw_path
    resolved = candidate.resolve()
    roots = tuple(Path(root).resolve() for root in allowed_roots)
    if not roots:
        raise PermissionError("evidence gate requires an explicit nonempty local-root allowlist")
    if not any(resolved == root or root in resolved.parents for root in roots):
        raise PermissionError(f"evidence file is outside explicitly allowed roots: {field}")
    if not resolved.is_file():
        raise PermissionError(f"evidence file does not exist locally: {field}")
    return resolved


def _validate_evidence_file(section: dict, *, field: str, allowed_roots: Iterable[Path]) -> str:
    path = _resolve_local_evidence(
        section.get("evidence_reference"), allowed_roots=allowed_roots, field=f"{field}.evidence_reference"
    )
    declared = _evidence_hash(section.get("evidence_sha256"), f"{field}.evidence_sha256")
    actual = _sha256_file(path)
    if actual != declared:
        raise PermissionError(f"evidence file hash mismatch: {field}")
    return actual


def _validate_evidence_gate(evidence_path: Path, selected_ids: list[str], manifest_hashes: dict,
                            allowed_roots: Iterable[Path]) -> dict:
    evidence_path = _resolve_local_evidence(
        str(evidence_path), allowed_roots=allowed_roots, field="evidence_document"
    )
    if not evidence_path.is_file():
        raise PermissionError(
            "patient mask loading blocked: provide an approved evidence JSON with source-backed registration "
            "QC and valid target-support evidence"
        )
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if evidence.get("schema_version") != 1:
        raise PermissionError("patient mask loading blocked: unsupported evidence schema")
    declared_hashes = evidence.get("manifest_sha256")
    if not isinstance(declared_hashes, dict):
        raise PermissionError("patient mask loading blocked: evidence manifest hashes are missing")
    for key, actual in manifest_hashes.items():
        if str(declared_hashes.get(key, "")).lower() != actual:
            raise PermissionError(f"patient mask loading blocked: evidence hash mismatch for {key}")
    expected_ids_hash = _selected_id_hash(selected_ids)
    if evidence.get("selected_patient_ids_sha256_sorted_nul_joined") != expected_ids_hash:
        raise PermissionError("patient mask loading blocked: evidence does not cover this selected cohort")

    registration = evidence.get("registration")
    if not isinstance(registration, dict):
        raise PermissionError("patient mask loading blocked: registration evidence section is missing")
    if (registration.get("source_backed") is not True
            or registration.get("all_selected_patients_covered") is not True
            or registration.get("qc_status") != "pass"
            or not str(registration.get("method", "")).strip()
            or not str(registration.get("reference", "")).strip()
            or not str(registration.get("evidence_reference", "")).strip()):
        raise PermissionError("patient mask loading blocked: source-backed longitudinal registration QC did not pass")
    registration_hash = _validate_evidence_file(
        registration, field="registration", allowed_roots=allowed_roots
    )

    support = evidence.get("target_support")
    if not isinstance(support, dict):
        raise PermissionError("patient mask loading blocked: target-support evidence section is missing")
    if (support.get("source_backed") is not True
            or support.get("all_selected_patients_covered") is not True
            or support.get("zero_is_valid_negative") is not True
            or support.get("coverage_scope") != "entire_image_grid"
            or not str(support.get("definition", "")).strip()
            or not str(support.get("evidence_reference", "")).strip()):
        raise PermissionError(
            "patient mask loading blocked: source-backed full-grid valid target-label support was not established"
        )
    support_hash = _validate_evidence_file(
        support, field="target_support", allowed_roots=allowed_roots
    )

    origin_support = evidence.get("origin_support")
    if (not isinstance(origin_support, dict)
            or origin_support.get("source") != "origin-only"
            or origin_support.get("not_target_derived") is not True):
        raise PermissionError("patient mask loading blocked: origin-only domain/support rule is not attested")
    return {"document": evidence, "registration_evidence_sha256": registration_hash,
            "target_support_evidence_sha256": support_hash,
            "evidence_document_sha256": _sha256_file(evidence_path)}


def _header(path: Path, data_root: Path) -> dict:
    try:
        path.resolve().relative_to(data_root.resolve())
    except ValueError as exc:
        raise ValueError("selected file path escapes the MU data root") from exc
    if not path.is_file():
        raise FileNotFoundError(path)
    import nibabel as nib

    image = nib.load(str(path), mmap=True)
    shape = tuple(int(value) for value in image.shape)
    affine = np.asarray(image.affine, dtype=float)
    spacing = tuple(float(value) for value in image.header.get_zooms()[:3])
    unit = image.header.get_xyzt_units()[0]
    if len(shape) != 3 or affine.shape != (4, 4) or not np.isfinite(affine).all():
        raise ValueError(f"invalid selected 3D NIfTI header: {path.name}")
    if unit != "mm" or not np.isfinite(spacing).all() or any(value <= 0 for value in spacing):
        raise ValueError(f"selected NIfTI must have finite millimeter spacing: {path.name}")
    if not np.allclose(spacing, (1.0, 1.0, 1.0), rtol=0.0, atol=GRID_TOLERANCE_MM):
        raise ValueError(f"selected MU files no longer match the reviewed 1-mm grid: {path.name}")
    if abs(float(np.linalg.det(affine[:3, :3]))) < 1e-12:
        raise ValueError(f"selected NIfTI has singular affine: {path.name}")
    return {"shape": shape, "affine": affine, "spacing_mm": spacing}


def _paths(patient_id: str, scans: list[int], data_root: Path) -> tuple[list[Path], Path]:
    patient_dir = data_root / patient_id
    mask_paths = [patient_dir / f"Timepoint_{n}" / f"{patient_id}_Timepoint_{n}_tumorMask.nii.gz"
                  for n in scans]
    origin_brain = (patient_dir / f"Timepoint_{scans[1]}" /
                    f"{patient_id}_Timepoint_{scans[1]}_brain_t1n.nii.gz")
    return mask_paths, origin_brain


def _metadata_candidates(scan_rows: list[dict[str, str]], patient_rows: list[dict[str, str]]) -> tuple[list[dict], Counter]:
    gbm_rows = [row for row in patient_rows if row["dataset"] == "MU" and _truth(row["is_gbm"])]
    patient_by_id = {row["patient_id"]: row for row in gbm_rows}
    if len(patient_by_id) != len(gbm_rows):
        raise ValueError("duplicate MU GBM patient rows in patient manifest")
    visits: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in scan_rows:
        if row["dataset"] == "MU" and _truth(row["mask_exists"]) and row["patient_id"] in patient_by_id:
            visits[row["patient_id"]].append(row)

    counts: Counter = Counter()
    candidates = []
    for patient_id, patient_visits in visits.items():
        try:
            ordered = sorted(patient_visits, key=lambda row: int(row["scan_number"]))
        except (TypeError, ValueError):
            counts["invalid_scan_number_rows"] += 1
            continue
        if len(ordered) < 3:
            continue
        counts["gbm_patients_with_three_mask_available_visits"] += 1
        first = ordered[:3]
        try:
            scans = [int(row["scan_number"]) for row in first]
            days = [float(row["day"]) for row in first]
        except (TypeError, ValueError):
            counts["first_three_dates_or_scan_numbers_invalid"] += 1
            continue
        if len(set(scans)) != 3 or not np.isfinite(days).all() or not (days[0] < days[1] < days[2]):
            counts["first_three_dates_or_scan_numbers_invalid"] += 1
            continue
        counts["gbm_with_valid_first_three_dates"] += 1
        gap = days[2] - days[1]
        if not TARGET_GAP_MIN_DAYS <= gap <= TARGET_GAP_MAX_DAYS:
            continue
        counts["gbm_in_76_104_day_target_window"] += 1
        fold_raw = patient_by_id[patient_id].get("fold", "").strip()
        try:
            fold = int(fold_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"selected patient is missing its existing patient_manifest fold: {patient_id}") from exc
        if str(fold) != fold_raw or not 0 <= fold < 5:
            raise ValueError(f"selected patient has an invalid existing fold: {patient_id}")
        candidates.append({"patient_id": patient_id, "scan_numbers": scans, "days": days, "gap_days": gap,
                           "fold": fold})
    return sorted(candidates, key=lambda candidate: candidate["patient_id"]), counts


def _verify_origin_supports(origin_support_by_patient: Mapping[str, OriginSupport], candidates: list[dict],
                            reference_headers: dict[str, dict]) -> dict[str, np.ndarray]:
    if origin_support_by_patient is None:
        raise PermissionError("scoring refused: caller must provide origin-only support for every selected patient")
    expected_ids = {candidate["patient_id"] for candidate in candidates}
    actual_ids = {str(patient_id) for patient_id in origin_support_by_patient}
    if actual_ids != expected_ids:
        raise ValueError("origin support must cover exactly the selected patient cohort")
    supports = {}
    for patient_id, support in origin_support_by_patient.items():
        if not isinstance(support, OriginSupport) or not support.source.strip():
            raise ValueError("each support entry must be an OriginSupport with a source description")
        mask = np.asarray(support.mask, dtype=bool)
        affine = np.asarray(support.affine, dtype=float)
        shape = reference_headers[patient_id]["shape"]
        if mask.ndim != 3 or mask.shape != shape or affine.shape != (4, 4) or not np.isfinite(affine).all():
            raise ValueError(f"origin-only support geometry is invalid for {patient_id}")
        if not np.allclose(affine, reference_headers[patient_id]["affine_2mm"], rtol=0.0, atol=GRID_TOLERANCE_MM):
            raise ValueError(f"origin-only support affine does not match the 2-mm reference for {patient_id}")
        if not mask.any():
            raise ValueError(f"origin-only support is empty for {patient_id}")
        frozen_mask = mask.copy()
        frozen_mask.setflags(write=False)
        supports[patient_id] = frozen_mask
    return supports


def _target_grid(header: dict) -> tuple[tuple[int, int, int], np.ndarray]:
    shape_native = np.asarray(header["shape"], dtype=int)
    shape_2mm = tuple(int(value) for value in ((shape_native + 1) // 2))
    affine_2mm = header["affine"].copy()
    affine_2mm[:3, :3] *= 2.0
    return shape_2mm, affine_2mm


def _selected_mask_hash_manifest(candidates: list[dict], selected_paths: dict[str, tuple[list[Path], Path]],
                                 data_root: Path) -> tuple[dict[Path, str], str]:
    """Hash every selected mask file and return only a cohort aggregate hash to callers."""
    entries = []
    hashes: dict[Path, str] = {}
    for candidate in candidates:
        patient_id = candidate["patient_id"]
        mask_paths, _ = selected_paths[patient_id]
        for scan_number, path in zip(candidate["scan_numbers"], mask_paths):
            digest = _sha256_file(path)
            relative_path = path.resolve().relative_to(data_root.resolve()).as_posix()
            entries.append({"path": relative_path, "scan_number": int(scan_number), "sha256": digest})
            hashes[path.resolve()] = digest
    entries.sort(key=lambda item: (item["path"], item["scan_number"]))
    payload = json.dumps(entries, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashes, hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_mu_first_two_to_third(
    *,
    evidence_path: Path,
    origin_support_by_patient: Mapping[str, OriginSupport] | None,
    scan_manifest_path: Path = DEFAULT_SCAN_MANIFEST,
    patient_manifest_path: Path = DEFAULT_PATIENT_MANIFEST,
    data_root: Path = DEFAULT_DATA_ROOT,
    metadata_report_path: Path = DEFAULT_METADATA_REPORT,
    mri_qc_report_path: Path = DEFAULT_MRI_QC_REPORT,
    allowed_evidence_roots: Iterable[Path] | None = None,
) -> MUDevelopmentCohort:
    """Load 18 reviewed MU GBM first-two-to-third histories only after metadata and evidence gates pass.

    Metadata selection is exactly: MU patients with ``is_gbm=true``; first three rows with
    ``mask_exists=true`` ordered by numeric scan number; three unique scans with finite strictly increasing days;
    second-to-third gap in [76,104] days; all three selected mask files plus second-visit T1n file present; and
    matching 1-mm millimeter headers. Target mask values are not read until all of these checks, snapshot checks,
    caller-domain checks, and the source-backed registration/target-support evidence gate have passed.
    """
    scan_path = Path(scan_manifest_path).resolve()
    patient_path = Path(patient_manifest_path).resolve()
    root = Path(data_root).resolve()
    metadata_path = Path(metadata_report_path).resolve()
    mri_path = Path(mri_qc_report_path).resolve()
    evidence_path = Path(evidence_path).resolve()

    snapshot_hashes = {"scan_manifest": _sha256_file(scan_path),
                       "patient_manifest": _sha256_file(patient_path),
                       "metadata_screen": _sha256_file(metadata_path),
                       "mri_only_qc": _sha256_file(mri_path)}
    metadata_report, mri_report = _validate_snapshot(scan_path, patient_path, metadata_path, mri_path)

    scan_selection_fields = {"patient_id", "dataset", "scan_number", "day", "mask_exists"}
    patient_selection_fields = {"patient_id", "dataset", "is_gbm"}
    scan_rows = _read_csv(scan_path, scan_selection_fields)
    patient_rows_selection = _read_csv(patient_path, patient_selection_fields)
    patient_rows_all = _read_csv(patient_path, patient_selection_fields | {"fold"})
    projected_hashes = {
        "scan_selection_columns": _projected_rows_sha256(scan_rows),
        "patient_selection_columns": _projected_rows_sha256(patient_rows_selection),
    }
    if projected_hashes != EXPECTED_PROJECTED_MANIFEST_SHA256:
        raise ValueError("current manifests do not reproduce the committed metadata-screen selection hashes")

    candidates, selection_counts = _metadata_candidates(scan_rows, patient_rows_all)
    if len(candidates) != EXPECTED_PATIENTS:
        raise ValueError(f"metadata selection found {len(candidates)} patients; reviewed expected count is 18")
    selected_ids = [candidate["patient_id"] for candidate in candidates]

    # Header-only preflight for every selected mask and the forecast-origin MRI. No mask voxel values are touched.
    file_counts = Counter()
    reference_headers: dict[str, dict] = {}
    selected_paths: dict[str, tuple[list[Path], Path]] = {}
    for candidate in candidates:
        patient_id = candidate["patient_id"]
        mask_paths, origin_brain_path = _paths(patient_id, candidate["scan_numbers"], root)
        selected_paths[patient_id] = (mask_paths, origin_brain_path)
        if not all(path.is_file() for path in mask_paths) or not origin_brain_path.is_file():
            raise FileNotFoundError(f"selected mask or origin T1n file is missing for a reviewed candidate")
        headers = [_header(path, root) for path in mask_paths]
        origin_header = _header(origin_brain_path, root)
        ref = headers[0]
        aligned = all(
            header["shape"] == ref["shape"]
            and np.allclose(header["affine"], ref["affine"], rtol=0.0, atol=GRID_TOLERANCE_MM)
            and np.allclose(header["spacing_mm"], ref["spacing_mm"], rtol=0.0, atol=GRID_TOLERANCE_MM)
            for header in headers[1:] + [origin_header]
        )
        if not aligned:
            raise ValueError(f"selected MU mask/T1n headers are not aligned for {patient_id}")
        shape_2mm, affine_2mm = _target_grid(ref)
        reference_headers[patient_id] = {"shape": shape_2mm, "affine_2mm": affine_2mm,
                                         "native_shape": ref["shape"], "native_affine": ref["affine"]}
        file_counts["mask_files_header_checked"] += len(mask_paths)
        file_counts["origin_t1n_headers_checked"] += 1

    if file_counts["mask_files_header_checked"] != 3 * EXPECTED_PATIENTS:
        raise ValueError("selected mask header count does not match the reviewed cohort")
    if file_counts["origin_t1n_headers_checked"] != EXPECTED_PATIENTS:
        raise ValueError("origin T1n header count does not match the reviewed cohort")

    # This gate precedes all reads of tumor-mask voxel arrays.
    if allowed_evidence_roots is None:
        raise PermissionError("provide explicit allowed_evidence_roots for local registration/support evidence")
    evidence_roots = tuple(Path(path) for path in allowed_evidence_roots)
    evidence = _validate_evidence_gate(evidence_path, selected_ids, snapshot_hashes, evidence_roots)
    origin_domains = _verify_origin_supports(origin_support_by_patient, candidates, reference_headers)
    # The evidence gate currently permits only source-attested full-grid coverage. Represent that
    # explicitly for scoring; partial or unknown target support is not silently treated as negative.
    target_valid_support = {
        candidate["patient_id"]: np.ones(reference_headers[candidate["patient_id"]]["shape"], dtype=bool)
        for candidate in candidates
    }
    for support_mask in target_valid_support.values():
        support_mask.setflags(write=False)

    # Hash full selected mask files immediately before any mask voxel array is opened.
    mask_hashes_before, mask_manifest_sha256 = _selected_mask_hash_manifest(candidates, selected_paths, root)

    import nibabel as nib
    from nibabel.processing import resample_from_to

    records = []
    loaded_counts = Counter()
    try:
        for candidate in candidates:
            patient_id = candidate["patient_id"]
            mask_paths, _ = selected_paths[patient_id]
            grid = reference_headers[patient_id]
            target = (grid["shape"], grid["affine_2mm"])
            masks_2mm = []
            for scan_index, path in enumerate(mask_paths):
                image = nib.load(str(path), mmap=True)
                # Extract only the binary core classes; other label classes are not used downstream.
                raw_core = np.isin(np.asarray(image.dataobj), (1, 3)).astype(np.uint8)
                core_image = nib.Nifti1Image(raw_core, np.asarray(image.affine, dtype=float))
                resampled = resample_from_to(core_image, target, order=0, mode="constant", cval=0)
                core_2mm = np.asarray(resampled.dataobj, dtype=np.uint8) > 0
                if core_2mm.shape != grid["shape"]:
                    raise ValueError(f"resampled core mask has wrong 2-mm shape for {patient_id}")
                core_2mm.setflags(write=False)
                masks_2mm.append(core_2mm)
                loaded_counts["mask_voxel_arrays_read"] += 1
                loaded_counts["mask_arrays_resampled_to_2mm"] += 1
                loaded_counts[f"visit_{scan_index + 1}_arrays_read"] += 1

            # Domain comes only from the caller's source-verified origin support. Never clip target labels.
            if np.any(masks_2mm[0] & ~origin_domains[patient_id]) or np.any(
                    masks_2mm[1] & ~origin_domains[patient_id]):
                raise ValueError(f"history core lies outside supplied origin-only support for {patient_id}")
            if np.any(masks_2mm[2] & ~origin_domains[patient_id]):
                raise ValueError(f"third-visit core lies outside origin-only support for {patient_id}")
            days_immutable = np.asarray(candidate["days"], dtype=float)
            days_immutable.setflags(write=False)
            records.append(PatientSeries(
                patient_id=patient_id,
                days=days_immutable,
                masks=tuple(masks_2mm),
            ))
    finally:
        # Detect mutation during any load, including a partial cohort load that raises.
        mask_hashes_after, post_manifest_sha256 = _selected_mask_hash_manifest(candidates, selected_paths, root)
        if mask_hashes_after != mask_hashes_before or post_manifest_sha256 != mask_manifest_sha256:
            raise RuntimeError("a selected MU mask file changed while the cohort was loading")

    folds = MappingProxyType({candidate["patient_id"]: int(candidate["fold"]) for candidate in candidates})
    audit = {
        "schema_version": 1,
        "study": "MU exploratory first-two-to-third PDE development cohort",
        "selection": {
            "rule": "is_gbm=true; first three mask-available visits by numeric scan number; finite strictly increasing days; second-to-third gap 76-104 days",
            "candidate_patients": len(candidates),
            "selection_counts": dict(sorted(selection_counts.items())),
            "target_mask_values_used_for_selection": False,
            "target_masks_used_for_domain": False,
            "target_visit_mri_used_for_selection": False,
            "target_support_coverage_scope": evidence["document"]["target_support"]["coverage_scope"],
            "patient_fold_source": "existing patient_manifest.csv fold column",
        },
        "source_hashes": {
            **snapshot_hashes,
            **projected_hashes,
            "selected_patient_ids_sha256_sorted_nul_joined": _selected_id_hash(selected_ids),
            "registration_evidence_sha256": evidence["registration_evidence_sha256"],
            "target_support_evidence_sha256": evidence["target_support_evidence_sha256"],
            "evidence_document_sha256": evidence["evidence_document_sha256"],
            "selected_mask_file_hash_manifest_sha256": mask_manifest_sha256,
        },
        "header_preflight": {
            "selected_mask_headers_checked": int(file_counts["mask_files_header_checked"]),
            "origin_t1n_headers_checked": int(file_counts["origin_t1n_headers_checked"]),
            "spatial_unit": "mm",
            "native_spacing_mm": [1.0, 1.0, 1.0],
            "2mm_resampling": "nearest neighbor, shared origin-reference grid per patient",
            "affine_tolerance_mm": GRID_TOLERANCE_MM,
        },
        "access": {
            "mask_voxel_arrays_read": int(loaded_counts["mask_voxel_arrays_read"]),
            "mask_arrays_resampled_to_2mm": int(loaded_counts["mask_arrays_resampled_to_2mm"]),
            "visit_1_mask_arrays_read": int(loaded_counts["visit_1_arrays_read"]),
            "visit_2_mask_arrays_read": int(loaded_counts["visit_2_arrays_read"]),
            "visit_3_target_mask_arrays_read": int(loaded_counts["visit_3_arrays_read"]),
            "selected_mask_file_hash_count": len(mask_hashes_before),
            "selected_mask_files_rehashed_after_load": True,
            "label_values_used": [1, 3],
            "origin_support_source": "caller-supplied, exact 2-mm origin-only grid",
            "target_valid_support_source": "source-attested entire image grid; explicit all-true scoring mask",
            "target_valid_support_voxels_per_patient": int(next(iter(target_valid_support.values())).sum()),
            "source_backed_registration_gate_passed": True,
            "source_backed_target_support_gate_passed": True,
        },
        "interpretation_limit": (
            "MU was previously explored. This loader does not create independent validation. The caller's evidence "
            "must justify registration and target support for the selected cohort; the earlier MRI-only screen is "
            "descriptive and does not itself clear either gate."
        ),
    }
    cohort = MUDevelopmentCohort(
        tuple(records), folds, MappingProxyType(origin_domains),
        MappingProxyType(target_valid_support), audit,
    )
    _LOADER_ISSUED_COHORTS.add(cohort)
    return cohort
