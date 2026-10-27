#!/usr/bin/env python3
"""Locked external patient test for next-visit annotated GBM-core additions.

MU-Glioma-Post is development data only. UCSD-PTGBM v3 is the external test.
One first chronological adjacent pair is used per unique GBM patient. The test
target is a change in source-provided labels, not biological growth, recurrence,
or treatment response. The script stores aggregate results only and refuses to
overwrite a prior locked-test output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import platform
import re
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import GroupKFold


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_UCSD_ROOT = ROOT / "data/external/UCSD-PTGBM-v3"
DEFAULT_UCSD_CLINICAL = ROOT / "research_notes/trackB_confirmatory_audit_ucsd/UCSD_PTGBM-clinical-information_v3_2026-12-Mar.xlsx"
DEFAULT_MU_ROOT = ROOT / "data/tcia/MU-Glioma-Post"
DEFAULT_PATIENT_MANIFEST = ROOT / "data/manifests/patient_manifest.csv"
DEFAULT_PAIR_MANIFEST = ROOT / "data/manifests/forecast_pairs_mu.csv"
DEFAULT_PROTOCOL = ROOT / "paper/track_b_ucsd_external_test_protocol.md"
DEFAULT_OUTPUT = ROOT / "output/track_b_ucsd_locked_test/aggregate_result.json"
EXPECTED_UCSD_CLINICAL_SHA256 = "bbf40bbdd289266924869052806de4d7e8fc9ee73286cdc0cc066672ef5bb51d"
LOCKED_INPUT_PATHS = (
    "paper/track_b_ucsd_external_test_protocol.md",
    "src/160_growth_zone_features.py",
    "src/161_growth_zone_model.py",
    "src/162_growth_zone_cv.py",
    "src/track_b_run_provenance.py",
    "src/164_signed_spatial_heads.py",
    "src/166_trackb_ucsd_external_test.py",
)


def _load_local_module(alias: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(alias, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[alias] = module
    spec.loader.exec_module(module)
    return module


BASE = _load_local_module("trackb_ucsd_base_cv", ROOT / "src/162_growth_zone_cv.py")
FEATURES = BASE.FEATURES
HEADS = _load_local_module("trackb_ucsd_signed_heads", ROOT / "src/164_signed_spatial_heads.py")

RIDGE = 1.0
PRIOR_PENALTY = 1e-3
MAX_PER_CLASS_PER_PATIENT = 4000
FIT_SEED = 20261026
N_DEV_FOLDS = 5
N_BOOTSTRAPS = 10_000
MIN_EFFECT_AUC = 0.02
GATE_STRENGTH = 1.0
BASELINE_TIE_ORDER = ("ungated_hybrid", "mri_only", "pde_only", "distance_only", "no_change")


@dataclass(frozen=True)
class ClinicalVisit:
    patient_id: str
    timepoint: int
    timepoint_code: str
    brats_id: str
    days_from_anchor: float | None


@dataclass(frozen=True)
class ScanFiles:
    paths: Mapping[str, Path]


@dataclass(frozen=True)
class PatientPair:
    patient_id: str
    first: ClinicalVisit
    second: ClinicalVisit
    first_files: ScanFiles
    second_files: ScanFiles
    interval_days: float | None


@dataclass(frozen=True)
class PatientSample:
    patient_key: str
    x: np.ndarray
    prior: np.ndarray
    instability: np.ndarray
    target: np.ndarray
    risk_voxels: int
    additions_total: int
    additions_outside_risk: int
    additions_outside_current_support: int
    interval_days: float | None


def _truth(value: Any) -> bool:
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def _safe_float(value: Any) -> float | None:
    if value is None or str(value).strip().casefold() in {"", "unknown", "na", "n/a", "none"}:
        return None
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if np.isfinite(result) else None


def read_ucsd_gbm_visits(clinical_path: Path) -> list[ClinicalVisit]:
    """Read current v3 clinical rows; never infer diagnosis from MRI labels."""
    try:
        import openpyxl
    except ImportError as exc:
        raise RuntimeError("openpyxl is required to read the TCIA v3 clinical workbook") from exc
    workbook = openpyxl.load_workbook(clinical_path, read_only=True, data_only=True)
    try:
        sheet = next((s for s in workbook.worksheets if "clinical" in s.title.casefold() or "cinical" in s.title.casefold()), None)
        if sheet is None:
            raise ValueError("UCSD v3 workbook has no clinical information worksheet")
        rows = sheet.iter_rows(values_only=True)
        header = next(rows, None)
        if header is None:
            raise ValueError("UCSD clinical worksheet is empty")
        columns = {str(value).strip(): index for index, value in enumerate(header) if value is not None}
        required = {"ID", "BraTS ID", "Primary Diagnosis"}
        if not required.issubset(columns):
            raise ValueError(f"UCSD v3 clinical worksheet lacks required columns: {sorted(required - columns.keys())}")
        anchor_column = next(
            (name for name in columns if name.startswith("Days from Acquisition to Date of initial surgery")),
            None,
        )
        result: list[ClinicalVisit] = []
        seen_ids: set[str] = set()
        seen_patient_timepoints: set[tuple[str, int]] = set()
        for row in rows:
            if str(row[columns["Primary Diagnosis"]]).strip().casefold() != "glioblastoma":
                continue
            tcia_id = str(row[columns["ID"]]).strip()
            match = re.fullmatch(r"(.+?)_(\d+)", tcia_id)
            if not match:
                raise ValueError("UCSD v3 clinical ID does not end in a numeric timepoint")
            patient_id, timepoint = match.group(1), int(match.group(2))
            if tcia_id.casefold() in seen_ids or (patient_id, timepoint) in seen_patient_timepoints:
                raise ValueError("UCSD v3 clinical workbook contains a duplicate patient/timepoint row")
            seen_ids.add(tcia_id.casefold())
            seen_patient_timepoints.add((patient_id, timepoint))
            brats_id = str(row[columns["BraTS ID"]]).strip()
            if not brats_id or brats_id.casefold() in {"none", "unknown"}:
                raise ValueError("UCSD v3 clinical row lacks its BraTS subject/timepoint crosswalk")
            anchor = _safe_float(row[columns[anchor_column]]) if anchor_column else None
            result.append(ClinicalVisit(patient_id, timepoint, match.group(2), brats_id, anchor))
    finally:
        workbook.close()
    if not result:
        raise ValueError("no exact-primary-diagnosis Glioblastoma rows found in UCSD v3")
    return result


def _nifti_files(root: Path) -> list[Path]:
    if not root.is_dir():
        raise NotADirectoryError(f"UCSD image root does not exist: {root}")
    return sorted(
        path for path in root.rglob("*")
        if path.is_file() and (path.name.casefold().endswith(".nii") or path.name.casefold().endswith(".nii.gz"))
    )


def _nifti_stem(path: Path) -> str:
    name = path.name.casefold()
    return name[:-7] if name.endswith(".nii.gz") else name[:-4]


def _modality(path: Path) -> str | None:
    stem = _nifti_stem(path)
    tokens = set(re.split(r"[^a-z0-9]+", stem))
    if tokens & {"t1n", "t1pre", "t1w"}:
        return "t1n"
    if tokens & {"t1c", "t1ce", "t1gd", "t1post"}:
        return "t1c"
    if tokens & {"t2", "t2w"}:
        return "t2w"
    if tokens & {"t2f", "flair"}:
        return "t2f"
    # UCSD-PTGBM also supplies cellular-tumor segmentations whose labels use
    # information across visits. The frozen endpoint requires the standard
    # BraTS compartment mask only.
    if stem.endswith("_brats_tumor_seg") or stem == "brats_tumor_seg":
        return "seg"
    return None


def _identifier_matches_path(identifier: str, path: Path) -> bool:
    """Match a complete ID token so visit 1 cannot collide with visit 10."""
    pattern = rf"(?<![a-z0-9]){re.escape(identifier.casefold())}(?![a-z0-9])"
    return re.search(pattern, str(path).casefold()) is not None


def _files_for_visit(files: Sequence[Path], visit: ClinicalVisit) -> ScanFiles | None:
    identifiers = (f"{visit.patient_id}_{visit.timepoint_code}", visit.brats_id)
    candidates = [
        path for path in files
        if any(_identifier_matches_path(identifier, path) for identifier in identifiers)
    ]
    found: dict[str, list[Path]] = defaultdict(list)
    for path in candidates:
        kind = _modality(path)
        if kind:
            found[kind].append(path)
    expected = (*FEATURES.MODALITIES, "seg")
    if any(len(found.get(kind, [])) != 1 for kind in expected):
        return None
    return ScanFiles({kind: found[kind][0] for kind in expected})


def _load_labels(path: Path, image: Any | None = None) -> np.ndarray:
    import nibabel as nib

    if image is None:
        image = nib.load(str(path))
    values = np.asanyarray(image.dataobj)
    if values.ndim != 3 or not np.all(np.isfinite(values)):
        raise ValueError(f"segmentation must be a finite 3-D label image: {path.name}")
    rounded = np.rint(values)
    if not np.allclose(values, rounded, rtol=0.0, atol=1e-5):
        raise ValueError(f"segmentation contains noninteger label values: {path.name}")
    if np.any((rounded < 0) | (rounded > 4)):
        raise ValueError(f"segmentation contains labels outside the documented 0-4 range: {path.name}")
    return rounded.astype(np.int16, copy=False)


def _valid_orthogonal_affine(affine: np.ndarray) -> bool:
    """Accept only orthogonal 1 mm grids used by the frozen perturbation scale."""
    value = np.asarray(affine, dtype=np.float64)
    if value.shape != (4, 4) or not np.all(np.isfinite(value)):
        return False
    spacing = np.linalg.norm(value[:3, :3], axis=0)
    if np.any(spacing <= 0) or not np.allclose(spacing, 1.0, rtol=0.0, atol=1e-3):
        return False
    directions = value[:3, :3] / spacing[np.newaxis, :]
    return bool(np.allclose(directions.T @ directions, np.eye(3), rtol=0.0, atol=1e-5))


def _image_grid_valid(image: Any) -> bool:
    units = image.header.get_xyzt_units()[0]
    return units == "mm" and _valid_orthogonal_affine(np.asarray(image.affine, dtype=np.float64))


def _validate_one_mm_file_headers(paths: Sequence[Path], cohort: str) -> None:
    import nibabel as nib

    for path in paths:
        image = nib.load(str(path))
        if not _image_grid_valid(image):
            raise ValueError(f"{cohort} input is not on an orthogonal 1 mm grid with millimetre units: {path.name}")


def _header_aligned(scan: ScanFiles) -> bool:
    import nibabel as nib

    images = [nib.load(str(scan.paths[name])) for name in (*FEATURES.MODALITIES, "seg")]
    shape = images[-1].shape
    affine = np.asarray(images[-1].affine, dtype=np.float64)
    return _image_grid_valid(images[-1]) and all(
        image.shape == shape and np.allclose(np.asarray(image.affine), affine, rtol=0.0, atol=1e-3)
        and _image_grid_valid(image)
        for image in images[:-1]
    )


def _visits_cross_aligned(first: ScanFiles, second: ScanFiles) -> bool:
    import nibabel as nib

    images = [
        nib.load(str(scan.paths[name]))
        for scan in (first, second)
        for name in (*FEATURES.MODALITIES, "seg")
    ]
    shape = images[0].shape
    affine = np.asarray(images[0].affine, dtype=np.float64)
    return _image_grid_valid(images[0]) and all(
        image.shape == shape and np.allclose(np.asarray(image.affine), affine, rtol=0.0, atol=1e-3)
        and _image_grid_valid(image)
        for image in images[1:]
    )


def select_ucsd_pairs(clinical_path: Path, data_root: Path) -> tuple[list[PatientPair], dict[str, int]]:
    """Use the first chronological adjacent pair with complete, aligned files."""
    visits = read_ucsd_gbm_visits(clinical_path)
    files = _nifti_files(data_root)
    by_patient: dict[str, list[ClinicalVisit]] = defaultdict(list)
    for visit in visits:
        by_patient[visit.patient_id].append(visit)
    summary: Counter[str] = Counter()
    summary["gbm_clinical_visit_rows"] = len(visits)
    summary["gbm_patients_in_clinical_metadata"] = len(by_patient)
    selected: list[PatientPair] = []
    for patient_id, patient_visits in sorted(by_patient.items()):
        patient_visits.sort(key=lambda row: row.timepoint)
        if len(patient_visits) < 2:
            summary["gbm_without_repeat"] += 1
            continue
        summary["gbm_patients_with_repeat"] += 1
        adjacent_pairs = list(zip(patient_visits, patient_visits[1:]))
        for first, second in adjacent_pairs:
            summary["adjacent_chronology_pairs_checked"] += 1
            if first.days_from_anchor is not None and second.days_from_anchor is not None:
                if second.days_from_anchor >= first.days_from_anchor:
                    raise ValueError("UCSD timepoint suffix conflicts with recorded acquisition-to-anchor chronology")
        selected_pair = None
        for first, second in adjacent_pairs:
            summary["adjacent_pairs_examined"] += 1
            first_files = _files_for_visit(files, first)
            second_files = _files_for_visit(files, second)
            if first_files is None or second_files is None:
                summary["adjacent_pairs_incomplete_files"] += 1
                continue
            if not _header_aligned(first_files) or not _header_aligned(second_files):
                summary["adjacent_pairs_within_visit_grid_failure"] += 1
                continue
            if not _visits_cross_aligned(first_files, second_files):
                summary["adjacent_pairs_between_visit_grid_failure"] += 1
                continue
            interval = (
                abs(second.days_from_anchor - first.days_from_anchor)
                if first.days_from_anchor is not None and second.days_from_anchor is not None
                else None
            )
            selected_pair = PatientPair(patient_id, first, second, first_files, second_files, interval)
            break
        if selected_pair is None:
            summary["gbm_no_complete_aligned_adjacent_pair"] += 1
            continue
        selected.append(selected_pair)
        summary["selected_gbm_pairs"] += 1
        if selected_pair.interval_days is None:
            summary["selected_interval_missing"] += 1
    return selected, dict(summary)


def _load_ucsd_visit(scan: ScanFiles) -> Any:
    import nibabel as nib

    images: dict[str, np.ndarray] = {}
    reference = nib.load(str(scan.paths["seg"]))
    labels = _load_labels(scan.paths["seg"], reference)
    affine = np.asarray(reference.affine, dtype=np.float64)
    for modality in FEATURES.MODALITIES:
        image = nib.load(str(scan.paths[modality]))
        if image.shape != labels.shape or not np.allclose(image.affine, affine, rtol=0.0, atol=1e-3):
            raise ValueError("UCSD modality is not aligned with its source segmentation")
        images[modality] = np.asarray(image.dataobj, dtype=np.float32)
    return FEATURES.VisitVolume(images=images, labels=labels, affine=affine)


def _patient_sample(
    patient_key: str,
    current: Any,
    future_labels: np.ndarray,
    future_affine: np.ndarray,
    parameters: Any,
    interval_days: float | None,
) -> PatientSample:
    FEATURES.validate_aligned_visits(
        current,
        FEATURES.VisitVolume(images=current.images, labels=future_labels, affine=future_affine),
    )
    visit = FEATURES.build_origin_features(current, parameters)
    log_draws = FEATURES.fixed_origin_perturbation_draws(current, visit, parameters)
    instability_map = FEATURES.origin_perturbation_instability(
        visit.log_pde_field, log_draws, visit.candidate
    )
    x, _ = FEATURES.origin_feature_matrix(visit)
    future_core = FEATURES.annotated_core_mask(future_labels)
    additions_all = future_core & ~visit.source
    target_full = additions_all & visit.candidate
    additions_total = int(np.count_nonzero(additions_all))
    additions_inside = int(np.count_nonzero(target_full))
    additions_outside_current_support = int(np.count_nonzero(additions_all & ~FEATURES.brain_support(current.images, current.labels)))
    return PatientSample(
        patient_key=patient_key,
        x=x,
        prior=visit.log_pde_field[visit.candidate].astype(np.float32, copy=False),
        instability=instability_map[visit.candidate].astype(np.float32, copy=False),
        target=target_full[visit.candidate].astype(np.uint8, copy=False),
        risk_voxels=int(np.count_nonzero(visit.candidate)),
        additions_total=additions_total,
        additions_outside_risk=max(0, additions_total - additions_inside),
        additions_outside_current_support=additions_outside_current_support,
        interval_days=interval_days,
    )


def load_mu_development_samples(
    data_root: Path,
    patient_manifest: Path,
    pair_manifest: Path,
    parameters: Any,
) -> tuple[list[PatientSample], list[Path]]:
    patients = {row["patient_id"]: row for row in BASE._read_csv(patient_manifest, BASE.PATIENT_MANIFEST_COLUMNS)}
    pairs = BASE._read_csv(pair_manifest, BASE.PAIR_MANIFEST_COLUMNS)
    selected: dict[str, dict[str, str]] = {}
    for pair in pairs:
        patient_id = pair["patient_id"]
        patient = patients.get(patient_id)
        # This repository stores MU-Glioma-Post as "MU" in patient_manifest.csv.
        if patient is None or str(patient.get("dataset", "")).strip().casefold() not in {"mu", "mu-glioma-post"}:
            continue
        if not (BASE._truth(patient["eligible"]) and BASE._truth(patient["is_gbm"])):
            continue
        if not (BASE._truth(pair["in_primary"]) and BASE._truth(pair["input_core_nonempty"])):
            continue
        old = selected.get(patient_id)
        sort_key = (BASE._finite_number(pair["day_in"], "day_in"), int(pair["tp_in"]), int(pair["tp_out"]))
        if old is None:
            selected[patient_id] = pair
        else:
            old_key = (BASE._finite_number(old["day_in"], "day_in"), int(old["tp_in"]), int(old["tp_out"]))
            if sort_key < old_key:
                selected[patient_id] = pair
    samples: list[PatientSample] = []
    source_paths: list[Path] = []
    for index, (patient_id, pair) in enumerate(sorted(selected.items()), start=1):
        timepoint_in, timepoint_out = int(pair["tp_in"]), int(pair["tp_out"])
        current_dir = data_root / patient_id / f"Timepoint_{timepoint_in}"
        future_dir = data_root / patient_id / f"Timepoint_{timepoint_out}"
        input_paths = [BASE._single_file(current_dir, "*_tumorMask.nii.gz")]
        input_paths.extend(
            BASE._single_file(current_dir, f"*_brain_{modality}.nii.gz")
            for modality in FEATURES.MODALITIES
        )
        input_paths.append(BASE._single_file(future_dir, "*_tumorMask.nii.gz"))
        source_paths.extend(input_paths)
        _validate_one_mm_file_headers(input_paths, "MU development")
        current = BASE._load_visit(data_root, patient_id, timepoint_in)
        future_labels, future_affine = BASE._load_mask(data_root, patient_id, timepoint_out)
        samples.append(
            _patient_sample(
                patient_id,
                current,
                future_labels,
                future_affine,
                parameters,
                BASE._finite_number(pair["dt_days"], "dt_days"),
            )
        )
        print(f"loaded development pair {index}/{len(selected)}", flush=True)
    if len(samples) < N_DEV_FOLDS:
        raise ValueError("too few unique MU development patients for the frozen grouped development folds")
    return samples, source_paths


def _patient_gate_scale(samples: Sequence[PatientSample]) -> float:
    values = [float(np.percentile(row.instability, 90.0)) for row in samples if row.instability.size]
    if not values:
        return 1.0
    scale = float(np.median(values))
    return scale if np.isfinite(scale) and scale > 1e-8 else 1.0


def _gate(sample: PatientSample, scale: float, *, active: bool = True) -> np.ndarray:
    if not active:
        return np.ones(sample.target.size, dtype=np.float32)
    return np.exp(-GATE_STRENGTH * sample.instability / scale).astype(np.float32)


def _head_samples(
    samples: Sequence[PatientSample],
    *,
    scale: float,
    use_prior: bool,
    use_gate: bool,
) -> list[Any]:
    result = []
    for row in samples:
        prior = row.prior if use_prior else np.zeros_like(row.prior)
        result.append(
            HEADS.SignedVoxelSample(
                patient_id=row.patient_key,
                features=row.x,
                spatial_prior_score=prior,
                stability_gate=_gate(row, scale, active=use_gate),
                target=row.target,
            )
        )
    return result


def _fit_head(samples: Sequence[PatientSample], *, scale: float, use_prior: bool, use_gate: bool, seed: int) -> Any:
    return HEADS.fit_signed_residual_head(
        _head_samples(samples, scale=scale, use_prior=use_prior, use_gate=use_gate),
        ridge=RIDGE,
        prior_penalty=PRIOR_PENALTY,
        max_per_class_per_interval=MAX_PER_CLASS_PER_PATIENT,
        seed=seed,
    )


def _scores(row: PatientSample, *, scale: float, models: Mapping[str, Any]) -> dict[str, np.ndarray]:
    gate = _gate(row, scale)
    ones = np.ones(row.target.size, dtype=np.float32)
    zeros = np.zeros(row.target.size, dtype=np.float32)
    return {
        "full_hybrid": models["full_hybrid"].predict_score(row.x, row.prior, gate),
        "ungated_hybrid": models["ungated_hybrid"].predict_score(row.x, row.prior, ones),
        "mri_only": models["mri_only"].predict_score(row.x, zeros, ones),
        "pde_only": row.prior.astype(np.float64),
        "distance_only": -row.x[:, 0].astype(np.float64),
        "no_change": np.zeros(row.target.size, dtype=np.float64),
    }


def _auc(target: np.ndarray, score: np.ndarray) -> float | None:
    if target.size == 0 or np.unique(target).size != 2:
        return None
    return float(roc_auc_score(target, score))


def _interval_stratum(interval_days: float | None) -> str:
    if interval_days is None:
        return "unknown"
    if interval_days <= 90.0:
        return "le_90_days"
    if interval_days <= 180.0:
        return "91_180_days"
    return "gt_180_days"


def develop_and_select_baseline(samples: Sequence[PatientSample]) -> tuple[dict[str, float], str]:
    groups = np.asarray([row.patient_key for row in samples])
    cv_scores: dict[str, list[float]] = {"full_hybrid": []}
    cv_scores.update({name: [] for name in BASELINE_TIE_ORDER})
    splitter = GroupKFold(n_splits=N_DEV_FOLDS)
    for fold, (train_idx, valid_idx) in enumerate(
        splitter.split(np.zeros(len(samples)), groups=groups), start=1
    ):
        train = [samples[index] for index in train_idx]
        valid = [samples[index] for index in valid_idx]
        scale = _patient_gate_scale(train)
        models = {
            "full_hybrid": _fit_head(train, scale=scale, use_prior=True, use_gate=True, seed=FIT_SEED + fold),
            "ungated_hybrid": _fit_head(train, scale=scale, use_prior=True, use_gate=False, seed=FIT_SEED + fold),
            "mri_only": _fit_head(train, scale=scale, use_prior=False, use_gate=False, seed=FIT_SEED + fold),
        }
        for row in valid:
            for name, score in _scores(row, scale=scale, models=models).items():
                value = _auc(row.target, score)
                if value is not None:
                    cv_scores[name].append(value)
    means = {name: float(np.mean(values)) for name, values in cv_scores.items() if values}
    missing = [name for name in (*BASELINE_TIE_ORDER, "full_hybrid") if name not in means]
    if missing:
        raise ValueError(f"development folds had no evaluable patient AUC for: {missing}")
    primary_comparator = max(BASELINE_TIE_ORDER, key=lambda name: means[name])
    return means, primary_comparator


def _hash_path_size_manifest(paths: Iterable[Path], root: Path) -> str:
    """Hash selected relative paths and sizes without reading image voxels."""
    digest = hashlib.sha256()
    rows = sorted({
        (str(path.resolve().relative_to(root.resolve())).replace("\\", "/"), path.stat().st_size)
        for path in paths
    })
    for relative, size in rows:
        digest.update(f"{relative}\t{size}\n".encode("utf-8"))
    return digest.hexdigest()


def _hash_files(paths: Iterable[Path], root: Path) -> str:
    """Hash selected file contents with relative names and sizes, but no IDs in output."""
    digest = hashlib.sha256()
    rows = sorted({
        (str(path.resolve().relative_to(root.resolve())).replace("\\", "/"), path)
        for path in paths
    })
    for relative, path in rows:
        before = path.stat()
        content_sha256 = _sha256(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError(f"source file changed while hashing: {relative}")
        digest.update(f"{relative}\t{before.st_size}\t".encode("utf-8"))
        digest.update(bytes.fromhex(content_sha256))
        digest.update(b"\n")
    return digest.hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _software_versions() -> dict[str, str]:
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for package in ("numpy", "scipy", "scikit-learn", "nibabel", "openpyxl"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "missing"
    return versions


def _git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, capture_output=True, text=True, check=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def _assert_locked_inputs_committed(protocol: Path) -> str:
    """Require the protocol and all scoring code to match a committed snapshot."""
    if protocol.resolve() != DEFAULT_PROTOCOL.resolve():
        raise ValueError("locked scoring requires the repository's frozen Track B protocol")
    status = subprocess.run(
        ["git", "status", "--porcelain", "--", *LOCKED_INPUT_PATHS],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=True,
    ).stdout.strip()
    if status:
        raise RuntimeError("locked protocol or scoring code is uncommitted; commit it before accessing UCSD outcomes")
    commit = _git_commit()
    if not commit:
        raise RuntimeError("locked scoring requires a Git commit for provenance")
    return commit


def _validate_locked_source_paths(args: argparse.Namespace) -> None:
    expected = (
        (args.ucsd_root, DEFAULT_UCSD_ROOT),
        (args.ucsd_clinical, DEFAULT_UCSD_CLINICAL),
        (args.mu_data_root, DEFAULT_MU_ROOT),
        (args.patient_manifest, DEFAULT_PATIENT_MANIFEST),
        (args.pair_manifest, DEFAULT_PAIR_MANIFEST),
        (args.protocol, DEFAULT_PROTOCOL),
    )
    if any(actual.resolve() != canonical.resolve() for actual, canonical in expected):
        raise ValueError("locked scoring only accepts the source roots and protocol frozen in the repository")
    clinical_sha256 = _sha256(args.ucsd_clinical)
    if clinical_sha256 != EXPECTED_UCSD_CLINICAL_SHA256:
        raise ValueError("UCSD clinical workbook checksum differs from the pinned v3 source file")


def _create_analysis_lock(
    *,
    protocol: Path,
    clinical_path: Path,
    patient_manifest: Path,
    pair_manifest: Path,
    ucsd_path_size_manifest_sha256: str,
    development_file_manifest_sha256: str,
) -> Path:
    """Create a one-shot marker immediately before reading UCSD label arrays."""
    if _sha256(clinical_path) != EXPECTED_UCSD_CLINICAL_SHA256:
        raise ValueError("UCSD clinical workbook changed before the locked run started")
    output_dir = DEFAULT_OUTPUT.parent
    output_dir.mkdir(parents=True, exist_ok=True)
    lock_path = output_dir / "analysis_started.json"
    record = {
        "status": "external_outcomes_access_started",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "git_commit": _assert_locked_inputs_committed(protocol),
        "locked_code_sha256": {
            relative: _sha256(ROOT / relative) for relative in LOCKED_INPUT_PATHS
        },
        "software_versions": _software_versions(),
        "protocol_sha256": _sha256(protocol),
        "runner_sha256": _sha256(Path(__file__)),
        "feature_module_sha256": _sha256(Path(FEATURES.__file__)),
        "development_cv_module_sha256": _sha256(Path(BASE.__file__)),
        "signed_head_module_sha256": _sha256(Path(HEADS.__file__)),
        "clinical_workbook_sha256": _sha256(clinical_path),
        "ucsd_selected_path_size_manifest_sha256": ucsd_path_size_manifest_sha256,
        "development_patient_manifest_sha256": _sha256(patient_manifest),
        "development_pair_manifest_sha256": _sha256(pair_manifest),
        "development_source_content_manifest_sha256": development_file_manifest_sha256,
        "note": "If scoring stops after this marker, inspect the failure before any amended rerun.",
    }
    try:
        with lock_path.open("x", encoding="utf-8") as handle:
            json.dump(record, handle, indent=2)
            handle.write("\n")
    except FileExistsError as exc:
        raise FileExistsError(f"external test already started; inspect {lock_path} before any rerun") from exc
    return lock_path


def _update_analysis_lock(lock_path: Path, **updates: Any) -> dict[str, Any]:
    """Atomically add post-lock input fingerprints to the one-shot record."""
    record = json.loads(lock_path.read_text(encoding="utf-8"))
    record.update(updates)
    temporary_path = lock_path.with_name(lock_path.name + ".tmp")
    with temporary_path.open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
    os.replace(temporary_path, lock_path)
    return record


def _preflight(args: argparse.Namespace) -> tuple[list[PatientPair], dict[str, int]]:
    pairs, counts = select_ucsd_pairs(args.ucsd_clinical, args.ucsd_root)
    if not pairs:
        raise ValueError("no GBM patient has a complete, aligned pair with all four standard modalities and masks")
    return pairs, counts


def _read_pair(pair: PatientPair) -> tuple[Any, np.ndarray, np.ndarray]:
    import nibabel as nib

    current = _load_ucsd_visit(pair.first_files)
    future_image = nib.load(str(pair.second_files.paths["seg"]))
    future_labels = _load_labels(pair.second_files.paths["seg"], future_image)
    return current, future_labels, np.asarray(future_image.affine, dtype=np.float64)


def _bootstrap_mean_ci(values: Sequence[float], *, seed: int) -> tuple[float, float]:
    if not values:
        return float("nan"), float("nan")
    array = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(seed)
    draws = rng.choice(array, size=(N_BOOTSTRAPS, array.size), replace=True).mean(axis=1)
    return float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def _added_value_criterion(delta: float, ci_low: float) -> bool:
    """Require evidence that the patient-average gain exceeds the preset margin."""
    return bool(np.isfinite(delta) and np.isfinite(ci_low) and ci_low > MIN_EFFECT_AUC)


def run_locked_test(args: argparse.Namespace) -> dict[str, Any]:
    if DEFAULT_OUTPUT.exists() or (DEFAULT_OUTPUT.parent / "analysis_started.json").exists():
        raise FileExistsError("locked UCSD analysis has already started; inspect its marker/output before rerunning")
    if not args.ucsd_clinical.is_file() or not args.protocol.is_file():
        raise FileNotFoundError("the current UCSD v3 clinical workbook and frozen protocol are required")
    if not args.patient_manifest.is_file() or not args.pair_manifest.is_file():
        raise FileNotFoundError("the frozen MU development manifests are required")
    _validate_locked_source_paths(args)
    pairs, preflight_counts = _preflight(args)
    ucsd_source_paths = [
        path for pair in pairs for bundle in (pair.first_files, pair.second_files) for path in bundle.paths.values()
    ]
    ucsd_path_size_manifest_sha256 = _hash_path_size_manifest(ucsd_source_paths, args.ucsd_root)
    parameters = FEATURES.OriginDiffusionParameters(
        edge_sensitivity=0.5,
        screening_length_mm=10.0,
        label2_log_permeability=0.69314718056,
        block_cavity=True,
    )

    print("loading MU development cohort", flush=True)
    development, development_source_paths = load_mu_development_samples(
        args.mu_data_root, args.patient_manifest, args.pair_manifest, parameters
    )
    development_cv, primary_comparator = develop_and_select_baseline(development)
    development_file_manifest_sha256 = _hash_files(development_source_paths, args.mu_data_root)
    gate_scale = _patient_gate_scale(development)
    final_models = {
        "full_hybrid": _fit_head(development, scale=gate_scale, use_prior=True, use_gate=True, seed=FIT_SEED),
        "ungated_hybrid": _fit_head(development, scale=gate_scale, use_prior=True, use_gate=False, seed=FIT_SEED),
        "mri_only": _fit_head(development, scale=gate_scale, use_prior=False, use_gate=False, seed=FIT_SEED),
    }

    # Preflight inspected only file names and NIfTI headers. Lock immediately
    # before reading the first UCSD segmentation voxel.
    analysis_lock = _create_analysis_lock(
        protocol=args.protocol,
        clinical_path=args.ucsd_clinical,
        patient_manifest=args.patient_manifest,
        pair_manifest=args.pair_manifest,
        ucsd_path_size_manifest_sha256=ucsd_path_size_manifest_sha256,
        development_file_manifest_sha256=development_file_manifest_sha256,
    )
    # Create the one-shot marker before hashing or loading UCSD array contents.
    # The content digest is added to the marker before any segmentation is read.
    ucsd_content_manifest_sha256 = _hash_files(ucsd_source_paths, args.ucsd_root)
    lock_record = _update_analysis_lock(
        analysis_lock,
        ucsd_selected_content_manifest_sha256=ucsd_content_manifest_sha256,
    )

    patient_aucs: dict[str, list[float]] = {"full_hybrid": []}
    patient_aucs.update({name: [] for name in BASELINE_TIE_ORDER})
    contrast: list[float] = []
    event_counts: list[int] = []
    outside_counts: list[int] = []
    outside_support_counts: list[int] = []
    intervals: list[float] = []
    interval_stratum_differences: dict[str, list[float]] = {
        "le_90_days": [],
        "91_180_days": [],
        "gt_180_days": [],
        "unknown": [],
    }
    zero_addition_pairs = 0
    print("scoring locked UCSD v3 patient pairs", flush=True)
    for index, pair in enumerate(pairs, start=1):
        current, future_labels, future_affine = _read_pair(pair)
        sample = _patient_sample(
            pair.patient_id,
            current,
            future_labels,
            future_affine,
            parameters,
            pair.interval_days,
        )
        if sample.target.sum() == 0:
            zero_addition_pairs += 1
        score_maps = _scores(sample, scale=gate_scale, models=final_models)
        current_aucs: dict[str, float] = {}
        for name, score in score_maps.items():
            value = _auc(sample.target, score)
            if value is not None:
                current_aucs[name] = value
                patient_aucs[name].append(value)
        if primary_comparator in current_aucs and "full_hybrid" in current_aucs:
            difference = current_aucs["full_hybrid"] - current_aucs[primary_comparator]
            contrast.append(difference)
            interval_stratum_differences[_interval_stratum(pair.interval_days)].append(difference)
        event_counts.append(sample.additions_total)
        outside_counts.append(sample.additions_outside_risk)
        outside_support_counts.append(sample.additions_outside_current_support)
        if pair.interval_days is not None:
            intervals.append(pair.interval_days)
        print(f"scored patient pair {index}/{len(pairs)}", flush=True)

    ci_low, ci_high = _bootstrap_mean_ci(contrast, seed=FIT_SEED + 100)
    if _hash_files(ucsd_source_paths, args.ucsd_root) != ucsd_content_manifest_sha256:
        raise RuntimeError("selected UCSD files changed during scoring; preserve the one-shot lock and investigate")
    if _sha256(args.ucsd_clinical) != lock_record["clinical_workbook_sha256"]:
        raise RuntimeError("UCSD clinical workbook changed during scoring; preserve the one-shot lock and investigate")
    if _sha256(args.patient_manifest) != lock_record["development_patient_manifest_sha256"]:
        raise RuntimeError("development patient manifest changed during scoring; preserve the one-shot lock and investigate")
    if _sha256(args.pair_manifest) != lock_record["development_pair_manifest_sha256"]:
        raise RuntimeError("development pair manifest changed during scoring; preserve the one-shot lock and investigate")
    if _hash_files(development_source_paths, args.mu_data_root) != development_file_manifest_sha256:
        raise RuntimeError("development source files changed during scoring; preserve the one-shot lock and investigate")
    delta = float(np.mean(contrast)) if contrast else float("nan")
    primary = {
        "contrast": f"full_hybrid minus development-selected {primary_comparator}",
        "n_evaluable_patients": len(contrast),
        "mean_patient_auc_difference": delta if np.isfinite(delta) else None,
        "paired_patient_bootstrap_95_percentile_ci": [ci_low, ci_high] if np.isfinite(ci_low) else None,
        "minimum_effect_auc": MIN_EFFECT_AUC,
        "criterion_met": _added_value_criterion(delta, ci_low),
    }
    means = {name: float(np.mean(values)) if values else None for name, values in patient_aucs.items()}
    result = {
        "status": "protocol_locked_external_patient_test_complete",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "interpretation": "next-visit change in source-provided NETC+ET annotations; not biological growth, recurrence, clinical progression, or treatment response",
        "source": {
            "dataset": "UCSD-PTGBM v3",
            "license": "CC BY 4.0 per current TCIA collection page",
            "dataset_doi": "10.7937/fwv2-dt74",
            "project_data_role": "external locked test; no UCSD outcomes used for model fitting or comparator selection",
            "clinical_workbook_sha256": lock_record["clinical_workbook_sha256"],
            "selected_source_path_size_manifest_sha256": ucsd_path_size_manifest_sha256,
            "selected_source_content_manifest_sha256": ucsd_content_manifest_sha256,
            "development_source_content_manifest_sha256": lock_record["development_source_content_manifest_sha256"],
            "development_patient_manifest_sha256": lock_record["development_patient_manifest_sha256"],
            "development_pair_manifest_sha256": lock_record["development_pair_manifest_sha256"],
        },
        "development": {
            "dataset": "MU-Glioma-Post",
            "n_patients": len(development),
            "patients_used_by_final_full_hybrid_fit": final_models["full_hybrid"].n_training_patients,
            "role": "previously explored development/training only; no confirmation claim",
            "patient_grouped_cv_auc_for_baseline_selection": development_cv,
            "primary_comparator_selected_without_UCSD_outcomes": primary_comparator,
            "gate_scale_from_development_only": gate_scale,
        },
        "cohort": {
            "metadata_and_file_eligible_unique_patients": len(pairs),
            "selected_pair_policy": "first chronological adjacent pair with all four standard modalities and standard segmentation, complete source registration, and common grid",
            "preflight_counts": preflight_counts,
            "auc_evaluable_patients": len(contrast),
            "pairs_with_no_in_risk_additions": zero_addition_pairs,
            "median_observed_interval_days": float(np.median(intervals)) if intervals else None,
            "min_observed_interval_days": float(np.min(intervals)) if intervals else None,
            "max_observed_interval_days": float(np.max(intervals)) if intervals else None,
            "pairs_with_unknown_interval": len(pairs) - len(intervals),
            "primary_difference_by_prespecified_interval_stratum": {
                name: {
                    "n_evaluable_patients": len(values),
                    "mean_patient_auc_difference": float(np.mean(values)) if values else None,
                }
                for name, values in interval_stratum_differences.items()
            },
            "total_annotated_additions": int(sum(event_counts)),
            "additions_outside_risk_set": int(sum(outside_counts)),
            "additions_outside_current_brain_support": int(sum(outside_support_counts)),
        },
        "primary_test": primary,
        "patient_macro_auc": means,
        "analysis": {
            "primary_metric": "mean of patient-specific ROC AUCs among patients with at least one in-risk addition and at least one in-risk non-addition voxel",
            "confidence_interval": "paired nonparametric bootstrap resampling patients, 10000 draws, percentile 95% interval",
            "training_ridge": RIDGE,
            "prior_penalty": PRIOR_PENALTY,
            "max_per_class_per_patient": MAX_PER_CLASS_PER_PATIENT,
            "fit_seed_base": FIT_SEED,
            "development_fold_seed_rule": "fit_seed_base plus one-based fold index",
            "final_fit_seed": FIT_SEED,
            "bootstrap_seed": FIT_SEED + 100,
            "bootstrap_draws": N_BOOTSTRAPS,
            "development_folds": N_DEV_FOLDS,
            "baseline_tie_order": list(BASELINE_TIE_ORDER),
            "gate_strength": GATE_STRENGTH,
            "diffusion_parameters": {
                "edge_sensitivity": parameters.edge_sensitivity,
                "screening_length_mm": parameters.screening_length_mm,
                "label2_log_permeability": parameters.label2_log_permeability,
                "block_cavity": parameters.block_cavity,
                "risk_set": "full current skull-stripped brain support, excluding current core and label-4 cavity",
            },
            "patient_ids_and_voxel_predictions_written": False,
            "one_shot_lock_file": str(analysis_lock.relative_to(ROOT)).replace("\\", "/"),
        },
        "provenance": {
            "git_commit_at_run": lock_record["git_commit"],
            "locked_code_sha256": lock_record["locked_code_sha256"],
            "protocol_sha256": lock_record["protocol_sha256"],
            "runner_sha256": lock_record["runner_sha256"],
            "feature_module_sha256": lock_record["feature_module_sha256"],
            "development_cv_module_sha256": lock_record["development_cv_module_sha256"],
            "signed_head_module_sha256": lock_record["signed_head_module_sha256"],
            "software_versions": lock_record["software_versions"],
            "ucsd_selected_content_manifest_sha256": ucsd_content_manifest_sha256,
        },
        "limitations": [
            "The endpoint measures change in supplied annotations and cannot separate biology, therapy, registration, or annotation effects.",
            "The next-visit interval varies and is not a predictor; results do not estimate a fixed-horizon forecast.",
            "The sensitivity gate is a deterministic boundary stress heuristic, not calibrated uncertainty.",
            "The primary AUC is conditional on an observed addition in the prespecified risk set; zero-addition pairs are reported separately.",
            "UCSD is publicly downloadable and methods may have seen the same challenge data; this is project-level held-out external evaluation, not a claim of novelty or clinical confirmation.",
        ],
    }
    DEFAULT_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with DEFAULT_OUTPUT.open("x", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, allow_nan=False)
        handle.write("\n")
    return result


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ucsd-root", type=Path, default=DEFAULT_UCSD_ROOT)
    parser.add_argument("--ucsd-clinical", type=Path, default=DEFAULT_UCSD_CLINICAL)
    parser.add_argument("--mu-data-root", type=Path, default=DEFAULT_MU_ROOT)
    parser.add_argument("--patient-manifest", type=Path, default=DEFAULT_PATIENT_MANIFEST)
    parser.add_argument("--pair-manifest", type=Path, default=DEFAULT_PAIR_MANIFEST)
    parser.add_argument("--protocol", type=Path, default=DEFAULT_PROTOCOL)
    parser.add_argument("--preflight-only", action="store_true", help="check source pair completeness and grid headers only")
    parser.add_argument("--score-locked-test", action="store_true", help="fit on MU and score UCSD once using the frozen protocol")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.preflight_only == args.score_locked_test:
        parser.error("choose exactly one: --preflight-only or --score-locked-test")
    if args.preflight_only:
        pairs, counts = _preflight(args)
        print(json.dumps({"ucsd_file_complete_pairs": len(pairs), "counts": counts}, indent=2))
        return 0
    result = run_locked_test(args)
    print(json.dumps({"status": result["status"], "n_ucsd_pairs": result["cohort"]["metadata_and_file_eligible_unique_patients"], "primary": result["primary_test"]}, indent=2))
    print(f"wrote aggregate locked-test result: {DEFAULT_OUTPUT}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
