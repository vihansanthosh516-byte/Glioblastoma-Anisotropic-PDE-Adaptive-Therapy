"""Patient-level bootstrap precision gate for the Track B external-test plan.

This estimates study resolution from development-cohort paired differences. It is not an effect estimate for
the external cohort and must receive exactly one paired difference per development patient.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from types import MappingProxyType
import weakref

import numpy as np

BOOTSTRAP_REPS = 10_000
BOOTSTRAP_SEED = 20261028
MAX_HALF_WIDTH_DICE = 0.05
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DEFAULT_ELIGIBILITY_REPORT_ROOT = PROJECT_ROOT / "output"

DUKE_ELIGIBILITY_RULE_V1 = (
    "Duke first-two-to-third patient cohort; use first three chronological eligible MRI/mask visits; "
    "require source-valid strictly increasing dates and immediate second-to-third gap 76-104 days inclusive; "
    "require source manifests, access and school-use approval, overlap, date, registration, support, and mask-provenance gates; "
    "select without reading target-mask voxel values; one eligible row per patient."
)
DUKE_ELIGIBILITY_RULE_SHA256 = hashlib.sha256(DUKE_ELIGIBILITY_RULE_V1.encode("utf-8")).hexdigest()


@dataclass(frozen=True, eq=False)
class ExternalEligibilityEvidence:
    """Validated pre-target-access eligibility report that supplies external N."""

    eligible_patient_count: int
    eligibility_rule_sha256: str
    report_path: str
    report_sha256: str
    selected_patient_ids_sha256_sorted_nul_joined: str
    source_manifest_sha256: Mapping[str, str]
    gate_evidence_sha256: Mapping[str, str]


_VALIDATED_ELIGIBILITY_REPORTS = weakref.WeakSet()


def is_validated_external_eligibility(value) -> bool:
    """Check that this exact evidence object came from report validation."""
    return isinstance(value, ExternalEligibilityEvidence) and value in _VALIDATED_ELIGIBILITY_REPORTS


def _valid_sha256(value) -> bool:
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value.lower()) is not None


def _inside_any_root(path: Path, roots: tuple[Path, ...]) -> bool:
    return any(path == root or root in path.parents for root in roots)


def _verify_report_artifact(
    report_path: Path, artifact, *, field: str, allowed_evidence_roots: tuple[Path, ...],
) -> tuple[Path, str]:
    if not isinstance(artifact, Mapping):
        raise ValueError(f"eligibility report requires a file-backed {field} record")
    raw_path = artifact.get("path")
    if not isinstance(raw_path, str) or not raw_path.strip():
        raise ValueError(f"eligibility report requires a path for {field}")
    candidate = Path(raw_path)
    if not candidate.is_absolute():
        candidate = report_path.parent / candidate
    path = candidate.resolve()
    if not _inside_any_root(path, allowed_evidence_roots):
        raise PermissionError(f"{field} file is outside explicitly allowed evidence roots")
    if not path.is_file():
        raise FileNotFoundError(path)
    expected_hash = str(artifact.get("sha256", "")).strip().lower()
    if not _valid_sha256(expected_hash):
        raise ValueError(f"eligibility report requires a valid {field}.sha256")
    actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_hash != expected_hash:
        raise PermissionError(f"{field} file hash does not match the eligibility report")
    return path, actual_hash


def _verify_reviewed_gate(
    report_path: Path, section, *, field: str, allowed_evidence_roots: tuple[Path, ...],
    status_required: bool = True,
) -> str:
    if not isinstance(section, Mapping):
        raise ValueError(f"{field} gate requires a reviewed evidence record")
    if status_required and section.get("status") != "passed":
        raise PermissionError(f"{field} gate is not recorded as passed")
    reviewer = section.get("reviewed_by")
    reviewed_at = section.get("reviewed_at")
    if not isinstance(reviewer, str) or not reviewer.strip():
        raise ValueError(f"{field} gate requires a named reviewer")
    if not isinstance(reviewed_at, str) or not reviewed_at.strip():
        raise ValueError(f"{field} gate requires a review timestamp")
    _, digest = _verify_report_artifact(
        report_path, section, field=field, allowed_evidence_roots=allowed_evidence_roots,
    )
    return digest


def validate_external_eligibility_report(
    report_path: str | Path, *, expected_report_sha256: str,
    allowed_report_roots: Iterable[str | Path] = (DEFAULT_ELIGIBILITY_REPORT_ROOT,),
    allowed_evidence_roots: Iterable[str | Path] = (PROJECT_ROOT,),
) -> ExternalEligibilityEvidence:
    """Validate the locked metadata-only Duke N report before precision can be calculated.

    The report must show access approval, completed source-backed cohort gates, and a closed target-mask-value
    state. The report, source manifests, selected-ID file, and reviewed gate evidence are hash-bound. A human must
    review the referenced evidence; this function verifies file identity and required report fields, not meaning.
    """
    path = Path(report_path).resolve()
    roots = tuple(Path(root).resolve() for root in allowed_report_roots)
    evidence_roots = tuple(Path(root).resolve() for root in allowed_evidence_roots)
    if not roots or not _inside_any_root(path, roots):
        raise PermissionError("external eligibility report must be inside an explicitly allowed report root")
    if not evidence_roots:
        raise PermissionError("at least one explicitly allowed evidence root is required")
    if not path.is_file():
        raise FileNotFoundError(path)
    expected_hash = str(expected_report_sha256 or "").strip().lower()
    if not _valid_sha256(expected_hash):
        raise ValueError("expected_report_sha256 must be a 64-character SHA-256 digest")
    actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
    if actual_hash != expected_hash:
        raise PermissionError("external eligibility report hash does not match the locked digest")
    report = json.loads(path.read_text(encoding="utf-8"))
    if (not isinstance(report, Mapping) or type(report.get("schema_version")) is not int
            or report.get("schema_version") != 1):
        raise ValueError("unsupported external eligibility report schema")
    if (report.get("dataset") != "Duke GBM-DIR-QA"
            or report.get("screen_scope") != "post_overlap_registration_support_pre_target_mask_access"
            or report.get("target_access_state") != "closed"
            or report.get("target_mask_values_used_for_selection") is not False
            or type(report.get("target_mask_voxel_arrays_read")) is not int
            or report.get("target_mask_voxel_arrays_read") != 0
            or report.get("eligibility_rule_sha256") != DUKE_ELIGIBILITY_RULE_SHA256):
        raise PermissionError("external eligibility report does not prove the locked pre-target metadata screen")

    permitted = report.get("permitted_use")
    if (not isinstance(permitted, Mapping) or permitted.get("access_approved") is not True
            or permitted.get("school_use_approved") is not True):
        raise PermissionError("access and permitted-use approval must pass before external N is used")
    permitted_hash = _verify_reviewed_gate(
        path, permitted, field="permitted_use",
        allowed_evidence_roots=evidence_roots,
        status_required=False,
    )

    source_artifacts = report.get("source_metadata_manifests")
    if not isinstance(source_artifacts, Mapping) or not source_artifacts:
        raise ValueError("external eligibility report must identify source metadata manifest files")
    source_hashes = {}
    for name, artifact in source_artifacts.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("source metadata manifest names must be nonempty strings")
        _, source_hashes[name] = _verify_report_artifact(
            path, artifact, field=f"source_metadata_manifests.{name}",
            allowed_evidence_roots=evidence_roots,
        )

    gate_hashes = {"permitted_use": permitted_hash}
    for field in ("overlap_check", "date_quality", "registration_qc", "target_support", "mask_provenance"):
        gate_hashes[field] = _verify_reviewed_gate(
            path, report.get(field), field=field, allowed_evidence_roots=evidence_roots,
        )

    support = report["target_support"]
    if (support.get("status") != "passed" or support.get("coverage_scope") != "entire_image_grid"
            or support.get("zero_is_valid_negative") is not True):
        raise PermissionError("full-grid target-support gate is not recorded as passed")

    selected_ids_path, _ = _verify_report_artifact(
        path, report.get("selected_patient_ids_file"), field="selected_patient_ids_file",
        allowed_evidence_roots=evidence_roots,
    )
    selected_ids = [line.strip() for line in selected_ids_path.read_text(encoding="utf-8").splitlines()]
    if any(not patient_id or "\0" in patient_id for patient_id in selected_ids):
        raise ValueError("selected patient-ID file must contain only nonempty one-per-line IDs")
    if len(set(selected_ids)) != len(selected_ids):
        raise ValueError("selected patient-ID file contains duplicate IDs")
    patient_id_hash = hashlib.sha256("\0".join(sorted(selected_ids)).encode("utf-8")).hexdigest()
    reported_patient_id_hash = str(report.get("selected_patient_ids_sha256_sorted_nul_joined", "")).lower()
    if not _valid_sha256(reported_patient_id_hash) or reported_patient_id_hash != patient_id_hash:
        raise ValueError("selected patient-ID set does not match its reported hash")

    eligibility = report.get("eligibility")
    if (not isinstance(eligibility, Mapping) or eligibility.get("selection_unit") != "patient"):
        raise ValueError("external eligibility report must define patient-level eligibility")
    count = eligibility.get("eligible_patient_count")
    if isinstance(count, bool) or not isinstance(count, int) or count < 2:
        raise ValueError("metadata-only external eligible-patient count must be an integer of at least two")
    if count != len(selected_ids):
        raise ValueError("external eligible-patient count does not match the selected patient-ID file")

    evidence = ExternalEligibilityEvidence(
        eligible_patient_count=int(count),
        eligibility_rule_sha256=DUKE_ELIGIBILITY_RULE_SHA256,
        report_path=str(path),
        report_sha256=actual_hash,
        selected_patient_ids_sha256_sorted_nul_joined=patient_id_hash,
        source_manifest_sha256=MappingProxyType(dict(source_hashes)),
        gate_evidence_sha256=MappingProxyType(dict(gate_hashes)),
    )
    _VALIDATED_ELIGIBILITY_REPORTS.add(evidence)
    return evidence


def precision_gate_from_patient_differences(
    patient_differences: Iterable[tuple[str, float]] | Mapping[str, float], *,
    expected_patient_ids: Iterable[str],
    external_eligibility: ExternalEligibilityEvidence,
) -> dict:
    """Bootstrap paired mean Dice difference at a planned external patient count.

    Rows must contain exactly one finite paired value for every expected patient. Patients are resampled with
    replacement; visit/pair rows are never treated as independent observations.
    """
    if not is_validated_external_eligibility(external_eligibility):
        raise PermissionError("precision calculation requires a validated metadata-only external eligibility report")
    target_n = external_eligibility.eligible_patient_count
    if isinstance(expected_patient_ids, (str, bytes)):
        raise ValueError("expected_patient_ids must be a collection of patient identifiers")
    expected = []
    for raw_patient_id in expected_patient_ids:
        patient_id = str(raw_patient_id).strip()
        if not patient_id or "\0" in patient_id or patient_id in expected:
            raise ValueError("expected patient IDs must be nonempty, unique, and contain no NUL characters")
        expected.append(patient_id)
    if len(expected) < 2:
        raise ValueError("at least two expected patients are required")
    expected_set = set(expected)
    items = patient_differences.items() if isinstance(patient_differences, Mapping) else patient_differences
    normalized: list[tuple[str, float]] = []
    seen: set[str] = set()
    for item in items:
        if not isinstance(item, (tuple, list)) or len(item) != 2:
            raise ValueError("patient_differences must contain (patient_id, difference) pairs")
        patient_id = str(item[0]).strip()
        if not patient_id or patient_id in seen:
            raise ValueError("patient IDs must be nonempty and unique")
        try:
            difference = float(item[1])
        except (TypeError, ValueError, OverflowError) as exc:
            raise ValueError("paired differences must be finite numbers") from exc
        if not np.isfinite(difference):
            raise ValueError("paired differences must be finite numbers")
        seen.add(patient_id)
        normalized.append((patient_id, difference))
    observed_set = seen
    if observed_set != expected_set:
        raise ValueError(
            "paired differences must cover exactly the expected patient cohort "
            f"(missing={len(expected_set - observed_set)}, extra={len(observed_set - expected_set)})"
        )

    normalized.sort(key=lambda item: item[0])
    values = np.asarray([value for _, value in normalized], dtype=float)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    indices = rng.integers(0, len(values), size=(BOOTSTRAP_REPS, int(target_n)))
    bootstrap_means = values[indices].mean(axis=1)
    lower, upper = np.quantile(bootstrap_means, [0.025, 0.975])
    half_width = float((upper - lower) / 2.0)
    expected_id_hash = hashlib.sha256("\0".join(sorted(expected)).encode("utf-8")).hexdigest()
    return {
        "n_development_patients": int(len(values)),
        "n_expected_patients": int(len(expected_set)),
        "expected_patient_ids_sha256_sorted_nul_joined": expected_id_hash,
        "target_external_patient_count": int(target_n),
        "external_eligibility_rule_sha256": external_eligibility.eligibility_rule_sha256,
        "external_eligibility_report_sha256": external_eligibility.report_sha256,
        "external_eligibility_report_path": external_eligibility.report_path,
        "external_selected_patient_ids_sha256_sorted_nul_joined": (
            external_eligibility.selected_patient_ids_sha256_sorted_nul_joined
        ),
        "external_source_manifest_sha256": dict(external_eligibility.source_manifest_sha256),
        "external_gate_evidence_sha256": dict(external_eligibility.gate_evidence_sha256),
        "n_bootstrap": BOOTSTRAP_REPS,
        "seed": BOOTSTRAP_SEED,
        "mean_paired_dice_difference_development_only": float(values.mean()),
        "simulated_precision_interval_95": [float(lower), float(upper)],
        "simulated_precision_half_width_dice": half_width,
        "required_max_half_width_dice": MAX_HALF_WIDTH_DICE,
        "precision_gate_pass": bool(half_width <= MAX_HALF_WIDTH_DICE),
        "interpretation": "development-only study-resolution estimate; not an external effect estimate or clinical threshold",
        "bootstrap_unit": "patient; one paired difference per unique patient",
    }
