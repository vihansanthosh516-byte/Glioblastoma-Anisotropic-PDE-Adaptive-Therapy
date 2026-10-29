from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_mu_patient_runner  # noqa: E402
from src.fk_precision_gate import (  # noqa: E402
    BOOTSTRAP_REPS,
    BOOTSTRAP_SEED,
    DUKE_ELIGIBILITY_RULE_SHA256,
    MAX_HALF_WIDTH_DICE,
    ExternalEligibilityEvidence,
    precision_gate_from_patient_differences,
    validate_external_eligibility_report,
)


def _write_artifact(root: Path, relative_path: str, contents: str) -> dict[str, str]:
    path = root / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents, encoding="utf-8")
    return {"path": relative_path, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def _write_report(root: Path, mutate=None) -> tuple[Path, str]:
    selected_ids = [f"synthetic-patient-{index:02d}" for index in range(7)]
    id_file = _write_artifact(root, "evidence/selected_ids.txt", "\n".join(selected_ids) + "\n")
    patient_id_hash = hashlib.sha256("\0".join(sorted(selected_ids)).encode("utf-8")).hexdigest()

    def reviewed_gate(name: str, **fields) -> dict:
        return {
            **_write_artifact(root, f"evidence/{name}.txt", f"synthetic evidence for {name}\n"),
            "status": "passed",
            "reviewed_by": "synthetic test reviewer",
            "reviewed_at": "2026-10-28T12:00:00Z",
            **fields,
        }

    report = {
        "schema_version": 1,
        "dataset": "Duke GBM-DIR-QA",
        "screen_scope": "post_overlap_registration_support_pre_target_mask_access",
        "target_access_state": "closed",
        "target_mask_values_used_for_selection": False,
        "target_mask_voxel_arrays_read": 0,
        "eligibility_rule_sha256": DUKE_ELIGIBILITY_RULE_SHA256,
        "permitted_use": reviewed_gate(
            "permitted_use", access_approved=True, school_use_approved=True,
        ),
        "source_metadata_manifests": {
            "visit_metadata": _write_artifact(
                root, "evidence/visit_metadata.csv", "patient_id,visit_day\nsynthetic,0\n",
            ),
        },
        "overlap_check": reviewed_gate("overlap_check"),
        "date_quality": reviewed_gate("date_quality"),
        "registration_qc": reviewed_gate("registration_qc"),
        "target_support": reviewed_gate(
            "target_support", coverage_scope="entire_image_grid", zero_is_valid_negative=True,
        ),
        "mask_provenance": reviewed_gate("mask_provenance"),
        "selected_patient_ids_file": id_file,
        "selected_patient_ids_sha256_sorted_nul_joined": patient_id_hash,
        "eligibility": {
            "selection_unit": "patient",
            "eligible_patient_count": len(selected_ids),
        },
    }
    if mutate is not None:
        mutate(report)
    path = root / "eligibility_report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, sort_keys=True, separators=(",", ":")), encoding="utf-8")
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return path, digest


def _validated_report(root: Path):
    path, digest = _write_report(root)
    return validate_external_eligibility_report(
        path,
        expected_report_sha256=digest,
        allowed_report_roots=(root,),
        allowed_evidence_roots=(root,),
    )


def test_negative_mean_can_pass_precision_without_being_positive(tmp_path):
    patient_ids = [f"development-patient-{index:02d}" for index in range(18)]
    differences = [(patient_id, -0.10) for patient_id in patient_ids]

    result = precision_gate_from_patient_differences(
        differences,
        expected_patient_ids=patient_ids,
        external_eligibility=_validated_report(tmp_path),
    )

    assert result["n_development_patients"] == 18
    assert result["target_external_patient_count"] == 7  # synthetic report fixture only
    assert result["n_bootstrap"] == BOOTSTRAP_REPS == 10_000
    assert result["seed"] == BOOTSTRAP_SEED == 20261028
    assert np.allclose(result["simulated_precision_interval_95"], [-0.10, -0.10])
    assert result["simulated_precision_half_width_dice"] == 0.0
    assert result["required_max_half_width_dice"] == MAX_HALF_WIDTH_DICE == 0.05
    assert result["precision_gate_pass"] is True
    assert np.isclose(result["mean_paired_dice_difference_development_only"], -0.10)
    assert result["external_eligibility_rule_sha256"] == DUKE_ELIGIBILITY_RULE_SHA256
    assert result["external_eligibility_report_sha256"]
    assert result["external_gate_evidence_sha256"].keys() == {
        "permitted_use", "overlap_check", "date_quality", "registration_qc",
        "target_support", "mask_provenance",
    }


def test_precision_result_is_deterministic_and_patient_order_invariant(tmp_path):
    eligibility = _validated_report(tmp_path)
    patient_ids = [f"patient-{index:02d}" for index in range(18)]
    differences = list(zip(patient_ids, np.linspace(-0.25, 0.20, len(patient_ids))))
    forward = precision_gate_from_patient_differences(
        differences, expected_patient_ids=patient_ids, external_eligibility=eligibility,
    )
    reversed_order = precision_gate_from_patient_differences(
        list(reversed(differences)), expected_patient_ids=list(reversed(patient_ids)),
        external_eligibility=eligibility,
    )

    assert forward == reversed_order
    assert forward["simulated_precision_half_width_dice"] > 0.05
    assert forward["precision_gate_pass"] is False


@pytest.mark.parametrize(
    "differences, expected_ids",
    [
        ([('p1', 0.1)], ['p1', 'p2']),
        ([('p1', 0.1), ('p2', 0.2), ('extra', 0.0)], ['p1', 'p2']),
        ([('p1', 0.1), ('p1', 0.2)], ['p1', 'p2']),
        ([('p1', float('nan')), ('p2', 0.2)], ['p1', 'p2']),
    ],
)
def test_precision_gate_requires_one_finite_difference_per_patient(
    tmp_path, differences, expected_ids,
):
    with pytest.raises(ValueError):
        precision_gate_from_patient_differences(
            differences,
            expected_patient_ids=expected_ids,
            external_eligibility=_validated_report(tmp_path),
        )


def test_precision_gate_rejects_fabricated_evidence_object():
    forged = ExternalEligibilityEvidence(
        eligible_patient_count=7,
        eligibility_rule_sha256=DUKE_ELIGIBILITY_RULE_SHA256,
        report_path="not-a-validated-report.json",
        report_sha256="a" * 64,
        selected_patient_ids_sha256_sorted_nul_joined="b" * 64,
        source_manifest_sha256={},
        gate_evidence_sha256={},
    )
    with pytest.raises(PermissionError, match="validated metadata-only"):
        precision_gate_from_patient_differences(
            [("p1", 0.0), ("p2", 0.0)],
            expected_patient_ids=["p1", "p2"],
            external_eligibility=forged,
        )


def test_report_hash_and_evidence_file_hashes_are_required(tmp_path):
    path, digest = _write_report(tmp_path)
    with pytest.raises(PermissionError, match="report hash"):
        validate_external_eligibility_report(
            path, expected_report_sha256="0" * 64,
            allowed_report_roots=(tmp_path,), allowed_evidence_roots=(tmp_path,),
        )

    (tmp_path / "evidence/visit_metadata.csv").write_text("changed\n", encoding="utf-8")
    with pytest.raises(PermissionError, match="file hash"):
        validate_external_eligibility_report(
            path, expected_report_sha256=digest,
            allowed_report_roots=(tmp_path,), allowed_evidence_roots=(tmp_path,),
        )


@pytest.mark.parametrize(
    "section, key, value",
    [
        ("permitted_use", "access_approved", False),
        ("permitted_use", "school_use_approved", False),
        ("overlap_check", "status", "unknown"),
        ("date_quality", "status", "unknown"),
        ("registration_qc", "status", "unknown"),
        ("target_support", "status", "unknown"),
        ("target_support", "coverage_scope", "tumor_bounding_box"),
        ("target_support", "zero_is_valid_negative", False),
        ("mask_provenance", "status", "unknown"),
        ("top", "target_access_state", "open"),
        ("top", "target_mask_values_used_for_selection", True),
        ("top", "target_mask_voxel_arrays_read", 1),
    ],
)
def test_report_rejects_any_unresolved_external_gate(tmp_path, section, key, value):
    def mutate(report):
        target = report if section == "top" else report[section]
        target[key] = value

    path, digest = _write_report(tmp_path, mutate=mutate)
    with pytest.raises((PermissionError, ValueError)):
        validate_external_eligibility_report(
            path, expected_report_sha256=digest,
            allowed_report_roots=(tmp_path,), allowed_evidence_roots=(tmp_path,),
        )


@pytest.mark.parametrize("field, invalid_value", [("schema_version", True), ("target_mask_voxel_arrays_read", False)])
def test_report_rejects_boolean_values_for_integer_fields(tmp_path, field, invalid_value):
    path, digest = _write_report(tmp_path, mutate=lambda report: report.__setitem__(field, invalid_value))
    with pytest.raises((PermissionError, ValueError)):
        validate_external_eligibility_report(
            path, expected_report_sha256=digest,
            allowed_report_roots=(tmp_path,), allowed_evidence_roots=(tmp_path,),
        )


def test_report_patient_count_must_match_the_hash_bound_id_file(tmp_path):
    path, _ = _write_report(
        tmp_path,
        mutate=lambda report: report["eligibility"].__setitem__("eligible_patient_count", 8),
    )
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    with pytest.raises(ValueError, match="does not match"):
        validate_external_eligibility_report(
            path, expected_report_sha256=digest,
            allowed_report_roots=(tmp_path,), allowed_evidence_roots=(tmp_path,),
        )


def test_report_must_be_inside_allowed_report_root(tmp_path):
    path, digest = _write_report(tmp_path)
    with pytest.raises(PermissionError, match="allowed report root"):
        validate_external_eligibility_report(
            path, expected_report_sha256=digest,
            allowed_report_roots=(tmp_path / "elsewhere",), allowed_evidence_roots=(tmp_path,),
        )


def test_sensitivity_wrapper_passes_precision_evidence_only_to_primary(tmp_path, monkeypatch):
    eligibility = _validated_report(tmp_path)
    calls = []

    def fake_run(_cohort, **kwargs):
        calls.append(kwargs)
        rows = [
            {"patient_id": patient_id, "scores": {"paired_delta_dice_hybrid_minus_no_pde": 0.01}}
            for patient_id in ("p1", "p2")
        ]
        return {
            "rows": rows,
            "cohort_audit": {"source_hashes": {"selected_patient_ids_sha256_sorted_nul_joined": "c" * 64}},
            "n_successful_hybrid_forecasts": 2,
            "n_persistence_fallbacks": 0,
        }

    monkeypatch.setattr(fk_mu_patient_runner, "run_mu_first_two_to_third_development_cv", fake_run)
    fk_mu_patient_runner.run_mu_first_two_to_third_with_prior_sensitivity(
        None,
        external_eligibility=eligibility,
    )

    assert len(calls) == 4
    assert all(dict(call["initial_boundary_covariance"]) == {
        "sigma_b_mm": 0.8, "tau_b_mm": 0.3, "ell_mm": 6.0,
    } for call in calls)
    assert calls[0]["external_eligibility"] is eligibility
    assert all(call["external_eligibility"] is None for call in calls[1:])
