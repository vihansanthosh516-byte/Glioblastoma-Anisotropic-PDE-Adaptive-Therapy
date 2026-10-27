"""Synthetic fail-closed tests for the Track A donor manifest gate."""
from __future__ import annotations

import json

import pandas as pd
import pytest

from src import trackA_identity as identity


def manifest():
    return pd.DataFrame([
        {"sample_id": "S1", "donor_id": "D1", "identity_status": "source_verified",
         "identity_evidence": "crosswalk.tsv#S1", "is_gbm": True,
         "n_core": 1, "n_periphery": 0, "n_healthy": 0},
        {"sample_id": "S2", "donor_id": "S2", "identity_status": "unresolved_missing_metadata",
         "identity_evidence": "", "is_gbm": True, "n_core": 0, "n_periphery": 1, "n_healthy": 0},
        {"sample_id": "S3", "donor_id": "S3", "identity_status": "unresolved_metadata_conflict",
         "identity_evidence": "", "is_gbm": False, "n_core": 0, "n_periphery": 0, "n_healthy": 1},
    ])


def test_unresolved_manifest_rows_block_default_donor_mapping():
    with pytest.raises(ValueError, match="unresolved sample identities"):
        identity.validate_manifest_mapping(manifest(), ["S1", "S2", "S3"])


def test_explicit_group_and_exclusion_policy_maps_synthetic_samples():
    donor_by_sample, excluded = identity.validate_manifest_mapping(
        manifest(), ["S1", "S2", "S3"],
        sensitivity_policy={"groups": {"S2": "D1"}, "exclude": ["S3"]},
    )
    assert donor_by_sample.to_dict() == {"S1": "D1", "S2": "D1", "S3": "S3"}
    assert excluded == frozenset({"S3"})
    counts = identity.validate_manifest_cell_counts(
        manifest(), ["S1", "S2", "S3"], ["Core", "Periphery", "Healthy"], donor_by_sample, excluded
    )
    assert counts == {"D1": {"Core": 1, "Periphery": 1, "Healthy": 0}}


def test_unverified_single_sample_id_blocks_default_scoring():
    data = manifest().iloc[[0]].copy()
    data.loc[data.index[0], "identity_status"] = "assumed_single_sample_id"
    data.loc[data.index[0], "identity_evidence"] = ""
    with pytest.raises(ValueError, match="unresolved sample identities"):
        identity.validate_manifest_mapping(data, ["S1"])


def test_source_verified_status_requires_evidence_reference():
    data = manifest().iloc[[0]].copy()
    data.loc[data.index[0], "identity_evidence"] = ""
    with pytest.raises(ValueError, match="identity evidence reference"):
        identity.validate_manifest_mapping(data, ["S1"])


def test_sensitivity_policy_must_account_for_each_unresolved_sample_once():
    with pytest.raises(ValueError, match="every unresolved sample"):
        identity.validate_manifest_mapping(
            manifest(), ["S1", "S2", "S3"], sensitivity_policy={"groups": {"S2": "D1"}, "exclude": []}
        )


def test_sensitivity_policy_rejects_duplicate_keys_after_normalization():
    with pytest.raises(ValueError, match="duplicate normalized sample key"):
        identity.validate_manifest_mapping(
            manifest(), ["S1", "S2", "S3"],
            sensitivity_policy={"groups": {"S2": "D1", " S2 ": "D1", "S3": "D1"}, "exclude": []},
        )


def test_sensitivity_policy_requires_exclude_to_be_a_list():
    with pytest.raises(ValueError, match="exclude must be a list"):
        identity.validate_manifest_mapping(
            manifest(), ["S1", "S2", "S3"],
            sensitivity_policy={"groups": {"S2": "D1", "S3": "D1"}, "exclude": {"unused": True}},
        )


def test_manifest_mapping_requires_exact_sample_coverage():
    resolved = manifest().iloc[[0]].copy()
    with pytest.raises(ValueError, match="coverage mismatch"):
        identity.validate_manifest_mapping(resolved, ["S1", "S2"])


def test_manifest_counts_must_match_the_observed_synthetic_labels():
    resolved = manifest().iloc[[0]].copy()
    donor_by_sample, excluded = identity.validate_manifest_mapping(resolved, ["S1"])
    with pytest.raises(ValueError, match="counts disagree"):
        identity.validate_manifest_cell_counts(resolved, ["S1"], ["Periphery"], donor_by_sample, excluded)


def test_manifest_count_validation_rejects_unrecognized_nonmissing_labels():
    resolved = manifest().iloc[[0]].copy()
    donor_by_sample, excluded = identity.validate_manifest_mapping(resolved, ["S1"])
    with pytest.raises(ValueError, match="absent from donor manifest schema"):
        identity.validate_manifest_cell_counts(resolved, ["S1"], ["UnknownRegion"], donor_by_sample, excluded)


def test_manifest_count_validation_checks_each_sample_before_donor_aggregation():
    data = pd.DataFrame([
        {"sample_id": "S1", "donor_id": "D", "identity_status": "source_verified",
         "identity_evidence": "crosswalk.tsv#S1", "is_gbm": True, "n_core": 1, "n_periphery": 0,
         "n_healthy": 0},
        {"sample_id": "S2", "donor_id": "D", "identity_status": "source_verified",
         "identity_evidence": "crosswalk.tsv#S2", "is_gbm": True, "n_core": 0, "n_periphery": 1,
         "n_healthy": 0},
    ])
    donor_by_sample, excluded = identity.validate_manifest_mapping(data, ["S1", "S2"])
    with pytest.raises(ValueError, match="per-sample cell counts"):
        identity.validate_manifest_cell_counts(data, ["S1", "S2"], ["Periphery", "Core"],
                                               donor_by_sample, excluded)


def test_sample_id_normalization_preserves_cell_order_and_missing_rows():
    assert identity.normalize_sample_ids([" S1 ", None, "unknown", b"S2"]) == ["S1", None, None, "S2"]


def test_eligible_non_gbm_core_periphery_donor_is_not_selected():
    data = pd.DataFrame([
        {"sample_id": "G", "donor_id": "GBM1", "identity_status": "source_verified",
         "identity_evidence": "crosswalk.tsv#G", "is_gbm": True,
         "n_core": 1, "n_periphery": 1, "n_healthy": 0},
        {"sample_id": "N", "donor_id": "NONGBM1", "identity_status": "source_verified",
         "identity_evidence": "crosswalk.tsv#N", "is_gbm": False,
         "n_core": 1, "n_periphery": 1, "n_healthy": 0},
    ])
    donor_by_sample, excluded = identity.validate_manifest_mapping(data, ["G", "N"])
    counts = identity.validate_manifest_cell_counts(data, ["G", "G", "N", "N"],
                                                     ["Core", "Periphery", "Core", "Periphery"],
                                                     donor_by_sample, excluded)
    gbm = identity.manifest_gbm_donors(data, donor_by_sample, excluded)
    eligible = sorted(d for d, values in counts.items()
                      if d in gbm and values["Core"] > 0 and values["Periphery"] > 0)
    assert eligible == ["GBM1"]


def test_metadata_only_candidate_status_blocks_without_policy():
    data = manifest().iloc[[0]].copy()
    data.loc[data.index[0], "identity_status"] = "unresolved_metadata_only_candidate"
    with pytest.raises(ValueError, match="unresolved sample identities"):
        identity.validate_manifest_mapping(data, ["S1"])


def test_manifest_provenance_binds_dataset_manifest_and_generator(tmp_path):
    manifest_path = tmp_path / "donors.csv"
    summary_path = tmp_path / "summary.json"
    input_path = tmp_path / "source.h5ad"
    script_path = tmp_path / "generator.py"
    dependencies = {"donor_cv": tmp_path / "donor_cv.py", "leakage_checks": tmp_path / "leakage_checks.py"}
    for path, content in ((manifest_path, b"sample_id\nS1\n"), (input_path, b"source-data"),
                          (script_path, b"generator-code"), (dependencies["donor_cv"], b"donor-map-code"),
                          (dependencies["leakage_checks"], b"leakage-check-code")):
        path.write_bytes(content)
    summary_path.write_text(json.dumps({
        "input_file": "data/source.h5ad",
        "input_sha256": identity.sha256_file(input_path),
        "manifest_sha256": identity.sha256_file(manifest_path),
        "source_script_sha256": identity.sha256_file(script_path),
        "dependency_sha256": {name: identity.sha256_file(path) for name, path in dependencies.items()},
    }), encoding="utf-8")

    identity.verify_manifest_provenance(manifest_path, summary_path, input_path, script_path, dependencies,
                                        expected_input_relative_path="data/source.h5ad")

    input_path.write_bytes(b"edited-source-data")
    with pytest.raises(ValueError, match="input_sha256"):
        identity.verify_manifest_provenance(manifest_path, summary_path, input_path, script_path, dependencies,
                                            expected_input_relative_path="data/source.h5ad")


def test_identity_policy_loader_rejects_duplicate_json_keys(tmp_path):
    path = tmp_path / "policy.json"
    path.write_text('{"groups": {"S2": "D1"}, "groups": {"S2": "D2"}, "exclude": []}',
                    encoding="utf-8")
    with pytest.raises(ValueError, match="repeated JSON key"):
        identity.load_identity_policy(path)


def test_identity_resolution_label_distinguishes_assumed_and_excluded_rows():
    counts = {"source_verified": 2, "assumed_single_sample_id": 3}
    assert "3 non-source-verified" in identity.identity_resolution_label(counts, 0)
    assert "all 3 non-source-verified" in identity.identity_resolution_label(counts, 3)
    assert "human review" in identity.identity_resolution_label({"source_verified": 2}, 0)
    with pytest.raises(ValueError, match="between zero"):
        identity.identity_resolution_label(counts, 4)
