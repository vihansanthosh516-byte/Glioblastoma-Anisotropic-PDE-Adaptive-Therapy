"""Fail-closed validation for the Track A sample-to-donor manifest."""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd

# A person must verify these rows against the cited source. The code checks only that a reference exists.
SOURCE_VERIFIED_STATUSES = frozenset({"source_verified"})
MISSING_ID_TOKENS = frozenset({"", "na", "n/a", "nan", "none", "null", "<na>", "unknown",
                               "not known", "not available", "not recorded", "not reported", "missing"})
COUNT_COLUMNS = {"Core": "n_core", "Periphery": "n_periphery", "Healthy": "n_healthy"}


def _clean_id(value) -> str | None:
    try:
        missing = pd.isna(value)
        if pd.api.types.is_scalar(missing) and bool(missing):
            return None
    except (TypeError, ValueError):
        pass
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    text = str(value).strip()
    return None if text.casefold() in MISSING_ID_TOKENS else text


def known_sample_ids(values: Sequence) -> list[str]:
    """Normalize observed sample IDs and omit rows with missing-ID sentinels."""
    return sorted({sample_id for value in values if (sample_id := _clean_id(value)) is not None})


def normalize_sample_ids(values: Sequence) -> list[str | None]:
    """Normalize IDs in-place order for per-cell mapping without repeated set/sort work."""
    return [_clean_id(value) for value in values]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_manifest_provenance(
    manifest_path: Path,
    summary_path: Path,
    input_path: Path,
    source_script_path: Path,
    dependency_paths: Mapping[str, Path],
    *,
    expected_input_relative_path: str,
) -> dict[str, object]:
    """Reject donor tables not tied to the current source data and generator code."""
    required_paths = (
        (manifest_path, "manifest"),
        (summary_path, "manifest summary"),
        (input_path, "source data"),
        (source_script_path, "manifest generator"),
    )
    for path, label in required_paths:
        if not Path(path).is_file():
            raise FileNotFoundError(f"Track A {label} file is missing: {path}")
    try:
        summary = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError("Track A manifest summary is unreadable") from exc

    actual_hashes = {
        "input_sha256": sha256_file(input_path),
        "manifest_sha256": sha256_file(manifest_path),
        "source_script_sha256": sha256_file(source_script_path),
    }
    actual_dependencies = {
        name: sha256_file(path) for name, path in sorted(dependency_paths.items())
    }
    if summary.get("input_file") != expected_input_relative_path:
        raise ValueError("Track A manifest summary names a different source dataset")
    for key, actual in actual_hashes.items():
        if summary.get(key) != actual:
            raise ValueError(f"Track A manifest provenance mismatch for {key}")
    if summary.get("dependency_sha256") != actual_dependencies:
        raise ValueError("Track A manifest provenance mismatch for dependency_sha256")
    return summary


def load_identity_policy(path: Path) -> dict:
    """Read an identity sensitivity JSON object and reject repeated JSON keys."""
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"identity policy contains a repeated JSON key: {key}")
            result[key] = value
        return result

    try:
        policy = json.loads(Path(path).read_text(encoding="utf-8"), object_pairs_hook=unique_object)
    except json.JSONDecodeError as exc:
        raise ValueError("identity policy is not valid JSON") from exc
    if not isinstance(policy, Mapping):
        raise ValueError("identity policy JSON must contain one object")
    return dict(policy)


def identity_resolution_label(status_counts: Mapping[str, int], excluded_sample_count: int) -> str:
    """Describe which non-source-verified manifest rows were actually included in scoring."""
    unverified = sum(
        int(count) for status, count in status_counts.items() if status not in SOURCE_VERIFIED_STATUSES
    )
    excluded = int(excluded_sample_count)
    if excluded < 0 or excluded > unverified:
        raise ValueError("excluded sample count must be between zero and the unverified manifest-row count")
    included = unverified - excluded
    if included:
        return f"{included} non-source-verified sample identities were included under explicit sensitivity assumptions"
    if unverified:
        return f"all {unverified} non-source-verified sample identities were excluded from scoring"
    return "identity rows were marked source_verified by human review; evidence references were not checked by code"


def validate_manifest_mapping(manifest: pd.DataFrame, expected_sample_ids: Sequence,
                              sensitivity_policy: Mapping | None = None) -> tuple[pd.Series, frozenset[str]]:
    """Validate exact sample coverage and return donor mapping plus explicitly excluded samples.

    Only ``source_verified`` rows with a nonempty human evidence reference pass by default. The code does not
    inspect that reference; all rule-derived identities need an explicit sensitivity policy. The policy must
    account for each such sample exactly once using ``{"groups": {sample_id: donor_id}, "exclude": [...]}``.
    Group targets must already occur as donor IDs in the manifest. Policy never silently affects verified rows.
    """
    required = {"sample_id", "donor_id", "identity_status", "identity_evidence"}
    missing_columns = sorted(required - set(manifest.columns))
    if missing_columns:
        raise ValueError(f"donor manifest lacks required identity columns: {missing_columns}")
    if manifest.empty:
        raise ValueError("donor manifest is empty")

    df = manifest[["sample_id", "donor_id", "identity_status", "identity_evidence"]].copy()
    df["sample_id"] = df["sample_id"].map(_clean_id)
    df["donor_id"] = df["donor_id"].map(_clean_id)
    if df[["sample_id", "donor_id"]].isna().any().any():
        raise ValueError("donor manifest contains missing sample or donor IDs")
    if df["sample_id"].duplicated().any():
        duplicates = sorted(df.loc[df["sample_id"].duplicated(keep=False), "sample_id"].unique())
        raise ValueError(f"donor manifest has duplicate sample IDs: {duplicates[:5]}")

    expected = [_clean_id(value) for value in expected_sample_ids]
    if any(value is None for value in expected):
        raise ValueError("observed sample IDs contain missing values; filter them before manifest mapping")
    expected_set = set(expected)
    if len(expected_set) != len(expected):
        raise ValueError("expected sample IDs must be unique")
    manifest_set = set(df["sample_id"])
    missing = sorted(expected_set - manifest_set)
    extra = sorted(manifest_set - expected_set)
    if missing or extra:
        raise ValueError(f"donor manifest sample coverage mismatch (missing={missing[:5]}, extra={extra[:5]})")

    df["identity_status"] = df["identity_status"].fillna("").astype(str).str.strip()
    df["identity_evidence"] = df["identity_evidence"].map(_clean_id)
    verified_without_evidence = df.loc[
        df["identity_status"].isin(SOURCE_VERIFIED_STATUSES) & df["identity_evidence"].isna(), "sample_id"
    ]
    if not verified_without_evidence.empty:
        raise ValueError("source-verified donor rows must include an identity evidence reference")
    unresolved = set(df.loc[~df["identity_status"].isin(SOURCE_VERIFIED_STATUSES), "sample_id"])
    groups: dict[str, str] = {}
    excluded: set[str] = set()
    if unresolved:
        if sensitivity_policy is None:
            examples = sorted(unresolved)[:5]
            raise ValueError(f"donor-level scoring blocked by unresolved sample identities: {examples}")
        if not isinstance(sensitivity_policy, Mapping):
            raise ValueError("sensitivity policy must be a mapping")
        if set(sensitivity_policy) != {"groups", "exclude"}:
            raise ValueError("sensitivity policy must contain exactly 'groups' and 'exclude'")
        raw_groups = sensitivity_policy["groups"]
        raw_excluded = sensitivity_policy["exclude"]
        if not isinstance(raw_groups, Mapping) or not isinstance(raw_excluded, list):
            raise ValueError("sensitivity policy groups must be a mapping and exclude must be a list")
        groups = {}
        for sample, target in raw_groups.items():
            sample_id, donor_id = _clean_id(sample), _clean_id(target)
            if sample_id is None or donor_id is None:
                raise ValueError("sensitivity policy contains a missing sample or donor ID")
            if sample_id in groups:
                raise ValueError(f"sensitivity policy has duplicate normalized sample key: {sample_id}")
            groups[sample_id] = donor_id
        try:
            excluded_values = list(raw_excluded)
        except TypeError as exc:
            raise ValueError("sensitivity policy exclude must be iterable") from exc
        normalized_excluded = [_clean_id(value) for value in excluded_values]
        if any(value is None for value in normalized_excluded):
            raise ValueError("sensitivity policy contains a missing excluded sample ID")
        if len(set(normalized_excluded)) != len(normalized_excluded):
            raise ValueError("sensitivity policy lists an excluded sample more than once")
        excluded = set(normalized_excluded)
        accounted = set(groups) | excluded
        if set(groups) & excluded:
            raise ValueError("a sample cannot be both grouped and excluded by sensitivity policy")
        if accounted != unresolved:
            raise ValueError("sensitivity policy must group or exclude every unresolved sample and no resolved sample")
        donor_ids = set(df["donor_id"])
        unknown_targets = sorted(set(groups.values()) - donor_ids)
        if unknown_targets:
            raise ValueError(f"sensitivity policy group targets are absent from donor manifest: {unknown_targets[:5]}")
    elif sensitivity_policy is not None:
        if not isinstance(sensitivity_policy, Mapping) or set(sensitivity_policy) != {"groups", "exclude"}:
            raise ValueError("sensitivity policy must contain exactly 'groups' and 'exclude'")
        if sensitivity_policy["groups"] or sensitivity_policy["exclude"]:
            raise ValueError("sensitivity policy may only address unresolved samples")

    donor_by_sample = df.set_index("sample_id")["donor_id"].copy()
    for sample_id, donor_id in groups.items():
        donor_by_sample.loc[sample_id] = donor_id
    return donor_by_sample.reindex(expected), frozenset(excluded)


def _count_value(value, column: str) -> int:
    try:
        numeric = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"manifest {column} must contain non-negative integer counts") from exc
    if not math.isfinite(numeric) or numeric < 0 or not numeric.is_integer():
        raise ValueError(f"manifest {column} must contain non-negative integer counts")
    return int(numeric)


def validate_manifest_cell_counts(manifest: pd.DataFrame, sample_ids: Sequence, class_labels: Sequence,
                                  donor_by_sample: pd.Series,
                                  excluded_samples: frozenset[str] = frozenset()) -> dict[str, dict[str, int]]:
    """Check manifest class totals against observed labels without inspecting expression data."""
    if len(sample_ids) != len(class_labels):
        raise ValueError("sample IDs and class labels must have equal length")
    needed = {"sample_id", *COUNT_COLUMNS.values()}
    if needed - set(manifest.columns):
        raise ValueError(f"donor manifest lacks count columns: {sorted(needed - set(manifest.columns))}")

    expected_by_sample: dict[str, Counter] = {}
    for row in manifest.to_dict("records"):
        sample_id = _clean_id(row["sample_id"])
        if sample_id in excluded_samples:
            continue
        if sample_id not in donor_by_sample.index:
            raise ValueError(f"manifest sample is absent from validated donor map: {sample_id}")
        counts = expected_by_sample.setdefault(sample_id, Counter())
        for class_label, column in COUNT_COLUMNS.items():
            counts[class_label] += _count_value(row[column], column)

    observed_by_sample: dict[str, Counter] = {}
    for raw_sample, raw_class in zip(sample_ids, class_labels):
        sample_id = _clean_id(raw_sample)
        if sample_id is None or sample_id in excluded_samples or raw_class is None:
            continue
        if raw_class not in COUNT_COLUMNS:
            raise ValueError(f"observed class label is absent from donor manifest schema: {raw_class!r}")
        if sample_id not in donor_by_sample.index:
            raise ValueError(f"observed sample is absent from validated donor map: {sample_id}")
        observed_by_sample.setdefault(sample_id, Counter())[raw_class] += 1

    sample_ids = set(expected_by_sample) | set(observed_by_sample)
    mismatches = []
    for sample_id in sorted(sample_ids):
        for class_label in COUNT_COLUMNS:
            wanted = expected_by_sample.get(sample_id, {}).get(class_label, 0)
            found = observed_by_sample.get(sample_id, {}).get(class_label, 0)
            if wanted != found:
                mismatches.append((sample_id, class_label, wanted, found))
    if mismatches:
        raise ValueError(f"donor manifest per-sample cell counts disagree with observed labels: {mismatches[:5]}")

    expected_by_donor: dict[str, Counter] = {}
    for sample_id, counts in expected_by_sample.items():
        donor_id = str(donor_by_sample.loc[sample_id])
        expected_by_donor.setdefault(donor_id, Counter()).update(counts)
    return {donor_id: {label: int(expected_by_donor.get(donor_id, {}).get(label, 0)) for label in COUNT_COLUMNS}
            for donor_id in sorted(expected_by_donor)}


def manifest_gbm_donors(manifest: pd.DataFrame, donor_by_sample: pd.Series,
                        excluded_samples: frozenset[str] = frozenset()) -> frozenset[str]:
    """Return donors whose included manifest rows are all explicitly marked GBM."""
    if "is_gbm" not in manifest.columns:
        raise ValueError("donor manifest lacks required diagnosis field 'is_gbm'")
    flags_by_donor: dict[str, list[bool]] = {}
    for row in manifest.to_dict("records"):
        sample_id = _clean_id(row["sample_id"])
        if sample_id in excluded_samples:
            continue
        if sample_id not in donor_by_sample.index:
            raise ValueError(f"manifest sample is absent from validated donor map: {sample_id}")
        flag = str(row["is_gbm"]).strip().casefold()
        if flag in {"true", "1", "yes"}:
            is_gbm = True
        elif flag in {"false", "0", "no"}:
            is_gbm = False
        else:
            raise ValueError(f"manifest is_gbm must be boolean for sample {sample_id}")
        donor_id = str(donor_by_sample.loc[sample_id])
        flags_by_donor.setdefault(donor_id, []).append(is_gbm)
    return frozenset(donor_id for donor_id, flags in flags_by_donor.items() if flags and all(flags))
