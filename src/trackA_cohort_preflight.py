"""Inventory Track A cohort eligibility from coded provenance CSV files only."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from .trackA_margin_evaluation import minimum_possible_adjusted_p


ROLES = ("development", "calibration", "evaluation", "shift")
STATUS = {"verified", "missing", "unknown"}
COMMON = {
    "study", "accession", "origin", "item", "role", "condition",
    "reuse_terms", "origin_lineage", "prior_exposure", "evidence",
}
LANE_COLUMNS = {
    "cell-spatial": {"section", "region", "label_basis", "modality", "registration",
                     "spatial_inputs", "independent_audit"},
    "spot-spatial": {"section", "region", "label_basis", "modality", "registration",
                     "spatial_inputs", "independent_audit"},
    "whole-tissue": {"region", "label_basis", "modality"},
    "tmz": {"untreated_rna", "measured_viability", "assay_definition", "independent_culture", "mgmt"},
}
MODALITY = {"cell-spatial": "cell_spatial", "spot-spatial": "spot_spatial", "whole-tissue": "bulk_rna"}


def _enum(row, key, allowed):
    value = row[key]
    if value not in allowed:
        raise ValueError(f"Unsupported {key}: {value!r}")
    return value


def _missing(row, *fields):
    return [field for field in fields if _enum(row, field, STATUS) != "verified"]


def audit(rows, *, lane, minimum_evaluation, minimum_shift, minimum_repeated_development,
          comparison_count=None, alpha=None):
    """Audit declared evidence. Confirm provenance and fit/power separately."""
    if lane not in LANE_COLUMNS or any(
        type(n) is not int or n < 1
        for n in (minimum_evaluation, minimum_shift)
    ):
        raise ValueError("Provide a valid lane and positive origin floors")
    if lane == "tmz":
        if comparison_count is not None or alpha is not None:
            raise ValueError("A paired region sign-test design does not apply to TMZ")
        if type(minimum_repeated_development) is not int or minimum_repeated_development < 1:
            raise ValueError("TMZ preflight needs a planned repeated-origin floor")
    elif (type(comparison_count) is not int or comparison_count < 1
          or alpha is None or not 0 < alpha < 1):
        raise ValueError("Regional preflight needs the planned comparison count and alpha")

    expected = COMMON | LANE_COLUMNS[lane]
    items, origins, cohorts, studies = set(), {}, {}, {}
    included, rejected = defaultdict(list), []
    observed = Counter()
    for row in rows:
        if set(row) != expected or any(
            not isinstance(value, str) or not value or value.strip() != value for value in row.values()
        ):
            raise ValueError("Provide exact metadata fields with nonempty, trimmed text values")
        study, accession, origin, item, role, condition = (
            row[name] for name in ("study", "accession", "origin", "item", "role", "condition")
        )
        _enum(row, "role", ROLES)
        observed[role] += 1
        if (study, item) in items:
            raise ValueError(f"Duplicate study and item: {study}/{item}")
        items.add((study, item))
        if origins.setdefault(origin, (role, study, condition)) != (role, study, condition):
            raise ValueError(f"Origin {origin} crosses a role, study or assay condition")
        # Calibration and model fitting can come from distinct donors in a
        # development cohort. An evaluation cohort stays outside both.
        family = "development" if role in {"development", "calibration"} else role
        if cohorts.setdefault(study, family) != family:
            raise ValueError(f"Study {study} crosses study roles")
        if studies.setdefault(accession, study) != study:
            raise ValueError(f"Accession {accession} belongs to multiple studies")

        reasons = _missing(row, "reuse_terms", "origin_lineage")
        reasons.extend(field for field in ("condition", "evidence") if row[field] == "unknown")
        exposure = _enum(row, "prior_exposure", {"untouched", "seen", "unknown"})
        if role in {"evaluation", "shift"} and exposure != "untouched":
            reasons.append("prior_exposure")
        if lane == "tmz":
            reasons.extend(_missing(row, "untreated_rna", "measured_viability",
                                    "assay_definition", "independent_culture"))
            _enum(row, "mgmt", STATUS)
        else:
            if _enum(row, "region", {"core", "infiltrating", "unknown"}) == "unknown":
                reasons.append("region")
            if _enum(row, "label_basis", {"independent_pathology", "other", "unknown"}) != "independent_pathology":
                reasons.append("label_basis")
            if _enum(row, "modality", set(MODALITY.values()) | {"dissociated_rna", "unknown"}) != MODALITY[lane]:
                reasons.append("modality")
            if lane.endswith("-spatial"):
                reasons.extend(_missing(row, "registration"))
                if row["section"] == "unknown":
                    reasons.append("section")
                _enum(row, "spatial_inputs", STATUS)
                _enum(row, "independent_audit", STATUS)
        if reasons:
            rejected.append({"accession": accession, "origin": origin, "item": item, "reasons": reasons})
        else:
            included[origin].append(row)

    report = {role: {"rows": observed[role], "eligible_origins": [], "extra_support_origins": [],
                     "conditions": {}} for role in ROLES}
    if lane == "tmz":
        for detail in report.values():
            detail["mgmt_complete_origins"] = []
    incomplete = []
    for origin, (role, _, condition) in sorted(origins.items()):
        records = included[origin]
        if lane == "tmz":
            eligible = bool(records)
            extra_support = len(records) >= 2
        else:
            eligible = {row["region"] for row in records} == {"core", "infiltrating"}
            spatial_labels = {
                row["region"] for row in records
                if row.get("spatial_inputs") == row.get("independent_audit") == "verified"
            }
            extra_support = spatial_labels == {"core", "infiltrating"}
        if not eligible:
            incomplete.append(origin)
            continue
        role_report = report[role]
        role_report["eligible_origins"].append(origin)
        role_report["conditions"][condition] = role_report["conditions"].get(condition, 0) + 1
        if lane == "tmz" and all(row["mgmt"] == "verified" for row in records):
            role_report["mgmt_complete_origins"].append(origin)
        if extra_support:
            role_report["extra_support_origins"].append(origin)

    evaluation = len(report["evaluation"]["eligible_origins"])
    prior_conditions = set().union(*(report[role]["conditions"] for role in ROLES[:-1]))
    new_shift_origins = sum(
        count for condition, count in report["shift"]["conditions"].items()
        if condition not in prior_conditions
    )
    checks = {
        "development_origins": bool(report["development"]["eligible_origins"]),
        "evaluation_origin_floor": evaluation >= minimum_evaluation,
        "shift_origin_floor": len(report["shift"]["eligible_origins"]) >= minimum_shift,
        "held_out_shift_condition": new_shift_origins >= minimum_shift,
    }
    best_p = None
    if lane == "tmz":
        checks["independent_sister_cultures"] = (
            len(report["development"]["extra_support_origins"]) >= minimum_repeated_development
        )
    else:
        best_p = minimum_possible_adjusted_p(evaluation, comparison_count) if evaluation else 1.0
        checks["sign_test_attainability"] = best_p <= alpha
        if lane.endswith("-spatial"):
            checks["independent_calibration_pairs"] = bool(report["calibration"]["extra_support_origins"])

    return {
        "scope": "metadata only", "lane": lane, "checks": checks,
        "passed_metadata_checks": all(checks.values()), "roles": report,
        "minimum_possible_first_holm_p": best_p, "incomplete_origins": incomplete,
        "excluded_items": rejected,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--lane", choices=LANE_COLUMNS, required=True)
    parser.add_argument("--minimum-evaluation", type=int, required=True)
    parser.add_argument("--minimum-shift", type=int, required=True)
    parser.add_argument("--minimum-repeated-development", type=int)
    parser.add_argument("--comparison-count", type=int)
    parser.add_argument("--alpha", type=float)
    options = vars(parser.parse_args())
    manifest = options.pop("manifest")
    with manifest.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        expected = COMMON | LANE_COLUMNS[options["lane"]]
        if len(reader.fieldnames or ()) != len(expected) or set(reader.fieldnames or ()) != expected:
            parser.error("Manifest header must match the metadata-only schema exactly")
        try:
            report = audit(reader, **options)
        except ValueError as exc:
            parser.error(str(exc))
    print(json.dumps(report, indent=2))
    if not report["passed_metadata_checks"]:
        parser.exit(1)


if __name__ == "__main__":
    main()
