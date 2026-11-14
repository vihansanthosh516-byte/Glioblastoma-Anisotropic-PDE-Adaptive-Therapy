"""Validate the metadata needed before opening Track A expression."""

from collections import defaultdict
from collections.abc import Iterable, Mapping

from .trackA_margin_evaluation import minimum_possible_adjusted_p


def audit_specimens(
    rows: Iterable[Mapping[str, str]], *, comparison_count: int, minimum_donors: int = 6,
    alpha: float = 0.05,
) -> dict:
    """Check a curated manifest of eligible specimens using study-global origin IDs.

    The curator resolves origin IDs and verifies independent region provenance
    against source evidence. This function checks the declared metadata as given.
    """
    if minimum_donors < 1 or not 0 < alpha < 1 or comparison_count < 1:
        raise ValueError("Invalid benchmark thresholds")
    specimens = set()
    origins = {}
    study_roles = {}
    regions = defaultdict(set)
    for row in rows:
        study = row["study"].strip()
        specimen = row["specimen"].strip()
        origin = row["origin"].strip()
        role = row["role"].strip()
        region = row["region"].strip()
        if not study or not specimen or not origin or not row["region_source"].strip():
            raise ValueError("Study, specimen, origin and region source are required")
        if role not in {"development", "test"} or region not in {"core", "margin"}:
            raise ValueError("Invalid role or region")
        if row["label_basis"] not in {"pathology", "surgical_sampling"}:
            raise ValueError("Region labels must come from independent sampling or pathology")
        if (study, specimen) in specimens:
            raise ValueError(f"Duplicate specimen in study {study}")
        specimens.add((study, specimen))
        if study_roles.setdefault(study, role) != role:
            raise ValueError(f"Study {study} crosses development and test")
        if origins.setdefault(origin, (study, role)) != (study, role):
            raise ValueError(f"Origin {origin} is repeated across studies or roles")
        regions[origin].add(region)

    if list(study_roles.values()).count("test") != 1:
        raise ValueError("The current protocol requires one external test study")
    if any(values != {"core", "margin"} for values in regions.values()):
        raise ValueError("Every eligible origin must have both regions")
    test = sorted(origin for origin, (_, role) in origins.items() if role == "test")
    development = sorted(origin for origin, (_, role) in origins.items() if role == "development")
    if not development:
        raise ValueError("Development origins are required")
    best_p = minimum_possible_adjusted_p(len(test), comparison_count)
    return {
        "test_origins": test,
        "development_origins": development,
        "minimum_possible_first_holm_p": best_p,
        "count_gate_passed": len(test) >= minimum_donors,
        "attainability_gate_passed": best_p <= alpha,
    }
