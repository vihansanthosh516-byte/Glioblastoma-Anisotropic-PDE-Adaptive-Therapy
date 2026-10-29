"""Patient-level bootstrap precision gate for the Track B external-test plan.

This estimates study resolution from development-cohort paired differences. It is not an effect estimate for
the external cohort and must receive exactly one paired difference per development patient.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
import hashlib

import numpy as np

BOOTSTRAP_REPS = 10_000
BOOTSTRAP_SEED = 20261028
MAX_HALF_WIDTH_DICE = 0.05


def precision_gate_from_patient_differences(
    patient_differences: Iterable[tuple[str, float]] | Mapping[str, float], *,
    expected_patient_ids: Iterable[str],
    target_n: int,
) -> dict:
    """Bootstrap paired mean Dice difference at a planned external patient count.

    Rows must contain exactly one finite paired value for every expected patient. Patients are resampled with
    replacement; visit/pair rows are never treated as independent observations.
    """
    if isinstance(target_n, (bool, np.bool_)) or int(target_n) != target_n or target_n < 1:
        raise ValueError("target_n must be a positive integer")
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
