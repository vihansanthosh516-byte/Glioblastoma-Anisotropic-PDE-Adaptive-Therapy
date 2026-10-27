"""Metadata-only regression tests for Track B forecast-triplet selection."""
from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "162_growth_zone_cv.py"
SPEC = importlib.util.spec_from_file_location("track_b_manifest_triplets_test", MODULE_PATH)
CV = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = CV
assert SPEC.loader is not None
SPEC.loader.exec_module(CV)


def _write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _patient(patient_id: str, *, eligible: bool = True, is_gbm: bool = True) -> dict[str, str]:
    return {
        "patient_id": patient_id,
        "dataset": "MU",
        "eligible": str(eligible).lower(),
        "is_gbm": str(is_gbm).lower(),
    }


def _pair(
    patient_id: str,
    tp_in: int,
    tp_out: int,
    day_in: float,
    day_out: float,
    *,
    in_primary: bool = True,
    input_core_nonempty: bool = True,
) -> dict[str, object]:
    return {
        "patient_id": patient_id,
        "tp_in": tp_in,
        "tp_out": tp_out,
        "day_in": day_in,
        "day_out": day_out,
        "in_primary": str(in_primary).lower(),
        "input_core_nonempty": str(input_core_nonempty).lower(),
    }


def test_manifest_triplets_apply_patient_pair_and_chain_filters(tmp_path: Path) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    _write_csv(
        patient_path,
        ["patient_id", "dataset", "eligible", "is_gbm"],
        [
            _patient("eligible_gbm"),
            _patient("eligible_non_gbm", is_gbm=False),
            _patient("ineligible_gbm", eligible=False),
            _patient("unchained_gbm"),
        ],
    )
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "in_primary", "input_core_nonempty"],
        [
            # Valid consecutive chain. Day mismatch is within the source tolerance.
            _pair("eligible_gbm", 1, 2, 0, 30),
            _pair("eligible_gbm", 2, 3, 30.0005, 60),
            # Consecutive timepoints, but day mismatch exceeds 1e-3 days.
            _pair("eligible_gbm", 3, 4, 60.002, 90),
            # These rows fail the pair-level filters.
            _pair("eligible_gbm", 4, 5, 90, 120, in_primary=False),
            _pair("eligible_gbm", 5, 6, 120, 150, input_core_nonempty=False),
            # Both rows chain, but patient-level eligibility/diagnosis excludes them.
            _pair("eligible_non_gbm", 1, 2, 0, 30),
            _pair("eligible_non_gbm", 2, 3, 30, 60),
            _pair("ineligible_gbm", 1, 2, 0, 30),
            _pair("ineligible_gbm", 2, 3, 30, 60),
            # Eligible GBM rows do not chain by timepoint.
            _pair("unchained_gbm", 1, 2, 0, 30),
            _pair("unchained_gbm", 3, 4, 30, 60),
        ],
    )

    triplets = CV._manifest_triplets(patient_path, pair_path)

    assert len(triplets) == 1
    patient_id, first, second = triplets[0]
    assert patient_id == "eligible_gbm"
    assert (int(first["tp_in"]), int(first["tp_out"])) == (1, 2)
    assert (int(second["tp_in"]), int(second["tp_out"])) == (2, 3)


def test_manifest_triplets_fail_when_no_valid_chain_exists(tmp_path: Path) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    _write_csv(
        patient_path,
        ["patient_id", "dataset", "eligible", "is_gbm"],
        [_patient("gbm")],
    )
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "in_primary", "input_core_nonempty"],
        [
            _pair("gbm", 1, 2, 0, 30),
            _pair("gbm", 2, 3, 31, 60),  # sequential scan IDs, but inconsistent visit day
            _pair("gbm", 4, 5, 60, 90),  # neither row chains to the previous pair
        ],
    )

    with pytest.raises(ValueError, match="no chained, eligible GBM scan triplets"):
        CV._manifest_triplets(patient_path, pair_path)
