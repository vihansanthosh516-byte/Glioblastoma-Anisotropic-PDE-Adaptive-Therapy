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
    dt_days: float | None = None,
    in_primary: bool = True,
    input_core_nonempty: bool = True,
) -> dict[str, object]:
    return {
        "patient_id": patient_id,
        "tp_in": tp_in,
        "tp_out": tp_out,
        "day_in": day_in,
        "day_out": day_out,
        "dt_days": day_out - day_in if dt_days is None else dt_days,
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
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
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
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [
            _pair("gbm", 1, 2, 0, 30),
            _pair("gbm", 2, 3, 31, 60),  # sequential scan IDs, but inconsistent visit day
            _pair("gbm", 4, 5, 60, 90),  # neither row chains to the previous pair
        ],
    )

    with pytest.raises(ValueError, match="no chained, eligible GBM scan triplets"):
        CV._manifest_triplets(patient_path, pair_path)


def test_manifest_rejects_duplicate_patient_ids(tmp_path: Path) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    _write_csv(
        patient_path,
        ["patient_id", "dataset", "eligible", "is_gbm"],
        [_patient("P1"), _patient(" P1 ")],
    )
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [_pair("P1", 1, 2, 0, 30), _pair("P1", 2, 3, 30, 60)],
    )

    with pytest.raises(ValueError, match="duplicate patient_id"):
        CV._manifest_triplets(patient_path, pair_path)


def test_manifest_rejects_duplicate_pairs_and_unknown_patient_ids(tmp_path: Path) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    _write_csv(patient_path, ["patient_id", "dataset", "eligible", "is_gbm"], [_patient("P1")])
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [_pair("P1", 1, 2, 0, 30), _pair("P1", 1, 2, 0, 30)],
    )
    with pytest.raises(ValueError, match="duplicate patient/timepoint pair"):
        CV._manifest_triplets(patient_path, pair_path)

    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [_pair("UNKNOWN", 1, 2, 0, 30)],
    )
    with pytest.raises(ValueError, match="absent from patient manifest"):
        CV._manifest_triplets(patient_path, pair_path)


@pytest.mark.parametrize(
    ("pair", "message"),
    [
        (_pair("P1", 1, 2, 0, 30, dt_days=29), "dt_days disagrees"),
        (_pair("P1", 1, 2, 0, 30, dt_days=float("nan")), "finite numeric"),
        (_pair("P1", 1, 2, 0, 30, dt_days=0), "interval must be positive"),
        (_pair("P1", 1.5, 2, 0, 30), "nonnegative integer"),
    ],
)
def test_manifest_rejects_invalid_pair_intervals_and_timepoints(
    tmp_path: Path, pair: dict[str, object], message: str
) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    _write_csv(patient_path, ["patient_id", "dataset", "eligible", "is_gbm"], [_patient("P1")])
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [pair],
    )

    with pytest.raises(ValueError, match=message):
        CV._manifest_triplets(patient_path, pair_path)


def test_manifest_rejects_ambiguous_flags_and_missing_required_columns(tmp_path: Path) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    _write_csv(
        patient_path,
        ["patient_id", "dataset", "eligible", "is_gbm"],
        [{"patient_id": "P1", "dataset": "MU", "eligible": "sometimes", "is_gbm": "true"}],
    )
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [_pair("P1", 1, 2, 0, 30)],
    )
    with pytest.raises(ValueError, match="boolean must be explicit"):
        CV._manifest_triplets(patient_path, pair_path)

    _write_csv(patient_path, ["patient_id", "eligible"], [{"patient_id": "P1", "eligible": "true"}])
    with pytest.raises(ValueError, match="lacks required columns"):
        CV._manifest_triplets(patient_path, pair_path)


def test_run_preflight_rejects_existing_output_before_manifest_or_image_access(tmp_path: Path, monkeypatch) -> None:
    output_path = tmp_path / "already-there.json"
    output_path.write_text("synthetic sentinel", encoding="utf-8")
    args = type("Args", (), {
        "outer_folds": 2,
        "inner_folds": 2,
        "ridge_grid": [0.1],
        "output": output_path,
        "patient_manifest": tmp_path / "missing-patients.csv",
        "pair_manifest": tmp_path / "missing-pairs.csv",
        "data_root": tmp_path / "missing-data",
        "overwrite": False,
    })()

    def forbidden(*_args, **_kwargs):
        pytest.fail("patient image loader must not run during preflight")

    monkeypatch.setattr(CV, "_load_raw_samples", forbidden)
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        CV.run(args)


def test_run_preflight_rejects_output_inside_data_root_even_with_overwrite(tmp_path: Path, monkeypatch) -> None:
    data_root = tmp_path / "data"
    data_root.mkdir()
    output_path = data_root / "result.json"
    output_path.write_text("synthetic sentinel", encoding="utf-8")
    args = type("Args", (), {
        "outer_folds": 2,
        "inner_folds": 2,
        "ridge_grid": [0.1],
        "output": output_path,
        "patient_manifest": tmp_path / "missing-patients.csv",
        "pair_manifest": tmp_path / "missing-pairs.csv",
        "data_root": data_root,
        "overwrite": True,
    })()

    def forbidden(*_args, **_kwargs):
        pytest.fail("patient image loader must not run during preflight")

    monkeypatch.setattr(CV, "_load_raw_samples", forbidden)
    with pytest.raises(ValueError, match="output path must not be inside the data root"):
        CV.run(args)


def test_run_preflight_checks_patient_count_before_image_access(tmp_path: Path, monkeypatch) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    data_root = tmp_path / "data"
    data_root.mkdir()
    output_path = tmp_path / "result.json"
    _write_csv(patient_path, ["patient_id", "dataset", "eligible", "is_gbm"], [_patient("P1")])
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [_pair("P1", 1, 2, 0, 30), _pair("P1", 2, 3, 30, 60)],
    )
    args = type("Args", (), {
        "outer_folds": 2,
        "inner_folds": 2,
        "ridge_grid": [0.1],
        "output": output_path,
        "patient_manifest": patient_path,
        "pair_manifest": pair_path,
        "data_root": data_root,
        "overwrite": False,
        "edge_sensitivity": 0.5,
        "screening_length_mm": 10.0,
        "label2_log_permeability": 0.69314718056,
        "candidate_radius_mm": 20.0,
        "prior_smoothing_mm": 2.0,
        "open_cavity_ablation": False,
    })()

    def forbidden(*_args, **_kwargs):
        pytest.fail("patient image loader must not run during preflight")

    monkeypatch.setattr(CV, "_load_raw_samples", forbidden)
    with pytest.raises(ValueError, match="fewer unique patients than outer folds"):
        CV.run(args)


def test_run_preflight_rejects_bad_interval_before_image_access(tmp_path: Path, monkeypatch) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    data_root = tmp_path / "data"
    data_root.mkdir()
    _write_csv(
        patient_path,
        ["patient_id", "dataset", "eligible", "is_gbm"],
        [_patient("P1"), _patient("P2")],
    )
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [
            _pair("P1", 1, 2, 0, 30),
            _pair("P1", 2, 3, 30, 60),
            _pair("P2", 1, 2, 0, 30, dt_days=29),
            _pair("P2", 2, 3, 30, 60),
        ],
    )
    args = type("Args", (), {
        "outer_folds": 2,
        "inner_folds": 2,
        "ridge_grid": [0.1],
        "output": tmp_path / "result.json",
        "patient_manifest": patient_path,
        "pair_manifest": pair_path,
        "data_root": data_root,
        "overwrite": False,
        "edge_sensitivity": 0.5,
        "screening_length_mm": 10.0,
        "label2_log_permeability": 0.69314718056,
        "candidate_radius_mm": 20.0,
        "prior_smoothing_mm": 2.0,
        "open_cavity_ablation": False,
    })()

    def forbidden(*_args, **_kwargs):
        pytest.fail("patient image loader must not run during preflight")

    monkeypatch.setattr(CV, "_load_raw_samples", forbidden)
    with pytest.raises(ValueError, match="dt_days disagrees"):
        CV.run(args)


def test_run_preflight_requires_inner_training_patients_before_image_access(tmp_path: Path, monkeypatch) -> None:
    patient_path = tmp_path / "patients.csv"
    pair_path = tmp_path / "pairs.csv"
    data_root = tmp_path / "data"
    data_root.mkdir()
    _write_csv(
        patient_path,
        ["patient_id", "dataset", "eligible", "is_gbm"],
        [_patient("P1"), _patient("P2")],
    )
    _write_csv(
        pair_path,
        ["patient_id", "tp_in", "tp_out", "day_in", "day_out", "dt_days", "in_primary", "input_core_nonempty"],
        [
            _pair("P1", 1, 2, 0, 30),
            _pair("P1", 2, 3, 30, 60),
            _pair("P2", 1, 2, 0, 30),
            _pair("P2", 2, 3, 30, 60),
        ],
    )
    args = type("Args", (), {
        "outer_folds": 2,
        "inner_folds": 2,
        "ridge_grid": [0.1],
        "output": tmp_path / "result.json",
        "patient_manifest": patient_path,
        "pair_manifest": pair_path,
        "data_root": data_root,
        "overwrite": False,
        "edge_sensitivity": 0.5,
        "screening_length_mm": 10.0,
        "label2_log_permeability": 0.69314718056,
        "candidate_radius_mm": 20.0,
        "prior_smoothing_mm": 2.0,
        "open_cavity_ablation": False,
    })()

    def forbidden(*_args, **_kwargs):
        pytest.fail("patient image loader must not run during preflight")

    monkeypatch.setattr(CV, "_load_raw_samples", forbidden)
    with pytest.raises(ValueError, match="outer training fold has too few patients"):
        CV.run(args)
