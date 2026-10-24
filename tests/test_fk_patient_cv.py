from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src.fk_patient_cv import (  # noqa: E402
    ForecastFailure,
    PatientSeries,
    run_outer_patient_cv,
)


def make_series(patient_id: str, *, n_scans: int = 5, empty_target: bool = False) -> PatientSeries:
    masks = []
    for scan_index in range(n_scans):
        mask = np.zeros((8, 8, 8), dtype=bool)
        mask[2:4, 2:4, 2:4] = True
        mask[2:4, 2:4, 4:4 + min(scan_index, 2)] = True
        masks.append(mask)
    if empty_target:
        masks[-1] = np.zeros_like(masks[-1])
    return PatientSeries(patient_id, np.arange(n_scans, dtype=float) * 30.0, tuple(masks))


def make_cohort():
    records = [make_series(f"patient-{index}") for index in range(6)]
    folds = {record.patient_id: index % 3 for index, record in enumerate(records)}
    return records, folds


def test_outer_patient_folds_keep_all_scans_of_test_patient_out_of_nuisance_fit():
    records, folds = make_cohort()
    fit_calls = []
    forecast_calls = []

    def fit_nuisance(train_records):
        ids = {record.patient_id for record in train_records}
        fit_calls.append(ids)
        return {"training_ids": ids}

    def forecast(history, nuisance):
        forecast_calls.append((history.patient_id, tuple(history.days), len(history.masks), nuisance["training_ids"]))
        return {"history_days": tuple(history.days), "history_scans": len(history.masks)}

    def score(prediction, target):
        return {"target_volume": int(target.mask.sum()), "history_scans": prediction["history_scans"]}

    result = run_outer_patient_cv(
        records, folds, fit_nuisance, forecast, score, n_folds=3, minimum_history=3,
    )

    assert len(result["rows"]) == 12
    assert len(result["fold_audit"]) == 3
    for call_ids in fit_calls:
        assert len(call_ids) == 4
    for row in result["rows"]:
        assert row["target_index"] in (3, 4)
        assert row["scores"]["history_scans"] == row["target_index"]
        assert row["forecast_status"] == "ok"
        assert row["patient_id"] not in fit_calls[row["fold"]]
    for patient_id, history_days, history_scans, nuisance_ids in forecast_calls:
        assert len(history_days) == history_scans
        assert patient_id not in nuisance_ids


def test_empty_heldout_target_is_scored_but_not_passed_to_forecast():
    records, folds = make_cohort()
    records[0] = make_series(records[0].patient_id, empty_target=True)
    forecast_target_ids = []

    def forecast(history, nuisance):
        forecast_target_ids.append(history.patient_id)
        assert not hasattr(history, "target")
        return {"history_scans": len(history.masks)}

    result = run_outer_patient_cv(
        records, folds, lambda train: {}, forecast,
        lambda prediction, target: {"target_volume": int(target.mask.sum())},
        n_folds=3, minimum_history=3,
    )

    assert len(forecast_target_ids) == 12
    empty_rows = [row for row in result["rows"] if row["patient_id"] == records[0].patient_id]
    assert empty_rows[-1]["scores"]["target_volume"] == 0


def test_forecast_failure_uses_history_only_fallback_and_is_counted():
    records, folds = make_cohort()
    fallback_calls = []

    def forecast(history, nuisance):
        if history.patient_id == "patient-0" and history.target_index == 3:
            raise ForecastFailure("synthetic numerical failure")
        return {"history_scans": len(history.masks)}

    def fallback(history, nuisance):
        fallback_calls.append((history.patient_id, history.target_index))
        return {"history_scans": len(history.masks), "fallback": True}

    result = run_outer_patient_cv(
        records, folds, lambda train: {}, forecast,
        lambda prediction, target: {"target_volume": int(target.mask.sum())},
        fallback_forecast=fallback, n_folds=3, minimum_history=3,
    )

    failed = [row for row in result["rows"] if row["forecast_status"] == "fallback"]
    assert len(failed) == len(fallback_calls) == 1
    assert failed[0]["patient_id"] == "patient-0"
    assert result["fold_audit"][0]["fallback_count"] == 1


def test_runner_rejects_duplicate_patients_missing_folds_and_bad_time_order():
    records, folds = make_cohort()
    with pytest.raises(ValueError, match="duplicate patient"):
        run_outer_patient_cv(records + [records[0]], folds, lambda train: {}, lambda h, n: {},
                             lambda p, t: {}, n_folds=3)
    incomplete = dict(folds)
    incomplete.pop("patient-0")
    with pytest.raises(ValueError, match="fold map"):
        run_outer_patient_cv(records, incomplete, lambda train: {}, lambda h, n: {},
                             lambda p, t: {}, n_folds=3)
    bad_days = PatientSeries(records[0].patient_id, np.asarray([0, 30, 20, 90, 120]), records[0].masks)
    with pytest.raises(ValueError, match="strictly increasing"):
        run_outer_patient_cv([bad_days, *records[1:]], folds, lambda train: {}, lambda h, n: {},
                             lambda p, t: {}, n_folds=3)


def test_runner_rejects_empty_historical_masks_and_missing_folds():
    records, folds = make_cohort()
    bad_masks = list(records[0].masks)
    bad_masks[2] = np.zeros_like(bad_masks[2])
    bad_record = PatientSeries(records[0].patient_id, records[0].days, tuple(bad_masks))
    with pytest.raises(ValueError, match="history masks must be nonempty"):
        run_outer_patient_cv([bad_record, *records[1:]], folds, lambda train: {}, lambda h, n: {},
                             lambda p, t: {}, n_folds=3)
    with pytest.raises(ValueError, match="every outer fold"):
        run_outer_patient_cv(records, {patient: 0 for patient in folds}, lambda train: {}, lambda h, n: {},
                             lambda p, t: {}, n_folds=3)


def test_runner_rejects_bad_records_masks_and_parameters():
    records, folds = make_cohort()
    with pytest.raises(ValueError, match="records must contain PatientSeries"):
        run_outer_patient_cv([{"patient_id": "p"}], folds, lambda train: {}, lambda h, n: {},
                             lambda p, t: {}, n_folds=3)
    with pytest.raises(ValueError, match="at least one patient"):
        run_outer_patient_cv([], {}, lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=3)
    with pytest.raises(ValueError, match="patient IDs must be nonempty"):
        run_outer_patient_cv([PatientSeries(" ", records[0].days, records[0].masks)], {},
                             lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=3)
    bad_binary = list(records[0].masks)
    bad_binary[0] = np.full((8, 8, 8), 2, dtype=np.uint8)
    with pytest.raises(ValueError, match="masks must be binary"):
        run_outer_patient_cv([PatientSeries(records[0].patient_id, records[0].days, tuple(bad_binary)),
                              *records[1:]], folds, lambda train: {}, lambda h, n: {}, lambda p, t: {},
                             n_folds=3)
    bad_shape = list(records[0].masks)
    bad_shape[0] = np.zeros((7, 8, 8), dtype=bool)
    bad_shape[0][2:4, 2:4, 2:4] = True
    with pytest.raises(ValueError, match="share one shape"):
        run_outer_patient_cv([PatientSeries(records[0].patient_id, records[0].days, tuple(bad_shape)),
                              *records[1:]], folds, lambda train: {}, lambda h, n: {}, lambda p, t: {},
                             n_folds=3)
    with pytest.raises(ValueError, match="n_folds"):
        run_outer_patient_cv(records, folds, lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=1)
    with pytest.raises(ValueError, match="minimum_history"):
        run_outer_patient_cv(records, folds, lambda train: {}, lambda h, n: {}, lambda p, t: {},
                             n_folds=3, minimum_history=2)


def test_fold_map_rejects_noninteger_boolean_out_of_range_and_normalized_duplicates():
    records, folds = make_cohort()
    invalid = dict(folds)
    invalid["patient-0"] = True
    with pytest.raises(ValueError, match="not booleans"):
        run_outer_patient_cv(records, invalid, lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=3)
    invalid = dict(folds)
    invalid["patient-0"] = "bad"
    with pytest.raises(ValueError, match="folds must be integers"):
        run_outer_patient_cv(records, invalid, lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=3)
    invalid = dict(folds)
    invalid["patient-0"] = 3
    with pytest.raises(ValueError, match="outside the declared"):
        run_outer_patient_cv(records, invalid, lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=3)
    invalid["patient-0"] = 0
    invalid[" patient-0"] = 0
    with pytest.raises(ValueError, match="duplicate normalized"):
        run_outer_patient_cv(records, invalid, lambda train: {}, lambda h, n: {}, lambda p, t: {}, n_folds=3)


def test_missing_fallback_and_invalid_callback_outputs_fail_loudly():
    records, folds = make_cohort()

    def failing_forecast(history, nuisance):
        raise ForecastFailure("expected synthetic failure")

    with pytest.raises(RuntimeError, match="no history-only fallback"):
        run_outer_patient_cv(records, folds, lambda train: {}, failing_forecast,
                             lambda prediction, target: {}, n_folds=3)
    with pytest.raises(ValueError, match="return a prediction"):
        run_outer_patient_cv(records, folds, lambda train: {}, lambda history, nuisance: None,
                             lambda prediction, target: {}, n_folds=3)
    with pytest.raises(ValueError, match="return a mapping"):
        run_outer_patient_cv(records, folds, lambda train: {}, lambda history, nuisance: {},
                             lambda prediction, target: [], n_folds=3)
