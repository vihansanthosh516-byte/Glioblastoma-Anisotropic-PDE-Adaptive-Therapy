"""Synthetic tests for the LUMIERE result writer. No patient data."""
import hashlib
import json
from dataclasses import replace

import pytest

from src import fk_lumiere_protocol as lp
from src import fk_lumiere_result_writer as rw
from src import fk_lumiere_runner as lr
from tests.test_fk_lumiere_protocol import _factory, _spread_cases
from tests.test_fk_lumiere_runner import TINY, FakeProvider, _folds


def _protocol(cases):
    return lp.run_lumiere_protocol(cases, _folds(cases), provider_factory=_factory([]), recipe=TINY, n_bootstrap=200)


def test_written_record_reads_back_with_all_fields(tmp_path):
    cases = _spread_cases()
    result = _protocol(cases)
    path = tmp_path / "out" / "lumiere_result.json"
    rw.write_lumiere_protocol_result(result, path, cases=cases, fold_by_patient=_folds(cases))
    rec = json.loads(path.read_text())
    assert rec["label"] == "synthetic protocol implementation" and "NOT adopted" in rec["amendment_status"]
    assert rec["overall_status"] == "completed (synthetic)"
    for f, h in rec["code_sha256"].items():
        assert h == hashlib.sha256((rw.PROJECT_ROOT / f).read_bytes()).hexdigest()
    assert all(rec["protocol_sha256"][f] for f in rw.PROTOCOL_FILES)
    assert rec["cohort"]["n_input"] == 26 and len(rec["cohort"]["fold_by_patient"]) == 26
    assert rec["gates"] == {g: False for g in lr.LUMIERE_GATES}
    assert set(rec["analyses"]) == set(rw.ANALYSES)
    for name, a in rec["analyses"].items():
        assert a["status"] == "scored (synthetic)"
        assert set(a["patient_dice"]) >= set(lr.COMPARATORS) | {"candidate", "matched_volume_no_pde_ablation"}
        assert a["pass_test"]["status"] in {"pass", "not passed", "incomplete"}
        fc = a["fallback_counts"]
        assert isinstance(fc["pde_forecast_failure_fallbacks"], int) and fc["folds_pde_not_ready"] == []
        if any(v is None for v in a["patient_dice"]["film_unet"].values()):
            assert a["missing_arm_reasons"]["film_unet"]  # the reason is recorded
        assert "count_arm" in a and "mean_poisson_deviance" in a["count_arm"]
    assert rec["settings"]["minimums"] == {"primary": 25, "central_gap": 10}
    assert rec["settings"]["seeds"] == {"forecast": 20261028, "boundary_calibration": 20261029, "bootstrap": 20261110}
    assert set(rec["analyses"]["central_gap"]["patient_dice"]["candidate"]) == {
        c.patient_id for c in cases if c.forecast_weeks in lp.CENTRAL_GAP_WEEKS}
    with pytest.raises(FileExistsError):
        rw.write_lumiere_protocol_result(result, path, cases=cases, fold_by_patient=_folds(cases))


def test_fallbacks_and_incomplete_protocol_are_recorded(tmp_path):
    cases = _spread_cases()
    def factory(rho_b, d_b):
        return FakeProvider(fail_ids={"s07"})
    result = lp.run_lumiere_protocol(cases, _folds(cases), provider_factory=factory, recipe=TINY, n_bootstrap=200)
    rec = rw.write_lumiere_protocol_result(result, tmp_path / "a.json", cases=cases, fold_by_patient=_folds(cases))
    assert rec["analyses"]["primary"]["fallback_counts"]["pde_forecast_failure_fallbacks"] == 1
    few = [replace(c, forecast_weeks=11) if c.forecast_weeks in lp.CENTRAL_GAP_WEEKS and i > 12 else c
           for i, c in enumerate(cases)]
    r2 = _protocol(few)
    rec2 = rw.write_lumiere_protocol_result(r2, tmp_path / "b.json", cases=few, fold_by_patient=_folds(few))
    assert rec2["overall_status"].startswith("incomplete") and rec2["analyses"]["central_gap"]["flow"]["stopped"]
    assert "patient_dice" not in rec2["analyses"]["central_gap"]
