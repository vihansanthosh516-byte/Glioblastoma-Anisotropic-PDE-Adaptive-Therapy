"""End-to-end synthetic test: coverage-aware scoring inside the LUMIERE forecast runner. No patient data."""
from dataclasses import replace

import numpy as np
import pytest

from src import fk_lumiere_runner as lr
from src.track_b_coverage_scoring import coverage_from_geometry
from tests.test_fk_lumiere_runner import SHAPE, TINY, FakeProvider, _ball, _cases, _folds


def _cov(shift_target=4, visits=("visit1", "visit2", "target")):
    out = []
    for v in visits:
        m = np.eye(4)
        if v == "target":
            m[0, 3] = shift_target          # target scan misses the last `shift_target` slices
        out.append(coverage_from_geometry(v, SHAPE, SHAPE, m))
    return tuple(out)


def _run(cases, **kw):
    return lr.run_lumiere_cv(cases, _folds(cases), pde_provider=FakeProvider(), recipe=TINY,
                             minimum=10, n_bootstrap=50, **kw)


def _with_outside_lesion(cases, loads=None):
    """Same targets inside coverage; an extra lesion only in uncovered slices 12-15."""
    out = []
    for c in cases:
        def load(c=c):
            t = np.asarray(c.load_target(), bool).copy()
            t[13:16, 1:4, 1:4] = True
            return t
        out.append(replace(c, load_target=load))
    return out


def test_coverage_applied_reported_and_identical_for_all_methods():
    cases = [replace(c, coverage=_cov()) for c in _cases(12)]
    res = _run(cases, require_coverage=True)
    cs = res["coverage_scoring"]
    assert cs["applied"] and cs["visits"] == ["visit1", "visit2", "target"]
    expected_voxels = 12 * 16 * 16
    for pid, rep in cs["per_patient"].items():
        assert rep["support_voxels"] == expected_voxels
        assert rep["support_fraction"] == pytest.approx(12 / 16)
    # one support per patient is shared by every scored method
    assert set(cs["per_patient"]) == {c.patient_id for c in cases}


def test_labels_outside_support_cannot_change_any_score():
    base_cases = [replace(c, coverage=_cov()) for c in _cases(12)]
    base = _run(base_cases, require_coverage=True)
    altered = _run(_with_outside_lesion(base_cases), require_coverage=True)
    for method, per_patient in base["dice"].items():
        assert per_patient == altered["dice"][method], method


def test_without_coverage_the_outside_lesion_does_change_scores():
    """Contrast: historical (no-coverage) scoring counts the uncovered lesion as missed."""
    cases = _cases(12)
    base = _run(cases)
    altered = _run(_with_outside_lesion(cases))
    assert base["coverage_scoring"] == {"applied": False}
    assert base["dice"]["copy_last"] != altered["dice"]["copy_last"]


def test_missing_coverage_stops_before_any_target_read():
    loads = []
    cases = [replace(c, coverage=_cov()) for c in _cases(12, loads=loads)]
    cases[3] = replace(cases[3], coverage=())
    with pytest.raises(PermissionError, match="no scan coverage"):
        _run(cases, require_coverage=True)
    assert loads == []


def test_invalid_coverage_stops_scoring():
    loads = []
    cases = [replace(c, coverage=_cov()) for c in _cases(12, loads=loads)]
    cases[0] = replace(cases[0], coverage=_cov(visits=("visit1", "visit2")))           # target missing
    with pytest.raises(ValueError, match="no scan coverage for required visits"):
        _run(cases, require_coverage=True)
    cases[0] = replace(cases[0], coverage=_cov() + (_cov()[0],))                        # duplicate visit
    with pytest.raises(ValueError, match="duplicate coverage"):
        _run(cases, require_coverage=True)
    cases[0] = replace(cases[0], coverage=(np.ones(SHAPE, bool),) * 3)                  # raw masks
    with pytest.raises(TypeError):
        _run(cases, require_coverage=True)
    assert loads == []
