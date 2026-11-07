"""Synthetic tests for source-agnostic Track B prep. No patient data."""
import inspect
from types import SimpleNamespace

import numpy as np
import pytest

from src import fk_longitudinal_prep as lp
from src.track_b_transition_estimand import geometric_coverage


def _v(pid, week, usable=True, order=0):
    return lp.Visit(pid, week, order, usable, usable, usable)


def test_primary_triplet_uses_first_three_usable_postop_visits_only():
    visits = [_v("a", 0), _v("a", 4), _v("a", 8), _v("a", 21), _v("a", 34),  # first gap 13: eligible
              _v("b", 4), _v("b", 8), _v("b", 30), _v("b", 43),              # first gap 22: excluded, never re-searched
              _v("c", 4), _v("c", 8, usable=False), _v("c", 12), _v("c", 25),  # unusable visit skipped
              _v("d", 4), _v("d", 8)]
    out = lp.primary_triplets(visits, min_gap_weeks=11, max_gap_weeks=15)
    assert [v.week for v in out["eligible"]["a"]] == [4, 8, 21]
    assert "b" in out["excluded"] and "22 weeks" in out["excluded"]["b"]
    assert [v.week for v in out["eligible"]["c"]] == [4, 12, 25]
    assert out["excluded"]["d"] == "fewer than 3 usable post-op visits"


def test_origin_support_has_no_target_input_and_keeps_history_core():
    assert list(inspect.signature(lp.origin_support).parameters) == ["origin_brain_mask", "history_label_maps"]
    brain = np.zeros((10, 10, 10), bool)
    brain[2:8, 2:8, 2:8] = True
    hist = np.zeros((10, 10, 10), np.uint8)
    hist[8, 5, 5] = 2  # history core just outside the brain mask
    s = lp.origin_support(brain, [hist])
    assert s[brain].all() and s[8, 5, 5] and s.sum() == brain.sum() + 1
    with pytest.raises(ValueError, match="empty"):
        lp.origin_support(np.zeros((4, 4, 4), bool))


def _phantom(shape=(40, 40, 32)):
    g = np.indices(shape).astype(float)
    c = (np.array(shape) - 1) / 2
    r = np.sqrt(sum(((g[i] - c[i]) / (shape[i] / 2.6)) ** 2 for i in range(3)))
    vol = (r < 1).astype(float) * (1 + 0.5 * np.sin(g[0] / 3) + 0.3 * np.cos(g[1] / 4) + 0.2 * (g[2] / shape[2]))
    return vol, r < 1.05


def test_registration_recovers_known_rigid_motion_and_passes_qc():
    fixed, sup = _phantom()
    truth = np.array([3.0, -2.0, 4.0, 4.0, -2.0, 2.0])
    moving = lp.apply_rigid(fixed, truth, spacing_mm=(2, 2, 2))
    out = lp.register_rigid_mri(fixed, moving, fixed_support=sup, spacing_mm=(2, 2, 2))
    assert out["ncc_after"] > out["ncc_before"]
    assert out["ncc_after"] > 0.95 and out["qc"]["passed"]


def test_registration_qc_fails_closed_and_reports_exclusions():
    bad = lp.registration_qc(np.array([0, 0, 25.0, 0, 0, 0]), 0.5)
    assert not bad["passed"] and len(bad["reasons"]) == 2
    good = lp.registration_qc(np.zeros(6), 0.9)
    rep = lp.exclusion_report({"p1": good, "p2": bad})
    assert rep == {"n_checked": 2, "n_excluded": 1, "excluded": {"p2": bad["reasons"]}}
    assert lp.registration_qc(np.zeros(6), float("nan"))["passed"] is False


def test_registration_qc_rejects_nonfinite_or_malformed_transform():
    for params in ([float('nan'), 0, 0, 0, 0, 0],
                   [0, 0, 0, float('nan'), 0, 0],
                   [0, 0, 0, float('inf'), 0, 0],
                   [0, 0, 0]):
        result = lp.registration_qc(np.asarray(params), 0.95)
        assert result['passed'] is False
        assert result['reasons']


@pytest.mark.parametrize('failure_stage', [0, 1])
def test_registration_optimizer_failure_never_passes_qc(monkeypatch, failure_stage):
    calls = []
    def fake_minimize(*args, **kwargs):
        stage = len(calls)
        calls.append(stage)
        return SimpleNamespace(x=np.zeros(6), success=stage != failure_stage)
    monkeypatch.setattr(lp.optimize, 'minimize', fake_minimize)
    fixed, support = _phantom((8, 8, 8))
    result = lp.register_rigid_mri(fixed, fixed, fixed_support=support, spacing_mm=(2, 2, 2))
    assert result['qc']['passed'] is False
    assert any('optimizer' in reason for reason in result['qc']['reasons'])


def test_registration_uses_no_label_input():
    params = list(inspect.signature(lp.register_rigid_mri).parameters)
    assert params == ["fixed", "moving", "fixed_support", "spacing_mm", "downsample"]


@pytest.mark.parametrize('params,spacing', [
    ([0, 0, 0, 0, 0, 0], (2, 2, 2)),
    ([0, 0, 0, 2, -2, 0], (2, 2, 2)),
    ([4, -6, 8, 1.2, -2.1, 0.4], (1, 2, 3)),
])
def test_coverage_matches_existing_rigid_label_resampler(params, spacing):
    shape = (12, 10, 8)
    center = (np.array(shape) - 1) / 2
    rotation = lp._rotation(params[:3])
    scale = np.diag(spacing)
    matrix = np.linalg.inv(scale) @ rotation @ scale
    offset = center + np.array(params[3:]) / spacing - matrix @ center
    affine = np.eye(4)
    affine[:3, :3], affine[:3, 3] = matrix, offset
    coverage = geometric_coverage(shape, shape, affine)
    reference = lp.apply_rigid(np.ones(shape), params, spacing_mm=spacing, order=0).astype(bool)
    np.testing.assert_array_equal(coverage, reference)
