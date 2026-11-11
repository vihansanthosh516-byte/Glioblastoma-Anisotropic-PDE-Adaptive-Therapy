"""Synthetic tests for coverage-aware Track B scoring. No patient data."""
import numpy as np
import pytest

from src.track_b_coverage_scoring import (ScanCoverage, common_support, coverage_from_geometry,
                                          forecast_dice_on_support, transition_window_on_support)
from src.track_b_transition_estimand import construct_history_states

SHAPE = (20, 20, 20)
SP = (2.0, 2.0, 2.0)


def shifted(visit, shift_z):
    """Coverage of a moving grid shifted by shift_z voxels along axis 0 (partial overlap)."""
    m = np.eye(4)
    m[0, 3] = shift_z
    return coverage_from_geometry(visit, SHAPE, SHAPE, m)


def history():
    previous = np.zeros(SHAPE, int)
    middle = np.zeros(SHAPE, int)
    middle[9:11, 9:11, 9:11] = 1            # core
    middle[:, 2:18, 2:18][middle[:, 2:18, 2:18] == 0] = 2   # FLAIR around
    previous[:, 2:10, 2:18] = 2              # previous FLAIR covers half -> persistent vs new
    return construct_history_states(previous, middle, voxel_spacing_mm=SP)


def test_partial_coverage_and_multi_visit_and():
    covs = [shifted("previous", 0), shifted("middle", 3), shifted("target", -4)]
    support = common_support(covs, ("previous", "middle", "target"))
    # shift +3 loses the last 3 slices, shift -4 loses the first 4
    assert support[:4].sum() == 0 and support[-3:].sum() == 0 and support[4:17].all()


def test_fail_closed_on_missing_or_invalid_coverage():
    covs = [shifted("previous", 0), shifted("middle", 0)]
    with pytest.raises(ValueError, match="no scan coverage"):
        common_support(covs, ("previous", "middle", "target"))
    with pytest.raises(TypeError):
        ScanCoverage(mask=np.ones(SHAPE, bool), visit="target")       # not from geometry
    with pytest.raises(TypeError):
        common_support([np.ones(SHAPE, bool)], ("target",))
    with pytest.raises(ValueError, match="empty"):
        common_support([shifted("a", 0), shifted("b", 25)], ("a", "b"))


def test_labels_outside_support_cannot_change_transition_score():
    state, bands = history()
    support = common_support([shifted("previous", 0), shifted("middle", 0), shifted("target", 6)],
                             ("previous", "middle", "target"))
    target = np.zeros(SHAPE, int)
    target[8:12, 8:12, 8:12] = 1
    base = transition_window_on_support(state, bands, target, support, min_voxels=5)
    assert base is not None
    altered = target.copy()
    altered[~support] = np.random.default_rng(1).choice([0, 1, 2, 3, 4], size=int((~support).sum()))
    again = transition_window_on_support(state, bands, altered, support, min_voxels=5)
    assert again["difference"] == base["difference"] and again["bands"] == base["bands"]
    assert base["support_voxels"] == int(support.sum())


def test_lesion_outside_coverage_is_unobserved_not_negative():
    target = np.zeros(SHAPE, bool)
    target[16:19, 5:8, 5:8] = True           # lesion only in slices 16-18
    support = common_support([shifted("target", 4), shifted("v2", 0)], ("target", "v2"))  # covers 0..15
    assert not support[16:].any()
    # A method that misses the uncovered lesion is not penalized; Dice is both-empty on support.
    out = forecast_dice_on_support({"miss": np.zeros(SHAPE, bool), "hit": target}, target, support)
    assert out["dice"]["miss"] == 1.0 and out["dice"]["hit"] == 1.0


def test_common_support_same_for_all_methods_and_outside_predictions_ignored():
    target = np.zeros(SHAPE, bool)
    target[5:9, 5:9, 5:9] = True
    support = common_support([shifted("target", 4), shifted("v2", 0)], ("target", "v2"))
    a = target.copy()
    b = target.copy()
    b[17:, :, :] = True                       # big false positive, all outside support
    out = forecast_dice_on_support({"a": a, "b": b}, target, support)
    assert out["dice"]["a"] == out["dice"]["b"] == 1.0
    c = target.copy()
    c[10:12, 10:12, 10:12] = True             # false positive inside support counts
    assert forecast_dice_on_support({"c": c}, target, support)["dice"]["c"] < 1.0


def test_no_eligible_band_after_restriction_returns_none():
    state, bands = history()
    support = np.zeros(SHAPE, bool)
    support[0, 0, 0] = True
    assert transition_window_on_support(state, bands, np.zeros(SHAPE, int), support) is None
