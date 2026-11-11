"""Synthetic tests for coverage-aware Track B scoring. No patient data."""
import numpy as np
import pytest

from src.track_b_coverage_scoring import (ObservableSupport, ScanCoverage, common_support, coverage_from_geometry,
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
    support = common_support(covs, ("previous", "middle", "target")).mask
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
    outside = ~support.mask
    altered[outside] = np.random.default_rng(1).choice([0, 1, 2, 3, 4], size=int(outside.sum()))
    again = transition_window_on_support(state, bands, altered, support, min_voxels=5)
    assert again["difference"] == base["difference"] and again["bands"] == base["bands"]
    assert base["support_voxels"] == int(support.mask.sum())


def test_lesion_outside_coverage_is_unobserved_not_negative():
    target = np.zeros(SHAPE, bool)
    target[16:19, 5:8, 5:8] = True           # lesion only in slices 16-18
    support = common_support([shifted("target", 4), shifted("v2", 0)], ("target", "v2"))  # covers 0..15
    assert not support.mask[16:].any()
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
    support = common_support([shifted("target", 0)], ("target",))
    assert transition_window_on_support(state, bands, np.zeros(SHAPE, int), support, min_voxels=100000) is None


# ---------------------------------------------------------------- boundary enforcement (Codex review)
def test_raw_boolean_masks_rejected_by_both_scorers():
    state, bands = history()
    raw = np.ones(SHAPE, bool)
    with pytest.raises(TypeError, match="raw masks are rejected"):
        transition_window_on_support(state, bands, np.zeros(SHAPE, int), raw)
    with pytest.raises(TypeError, match="raw masks are rejected"):
        forecast_dice_on_support({"a": raw}, raw, raw)
    with pytest.raises(TypeError):
        ObservableSupport(mask=raw, visits=("t",), provenance="scan_geometry")


def test_coverage_and_support_are_copied_and_read_only():
    cov = shifted("target", 0)
    assert not cov.mask.flags.writeable
    with pytest.raises(ValueError):
        cov.mask[0, 0, 0] = False
    support = common_support([cov], ("target",))
    assert not support.mask.flags.writeable
    with pytest.raises(ValueError):
        support.mask[0, 0, 0] = False
    with pytest.raises(Exception):
        support.mask = np.zeros(SHAPE, bool)      # frozen dataclass


def test_duplicate_visits_rejected():
    with pytest.raises(ValueError, match="duplicate coverage"):
        common_support([shifted("target", 0), shifted("target", 2)], ("target",))
    with pytest.raises(ValueError, match="distinct"):
        common_support([shifted("target", 0)], ("target", "target"))


def test_support_records_provenance_and_visits():
    support = common_support([shifted("a", 0), shifted("b", 1)], ("a", "b"))
    assert support.provenance == "scan_geometry" and support.visits == ("a", "b")
    out = forecast_dice_on_support({"m": np.zeros(SHAPE, bool)}, np.zeros(SHAPE, bool), support)
    assert out["support_visits"] == ["a", "b"] and out["dice"]["m"] == 1.0
