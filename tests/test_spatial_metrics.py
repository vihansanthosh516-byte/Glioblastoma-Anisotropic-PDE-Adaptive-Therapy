#!/usr/bin/env python3
"""Test suite for spatial validation metrics (Phase 3).

Tests src/53_spatial_metrics.py:
  - Dice Similarity Coefficient correctness
  - Hausdorff Distance correctness
  - Mean Surface Distance
  - Normalized Volume Difference
  - Boundary extraction
  - Clinical threshold gating
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
MODULE_PATH = REPO_ROOT / "src" / "53_spatial_metrics.py"

spec = importlib.util.spec_from_file_location("spatial_metrics", MODULE_PATH)
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


def _disk(shape=(50, 50), center=(25, 25), radius=10, value=1.0):
    mask = np.zeros(shape)
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    mask[(xx - center[0]) ** 2 + (yy - center[1]) ** 2 <= radius ** 2] = value
    return mask


def test_dice_identical_is_one():
    mask = _disk()
    assert m.dice_similarity_coefficient(mask, mask) == 1.0


def test_dice_disjoint_is_zero():
    a = _disk(center=(25, 25), radius=8)
    b = _disk(center=(2, 25), radius=6)  # clearly disjoint
    assert m.dice_similarity_coefficient(a, b) == 0.0


def test_dice_both_empty_is_one():
    empty = np.zeros((50, 50))
    assert m.dice_similarity_coefficient(empty, empty) == 1.0


def test_dice_partial_overlap_between_zero_and_one():
    a = _disk(center=(25, 25), radius=10)
    b = _disk(center=(27, 25), radius=10)  # shifted 2px
    dsc = m.dice_similarity_coefficient(a, b)
    assert 0.0 < dsc < 1.0, f"DSC={dsc} should be strictly between 0 and 1"


def test_dice_thresholds_continuous_mask():
    a = _disk(value=0.8)
    b = _disk(center=(25, 25), radius=10, value=0.5)
    assert m.dice_similarity_coefficient(a, b, threshold=0.1) == 1.0


def test_hausdorff_identical_is_zero():
    mask = _disk()
    assert m.hausdorff_distance(mask, mask) == 0.0


def test_hausdorff_disjoint_positive():
    a = _disk(center=(25, 25), radius=8)
    b = _disk(center=(2, 25), radius=6)
    hd = m.hausdorff_distance(a, b)
    assert hd > 0.0, "Disjoint masks must have HD>0"


def test_hausdorff_empty_returns_zero():
    a = _disk()
    empty = np.zeros((50, 50))
    assert m.hausdorff_distance(a, empty) == 0.0


def test_hausdorff95_le_max_hd():
    a = _disk(center=(25, 25), radius=10)
    b = _disk(center=(27, 25), radius=10)  # shifted 2px
    hd95 = m.hausdorff_distance(a, b, percentile=95)
    hd_max = m.hausdorff_distance(a, b)
    assert hd95 <= hd_max, "HD95 must be <= max HD"


def test_mean_surface_distance_identical_is_zero():
    mask = _disk()
    assert m.mean_surface_distance(mask, mask) == 0.0


def test_normalized_volume_difference_identical():
    assert m.normalized_volume_difference(100.0, 100.0) == 0.0


def test_normalized_volume_difference_half():
    assert m.normalized_volume_difference(50.0, 100.0) == 0.5


def test_normalized_volume_difference_both_zero():
    assert m.normalized_volume_difference(0.0, 0.0) == 0.0


def test_extract_boundary_nonempty():
    mask = _disk(radius=10)
    bnd = m.extract_boundary(mask)
    assert bnd.shape[0] > 0, "Boundary should be non-empty for a non-empty mask"
    assert bnd.shape[1] == 2, "2D mask boundary should have 2 columns"


def test_extract_boundary_empty_mask():
    empty = np.zeros((20, 20))
    bnd = m.extract_boundary(empty)
    assert bnd.shape[0] == 0, "Empty mask boundary should be empty"


def test_compute_spatial_metrics_returns_full_schema():
    a = _disk(center=(25, 25), radius=10)
    b = _disk(center=(27, 25), radius=10)
    out = m.compute_spatial_metrics(a, b)
    required = [
        "dice_coefficient", "hausdorff_distance_mm",
        "mean_surface_distance_mm", "normalized_volume_difference",
        "meets_clinical_threshold",
    ]
    for k in required:
        assert k in out, f"Missing key '{k}' in compute_spatial_metrics output"


def test_meets_clinical_threshold_pass():
    a = _disk(center=(25, 25), radius=10)
    out = m.compute_spatial_metrics(a, a)  # identical -> DSC=1, HD=0
    assert out["meets_clinical_threshold"] is True


def run_all():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed, failed = 0, 0
    for test in tests:
        try:
            test()
            print(f"  PASS  {test.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  FAIL  {test.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ERROR {test.__name__}: {type(e).__name__}: {e}")
            failed += 1
    print(f"\nResults: {passed} passed, {failed} failed, {passed + failed} total")
    return failed == 0


if __name__ == "__main__":
    sys.exit(0 if run_all() else 1)
