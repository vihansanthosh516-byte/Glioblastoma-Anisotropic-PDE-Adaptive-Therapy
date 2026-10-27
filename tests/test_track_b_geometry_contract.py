import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest


MODULE_PATH = Path(__file__).resolve().parents[1] / "src" / "160_growth_zone_features.py"
SPEC = importlib.util.spec_from_file_location("growth_zone_features_geometry_test", MODULE_PATH)
FEATURES = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = FEATURES
SPEC.loader.exec_module(FEATURES)
SIGNED_PATH = Path(__file__).resolve().parents[1] / "src" / "163_signed_transition.py"
SIGNED_SPEC = importlib.util.spec_from_file_location("signed_transition_geometry_test", SIGNED_PATH)
SIGNED = importlib.util.module_from_spec(SIGNED_SPEC)
sys.modules[SIGNED_SPEC.name] = SIGNED
SIGNED_SPEC.loader.exec_module(SIGNED)


def _visit(affine):
    shape = (3, 4, 5)
    images = {name: np.arange(np.prod(shape), dtype=np.float32).reshape(shape) for name in FEATURES.MODALITIES}
    labels = np.zeros(shape, dtype=np.uint8)
    labels[1, 1, 1] = 1
    return FEATURES.VisitVolume(images=images, labels=labels, affine=np.asarray(affine, dtype=np.float64))


def test_validate_visit_accepts_orthogonal_oblique_affine_and_returns_spacing():
    angle = np.deg2rad(31.0)
    rotation = np.array(
        [[np.cos(angle), -np.sin(angle), 0.0],
         [np.sin(angle), np.cos(angle), 0.0],
         [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    affine = np.eye(4, dtype=np.float64)
    affine[:3, :3] = rotation @ np.diag([1.0, 1.5, 2.0])

    *_, spacing = FEATURES.validate_visit(_visit(affine))

    assert spacing == pytest.approx([1.0, 1.5, 2.0])


def test_validate_visit_rejects_sheared_affine_before_distance_calculation():
    affine = np.eye(4, dtype=np.float64)
    affine[0, 1] = 0.2

    with pytest.raises(ValueError, match="non-orthogonal voxel axes"):
        FEATURES.validate_visit(_visit(affine))


def test_signed_transition_targets_use_disjoint_forecast_origin_risk_sets():
    current = np.zeros((2, 2, 2), dtype=np.uint8)
    current[0, 0, 0] = 1  # Current tumor: eligible for removal, not addition.
    current[0, 0, 1] = 4  # Current cavity: excluded from addition risk.
    future = np.zeros_like(current)
    future[0, 0, 1] = 3  # Future label cannot make current cavity an at-risk voxel.
    future[0, 1, 0] = 1  # New annotation in addition risk.
    domain = np.ones_like(current, dtype=bool)

    targets = SIGNED.make_signed_transition_targets(current, future, domain)

    assert targets.addition_risk.sum() == 6
    assert targets.removal_risk.sum() == 1
    assert targets.additions.sum() == 1
    assert targets.removals.sum() == 1
    assert not np.any(targets.addition_risk & targets.removal_risk)
    assert not targets.addition_risk[0, 0, 0]
    assert not targets.addition_risk[0, 0, 1]
