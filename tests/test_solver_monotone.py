"""Tests for the positivity-preserving solver (src/solver_monotone.py; GRAND_PLAN #15, weak point W21)."""
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "src"))

import run_improved_aniso as ria  # noqa: E402
from solver_monotone import TensorFKMonotone, lattice_directions  # noqa: E402


def _mass(a):
    return float(np.asarray(a).sum(dtype=np.float64))


def test_direction_sets():
    assert len(lattice_directions(False)) == 13
    assert len(lattice_directions(True)) == 37


def test_isotropic_equals_tensorfk():
    shape = (16, 16, 16)
    iso = np.zeros(shape + (6,), np.float32)
    iso[..., :3] = 0.3
    mask = np.ones(shape, bool)
    u0 = np.zeros(shape, np.float32)
    u0[6:10, 6:10, 6:10] = 1
    a, _ = ria.TensorFK(iso, mask, h=2.0).run(u0, [0.1], 30.0)
    b, _ = TensorFKMonotone(iso, mask, h=2.0).run(u0, [0.1], 30.0)
    assert np.abs(a - b).max() < 1e-6


def test_strong_anisotropy_adds_no_mass():
    # script 97 V5 stress case: TensorFK's clamp adds 2.8% mass here
    T = ria.mat_to_six(np.array([[0.6, 0.5, 0.0], [0.5, 0.6, 0.0], [0.0, 0.0, 0.05]])).astype(np.float32)
    d6 = np.tile(T, (20, 20, 20, 1)) * 0.3
    u0 = np.zeros((20, 20, 20), np.float32)
    u0[8:12, 8:12, 8:12] = 1.0
    s = TensorFKMonotone(d6, np.ones((20, 20, 20), bool), h=2.0)
    u, _ = s.run(u0, [0.0], 100.0)
    assert float(s.clamp_added[0]) == 0.0
    assert abs(_mass(u[0]) - _mass(u0)) / _mass(u0) < 1e-6
    assert s.residual_summary["rel_residual_max"] < 1e-6


def test_random_spd_conserves_mass_and_stays_bounded():
    shape = (20, 18, 16)
    zz, yy, xx = np.meshgrid(*[np.arange(n) for n in shape], indexing="ij")
    mask = ((xx - 8) ** 2 + (yy - 9) ** 2 + (zz - 10) ** 2) < 7 ** 2
    rng = np.random.default_rng(1)
    A = rng.normal(size=shape + (3, 3)) * 0.3
    d6 = ria.mat_to_six(A @ np.swapaxes(A, -1, -2) + 0.1 * np.eye(3)) * 0.3
    u0 = (rng.random(shape) * mask).astype(np.float32)
    s = TensorFKMonotone(d6, mask, h=2.0)
    u, _ = s.run(u0, [0.0], 60.0)
    assert float(s.clamp_added[0]) == 0.0
    assert abs(_mass(u[0]) - _mass(u0)) / _mass(u0) < 1e-6
    assert float(u.min()) >= 0.0 and float(u.max()) <= 1.0
    assert _mass(u[0][~mask]) == 0.0


def test_growth_stays_in_unit_interval_without_clamp():
    T = ria.mat_to_six(np.array([[0.6, 0.5, 0.0], [0.5, 0.6, 0.0], [0.0, 0.0, 0.05]])).astype(np.float32)
    d6 = np.tile(T, (16, 16, 16, 1)) * 0.3
    u0 = np.zeros((16, 16, 16), np.float32)
    u0[6:10, 6:10, 6:10] = 1.0
    s = TensorFKMonotone(d6, np.ones((16, 16, 16), bool), h=2.0)
    u, _ = s.run(u0, [0.2], 60.0)
    assert float(s.clamp_added[0]) == 0.0 and float(s.clamp_removed[0]) == 0.0
    assert float(u.min()) >= 0.0 and float(u.max()) <= 1.0
