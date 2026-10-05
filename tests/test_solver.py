"""Solver verification tests (analysis_plan_v1.md A1.19). Wraps the fast checks in src/97_solver_verification.py.

V7 and V9 (production-grid resolution) are slow and informational; run `python src/97_solver_verification.py`.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

spec = importlib.util.spec_from_file_location("solver_verification", REPO / "src" / "97_solver_verification.py")
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)


def test_spatial_order_is_second():
    r = sv.v1_spatial_order()
    assert r["pass"], r["orders_used"]


def test_temporal_order_is_first():
    r = sv.v2_temporal_order()
    assert r["pass"], r["orders_used"]


def test_front_converges_with_grid():
    assert sv.v3_grid_convergence()["pass"]


def test_isotropic_run_conserves_mass():
    r = sv.v4_conservation()
    assert r["rel_mass_change"] < 1e-5


def test_solution_bounded_and_finite():
    assert sv.v5_boundedness()["pass"]


def test_no_flux_outside_domain():
    r = sv.v6_boundary()
    assert r["mass_outside_mask"] == 0.0 and r["interior_rel_mass_change"] < 1e-5


def test_zero_growth_zero_diffusion_is_identity():
    import numpy as np
    shape = (8, 8, 8)
    s = sv.ria.TensorFK(np.zeros(shape + (6,), np.float32), np.ones(shape, bool), h=2.0)
    u0 = np.random.default_rng(0).random(shape).astype(np.float32)
    u, _ = s.run(u0, [0.0], 50.0)
    assert np.allclose(u[0], u0, atol=1e-7)


def test_pure_growth_matches_logistic():
    import math
    import numpy as np
    shape = (6, 6, 6)
    s = sv.ria.TensorFK(np.zeros(shape + (6,), np.float32), np.ones(shape, bool), h=2.0)
    u0 = np.full(shape, 0.1, np.float32)
    rho, t = 0.05, 40.0
    u, _ = s.run(u0, [rho], t)
    exact = 1 / (1 + (1 / 0.1 - 1) * math.exp(-rho * t))
    assert abs(float(u[0].mean()) / exact - 1) < 0.02
