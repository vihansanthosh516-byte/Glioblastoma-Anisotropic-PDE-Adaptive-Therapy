"""Track B synthetic parameter recovery (src/fk_mask_inverse.py, src/51_inverse_parameter_estimation.py).

Reference data come from the same TensorFK equation at h = 1 mm, block-averaged to 2 mm; the estimator fits at
h = 2 mm. So the estimator never sees its own numerical path, and the 1 mm vs 2 mm discretisation mismatch is
real and measured. Synthetic only: this does not show recovery from patient scans (unknown model error,
segmentation noise, treatment, anisotropy, irregular shapes are all absent).

Geometry: one octant of an ellipsoid (semi-axes 11.3, 8.7, 7.1 mm) at the grid corner; the no-flux faces act
as mirror planes. An ellipsoid, not a sphere, because a sphere flips whole voxel shells at once and makes the
thresholded volume a coarse staircase in the parameters. Scans at days 0, 60, 120, 180.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import fk_mask_inverse as fi  # noqa: E402

_spec = importlib.util.spec_from_file_location("s51_t", REPO / "src" / "51_inverse_parameter_estimation.py")
s51 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s51)

L_MM = 40.0
AXES = (11.3, 8.7, 7.1)
DAYS = [0, 60, 120, 180]
TRUTHS = [(0.03, 0.1), (0.05, 0.03)]   # (rho 1/day, D mm^2/day); front speeds 0.11 and 0.077 mm/day


def ellipsoid(h):
    n = int(round(L_MM / h))
    c = (np.arange(n) + 0.5) * h
    q = (c[:, None, None] / AXES[0]) ** 2 + (c[None, :, None] / AXES[1]) ** 2 + (c[None, None, :] / AXES[2]) ** 2
    return q <= 1.0


def down2(m):
    n = m.shape[0] // 2
    return m.reshape(n, 2, n, 2, n, 2).mean((1, 3, 5)) >= 0.5


@pytest.fixture(scope="module", params=TRUTHS, ids=lambda t: f"rho{t[0]}_D{t[1]}")
def reference(request):
    rho, D = request.param
    m0 = ellipsoid(1.0)
    fine = [m0] + fi.forward_masks(m0, np.ones(m0.shape, bool), rho, D, DAYS, 1.0)
    return {"truth": (rho, D), "fine": fine, "coarse": [down2(m) for m in fine]}


# ---- physical units: dimensionless groups of the equation ----

def _u(rho, D, t, h, n=12):
    m = np.zeros((n, n, n), bool)
    m[:4, :4, :4] = True
    s = fi.ria.TensorFK(np.tile(np.array([D, D, D, 0, 0, 0], np.float32), (n, n, n, 1)), np.ones((n,) * 3, bool), h=h)
    return s.run(m.astype(np.float32), [rho], t)[0][0]


def test_time_rescaling_leaves_solution_unchanged():
    """u depends on rho*t and D*t/h^2 only: (rho, D, t) -> (2 rho, 2 D, t/2) gives the same field."""
    assert np.allclose(_u(0.1, 0.1, 60.0, 2.0), _u(0.2, 0.2, 30.0, 2.0), atol=1e-5)


def test_space_rescaling_leaves_solution_unchanged():
    """(h, D) -> (2h, 4D) gives the same voxel field: D is in mm^2/day and h in mm."""
    assert np.allclose(_u(0.05, 0.1, 60.0, 1.0), _u(0.05, 0.4, 60.0, 2.0), atol=1e-5)


# ---- identifiability of the observation design ----

def test_equal_front_speed_pairs_give_different_scans():
    """Asymptotic Fisher speed 2 sqrt(rho D) is the same for both pairs, yet over 180 days from a sharp
    initial mask the volumes differ: the transient carries separate information on rho and D."""
    m0 = down2(ellipsoid(1.0))
    dom = np.ones(m0.shape, bool)
    a = fi.forward_masks(m0, dom, 0.03, 0.1, DAYS, 2.0)
    b = fi.forward_masks(m0, dom, 0.075, 0.04, DAYS, 2.0)
    assert all(abs(np.log(x.sum() / y.sum())) > 0.2 for x, y in zip(a[1:], b[1:]))


def test_one_followup_leaves_a_ridge():
    """With a single follow-up, a second (rho, D) pair matches that scan's volume to within 2% yet
    disagrees later: one follow-up constrains only a combination of rho and D."""
    m0 = down2(ellipsoid(1.0))
    dom = np.ones(m0.shape, bool)
    target = fi.forward_masks(m0, dom, 0.03, 0.1, DAYS, 2.0)
    lo, hi = 1e-3, 0.1   # bisection on D at rho = 0.06 for the day-60 volume
    for _ in range(14):
        mid = np.sqrt(lo * hi)
        v = fi.forward_masks(m0, dom, 0.06, mid, DAYS[:2], 2.0)[0].sum()
        lo, hi = (mid, hi) if v < target[0].sum() else (lo, mid)
    alt = fi.forward_masks(m0, dom, 0.06, np.sqrt(lo * hi), DAYS, 2.0)
    assert abs(np.log(alt[0].sum() / target[0].sum())) < 0.02
    assert abs(np.log(alt[2].sum() / target[2].sum())) > 0.1


def test_three_followups_give_full_rank(reference):
    m = reference["coarse"]
    J = fi.jacobian(np.log(reference["truth"]), m[0], np.ones(m[0].shape, bool), m[1:], DAYS, 2.0)
    sv = np.linalg.svd(J, compute_uv=False)
    assert sv[1] > 0 and sv[0] / sv[1] < 100


# ---- recovery against the finer reference ----

def test_recovery_within_discretisation_bias(reference):
    """Tolerance is derived, not chosen: r0 = residual of the 2 mm model AT the true parameters against the
    1 mm reference (pure discretisation mismatch). Linearised, it shifts the estimate by J^+ r0. The test
    allows 2 |J^+ r0| + 0.05 in log parameters (factor 2 for nonlinearity and the volume staircase, 0.05 for
    the optimiser's 1e-3 tolerance and finite-difference J), and requires the fit to be at least as good as
    the truth, to stay off the bounds and to match each follow-up mask with Dice >= 0.95."""
    m = reference["coarse"]
    dom = np.ones(m[0].shape, bool)
    xt = np.log(reference["truth"])
    r0 = fi.log_volume_residuals(xt, m[0], dom, m[1:], DAYS, 2.0)
    J = fi.jacobian(xt, m[0], dom, m[1:], DAYS, 2.0)
    bias = np.abs(np.linalg.pinv(J) @ r0)
    fit = fi.fit_fk_masks(m, DAYS, 2.0)
    err = np.abs(np.log([fit["rho"], fit["D"]]) - xt)
    print(reference["truth"], {k: fit[k] for k in ("rho", "D", "sse_log_volume", "dice", "jacobian_condition")},
          "r0", r0.round(4), "bias bound", bias.round(3), "log err", err.round(3))
    assert not fit["rho_at_bound"] and not fit["D_at_bound"]
    assert fit["sse_log_volume"] <= float(r0 @ r0) + 1e-12
    assert min(fit["dice"]) >= 0.95
    assert (err <= 2 * bias + 0.05).all()


# ---- the existing volume surrogate (script 51) on the same data ----

def test_volume_surrogate_cannot_recover_pde_parameters(reference):
    """Known limitation, kept as a test so a change is noticed. Script 51 fits dV/dt = rho V + (36 pi)^(1/3) D V^(1/3)
    (radial speed D/R), not the Fisher front. On PDE volumes it puts D on its lower bound and misses rho."""
    vols = [8.0 * float(x.sum()) for x in reference["fine"]]   # octant x 8, 1 mm^3 voxels
    fit = s51.fit_trajectory(DAYS, vols, [(s51.RHO_MIN, s51.RHO_MAX), (s51.D_MIN, s51.D_MAX)])
    rho, D = reference["truth"]
    assert fit["D_at_lower_bound"]
    assert abs(np.log(fit["rho"] / rho)) > np.log(2)
