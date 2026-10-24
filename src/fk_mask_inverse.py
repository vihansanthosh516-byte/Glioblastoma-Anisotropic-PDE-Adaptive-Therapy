"""Fit isotropic Fisher-Kolmogorov (rho, D) to a series of binary tumour masks (Track B; no data I/O).

Forward model: run_improved_aniso.TensorFK with an isotropic tensor D * I on a masked grid of spacing h (mm),
seeded with the first mask as u = 1, thresholded at run_improved_aniso.U_VISIBLE at each later scan day.
Objective: sum over follow-up scans of log(simulated voxels / observed voxels)^2, minimised over
(log rho, log D) by bounded Nelder-Mead from two starts (the thresholded volume is piecewise constant in the
parameters, so gradient methods are not used). Units: rho 1/day, D mm^2/day, days, h mm.

Separate recovery of rho and D needs at least two follow-up scans; with one there is a single residual for
two parameters and only a combination is constrained. fit_fk_masks reports n_followups and the condition
number of the finite-difference Jacobian so callers can see this; it does not decide identifiability for them.
"""
from __future__ import annotations

import sys
from pathlib import Path as _Path

import numpy as np
from scipy.optimize import minimize

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import run_improved_aniso as ria  # noqa: E402

BOUNDS = ((1e-3, 0.5), (1e-3, 1.0))   # (rho 1/day, D mm^2/day): wide, so "interior" is a real check
STARTS = ((0.01, 0.03), (0.1, 0.3))


def forward_masks(mask0: np.ndarray, domain: np.ndarray, rho: float, D: float, days, h: float) -> list:
    """Binary masks at days[1:], integrating continuously from days[0]."""
    d6 = np.zeros(np.asarray(mask0).shape + (6,), np.float32)
    d6[..., :3] = float(D)
    return forward_masks_tensor(mask0, domain, rho, d6, days, h)


def forward_masks_tensor(mask0: np.ndarray, domain: np.ndarray, rho: float, diffusion_tensor: np.ndarray,
                         days, h: float) -> list:
    """Binary masks at days[1:] under a supplied six-component diffusion tensor field.

    This is used to generate synthetic anisotropic references. The fitter remains isotropic, so comparing
    those references with ``fit_fk_masks`` is a deliberate model-mismatch stress test, not patient evidence.
    Tensor component order is xx, yy, zz, xy, xz, yz, as required by ``run_improved_aniso.TensorFK``.
    """
    mask0 = np.asarray(mask0, dtype=bool)
    domain = np.asarray(domain, dtype=bool)
    diffusion_tensor = np.asarray(diffusion_tensor, dtype=np.float32)
    days = np.asarray(days, dtype=float)
    if mask0.ndim != 3 or domain.shape != mask0.shape:
        raise ValueError("mask0 and domain must have the same 3D shape")
    if diffusion_tensor.shape != mask0.shape + (6,) or not np.isfinite(diffusion_tensor).all():
        raise ValueError("diffusion_tensor must be finite with shape mask0.shape + (6,)")
    if days.ndim != 1 or len(days) < 2 or not np.isfinite(days).all() or np.any(np.diff(days) <= 0):
        raise ValueError("days must contain at least two finite, strictly increasing values")
    if not np.isfinite(rho) or rho < 0 or not np.isfinite(h) or h <= 0:
        raise ValueError("rho must be nonnegative and h must be positive")

    s = ria.TensorFK(diffusion_tensor, domain, h=h)
    u = mask0.astype(np.float32)
    out, prev = [], days[0]
    for t in days[1:]:
        u = s.run(u, [rho], float(t - prev))[0][0]
        prev = t
        out.append(u > ria.U_VISIBLE)
    return out


def log_volume_residuals(x, mask0, domain, observed, days, h) -> np.ndarray:
    sim = forward_masks(mask0, domain, float(np.exp(x[0])), float(np.exp(x[1])), days, h)
    return np.array([np.log(max(int(a.sum()), 1) / max(int(b.sum()), 1)) for a, b in zip(sim, observed)])


def jacobian(x, mask0, domain, observed, days, h, step: float = 0.05) -> np.ndarray:
    """Central differences in log parameters; step 0.05 (5%) spans several voxel flips, so the staircase in
    the thresholded volume does not dominate."""
    cols = [(log_volume_residuals(x + step * e, mask0, domain, observed, days, h)
             - log_volume_residuals(x - step * e, mask0, domain, observed, days, h)) / (2 * step) for e in np.eye(2)]
    return np.stack(cols, 1)


def fit_fk_masks(masks: list, days, h: float, domain: np.ndarray | None = None, bounds=BOUNDS,
                 starts=STARTS, maxfev: int = 120) -> dict:
    """masks[0] seeds the model; masks[1:] are the observations at days[1:]."""
    if len(masks) != len(days) or len(masks) < 2:
        raise ValueError("need one mask per day and at least one follow-up")
    if any(b <= a for a, b in zip(days, days[1:])):
        raise ValueError("scan days must increase")
    domain = np.ones(masks[0].shape, bool) if domain is None else domain
    mask0, observed = masks[0], masks[1:]
    lb = [(np.log(lo), np.log(hi)) for lo, hi in bounds]

    def f(x):
        return float(np.sum(log_volume_residuals(x, mask0, domain, observed, days, h) ** 2))

    best = None
    for s in starts:
        r = minimize(f, np.log(s), method="Nelder-Mead", bounds=lb,
                     options={"xatol": 1e-3, "fatol": 1e-8, "maxfev": maxfev})
        if best is None or r.fun < best.fun:
            best = r
    x = best.x
    J = jacobian(x, mask0, domain, observed, days, h)
    sv = np.linalg.svd(J, compute_uv=False)
    rho, D = float(np.exp(x[0])), float(np.exp(x[1]))
    sim = forward_masks(mask0, domain, rho, D, days, h)
    return {"rho": rho, "D": D, "front_speed_mm_per_day": 2.0 * float(np.sqrt(rho * D)),
            "sse_log_volume": float(best.fun), "n_followups": len(observed), "n_evals": int(best.nfev),
            "dice": [float(ria.dice(a, b)) for a, b in zip(sim, observed)],
            "jacobian_condition": float(sv[0] / sv[-1]) if len(sv) == 2 and sv[-1] > 0 else float("inf"),
            "rho_at_bound": bool(min(abs(x[0] - lb[0][0]), abs(x[0] - lb[0][1])) < 1e-3),
            "D_at_bound": bool(min(abs(x[1] - lb[1][0]), abs(x[1] - lb[1][1])) < 1e-3)}
