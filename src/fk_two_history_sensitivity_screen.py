"""Synthetic-only sensitivity screen of the two-history shape control (Track B development).

Synthetic 3D ellipsoid-like tumours inside a fixed ellipsoidal support S, three truth regimes (growth with
anisotropic drift, shrinkage, plateau), varied history/forecast gaps. The target mask is used ONLY to score
and to supply the TRUE volume draw; every arm receives the same volume draws. No patient data, no PDE.

Arms: copy_last (h2 unchanged, NOT volume matched), last_mask_sd (control with f=0, volume matched; this is also
the volume-matched copy-last, so no separate arm), control_default (module defaults) and a grid over
SHRINK_TAU_DAYS x MAX_FACTOR x MAX_SD_DELTA_MM. Pre-boundary (apply_boundary=False quantities) and
boundary-perturbed (shared candidate seeds, fixed synthetic covariance) scores are separate fields, never merged.
Results are screening evidence only; they do not validate the control as a patient baseline.
"""
from __future__ import annotations

import argparse
import itertools
import json
import time
from pathlib import Path as _Path

import numpy as np

from src import fk_two_history_shape as th

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
DEFAULT_OUT = PROJECT_ROOT / "output" / "track_b_two_history_screen_20261102" / "screen.json"

VOXEL_MM = 2.0
BOUNDARY_COV = {"sigma_b_mm": 2.0, "tau_b_mm": 0.5, "ell_mm": 8.0}  # fixed synthetic, not calibrated
REGIMES = ("growth", "shrinkage", "plateau")
HISTORY_WEEKS = (4, 8, 16)
FORECAST_WEEKS = (4, 8, 16)
TAUS = (0.0, 7.0, 28.0)
CAPS = (1.0, 4.0, 8.0)
CLIPS = (10.0, 30.0)
NOISE_SD_LOG = 0.15
N_NOISY = 3
T0_WEEKS = 10.0  # process time of h1
DICE_FIELDS = ("dice_true_vol_pre", "dice_noisy_vol_pre", "dice_true_vol_boundary", "dice_noisy_vol_boundary")


def dice(a, b) -> float:
    a, b = np.asarray(a, bool), np.asarray(b, bool)
    s = a.sum() + b.sum()
    return 1.0 if s == 0 else float(2.0 * (a & b).sum() / s)


def support(n: int) -> np.ndarray:
    g = (np.indices((n,) * 3) - (n - 1) / 2.0) / (n / 2.0)
    return (g[0] / 0.80) ** 2 + (g[1] / 0.68) ** 2 + (g[2] / 0.72) ** 2 <= 1.0


def _smooth_field(rng, n):
    f = ndimage_gauss(rng.standard_normal((n,) * 3), n / 8.0)
    return f / (f.std() + 1e-12)


def ndimage_gauss(a, s):
    from scipy import ndimage
    return ndimage.gaussian_filter(a, s, mode="wrap")


def truth_mask(regime: str, t_weeks: float, case_seed: int, n: int, S: np.ndarray) -> np.ndarray:
    """Known synthetic process: ellipsoid whose radii/centre follow the regime, plus a persistent irregular field."""
    rng = np.random.default_rng(1000 + case_seed)
    scale = n / 48.0
    c0 = (n - 1) / 2.0 + rng.uniform(-3, 3, 3) * scale
    drift = rng.uniform(0.6, 1.0) * np.array([0.10, 0.04, 0.0]) * scale * rng.choice([-1, 1], 3)
    axis_rate = rng.uniform(0.8, 1.2) * np.array([0.22, 0.08, 0.06]) * scale  # anisotropic voxel/week
    irr = _smooth_field(rng, n)
    t = float(t_weeks)
    if regime == "growth":
        r0, r, ctr = np.array([4.0, 4.0, 4.0]) * scale, None, None
        r = r0 + axis_rate * t
        ctr = c0 + drift * t
    elif regime == "shrinkage":
        r0 = np.array([9.0, 8.0, 8.0]) * scale
        r = np.maximum(r0 - 0.55 * axis_rate * t * 2.0, 1.5 * scale)
        ctr = c0 + 0.3 * drift * t
    elif regime == "plateau":
        r0 = np.array([4.0, 4.0, 4.0]) * scale
        r = r0 + axis_rate * 3.0 * 3.0 * (1.0 - np.exp(-t / 3.0))
        ctr = c0 + drift * 3.0 * (1.0 - np.exp(-t / 3.0))
    else:
        raise ValueError(regime)
    g = np.indices((n,) * 3).astype(float)
    q = sum(((g[i] - ctr[i]) / r[i]) ** 2 for i in range(3))
    m = (q + 0.12 * irr) <= 1.0
    return m & S


def volume_draws(true_vol: int, cap: int, case_seed: int) -> np.ndarray:
    rng = np.random.default_rng(5000 + case_seed)
    noisy = np.rint(true_vol * np.exp(rng.normal(0, NOISE_SD_LOG, N_NOISY))).astype(np.int64)
    d = np.concatenate([[true_vol], noisy]).astype(np.int64)
    return np.clip(d, 1, cap - 1)


def control_settings(grid=True):
    out = {"control_default": (th.SHRINK_TAU_DAYS, th.MAX_FACTOR, th.MAX_SD_DELTA_MM)}
    if grid:
        for tau, cap, clip in itertools.product(TAUS, CAPS, CLIPS):
            out[f"grid_tau{tau:g}_cap{cap:g}_clip{clip:g}"] = (tau, cap, clip)
    return out


def _score(draws, target):
    return [dice(d, target) for d in draws]


def _fields(pre, boundary, target):
    p = _score(pre, target)
    f = {"dice_true_vol_pre": p[0], "dice_noisy_vol_pre": float(np.mean(p[1:]))}
    if boundary is None:
        f.update(dice_true_vol_boundary=None, dice_noisy_vol_boundary=None)
    else:
        b = _score(boundary, target)
        f.update(dice_true_vol_boundary=b[0], dice_noisy_vol_boundary=float(np.mean(b[1:])))
    return f


def run_case(regime, hist_w, fc_w, seed, n=48, grid=True, n_anchors=64) -> dict:
    S = support(n)
    case_seed = seed + 17 * hist_w + 131 * fc_w + 7919 * REGIMES.index(regime)
    h1 = truth_mask(regime, T0_WEEKS, case_seed, n, S)
    h2 = truth_mask(regime, T0_WEEKS + hist_w, case_seed, n, S)
    target = truth_mask(regime, T0_WEEKS + hist_w + fc_w, case_seed, n, S)
    cap = int(S.sum())
    if not (h1.any() and h2.any() and target.any()):
        return {"regime": regime, "history_weeks": hist_w, "forecast_weeks": fc_w, "seed": seed,
                "skipped": "empty mask"}
    vols = volume_draws(int(target.sum()), cap, case_seed)
    gap, fc = 7.0 * hist_w, 7.0 * fc_w
    bseeds = th.candidate_boundary_seeds(case_seed, len(vols))
    spacing = VOXEL_MM
    arms = {"copy_last": _fields([h2] * len(vols), None, target)}  # not volume matched; no boundary variant

    def run_arm(g_days, h1_):
        r = th.two_history_shape_draws(h1_, h2, g_days, fc, vols, domain=S, spacing_mm=spacing,
                                       boundary_covariance=BOUNDARY_COV, boundary_seeds=bseeds,
                                       n_anchors=n_anchors, apply_boundary=True)
        return _fields(r["pre_boundary_mask_draws"], r["mask_draws"], target) | {"factor": r["factor"]}

    arms["last_mask_sd"] = run_arm(0.0, h1)  # gap 0 -> f=0 (last-mask rule)
    saved = (th.SHRINK_TAU_DAYS, th.MAX_FACTOR, th.MAX_SD_DELTA_MM)
    try:
        for name, (tau, capf, clip) in control_settings(grid).items():
            th.SHRINK_TAU_DAYS, th.MAX_FACTOR, th.MAX_SD_DELTA_MM = tau, capf, clip
            arms[name] = run_arm(gap, h1)
    finally:
        th.SHRINK_TAU_DAYS, th.MAX_FACTOR, th.MAX_SD_DELTA_MM = saved
    return {"regime": regime, "history_weeks": hist_w, "forecast_weeks": fc_w, "seed": seed,
            "volumes_voxels": [int(v) for v in vols], "h2_voxels": int(h2.sum()), "arms": arms}


def _agg(cases, key_fn):
    groups: dict = {}
    for c in cases:
        if "arms" in c:
            groups.setdefault(key_fn(c), []).append(c)
    out = {}
    for key, cs in sorted(groups.items()):
        entry = {"n_cases": len(cs), "arms": {}}
        for arm in cs[0]["arms"]:
            a = {}
            for fld in DICE_FIELDS:
                vals = [c["arms"][arm][fld] for c in cs]
                if any(v is None for v in vals):
                    a[fld] = None
                    continue
                base = [c["arms"]["last_mask_sd"][fld] for c in cs]
                d = np.array(vals) - np.array(base)
                se = float(d.std(ddof=1) / np.sqrt(len(d))) if len(d) > 1 else None
                a[fld] = {"mean": float(np.mean(vals)), "paired_diff_vs_last_mask_sd_mean": float(d.mean()),
                          "paired_diff_se": se, "frac_cases_better": float((d > 0).mean())}
            entry["arms"][arm] = a
        out[key] = entry
    return out


def run_screen(regimes=REGIMES, hist=HISTORY_WEEKS, fcs=FORECAST_WEEKS, seeds=(0, 1), n=48, grid=True) -> dict:
    t0 = time.time()
    cases = [run_case(r, h, f, s, n=n, grid=grid)
             for r in regimes for h in hist for f in fcs for s in seeds]
    return {
        "description": "synthetic-only sensitivity screen; pre_boundary and boundary_perturbed are separate fields",
        "config": {"grid_n": n, "voxel_mm": VOXEL_MM, "boundary_covariance": BOUNDARY_COV, "regimes": list(regimes),
                   "history_weeks": list(hist), "forecast_weeks": list(fcs), "seeds": list(seeds),
                   "noisy_volume_sd_log": NOISE_SD_LOG, "n_noisy_draws": N_NOISY, "t0_weeks": T0_WEEKS,
                   "control_defaults": {"SHRINK_TAU_DAYS": th.SHRINK_TAU_DAYS, "MAX_FACTOR": th.MAX_FACTOR,
                                        "MAX_SD_DELTA_MM": th.MAX_SD_DELTA_MM},
                   "note": "volume-matched copy-last is identical to last_mask_sd; copy_last has no boundary field"},
        "runtime_seconds": time.time() - t0,
        "by_regime": _agg(cases, lambda c: c["regime"]),
        "by_regime_history_forecast": _agg(
            cases, lambda c: f"{c['regime']}|hist{c['history_weeks']}w|fc{c['forecast_weeks']}w"),
        "cases": cases,
    }


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--grid-n", type=int, default=48)
    a = ap.parse_args(argv)
    res = run_screen(seeds=tuple(range(a.seeds)), n=a.grid_n)
    out = _Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=1), encoding="utf-8")
    print(f"wrote {out} ({res['runtime_seconds']:.0f}s)")


if __name__ == "__main__":
    main()
