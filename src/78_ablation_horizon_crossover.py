#!/usr/bin/env python3
"""
Script 78: Horizon-crossover ablation. At what horizon / diffusivity does DTI
orientation start to matter? (Track B negative 2, follow-up to script 77.)
============================================================================
Script 77 showed that at 2 mm voxels and 90 d the diffusion length (0.46-1.18 mm)
is below one voxel, so removing DTI cannot change anything (mask Dice 0.994).
Script 74 showed that anisotropy changes shape at 600 d. This script maps the
crossover instead of asserting it.

PRE-SPECIFIED DESIGN (fixed before the first run; not changed after seeing data)
--------------------------------------------------------------------------------
Model     script-60 FastPDESolver tensor (synthetic tracts), natural growth (no kill),
          seed = script-60 seed x 0.1 (peak 0.08 K, sigma 5 mm), dt 0.2 d, K = 1.
          The tensor is built on the 64^3 (2 mm) grid and repeated 2x2x2 onto the
          128^3 (1 mm) grid, so both grids see the SAME physical tensor field.
Arms      full         script-60 DTI tensor
          iso_matched  D_iso(x) = trace(D(x))/3 * I   (same mean diffusivity, no orientation)
                       -> PRIMARY comparison: full vs iso_matched
          x_aligned    script-60 "No DTI" arm (uniform x-aligned tensor), for continuity
                       with script 60/77 -> secondary
Sweep     horizon  {90, 180, 365, 600, 900} d (snapshots of one run)
          D_white  {0.01, 0.03, 0.1, 0.3} mm^2/d; D_gray = D_white / 10
                   (Swanson ratio; script 60's own constants are 0.013 / 0.0013)
          rho      0.02 /d primary; 0.01 and 0.04 on the 2 mm grid only (secondary)
          grid     2 mm (64^3) and 1 mm (128^3), rho = 0.02 only on the 1 mm grid
Endpoints (all vs the full arm, same D / rho / grid / horizon)
          primary    mask Dice at u >= 0.16 K (T2/FLAIR threshold, Swanson 2008)
          secondary  mask Dice at u >= 0.80 K (T1Gd threshold)
                     elongation sqrt(lambda_max / lambda_min) of the 0.16 mask (3-D PCA)
                     extent: 90th percentile radial distance of the 0.16 mask (mm)
                     L_diff/dx, L_diff = sqrt(2 D_white t)
Crossover horizon (pre-specified): smallest horizon at which primary Dice (full vs
          iso_matched) < 0.95. "Not reached" is reported as such.
Validity  (a) mask touching the domain boundary is flagged (result censored);
          (b) grid convergence: 1 mm full arm downsampled to 2 mm vs 2 mm full arm (Dice 0.16).
Literature audit: script-60 D_white 0.013 and D_gray 0.0013 mm^2/d vs Swanson-type
          values 0.13 and 0.013 mm^2/d (arXiv 2402.02273 and similar use these).
No statistics: the model is deterministic per (D, rho, grid); there is no n.

Output: output/ablation_horizon_crossover.json
"""
from __future__ import annotations

import json
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import torch

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUT_JSON = PROJECT_ROOT / "output" / "ablation_horizon_crossover.json"

from rl.equal_budget_arms import BatchedSolver  # noqa: E402

_spec = spec_from_file_location("s60", PROJECT_ROOT / "src" / "60_baselines_and_ablation.py")
s60 = module_from_spec(_spec)
_spec.loader.exec_module(s60)

HORIZONS = (90, 180, 365, 600, 900)
D_LIST = (0.01, 0.03, 0.1, 0.3)
THR_FLAIR, THR_T1GD = 0.16, 0.80
CROSSOVER_DICE = 0.95
ARMS = ("full", "iso_matched", "x_aligned")
# (voxels per side, rho, role)
CONFIGS = ((64, 0.02, "primary"), (128, 0.02, "primary"), (64, 0.01, "secondary"), (64, 0.04, "secondary"))
LIT_AUDIT = {"script60_D_white_mm2_per_day": s60.D_WHITE_BASE, "script60_D_gray_mm2_per_day": s60.D_GRAY_BASE,
             "swanson_type_D_white_mm2_per_day": 0.13, "swanson_type_D_gray_mm2_per_day": 0.013,
             "swanson_type_rho_per_day": 0.025,
             "ratio_script60_to_literature_white": s60.D_WHITE_BASE / 0.13,
             "ratio_script60_to_literature_gray": s60.D_GRAY_BASE / 0.013}


def make_solver(arm: str, n: int, d_white: float, rho: float):
    s60.D_GRAY_BASE = d_white / 10.0  # read by _build_tensor_field at construction
    base = s60.FastPDESolver(grid_size=(64, 64, 64), dt_pde=s60.DT_PDE_EVAL, rho=rho, D_white=d_white,
                             alpha_sens=1.0, use_dti=(arm != "x_aligned"), use_mechanics=True)
    s = base
    if n != 64:
        s = s60.FastPDESolver(grid_size=(n, n, n), dt_pde=s60.DT_PDE_EVAL, rho=rho, D_white=d_white,
                              alpha_sens=1.0, use_dti=(arm != "x_aligned"), use_mechanics=True)
        rep = n // 64
        for f in ("D_xx", "D_yy", "D_zz", "D_xy", "D_xz", "D_yz"):  # same physical tensor field
            setattr(s, f, np.repeat(np.repeat(np.repeat(getattr(base, f), rep, 0), rep, 1), rep, 2))
    if arm == "iso_matched":
        md = (s.D_xx + s.D_yy + s.D_zz) / 3.0
        s.D_xx = s.D_yy = s.D_zz = md
        s.D_xy = np.zeros_like(md)
        s.D_xz = np.zeros_like(md)
        s.D_yz = np.zeros_like(md)
    s._compute_face_diffusivities()
    return s


def dice(a: np.ndarray, b: np.ndarray) -> float:
    tot = a.sum() + b.sum()
    return float(2 * (a & b).sum() / tot) if tot > 0 else float("nan")


def elongation(mask: np.ndarray) -> float:
    pts = np.argwhere(mask).astype(float)
    if len(pts) < 10:
        return float("nan")
    ev = np.linalg.eigvalsh(np.cov(pts.T))
    return float(np.sqrt(ev[-1] / max(ev[0], 1e-12)))


def extent_p90(mask: np.ndarray, dx: float) -> float:
    if not mask.any():
        return float("nan")
    c = (mask.shape[0] - 1) / 2.0
    pts = np.argwhere(mask) - c
    return float(np.percentile(np.sqrt((pts ** 2).sum(1)), 90) * dx)


def touches_boundary(mask: np.ndarray) -> bool:
    return bool(mask[0].any() or mask[-1].any() or mask[:, 0].any() or mask[:, -1].any()
                or mask[:, :, 0].any() or mask[:, :, -1].any())


def run_cell(n: int, d_white: float, rho: float):
    solvers = [make_solver(a, n, d_white, rho) for a in ARMS]
    bs = BatchedSolver(solvers, np.zeros((len(ARMS), 4)))
    zero = torch.zeros(len(ARMS), dtype=torch.long)
    dx = bs.dx
    snaps, ds_full = {}, {}
    for day in range(1, max(HORIZONS) + 1):
        bs.step(zero)
        if day in HORIZONS:
            u = bs.u.numpy()
            m = u >= THR_FLAIR
            h = u >= THR_T1GD
            rec = {"L_diff_mm": float(np.sqrt(2 * d_white * day)), "dx_mm": dx}
            rec["L_diff_over_dx"] = rec["L_diff_mm"] / dx
            for i, a in enumerate(ARMS):
                rec[a] = {"volume_flair_mm3": float(m[i].sum() * dx ** 3), "elongation": elongation(m[i]),
                          "extent_p90_mm": extent_p90(m[i], dx), "touches_boundary": touches_boundary(m[i]),
                          "n_flair_voxels": int(m[i].sum()), "n_t1gd_voxels": int(h[i].sum())}
                if a != "full":
                    rec[a]["dice_flair_vs_full"] = dice(m[i], m[0])
                    rec[a]["dice_t1gd_vs_full"] = dice(h[i], h[0])
                    rec[a]["elongation_diff_vs_full"] = rec[a]["elongation"] - rec["full"]["elongation"]
                    rec[a]["extent_diff_vs_full_mm"] = rec[a]["extent_p90_mm"] - rec["full"]["extent_p90_mm"]
            snaps[day] = rec
            r = n // 64
            ds_full[day] = u[0].reshape(64, r, 64, r, 64, r).mean(axis=(1, 3, 5)) if r > 1 else u[0].copy()
    return snaps, ds_full


def main():
    t0 = time.time()
    res = {"script": "78_ablation_horizon_crossover", "horizons_d": HORIZONS, "D_white_list": D_LIST,
           "arms": ARMS, "thresholds_K": {"flair": THR_FLAIR, "t1gd": THR_T1GD},
           "crossover_dice_threshold": CROSSOVER_DICE, "literature_audit": LIT_AUDIT, "cells": {}}
    full_fields = {}
    for n, rho, role in CONFIGS:
        for d in D_LIST:
            key = f"n{n}_rho{rho}_D{d}"
            snaps, ds = run_cell(n, d, rho)
            res["cells"][key] = {"grid_voxels": n, "rho": rho, "D_white": d, "role": role,
                                 "snapshots": {str(k): v for k, v in snaps.items()}}
            full_fields[(n, rho, d)] = ds
            cross = next((h for h in HORIZONS if snaps[h]["iso_matched"]["dice_flair_vs_full"] < CROSSOVER_DICE), None)
            res["cells"][key]["crossover_horizon_d"] = cross
            print(f"[{time.time() - t0:7.0f}s] {key}: crossover {cross} | "
                  + " ".join(f"{h}d dice {snaps[h]['iso_matched']['dice_flair_vs_full']:.3f}" for h in HORIZONS),
                  flush=True)
            OUT_JSON.write_text(json.dumps(res, indent=2))
    conv = {}
    for rho in (0.02,):
        for d in D_LIST:
            a, b = full_fields[(64, rho, d)], full_fields[(128, rho, d)]
            conv[f"rho{rho}_D{d}"] = {str(h): dice(a[h] >= THR_FLAIR, b[h] >= THR_FLAIR) for h in HORIZONS}
    res["grid_convergence_dice_1mm_downsampled_vs_2mm"] = conv
    res["elapsed_s"] = time.time() - t0
    OUT_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "cells"}, indent=1))
    print(f"[saved] {OUT_JSON}")


if __name__ == "__main__":
    main()
