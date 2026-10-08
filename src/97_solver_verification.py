#!/usr/bin/env python3
"""Script 97: numerical verification of the tumour PDE solver (Phase 2.1; plan A1.19, masterplan §9, §105).

Verification asks "are the equations solved correctly?" It does not ask whether they describe a tumour.
Solver under test: run_improved_aniso.TensorFK (finite volume, explicit Euler, zero flux outside the domain mask).
It already has a short selftest (output/forecast_validation/solver_selftest.json). This script adds what
the selftest lacks: convergence orders, boundedness, boundary behaviour, init sensitivity, and the
resolution of the production grid.

PRE-DECLARED CRITERIA (fixed before the first run of this script)
  V1 spatial order      diffusion of a cosine mode, exact solution, N = 8, 16, 32 cells on 64 mm.
                        Observed order from successive refinements must be in [1.8, 2.2].
  V2 temporal order     Fisher-KPP, max_dt = 0.4, 0.2, 0.1, 0.05 vs reference 0.0125.
                        Observed order must be in [0.8, 1.2] (explicit Euler).
  V3 grid convergence   front position at t = 300 d vs finest grid (h = 0.75 mm): error must fall
                        monotonically over h = 3, 1.5. Reported, no order claimed (pulled fronts).
  V4 conservation       isotropic, rho = 0: relative mass change < 1e-5.
  V5 boundedness        random u0 in [0, 1], rho = 0.1, 0 <= u <= 1 and finite.
  V6 boundary           mass outside the domain mask exactly 0; interior mass conserved (isotropic, rho = 0).
  V7 production grid    h = 2 mm, D in {0.01, 0.03, 0.1, 0.3} mm2/d, rho in {0.01, 0.03, 0.1} /d.
                        V7a (declared): Fisher front speed vs 2 sqrt(D rho), verified if |error| < 10%.
                        FOUND INVALID on first run for thick fronts (see v7_production_grid docstring).
                        V7b (added after the first run, before its results were read): speed at h = 2 mm vs
                        h = 0.25 mm, verified if |error| < 10%. Both are reported.
                        Consequence rule: forecast fits that land in a V7b-unverified cell are reported separately.
  V9 forecast horizon   what the forecast actually uses: a sphere (R = 10 mm) grown for 70 d (the median MU gap), on the
                        production grid h = 2 mm vs h = 0.5 mm, same (D, rho) grid as V7. Compared by radius gain from the
                        u > 0.5 volume. Declared criterion: |radius(h = 2) - radius(h = 0.5)| < 1 mm (half a voxel).
                        Octant of the sphere with the corner at the centre; the array edges are the symmetry planes.
  V8 init sensitivity   +/- 1 voxel (2 mm) change of the initial radius: Dice of outputs, reported (no pass/fail).
V1-V6 failing means STOP: no biological number is reported until fixed (analysis_plan_v1.md A1.19).

Output: output/solver_verification.json, output/solver_verification.manifest.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path as _Path

import numpy as np
import torch

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
sys.path.insert(0, str(PROJECT_ROOT))

import run_improved_aniso as ria  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

SEED = 20261004
torch.set_num_threads(2)


def iso6(shape, D):
    return np.tile(np.array([D, D, D, 0, 0, 0], np.float32), tuple(shape) + (1,))


def front_x(u_line, h, start):
    """Position (mm) where the centreline profile first falls through 0.5 beyond index `start`, linear interpolation."""
    idx = np.where(u_line[start:] < 0.5)[0]
    if len(idx) == 0:
        return float("nan")
    i = start + int(idx[0])
    if i == 0:
        return 0.5 * h
    u0, u1 = float(u_line[i - 1]), float(u_line[i])
    frac = (u0 - 0.5) / max(u0 - u1, 1e-12)
    return ((i - 1) + 0.5 + frac) * h


def v1_spatial_order():
    L, D, t_end, amp = 64.0, 0.1, 2000.0, 0.4   # base 0.5 +/- 0.4 keeps u in [0.1, 0.9]; the clamp never acts
    rows = []
    for N in (8, 16, 32, 64):
        h = L / N
        shape = (N, 4, 4)
        solver = ria.TensorFK(iso6(shape, D), np.ones(shape, bool), h=h)
        x = (np.arange(N) + 0.5) * h
        u0 = np.broadcast_to((0.5 + amp * np.cos(math.pi * x / L))[:, None, None], shape).astype(np.float32)
        u, n = solver.run(u0, [0.0], t_end)
        exact = 0.5 + amp * np.cos(math.pi * x / L) * math.exp(-D * (math.pi / L) ** 2 * t_end)
        err = float(np.abs(u[0, :, 1, 1] - exact).max())
        rows.append({"N": N, "h_mm": h, "max_abs_err": err, "steps": n})
    for a, b in zip(rows, rows[1:]):
        b["observed_order_vs_previous"] = math.log2(a["max_abs_err"] / b["max_abs_err"])
    orders = [r["observed_order_vs_previous"] for r in rows[1:3]]
    return {"rows": rows, "orders_used": orders, "pass": all(1.8 <= o <= 2.2 for o in orders),
            "note": "N = 64 is reported but not used: float32 round-off is comparable to the truncation error there."}


def v2_temporal_order():
    shape, h, D, rho, t_end = (60, 4, 4), 2.0, 0.1, 0.05, 100.0
    x = (np.arange(shape[0]) + 0.5) * h
    u0 = np.broadcast_to((0.5 * (1 - np.tanh((x - 40) / 6)))[:, None, None], shape).astype(np.float32)
    solver = ria.TensorFK(iso6(shape, D), np.ones(shape, bool), h=h)
    ref, _ = solver.run(u0, [rho], t_end, max_dt=0.0125)
    rows = []
    for dt in (0.4, 0.2, 0.1, 0.05):
        u, n = solver.run(u0, [rho], t_end, max_dt=dt)
        rows.append({"max_dt": dt, "steps": n, "max_abs_err_vs_ref": float(np.abs(u - ref).max())})
    for a, b in zip(rows, rows[1:]):
        b["observed_order_vs_previous"] = math.log2(a["max_abs_err_vs_ref"] / b["max_abs_err_vs_ref"])
    orders = [r["observed_order_vs_previous"] for r in rows[1:3]]
    return {"rows": rows, "orders_used": orders, "pass": all(0.8 <= o <= 1.2 for o in orders)}


def v3_grid_convergence():
    L, D, rho, t_end = 96.0, 0.1, 0.05, 300.0
    pos = {}
    for N in (16, 32, 64, 128):
        h = L / N
        shape = (N, 4, 4)
        x = (np.arange(N) + 0.5) * h
        u0 = np.broadcast_to((0.5 * (1 - np.tanh((x - 24) / 6)))[:, None, None], shape).astype(np.float32)
        solver = ria.TensorFK(iso6(shape, D), np.ones(shape, bool), h=h)
        u, _ = solver.run(u0, [rho], t_end)
        pos[h] = front_x(u[0, :, 1, 1], h, 0)
    fine = pos[L / 128]
    rows = [{"h_mm": h, "front_mm": p, "abs_err_vs_finest": abs(p - fine)} for h, p in pos.items()]
    e = [r["abs_err_vs_finest"] for r in rows[:3]]
    return {"xi_mm": math.sqrt(D / rho), "rows": rows, "pass": e[1] < e[0] and e[2] < e[1]}


def v4_conservation():
    shape = (30, 28, 26)
    zz, yy, xx = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    mask = ((xx - 14) ** 2 + (yy - 14) ** 2 + (zz - 13) ** 2) < 12 ** 2
    solver = ria.TensorFK(iso6(shape, 0.3), mask, h=2.0)
    rng = np.random.default_rng(SEED)
    u0 = (rng.random(shape) * mask).astype(np.float32)
    u, _ = solver.run(u0, [0.0], 200.0)
    rel = float(abs(u[0].sum(dtype=np.float64) - u0.sum(dtype=np.float64)) / u0.sum(dtype=np.float64))
    return {"rel_mass_change": rel, "pass": rel < 1e-5}


def v5_boundedness():
    shape = (30, 28, 26)
    zz, yy, xx = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    mask = ((xx - 14) ** 2 + (yy - 14) ** 2 + (zz - 13) ** 2) < 12 ** 2
    rng = np.random.default_rng(SEED)
    A = rng.normal(size=shape + (3, 3)) * 0.3
    T = A @ np.swapaxes(A, -1, -2) + 0.1 * np.eye(3)
    d6 = ria.mat_to_six(T) * 0.3
    out = {}
    for name, d in (("isotropic", iso6(shape, 0.3)), ("random_spd_tensor", d6)):
        solver = ria.TensorFK(d, mask, h=2.0)
        u0 = (rng.random(shape) * mask).astype(np.float32)
        u, _ = solver.run(u0, [0.1], 100.0)
        out[name] = {"min": float(u.min()), "max": float(u.max()), "finite": bool(np.isfinite(u).all())}
    out["pass"] = all(v["min"] >= 0 and v["max"] <= 1 and v["finite"] for v in out.values() if isinstance(v, dict))
    # How much mass does the clamp add? Strong anisotropy, sharp edge, rho = 0.
    Tm = np.tile(ria.mat_to_six(np.array([[0.6, 0.5, 0.0], [0.5, 0.6, 0.0], [0.0, 0.0, 0.05]])).astype(np.float32),
                 (24, 24, 24, 1))
    solver = ria.TensorFK(Tm * 0.3, np.ones((24, 24, 24), bool), h=2.0)
    u0 = np.zeros((24, 24, 24), np.float32)
    u0[10:14, 10:14, 10:14] = 1.0
    u, _ = solver.run(u0, [0.0], 100.0)
    out["clamp_mass_gain_strong_anisotropy"] = float((u[0].sum(dtype=np.float64) - u0.sum(dtype=np.float64)) / u0.sum(dtype=np.float64))
    return out


def v6_boundary():
    shape = (30, 28, 26)
    zz, yy, xx = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    mask = ((xx - 14) ** 2 + (yy - 14) ** 2 + (zz - 13) ** 2) < 10 ** 2
    solver = ria.TensorFK(iso6(shape, 0.3), mask, h=2.0)
    u0 = np.zeros(shape, np.float32)
    u0[10:19, 10:19, 9:18] = 1.0
    u0 *= mask
    u, _ = solver.run(u0, [0.0], 400.0)
    outside = float(np.abs(u[0][~mask]).sum(dtype=np.float64))
    rel = float(abs(u[0].sum(dtype=np.float64) - u0.sum(dtype=np.float64)) / u0.sum(dtype=np.float64))
    # mass reaches the wall: fraction of mass in the outermost interior shell grows but never leaves
    return {"mass_outside_mask": outside, "interior_rel_mass_change": rel, "pass": outside == 0.0 and rel < 1e-5}


def dice(a, b):
    return float(2 * (a & b).sum() / max(a.sum() + b.sum(), 1))


def v8_init_sensitivity():
    shape = (40, 40, 40)
    zz, yy, xx = np.meshgrid(*[np.arange(s) for s in shape], indexing="ij")
    r2 = (xx - 20) ** 2 + (yy - 20) ** 2 + (zz - 20) ** 2
    solver = ria.TensorFK(iso6(shape, 0.1), np.ones(shape, bool), h=2.0)
    outs = {}
    for r in (4, 5, 6):
        u0 = (r2 < r * r).astype(np.float32)
        u, _ = solver.run(u0, [0.03], 60.0)
        outs[r] = (u[0] > ria.U_VISIBLE, u0 > 0)
    return {"input_dice_r4_vs_r5": dice(outs[4][1], outs[5][1]), "input_dice_r6_vs_r5": dice(outs[6][1], outs[5][1]),
            "output_dice_r4_vs_r5": dice(outs[4][0], outs[5][0]), "output_dice_r6_vs_r5": dice(outs[6][0], outs[5][0]),
            "note": "Reported only. A 2 mm change in the initial radius changes the output by this much."}


def speed(D, rho, h, length_mm=90.0):
    N = int(round(length_mm / h))
    shape = (N, 4, 4)
    solver = ria.TensorFK(iso6(shape, D), np.ones(shape, bool), h=h)
    c = 2 * math.sqrt(D * rho)
    block = int(round(8.0 / h))
    u0 = np.zeros(shape, np.float32)
    u0[:block] = 1.0
    t1, t2 = 15.0 / c, 45.0 / c
    u1, _ = solver.run(u0, [rho], t1)
    u2, _ = solver.run(u0, [rho], t2)
    p1, p2 = front_x(u1[0, :, 1, 1], h, block), front_x(u2[0, :, 1, 1], h, block)
    return (p2 - p1) / (t2 - t1), c


H_REF = 0.25


def v7_production_grid():
    """Two criteria. V7a (as pre-declared): speed vs the asymptotic 2 sqrt(D rho). First run showed V7a is the wrong
    test for thick fronts: a pulled front reaches 2 sqrt(D rho) from below (Bramson delay, about 3/(2 lambda) ln t with
    lambda = sqrt(rho / D)), so a 90 mm slab and 45 mm of travel give -10 to -35% at xi = 3-5 mm on EVERY grid.
    V7b (added after that first run, before looking at its results): discretisation error only,
    speed(h = 2 mm) vs speed(h = 0.25 mm) over the same time windows; verified if |error| < 10%."""
    rows = []
    for D in (0.01, 0.03, 0.1, 0.3):
        for rho in (0.01, 0.03, 0.1):
            s2, c = speed(D, rho, 2.0)
            xi = math.sqrt(D / rho)
            s1, _ = speed(D, rho, 1.0)
            sr, _ = speed(D, rho, H_REF)
            row = {"D": D, "rho": rho, "xi_mm": xi, "h_over_xi": 2.0 / xi, "speed_analytic": c,
                   "speed_h2": s2, "speed_h1": s1, f"speed_h{H_REF:g}": sr,
                   "V7a_rel_err_h2_vs_analytic": s2 / c - 1, "V7a_verified": bool(abs(s2 / c - 1) < 0.10),
                   "V7b_rel_err_h2_vs_ref": s2 / sr - 1, "V7b_rel_err_h1_vs_ref": s1 / sr - 1,
                   "V7b_verified_h2": bool(abs(s2 / sr - 1) < 0.10)}
            rows.append(row)
    return {"rows": rows, "n_cells": len(rows),
            "V7a_n_verified": int(sum(r["V7a_verified"] for r in rows)),
            "V7b_n_verified_h2": int(sum(r["V7b_verified_h2"] for r in rows)),
            "V7b_unverified_h2": [(r["D"], r["rho"]) for r in rows if not r["V7b_verified_h2"]]}


def sphere_radius(D, rho, h, R0=10.0, t_end=70.0, half=40.0):
    N = int(round(half / h))
    shape = (N, N, N)
    c = (np.arange(N) + 0.5) * h
    zz, yy, xx = np.meshgrid(c, c, c, indexing="ij")
    u0 = ((xx ** 2 + yy ** 2 + zz ** 2) < R0 ** 2).astype(np.float32)
    solver = ria.TensorFK(iso6(shape, D), np.ones(shape, bool), h=h)
    u, _ = solver.run(u0, [rho], t_end)
    vol = 8.0 * float((u[0] > ria.U_VISIBLE).sum()) * h ** 3
    return (3 * vol / (4 * math.pi)) ** (1 / 3)


def v9_forecast_horizon():
    rows = []
    for D in (0.01, 0.03, 0.1, 0.3):
        for rho in (0.01, 0.03, 0.1):
            r2, r1, r05 = (sphere_radius(D, rho, h) for h in (2.0, 1.0, 0.5))
            rows.append({"D": D, "rho": rho, "radius_h2_mm": r2, "radius_h1_mm": r1, "radius_h0.5_mm": r05,
                         "err_h2_mm": r2 - r05, "err_h1_mm": r1 - r05, "verified_h2": bool(abs(r2 - r05) < 1.0)})
    return {"rows": rows, "n_verified_h2": int(sum(r["verified_h2"] for r in rows)), "n_cells": len(rows),
            "unverified_h2": [(r["D"], r["rho"]) for r in rows if not r["verified_h2"]]}


def main(solver_name="tensorfk"):
    t0 = time.time()
    suffix = ""
    if solver_name == "monotone":   # GRAND_PLAN #15: same checks on the positivity-preserving solver
        sys.path.insert(0, str(PROJECT_ROOT / "src"))
        from solver_monotone import TensorFKMonotone
        ria.TensorFK = TensorFKMonotone
        suffix = "_monotone"
    res = {"script": "97_solver_verification", "solver": f"run_improved_aniso.TensorFK as {ria.TensorFK.__name__}", "seed": SEED}
    for name, fn in [("V1_spatial_order", v1_spatial_order), ("V2_temporal_order", v2_temporal_order),
                     ("V3_grid_convergence_front", v3_grid_convergence), ("V4_conservation", v4_conservation),
                     ("V5_boundedness", v5_boundedness), ("V6_boundary", v6_boundary),
                     ("V7_production_grid", v7_production_grid), ("V9_forecast_horizon", v9_forecast_horizon),
                     ("V8_init_sensitivity", v8_init_sensitivity)]:
        res[name] = fn()
        print(f"[{time.time() - t0:5.0f}s] {name}: {res[name].get('pass', 'reported')}", flush=True)
    core = ["V1_spatial_order", "V2_temporal_order", "V3_grid_convergence_front", "V4_conservation",
            "V5_boundedness", "V6_boundary"]
    res["core_pass"] = all(res[k]["pass"] for k in core)
    res["production_grid_all_verified"] = res["V7_production_grid"]["V7b_n_verified_h2"] == res["V7_production_grid"]["n_cells"]
    (OUT / f"solver_verification{suffix}.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("solver_verification_97" + suffix, OUT / f"solver_verification{suffix}.manifest.json",
                       script="src/97_solver_verification.py", seed=SEED, config={"h_production_mm": 2.0, "solver": solver_name},
                       inputs=[PROJECT_ROOT / "run_improved_aniso.py", PROJECT_ROOT / "src" / "solver_monotone.py"], dataset="synthetic",
                       primary_endpoint="n/a (verification)")
    print(json.dumps({k: res[k] for k in ("core_pass", "production_grid_all_verified")}))


if __name__ == "__main__":
    main("monotone" if "--solver" in sys.argv and sys.argv[sys.argv.index("--solver") + 1] == "monotone" else "tensorfk")
