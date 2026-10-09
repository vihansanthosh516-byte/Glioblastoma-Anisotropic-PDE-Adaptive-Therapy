#!/usr/bin/env python3
"""Script 107: grid front-speed correction at h = 2 mm (GRAND_PLAN L2b; weak point W22). Numerics only, no patient data.

Problem (script 97 V7b/V9): at h = 2 mm the front width xi = sqrt(D/rho) is often below one voxel, and the discrete
front then moves faster than the converged one (up to +50%). So a (D, rho) cell at 2 mm does not mean the same
physical tumour as at fine h.
Fix: for each (D, rho), find s* such that the 2 mm front with D_eff = s* D moves at the speed of the h = 0.25 mm front
with D (bisection on s in [0.05, 2]; the discrete speed rises with D). Same slab test and time windows as V7b.
Declared before the run: corrected cells pass if |speed error| < 10% (V7b rule) AND |sphere radius error at 70 d| < 1 mm
against the h = 0.5 mm radius (V9 rule; reference radii read from output/solver_verification.json, same code path).
Use: the table maps a fitted 2 mm cell to its physical (D, rho). Production forecasts are NOT changed here: grid search
picks cells by Dice, so the correction changes what fitted parameters mean, not forecast skill. Using D_eff inside the
fit would need an amendment before any test-set use.
Output: output/front_speed_correction.json (+ manifest)
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
import time
from pathlib import Path as _Path

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from src.run_manifest import write_run_manifest  # noqa: E402

spec = importlib.util.spec_from_file_location("m97", PROJECT_ROOT / "src" / "97_solver_verification.py")
m97 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m97)

H = 2.0


def speed_at(D, rho, h, c_ref):
    """m97.speed with the time windows of the TARGET (D, rho) kept fixed, so s changes only the diffusion."""
    import numpy as np
    N = int(round(90.0 / h))
    shape = (N, 4, 4)
    solver = m97.ria.TensorFK(m97.iso6(shape, D), np.ones(shape, bool), h=h)
    block = int(round(8.0 / h))
    u0 = np.zeros(shape, np.float32)
    u0[:block] = 1.0
    t1, t2 = 15.0 / c_ref, 45.0 / c_ref
    u1, _ = solver.run(u0, [rho], t1)
    u2, _ = solver.run(u0, [rho], t2)
    p1, p2 = m97.front_x(u1[0, :, 1, 1], h, block), m97.front_x(u2[0, :, 1, 1], h, block)
    return (p2 - p1) / (t2 - t1)


def main():
    t0 = time.time()
    ref = json.loads((OUT / "solver_verification.json").read_text())
    v7 = {(r["D"], r["rho"]): r for r in ref["V7_production_grid"]["rows"]}
    v9 = {(r["D"], r["rho"]): r for r in ref["V9_forecast_horizon"]["rows"]}
    rows = []
    for (D, rho), r7 in v7.items():
        c = 2 * math.sqrt(D * rho)
        target = r7[f"speed_h{m97.H_REF:g}"]
        lo, hi = 0.05, 2.0
        for _ in range(18):
            mid = math.sqrt(lo * hi)
            v = speed_at(mid * D, rho, H, c)
            if math.isnan(v) or v > target:   # nan = front left the slab, i.e. too fast
                hi = mid
            else:
                lo = mid
        s = math.sqrt(lo * hi)
        sp = speed_at(s * D, rho, H, c)
        rad = m97.sphere_radius(s * D, rho, H)
        r05 = v9[(D, rho)]["radius_h0.5_mm"]
        row = {"D": D, "rho": rho, "h_over_xi": H / math.sqrt(D / rho), "s_star": s, "D_eff": s * D,
               "speed_err_before": r7["V7b_rel_err_h2_vs_ref"], "speed_err_after": sp / target - 1,
               "radius_err_before_mm": v9[(D, rho)]["err_h2_mm"], "radius_err_after_mm": rad - r05}
        row["pass_after"] = bool(abs(row["speed_err_after"]) < 0.10 and abs(row["radius_err_after_mm"]) < 1.0)
        rows.append(row)
        print(f"[{time.time() - t0:5.0f}s] D {D} rho {rho}: s* {s:.3f} speed {row['speed_err_before']:+.3f} -> "
              f"{row['speed_err_after']:+.3f}  radius {row['radius_err_before_mm']:+.2f} -> {row['radius_err_after_mm']:+.2f} mm", flush=True)
    res = {"script": "107_front_speed_correction", "h_mm": H, "reference": "h = 0.25 mm speed (V7b), h = 0.5 mm radius (V9)",
           "rows": rows, "n_cells": len(rows),
           "n_pass_before": int(sum(abs(r_["speed_err_before"]) < 0.10 and abs(r_["radius_err_before_mm"]) < 1.0 for r_ in rows)),
           "n_pass_after": int(sum(r_["pass_after"] for r_ in rows)),
           "weakest_points": [
               "Correction is fitted on a planar front (speed) and checked on a sphere (radius); curved, anisotropic fronts in real anatomy are not checked.",
               "One scalar per cell corrects speed, not front shape; the 2 mm front is still steeper than the converged one.",
               "Isotropic tensor only; it does not test the off-diagonal (Selling) stencil."]}
    (OUT / "front_speed_correction.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("front_speed_correction_107", OUT / "front_speed_correction.manifest.json",
                       script="src/107_front_speed_correction.py", seed=0, config={"h_mm": H, "bisection_steps": 18},
                       inputs=[OUT / "solver_verification.json"], dataset="none (numerics)", patient_split="none",
                       primary_endpoint="cells passing V7b and V9 after correction")
    print(json.dumps({k: res[k] for k in ("n_cells", "n_pass_before", "n_pass_after")}))


if __name__ == "__main__":
    main()
