#!/usr/bin/env python3
"""
Script 77: Is the script-60 "DTI/mechanics/diffusion < 0.01%" ablation a null
finding, or a property of the endpoint?
============================================================================
Track B negative 2. Script 60 scores every ablation by total tumour volume
V = sum(u) dx^3. With no-flux boundaries the diffusion term is conservative,
d/dt sum(u) = sum(rho u (1 - u/K) - k u), so diffusion (and therefore the tensor
D and its DTI orientation) can change V only through the logistic term
u/K. The script-59 seed peaks at u = 0.08 K, so that coupling is weak and the
ablation is close to invariant by construction (script 66 validate.json: the
largest diffusion effect on Stupp's day-90 V is 0.029%).

A test that can see D needs an endpoint that depends on where tumour cells
are. Here the script-60 solver (its own physics and alpha_sens kill rates, 30
LHS scenarios) is run under three configurations:
    full        script-60 full model (DTI tensor)
    no_dti      script-60 use_dti=False (uniform x-aligned tensor, as in script 60)
    no_diff     diffusion switched off (reaction only)
and two fixed equal-budget schedules (Stupp; combo on days 56-90, the script-66
oracle). Endpoints on day 90:
    V_total     script 60's endpoint (sum u)
    V_vis(t)    "visible" volume, voxels with u >= t, t in {0.01, 0.04} K
                (imaging sees a density threshold, Swanson et al. 2008)
    extent      90th-percentile distance from the seed centre of voxels with u >= 0.01 K
    dice_vs_full Dice of the u >= 0.01 K mask against the full model
Also: does the Stupp-vs-oracle ranking change under V_vis?

Also reported: the diffusion length sqrt(2 D t) and Fisher front travel
2 sqrt(D rho) t over the 90 days, against the 2 mm voxel.

Output: output/ablation_spatial_endpoints.json
"""
from __future__ import annotations

import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
import torch

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUT_JSON = PROJECT_ROOT / "output" / ("ablation_spatial_endpoints.json" if __import__("os").environ.get("GBM_D_REGIME") == "legacy"
                                     else "ablation_spatial_endpoints_swanson.json")  # default regime: Swanson D (script 59)

from rl.equal_budget_arms import BatchedSolver, s66  # noqa: E402

_spec = spec_from_file_location("s60", PROJECT_ROOT / "src" / "60_baselines_and_ablation.py")
s60 = module_from_spec(_spec)
_spec.loader.exec_module(s60)

THRESHOLDS = (0.01, 0.04)
CHUNK = 10
SCHEDULES = {"stupp": s66.STUPP_SCHEDULE, "combo_days56_90": s66.block_schedule(3, 56, 35)}
CONFIGS = {"full": (s60.ABLATIONS[s60.FULL], True),
           "no_dti": (s60.ABLATIONS["No DTI (x-aligned tensor)"], True),
           "no_diff": (s60.ABLATIONS[s60.FULL], False)}


def final_fields(scenarios, flags, diffuse: bool, sched) -> np.ndarray:
    out = []
    for i in range(0, len(scenarios), CHUNK):
        sc = scenarios[i:i + CHUNK]
        solvers = [s60.FastPDESolver(grid_size=s60.EVAL_GRID, dt_pde=s60.DT_PDE_EVAL, rho=p["rho"],
                                     D_white=p["D_w"], alpha_sens=p["alpha_sens"], **flags) for p in sc]
        bs = BatchedSolver(solvers, [s60.kill_row(p["alpha_sens"]) for p in sc])
        if not diffuse:
            bs._div = lambda u: torch.zeros_like(u)
        s66.run_schedules(bs, torch.tensor([sched] * len(sc)))
        out.append(bs.u.numpy().copy())
        dx = bs.dx
    return np.concatenate(out), dx


def endpoints(u: np.ndarray, dx: float) -> dict:
    n = u.shape[1]
    c = n // 2
    zz, yy, xx = np.meshgrid(*[np.arange(n)] * 3, indexing="ij")
    r = np.sqrt((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2) * dx
    out = {"V_total": (u.sum(axis=(1, 2, 3)) * dx ** 3).tolist()}
    for t in THRESHOLDS:
        out[f"V_vis_{t}"] = ((u >= t).sum(axis=(1, 2, 3)) * dx ** 3).tolist()
    ext = []
    for b in range(u.shape[0]):
        m = u[b] >= THRESHOLDS[0]
        ext.append(float(np.percentile(r[m], 90)) if m.any() else 0.0)
    out["extent_p90_mm"] = ext
    return out


def rel(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = b > 0
    return float(np.max(np.abs(a[ok] - b[ok]) / b[ok])) if ok.any() else None


def main():
    scenarios = s66.lhs_scenarios()
    res, masks = {}, {}
    for sname, sched in SCHEDULES.items():
        for cname, (flags, diffuse) in CONFIGS.items():
            u, dx = final_fields(scenarios, flags, diffuse, sched)
            res[(sname, cname)] = endpoints(u, dx)
            masks[(sname, cname)] = u >= THRESHOLDS[0]
            print(f"{sname:16s} {cname:8s} V_total {np.mean(res[(sname, cname)]['V_total']):8.2f}  "
                  f"V_vis0.01 {np.mean(res[(sname, cname)]['V_vis_0.01']):8.1f}  "
                  f"extent {np.mean(res[(sname, cname)]['extent_p90_mm']):5.1f} mm")

    out = {"n_scenarios": len(scenarios), "thresholds_K": THRESHOLDS, "schedules": list(SCHEDULES),
           "configs": {k: {"flags": v[0], "diffusion": v[1]} for k, v in CONFIGS.items()},
           "per_schedule": {}}
    for sname in SCHEDULES:
        full = res[(sname, "full")]
        blk = {}
        for cname in CONFIGS:
            e = res[(sname, cname)]
            mf, mc = masks[(sname, "full")], masks[(sname, cname)]
            inter = (mf & mc).sum(axis=(1, 2, 3))
            tot = mf.sum(axis=(1, 2, 3)) + mc.sum(axis=(1, 2, 3))
            dice = np.where(tot > 0, 2 * inter / np.maximum(tot, 1), 1.0)
            blk[cname] = {
                "mean": {k: float(np.mean(v)) for k, v in e.items()},
                "max_rel_change_vs_full": {k: rel(e[k], full[k]) for k in e} if cname != "full" else None,
                "mask_dice_vs_full_mean": float(dice.mean()), "mask_dice_vs_full_min": float(dice.min()),
                "per_scenario": e}
        out["per_schedule"][sname] = blk
    rank = {}
    for cname in CONFIGS:
        for ep in ["V_total"] + [f"V_vis_{t}" for t in THRESHOLDS]:
            a = np.array(res[("combo_days56_90", cname)][ep])
            b = np.array(res[("stupp", cname)][ep])
            rank[f"{cname}:{ep}"] = {"oracle_smaller_pct": float((a < b).mean() * 100),
                                     "mean_log_ratio_oracle_vs_stupp":
                                         float(np.mean(np.log(np.maximum(a, 1e-9) / np.maximum(b, 1e-9))))}
    out["ranking_oracle_vs_stupp"] = rank
    # Length scales: can diffusion move cells across even one voxel in 90 days?
    Dw = np.array([p["D_w"] for p in scenarios])
    rho = np.array([p["rho"] for p in scenarios])
    days = len(s66.STUPP_SCHEDULE)
    out["length_scales_mm"] = {
        "dx": float(128.0 / s60.EVAL_GRID[0]),
        "diffusion_length_sqrt_2Dt_min_max": [float(np.sqrt(2 * Dw * days).min()),
                                              float(np.sqrt(2 * Dw * days).max())],
        "fisher_front_travel_2sqrt_Drho_t_min_max": [float((2 * np.sqrt(Dw * rho) * days).min()),
                                                     float((2 * np.sqrt(Dw * rho) * days).max())],
        "days": days}
    print(json.dumps(out["length_scales_mm"], indent=1))
    OUT_JSON.write_text(json.dumps(out, indent=2))
    print(json.dumps({s: {c: {k: v for k, v in d.items() if k != "per_scenario"} for c, d in b.items()}
                      for s, b in out["per_schedule"].items()}, indent=1))
    print(json.dumps(rank, indent=1))
    print(f"[saved] {OUT_JSON}")


if __name__ == "__main__":
    main()
