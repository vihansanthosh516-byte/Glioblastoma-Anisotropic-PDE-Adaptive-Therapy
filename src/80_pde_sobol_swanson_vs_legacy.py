#!/usr/bin/env python3
"""
Script 80: PDE-based Sobol sensitivity at legacy vs Swanson-type diffusivities.
==============================================================================
Replaces the archived script-46 number (S1(rho) = 0.998), which came from a reduced ODE with an
arbitrary coupling k_diff = 15 and D_white bounded to +/-20% around 0.013 mm^2/d. Here the
script-60 3-D anisotropic PDE solver (64^3, 2 mm) is run for every Sobol sample.

PRE-SPECIFIED DESIGN (fixed before the first run; see vault/Script-60-66-Swanson-D.md)
--------------------------------------------------------------------------------------
Parameters  rho in (0.005, 0.035) /d, D_w, alpha_sens in (0.5, 1.5) (script-59 PARAM_RANGES)
Regimes     legacy   D_w in (0.001, 0.008) mm^2/d, D_gray = 0.0013
            swanson  D_w in (0.1,   0.8)   mm^2/d, D_gray = 0.013
Sampling    SALib Saltelli (deterministic Sobol sequence), N = 128 base samples, first-/total-order,
            k = 3 -> N (k + 2) = 640 PDE runs per regime. Same sample design in both regimes.
Schedule    Stupp (script 66 STUPP_SCHEDULE), 90 days, script-60 full model (DTI + mechanics flag)
Outputs Y   V_total   day-90 total volume (script 60's primary endpoint)   <- PRIMARY
            V_vis_0.04  voxels with u >= 0.04 K (imaging-visible volume)
            extent_p90_mm  90th percentile radial distance of voxels with u >= 0.01 K
Report      S1 and ST with bootstrap CI (SALib default, 95%) per parameter, per output, per regime.
Output: output/sobol_pde_<regime>.json, output/sobol_pde_swanson_vs_legacy.json (compare)
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUT_DIR = PROJECT_ROOT / "output"

_spec = spec_from_file_location("s77", PROJECT_ROOT / "src" / "77_ablation_spatial_endpoints.py")
s77 = module_from_spec(_spec)
_spec.loader.exec_module(s77)
s60, s66 = s77.s60, s77.s66

N_BASE = 128
SEED = 42
REGIMES = {"legacy": {"D_w": (0.001, 0.008), "D_gray": 0.0013},
           "swanson": {"D_w": (0.1, 0.8), "D_gray": 0.013}}
OUTPUTS = ("V_total", "V_vis_0.04", "extent_p90_mm")


def problem(regime: str) -> dict:
    return {"num_vars": 3, "names": ["rho", "D_w", "alpha_sens"],
            "bounds": [[0.005, 0.035], list(REGIMES[regime]["D_w"]), [0.5, 1.5]]}


def run_regime(regime: str) -> None:
    from SALib.analyze import sobol
    from SALib.sample import saltelli
    prob = problem(regime)
    X = saltelli.sample(prob, N_BASE, calc_second_order=False)
    s60.D_GRAY_BASE = REGIMES[regime]["D_gray"]          # read by FastPDESolver at construction
    flags = s60.ABLATIONS[s60.FULL]
    scenarios = [{"rho": float(r), "D_w": float(d), "alpha_sens": float(a)} for r, d, a in X]
    Y = {k: np.zeros(len(X)) for k in OUTPUTS}
    t0 = time.time()
    for i in range(0, len(X), s77.CHUNK):
        sc = scenarios[i:i + s77.CHUNK]
        u, dx = s77.final_fields(sc, flags, True, s66.STUPP_SCHEDULE)
        e = s77.endpoints(u, dx)
        Y["V_total"][i:i + len(sc)] = e["V_total"]
        Y["V_vis_0.04"][i:i + len(sc)] = e["V_vis_0.04"]
        Y["extent_p90_mm"][i:i + len(sc)] = e["extent_p90_mm"]
        print(f"[{time.time() - t0:6.0f}s] {regime} {i + len(sc)}/{len(X)}", flush=True)
    res = {"regime": regime, "n_base": N_BASE, "n_runs": len(X), "problem": prob, "seed": SEED,
           "D_gray": REGIMES[regime]["D_gray"], "stupp_budget": s66.BUDGET, "outputs": {}}
    for k in OUTPUTS:
        si = sobol.analyze(prob, Y[k], calc_second_order=False, seed=SEED, print_to_console=False)
        res["outputs"][k] = {
            "S1": dict(zip(prob["names"], map(float, si["S1"]))),
            "S1_conf": dict(zip(prob["names"], map(float, si["S1_conf"]))),
            "ST": dict(zip(prob["names"], map(float, si["ST"]))),
            "ST_conf": dict(zip(prob["names"], map(float, si["ST_conf"]))),
            "mean": float(Y[k].mean()), "std": float(Y[k].std()),
            "min": float(Y[k].min()), "max": float(Y[k].max())}
    res["samples"] = {"X": X.tolist(), **{k: Y[k].tolist() for k in OUTPUTS}}
    (OUT_DIR / f"sobol_pde_{regime}.json").write_text(json.dumps(res, indent=2))
    print(json.dumps({k: res[k] for k in ("regime", "n_runs", "outputs")}, indent=1))
    print(f"[saved] output/sobol_pde_{regime}.json")


def compare() -> None:
    a = json.loads((OUT_DIR / "sobol_pde_legacy.json").read_text())
    b = json.loads((OUT_DIR / "sobol_pde_swanson.json").read_text())
    out = {"archived_script46_reduced_ode_S1": json.loads(
        (OUT_DIR / "sobol_sensitivity_results.json").read_text())["S1"],
        "legacy": {k: v for k, v in a["outputs"].items()}, "swanson": {k: v for k, v in b["outputs"].items()}}
    (OUT_DIR / "sobol_pde_swanson_vs_legacy.json").write_text(json.dumps(out, indent=2))
    for k in OUTPUTS:
        print(k)
        for name in ("rho", "D_w", "alpha_sens"):
            print(f"  {name:11s} legacy S1 {a['outputs'][k]['S1'][name]:+.3f} ST {a['outputs'][k]['ST'][name]:+.3f} | "
                  f"swanson S1 {b['outputs'][k]['S1'][name]:+.3f} ST {b['outputs'][k]['ST'][name]:+.3f}")
    print("[saved] output/sobol_pde_swanson_vs_legacy.json")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("regime", choices=["legacy", "swanson", "compare"])
    r = ap.parse_args().regime
    compare() if r == "compare" else run_regime(r)
