#!/usr/bin/env python3
"""
Script 155: Track C wrong-model suites for the belief-containment rule.
=====================================================================
SIMULATION ONLY. The policies keep the specs chosen in script 154 (tuned on the nominal
training draw). The filter keeps its nominal model and priors. Only the TRUTH changes:

  M0 nominal        same generator as script 154 (new seed; reference)
  M1 competition    asymmetric competition a_sr = 2, a_rs = 0.5 (filter assumes 1, 1)
  M2 fast growth    r in [0.05, 0.10], outside the filter prior [0.001, 0.05]
  M3 weak drug      m = k/r in [3, 6], outside the prior [6, 20]
  M4 noisy readings true noise sd = 2 x the sd the filter is told
  M5 heavy tails    Student-t readings (df = 3), same sd
  M6 reading bias   every reading 10% high (log offset +0.095); ratio rules are unaffected,
                    the filter's absolute burden and level are not

Primary per suite: mean capped TTP, containment minus the script-154 primary comparator
(raw on/off 0.6/0.5), paired bootstrap. Also vs front-loaded and the other arms.
Seeds: suite k uses scenario seed 20261060 + k and noise seeds 1_060_00k*, all new.

Run:  python src/155_trackc_wrong_models.py --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_wrong_models"
S154 = PROJECT_ROOT / "output" / "trackc_containment" / "containment_results.json"
ARMS = ("risk_by_mean", "base_best_by_mean", "base_frontloaded", "base_raw_single", "pilot2_hybrid", "base_zhang")


def suite(k, n):
    seed = 20261060 + k
    if k == 0:
        return "M0_nominal", tc.draw_scenarios(n, seed)
    if k == 1:
        return "M1_competition", tc.draw_scenarios(n, seed, a_sr=2.0, a_rs=0.5)
    if k == 2:
        return "M2_fast_growth", tc.draw_scenarios(n, seed, ranges=dict(tc.RANGES, r=(0.05, 0.10, "log")))
    if k == 3:
        return "M3_weak_drug", tc.draw_scenarios(n, seed, ranges=dict(tc.RANGES, m=(3.0, 6.0, "log")))
    s = tc.draw_scenarios(n, seed)
    if k == 4:
        s["sigma_true"] = 2.0 * s["sigma"]
        return "M4_noise_2x", s
    if k == 5:
        s["noise_df"] = np.full(n, 3.0)
        return "M5_heavy_tails", s
    s["bias"] = np.full(n, np.log(1.10))
    return "M6_reading_bias", s


def boot_mean_diff(a, b, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    d = a - b
    idx = rng.integers(0, d.size, (n_boot, d.size))
    return float(d.mean()), np.percentile(d[idx].mean(1), [2.5, 97.5]).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    s154 = json.loads(S154.read_text())["arms"]
    specs = {a: tuple(s154[a]["spec"]) for a in ARMS}

    suites = [suite(k, args.n) for k in range(7)]
    jobs, keys = [], []
    for k, (name, scen) in enumerate(suites):
        noise = [1_060_000 + 10 * k + 1, 1_060_000 + 10 * k + 2]
        for spec in dict.fromkeys(specs.values()):
            jobs.append((spec, scen, noise))
            keys.append((name, spec))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(tc.score_spec, jobs))
    ttp = {key: r[1] for key, r in zip(keys, res)}

    out = {}
    for name, _ in suites:
        prim = ttp[(name, specs["base_best_by_mean"])]
        front = ttp[(name, specs["base_frontloaded"])]
        rows = {}
        for a in ARMS:
            x = ttp[(name, specs[a])]
            dm, cim = boot_mean_diff(x, prim, seed=41)
            df, cif = boot_mean_diff(x, front, seed=42)
            dc, cic = tc.paired_bootstrap_cvar(x, prim, seed=43)
            rows[a] = {"spec": list(specs[a]), "mean": float(x.mean()), "cvar10": tc.cvar_lower(x),
                       "mean_diff_vs_primary": dm, "mean_ci95": cim,
                       "mean_diff_vs_frontloaded": df, "frontloaded_ci95": cif,
                       "cvar10_diff_vs_primary": dc, "cvar10_ci95": cic}
        out[name] = rows

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {
        "script": "155_trackc_wrong_models",
        "label": "SIMULATION, assumed parameters; truth deliberately differs from the filter's model in M1-M6",
        "git_commit_at_run": commit, "n_per_suite": args.n,
        "specs_from": str(S154.relative_to(PROJECT_ROOT)),
        "primary_comparator": "base_best_by_mean (script 154)",
        "suites": out,
        "runtime_s": time.time() - t0,
    }
    (out_dir / "wrong_model_results.json").write_text(json.dumps(result, indent=2, default=float))
    for name, rows in out.items():
        r = rows["risk_by_mean"]
        print(f"{name:18s} risk {r['mean']:.1f}  vs primary {r['mean_diff_vs_primary']:+.1f} "
              f"{np.round(r['mean_ci95'], 1)}  vs front {r['mean_diff_vs_frontloaded']:+.1f}")


if __name__ == "__main__":
    main()
