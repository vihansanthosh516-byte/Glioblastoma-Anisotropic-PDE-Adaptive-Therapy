#!/usr/bin/env python3
"""
Script 176: Track C DEVELOPMENT-ONLY exploration of reversible drug-induced tolerance.
=====================================================================================
SIMULATION ONLY; exploratory; development seeds only (never a test). Mechanism (assumed rates):
drug-induced switching S -> R while on drug (mu, existing script-169 term) plus reversion
R -> S while off drug (nu, new). Qualitative support only: transient TMZ-tolerant states
(Rabe et al. 2020 Cell Death Dis, doi 10.1038/s41419-019-2200-2) and reversible GSC persister
states (Liau et al. 2017 Cell Stem Cell, doi 10.1016/j.stem.2016.11.003); no source gives
day-scale rates, so all mu/nu ranges are ASSUMED.

Question: does this mechanism give a break-based containment rule ("riskh": forced holiday
after max_run treated cycles) room to beat front-loaded dosing in the worst 10%, and which
simple baseline is strongest? Output feeds the choice of the method and comparators that a
later sealed script freezes before any test seed is opened.

Dev seeds: scenario 20262100 + k; noise [1_300_000 + 10k + rep, chunk]. Unused before.

Run:  python src/176_trackc_dev_plasticity.py --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_dev_plasticity"
N, CHUNK, N_REP = 150, 75, 2
PLASTIC = {"nominal": None,                                   # mechanism absent (control)
           "plastic_mod": ((0.005, 0.02), (0.005, 0.02)),     # (mu range, nu range) per day, log-uniform
           "plastic_strong": ((0.02, 0.05), (0.02, 0.05))}
SPECS = [("periodic", 1, 0), ("periodic", 2, 0), ("periodic", 3, 0), ("periodic", 4, 0),
         ("raw", 0.6, 0.6), ("riskf", 0.1, 2, 0), ("riskf", 0.05, 2, 4),
         ("riskh", 0.1, 2, 0, 2), ("riskh", 0.1, 2, 0, 3), ("riskh", 0.05, 2, 0, 2),
         ("riskh", 0.05, 2, 2, 2)]


def build(k, name, n):
    s = tc.draw_scenarios(n, 20262100 + k)
    if PLASTIC[name]:
        rng = np.random.default_rng([20262100 + k, 3])
        (m0, m1), (n0, n1) = PLASTIC[name]
        s["mu"] = np.exp(rng.uniform(np.log(m0), np.log(m1), n))
        s["nu"] = np.exp(rng.uniform(np.log(n0), np.log(n1), n))
    return s


def _job(args):
    spec, k, name, c = args
    sub = tc.subset(build(k, name, N), np.arange(c * CHUNK, min(N, (c + 1) * CHUNK)))
    pol = tc.make_policy(tuple(spec))
    reps = np.stack([tc.simulate(sub, pol, [1_300_000 + 10 * k + r, c])["ttp"] for r in range(1, N_REP + 1)])
    return (name, tuple(spec), c), reps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fp = OUT_DIR / "dev_results.json"
    if fp.exists():
        sys.exit(f"refusing to overwrite {fp}")
    t0 = time.time()
    n_chunks = -(-N // CHUNK)
    jobs = [(s, k, nm, c) for k, nm in enumerate(PLASTIC) for c in range(n_chunks) for s in SPECS]
    jobs.sort(key=lambda j: {"riskf": 0, "riskh": 0}.get(j[0][0], 1))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = dict(ex.map(_job, jobs))
    out = {}
    for nm in PLASTIC:
        out[nm] = {}
        for s in SPECS:
            r = np.concatenate([parts[(nm, s, c)] for c in range(n_chunks)], 1)
            out[nm][str(list(s))] = {"mean": float(r.mean()), "cvar10_scenario_mean": tc.cvar_scenario_mean(r)}
    fp.write_text(json.dumps({"script": "176_trackc_dev_plasticity",
                              "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; mu/nu ASSUMED",
                              "plastic": PLASTIC, "n": N, "n_rep": N_REP, "specs": SPECS,
                              "seed_rule": "scenario 20262100 + k; noise [1_300_000 + 10k + rep, chunk]",
                              "results": out, "runtime_s": time.time() - t0}, indent=2, default=float))
    for nm in PLASTIC:
        print("==", nm)
        for s in SPECS:
            v = out[nm][str(list(s))]
            print(f"  {str(s):32s} mean {v['mean']:7.1f}  cvar10 {v['cvar10_scenario_mean']:6.1f}")


if __name__ == "__main__":
    main()
