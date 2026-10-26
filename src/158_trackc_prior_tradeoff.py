#!/usr/bin/env python3
"""
Script 158: Track C narrow- vs wide-prior containment, head to head on the SAME new scenarios.
=============================================================================================
SIMULATION ONLY. Scripts 155 and 156 used different seeds, so their gains cannot be compared
directly. This script fixes the specs chosen earlier and scores them on identical scenarios:
  narrow  ("risk", 0.05, 2)               script 154 choice (prior m in [6, 20])
  wide    ("risk", 0.02, 1, 0.0, True)    script 156 unguarded choice (prior m in [3, 20])
  base    ("raw", 0.6, 0.5)               best tuned simple rule in scripts 154 and 156
  front   ("periodic", 1)
Suites M0-M7 as in script 156, new seeds 20261090 + k, noise 1_090_0k1 and 1_090_0k2, n = 400.
Primary: wide minus narrow per suite (paired bootstrap of the mean).

Run:  python src/158_trackc_prior_tradeoff.py --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import importlib.util
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

_spec = importlib.util.spec_from_file_location("s156", PROJECT_ROOT / "src" / "156_trackc_guarded_containment.py")
s156 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s156)

OUT_DIR = PROJECT_ROOT / "output" / "trackc_prior_tradeoff"
ARMS = {"narrow": ("risk", 0.05, 2), "wide": ("risk", 0.02, 1, 0.0, True),
        "base": ("raw", 0.6, 0.5), "front": ("periodic", 1)}
PAIRS = [("wide", "narrow"), ("narrow", "base"), ("wide", "base"), ("narrow", "front"), ("wide", "front")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    suites = [s156.suite(k, args.n, base_seed=20261090) for k in range(8)]
    jobs, keys = [], []
    for k, (name, scen) in enumerate(suites):
        for a, spec in ARMS.items():
            jobs.append((spec, scen, [1_090_000 + 10 * k + 1, 1_090_000 + 10 * k + 2]))
            keys.append((name, a))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        ttp = {key: r[1] for key, r in zip(keys, ex.map(tc.score_spec, jobs))}

    out = {}
    for name, _ in suites:
        row = {a: {"mean": float(ttp[(name, a)].mean()), "cvar10": tc.cvar_lower(ttp[(name, a)])} for a in ARMS}
        for a, b in PAIRS:
            d, ci = s156.boot_mean_diff(ttp[(name, a)], ttp[(name, b)], seed=61)
            row[f"{a}_minus_{b}"] = {"mean": d, "ci95": ci}
        out[name] = row
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {"script": "158_trackc_prior_tradeoff",
              "label": "SIMULATION, assumed parameters; fixed specs from scripts 154 and 156; same scenarios for all arms",
              "git_commit_at_run": commit, "arms": {a: list(s) for a, s in ARMS.items()},
              "seeds": {"scenario": "20261090 + k", "noise": "1_090_0k1, 1_090_0k2"}, "n_per_suite": args.n,
              "suites": out, "runtime_s": time.time() - t0}
    (out_dir / "tradeoff_results.json").write_text(json.dumps(result, indent=2, default=float))
    for name, row in out.items():
        print(name, {k: round(v["mean"], 1) for k, v in row.items() if "_minus_" in k})


if __name__ == "__main__":
    main()
