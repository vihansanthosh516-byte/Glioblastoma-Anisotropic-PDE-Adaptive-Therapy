#!/usr/bin/env python3
"""
Script 153: Track C headroom map on TRAINING pilot seeds (plan trackAC_new_results_plan_claude.md).
=================================================================================================
SIMULATION ONLY, assumed parameters. Question: how many days could ANY policy gain over the
best simple rule, if it knew the true parameters? This is an upper reference, not a policy.

Per training scenario (same draw as script 150: seed 20261031, n = 200, noise 1_031_001):
  oracle      best schedule with known parameters (trackc_oracle.best_known; an achievable
              lower bound on the true optimum, exact on short-horizon brute-force tests)
  rules       periodic q = 1..6, raw on/off grid, Zhang 1.0/0.5, and the pilot-2 hybrid
              if output/trackc_pilot2/pilot_results.json exists
  best_rule   the one rule with the best training CVaR_0.10 (fixed across scenarios)
  headroom    oracle - best_rule, per scenario

Not yet included: the readings-only bound (belief look-ahead), so headroom here mixes
"value of knowing the parameters" with "value of a smarter policy".
Training seeds only: nothing here touches the evaluation seeds.

Run:  python src/153_trackc_headroom.py --workers 7
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
import trackc_oracle as to  # noqa: E402
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_headroom"
PILOT2 = PROJECT_ROOT / "output" / "trackc_pilot2" / "pilot_results.json"
SEED_TRAIN, NOISE_TRAIN = 20261031, [1_031_001]          # = script 150 training draw
ONOFF = [(0.6, 0.3), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5), (0.8, 0.7), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7)]
MIN_GAIN = 14.0                                          # days; the plan's win margin


def _oracle(args):
    scen, i, G = args
    out = to.best_known(to._scalar(scen, i), G=G)
    return i, out["ttp"], out["ttp_dp"], int(out["schedule"].sum())


def strata(x, headroom, n_bins=4):
    if np.unique(x).size <= n_bins:                      # discrete (sigma): one row per value
        return [{"value": float(v), "n": int((x == v).sum()),
                 "mean_headroom": float(headroom[x == v].mean()),
                 "frac_ge_min_gain": float(np.mean(headroom[x == v] >= MIN_GAIN))} for v in np.unique(x)]
    edges = np.quantile(x, np.linspace(0, 1, n_bins + 1))
    rows = []
    for b in range(n_bins):
        m = (x >= edges[b]) & ((x < edges[b + 1]) if b < n_bins - 1 else (x <= edges[b + 1]))
        rows.append({"lo": float(edges[b]), "hi": float(edges[b + 1]), "n": int(m.sum()),
                     "mean_headroom": float(headroom[m].mean()),
                     "frac_ge_min_gain": float(np.mean(headroom[m] >= MIN_GAIN))})
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-train", type=int, default=200)
    ap.add_argument("--grid", type=int, default=80)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    scen = tc.draw_scenarios(args.n_train, SEED_TRAIN)
    n = args.n_train

    # ---- rules on the same training scenarios and noise as script 150
    specs = [("periodic", q) for q in range(1, 7)] + [("raw", a, o) for a, o in ONOFF] + [("raw", 1.0, 0.5)]
    hybrid_source = None
    if PILOT2.exists():
        chosen = json.loads(PILOT2.read_text())["arms"]["hybrid"]["spec"]
        specs.append(tuple(chosen))
        hybrid_source = str(PILOT2.relative_to(PROJECT_ROOT))
    specs = list(dict.fromkeys(specs))
    jobs = [(s, scen, NOISE_TRAIN) for s in specs]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rule_ttp = dict(ex.map(tc.score_spec, jobs))
    t_rules = time.time() - t0

    # ---- oracle per scenario
    t1 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(_oracle, [(scen, i, args.grid) for i in range(n)]))
    t_oracle = time.time() - t1
    oracle = np.empty(n)
    oracle_dp = np.empty(n)
    oracle_courses = np.empty(n, int)
    for i, t, tdp, c in res:
        oracle[i], oracle_dp[i], oracle_courses[i] = t, tdp, c

    cv = {s: tc.cvar_lower(v) for s, v in rule_ttp.items()}
    best = max(cv, key=cv.get)
    best_ttp = rule_ttp[best]
    per_scen_best = np.max(np.stack(list(rule_ttp.values())), 0)
    headroom = oracle - best_ttp
    d, ci = tc.paired_bootstrap_cvar(oracle, best_ttp, seed=11)

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {
        "script": "153_trackc_headroom",
        "label": "SIMULATION, training pilot seeds only; oracle knows true parameters (upper reference)",
        "git_commit_at_run": commit,
        "seeds": {"train": SEED_TRAIN, "noise_train": NOISE_TRAIN},
        "n": n, "grid": args.grid,
        "hybrid_spec_source": hybrid_source,
        "rules": {str(list(s)): {"cvar10": cv[s], "mean": float(v.mean())} for s, v in rule_ttp.items()},
        "best_rule_by_train_cvar": list(best),
        "oracle": {"cvar10": tc.cvar_lower(oracle), "mean": float(oracle.mean()),
                   "capped_frac": float(np.mean(oracle >= tc.HORIZON - 1e-9)),
                   "mean_courses": float(oracle_courses.mean()),
                   "dp_before_polish_mean": float(oracle_dp.mean())},
        "oracle_minus_best_rule": {
            "cvar10_diff": d, "ci95": ci,
            "mean_headroom": float(headroom.mean()), "median_headroom": float(np.median(headroom)),
            "frac_ge_min_gain": float(np.mean(headroom >= MIN_GAIN)),
            "frac_negative": float(np.mean(headroom < -1e-9)),
            "min_gain_days": MIN_GAIN},
        "oracle_minus_per_scenario_best_rule_mean": float((oracle - per_scen_best).mean()),
        "strata": {k: strata(scen[k], headroom) for k in ("r", "m", "f0", "c", "eps", "B0", "sigma")},
        "not_included": ["readings-only belief look-ahead bound"],
        "runtime_s": {"rules": t_rules, "oracle": t_oracle, "total": time.time() - t0},
    }
    (out_dir / "headroom_results.json").write_text(json.dumps(result, indent=2, default=float))
    cols = list(tc.PARAMS) + ["sigma"]
    rows = ["scenario," + ",".join(cols) + ",oracle,best_rule,headroom,oracle_courses"]
    rows += [f"{i}," + ",".join(f"{scen[c][i]:.6g}" for c in cols)
             + f",{oracle[i]:.2f},{best_ttp[i]:.2f},{headroom[i]:.2f},{oracle_courses[i]}" for i in range(n)]
    (out_dir / "headroom_per_scenario.csv").write_text("\n".join(rows) + "\n")
    print(json.dumps({k: result[k] for k in ("best_rule_by_train_cvar", "oracle", "oracle_minus_best_rule", "runtime_s")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
