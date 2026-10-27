#!/usr/bin/env python3
"""
Script 159: Track C numerical check of progression times and policy rankings (audit C1).
=======================================================================================
SIMULATION ONLY. Numerics, not policy evaluation: equations unchanged, fresh numerics-only
seeds (20261101-20261103) never used for tuning or scoring.

Part A, fixed schedules (none, periodic q = 1, 2, 4): time to 1.2 x B0 from
  - trackc_sim RK4, dt = 0.5 / 0.25 / 0.1 / 0.05, EVENT_MODE "refine"
  - trackc_sim RK4, dt = 0.25, EVENT_MODE "grid" (the rule used by scripts 150-158)
  - SciPy solve_ivp DOP853 (rtol 1e-10) with an event function (trackc_reference)
  reported as absolute differences to the SciPy reference.
Part B, policy rankings: the cheap rule families (periodic, raw on/off, raw single, Zhang)
  scored at each dt / event mode on the same scenarios and noise; ranking agreement with the
  finest setting (Kendall tau, top-1, max change in mean and CVaR_0.10).
Part C, two filter policies (containment narrow and wide) at dt = 0.25 refine / grid / 0.1.

Run:  python src/159_trackc_numerics_check.py --workers 7
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
from scipy.stats import kendalltau

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_oracle as to  # noqa: E402
import trackc_reference as ref  # noqa: E402
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_numerics_check"
SEED_A, SEED_B, SEED_C = 20261101, 20261102, 20261103
NOISE = 1_101_001
SETTINGS = [("dt0.5", 0.5, "refine"), ("dt0.25_grid", 0.25, "grid"), ("dt0.25", 0.25, "refine"),
            ("dt0.1", 0.1, "refine"), ("dt0.05", 0.05, "refine")]
ONOFF = [(0.6, 0.3), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5), (0.8, 0.7), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7)]
CHEAP = ([("periodic", q) for q in range(1, 7)] + [("raw", a, o) for a, o in ONOFF]
         + [("raw", a, a) for a in (0.8, 0.9, 1.0, 1.05, 1.1, 1.15)])
FILTER = [("risk", 0.05, 2), ("risk", 0.02, 1, 0.0, True)]


def _set(dt, mode):
    tc.DT, tc.EVENT_MODE = dt, mode


def _part_a(args):
    scen, i = args
    p = to._scalar(scen, i)
    sched = np.vstack([np.zeros(tc.N_CYCLES, bool), to.periodic_schedules(qs=[1, 2, 4])])
    row = {"ref": [ref.reference_ttp(p, s) for s in sched]}
    for name, dt, mode in SETTINGS:
        _set(dt, mode)
        row[name] = to.rollout(p, sched).tolist()
    _set(0.25, "refine")
    return i, row


def _score(args):
    spec, scen, setting = args
    _, dt, mode = setting
    _set(dt, mode)
    out = tc.simulate(scen, tc.make_policy(spec), NOISE)["ttp"]
    _set(0.25, "refine")
    return spec, setting[0], out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-a", type=int, default=150)
    ap.add_argument("--n-b", type=int, default=200)
    ap.add_argument("--n-c", type=int, default=60)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # ---- A: event times vs SciPy
    sa = tc.draw_scenarios(args.n_a, SEED_A)
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = dict(ex.map(_part_a, [(sa, i) for i in range(args.n_a)]))
    refm = np.array([rows[i]["ref"] for i in range(args.n_a)])
    part_a = {}
    for name, _, _ in SETTINGS:
        m = np.array([rows[i][name] for i in range(args.n_a)])
        d = np.abs(m - refm)
        part_a[name] = {"max_abs_days": float(d.max()), "p99_abs_days": float(np.quantile(d, 0.99)),
                        "median_abs_days": float(np.median(d)), "mean_signed_days": float((m - refm).mean())}
    part_a["n_event_times"] = int(refm.size)
    part_a["frac_censored_ref"] = float(np.mean(refm >= tc.HORIZON - 1e-9))

    # ---- B: cheap policy rankings
    sb = tc.draw_scenarios(args.n_b, SEED_B)
    jobs = [(s, sb, st) for s in CHEAP for st in SETTINGS]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(_score, jobs))
    tb = {(s, n): v for s, n, v in res}
    part_b = {}
    finest = "dt0.05"
    for stat_name, stat in (("mean", np.mean), ("cvar10", tc.cvar_lower)):
        f = np.array([stat(tb[(s, finest)]) for s in CHEAP])
        rows_b = {}
        for name, _, _ in SETTINGS:
            v = np.array([stat(tb[(s, name)]) for s in CHEAP])
            tau = kendalltau(v, f).statistic
            rows_b[name] = {"kendall_tau_vs_dt0.05": float(tau),
                            "top1_same": bool(np.argmax(v) == np.argmax(f)),
                            "max_abs_change_days": float(np.max(np.abs(v - f)))}
        part_b[stat_name] = rows_b
    part_b["max_per_scenario_abs_change_dt0.25_vs_dt0.05_days"] = float(
        max(np.max(np.abs(tb[(s, "dt0.25")] - tb[(s, finest)])) for s in CHEAP))
    part_b["frac_scenarios_changed_gt_1day_dt0.25_vs_dt0.05"] = float(
        np.mean([np.mean(np.abs(tb[(s, "dt0.25")] - tb[(s, finest)]) > 1.0) for s in CHEAP]))

    # ---- C: filter policies
    sc = tc.draw_scenarios(args.n_c, SEED_C)
    c_set = [st for st in SETTINGS if st[0] in ("dt0.25_grid", "dt0.25", "dt0.1")]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        resc = list(ex.map(_score, [(s, sc, st) for s in FILTER for st in c_set]))
    tcc = {(s, n): v for s, n, v in resc}
    part_c = {}
    for s in FILTER:
        base = tcc[(s, "dt0.1")]
        part_c[str(list(s))] = {n: {"mean": float(tcc[(s, n)].mean()),
                                    "max_abs_change_vs_dt0.1_days": float(np.max(np.abs(tcc[(s, n)] - base))),
                                    "frac_scenarios_changed_gt_1day": float(np.mean(np.abs(tcc[(s, n)] - base) > 1))}
                                for n, _, _ in c_set}

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {"script": "159_trackc_numerics_check",
              "label": "SIMULATION numerics check; equations unchanged; fresh numerics-only seeds",
              "git_commit_at_run": commit,
              "seeds": {"part_a": SEED_A, "part_b": SEED_B, "part_c": SEED_C, "noise": NOISE},
              "n": {"part_a_scenarios": args.n_a, "part_b_scenarios": args.n_b, "part_c_scenarios": args.n_c},
              "reference": "scipy.integrate.solve_ivp DOP853 rtol 1e-10 atol 1e-13, piecewise on/off, event at 1.2 B0",
              "part_a_event_times_vs_reference": part_a,
              "part_b_cheap_policy_rankings": part_b,
              "part_c_filter_policies": part_c,
              "runtime_s": time.time() - t0}
    (out_dir / "numerics_results.json").write_text(json.dumps(result, indent=2, default=float))
    print(json.dumps({"A": part_a, "B": part_b["mean"], "C": part_c}, indent=1, default=float))


if __name__ == "__main__":
    main()
