#!/usr/bin/env python3
"""
Script 154: Track C belief containment vs the best rules, on NEW evaluation seeds.
================================================================================
SIMULATION ONLY, assumed parameters. Exploratory, post hoc design:
- Pilot 2 (script 150) showed the hybrid losing to front-loaded dosing in CVaR_0.10.
- The headroom map (script 153, training seeds) showed the known-parameter oracle gains
  little in CVaR_0.10 (+3.4 d) but much in the mean (+83.5 d), and mostly doses near
  1.1-1.19 x B0 ("contain just below the progression line").
So this script makes MEAN capped TTP the primary endpoint (chosen after seeing pilot data;
label it post hoc) and tests a belief-containment rule built from that oracle pattern.

Arms, all tuned on the script-150 training draw (seed 20261031, n = 200, noise 1_031_001):
  risk        RiskContainment: treat when P(progress within h untreated cycles) >= p*
  raw_single  single raw threshold m >= a, a in 0.8..1.15 (gives the raw rule the same
              near-the-line freedom; this is the fair no-filter comparator)
  raw_onoff   script-150 raw on/off grid; Zhang 1.0/0.5; periodic q = 1..6;
  hybrid      the pilot-2 hybrid spec
Primary comparator: the non-filter baseline with the best TRAINING mean.
Evaluation: new seeds 20261042 (n = 400), noise 1_042_001 and 1_042_002, never used before.
Also: known-parameter oracle on the evaluation scenarios, to report the share of headroom
each arm captures.

Run:  python src/154_trackc_containment.py --workers 7
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

OUT_DIR = PROJECT_ROOT / "output" / "trackc_containment"
PILOT2 = PROJECT_ROOT / "output" / "trackc_pilot2" / "pilot_results.json"
SEED_TRAIN, NOISE_TRAIN = 20261031, [1_031_001]
SEED_EVAL, NOISE_EVAL = 20261042, [1_042_001, 1_042_002]
ONOFF = [(0.6, 0.3), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5), (0.8, 0.7), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7)]
P_STAR = [0.02, 0.05, 0.1, 0.2, 0.35, 0.5]
HORIZONS = [1, 2]
SINGLE = [0.8, 0.9, 1.0, 1.05, 1.1, 1.15]
MIN_GAIN = 14.0


def _oracle(args):
    scen, i = args
    return i, to.best_known(to._scalar(scen, i))["ttp"]


def boot_mean_diff(a, b, n_boot=2000, seed=0):
    rng = np.random.default_rng(seed)
    d = a - b
    idx = rng.integers(0, d.size, (n_boot, d.size))
    return float(d.mean()), np.percentile(d[idx].mean(1), [2.5, 97.5]).tolist()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-train", type=int, default=200)
    ap.add_argument("--n-eval", type=int, default=400)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    train = tc.draw_scenarios(args.n_train, SEED_TRAIN)
    ev = tc.draw_scenarios(args.n_eval, SEED_EVAL)

    risk_specs = [("risk", p, h) for p in P_STAR for h in HORIZONS]
    single_specs = [("raw", a, a) for a in SINGLE]
    onoff_specs = [("raw", a, o) for a, o in ONOFF]
    per_specs = [("periodic", q) for q in range(1, 7)]
    zhang = ("raw", 1.0, 0.5)
    hybrid = tuple(json.loads(PILOT2.read_text())["arms"]["hybrid"]["spec"])
    base_specs = list(dict.fromkeys(single_specs + onoff_specs + per_specs + [zhang]))
    specs = list(dict.fromkeys(risk_specs + base_specs + [hybrid]))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        tr = dict(ex.map(tc.score_spec, [(s, train, NOISE_TRAIN) for s in specs]))
    t_train = time.time() - t0

    mean = lambda s: float(tr[s].mean())  # noqa: E731
    chosen = {
        "risk_by_mean": max(risk_specs, key=mean),
        "risk_by_cvar": max(risk_specs, key=lambda s: tc.cvar_lower(tr[s])),
        "base_best_by_mean": max(base_specs, key=mean),
        "base_raw_single": max(single_specs, key=mean),
        "base_raw_onoff": max(onoff_specs, key=mean),
        "base_best_periodic": max(per_specs, key=mean),
        "base_frontloaded": ("periodic", 1),
        "base_zhang": zhang,
        "pilot2_hybrid": hybrid,
    }
    t1 = time.time()
    eval_specs = list(dict.fromkeys(chosen.values()))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        es = dict(ex.map(tc.score_spec, [(s, ev, NOISE_EVAL) for s in eval_specs]))
        orc = np.empty(args.n_eval)
        for i, t in ex.map(_oracle, [(ev, i) for i in range(args.n_eval)]):
            orc[i] = t
    t_eval = time.time() - t1

    prim = es[chosen["base_best_by_mean"]]
    gap = orc - prim
    arms = {}
    for name, spec in chosen.items():
        x = es[spec]
        dm, cim = boot_mean_diff(x, prim, seed=21)
        dc, cic = tc.paired_bootstrap_cvar(x, prim, seed=22)
        arms[name] = {"spec": list(spec), "train_mean": mean(spec), "train_cvar10": tc.cvar_lower(tr[spec]),
                      "eval_mean": float(x.mean()), "eval_cvar10": tc.cvar_lower(x),
                      "eval_median": float(np.median(x)),
                      "mean_diff_vs_primary": dm, "mean_ci95": cim,
                      "cvar10_diff_vs_primary": dc, "cvar10_ci95": cic,
                      "share_of_oracle_gap_captured": float((x - prim).mean() / gap.mean()) if gap.mean() > 0 else None}

    strata = {}
    for k in ("r", "B0", "sigma", "f0"):
        q = np.unique(ev[k]) if k == "sigma" else np.quantile(ev[k], [0, 0.25, 0.5, 0.75, 1])
        groups = ([ev[k] == v for v in q] if k == "sigma"
                  else [(ev[k] >= q[b]) & ((ev[k] < q[b + 1]) if b < 3 else (ev[k] <= q[b + 1])) for b in range(4)])
        rows = []
        for gi, m in enumerate(groups):
            d, ci = boot_mean_diff(es[chosen["risk_by_mean"]][m], prim[m], seed=30 + gi)
            rows.append({"group": float(q[gi]), "n": int(m.sum()), "mean_diff": d, "ci95": ci,
                         "oracle_gap": float(gap[m].mean())})
        strata[k] = rows

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {
        "script": "154_trackc_containment",
        "label": ("SIMULATION, assumed parameters. POST HOC: mean endpoint and the containment arm were "
                  "chosen after pilot 2 and the training headroom map. Evaluation seeds are new."),
        "git_commit_at_run": commit,
        "seeds": {"train": SEED_TRAIN, "noise_train": NOISE_TRAIN, "eval": SEED_EVAL, "noise_eval": NOISE_EVAL},
        "n": {"train": args.n_train, "eval": args.n_eval},
        "primary_endpoint": "mean capped TTP (days)", "win_rule": {"mean_ci_above_0": True, "min_gain_days": MIN_GAIN},
        "primary_comparator": "base_best_by_mean",
        "train_all": {str(list(s)): {"mean": mean(s), "cvar10": tc.cvar_lower(tr[s])} for s in specs},
        "arms": arms,
        "oracle_eval": {"mean": float(orc.mean()), "cvar10": tc.cvar_lower(orc),
                        "mean_gap_over_primary": float(gap.mean())},
        "strata_risk_by_mean_vs_primary": strata,
        "runtime_s": {"train": t_train, "eval_and_oracle": t_eval, "total": time.time() - t0},
    }
    (out_dir / "containment_results.json").write_text(json.dumps(result, indent=2, default=float))
    names = list(chosen)
    rows = ["scenario," + ",".join(names) + ",oracle"]
    rows += [f"{i}," + ",".join(f"{es[chosen[k]][i]:.2f}" for k in names) + f",{orc[i]:.2f}" for i in range(args.n_eval)]
    (out_dir / "eval_ttp_per_scenario.csv").write_text("\n".join(rows) + "\n")
    print(json.dumps({k: {a: v[a] for a in ("spec", "eval_mean", "mean_diff_vs_primary", "mean_ci95",
                                                "eval_cvar10", "share_of_oracle_gap_captured")}
                      for k, v in arms.items()}, indent=1, default=float))


if __name__ == "__main__":
    main()
