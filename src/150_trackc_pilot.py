#!/usr/bin/env python3
"""
Script 150: Track C PILOT arena (Claude final plan, appendix Track C; exploratory).
===================================================================================
SIMULATION ONLY, assumed parameters, PILOT seeds. This is not the frozen evaluation:
the final plan must draw new training, validation and held-out seeds. Nothing here is
evidence about patients or drugs.

Question: under one shared 12-course budget, does the hybrid (particle-filter burden
estimate + on/off thresholds + resistance mode, thresholds chosen for the worst 10% of
scenarios) beat the strongest existing-style rule chosen on training scenarios?

Arms
  hybrid      FilteredThreshold, 128 configs (a_on, a_off) x a_R x b, chosen by training CVaR_0.10
  ablations   no resistance mode (a_R = a_on), mean objective (same grid, chosen by mean),
              no filter (hybrid's a_on/a_off on raw readings), no hysteresis (single threshold)
  baselines   periodic q = 1..6 (q = 1 front-loaded, q = 4 paced), Zhang 1.0/0.5 on raw
              readings, raw on/off rule re-tuned by training CVaR
  not in the pilot: PPO, belief look-ahead planner, dynamic-programming oracle
Primary comparator: the baseline with the best training CVaR_0.10.
Endpoint: first day true burden >= 1.2 B0, capped at 1,456 days; CVaR_0.10 = mean of the
lowest 10% of per-scenario TTP (noise replicates averaged first).

Run:  python src/150_trackc_pilot.py --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")   # one BLAS thread per worker process

import argparse
import itertools
import json
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_pilot2"

# PILOT seeds only. The frozen evaluation must use different, newly recorded seeds.
# Pilot 2 (2026-10-24): new seeds after the kill-ratio and h fixes.
SEED_TRAIN, SEED_EVAL, SEED_POS = 20261031, 20261032, 20261033
NOISE_TRAIN = [1_031_001]
NOISE_EVAL = [1_032_001, 1_032_002]
NOISE_POS = [1_033_001]

ONOFF = [(0.6, 0.3), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5), (0.8, 0.7), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7)]
A_R = [0.6, 0.8, 1.0, 1.1]
B_TRIG = [0.05, 0.10, 0.20, 0.35]
POS_RANGES = dict(tc.RANGES, B0=(0.6, 0.8, "lin"), f0=(0.001, 0.01, "log"), c=(0.1, 0.25, "lin"))


def score_all(specs, scen, seeds, workers):
    jobs = [(s, scen, seeds) for s in specs]
    if workers <= 1:
        return dict(map(tc.score_spec, jobs))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        return dict(ex.map(tc.score_spec, jobs))


def best(scores, key):
    return max(scores, key=lambda s: key(scores[s]))


def summary(ttp):
    return {"cvar10": tc.cvar_lower(ttp), "mean": float(ttp.mean()), "median": float(np.median(ttp)),
            "capped_frac_all_replicates": float(np.mean(ttp >= tc.HORIZON - 1e-9))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-train", type=int, default=200)
    ap.add_argument("--n-eval", type=int, default=400)
    ap.add_argument("--n-pos", type=int, default=200)
    args = ap.parse_args()
    t_start = time.time()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    train = tc.draw_scenarios(args.n_train, SEED_TRAIN)
    ev = tc.draw_scenarios(args.n_eval, SEED_EVAL)

    # ---- stage 1: choose every tuned arm on training scenarios only
    hyb_specs = [("filt", a, o, r, b) for (a, o), r, b in itertools.product(ONOFF, A_R, B_TRIG)]
    single_specs = [("filt", a, a, r, b, True) for a, r, b in itertools.product([0.6, 0.8, 1.0], A_R, B_TRIG)]
    nores_specs = [("filt", a, o, a, 0.2) for a, o in ONOFF]
    raw_specs = [("raw", a, o) for a, o in ONOFF]
    per_specs = [("periodic", q) for q in range(1, 7)]
    zhang = ("raw", 1.0, 0.5)
    train_specs = list(dict.fromkeys(hyb_specs + single_specs + nores_specs + raw_specs + per_specs + [zhang]))
    t0 = time.time()
    tr = score_all(train_specs, train, NOISE_TRAIN, args.workers)
    t_train = time.time() - t0
    cv = lambda s: tc.cvar_lower(tr[s])  # noqa: E731

    chosen = {
        "hybrid": best({s: tr[s] for s in hyb_specs}, tc.cvar_lower),
        "abl_mean_objective": best({s: tr[s] for s in hyb_specs}, np.mean),
        "abl_no_resistance_mode": best({s: tr[s] for s in nores_specs}, tc.cvar_lower),
        "abl_no_hysteresis": best({s: tr[s] for s in single_specs}, tc.cvar_lower),
        "base_raw_tuned": best({s: tr[s] for s in raw_specs}, tc.cvar_lower),
        "base_best_periodic": best({s: tr[s] for s in per_specs}, tc.cvar_lower),
        "base_fixed_frontloaded": ("periodic", 1),
        "base_paced": ("periodic", 4),
        "base_zhang": zhang,
    }
    h = chosen["hybrid"]
    chosen["abl_no_filter"] = ("raw", h[1], h[2])
    baselines = [k for k in chosen if k.startswith("base_")]
    primary = max(baselines, key=lambda k: cv(chosen[k]) if chosen[k] in tr else -1)

    # ---- stage 2: score the chosen arms once on held-out pilot scenarios
    eval_specs = list(dict.fromkeys(chosen.values()))
    t0 = time.time()
    ev_scores = score_all(eval_specs, ev, NOISE_EVAL, args.workers)
    t_eval = time.time() - t0

    arms = {}
    for name, spec in chosen.items():
        ttp = ev_scores[spec]
        arms[name] = {"spec": list(spec), "train_cvar10": cv(spec) if spec in tr else None,
                      "eval": summary(ttp)}
    hyb = ev_scores[h]
    comparisons = {}
    for name, spec in chosen.items():
        if name == "hybrid":
            continue
        d, ci = tc.paired_bootstrap_cvar(hyb, ev_scores[spec], seed=7)
        comparisons[name] = {"cvar10_diff_hybrid_minus_arm": d, "ci95": ci}

    strata = {}
    f0_lo, c_lo = ev["f0"] < 0.01, ev["c"] < 0.125
    for lab, m in {"f0<0.01,c<0.125": f0_lo & c_lo, "f0<0.01,c>=0.125": f0_lo & ~c_lo,
                   "f0>=0.01,c<0.125": ~f0_lo & c_lo, "f0>=0.01,c>=0.125": ~f0_lo & ~c_lo}.items():
        d, ci = tc.paired_bootstrap_cvar(hyb[m], ev_scores[chosen[primary]][m], seed=8)
        strata[lab] = {"n": int(m.sum()), "cvar10_diff": d, "ci95": ci}

    # ---- filter gate (pilot: on training scenarios, selected hybrid)
    pol = tc.make_policy(h)
    cov = tc.simulate(train, pol, NOISE_TRAIN[0], track_fr=True)["fr_cover"]

    # ---- positive control: regime where theory expects adaptive dosing to win
    pos = tc.draw_scenarios(args.n_pos, SEED_POS, ranges=POS_RANGES, a_sr=2.0)
    pos_arms = {"paced": ("periodic", 4), "zhang": zhang, "raw_tuned": chosen["base_raw_tuned"]}
    pos_scores = {k: tc.evaluate(pos, tc.make_policy(s), NOISE_POS)["ttp"] for k, s in pos_arms.items()}
    hp = tc.FilteredThreshold(*h[1:], ranges=POS_RANGES)
    pos_scores["hybrid_prior_matched_B0"] = tc.evaluate(pos, hp, NOISE_POS)["ttp"]
    pos_ctrl = {}
    for k, v in pos_scores.items():
        entry = summary(v)
        if k != "paced":
            d, ci = tc.paired_bootstrap_cvar(v, pos_scores["paced"], seed=9)
            entry.update(cvar10_diff_vs_paced=d, ci95=ci)
        pos_ctrl[k] = entry

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {
        "script": "150_trackc_pilot",
        "label": "PILOT, exploratory simulation with assumed parameters; not the frozen evaluation",
        "git_commit_at_run": commit, "python": platform.python_version(), "numpy": np.__version__,
        "seeds": {"train": SEED_TRAIN, "eval": SEED_EVAL, "pos": SEED_POS,
                  "noise_train": NOISE_TRAIN, "noise_eval": NOISE_EVAL, "noise_pos": NOISE_POS},
        "n": {"train": args.n_train, "eval": args.n_eval, "pos": args.n_pos, "particles": 500},
        "settings": {"dt_days": tc.DT, "budget_courses": tc.BUDGET, "horizon_days": tc.HORIZON,
                     "progression": "true burden >= 1.2 x B0", "ranges": tc.RANGES, "pos_ranges": POS_RANGES,
                     "pos_competition": {"a_sr": 2.0, "a_rs": 1.0}},
        "n_configs_scored_on_train": len(train_specs),
        "primary_comparator": primary,
        "arms": arms,
        "hybrid_vs_arms_eval": comparisons,
        "strata_hybrid_vs_primary": strata,
        "filter_fr_90pct_interval_coverage_train": {"hits": cov[0], "total": cov[1],
                                                    "coverage": cov[0] / max(cov[1], 1)},
        "positive_control": pos_ctrl,
        "not_run": ["PPO", "belief look-ahead planner", "dynamic-programming oracle", "misspecified suites M1-M5"],
        "runtime_s": {"train": t_train, "eval": t_eval, "total": time.time() - t_start},
    }
    (OUT_DIR / "pilot_results.json").write_text(json.dumps(result, indent=2, default=float))
    rows = ["scenario," + ",".join(chosen)]
    rows += [f"{i}," + ",".join(f"{ev_scores[chosen[k]][i]:.2f}" for k in chosen) for i in range(args.n_eval)]
    (OUT_DIR / "eval_ttp_per_scenario.csv").write_text("\n".join(rows) + "\n")
    print(json.dumps({k: result[k] for k in ("primary_comparator", "hybrid_vs_arms_eval", "runtime_s")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
