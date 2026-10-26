#!/usr/bin/env python3
"""
Script 156: Track C containment with a weak-drug guard; retune and rerun the wrong-model suites.
================================================================================================
SIMULATION ONLY, assumed parameters. Post hoc follow-up to script 155, where containment
failed under a weak drug (M3, m in [3, 6]: -76.8 days vs the best baseline).

Change: the filter prior on m is widened to [3, 20] and a guard treats every cycle while
P(m < guard) >= 0.5. All arms (containment variants AND baselines) are retuned by mean
capped TTP on one MIXED training set:
  nominal  seed 20261031, n = 200, noise 1_031_001 (script-150 training draw)
  weak     seed 20261071, n = 100, m in [3, 6], noise 1_071_001 (new)
So weak drugs are now anticipated in training: M3 below is no longer an unseen wrong model.
M7 (m in [4, 8], straddling the nominal bound) is added as a partly unseen check.

Evaluation: suites M0-M7 on new seeds 20261080 + k, noise 1_080_0k1 and 1_080_0k2, n = 400.
Primary per suite: guarded containment minus the best retuned baseline (paired bootstrap of the mean).

Run:  python src/156_trackc_guarded_containment.py --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import itertools
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

OUT_DIR = PROJECT_ROOT / "output" / "trackc_guarded"
ONOFF = [(0.6, 0.3), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5), (0.8, 0.7), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7)]
SINGLE = [0.8, 0.9, 1.0, 1.05, 1.1, 1.15]
P_STAR, HORIZONS, GUARDS = [0.02, 0.05, 0.1], [1, 2], [5.0, 6.0, 7.0]


def suite(k, n, base_seed=20261080):
    seed = base_seed + k
    if k == 0:
        return "M0_nominal", tc.draw_scenarios(n, seed)
    if k == 1:
        return "M1_competition", tc.draw_scenarios(n, seed, a_sr=2.0, a_rs=0.5)
    if k == 2:
        return "M2_fast_growth", tc.draw_scenarios(n, seed, ranges=dict(tc.RANGES, r=(0.05, 0.10, "log")))
    if k == 3:
        return "M3_weak_drug", tc.draw_scenarios(n, seed, ranges=dict(tc.RANGES, m=(3.0, 6.0, "log")))
    if k == 7:
        return "M7_drug_4_to_8", tc.draw_scenarios(n, seed, ranges=dict(tc.RANGES, m=(4.0, 8.0, "log")))
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
    ap.add_argument("--n-nominal", type=int, default=200)
    ap.add_argument("--n-weak", type=int, default=100)
    ap.add_argument("--n", type=int, default=400)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    nominal = tc.draw_scenarios(args.n_nominal, 20261031)
    weak = tc.draw_scenarios(args.n_weak, 20261071, ranges=dict(tc.RANGES, m=(3.0, 6.0, "log")))
    plain = [("risk", p, h, 0.0, w) for p, h, w in itertools.product(P_STAR, HORIZONS, (False, True))]
    guarded = [("risk", p, h, g, True) for p, h, g in itertools.product(P_STAR, HORIZONS, GUARDS)]
    base = list(dict.fromkeys([("raw", a, a) for a in SINGLE] + [("raw", a, o) for a, o in ONOFF]
                              + [("periodic", q) for q in range(1, 7)] + [("raw", 1.0, 0.5)]))
    specs = plain + guarded + base
    jobs = [(s, nominal, [1_031_001]) for s in specs] + [(s, weak, [1_071_001]) for s in specs]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = list(ex.map(tc.score_spec, jobs))
    n_s = len(specs)
    tr = {s: np.concatenate([res[i][1], res[n_s + i][1]]) for i, s in enumerate(specs)}
    tr_parts = {s: (float(res[i][1].mean()), float(res[n_s + i][1].mean())) for i, s in enumerate(specs)}
    t_train = time.time() - t0

    mean = lambda s: float(tr[s].mean())  # noqa: E731
    chosen = {
        "guarded": max(guarded, key=mean),
        "unguarded": max(plain, key=mean),
        "base_best": max(base, key=mean),
        "base_frontloaded": ("periodic", 1),
    }

    t1 = time.time()
    suites = [suite(k, args.n) for k in range(8)]
    jobs, keys = [], []
    for k, (name, scen) in enumerate(suites):
        for spec in dict.fromkeys(chosen.values()):
            jobs.append((spec, scen, [1_080_000 + 10 * k + 1, 1_080_000 + 10 * k + 2]))
            keys.append((name, spec))
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        ttp = {key: r[1] for key, r in zip(keys, ex.map(tc.score_spec, jobs))}
    t_eval = time.time() - t1

    out = {}
    for name, _ in suites:
        prim = ttp[(name, chosen["base_best"])]
        front = ttp[(name, chosen["base_frontloaded"])]
        rows = {}
        for a, spec in chosen.items():
            x = ttp[(name, spec)]
            dm, cim = boot_mean_diff(x, prim, seed=51)
            dfr, cif = boot_mean_diff(x, front, seed=52)
            rows[a] = {"spec": list(spec), "mean": float(x.mean()), "cvar10": tc.cvar_lower(x),
                       "mean_diff_vs_base_best": dm, "ci95": cim,
                       "mean_diff_vs_frontloaded": dfr, "frontloaded_ci95": cif}
        out[name] = rows

    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {
        "script": "156_trackc_guarded_containment",
        "label": ("SIMULATION, assumed parameters. POST HOC guard after script 155. Weak drugs are in the "
                  "training mix, so M3 is anticipated; M7 straddles the nominal m bound."),
        "git_commit_at_run": commit,
        "train": {"nominal": {"seed": 20261031, "n": args.n_nominal, "noise": [1_031_001]},
                  "weak": {"seed": 20261071, "n": args.n_weak, "m_range": [3.0, 6.0], "noise": [1_071_001]},
                  "objective": "mean capped TTP over the pooled training scenarios"},
        "train_all": {str(list(s)): {"pooled_mean": mean(s), "nominal_mean": tr_parts[s][0],
                                     "weak_mean": tr_parts[s][1]} for s in specs},
        "chosen": {k: list(v) for k, v in chosen.items()},
        "eval_seeds": {"scenario": "20261080 + k", "noise": "1_080_0k1, 1_080_0k2"}, "n_per_suite": args.n,
        "suites": out,
        "runtime_s": {"train": t_train, "eval": t_eval, "total": time.time() - t0},
    }
    (out_dir / "guarded_results.json").write_text(json.dumps(result, indent=2, default=float))
    print("chosen", chosen)
    for name, rows in out.items():
        g, u = rows["guarded"], rows["unguarded"]
        print(f"{name:16s} guarded {g['mean_diff_vs_base_best']:+.1f} {np.round(g['ci95'], 1)} "
              f"front {g['mean_diff_vs_frontloaded']:+.1f} | unguarded {u['mean_diff_vs_base_best']:+.1f}")


if __name__ == "__main__":
    main()
