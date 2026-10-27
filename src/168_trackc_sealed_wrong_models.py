#!/usr/bin/env python3
"""
Script 168: Track C sealed wrong-model grid for the frozen script-167 arms.
=========================================================================
SIMULATION ONLY. Arms are read from output/trackc_sealed_containment/frozen_protocol.json
(no re-tuning). The filter keeps its own model and prior; only the truth changes.

PROTOCOL (fixed in code and committed before any result; SHA-256 written to the output)
  Arms: primary_containment (narrow), winner_containment_wide, primary_nonfilter (raw single),
        winner_periodic (front-loaded), reference Zhang.
  Severity grid (moderate, severe), one axis at a time, plus 2 combinations and a nominal
  reference (17 suites):
    weak drug        m in [4, 6] / [3, 4]                  (prior [6, 20])
    fast growth      r in [0.05, 0.10] / [0.08, 0.15]      (prior [0.001, 0.05])
    noisy readings   true sd x 2 / x 3                     (filter told the nominal sd)
    reading bias     +10% / +20% on every reading
    heavy tails      Student-t df 5 / df 3, same sd
    competition      (a_sr, a_rs) = (2, 0.5) / (3, 0.33)   (filter assumes 1, 1)
    high resistance  f0 in [0.1, 0.3] / [0.2, 0.5]          (prior [0.001, 0.1])
    combo A          reading bias +10% and noise x 2
    combo B          r in [0.05, 0.10] and m in [4, 6]
  Seeds: suite k scenario seed 20261200 + k, noise 1_200_0k1 and 1_200_0k2; n = 150.
  Estimand per suite: mean TTP difference vs primary_nonfilter, two-stage paired bootstrap
    (2,000 draws). Also both CVaR_0.10 estimands (descriptive).
  Pre-set robustness reading, per containment arm:
    "harm"    95% CI upper bound < 0 in a suite
    "benefit" 95% CI lower bound > 0
    "no harm across the grid" = no suite with harm. Every suite is reported.
  Scored once: refuses to overwrite the sealed result.

Run:  python src/168_trackc_sealed_wrong_models.py --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import hashlib
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

FROZEN = PROJECT_ROOT / "output" / "trackc_sealed_containment" / "frozen_protocol.json"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_sealed_wrong_models"
N, N_BOOT = 150, 2000
R = tc.RANGES
SUITES = {  # name: (range overrides, a_sr, a_rs, sd multiplier, bias fraction, student-t df)
    "nominal":            ({}, 1.0, 1.0, 1.0, 0.0, None),
    "weak_drug_mod":      ({"m": (4.0, 6.0, "log")}, 1.0, 1.0, 1.0, 0.0, None),
    "weak_drug_sev":      ({"m": (3.0, 4.0, "log")}, 1.0, 1.0, 1.0, 0.0, None),
    "fast_growth_mod":    ({"r": (0.05, 0.10, "log")}, 1.0, 1.0, 1.0, 0.0, None),
    "fast_growth_sev":    ({"r": (0.08, 0.15, "log")}, 1.0, 1.0, 1.0, 0.0, None),
    "noise_x2":           ({}, 1.0, 1.0, 2.0, 0.0, None),
    "noise_x3":           ({}, 1.0, 1.0, 3.0, 0.0, None),
    "bias_10pct":         ({}, 1.0, 1.0, 1.0, 0.10, None),
    "bias_20pct":         ({}, 1.0, 1.0, 1.0, 0.20, None),
    "heavy_tails_df5":    ({}, 1.0, 1.0, 1.0, 0.0, 5.0),
    "heavy_tails_df3":    ({}, 1.0, 1.0, 1.0, 0.0, 3.0),
    "competition_mod":    ({}, 2.0, 0.5, 1.0, 0.0, None),
    "competition_sev":    ({}, 3.0, 0.33, 1.0, 0.0, None),
    "high_resist_mod":    ({"f0": (0.1, 0.3, "log")}, 1.0, 1.0, 1.0, 0.0, None),
    "high_resist_sev":    ({"f0": (0.2, 0.5, "log")}, 1.0, 1.0, 1.0, 0.0, None),
    "combo_bias_noise":   ({}, 1.0, 1.0, 2.0, 0.10, None),
    "combo_fast_weak":    ({"r": (0.05, 0.10, "log"), "m": (4.0, 6.0, "log")}, 1.0, 1.0, 1.0, 0.0, None),
}
ARM_KEYS = ("primary_containment", "winner_containment_wide", "primary_nonfilter", "winner_periodic", "reference_zhang")
PROTOCOL = {"suites": {k: [v[0], v[1], v[2], v[3], v[4], v[5]] for k, v in SUITES.items()}, "n": N,
            "n_boot": N_BOOT, "arms": ARM_KEYS, "seed_rule": "20261200 + k; noise 1_200_000 + 10k + {1, 2}",
            "harm": "CI95 upper < 0 vs primary_nonfilter", "benefit": "CI95 lower > 0", "event_mode": "refine"}
PROTOCOL_HASH = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True, default=str).encode()).hexdigest()


def build(k, name, n):
    ov, a_sr, a_rs, sdm, bias, df = SUITES[name]
    s = tc.draw_scenarios(n, 20261200 + k, ranges=dict(R, **ov), a_sr=a_sr, a_rs=a_rs)
    if sdm != 1.0:
        s["sigma_true"] = sdm * s["sigma"]
    if bias:
        s["bias"] = np.full(n, np.log(1.0 + bias))
    if df:
        s["noise_df"] = np.full(n, df)
    return s


def _job(args):
    spec, k, name, n = args
    noise = [1_200_000 + 10 * k + 1, 1_200_000 + 10 * k + 2]
    return (name, tuple(spec)), tc.evaluate(build(k, name, n), tc.make_policy(tuple(spec)), noise)["ttp_reps"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-scale", type=float, default=1.0)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    if args.n_scale != 1.0 and out_dir == OUT_DIR:
        sys.exit("reduced-size runs must write elsewhere (--out)")
    out_dir.mkdir(parents=True, exist_ok=True)
    sp = out_dir / "sealed_results.json"
    if sp.exists():
        sys.exit(f"refusing to overwrite {sp}: the sealed grid is scored once")
    frozen = json.loads(FROZEN.read_text())
    arms = {"primary_containment": frozen["primary_containment"],
            "winner_containment_wide": frozen["family_winners"]["containment_wide"],
            "primary_nonfilter": frozen["primary_nonfilter"],
            "winner_periodic": frozen["family_winners"]["periodic"],
            "reference_zhang": ["raw", 1.0, 0.5]}
    n = max(7, int(N * args.n_scale))
    t0 = time.time()
    jobs = [(spec, k, name, n) for k, name in enumerate(SUITES) for spec in
            dict.fromkeys(tuple(v) for v in arms.values())]
    jobs.sort(key=lambda j: j[0][0] != "risk")              # start the slow filter jobs first
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        reps = dict(ex.map(_job, jobs))
    stats = {"mean": lambda x: float(np.mean(x)), "cvar10_scenario_mean": tc.cvar_scenario_mean,
             "cvar10_joint": tc.cvar_joint}
    out = {}
    for name in SUITES:
        base = reps[(name, tuple(arms["primary_nonfilter"]))]
        row = {}
        for a, spec in arms.items():
            x = reps[(name, tuple(spec))]
            r = {k: f(x) for k, f in stats.items()}
            if a != "primary_nonfilter":
                d, ci = tc.hier_bootstrap_diff(x, base, stats["mean"], n_boot=N_BOOT, seed=168)
                r["mean_diff_vs_primary_nonfilter"] = {"diff": d, "ci95": ci,
                                                       "harm": ci[1] < 0, "benefit": ci[0] > 0}
            row[a] = r
        out[name] = row
    summary = {}
    for a in ("primary_containment", "winner_containment_wide", "winner_periodic"):
        rows = {s: out[s][a]["mean_diff_vs_primary_nonfilter"] for s in SUITES}
        worst = min(rows, key=lambda s: rows[s]["diff"])
        summary[a] = {"suites_benefit": [s for s in SUITES if rows[s]["benefit"]],
                      "suites_harm": [s for s in SUITES if rows[s]["harm"]],
                      "no_harm_across_grid": not any(rows[s]["harm"] for s in SUITES),
                      "worst_suite": worst, "worst_diff": rows[worst]["diff"], "worst_ci95": rows[worst]["ci95"]}
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {"script": "168_trackc_sealed_wrong_models",
              "label": ("SIMULATION; frozen script-167 arms; pre-set severity grid; sealed, scored once. "
                        "Suites are chosen stress tests, not a sampled uncertainty distribution"),
              "protocol_sha256": PROTOCOL_HASH, "protocol": PROTOCOL, "arms": arms,
              "frozen_protocol_sha256": frozen["protocol_sha256"], "git_commit_at_run": commit,
              "n_scale": args.n_scale, "summary": summary, "suites": out, "runtime_s": time.time() - t0}
    sp.write_text(json.dumps(result, indent=2, default=float))
    for name in SUITES:
        r = out[name]
        print(f"{name:18s} narrow {r['primary_containment']['mean_diff_vs_primary_nonfilter']['diff']:+7.1f} "
              f"wide {r['winner_containment_wide']['mean_diff_vs_primary_nonfilter']['diff']:+7.1f} "
              f"front {r['winner_periodic']['mean_diff_vs_primary_nonfilter']['diff']:+7.1f}")
    print(json.dumps(summary, indent=1, default=float))


if __name__ == "__main__":
    main()
