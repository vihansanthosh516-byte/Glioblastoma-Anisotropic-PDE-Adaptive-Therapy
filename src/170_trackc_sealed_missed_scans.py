#!/usr/bin/env python3
"""
Script 170: Track C sealed missed-scan test with a controller that knows a scan is missing.
=========================================================================================
SIMULATION ONLY. Follow-up to script 169, where wide containment failed H1 with 40% missed
readings (-4.7 days, CI -16.6 to 7.0): there the controller was handed the last reading again
as if it were new. Here a missed scan reaches the policy as NaN (`miss_visible`): the particle
filter skips its weight update; the raw rules reuse their last reading (so they act exactly as
in script 169). No arm is retuned: specs come from scripts 167 and 169.

PROTOCOL (fixed in code and committed before scoring; SHA-256 written to the output)
  Suites: scan lost with probability 0.2 / 0.4 / 0.6 (cycle 0 is never lost); nominal ranges.
  n 600 per suite; scenario seed 20261400 + k; noise seeds [1_400_000 + 10k + rep, chunk],
  rep 1-2, chunks of 150 (same chunks and seeds for every arm). All seeds unused before.
  Arms: wide containment told about missing scans (main), the same rule not told (script-169
    behaviour), raw single 0.6 (comparator), raw on/off 0.6/0.5, front-loaded periodic,
    ablation filtered single threshold 0.6 (wide prior, told).
  H: main arm minus comparator, mean capped TTP, 95% lower bound > -7 days in all 3 suites.
  Two-stage paired bootstrap (scenarios, then replicates), 4,000 draws, seed 170.
  Descriptive: main minus not-told, main minus front-loaded, CVaR_0.10 (both estimands).
  Scored once: refuses to overwrite the sealed result.

Run:  python src/170_trackc_sealed_missed_scans.py --workers 7
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

FROZEN_167 = PROJECT_ROOT / "output" / "trackc_sealed_containment" / "frozen_protocol.json"
FROZEN_169 = PROJECT_ROOT / "output" / "trackc_sealed_new_axes" / "frozen_ablation.json"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_sealed_missed_scans"
N, CHUNK, N_REP, N_BOOT, MARGIN = 600, 150, 2, 4000, -7.0
SUITES = {"missed_20": 0.2, "missed_40": 0.4, "missed_60": 0.6}
PROTOCOL = {"suites": SUITES, "n": N, "chunk": CHUNK, "n_rep": N_REP, "n_boot": N_BOOT, "boot_seed": 170,
            "margin": MARGIN, "event_mode": "refine", "dt": 0.25,
            "seed_rule": "scenario 20261400 + k; noise [1_400_000 + 10k + rep, chunk]",
            "H": "wide containment (told) - raw single 0.6, mean TTP, CI95 lower > -7 in all 3 suites"}
PROTOCOL_HASH = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True).encode()).hexdigest()


def build(k, name, n, told):
    s = tc.draw_scenarios(n, 20261400 + k)
    s["miss_p"] = np.full(n, SUITES[name])
    if told:
        s["miss_visible"] = np.ones(n, bool)
    return s


def _job(args):
    arm, spec, told, k, name, n, c = args
    sub = tc.subset(build(k, name, n, told), np.arange(c * CHUNK, min(n, (c + 1) * CHUNK)))
    pol = tc.make_policy(tuple(spec))
    reps = np.stack([tc.simulate(sub, pol, [1_400_000 + 10 * k + r, c])["ttp"] for r in range(1, N_REP + 1)])
    return (name, arm, c), reps


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
        sys.exit(f"refusing to overwrite {sp}: the sealed test is scored once")
    f167, f169 = json.loads(FROZEN_167.read_text()), json.loads(FROZEN_169.read_text())
    wide = f167["family_winners"]["containment_wide"]
    arms = {"wide_told": (wide, True), "wide_not_told": (wide, False),
            "comparator_raw_single": (f167["primary_nonfilter"], True),
            "raw_onoff": (f167["family_winners"]["raw_onoff"], True),
            "front_loaded": (f167["family_winners"]["periodic"], True),
            "ablation_filtered_threshold": (f169["ablation_winner"], True)}
    n = max(7, int(N * args.n_scale))
    n_chunks = -(-n // CHUNK)
    t0 = time.time()
    jobs = [(a, spec, told, k, name, n, c) for k, name in enumerate(SUITES) for c in range(n_chunks)
            for a, (spec, told) in arms.items()]
    jobs.sort(key=lambda j: {"risk": 0, "filtw": 1}.get(j[1][0], 2))      # slow filter jobs first
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = dict(ex.map(_job, jobs))
    reps = {(name, a): np.concatenate([parts[(name, a, c)] for c in range(n_chunks)], 1)
            for name in SUITES for a in arms}
    stats = {"mean": lambda x: float(np.mean(x)), "cvar10_scenario_mean": tc.cvar_scenario_mean,
             "cvar10_joint": tc.cvar_joint}
    out = {}
    for name in SUITES:
        base = reps[(name, "comparator_raw_single")]
        row = {"n": n}
        for a in arms:
            x = reps[(name, a)]
            r = {k: f(x) for k, f in stats.items()}
            if a != "comparator_raw_single":
                d, ci = tc.hier_bootstrap_diff(x, base, stats["mean"], n_boot=N_BOOT, seed=170)
                r["mean_diff_vs_comparator"] = {"diff": d, "ci95": ci}
            row[a] = r
        w = reps[(name, "wide_told")]
        for other in ("wide_not_told", "front_loaded", "ablation_filtered_threshold"):
            d, ci = tc.hier_bootstrap_diff(w, reps[(name, other)], stats["mean"], n_boot=N_BOOT, seed=170)
            row[f"wide_told_minus_{other}"] = {"diff": d, "ci95": ci}
        out[name] = row
    h = {name: out[name]["wide_told"]["mean_diff_vs_comparator"] for name in SUITES}
    verdict = {"H_pass": all(r["ci95"][0] > MARGIN for r in h.values()),
               "H_failing_suites": [name for name, r in h.items() if r["ci95"][0] <= MARGIN],
               "H_rows": h}
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    result = {"script": "170_trackc_sealed_missed_scans",
              "label": ("SIMULATION; frozen script-167/169 arms; controller told when a scan is missing; "
                        "new seeds; sealed, scored once. Follow-up written after seeing the script-169 failure"),
              "protocol_sha256": PROTOCOL_HASH, "protocol": PROTOCOL,
              "arms": {a: {"spec": s, "told": t} for a, (s, t) in arms.items()},
              "git_commit_at_run": commit, "n_scale": args.n_scale, "verdict": verdict, "suites": out,
              "runtime_s": time.time() - t0}
    sp.write_text(json.dumps(result, indent=2, default=float))
    for name in SUITES:
        r = out[name]
        print(f"{name:10s} " + " ".join(f"{a} {r[a]['mean_diff_vs_comparator']['diff']:+6.1f}"
                                        f" {np.round(r[a]['mean_diff_vs_comparator']['ci95'], 1)}"
                                        for a in arms if a != "comparator_raw_single"))
    print(json.dumps({k: v for k, v in verdict.items() if k != "H_rows"}, default=float))


if __name__ == "__main__":
    main()
