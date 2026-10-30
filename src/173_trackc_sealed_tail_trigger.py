#!/usr/bin/env python3
"""
Script 173: Track C sealed test of an earlier-acting (tail-aware) containment trigger.
=====================================================================================
SIMULATION ONLY. Follow-up to script 171, whose H1 failed: hybrid minus front-loaded CVaR_0.10
was -4.2 days (CI lower -8.0 < -7) under state-linked missed scans. Forced courses K = 1, 2 changed
nothing, so the tail gap comes later, in the waiting cycles. This test lets the trigger look 3-4
cycles ahead (horizon) so it fires earlier in fast-growing or weak-drug scenarios. Written after
seeing script 171. Missing scans are always told to the controller.

PROTOCOL (fixed in code and committed before tuning or scoring; SHA-256 written to the outputs)
  Stage "tune": grid p* {0.05, 0.10} x horizon {3, 4} x K {0, 4}, plus the script-171 winner
    (0.05, 2, 4) and wide containment (0.10, 2, 0) as references (10 configs). Tuning draws only:
    script-167 tuning draw (seed 20261121, n 300, noise 1_121_001/2) scored twice, nominal and
    with state-linked missed scans 0.1 / 0.4 (told). Objective: highest worst-of-the-two CVaR_0.10
    (scenario mean) among configs whose mean TTP is within 7 days of wide containment in both.
    Winner frozen.
  Stage "seal": 4 suites, n 600 each: nominal, state_linked_mod (0.1/0.4), state_linked_sev
    (0.2/0.8), mixture (script-169 rules). Scenario seed 20261700 + k; noise seeds
    [1_700_000 + 10k + rep, chunk], rep 1-2, chunks of 150. All seeds unused before.
  Arms: trigger (tuned), wide containment (167), raw single 0.6 (comparator), front-loaded (167).
  H1 (tail): trigger minus front-loaded, CVaR_0.10 (scenario mean), 95% lower bound > -7 days in
    all 4 suites (same margin as script 171; not relaxed).
  H2 (mean kept): trigger minus comparator, mean TTP, 95% lower bound > -7 days in all 4 suites.
  Two-stage paired bootstrap (scenarios, then replicates), 4,000 draws, seed 173.
  Descriptive: trigger minus wide (mean, CVaR), trigger minus front-loaded (mean), wide minus
  comparator (mean), wide minus front-loaded (CVaR).
  Scored once: refuses to overwrite the sealed result.

Run:  python src/173_trackc_sealed_tail_trigger.py tune --workers 7
      python src/173_trackc_sealed_tail_trigger.py seal --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import hashlib
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

_spec = importlib.util.spec_from_file_location("s171", PROJECT_ROOT / "src" / "171_trackc_sealed_tail_hybrid.py")
s171 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s171)

FROZEN = PROJECT_ROOT / "output" / "trackc_sealed_containment" / "frozen_protocol.json"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_sealed_tail_trigger"
N, CHUNK, N_REP, N_BOOT, MARGIN = 600, 150, 2, 4000, -7.0
STATE_LINKED = s171.STATE_LINKED
SUITES = s171.SUITES
GRID = [("riskf", p, h, k) for p in (0.05, 0.1) for h in (3, 4) for k in (0, 4)] + \
    [("riskf", 0.05, 2, 4), ("riskf", 0.1, 2, 0)]
BASE = ("riskf", 0.1, 2, 0)
TUNE = {"seed": 20261121, "n": 300, "noise": [1_121_001, 1_121_002]}
TUNE_MISS = (0.1, 0.4)
SCEN_SEED, NOISE_SEED = 20261700, 1_700_000
PROTOCOL = {"suites": SUITES, "state_linked": STATE_LINKED, "n": N, "chunk": CHUNK, "n_rep": N_REP,
            "n_boot": N_BOOT, "boot_seed": 173, "margin": MARGIN, "grid": GRID, "tune": TUNE,
            "tune_miss": TUNE_MISS,
            "tune_rule": "max min(CVaR10 nominal, CVaR10 state-linked) s.t. mean >= mean(wide) - 7 in both",
            "event_mode": "refine", "dt": 0.25, "mixture": "script-169 rules, new seed, missing scans told",
            "seed_rule": "scenario 20261700 + k; noise [1_700_000 + 10k + rep, chunk]",
            "H1": "trigger - front_loaded, CVaR10 scenario mean, CI95 lower > -7 in all 4 suites",
            "H2": "trigger - comparator, mean TTP, CI95 lower > -7 in all 4 suites"}
PROTOCOL_HASH = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True, default=str).encode()).hexdigest()


def build(k, name, n):
    """Script-171 suite builder on new scenario seeds (s171.build adds 20261500)."""
    return s171.build(k + SCEN_SEED - 20261500, name, n)


def _job(args):
    spec, k, name, n, c = args
    sub = tc.subset(build(k, name, n), np.arange(c * CHUNK, min(n, (c + 1) * CHUNK)))
    pol = tc.make_policy(tuple(spec))
    reps = np.stack([tc.simulate(sub, pol, [NOISE_SEED + 10 * k + r, c])["ttp"] for r in range(1, N_REP + 1)])
    return (name, tuple(spec), c), reps


def tune(out_dir, workers, n_scale):
    fp = out_dir / "frozen_trigger.json"
    if fp.exists():
        sys.exit(f"refusing to overwrite {fp}")
    n_t = max(7, int(TUNE["n"] * n_scale))
    nom = tc.draw_scenarios(n_t, TUNE["seed"])
    nom["miss_visible"] = np.ones(n_t, bool)
    draws = {"nominal": nom,
             "state_linked": dict(nom, miss_p=np.full(n_t, TUNE_MISS[0]), miss_p_hi=np.full(n_t, TUNE_MISS[1]))}
    jobs = [(s, d, TUNE["noise"]) for d in draws.values() for s in GRID]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        res = list(ex.map(tc.score_spec_reps, jobs))
    reps = {(dn, s): res[i * len(GRID) + g][1] for i, dn in enumerate(draws) for g, s in enumerate(GRID)}
    mean = {key: float(v.mean()) for key, v in reps.items()}
    cvar = {key: float(tc.cvar_scenario_mean(v)) for key, v in reps.items()}
    ok = [s for s in GRID if all(mean[(d, s)] >= mean[(d, BASE)] - 7.0 for d in draws)]
    worst = {s: min(cvar[(d, s)] for d in draws) for s in GRID}
    win = max(ok, key=worst.get)
    fp.write_text(json.dumps({"script": "173_trackc_sealed_tail_trigger", "stage": "tune",
                              "protocol_sha256": PROTOCOL_HASH, "git_commit_at_run": s171._commit(),
                              "n_scale": n_scale,
                              "tuning_mean_ttp": {f"{d} {list(s)}": mean[(d, s)] for d in draws for s in GRID},
                              "tuning_cvar10_scenario_mean": {f"{d} {list(s)}": cvar[(d, s)]
                                                              for d in draws for s in GRID},
                              "eligible": [list(s) for s in ok], "trigger_winner": list(win)}, indent=2))
    for s in GRID:
        print(s, *(f"{d} {mean[(d, s)]:.1f} {cvar[(d, s)]:.1f}" for d in draws), "eligible" if s in ok else "")
    print("trigger winner", win)


def seal(out_dir, workers, n_scale):
    sp = out_dir / "sealed_results.json"
    if sp.exists():
        sys.exit(f"refusing to overwrite {sp}: the sealed test is scored once")
    trg = json.loads((out_dir / "frozen_trigger.json").read_text())
    if trg["protocol_sha256"] != PROTOCOL_HASH or trg["n_scale"] != n_scale:
        sys.exit("protocol or n_scale differs from the tuning stage; refusing to score")
    frozen = json.loads(FROZEN.read_text())
    arms = {"trigger": trg["trigger_winner"], "wide_containment": list(BASE),
            "comparator_raw_single": frozen["primary_nonfilter"],
            "front_loaded": frozen["family_winners"]["periodic"]}
    n = max(7, int(N * n_scale))
    n_chunks = -(-n // CHUNK)
    t0 = time.time()
    specs = list(dict.fromkeys(tuple(v) for v in arms.values()))
    jobs = [(spec, k, nm, n, c) for k, nm in enumerate(SUITES) for c in range(n_chunks) for spec in specs]
    jobs.sort(key=lambda j: {"riskf": 0}.get(j[0][0], 1))      # slow filter jobs first
    with ProcessPoolExecutor(max_workers=workers) as ex:
        parts = dict(ex.map(_job, jobs))
    reps = {(nm, a): np.concatenate([parts[(nm, tuple(spec), c)] for c in range(n_chunks)], 1)
            for nm in SUITES for a, spec in arms.items()}
    stats = {"mean": lambda x: float(np.mean(x)), "cvar10_scenario_mean": tc.cvar_scenario_mean,
             "cvar10_joint": tc.cvar_joint}

    def diff(nm, a, b, stat):
        d, ci = tc.hier_bootstrap_diff(reps[(nm, a)], reps[(nm, b)], stats[stat], n_boot=N_BOOT, seed=173)
        return {"diff": d, "ci95": ci}

    out = {}
    for nm in SUITES:
        row = {"n": n}
        for a in arms:
            row[a] = {k: f(reps[(nm, a)]) for k, f in stats.items()}
        row["H1_trigger_minus_front_cvar"] = diff(nm, "trigger", "front_loaded", "cvar10_scenario_mean")
        row["H2_trigger_minus_comparator_mean"] = diff(nm, "trigger", "comparator_raw_single", "mean")
        row["wide_minus_comparator_mean"] = diff(nm, "wide_containment", "comparator_raw_single", "mean")
        row["trigger_minus_wide_mean"] = diff(nm, "trigger", "wide_containment", "mean")
        row["trigger_minus_wide_cvar"] = diff(nm, "trigger", "wide_containment", "cvar10_scenario_mean")
        row["trigger_minus_front_mean"] = diff(nm, "trigger", "front_loaded", "mean")
        row["wide_minus_front_cvar"] = diff(nm, "wide_containment", "front_loaded", "cvar10_scenario_mean")
        out[nm] = row

    def rule(key):
        return {"pass": all(out[nm][key]["ci95"][0] > MARGIN for nm in SUITES),
                "failing": [nm for nm in SUITES if out[nm][key]["ci95"][0] <= MARGIN]}

    verdict = {"H1": rule("H1_trigger_minus_front_cvar"), "H2": rule("H2_trigger_minus_comparator_mean")}
    result = {"script": "173_trackc_sealed_tail_trigger",
              "label": ("SIMULATION; trigger tuned on the script-167 tuning draw (nominal and "
                        "state-linked); frozen 167 arms; new seeds; sealed, scored once. Written "
                        "after seeing scripts 167-171"),
              "protocol_sha256": PROTOCOL_HASH, "protocol": PROTOCOL, "arms": arms,
              "git_commit_at_run": s171._commit(), "n_scale": n_scale, "verdict": verdict, "suites": out,
              "runtime_s": time.time() - t0}
    sp.write_text(json.dumps(result, indent=2, default=float))
    for nm in SUITES:
        r = out[nm]
        print(f"{nm:18s} " + " ".join(f"{k} {r[k]['diff']:+6.1f} {np.round(r[k]['ci95'], 1)}"
                                      for k in r if isinstance(r[k], dict) and "diff" in r[k]))
    print(json.dumps(verdict, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("tune", "seal"))
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-scale", type=float, default=1.0)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    if args.n_scale != 1.0 and out_dir == OUT_DIR:
        sys.exit("reduced-size runs must write elsewhere (--out)")
    out_dir.mkdir(parents=True, exist_ok=True)
    (tune if args.stage == "tune" else seal)(out_dir, args.workers, args.n_scale)


if __name__ == "__main__":
    main()
