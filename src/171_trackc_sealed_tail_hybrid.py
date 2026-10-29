#!/usr/bin/env python3
"""
Script 171: Track C sealed test of a front-load-then-contain hybrid and of state-linked missed scans.
====================================================================================================
SIMULATION ONLY. Follow-up to scripts 167-170, written after seeing that front-loaded dosing keeps
the better lower tail (CVaR_0.10) while wide containment wins on mean TTP. The hybrid gives K
courses in cycles 0..K-1, then follows wide containment (trackc_sim "riskf"). Missing scans are
always visible to the controller here (script-170 design).

PROTOCOL (fixed in code and committed before tuning or scoring; SHA-256 written to the outputs)
  Stage "tune": hybrid grid p* {0.05, 0.10} x K {0, 1, 2, 3, 4, 6} (12 configs, the same budget
    as each script-167 family; K = 0, p* 0.10 is frozen wide containment), horizon 2, wide prior,
    on the script-167 tuning draw (seed 20261121, n 300, noise 1_121_001/2). Objective: highest
    CVaR_0.10 (scenario mean) among configs whose mean TTP is within 7 days of K = 0, p* 0.10.
    Winner frozen.
  Stage "seal": 4 suites, n 600 each.
    nominal           script-167 generator
    state_linked_mod  scan lost w.p. 0.1 while true burden < B0, 0.4 while >= B0 (told)
    state_linked_sev  0.2 / 0.8 (told)
    mixture           script-169 mixture rules (noise, bias, m, r, competition, one new axis
                      w.p. 0.5), new draw; missing scans told
    Scenario seed 20261500 + k; noise seeds [1_500_000 + 10k + rep, chunk], rep 1-2, chunks of
    150 (same chunks and seeds for every arm). All seeds unused before.
  Arms: hybrid (tuned), wide containment (167), raw single 0.6 (comparator), front-loaded (167).
  H1 (tail): hybrid minus front-loaded, CVaR_0.10 (scenario mean), 95% lower bound > -7 days in
    all 4 suites.
  H2 (mean kept): hybrid minus comparator, mean TTP, 95% lower bound > -7 days in all 4 suites.
  H3 (state-linked loss): wide containment minus comparator, mean TTP, 95% lower bound > -7 days
    in both state-linked suites.
  Two-stage paired bootstrap (scenarios, then replicates), 4,000 draws, seed 171.
  Descriptive: hybrid minus wide (mean, CVaR), hybrid minus front-loaded (mean), CVaR joint.
  Scored once: refuses to overwrite the sealed result.

Run:  python src/171_trackc_sealed_tail_hybrid.py tune --workers 7
      python src/171_trackc_sealed_tail_hybrid.py seal --workers 7
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

_spec = importlib.util.spec_from_file_location("s169", PROJECT_ROOT / "src" / "169_trackc_sealed_new_axes.py")
s169 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s169)

FROZEN = PROJECT_ROOT / "output" / "trackc_sealed_containment" / "frozen_protocol.json"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_sealed_tail_hybrid"
N, CHUNK, N_REP, N_BOOT, MARGIN = 600, 150, 2, 4000, -7.0
STATE_LINKED = {"state_linked_mod": (0.1, 0.4), "state_linked_sev": (0.2, 0.8)}
SUITES = ["nominal", *STATE_LINKED, "mixture"]
HYBRID_GRID = [("riskf", p, 2, k) for p in (0.05, 0.1) for k in (0, 1, 2, 3, 4, 6)]
BASE = ("riskf", 0.1, 2, 0)
TUNE = {"seed": 20261121, "n": 300, "noise": [1_121_001, 1_121_002]}
PROTOCOL = {"suites": SUITES, "state_linked": STATE_LINKED, "n": N, "chunk": CHUNK, "n_rep": N_REP,
            "n_boot": N_BOOT, "boot_seed": 171, "margin": MARGIN, "hybrid_grid": HYBRID_GRID,
            "tune": TUNE, "tune_rule": "max CVaR10 scenario mean s.t. mean >= mean(K=0, p*=0.1) - 7",
            "event_mode": "refine", "dt": 0.25, "mixture": "script-169 rules, new seed, missing scans told",
            "seed_rule": "scenario 20261500 + k; noise [1_500_000 + 10k + rep, chunk]",
            "H1": "hybrid - front_loaded, CVaR10 scenario mean, CI95 lower > -7 in all 4 suites",
            "H2": "hybrid - comparator, mean TTP, CI95 lower > -7 in all 4 suites",
            "H3": "wide - comparator, mean TTP, CI95 lower > -7 in both state-linked suites"}
PROTOCOL_HASH = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True, default=str).encode()).hexdigest()


def build(k, name, n):
    seed = 20261500 + k
    if name == "mixture":
        rng = np.random.default_rng([seed, 2])
        comp = rng.integers(0, 2, n)
        s = s169._blank(tc.draw_scenarios(n, seed, ranges=dict(tc.RANGES, m=(3.0, 20.0, "log"),
                                                                r=(0.001, 0.10, "log")),
                                          a_sr=np.where(comp, 2.0, 1.0), a_rs=np.where(comp, 0.5, 1.0)), n)
        s["sigma_true"] = s["sigma"] * rng.choice([1.0, 2.0, 3.0], n)
        s["bias"] = np.log1p(rng.choice([0.0, 0.1, 0.2], n))
        has = rng.random(n) < 0.5
        axis = rng.integers(0, len(s169.AXES), n)
        sev = rng.integers(0, 2, n)
        for a_i, a in enumerate(s169.AXES):
            for v in (0, 1):
                s169._apply_axis(s, np.flatnonzero(has & (axis == a_i) & (sev == v)), a, v, rng)
    else:
        s = tc.draw_scenarios(n, seed)
        if name in STATE_LINKED:
            lo, hi = STATE_LINKED[name]
            s["miss_p"], s["miss_p_hi"] = np.full(n, lo), np.full(n, hi)
    s["miss_visible"] = np.ones(n, bool)
    return s


def _job(args):
    spec, k, name, n, c = args
    sub = tc.subset(build(k, name, n), np.arange(c * CHUNK, min(n, (c + 1) * CHUNK)))
    pol = tc.make_policy(tuple(spec))
    reps = np.stack([tc.simulate(sub, pol, [1_500_000 + 10 * k + r, c])["ttp"] for r in range(1, N_REP + 1)])
    return (name, tuple(spec), c), reps


def _commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        return "unknown"


def tune(out_dir, workers, n_scale):
    fp = out_dir / "frozen_hybrid.json"
    if fp.exists():
        sys.exit(f"refusing to overwrite {fp}")
    scen = tc.draw_scenarios(max(7, int(TUNE["n"] * n_scale)), TUNE["seed"])
    with ProcessPoolExecutor(max_workers=workers) as ex:
        reps = dict(ex.map(tc.score_spec_reps, [(s, scen, TUNE["noise"]) for s in HYBRID_GRID]))
    mean = {s: float(reps[s].mean()) for s in HYBRID_GRID}
    cvar = {s: float(tc.cvar_scenario_mean(reps[s])) for s in HYBRID_GRID}
    ok = [s for s in HYBRID_GRID if mean[s] >= mean[BASE] - 7.0]
    win = max(ok, key=cvar.get)
    fp.write_text(json.dumps({"script": "171_trackc_sealed_tail_hybrid", "stage": "tune",
                              "protocol_sha256": PROTOCOL_HASH, "git_commit_at_run": _commit(),
                              "n_scale": n_scale,
                              "tuning_mean_ttp": {str(list(s)): mean[s] for s in HYBRID_GRID},
                              "tuning_cvar10_scenario_mean": {str(list(s)): cvar[s] for s in HYBRID_GRID},
                              "eligible": [list(s) for s in ok], "hybrid_winner": list(win)}, indent=2))
    for s in HYBRID_GRID:
        print(s, round(mean[s], 1), round(cvar[s], 1), "eligible" if s in ok else "")
    print("hybrid winner", win)


def seal(out_dir, workers, n_scale):
    sp = out_dir / "sealed_results.json"
    if sp.exists():
        sys.exit(f"refusing to overwrite {sp}: the sealed test is scored once")
    hyb = json.loads((out_dir / "frozen_hybrid.json").read_text())
    if hyb["protocol_sha256"] != PROTOCOL_HASH or hyb["n_scale"] != n_scale:
        sys.exit("protocol or n_scale differs from the tuning stage; refusing to score")
    frozen = json.loads(FROZEN.read_text())
    arms = {"hybrid": hyb["hybrid_winner"], "wide_containment": list(BASE),
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
        d, ci = tc.hier_bootstrap_diff(reps[(nm, a)], reps[(nm, b)], stats[stat], n_boot=N_BOOT, seed=171)
        return {"diff": d, "ci95": ci}

    out = {}
    for nm in SUITES:
        row = {"n": n}
        for a in arms:
            row[a] = {k: f(reps[(nm, a)]) for k, f in stats.items()}
        row["H1_hybrid_minus_front_cvar"] = diff(nm, "hybrid", "front_loaded", "cvar10_scenario_mean")
        row["H2_hybrid_minus_comparator_mean"] = diff(nm, "hybrid", "comparator_raw_single", "mean")
        row["wide_minus_comparator_mean"] = diff(nm, "wide_containment", "comparator_raw_single", "mean")
        row["hybrid_minus_wide_mean"] = diff(nm, "hybrid", "wide_containment", "mean")
        row["hybrid_minus_wide_cvar"] = diff(nm, "hybrid", "wide_containment", "cvar10_scenario_mean")
        row["hybrid_minus_front_mean"] = diff(nm, "hybrid", "front_loaded", "mean")
        row["wide_minus_front_cvar"] = diff(nm, "wide_containment", "front_loaded", "cvar10_scenario_mean")
        out[nm] = row

    def rule(key, names):
        rows = {nm: out[nm][key] for nm in names}
        return {"pass": all(r["ci95"][0] > MARGIN for r in rows.values()),
                "failing": [nm for nm, r in rows.items() if r["ci95"][0] <= MARGIN]}

    verdict = {"H1": rule("H1_hybrid_minus_front_cvar", SUITES),
               "H2": rule("H2_hybrid_minus_comparator_mean", SUITES),
               "H3": rule("wide_minus_comparator_mean", list(STATE_LINKED))}
    result = {"script": "171_trackc_sealed_tail_hybrid",
              "label": ("SIMULATION; hybrid tuned on the script-167 tuning draw; frozen 167 arms; "
                        "state-linked missed scans; new seeds; sealed, scored once. Written after "
                        "seeing scripts 167-170"),
              "protocol_sha256": PROTOCOL_HASH, "protocol": PROTOCOL, "arms": arms,
              "git_commit_at_run": _commit(), "n_scale": n_scale, "verdict": verdict, "suites": out,
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
