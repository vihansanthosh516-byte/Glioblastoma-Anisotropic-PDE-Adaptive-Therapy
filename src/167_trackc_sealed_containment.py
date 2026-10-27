#!/usr/bin/env python3
"""
Script 167: Track C containment, frozen protocol with equal tuning budgets, one sealed score.
===========================================================================================
SIMULATION ONLY, assumed parameters, nominal generator (trackc_sim.RANGES). The design came
from earlier exploratory scripts (150-158), so this is a sealed re-test of a post hoc
hypothesis, not a pre-registered confirmation of the original question.

PROTOCOL (frozen in code; its SHA-256 is written before any tuning outcome is seen)
  Families, 12 configs each (equal search budget):
    containment_narrow  RiskContainment, prior m in [6, 20]: p* x horizon
    containment_wide    RiskContainment, prior m in [3, 20]: same grid
    raw_onoff           on/off thresholds on raw readings
    raw_single          single threshold on raw readings
    periodic            a course every q cycles, first course at cycle 0 or 2
  Reference (untuned, not a family): Zhang 1.0/0.5.
  Tuning: one shared draw for every family, seed 20261121, n = 300, noise 1_121_001/002.
  Selection objective: mean capped TTP (exact crossing times, EVENT_MODE "refine").
  Primary contrast: best containment arm (over both containment families) minus best
    non-filter arm (over raw_onoff, raw_single, periodic), each chosen on tuning only.
    Estimand: mean TTP over the nominal scenario x noise distribution.
    Uncertainty: two-stage paired bootstrap (scenarios, then noise replicates), 4,000 draws.
    Success rule: 95% CI lower bound > 0 AND point estimate >= 14 days.
  Secondary (descriptive, no success rule): CVaR_0.10 of scenario-mean TTP and CVaR_0.10
    over scenario x noise outcomes, same bootstrap; every family winner vs the primary
    non-filter arm; max courses per arm.
  Sealed evaluation: seed 20261122, n = 400, noise 1_122_001..004. Scored once.

Stages:  --stage tune   writes frozen_protocol.json (protocol hash + selections)
         --stage seal   requires frozen_protocol.json with a matching hash; refuses to
                        overwrite sealed_results.json
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import hashlib
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

OUT_DIR = PROJECT_ROOT / "output" / "trackc_sealed_containment"
P_GRID = [0.01, 0.02, 0.05, 0.1, 0.2, 0.35]
FAMILIES = {
    "containment_narrow": [("risk", p, h) for p, h in itertools.product(P_GRID, (1, 2))],
    "containment_wide": [("risk", p, h, 0.0, True) for p, h in itertools.product(P_GRID, (1, 2))],
    "raw_onoff": [("raw", a, o) for a, o in [(0.6, 0.3), (0.6, 0.4), (0.6, 0.5), (0.8, 0.3), (0.8, 0.5),
                                              (0.8, 0.7), (0.9, 0.5), (1.0, 0.3), (1.0, 0.5), (1.0, 0.7),
                                              (1.1, 0.7), (1.2, 0.8)]],
    "raw_single": [("raw", a, a) for a in (0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.18)],
    "periodic": [("periodic", q, s) for s in (0, 2) for q in range(1, 7)],
}
CONTAINMENT = ("containment_narrow", "containment_wide")
NONFILTER = ("raw_onoff", "raw_single", "periodic")
ZHANG = ("raw", 1.0, 0.5)
TUNE = {"seed": 20261121, "n": 300, "noise": [1_121_001, 1_121_002]}
SEAL = {"seed": 20261122, "n": 400, "noise": [1_122_001, 1_122_002, 1_122_003, 1_122_004]}
RULE = {"min_gain_days": 14.0, "ci": 0.95, "n_boot": 4000, "boot_seed": 167}
PROTOCOL = {"families": {k: [list(s) for s in v] for k, v in FAMILIES.items()}, "reference": list(ZHANG),
            "tune": TUNE, "seal": SEAL, "rule": RULE, "objective": "mean capped TTP",
            "event_mode": "refine", "dt": 0.25, "generator": "trackc_sim.RANGES (nominal)",
            "primary": "best containment (both families) minus best non-filter (raw_onoff, raw_single, periodic)"}
PROTOCOL_HASH = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True).encode()).hexdigest()


def _commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        return "unknown"


def tune(out_dir, workers, n_scale):
    scen = tc.draw_scenarios(max(7, int(TUNE["n"] * n_scale)), TUNE["seed"])
    specs = [s for v in FAMILIES.values() for s in v]
    with ProcessPoolExecutor(max_workers=workers) as ex:
        reps = dict(ex.map(tc.score_spec_reps, [(s, scen, TUNE["noise"]) for s in specs]))
    mean = {s: float(reps[s].mean()) for s in specs}
    winners = {f: max(v, key=mean.get) for f, v in FAMILIES.items()}
    frozen = {"script": "167_trackc_sealed_containment", "stage": "tune", "protocol_sha256": PROTOCOL_HASH,
              "protocol": PROTOCOL, "git_commit_at_run": _commit(), "n_scale": n_scale,
              "tuning_mean_ttp": {str(list(s)): mean[s] for s in specs},
              "tuning_cvar10_scenario_mean": {str(list(s)): tc.cvar_scenario_mean(reps[s]) for s in specs},
              "family_winners": {f: list(s) for f, s in winners.items()},
              "primary_containment": list(winners[max(CONTAINMENT, key=lambda f: mean[winners[f]])]),
              "primary_nonfilter": list(winners[max(NONFILTER, key=lambda f: mean[winners[f]])])}
    (out_dir / "frozen_protocol.json").write_text(json.dumps(frozen, indent=2))
    print(json.dumps({k: frozen[k] for k in ("family_winners", "primary_containment", "primary_nonfilter")}, indent=1))


def seal(out_dir, workers, n_scale):
    fp, sp = out_dir / "frozen_protocol.json", out_dir / "sealed_results.json"
    if sp.exists():
        sys.exit(f"refusing to overwrite {sp}: the sealed draw is scored once")
    frozen = json.loads(fp.read_text())
    if frozen["protocol_sha256"] != PROTOCOL_HASH:
        sys.exit("protocol in code differs from the frozen protocol; refusing to score the sealed draw")
    if n_scale != frozen["n_scale"]:
        sys.exit("n_scale differs from the tuning stage")
    arms = {f"winner_{f}": tuple(s) for f, s in frozen["family_winners"].items()}
    arms.update(primary_containment=tuple(frozen["primary_containment"]),
                primary_nonfilter=tuple(frozen["primary_nonfilter"]), reference_zhang=ZHANG)
    scen = tc.draw_scenarios(max(7, int(SEAL["n"] * n_scale)), SEAL["seed"])
    specs = list(dict.fromkeys(arms.values()))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        reps = dict(ex.map(tc.score_spec_reps, [(s, scen, SEAL["noise"]) for s in specs]))
    courses = {}
    for s in specs:                                        # budget check, one replicate is enough
        courses[s] = int(tc.simulate(scen, tc.make_policy(s), SEAL["noise"][0])["courses"].max())
        assert courses[s] <= tc.BUDGET
    base = reps[arms["primary_nonfilter"]]
    stats = {"mean": lambda x: float(np.mean(x)), "cvar10_scenario_mean": tc.cvar_scenario_mean,
             "cvar10_joint": tc.cvar_joint}
    out = {}
    for name, s in arms.items():
        row = {"spec": list(s), "max_courses": courses[s]}
        for k, f in stats.items():
            row[k] = f(reps[s])
            if s != arms["primary_nonfilter"]:
                d, ci = tc.hier_bootstrap_diff(reps[s], base, f, n_boot=RULE["n_boot"], seed=RULE["boot_seed"])
                row[f"{k}_diff_vs_primary_nonfilter"] = {"diff": d, "ci95": ci}
        out[name] = row
    prim = out["primary_containment"]["mean_diff_vs_primary_nonfilter"]
    success = bool(prim["ci95"][0] > 0 and prim["diff"] >= RULE["min_gain_days"])
    result = {"script": "167_trackc_sealed_containment", "stage": "seal",
              "label": ("SIMULATION, nominal assumed model. Sealed re-test of a post hoc hypothesis from scripts "
                        "150-158; secondary rows descriptive"),
              "protocol_sha256": PROTOCOL_HASH, "git_commit_at_run": _commit(), "n_scale": n_scale,
              "primary": {"contrast": "primary_containment minus primary_nonfilter, mean TTP", **prim,
                          "success_rule": RULE, "success": success},
              "arms": out, "runtime_s": time.time() - t0}
    sp.write_text(json.dumps(result, indent=2, default=float))
    rows = ["scenario,replicate," + ",".join(arms)]
    for r in range(len(SEAL["noise"])):
        for i in range(base.shape[1]):
            rows.append(f"{i},{r}," + ",".join(f"{reps[arms[a]][r, i]:.4f}" for a in arms))
    (out_dir / "sealed_ttp_by_scenario_replicate.csv").write_text("\n".join(rows) + "\n")
    print(json.dumps(result["primary"], indent=1, default=float))
    for name, row in out.items():
        print(name, row["spec"], round(row["mean"], 1), round(row["cvar10_scenario_mean"], 1), round(row["cvar10_joint"], 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["tune", "seal"], required=True)
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-scale", type=float, default=1.0, help="smoke runs only, with --out elsewhere")
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    if args.n_scale != 1.0 and out_dir == OUT_DIR:
        sys.exit("reduced-size runs must write elsewhere (--out)")
    out_dir.mkdir(parents=True, exist_ok=True)
    (tune if args.stage == "tune" else seal)(out_dir, args.workers, args.n_scale)


if __name__ == "__main__":
    main()
