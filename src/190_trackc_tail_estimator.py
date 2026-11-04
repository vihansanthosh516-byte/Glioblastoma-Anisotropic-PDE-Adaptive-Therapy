#!/usr/bin/env python3
"""
Script 190: Track C Stage 2, test A: how the worst-decile score moves with the number of noise repeats.
=======================================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Plan: paper/trackc_restart_plan_2026-11-03.md.

Question. In the sealed 167 per-scenario file, the worst-decile score of the per-scenario mean (CVaR10) of
wide containment rises with the number of noise repeats (171.1, 174.5, 176.0 at 1, 2, 4 repeats) and stays
flat for front-loaded dosing, which ignores the scans. Scripts 171 and 173 used 2 repeats. This script
re-measures that effect on fresh dev seeds with up to 16 repeats, for five frozen rules, and compares:
plug-in CVaR10 of scenario means, a split-sample version, the joint CVaR10 of single outcomes, and the
paired difference to front-loaded dosing. It also saves every single outcome so that nothing has to be
re-run to rebuild a confidence interval, and it saves the number of courses used (the cost of each rule).

Nothing is tuned here. The five rules are the ones frozen in scripts 167 and 171. The simulator,
policies and noise model are unchanged (src/trackc_sim.py). This is not a test of any hypothesis.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev):
    scenarios  seed 20261800, n = 200 (Latin hypercube of the assumed nominal ranges)
    noise      keys [1_800_000 + r, chunk], r = 1..16, chunks of 50 scenarios; the same key in every rule
Written to output/trackc_stage2_A_tail_estimator/ (manifest.json first, then the data, then results.json).

Run:  python src/190_trackc_tail_estimator.py --workers 7
      python src/190_trackc_tail_estimator.py --smoke        (tiny run into a *_smoke folder)
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_manifest as man  # noqa: E402
import trackc_sim as tc  # noqa: E402
import trackc_tail_estimators as te  # noqa: E402

SCRIPT = "src/190_trackc_tail_estimator.py"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_stage2_A_tail_estimator"
SCEN_SEED, NOISE_BASE = 20261800, 1_800_000
N, CHUNK, R_MAX = 200, 50, 16
BASELINE = "front_loaded_q1"
ARMS = {
    "front_loaded_q1": ("periodic", 1, 0),                 # open loop: ignores the scans
    "raw_single_0.6": ("raw", 0.6, 0.6),                   # 167 non-filter comparator
    "wide_containment_p0.10_h2": ("risk", 0.1, 2, 0.0, True),   # 167 wide-prior winner
    "narrow_containment_p0.05_h2": ("risk", 0.05, 2),      # 167 primary arm
    "hybrid_front4_wide_p0.05_h2": ("riskf", 0.05, 2, 4),  # 171 / 173 arm
}
ALPHA = 0.10
N_BOOT, N_SUBSETS, BOOT_SEED = 1000, 200, 190
WEAKEST = [
    "n = 200 scenarios, so the worst decile is 20 scenarios; levels are noisy, the shape of the curve is the point",
    "assumed model ranges, one ODE family, toy time scale (mean TTP 600-700 d against about 210 d median PFS in patients)",
    "the five rules are frozen from earlier tuning on the 167 draw; none is tuned here",
    "the inverse-R extrapolation is a model of the bias, not a theorem; read its residual before using it",
    "joint CVaR10 measures the realised tail of one patient and includes noise-driven misdecisions; the plug-in CVaR10 of scenario means averages them away",
]


def _smoke(args):
    global N, CHUNK, R_MAX, N_BOOT, N_SUBSETS, OUT_DIR
    if args.smoke:
        N, CHUNK, R_MAX, N_BOOT, N_SUBSETS = 12, 6, 4, 50, 20
        OUT_DIR = PROJECT_ROOT / "output" / "trackc_stage2_A_tail_estimator_smoke"


def _job(args):
    arm, c, n, chunk, r_max = args
    scen = tc.draw_scenarios(n, SCEN_SEED)
    sub = tc.subset(scen, np.arange(c * chunk, min(n, (c + 1) * chunk)))
    pol = tc.make_policy(ARMS[arm])
    ttp, courses = [], []
    for r in range(1, r_max + 1):
        out = tc.simulate(sub, pol, [NOISE_BASE + r, c])
        ttp.append(out["ttp"])
        courses.append(out["courses"])
    return (arm, c), (np.stack(ttp), np.stack(courses))


def analyse(ttp, courses, r_max):
    """ttp, courses: {arm: (r_max, n)}. Returns the results dict (no I/O)."""
    Rs = [r for r in (1, 2, 4, 8, 16) if r <= r_max]
    rng = np.random.default_rng(BOOT_SEED)
    res = {"Rs": Rs, "alpha": ALPHA, "arms": {}}
    base = ttp[BASELINE]
    for arm, X in ttp.items():
        plug = te.curve(X, Rs, te.plugin, rng, N_SUBSETS)
        joint = te.curve(X, Rs, te.joint, rng, N_SUBSETS)
        split = te.curve(X, [r for r in Rs if r >= 2], te.split_sample, rng, max(5, N_SUBSETS // 10))
        fit_R = [r for r in Rs if r >= 2]
        c_inf, a, rms = te.fit_inverse_R(fit_R, [plug[r][0] for r in fit_R]) if len(fit_R) >= 2 \
            else (float("nan"),) * 3
        entry = {"mean_ttp": float(X.mean()), "mean_courses": float(courses[arm].mean()),
                 "capped_fraction": float((X >= tc.HORIZON - 0.01).mean()),
                 "plugin_cvar10": {str(r): plug[r] for r in Rs},
                 "joint_cvar10": {str(r): joint[r] for r in Rs},
                 "split_cvar10": {str(r): split[r] for r in split},
                 "inverse_R_fit": {"c_at_infinite_repeats": c_inf, "a": a, "rms_residual": rms,
                                   "fitted_on_R": fit_R}}
        if arm != BASELINE:
            pr = te.paired_regret(X, base)
            d_plug, ci_plug = te.scenario_bootstrap(X, base, te.plugin, N_BOOT, BOOT_SEED)
            d_joint, ci_joint = te.scenario_bootstrap(X, base, te.joint, N_BOOT, BOOT_SEED)
            d_mean, ci_mean = te.scenario_bootstrap(X, base, lambda Z: float(Z.mean()), N_BOOT, BOOT_SEED)
            entry["vs_front_loaded"] = {
                "mean_diff": {"point": d_mean, "ci95": ci_mean},
                "plugin_cvar10_diff_all_repeats": {"point": d_plug, "ci95": ci_plug},
                "joint_cvar10_diff_all_repeats": {"point": d_joint, "ci95": ci_joint},
                "plugin_cvar10_diff_by_R": {str(r): plug[r][0] - res["arms"][BASELINE]["plugin_cvar10"][str(r)][0]
                                            for r in Rs} if BASELINE in res["arms"] else None,
                "paired_regret_all_repeats": pr,
                "courses_used_diff": float(courses[arm].mean() - courses[BASELINE].mean())}
        res["arms"][arm] = entry
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    _smoke(args)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    protocol = {"script": "190_trackc_tail_estimator", "stage": "2A", "smoke": bool(args.smoke), "n": N,
                "chunk": CHUNK, "r_max": R_MAX, "scenario_seed": SCEN_SEED, "noise_base": NOISE_BASE,
                "arms": {k: list(v) for k, v in ARMS.items()}, "baseline": BASELINE, "alpha": ALPHA,
                "n_boot": N_BOOT, "n_subsets": N_SUBSETS, "boot_seed": BOOT_SEED}
    manifest = man.build_manifest("190_trackc_tail_estimator", protocol, purpose="stage2_dev",
                                  scenario=(SCEN_SEED, SCEN_SEED), noise=(NOISE_BASE + 1, NOISE_BASE + R_MAX),
                                  script=SCRIPT, require_clean=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_tail_estimators.py",))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; "
          f"git clean {manifest['git']['clean']}", flush=True)

    t0 = time.time()
    n_chunks = -(-N // CHUNK)
    jobs = [(arm, c, N, CHUNK, R_MAX) for arm in ARMS for c in range(n_chunks)]
    jobs.sort(key=lambda j: {"periodic": 2, "raw": 1}.get(ARMS[j[0]][0], 0))      # slow filter rules first
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = dict(ex.map(_job, jobs))
    ttp = {arm: np.concatenate([parts[(arm, c)][0] for c in range(n_chunks)], axis=1) for arm in ARMS}
    courses = {arm: np.concatenate([parts[(arm, c)][1] for c in range(n_chunks)], axis=1) for arm in ARMS}

    scen = tc.draw_scenarios(N, SCEN_SEED)
    cols = [k for k in ("r", "m", "f0", "c", "eps", "B0", "sigma") if k in scen]
    with open(OUT_DIR / "scenarios.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["scenario"] + cols)
        for i in range(N):
            w.writerow([i] + [f"{scen[k][i]:.10g}" for k in cols])
    with open(OUT_DIR / "ttp_by_scenario_replicate.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["arm", "scenario", "replicate", "ttp_days", "courses"])
        for arm in ARMS:
            for r in range(R_MAX):
                for i in range(N):
                    w.writerow([arm, i, r + 1, f"{ttp[arm][r, i]:.6f}", int(courses[arm][r, i])])

    res = analyse(ttp, courses, R_MAX)
    res.update({"script": "190_trackc_tail_estimator",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; not a test of any hypothesis",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for arm, e in res["arms"].items():
        print(f"== {arm}: mean TTP {e['mean_ttp']:.1f} d, mean courses {e['mean_courses']:.2f}")
        print("   plug-in CVaR10 by R :", {r: round(v[0], 1) for r, v in e["plugin_cvar10"].items()})
        print("   joint   CVaR10 by R :", {r: round(v[0], 1) for r, v in e["joint_cvar10"].items()})
        print("   split   CVaR10 by R :", {r: round(v[0], 1) for r, v in e["split_cvar10"].items()})
        f = e["inverse_R_fit"]
        print(f"   inverse-R fit: c = {f['c_at_infinite_repeats']:.1f}, a = {f['a']:.1f}, rms = {f['rms_residual']:.2f}")


if __name__ == "__main__":
    main()
