#!/usr/bin/env python3
"""
Script 192: Track C Stage 2, test C: which assumed worlds reproduce the published standard arm?
===========================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Plan: paper/trackc_restart_plan_2026-11-03.md,
design: paper/trackc_anchored_spec_2026-11-04.md section 4.

Question. The old simulated world has a mean time to progression of 596 to 702 days. The published standard arm
(RTOG 0525, Gilbert 2013) has a median PFS of 5.5 months (95% CI 4.7 to 6.1, i.e. 143 to 186 days), with scans about
every 3 months. For which settings of the parameter ranges ("worlds", src/trackc_calibration.py) does a simulated
standard arm (6 courses, one per 28-day cycle, open loop, scan every 84 days, RANO-like progression) land inside that
interval? And which of the world axes does the anchor constrain?

Method. Draw worlds by Latin hypercube over the axes of src/trackc_calibration.py (r shift, kill-strength scale, f0
shift, B0 centre, cost c, switching mu, reversion nu). Draw N_PAT patients per world. Run the standard arm. A world
passes if the interpolated Kaplan-Meier median PFS lies in 143 to 186 days. The same worlds are re-run under
changed endpoint rules (RANO floor, RANO factor, scan interval, the old 1.2 x baseline rule) to show how much the
answer depends on those assumptions.

Nothing here looks at a dosing policy contrast. The only arm is the open-loop standard arm. Passing is a necessary
check, not evidence that a world is true.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev):
    worlds     LHS seed 20261801
    patients   seed [20261802, world index]
    noise      key [1_810_000, chunk] (the arm ignores readings; the key is only passed to the simulator)
Written to output/trackc_stage2_C_calibration/ (manifest.json first, then worlds.csv, kept_worlds.csv, results.json).

Run:  python src/192_trackc_calibration.py --workers 3
      python src/192_trackc_calibration.py --smoke        (tiny run into a *_smoke folder)
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
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_manifest as man  # noqa: E402
import trackc_sim as tc  # noqa: E402

SCRIPT = "src/192_trackc_calibration.py"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_stage2_C_calibration"
WORLD_SEED, PATIENT_SEED, NOISE_BASE = 20261801, 20261802, 1_810_000
N_WORLDS, N_PAT, CHUNK_WORLDS = 500, 400, 25
OLD_MODEL_N = 2000
VARIANTS = {
    "base": {"endpoint": "rano", "factor": 1.4, "floor": 0.1, "scan": 3},
    "floor0": {"endpoint": "rano", "factor": 1.4, "floor": 0.0, "scan": 3},
    "floor0.2": {"endpoint": "rano", "factor": 1.4, "floor": 0.2, "scan": 3},
    "factor1.25": {"endpoint": "rano", "factor": 1.25, "floor": 0.1, "scan": 3},
    "scan2": {"endpoint": "rano", "factor": 1.4, "floor": 0.1, "scan": 2},
    "scan4": {"endpoint": "rano", "factor": 1.4, "floor": 0.1, "scan": 4},
    "old_rule_1.2_baseline": {"endpoint": "baseline", "factor": 1.4, "floor": 0.1, "scan": None},
}
STAT_COLS = ("median_interp", "median_step", "S84", "S168", "S252", "event_fraction")
WEAKEST = [
    "one published number (a median with its interval) fixes one condition per world; many worlds pass and the anchor does not say which is true",
    "kill strength is not anchored (no volume change per TMZ cycle was found), so the m-scale axis is free by construction",
    "the RANO-like rule (1.4 x the smallest earlier scan, floor 0.1 x baseline) is an assumed translation of the trial rule; the sensitivity variants show how much it matters",
    "the progression call uses the true burden, with no measurement error in the volume, so false calls and missed calls are not modelled",
    "the model has no patients who never respond to TMZ and no early deaths or drop-outs; the early part of the PFS curve is not checked",
    "world classification uses 400 patients per world, so a world near the edge of the interval can flip; script 193 re-draws the kept worlds",
]


def _smoke(args):
    global N_WORLDS, N_PAT, CHUNK_WORLDS, OLD_MODEL_N, OUT_DIR
    if args.smoke:
        N_WORLDS, N_PAT, CHUNK_WORLDS, OLD_MODEL_N = 12, 100, 6, 300
        OUT_DIR = PROJECT_ROOT / "output" / "trackc_stage2_C_calibration_smoke"


def run_variant(batch, v, chunk_key):
    """Standard arm for one variant on a batch of patients; restores the RANO constants afterwards."""
    old = (tc.RANO_FACTOR, tc.RANO_FLOOR)
    tc.RANO_FACTOR, tc.RANO_FLOOR = v["factor"], v["floor"]
    try:
        return cal.run_arm(batch, scan_every=v["scan"], endpoint=v["endpoint"], seed=chunk_key)
    finally:
        tc.RANO_FACTOR, tc.RANO_FLOOR = old


def _job(args):
    chunk, worlds, n_pat = args
    batch = cal.concat([cal.draw_patients(w, n_pat, [PATIENT_SEED, w["world"]]) for w in worlds])
    rows = [dict(w) for w in worlds]
    for vname, v in VARIANTS.items():
        ttp, ev, courses = run_variant(batch, v, [NOISE_BASE, chunk])
        for k, row in enumerate(rows):
            sl = slice(k * n_pat, (k + 1) * n_pat)
            st = cal.pfs_stats(ttp[sl], ev[sl])
            row[vname] = {c: st[c] for c in STAT_COLS} | {"mean_courses": float(courses[sl].mean())}
    return rows


def old_model_row(n):
    """The earlier world (script 150 ranges, c uniform 0 to 0.25, no switching) under every variant."""
    scen = tc.draw_scenarios(n, [PATIENT_SEED, 999_999])
    out = {}
    for vname, v in VARIANTS.items():
        ttp, ev, courses = run_variant(scen, v, [NOISE_BASE, 999_999])
        st = cal.pfs_stats(ttp, ev)
        out[vname] = {c: st[c] for c in STAT_COLS} | {"mean_courses": float(courses.mean())}
    return out


def _finite(x):
    x = np.asarray(x, float)
    return x[np.isfinite(x)]


def analyse(rows, old_model=None):
    """rows: list of world dicts as returned by _job. Pure function, no I/O."""
    n = len(rows)
    ok = {v: np.array([cal.in_target(r[v]["median_interp"]) for r in rows]) for v in VARIANTS}
    med = np.array([r["base"]["median_interp"] for r in rows])
    res = {"n_worlds": n, "target_days": list(cal.TARGET_PFS_DAYS),
           "base": {"n_pass": int(ok["base"].sum()), "frac_pass": float(ok["base"].mean()),
                    "median_interp_quantiles_5_25_50_75_95": [float(q) for q in np.percentile(_finite(med), [5, 25, 50, 75, 95])],
                    "n_median_never_reached": int((~np.isfinite(med)).sum()),
                    "frac_below_target": float((med < cal.TARGET_PFS_DAYS[0]).mean()),
                    "frac_above_target": float((med > cal.TARGET_PFS_DAYS[1]).mean())},
           "variants": {}, "axes": {}}
    for v in VARIANTS:
        both, either = (ok[v] & ok["base"]).sum(), (ok[v] | ok["base"]).sum()
        mv = np.array([r[v]["median_interp"] for r in rows])
        res["variants"][v] = {"n_pass": int(ok[v].sum()), "frac_pass": float(ok[v].mean()),
                              "jaccard_with_base_pass": float(both / either) if either else None,
                              "median_interp_quantiles_5_25_50_75_95": [float(q) for q in np.percentile(_finite(mv), [5, 25, 50, 75, 95])] if np.isfinite(mv).any() else None,
                              "frac_below_target": float((mv < cal.TARGET_PFS_DAYS[0]).mean()),
                              "frac_above_target": float((mv > cal.TARGET_PFS_DAYS[1]).mean())}
    for name in cal.WORLD_KEYS:
        x = np.array([r[name] for r in rows], float)
        entry = {"spearman_with_median": None, "pass_fraction_by_bin": {}}
        good = np.isfinite(med)
        if good.sum() > 3 and np.ptp(x[good]) > 0:
            rho = spearmanr(x[good], med[good]).statistic
            entry["spearman_with_median"] = float(rho)
        if name in ("mu", "nu"):
            bins = {f"{lvl:g}": x == lvl for lvl in cal.MU_NU_LEVELS}
        else:
            q = np.percentile(x, [100 / 3, 200 / 3])
            bins = {"low_third": x <= q[0], "middle_third": (x > q[0]) & (x <= q[1]), "high_third": x > q[1]}
        for b, m in bins.items():
            entry["pass_fraction_by_bin"][b] = {"n": int(m.sum()), "frac_pass": float(ok["base"][m].mean()) if m.any() else None}
        kept = x[ok["base"]]
        entry["kept_quartiles_25_50_75"] = [float(q) for q in np.percentile(kept, [25, 50, 75])] if kept.size else None
        entry["all_quartiles_25_50_75"] = [float(q) for q in np.percentile(x, [25, 50, 75])]
        res["axes"][name] = entry
    if old_model is not None:
        res["old_model"] = old_model
    return res


def write_csvs(rows, out_dir):
    cols = list(cal.WORLD_KEYS)
    header = ["world"] + cols + [f"{v}_{c}" for v in VARIANTS for c in STAT_COLS] + [f"{v}_pass" for v in VARIANTS]
    with open(out_dir / "worlds.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow([r["world"]] + [f"{r[c]:.10g}" for c in cols]
                       + [f"{r[v][c]:.6g}" for v in VARIANTS for c in STAT_COLS]
                       + [int(cal.in_target(r[v]["median_interp"])) for v in VARIANTS])
    with open(out_dir / "kept_worlds.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["world"] + cols + ["base_median_interp"])
        for r in rows:
            if cal.in_target(r["base"]["median_interp"]):
                w.writerow([r["world"]] + [f"{r[c]:.10g}" for c in cols] + [f"{r['base']['median_interp']:.6g}"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    _smoke(args)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    protocol = {"script": "192_trackc_calibration", "stage": "2C", "smoke": bool(args.smoke), "n_worlds": N_WORLDS,
                "n_patients_per_world": N_PAT, "chunk_worlds": CHUNK_WORLDS, "old_model_n": OLD_MODEL_N,
                "world_seed": WORLD_SEED, "patient_seed": PATIENT_SEED, "noise_base": NOISE_BASE,
                "axes": [list(a) for a in cal.AXES], "mu_nu_levels": list(cal.MU_NU_LEVELS),
                "b0_half_width": cal.B0_HALF_WIDTH, "f0_cap": cal.F0_CAP, "std_courses": cal.STD_COURSES,
                "target_days": list(cal.TARGET_PFS_DAYS), "variants": VARIANTS, "stat": "interpolated KM median"}
    manifest = man.build_manifest("192_trackc_calibration", protocol, purpose="stage2_dev",
                                  scenario=(WORLD_SEED, PATIENT_SEED), noise=(NOISE_BASE, NOISE_BASE),
                                  script=SCRIPT, require_clean=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; "
          f"git clean {manifest['git']['clean']}", flush=True)

    t0 = time.time()
    worlds = cal.draw_worlds(N_WORLDS, WORLD_SEED)
    jobs = [(c, worlds[i:i + CHUNK_WORLDS], N_PAT) for c, i in enumerate(range(0, N_WORLDS, CHUNK_WORLDS))]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = [row for part in ex.map(_job, jobs) for row in part]
    rows.sort(key=lambda r: r["world"])
    old = old_model_row(OLD_MODEL_N)
    write_csvs(rows, OUT_DIR)

    res = analyse(rows, old)
    res.update({"script": "192_trackc_calibration",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; calibration of the standard arm only",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    b = res["base"]
    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    print(f"base: {b['n_pass']} of {res['n_worlds']} worlds pass ({100 * b['frac_pass']:.1f}%); "
          f"median-of-medians {b['median_interp_quantiles_5_25_50_75_95'][2]:.0f} d; "
          f"{100 * b['frac_below_target']:.0f}% below, {100 * b['frac_above_target']:.0f}% above the interval")
    print(f"old model (base rule): median PFS {old['base']['median_interp']:.0f} d (interpolated), "
          f"{old['base']['median_step']:.0f} d (step)")
    for v, e in res["variants"].items():
        print(f"  {v:24s} pass {e['n_pass']:4d}  jaccard with base {e['jaccard_with_base_pass']}")


if __name__ == "__main__":
    main()
