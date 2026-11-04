#!/usr/bin/env python3
"""
Script 199: Track C Stage 2, prediction check: how deep a response does a controlled world predict?
===================================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Plan: paper/trackc_restart_plan_2026-11-03.md.

Question. Test C found that matching the published median PFS of the standard arm needs a weak kill strength (the
median kept world has m about 0.8 to 2.7, against 6 to 20 in the old ranges). That implies few deep responses. The
calibration did not use response depth, so published response rates and volume changes are an independent check.
This script only computes the model's prediction; the published numbers are compared in the report once their sources
are verified (paper/trackc_stage2_report_2026-11-04.md).

Method. For each controlled world of script 193, draw 1,000 fresh patients and give the standard arm (six courses, one
per cycle, open loop). Read the true burden at the scans (every 84 days, days 84 to 336, the first four scans after
baseline) and record the best response: the smallest ratio of burden to baseline. Report, per world and over worlds,
the share with a best response of at most 0.35 (a 65% volume reduction, the volumetric-RANO partial-response line), at
most 0.5 (a 50% reduction), at most 0.75, and the median best ratio. The same is computed for the earlier ranges (old
model) for contrast.
It uses the true burden, so measurement error in real volumes is not modelled.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev):
    patients   seed [20261807, world index], n = 1000 per world; the old model uses [20261807, 999999]
    noise      none (open loop, no readings)
Written to output/trackc_stage2_J_response_depth/ (manifest.json first, then depth_by_world.csv, results.json).

The seed above is that of configuration "single"; --config single_rep | mixture takes the seed and folder suffix
defined in src/trackc_calibration.py and reads the matching folder of script 193.

Run:  python src/199_trackc_response_depth.py --workers 4
      python src/199_trackc_response_depth.py --smoke        (reads the smoke output of script 193)
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import hashlib
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_manifest as man  # noqa: E402
import trackc_sim as tc  # noqa: E402

SCRIPT = "src/199_trackc_response_depth.py"
N_PAT = 1000
SCANS_USED = 4                                         # scans at cycles 3, 6, 9, 12 (days 84, 168, 252, 336)
WEAKEST = [
    "the response is read from the true burden at the scan, with no measurement error and without the minimum measurable lesion size, so the model predicts more small responses than a clinic could record",
    "the model has no patients without measurable disease at the start; real response rates are computed among patients with measurable disease and the share differs between studies",
    "published response criteria (Macdonald, RANO) use area or volume thresholds that are not the same as a burden ratio, so the comparison is approximate",
    "one set of controlled worlds on development seeds; the share of deep responders varies a lot between worlds",
]


def _smoke(args):
    global N_PAT
    if args.smoke:
        N_PAT = 200


def read_controlled(d_dir):
    with open(d_dir / "controlled_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


def best_response(scen, courses=cal.STD_COURSES, scans=SCANS_USED, scan_every=cal.SCAN_EVERY):
    """Smallest ratio of true burden to baseline over the first `scans` scans, standard arm, open loop.
    courses: an int, or one course cap per patient (see cal.std_caps)."""
    S = (1 - scen["f0"]) * scen["B0"]
    R = scen["f0"] * scen["B0"]
    base = S + R
    best = np.ones_like(base)
    for j in range(scans * scan_every):
        treat = np.where(j < np.asarray(courses), 1.0, 0.0) * np.ones(base.shape)
        S, R = tc.run_cycle(S, R, scen, treat, j * tc.CYCLE_DAYS)
        if (j + 1) % scan_every == 0:                  # the scan at the start of cycle j + 1
            best = np.minimum(best, (S + R) / base)
    return best


def summarise(ratio):
    return {"frac_best_le_0.35": float(np.mean(ratio <= 0.35)), "frac_best_le_0.5": float(np.mean(ratio <= 0.5)),
            "frac_best_le_0.75": float(np.mean(ratio <= 0.75)),
            "frac_any_shrinkage": float(np.mean(ratio < 1.0)), "median_best_ratio": float(np.median(ratio))}


def evaluate_world(w, n_pat, seed, std_courses=cal.STD_COURSES):
    ratio = best_response(cal.draw_patients(w, n_pat, [seed, int(w["world"])]), courses=cal.std_caps(std_courses, n_pat))
    return {"world": int(w["world"]), **summarise(ratio)}


def _job(args):
    return evaluate_world(*args)


def analyse(rows, old):
    out = {"n_controlled_worlds": len(rows), "scans_used_days": [tc.CYCLE_DAYS * cal.SCAN_EVERY * (i + 1) for i in range(SCANS_USED)],
           "old_model": old}
    for key in ("frac_best_le_0.35", "frac_best_le_0.5", "frac_best_le_0.75", "frac_any_shrinkage", "median_best_ratio"):
        x = np.array([r[key] for r in rows])
        out[key] = {"mean_over_worlds": float(x.mean()), "q05_q25_q50_q75_q95": [float(q) for q in np.percentile(x, [5, 25, 50, 75, 95])],
                    "min_max": [float(x.min()), float(x.max())]}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py (single, single_rep, mixture, single_mix, mixture_mix)")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    OUT_DIR, D_DIR = cal.stage_dir("J", cfg, args.smoke), cal.stage_dir("D", cfg, args.smoke)
    seed = cfg["j_seed"]
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (D_DIR / "results.json").exists():
        sys.exit(f"script 193 has not finished: {D_DIR / 'results.json'} is missing")
    worlds = read_controlled(D_DIR)
    if not worlds:
        sys.exit("script 193 found no controlled world")
    protocol = {"script": "199_trackc_response_depth", "stage": "2J", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_patients_per_world": N_PAT,
                "scenario_seed": seed, "scans_used": SCANS_USED, **cal.exposure_protocol(cfg, always=True),
                "controlled_worlds": [w["world"] for w in worlds],
                "input": str((D_DIR / "controlled_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((D_DIR / "results.json").read_bytes()).hexdigest()}
    manifest = man.build_manifest("199_trackc_response_depth", protocol, purpose="stage2_dev",
                                  scenario=(seed, seed), noise=None,
                                  script=SCRIPT, require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(worlds)} controlled worlds", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(_job, [(w, N_PAT, seed, cfg["std_courses"]) for w in worlds]))
    old = summarise(best_response(tc.draw_scenarios(2 * N_PAT, [seed, 999_999]),
                                  courses=cal.std_caps(cfg["std_courses"], 2 * N_PAT)))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with open(OUT_DIR / "depth_by_world.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(list(rows[0].keys()))
        for r in rows:
            w.writerow([f"{v:.6g}" if isinstance(v, float) else v for v in r.values()])
    res = analyse(rows, old)
    res.update({"script": "199_trackc_response_depth",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; prediction check, not a fit",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for key in ("frac_best_le_0.35", "frac_best_le_0.5", "frac_best_le_0.75", "median_best_ratio"):
        e = res[key]
        print(f"{key:20s} mean over worlds {e['mean_over_worlds']:.3f}  q05-q95 {e['q05_q25_q50_q75_q95'][0]:.3f}-{e['q05_q25_q50_q75_q95'][4]:.3f}")
    print("old model:", {k: round(v, 3) for k, v in old.items()})


if __name__ == "__main__":
    main()
