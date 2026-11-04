#!/usr/bin/env python3
"""
Script 198: Track C Stage 2, test I: do the calibration and the negative controls depend on the integration step?
==================================================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Plan: paper/trackc_restart_plan_2026-11-03.md.

Question. Scripts 192 to 194 integrate the equations with a Runge-Kutta step of 0.25 days (tc.DT). The RANO-like call
reads the burden at scan times, so a small integration error could move a patient across the 1.4 x nadir line at one
scan and change the called time by a whole scan interval. Do the standard-arm median, the two negative-control hazard
ratios and their pass flags change when the step is 0.10 days?

Method. For each controlled world of script 193, draw 1,000 fresh patients (new seed), run the same three open-loop
quantities at DT = 0.25 and at DT = 0.10 on the same patients: the interpolated median PFS of the standard arm, the
dose-dense hazard ratio (multiplier 2.0), and the 12-versus-6 landmark hazard ratio. Report the share of patients whose
called time changes, the largest changes in the median and in the hazard ratios, and every flag that flips.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev):
    patients   seed [20261806, world index], n = 1000 per world
    noise      key [1_840_000, world index] (the arms ignore readings; the key is only passed to the simulator)
Written to output/trackc_stage2_I_numerics/ (manifest.json first, then numerics_by_world.csv, results.json).

Seeds above are those of configuration "single"; --config single_rep | mixture takes the seeds and folder suffix
defined in src/trackc_calibration.py and reads the matching folder of script 193.

Run:  python src/198_trackc_numerics.py --workers 4
      python src/198_trackc_numerics.py --smoke        (reads the smoke output of script 193)
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

SCRIPT = "src/198_trackc_numerics.py"
N_PAT = 1000
DTS = (0.25, 0.10)
WEAKEST = [
    "two step sizes only; agreement between them does not prove convergence, although the earlier check (script 159) found a relative error below 1e-4 for the fixed-step solver",
    "1,000 patients per world; a world near an edge of a published interval can flip for sampling reasons that this test cannot separate from numerical ones, so flips are reported with the distance to the edge",
    "the equations and the RANO-like rule are unchanged; this test says nothing about whether they describe tumours",
]


def _smoke(args):
    global N_PAT
    if args.smoke:
        N_PAT = 200


def read_controlled(d_dir):
    with open(d_dir / "controlled_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


def measure(w, n_pat, dt, cfg):
    """The calibration and control quantities for one world at step `dt`. The module step is restored afterwards."""
    old = tc.DT
    tc.DT = dt
    try:
        scen = cal.draw_patients(w, n_pat, [cfg["i_seed"], int(w["world"])])
        seed = [cfg["noise_i"], int(w["world"])]
        courses = cal.std_caps(cfg["std_courses"], n_pat)
        ref = cal.run_arm(scen, courses=courses, seed=seed)
        st = cal.pfs_stats(ref[0], ref[1])
        dd = cal.dose_dense_control(scen, cal.DOSE_DENSE_MULT, seed=seed, ref=ref, courses=courses)
        dd_ttp = cal.run_arm(scen, courses=courses, dose=cal.DOSE_DENSE_MULT, seed=seed)[0]
        lm = cal.landmark_control(scen, seed=seed, ref=ref, courses=courses)
    finally:
        tc.DT = old
    return {"median": st["median_interp"], "dd_hr": dd["hr"], "lm_hr": lm["hr"], "ttp_std": ref[0], "ttp_dd": dd_ttp}


def evaluate_world(w, n_pat, cfg):
    a, b = (measure(w, n_pat, dt, cfg) for dt in DTS)
    flags = lambda m: (cal.in_target(m["median"]),  # noqa: E731
                       cal.DOSE_DENSE_HR[1] <= m["dd_hr"] <= cal.DOSE_DENSE_HR[2],
                       cal.TWELVE_VS_SIX_HR[1] <= m["lm_hr"] <= cal.TWELVE_VS_SIX_HR[2])
    fa, fb = flags(a), flags(b)
    return {"world": int(w["world"]),
            "median_dt0.25": a["median"], "median_dt0.10": b["median"],
            "dd_hr_dt0.25": a["dd_hr"], "dd_hr_dt0.10": b["dd_hr"],
            "lm_hr_dt0.25": a["lm_hr"], "lm_hr_dt0.10": b["lm_hr"],
            "frac_patients_changed_std": float(np.mean(a["ttp_std"] != b["ttp_std"])),
            "frac_patients_changed_dd": float(np.mean(a["ttp_dd"] != b["ttp_dd"])),
            "flip_target": fa[0] != fb[0], "flip_dd": fa[1] != fb[1], "flip_lm": fa[2] != fb[2],
            "edge_distance_median": float(min(abs(b["median"] - cal.TARGET_PFS_DAYS[0]), abs(b["median"] - cal.TARGET_PFS_DAYS[1])))}


def _job(args):
    return evaluate_world(*args)


def analyse(rows):
    """Summary over worlds. Pure function, no I/O."""
    n = len(rows)
    d_med = np.array([abs(r["median_dt0.25"] - r["median_dt0.10"]) for r in rows])
    d_dd = np.array([abs(r["dd_hr_dt0.25"] - r["dd_hr_dt0.10"]) for r in rows])
    d_lm = np.array([abs(r["lm_hr_dt0.25"] - r["lm_hr_dt0.10"]) for r in rows])
    ch_std = np.array([r["frac_patients_changed_std"] for r in rows])
    ch_dd = np.array([r["frac_patients_changed_dd"] for r in rows])
    flips = [r for r in rows if r["flip_target"] or r["flip_dd"] or r["flip_lm"]]
    return {"n_worlds": n, "steps_days": list(DTS),
            "abs_change_median_days": {"mean": float(d_med.mean()), "max": float(d_med.max())},
            "abs_change_dose_dense_hr": {"mean": float(d_dd.mean()), "max": float(d_dd.max())},
            "abs_change_landmark_hr": {"mean": float(d_lm.mean()), "max": float(d_lm.max())},
            "frac_patients_with_changed_call_standard_arm": {"mean": float(ch_std.mean()), "max": float(ch_std.max())},
            "frac_patients_with_changed_call_dose_dense_arm": {"mean": float(ch_dd.mean()), "max": float(ch_dd.max())},
            "n_worlds_with_any_flag_flip": len(flips),
            "flips": [{"world": r["world"], "target": r["flip_target"], "dose_dense": r["flip_dd"], "landmark": r["flip_lm"],
                       "edge_distance_median_days": r["edge_distance_median"]} for r in flips]}


def write_csv(rows, out_dir):
    cols = list(rows[0].keys())
    with open(out_dir / "numerics_by_world.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else int(r[c]) if isinstance(r[c], (bool, np.bool_)) else r[c] for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py (single, single_rep, mixture, single_mix, mixture_mix)")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    OUT_DIR, D_DIR = cal.stage_dir("I", cfg, args.smoke), cal.stage_dir("D", cfg, args.smoke)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (D_DIR / "results.json").exists():
        sys.exit(f"script 193 has not finished: {D_DIR / 'results.json'} is missing")
    worlds = read_controlled(D_DIR)
    if not worlds:
        sys.exit("script 193 found no controlled world")
    protocol = {"script": "198_trackc_numerics", "stage": "2I", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_patients_per_world": N_PAT,
                "scenario_seed": cfg["i_seed"], "noise_base": cfg["noise_i"], "steps_days": list(DTS),
                **cal.exposure_protocol(cfg),
                "dose_dense_multiplier": cal.DOSE_DENSE_MULT, "controlled_worlds": [w["world"] for w in worlds],
                "input": str((D_DIR / "controlled_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((D_DIR / "results.json").read_bytes()).hexdigest(),
                "rano_factor": tc.RANO_FACTOR, "rano_floor": tc.RANO_FLOOR, "event_mode": tc.EVENT_MODE}
    manifest = man.build_manifest("198_trackc_numerics", protocol, purpose="stage2_dev",
                                  scenario=(cfg["i_seed"], cfg["i_seed"]), noise=(cfg["noise_i"], cfg["noise_i"]),
                                  script=SCRIPT, require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(worlds)} controlled worlds", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(_job, [(w, N_PAT, cfg) for w in worlds]))
    write_csv(rows, OUT_DIR)
    res = analyse(rows)
    res.update({"script": "198_trackc_numerics",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; integration-step check",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    print(f"largest change in median PFS {res['abs_change_median_days']['max']:.2f} d; in dose-dense HR "
          f"{res['abs_change_dose_dense_hr']['max']:.4f}; in landmark HR {res['abs_change_landmark_hr']['max']:.4f}")
    print(f"patients with a changed call: standard arm mean {100 * res['frac_patients_with_changed_call_standard_arm']['mean']:.2f}%, "
          f"max {100 * res['frac_patients_with_changed_call_standard_arm']['max']:.2f}%; worlds with a flag flip: {res['n_worlds_with_any_flag_flip']}")


if __name__ == "__main__":
    main()
