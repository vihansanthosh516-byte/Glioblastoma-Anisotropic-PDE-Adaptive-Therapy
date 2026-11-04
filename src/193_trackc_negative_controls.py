#!/usr/bin/env python3
"""
Script 193: Track C Stage 2, test D: do the kept worlds also pass the two negative controls?
===========================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Plan: paper/trackc_restart_plan_2026-11-03.md,
design: paper/trackc_anchored_spec_2026-11-04.md section 4 step 3.

Question. Script 192 found the worlds whose simulated standard arm matches the published median PFS. Real trials found
little or no gain from more drug: dose-dense TMZ against standard TMZ (RTOG 0525, PFS HR 0.87, 95% CI 0.75 to 1.00) and
12 against 6 cycles among patients still progression-free after cycle 6 (GEINO 14-01, PFS HR 0.99, 95% CI 0.71 to
1.38). A world in which more drug helps far more than that is not consistent with the trials, even if its standard
arm matches. Which of the kept worlds survive both controls?

Method. For each kept world of script 192, draw 3,000 NEW patients (fresh seed), re-check that the standard arm still
lands in 143 to 186 days, then run: (1) the standard arm with the drug effect multiplied by 2.0 (dose-dense; the
multiplier is an ASSUMED ratio of drug exposure per cycle, 1.5 and 3.0 are sensitivity runs), and (2) 12 against 6
courses among patients without a progression call at the day-168 scan, measured from that scan. Hazard ratios come
from a Cox model with Efron ties (src/trackc_survival.py). A control passes if the model hazard ratio lies inside
the published 95% interval. No dosing policy is compared here.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev), configuration "single" (the run of 2026-11-04):
    patients   seed [20261803, world index]
    noise      key [1_820_000, world index]   (the arms ignore readings; the key is only passed to the simulator)
Written to output/trackc_stage2_D_negative_controls/ (manifest.json first, then controls.csv, controlled_worlds.csv,
results.json). --config single_rep | mixture reads the matching folder of script 192 and uses the seeds and the folder
suffix defined in src/trackc_calibration.py.

Run:  python src/193_trackc_negative_controls.py --workers 3
      python src/193_trackc_negative_controls.py --config mixture --workers 3
      python src/193_trackc_negative_controls.py --smoke        (reads the smoke output of script 192)
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

SCRIPT = "src/193_trackc_negative_controls.py"
N_PAT = 3000
DD_MULTS = (1.5, 2.0, 3.0)
PRIMARY_MULT = cal.DOSE_DENSE_MULT
WEAKEST = [
    "the dose-dense multiplier is an assumed ratio of drug exposure per cycle; the label doses of the dose-dense arm were not verified, so 1.5 and 3.0 are reported next to the primary 2.0",
    "both controls compare open-loop arms; the model hazard ratio is computed on the same patients in both arms, so its Cox interval is not a valid confidence interval; the flags use the point estimate",
    "passing the controls is necessary, not sufficient: the controls say nothing about how a feedback rule behaves",
    "the landmark control uses 3,000 simulated patients per world, about 1,500 of them past the landmark; a world near an interval edge can flip",
    "the real trials include toxicity, drop-out and deaths before progression, none of which are modelled",
]


def _smoke(args):
    global N_PAT
    if args.smoke:
        N_PAT = 400


def read_kept(c_dir):
    with open(c_dir / "kept_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


def evaluate_world(w, n_pat, cfg, keys):
    """One kept world: re-draw patients, re-check the calibration, run both controls. Pure function."""
    scen = cal.draw_patients(w, n_pat, [cfg["d_seed"], int(w["world"])])
    seed = [cfg["noise_d"], int(w["world"])]
    ref = cal.run_arm(scen, seed=seed)
    st = cal.pfs_stats(ref[0], ref[1])
    row = {k: w[k] for k in ["world"] + list(keys)}
    row.update({"c_median_interp": w["base_median_interp"], "re_median_interp": st["median_interp"],
                "re_pass": bool(cal.in_target(st["median_interp"])), "re_S84": st["S84"], "re_S168": st["S168"]})
    for m in DD_MULTS:
        d = cal.dose_dense_control(scen, mult=m, seed=seed, ref=ref)
        tag = f"dd{m:g}"
        row[f"{tag}_hr"], row[f"{tag}_median"] = d["hr"], d["median_dose_dense"]
        row[f"{tag}_ok"] = bool(cal.DOSE_DENSE_HR[1] <= d["hr"] <= cal.DOSE_DENSE_HR[2])
    lm = cal.landmark_control(scen, seed=seed, ref=ref)
    row.update({"lm_hr": lm["hr"], "lm_ci_lo": lm["ci"][0], "lm_ci_hi": lm["ci"][1], "lm_n": lm["n_landmark"]})
    row["lm_ok"] = bool(cal.TWELVE_VS_SIX_HR[1] <= lm["hr"] <= cal.TWELVE_VS_SIX_HR[2])
    row["controlled"] = bool(row["re_pass"] and row[f"dd{PRIMARY_MULT:g}_ok"] and row["lm_ok"])
    return row


def _job(args):
    return evaluate_world(*args)


def analyse(rows, keys=None):
    """Summary of the control results. Pure function, no I/O."""
    keys = cal.WORLD_KEYS if keys is None else keys
    n = len(rows)
    res = {"n_kept_worlds_from_script_192": n, "primary_dose_multiplier": PRIMARY_MULT,
           "published": {"dose_dense_pfs_hr": list(cal.DOSE_DENSE_HR), "twelve_vs_six_pfs_hr": list(cal.TWELVE_VS_SIX_HR)}}
    if n == 0:
        res.update({"n_still_in_target": 0, "n_controlled": 0})
        return res
    flag = lambda key: np.array([bool(r[key]) for r in rows])  # noqa: E731
    p = f"dd{PRIMARY_MULT:g}"
    res["n_still_in_target"] = int(flag("re_pass").sum())
    res["n_dose_dense_ok"] = {f"dd{m:g}": int(flag(f"dd{m:g}_ok").sum()) for m in DD_MULTS}
    res["n_landmark_ok"] = int(flag("lm_ok").sum())
    res["n_both_controls_ok_primary"] = int((flag(f"{p}_ok") & flag("lm_ok")).sum())
    res["n_controlled"] = int(flag("controlled").sum())
    for key, label in ((f"{p}_hr", "dose_dense_hr_primary"), ("lm_hr", "landmark_hr")):
        x = np.array([r[key] for r in rows], float)
        res[label + "_quantiles_5_25_50_75_95"] = [float(q) for q in np.nanpercentile(x, [5, 25, 50, 75, 95])]
    res["frac_dose_dense_hr_below_published_lower_bound"] = float(np.mean([r[f"{p}_hr"] < cal.DOSE_DENSE_HR[1] for r in rows]))
    res["frac_landmark_hr_below_published_lower_bound"] = float(np.mean([r["lm_hr"] < cal.TWELVE_VS_SIX_HR[1] for r in rows]))
    res["axes"] = {}
    ctl = flag("controlled")
    for name in keys:
        x = np.array([r[name] for r in rows], float)
        entry = {"kept_mean": float(x.mean()), "controlled_mean": float(x[ctl].mean()) if ctl.any() else None,
                 "landmark_ok_mean": float(x[flag("lm_ok")].mean()) if flag("lm_ok").any() else None,
                 "landmark_fail_mean": float(x[~flag("lm_ok")].mean()) if (~flag("lm_ok")).any() else None}
        res["axes"][name] = entry
    if res["n_controlled"] == 0:
        res["kill_rule_fired"] = ("no kept world passes both controls: within the assumed world family, matching the standard arm "
                                  "forces a model benefit of more drug that the published trials exclude")
    return res


def write_csvs(rows, out_dir, keys=None):
    keys = list(cal.WORLD_KEYS if keys is None else keys)
    if not rows:
        return
    cols = list(rows[0].keys())
    with open(out_dir / "controls.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else int(r[c]) if isinstance(r[c], (bool, np.bool_)) else r[c] for c in cols])
    with open(out_dir / "controlled_worlds.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["world"] + keys + ["re_median_interp", f"dd{PRIMARY_MULT:g}_hr", "lm_hr"])
        for r in rows:
            if r["controlled"]:
                w.writerow([r["world"]] + [f"{r[c]:.10g}" for c in keys]
                           + [f"{r['re_median_interp']:.6g}", f"{r[f'dd{PRIMARY_MULT:g}_hr']:.4f}", f"{r['lm_hr']:.4f}"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py (single, single_rep, mixture)")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    OUT_DIR, C_DIR = cal.stage_dir("D", cfg, args.smoke), cal.stage_dir("C", cfg, args.smoke)
    keys = cal.world_keys(cfg["family"])
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (C_DIR / "results.json").exists():
        sys.exit(f"script 192 has not finished: {C_DIR / 'results.json'} is missing")
    kept = read_kept(C_DIR)
    protocol = {"script": "193_trackc_negative_controls", "stage": "2D", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_patients_per_world": N_PAT,
                "scenario_seed": cfg["d_seed"], "noise_base": cfg["noise_d"], "dose_multipliers": list(DD_MULTS),
                "primary_multiplier": PRIMARY_MULT, "landmark_days": cal.LANDMARK_DAYS,
                "published_dose_dense_hr": list(cal.DOSE_DENSE_HR), "published_12_vs_6_hr": list(cal.TWELVE_VS_SIX_HR),
                "target_days": list(cal.TARGET_PFS_DAYS), "kept_worlds": [k["world"] for k in kept],
                "input": str((C_DIR / "kept_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((C_DIR / "results.json").read_bytes()).hexdigest()}
    last = max([k["world"] for k in kept], default=0)
    manifest = man.build_manifest("193_trackc_negative_controls", protocol, purpose="stage2_dev",
                                  scenario=(cfg["d_seed"], cfg["d_seed"]), noise=(cfg["noise_d"], cfg["noise_d"]),
                                  script=SCRIPT, require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(kept)} kept worlds (largest index {last}); config {cfg['name']}", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(_job, [(w, N_PAT, cfg, keys) for w in kept]))
    write_csvs(rows, OUT_DIR, keys)
    res = analyse(rows, keys)
    res.update({"script": "193_trackc_negative_controls",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; negative controls on open-loop arms",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    print(f"kept worlds {res['n_kept_worlds_from_script_192']}; still in target {res.get('n_still_in_target')}; "
          f"dose-dense ok {res.get('n_dose_dense_ok')}; landmark ok {res.get('n_landmark_ok')}; controlled {res['n_controlled']}")
    if "kill_rule_fired" in res:
        print("KILL RULE:", res["kill_rule_fired"])


if __name__ == "__main__":
    main()
