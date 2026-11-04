#!/usr/bin/env python3
"""
Script 194: Track C Stage 2, test E: policy contrasts across the controlled worlds (regime map, first pass).
===========================================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Plan: paper/trackc_restart_plan_2026-11-03.md,
design: paper/trackc_anchored_spec_2026-11-04.md section 4 steps 4 and 5.

Question. In the worlds that match the published standard arm (script 192) and pass both negative controls
(script 193), how much can any dosing rule gain over the standard six courses, and which world axis decides it?

Arms (all scan every 84 days; rules that read scans see one noisy reading per scan, noise sd 0.05, 0.10 or 0.20):
    soc6                  six courses, one per cycle (the reference)
    front3, front9, front12   front-loaded k courses (the dose-response curve)
    paced2, paced3        a course every 2nd or 3rd cycle, up to 12 courses
    raw_zhang_1.0_0.5, raw_0.6_0.6   on/off rules on the raw reading against the first reading
    nadir_k3_on1.25, nadir_k3_on1.15, nadir_k6_on1.25   induction courses, then treat while the scan shows the burden
                          at least a_on times the smallest earlier reading (src/trackc_policies.py)
    best_fixed            the best open-loop schedule of a family of 56, chosen on separate patients
    informed_family_max   per patient, the best schedule of that family (a controller that knows the patient's
                          parameters; a family maximum, not a proven optimum, and chosen for 730-day restricted time;
                          ties go to fewer courses)
    best_fixed_lamL, informed_lamL   the same two arms, chosen to maximise restricted time minus L x courses, with
                          L = 2, 5, 10, 20 days per course (an assumed grid; the cost of a course is not measured)
Endpoints: 3b = RANO-like (1.4 x the smallest earlier scan, floor 0.1 x baseline, called at scans), the endpoint
the trials used; 3a = the old rule, true burden reaches 1.2 x baseline (the adaptive-therapy convention). Both are
reported for every arm because they can disagree.
Scores: restricted mean time to progression at 365 and 730 days, paired against soc6 on the same patients, the
paired-regret tail (mean of the worst 10% paired differences), and courses used.

Nothing is tuned on a sealed seed. The rules are simple on purpose; the particle-filter rules of scripts 154-173 are
NOT run here (cost), so this test says nothing about them.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev), configuration "single" (the run of 2026-11-04):
    evaluation patients   seed [20261804, world index], n = 1000 per world
    selection patients    seed [20261805, world index], n = 400 per world (best_fixed is chosen here)
    noise                 key [1_830_000, world index], the same in every arm
Written to output/trackc_stage2_E_regime_map/ (manifest.json first, then arms_by_world.csv, results.json).
--config single_rep | mixture reads the matching folder of script 193 and uses the seeds and folder suffix defined in
src/trackc_calibration.py.

Run:  python src/194_trackc_regime_map.py --workers 4
      python src/194_trackc_regime_map.py --config mixture --workers 4
      python src/194_trackc_regime_map.py --smoke        (reads the smoke output of script 193)
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
from scipy.stats import spearmanr

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_calibration as cal  # noqa: E402
import trackc_manifest as man  # noqa: E402
import trackc_regime as rg  # noqa: E402

SCRIPT = "src/194_trackc_regime_map.py"
N_EVAL, N_SEL = 1000, 400
LAMBDAS = (0.0, 2.0, 5.0, 10.0, 20.0)      # days of restricted time that one course is worth in the cost-aware arms (ASSUMED grid)
ENDPOINTS = {"3b_rano": "rano", "3a_baseline": "baseline"}
REFERENCE = "soc6"
WEAKEST = [
    "the controlled worlds are a thin slice of an assumed family that passed one published median and two published hazard ratios; the map describes the model, not tumours",
    "only simple rules are run; the particle-filter rules of scripts 154-173 are not, so nothing here says whether they would gain",
    "best_fixed and informed_family_max come from a family of 56 schedules; the maximum over a family is a lower bound on the true optimum",
    "the two endpoints disagree in places; 3b rests on an assumed translation of the RANO rule (factor 1.4, floor 0.1 x baseline) and uses the true burden without measurement error",
    "each world has 1,000 evaluation patients, so a per-world difference carries a sampling error of a few days",
    "nothing here is a clinical or patient-specific claim",
]
AXES = list(cal.WORLD_KEYS)                # default world columns ('single' family); other families pass their own


def _smoke(args):
    global N_EVAL, N_SEL
    if args.smoke:
        N_EVAL, N_SEL = 100, 60


def read_controlled(d_dir):
    with open(d_dir / "controlled_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


def evaluate_world(w, n_eval, n_sel, cfg):
    """All arms and both endpoints for one world. Pure function; returns a list of long-format rows."""
    world = int(w["world"])
    ev = cal.draw_patients(w, n_eval, [cfg["e_eval_seed"], world])
    sel = cal.draw_patients(w, n_sel, [cfg["e_sel_seed"], world])
    sched, names = rg.candidate_schedules()
    arms = rg.policy_arms()
    rows = []
    for ep_label, endpoint in ENDPOINTS.items():
        res = {name: rg.run_policy(ev, factory, budget, endpoint, [cfg["noise_e"], world])
               for name, (factory, budget) in arms.items()}
        ttp_sel, courses_sel = rg.evaluate_schedules(sel, sched, endpoint)
        ttp_ev, courses_ev = rg.evaluate_schedules(ev, sched, endpoint)
        cols = np.arange(ttp_ev.shape[1])
        notes = {}
        # lambda = days of restricted time one course is worth; 0 keeps the pure time score (ties go to fewer courses)
        for lam in LAMBDAS:
            tag = "" if lam == 0 else f"_lam{lam:g}"
            sel_score = [rg.rmst(t, 730.0) - lam * float(np.mean(c)) - 1e-6 * float(np.mean(c))
                         for t, c in zip(ttp_sel, courses_sel)]
            best = int(np.argmax(sel_score))
            res[f"best_fixed{tag}"] = (ttp_ev[best], courses_ev[best])
            notes[f"best_fixed{tag}"] = names[best]
            score = np.minimum(ttp_ev, 730.0) - lam * courses_ev - 1e-6 * courses_ev
            pick = score.argmax(0)
            res["informed_family_max" if lam == 0 else f"informed{tag}"] = (ttp_ev[pick, cols], courses_ev[pick, cols])
        ref_ttp, ref_courses = res[REFERENCE]
        inf_capped = np.minimum(res["informed_family_max"][0], 730.0)
        for name, (ttp, courses) in res.items():
            # share of patients on whom this arm beats the informed family maximum (it should be near 0 if that
            # maximum were a bound; a positive share means the family of 56 schedules is not exhaustive)
            beats = float(np.mean(np.minimum(ttp, 730.0) > inf_capped + 1e-6))
            for tau in rg.TAUS:
                st = rg.paired_stats(ttp, ref_ttp, tau)
                rows.append({"world": world, "endpoint": ep_label, "arm": name, "tau": int(tau),
                             "rmst": rg.rmst(ttp, tau), "diff_vs_soc6": st["diff"], "ci_lo": st["ci_lo"],
                             "ci_hi": st["ci_hi"], "regret_cvar10": st["regret_cvar"], "frac_worse": st["frac_worse"],
                             "mean_courses": float(np.mean(courses)),
                             "courses_diff_vs_soc6": float(np.mean(courses) - np.mean(ref_courses)),
                             "frac_beats_informed_730": beats,
                             "note": notes.get(name, "")})
    return rows


def _job(args):
    return evaluate_world(*args)


def _ols_standardised(X, y, names=None):
    """Standardised least-squares coefficients and R-squared; columns with no spread are dropped."""
    names = AXES if names is None else names
    keep = [i for i in range(X.shape[1]) if X[:, i].std() > 0]
    Z = (X[:, keep] - X[:, keep].mean(0)) / X[:, keep].std(0)
    yc = y - y.mean()
    beta, *_ = np.linalg.lstsq(Z, yc, rcond=None)
    ss = float(((yc - Z @ beta) ** 2).sum())
    r2 = 1.0 - ss / float((yc ** 2).sum()) if (yc ** 2).sum() > 0 else float("nan")
    return {names[k]: float(b) for k, b in zip(keep, beta)}, r2


def analyse(rows, worlds, keys=None):
    """Summary over worlds. Pure function, no I/O. keys: the world-column names of the family (default 'single')."""
    axes = AXES if keys is None else list(keys)
    w_by = {w["world"]: w for w in worlds}
    res = {"n_controlled_worlds": len(worlds), "reference": REFERENCE, "endpoints": {}}
    arms = sorted({r["arm"] for r in rows})
    for ep in ENDPOINTS:
        res["endpoints"][ep] = {}
        for tau in rg.TAUS:
            tau = int(tau)
            block = {}
            for arm in arms:
                sel = [r for r in rows if r["endpoint"] == ep and r["arm"] == arm and r["tau"] == tau]
                d = np.array([r["diff_vs_soc6"] for r in sel])
                if d.size == 0:
                    continue
                block[arm] = {"n_worlds": int(d.size), "median_diff_days": float(np.median(d)),
                              "mean_diff_days": float(d.mean()),
                              "q25_q75": [float(q) for q in np.percentile(d, [25, 75])],
                              "min_max": [float(d.min()), float(d.max())],
                              "frac_worlds_clear_gain": float(np.mean([r["ci_lo"] > 0 for r in sel])),
                              "frac_worlds_clear_loss": float(np.mean([r["ci_hi"] < 0 for r in sel])),
                              "median_regret_cvar10": float(np.median([r["regret_cvar10"] for r in sel])),
                              "median_courses_diff": float(np.median([r["courses_diff_vs_soc6"] for r in sel])),
                              "median_frac_patients_beating_informed_family_max": float(np.median([r["frac_beats_informed_730"] for r in sel])),
                              "median_mean_courses": float(np.median([r["mean_courses"] for r in sel]))}
            res["endpoints"][ep][str(tau)] = {"arms": block}
            # which axis decides the gain of the two upper-bound-style arms and of the best simple rule
            realizable = [a for a in block if a != "soc6" and not a.startswith(("best_fixed", "informed"))]
            best_rule = max(realizable, key=lambda a: block[a]["median_diff_days"]) if realizable else None
            # cost-aware view: restricted time minus lambda x courses, relative to soc6, median over worlds
            net = {}
            for lam in LAMBDAS:
                net[f"{lam:g}"] = {a: float(np.median([r["diff_vs_soc6"] - lam * r["courses_diff_vs_soc6"]
                                                      for r in rows if r["endpoint"] == ep and r["arm"] == a and r["tau"] == tau]))
                                   for a in block}
            axes_block = {}
            for arm in ("best_fixed", "informed_family_max", best_rule):
                if arm is None or len(worlds) < 4:
                    continue
                sel = sorted([r for r in rows if r["endpoint"] == ep and r["arm"] == arm and r["tau"] == tau],
                             key=lambda r: r["world"])
                y = np.array([r["diff_vs_soc6"] for r in sel])
                X = np.array([[w_by[r["world"]][a] for a in axes] for r in sel], float)
                coef, r2 = _ols_standardised(X, y, axes)
                sp = {}
                for i, a in enumerate(axes):
                    if X[:, i].std() > 0 and y.std() > 0:
                        sp[a] = float(spearmanr(X[:, i], y).statistic)
                axes_block[arm] = {"spearman": sp, "standardised_ols_coef": coef, "ols_r2": r2}
            res["endpoints"][ep][str(tau)]["net_score_median_by_lambda"] = net
            res["endpoints"][ep][str(tau)]["best_simple_rule_by_median_diff"] = best_rule
            res["endpoints"][ep][str(tau)]["what_decides_the_gain"] = axes_block
    return res


def write_csv(rows, out_dir):
    cols = list(rows[0].keys())
    with open(out_dir / "arms_by_world.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else r[c] for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py (single, single_rep, mixture)")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    OUT_DIR, D_DIR = cal.stage_dir("E", cfg, args.smoke), cal.stage_dir("D", cfg, args.smoke)
    keys = cal.world_keys(cfg["family"])
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (D_DIR / "results.json").exists():
        sys.exit(f"script 193 has not finished: {D_DIR / 'results.json'} is missing")
    worlds = read_controlled(D_DIR)
    if not worlds:
        sys.exit("script 193 found no controlled world; see its kill rule")
    sched, names = rg.candidate_schedules()
    protocol = {"script": "194_trackc_regime_map", "stage": "2E", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_eval": N_EVAL, "n_sel": N_SEL,
                "eval_seed": cfg["e_eval_seed"], "sel_seed": cfg["e_sel_seed"], "noise_base": cfg["noise_e"],
                "scan_every": rg.SCAN_EVERY,
                "taus": list(rg.TAUS), "endpoints": ENDPOINTS, "reference": REFERENCE, "lambdas": list(LAMBDAS),
                "arms": {k: v[1] for k, v in rg.policy_arms().items()}, "candidate_schedules": names,
                "controlled_worlds": [w["world"] for w in worlds],
                "input": str((D_DIR / "controlled_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((D_DIR / "results.json").read_bytes()).hexdigest(),
                "rano_factor": cal.tc.RANO_FACTOR, "rano_floor": cal.tc.RANO_FLOOR}
    manifest = man.build_manifest("194_trackc_regime_map", protocol, purpose="stage2_dev",
                                  scenario=(cfg["e_eval_seed"], cfg["e_sel_seed"]),
                                  noise=(cfg["noise_e"], cfg["noise_e"]),
                                  script=SCRIPT, require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py",
                                                         "src/trackc_policies.py", "src/trackc_regime.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(worlds)} controlled worlds", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = list(ex.map(_job, [(w, N_EVAL, N_SEL, cfg) for w in worlds]))
    rows = [r for part in parts for r in part]
    write_csv(rows, OUT_DIR)
    res = analyse(rows, worlds, keys)
    res.update({"script": "194_trackc_regime_map",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; simple rules only, no particle filter",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for ep in ENDPOINTS:
        for tau in (365, 730):
            blk = res["endpoints"][ep][str(tau)]["arms"]
            print(f"== {ep}, restricted mean {tau} d, median difference to soc6 over {len(worlds)} worlds (days)")
            for arm, e in sorted(blk.items(), key=lambda kv: -kv[1]["median_diff_days"]):
                print(f"   {arm:22s} {e['median_diff_days']:+7.1f}  [{e['min_max'][0]:+.0f}, {e['min_max'][1]:+.0f}]  "
                      f"courses {e['median_mean_courses']:.1f}  clear gain in {100 * e['frac_worlds_clear_gain']:.0f}% of worlds")


if __name__ == "__main__":
    main()
