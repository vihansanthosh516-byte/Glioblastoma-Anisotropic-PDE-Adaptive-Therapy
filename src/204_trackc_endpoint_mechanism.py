#!/usr/bin/env python3
"""
Script 204: Track C Stage 2, test M: what does the RANO-like headroom consist of?
================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. The design below was written and committed BEFORE any
result of this script existed. A pilot of the code on three worlds of the single family, 60 patients and throw-away seeds
(not the registered ones) was looked at once to check that the pipeline runs; it showed that no scan ever shows a fall to
below 0.9 of baseline in that family and that every arm is called 120 to 150 days after the burden first reaches 1.2 x
baseline. That observation is why the scan-interval sweep (part M2) and the 1.4 x baseline crossing (part M1) were added.

Question. Scripts 194 and 202 show that a controller that knows every patient's parameters (the "informed" arm: the best of
56 open-loop schedules for each patient) lies far above the best fixed schedule at the same mean number of courses when
progression is called by the RANO-like rule (+29.5 d in the single family, +40.9 d in the mixture family) but almost not at
all when progression is called by the baseline-referenced rule (+0.15 d and +5.1 d; ledger C21, C22). The RANO-like rule
compares each scan with the smallest earlier scan (floor 0.1 x baseline) and is checked only at scans, every 84 days. Is
that headroom better control of the disease, or a later call of a similar disease course? Two candidate causes are
looked at: (a) the scan grid (a burden path can cross the threshold between two scans and be called at the next scan, up
to 84 days later; a controller that knows the patient can place its courses so that the crossing falls between scans), and
(b) the nadir reference (a scanned fall lowers the reference, so a good response brings the call forward).

Method (open loop, so no measurement noise). For every controlled world of a cap-6 configuration (script 193) draw
evaluation patients and selection patients. Every one of the 56 candidate schedules of src/trackc_regime.py is run with the
RANO-like rule (a re-implementation of the RANO branch of src/trackc_sim.py that also keeps the true burden at every cycle
boundary, days 0, 28, ..., 728; a test checks that it reproduces the progression times and course counts of
trackc_sim.simulate). Arms: front6 (the standard arm), front3, front12, best_fixed (best mean restricted time at 730 d on the
selection patients; ties go to fewer courses), informed (per patient the schedule with the largest min(ttp, 730 d); ties go
to fewer courses) and informed_lam10 (the same with a cost of 10 days of restricted time per course, an ASSUMED value).
M1, mechanism at the primary scan interval (84 d). Per patient and arm: the RANO-like call time; the first cycle boundary at
which the true burden is at least 1.4 x baseline (the level of the RANO-like threshold when no scan has shown a fall; 730 d if
never) and at least 1.2 x baseline (the baseline-referenced level; 730 d if never), both on the trajectory of the RANO-like run;
the true burden relative to baseline at the call; whether the call was made while the burden was still below 1.2 x baseline;
whether a scan ever showed a fall to below 0.9 of baseline; the days spent at or above 1.2 x baseline; the mean log burden ratio
over days 0 to 728. Averages over patients are taken within a world.
M2, scan-interval sweep. For scan intervals of 28, 56, 84 and 112 days the RANO-like progression time of all 56 schedules is
recomputed on the same evaluation patients. The frontier is the upper concave envelope of the (mean courses, restricted mean
time to 730 d) points of ALL 56 fixed schedules (richer than the dozen fixed arms of script 202, so the gaps are less inflated);
the gap of the informed arm is its restricted mean time minus the frontier at its own mean number of courses (undefined when
no fixed schedule uses so few courses).

Readout fixed before any result (median over worlds, informed against front3, the fixed arm nearest in courses; d_X is the
gain of the informed arm in the mean of measure X; rano = RANO-like time, c14 = first crossing of 1.4 x baseline, base = first
crossing of 1.2 x baseline):
  "no gain"                      d_rano <= 0;
  "mostly disease control"       d_c14 > 0.75 d_rano and d_base > 0.75 d_rano and the mean log burden of the informed arm is lower;
  "mostly the scan grid"         d_c14 < 0.25 d_rano and d_base < 0.25 d_rano and the share of calls made below 1.2 x baseline
                                 differs by less than 0.20 (the burden path reaches the levels no later, only the call is later);
  "mostly the nadir reference"   d_base < 0.25 d_rano and the share of calls made below 1.2 x baseline is lower by at least 0.20;
  "mixed"                        otherwise (checked in this order).
For M2 the headroom is called "scan-grid dependent" if the median gap at 28 d is less than half the median gap at 84 d.
These thresholds are arbitrary and are stated only so that the labels cannot be chosen after seeing the numbers.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev): patients seed [selection seed, world index] and
[evaluation seed, world index] per configuration, see M_SEEDS. No noise key (open loop).
Written to output/trackc_stage2_M_endpoint_mechanism<suffix>/ (manifest.json first, then arms_by_world.csv,
sweep_by_world.csv, results.json).

Run:  python src/204_trackc_endpoint_mechanism.py --config single --workers 4
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
import trackc_regime as rg  # noqa: E402

SCRIPT = "src/204_trackc_endpoint_mechanism.py"
OUT_NAME = "trackc_stage2_M_endpoint_mechanism"
M_SEEDS = {"single": (20261860, 20261861), "single_rep": (20261862, 20261863), "mixture": (20261864, 20261865)}   # (selection, evaluation)
N_SEL = 200
N_EVAL = 500
TAU = 730.0
N_REC = 27                           # cycle boundaries 0..26, i.e. days 0..728
LIMIT = 1.2                          # level of the baseline-referenced rule, relative to baseline
RANO_LEVEL = 1.4                     # level of the RANO-like rule when no scan has shown a fall (cal.tc.RANO_FACTOR)
DIP = 0.9                            # a scan "shows a response" if the burden is below this fraction of baseline
LAM_COST = 10.0                      # ASSUMED cost of one course, in days of restricted time, for informed_lam10
SCAN_SWEEP = (1, 2, 3, 4)            # scan interval in cycles: 28, 56, 84 (primary), 112 days
ARMS = ("front6", "front3", "front12", "best_fixed", "informed", "informed_lam10")
METRICS = ("mean_courses", "rmst730_rano", "rmst730_c14", "rmst730_base", "mean_days_above", "mean_log_ratio", "frac_called",
           "frac_called_below_B0", "frac_called_below_1p2", "median_ratio_call", "frac_scan_dip", "median_ratio_728")
CONTRASTS = (("informed", "front3"), ("informed", "best_fixed"), ("informed", "front12"), ("informed_lam10", "front3"))
WEAKEST = [
    "open loop and noise-free: this is the disease course under each schedule, not a rule that has to read noisy scans",
    "the informed arms are maxima over the same 56 schedules and use the true parameters; they are not realisable and are a lower bound of the true optimum",
    "the burden path is sampled at the 28-day cycle boundaries, so times above a level and first crossings are known to within 28 days",
    "after a RANO-like call the drug stops, as in the other scripts, so the burden after the call is the untreated regrowth; crossings and times above a level include that regrowth",
    "the 1.2 x baseline crossing is measured on the RANO-run trajectory, not on a run that stops the drug at the baseline-referenced call, so it is not identical to endpoint 3a of script 194",
    "the frontier of part M2 is built in-sample on the evaluation patients from the 56 fixed schedules, so it is the best fixed dosing can do on those patients; it still omits schedules outside the family",
    "the labels use arbitrary thresholds; only the numbers carry information",
    "controlled worlds only (weak kill by construction, tests C and D); a sample of an assumed family; 500 evaluation patients per world",
]


def _smoke(args):
    global N_SEL, N_EVAL
    if args.smoke:
        N_SEL, N_EVAL = 40, 60


def read_controlled(d_dir):
    with open(d_dir / "controlled_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


# ---------------------------------------------------------------- simulation with the burden path
def rano_run_with_burden(scen, mask, scan_every=rg.SCAN_EVERY):
    """Open-loop run with the RANO-like rule of trackc_sim.simulate that also keeps the true burden at the start of every
    cycle j = 0..N_REC-1 (days 28 j). mask: (N, N_CYCLES) bool, the courses scheduled for each patient; as in simulate, no
    course is given after the call or beyond the budget. Returns ttp in days (HORIZON if never called, as simulate does),
    courses used, the burden array (N, N_REC) and the reference value (scanned nadir, floored) that the rule used at the
    call (NaN if not called)."""
    tc = cal.tc
    n = len(scen["r"])
    S = (1 - scen["f0"]) * scen["B0"]
    R = scen["f0"] * scen["B0"]
    ttp = np.full(n, np.nan)
    used = np.zeros(n, int)
    burden = np.empty((n, N_REC))
    ref_at_call = np.full(n, np.nan)
    base0 = nadir = None
    for j in range(tc.N_CYCLES):                 # all 52 cycles, so ttp and courses equal those of trackc_sim.simulate
        scan = (j % scan_every == 0)
        B_true = S + R
        if j < N_REC:
            burden[:, j] = B_true
        if j == 0:
            base0, nadir = B_true.copy(), B_true.copy()
        else:
            ref = tc.RANO_FACTOR * np.maximum(nadir, tc.RANO_FLOOR * base0)
            hit = scan & np.isnan(ttp) & (B_true >= ref)
            ttp[hit] = j * tc.CYCLE_DAYS
            ref_at_call[hit] = (ref / tc.RANO_FACTOR)[hit]
            nadir = np.where(scan, np.minimum(nadir, B_true), nadir)
        treat = mask[:, j] & (used < tc.BUDGET) & np.isnan(ttp)
        used += treat
        S, R = tc.run_cycle(S, R, scen, treat, j * tc.CYCLE_DAYS)
    ttp = np.where(np.isnan(ttp), tc.HORIZON, ttp)           # never called within the 52 cycles: capped at the horizon, as in simulate
    return {"ttp": ttp, "courses": used, "burden": burden, "ref_at_call": ref_at_call}


def evaluate_schedules_burden(scen, schedules, scan_every=rg.SCAN_EVERY):
    """Every schedule on every patient: arrays shaped (n_schedules, n_patients, ...). Same layout as rg.evaluate_schedules."""
    n_s, n_p = len(schedules), len(scen["r"])
    tiled = {k: np.tile(v, n_s) for k, v in scen.items()}
    mask = np.repeat(schedules, n_p, axis=0)
    out = rano_run_with_burden(tiled, mask, scan_every)
    return {k: v.reshape((n_s, n_p) + v.shape[1:]) for k, v in out.items()}


# ---------------------------------------------------------------- pure helpers (tested with synthetic data)
def first_crossing(ratio, level):
    """First cycle boundary (days, 1..N_REC-1 cycles) at which burden/baseline >= level; TAU if never."""
    above = ratio[:, 1:] >= level
    return np.where(above.any(axis=1), (above.argmax(axis=1) + 1) * 28.0, TAU)


def patient_metrics(burden, ttp, courses):
    """Per-patient measures from a burden path (N, N_REC), call times (days; at least 730 means not called) and courses."""
    burden = np.asarray(burden, float)
    n, k = burden.shape
    ratio = burden / burden[:, [0]]
    ttp = np.asarray(ttp, float)
    called = ttp <= (k - 1) * 28.0 + 1e-9
    idx = np.clip(np.rint(ttp / 28.0).astype(int), 0, k - 1)
    ratio_call = np.where(called, ratio[np.arange(n), idx], np.nan)
    return {"courses": np.asarray(courses, float), "rano": np.minimum(ttp, TAU), "c14": first_crossing(ratio, RANO_LEVEL),
            "base": first_crossing(ratio, LIMIT), "days_above": (ratio[:, 1:] >= LIMIT).sum(axis=1) * 28.0,
            "log_ratio": np.log(ratio).mean(axis=1), "called": called, "ratio_call": ratio_call,
            "scan_dip": ratio[:, ::rg.SCAN_EVERY].min(axis=1) < DIP, "ratio_end": ratio[:, -1]}


def summarise_arm(m):
    """World-level summary of the per-patient measures of one arm."""
    called = m["called"]
    rc = m["ratio_call"][called]
    return {"mean_courses": float(m["courses"].mean()), "rmst730_rano": float(m["rano"].mean()),
            "rmst730_c14": float(m["c14"].mean()), "rmst730_base": float(m["base"].mean()),
            "mean_days_above": float(m["days_above"].mean()), "mean_log_ratio": float(m["log_ratio"].mean()),
            "frac_called": float(called.mean()),
            "frac_called_below_B0": float((rc < 1.0).mean()) if rc.size else float("nan"),
            "frac_called_below_1p2": float((rc < LIMIT).mean()) if rc.size else float("nan"),
            "median_ratio_call": float(np.median(rc)) if rc.size else float("nan"),
            "frac_scan_dip": float(m["scan_dip"].mean()), "median_ratio_728": float(np.median(m["ratio_end"]))}


def upper_envelope(x, y):
    """Vertices (sorted by x) of the upper concave envelope of the points (x, y)."""
    best = {}
    for a, b in zip(np.asarray(x, float).tolist(), np.asarray(y, float).tolist()):
        best[a] = max(best.get(a, -np.inf), b)
    pts = sorted(best.items())
    hull = []
    for p in pts:
        while len(hull) >= 2 and ((hull[-1][0] - hull[-2][0]) * (p[1] - hull[-2][1]) - (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0])) >= 0:
            hull.pop()
        hull.append(p)
    return hull


def frontier_at(hull, x):
    """Value of the envelope at x (linear between vertices); NaN outside the range of the vertices."""
    xs = [h[0] for h in hull]
    ys = [h[1] for h in hull]
    if x < xs[0] - 1e-9 or x > xs[-1] + 1e-9:
        return float("nan")
    return float(np.interp(x, xs, ys))


def pick_arms(res_sel, res_ev, names):
    """Per-patient arrays of every arm: returns {arm: dict of burden, ttp, courses} and the name of best_fixed."""
    ttp_sel, c_sel = res_sel["ttp"], res_sel["courses"]
    ttp_ev, c_ev = res_ev["ttp"], res_ev["courses"]
    n_s, n_p = ttp_ev.shape
    cols = np.arange(n_p)
    score_sel = [float(np.minimum(t, TAU).mean()) - 1e-6 * float(c.mean()) for t, c in zip(ttp_sel, c_sel)]
    best = int(np.argmax(score_sel))
    picks = {"front6": np.full(n_p, names.index("front6")), "front3": np.full(n_p, names.index("front3")),
             "front12": np.full(n_p, names.index("front12")), "best_fixed": np.full(n_p, best)}
    picks["informed"] = (np.minimum(ttp_ev, TAU) - 1e-6 * c_ev).argmax(0)
    picks["informed_lam10"] = (np.minimum(ttp_ev, TAU) - LAM_COST * c_ev - 1e-6 * c_ev).argmax(0)
    out = {}
    for arm, pick in picks.items():
        out[arm] = {"burden": res_ev["burden"][pick, cols], "ttp": ttp_ev[pick, cols], "courses": c_ev[pick, cols]}
    return out, names[best]


def sweep_gaps(ttp, courses):
    """Gap of the informed arms above the fixed frontier for one scan interval. ttp, courses: (n_schedules, n_patients)."""
    cols = np.arange(ttp.shape[1])
    capped = np.minimum(ttp, TAU)
    hull = upper_envelope(courses.mean(axis=1), capped.mean(axis=1))
    res = {}
    for tag, lam in (("informed", 0.0), ("informed_lam10", LAM_COST)):
        pick = (capped - lam * courses - 1e-6 * courses).argmax(0)
        x, y = float(courses[pick, cols].mean()), float(capped[pick, cols].mean())
        f = frontier_at(hull, x)
        res[tag] = {"courses": x, "rmst": y, "frontier": f, "gap": y - f if np.isfinite(f) else float("nan")}
    return res, hull


def evaluate_world(w, n_sel, n_ev, seeds):
    """All arms (M1) and the scan-interval sweep (M2) for one world. Pure function of the world and the seeds."""
    world = int(w["world"])
    sel = cal.draw_patients(w, n_sel, [seeds[0], world])
    ev = cal.draw_patients(w, n_ev, [seeds[1], world])
    sched, names = rg.candidate_schedules()
    res_ev = evaluate_schedules_burden(ev, sched)
    arms, best_name = pick_arms(evaluate_schedules_burden(sel, sched), res_ev, names)
    arm_rows, sweep_rows = [], []
    for arm in ARMS:
        a = arms[arm]
        row = {"world": world, "arm": arm, "note": best_name if arm == "best_fixed" else ""}
        row.update(summarise_arm(patient_metrics(a["burden"], a["ttp"], a["courses"])))
        arm_rows.append(row)
    for s in SCAN_SWEEP:
        if s == rg.SCAN_EVERY:
            ttp, courses = res_ev["ttp"], res_ev["courses"]
        else:
            ttp, courses = rg.evaluate_schedules(ev, sched, "rano", scan_every=s)
        gaps, _ = sweep_gaps(ttp, courses)
        capped = np.minimum(ttp, TAU)
        for tag, g in gaps.items():
            sweep_rows.append({"world": world, "scan_days": 28 * s, "arm": tag, "courses": g["courses"], "rmst730": g["rmst"],
                               "frontier": g["frontier"], "gap": g["gap"],
                               "front12_rmst730": float(capped[names.index("front12")].mean()),
                               "front6_rmst730": float(capped[names.index("front6")].mean())})
    return {"arms": arm_rows, "sweep": sweep_rows}


def _job(args):
    return evaluate_world(*args)


def _q(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    return [float(v) for v in np.percentile(x, [25, 50, 75])] if x.size else None


def verdict(contrast, early_key="frac_called_below_1p2"):
    """Label of the headroom from the median-world contrast of informed against front3 (rule fixed in the docstring)."""
    d_rano, d_c14, d_base = (contrast[k]["median"] for k in ("rmst730_rano", "rmst730_c14", "rmst730_base"))
    d_early, d_log = contrast[early_key]["median"], contrast["mean_log_ratio"]["median"]
    if d_rano <= 0:
        return "no gain"
    if d_c14 > 0.75 * d_rano and d_base > 0.75 * d_rano and d_log < 0:
        return "mostly disease control"
    if d_c14 < 0.25 * d_rano and d_base < 0.25 * d_rano and abs(d_early) < 0.20:
        return "mostly the scan grid"
    if d_base < 0.25 * d_rano and d_early <= -0.20:
        return "mostly the nadir reference"
    return "mixed"


def analyse(rows):
    """Summary over worlds of the arm rows (M1). Pure function, no I/O."""
    worlds = sorted({r["world"] for r in rows})
    by = {(r["world"], r["arm"]): r for r in rows}
    out = {"n_worlds": len(worlds), "arms": {}, "contrasts": {}}
    for arm in ARMS:
        out["arms"][arm] = {m: {"q25_q50_q75": _q([by[(w, arm)][m] for w in worlds])} for m in METRICS}
    for a, b in CONTRASTS:
        block = {}
        for m in METRICS:
            d = np.array([by[(w, a)][m] - by[(w, b)][m] for w in worlds], float)
            ok = np.isfinite(d)
            block[m] = {"median": float(np.median(d[ok])) if ok.any() else float("nan"),
                        "q25_q75": _q(d)[::2] if ok.any() else None,
                        "n_positive": int((d[ok] > 0).sum()), "n_negative": int((d[ok] < 0).sum())}
        block["verdict"] = verdict(block) if (a, b) == ("informed", "front3") else None
        out["contrasts"][f"{a}_minus_{b}"] = block
    return out


def analyse_sweep(rows):
    """Summary over worlds of the scan-interval sweep (M2): median gap and the worlds above the frontier, per interval and arm."""
    out = {}
    for arm in ("informed", "informed_lam10"):
        out[arm] = {}
        for days in sorted({r["scan_days"] for r in rows}):
            sel = [r for r in rows if r["arm"] == arm and r["scan_days"] == days]
            gap = np.array([r["gap"] for r in sel], float)
            ok = np.isfinite(gap)
            out[arm][str(days)] = {"n_worlds": int(ok.sum()), "n_outside_fixed_range": int((~ok).sum()),
                                   "gap_q25_q50_q75": _q(gap), "n_above": int((gap[ok] > 0.05).sum()), "n_below": int((gap[ok] < -0.05).sum()),
                                   "median_courses": float(np.median([r["courses"] for r in sel])),
                                   "median_front12_rmst730": float(np.median([r["front12_rmst730"] for r in sel])),
                                   "median_front6_rmst730": float(np.median([r["front6_rmst730"] for r in sel]))}
    g = out["informed"]
    g28, g84 = g["28"]["gap_q25_q50_q75"], g["84"]["gap_q25_q50_q75"]
    out["verdict_informed"] = ("scan-grid dependent" if (g28 and g84 and g84[1] > 0 and g28[1] < 0.5 * g84[1]) else "not scan-grid dependent")
    return out


def write_csv(rows, out_dir, name):
    cols = list(rows[0].keys())
    with open(out_dir / name, "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(cols)
        for r in rows:
            wr.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else r[c] for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    if cfg["name"] not in M_SEEDS:
        sys.exit(f"no seeds for config {cfg['name']}")
    seeds = M_SEEDS[cfg["name"]]
    OUT_DIR, D_DIR = PROJECT_ROOT / "output" / (OUT_NAME + cfg["suffix"] + ("_smoke" if args.smoke else "")), cal.stage_dir("D", cfg, args.smoke)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (D_DIR / "results.json").exists():
        sys.exit(f"script 193 has not finished: {D_DIR / 'results.json'} is missing")
    worlds = read_controlled(D_DIR)
    if not worlds:
        sys.exit("script 193 found no controlled world")
    sched, names = rg.candidate_schedules()
    protocol = {"script": "204_trackc_endpoint_mechanism", "stage": "2M", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_sel": N_SEL, "n_eval": N_EVAL, "selection_seed": seeds[0], "evaluation_seed": seeds[1],
                "tau": TAU, "endpoint": "rano", "n_record_boundaries": N_REC, "limit": LIMIT, "rano_level": RANO_LEVEL, "dip": DIP,
                "lam_cost": LAM_COST, "scan_sweep_cycles": list(SCAN_SWEEP), "arms": list(ARMS), "contrasts": [list(c) for c in CONTRASTS],
                "candidate_schedules": names, "scan_every": rg.SCAN_EVERY, "controlled_worlds": [w["world"] for w in worlds],
                "input": str((D_DIR / "controlled_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((D_DIR / "results.json").read_bytes()).hexdigest(),
                "rano_factor": cal.tc.RANO_FACTOR, "rano_floor": cal.tc.RANO_FLOOR}
    manifest = man.build_manifest("204_trackc_endpoint_mechanism", protocol, purpose="stage2_dev",
                                  scenario=seeds, noise=None, script=SCRIPT,
                                  require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py", "src/trackc_regime.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(worlds)} controlled worlds", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = list(ex.map(_job, [(w, N_SEL, N_EVAL, seeds) for w in worlds]))
    arm_rows = [r for part in parts for r in part["arms"]]
    sweep_rows = [r for part in parts for r in part["sweep"]]
    write_csv(arm_rows, OUT_DIR, "arms_by_world.csv")
    write_csv(sweep_rows, OUT_DIR, "sweep_by_world.csv")
    res = analyse(arm_rows)
    res["sweep"] = analyse_sweep(sweep_rows)
    res.update({"script": "204_trackc_endpoint_mechanism",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; what the RANO-like headroom consists of",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"], "weakest_points": WEAKEST,
                "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for arm in ARMS:
        a = res["arms"][arm]
        print(f"{arm:15s} courses {a['mean_courses']['q25_q50_q75'][1]:5.2f}  rano730 {a['rmst730_rano']['q25_q50_q75'][1]:6.1f}  c14 {a['rmst730_c14']['q25_q50_q75'][1]:6.1f}  "
              f"base730 {a['rmst730_base']['q25_q50_q75'][1]:6.1f}  log burden {a['mean_log_ratio']['q25_q50_q75'][1]:+.3f}  "
              f"called below 1.2 {a['frac_called_below_1p2']['q25_q50_q75'][1]:.2f}  dip {a['frac_scan_dip']['q25_q50_q75'][1]:.2f}")
    for name, c in res["contrasts"].items():
        print(f"{name}: rano {c['rmst730_rano']['median']:+.1f} d, c14 {c['rmst730_c14']['median']:+.1f} d, base {c['rmst730_base']['median']:+.1f} d, "
              f"log burden {c['mean_log_ratio']['median']:+.3f}, called below 1.2 {c['frac_called_below_1p2']['median']:+.2f}" + (f"  -> {c['verdict']}" if c["verdict"] else ""))
    for arm in ("informed", "informed_lam10"):
        for days, e in res["sweep"][arm].items():
            print(f"sweep {arm:15s} scan every {days:>3s} d: gap {e['gap_q25_q50_q75'][1]:+6.1f} d (q25 {e['gap_q25_q50_q75'][0]:+.1f}, q75 {e['gap_q25_q50_q75'][2]:+.1f}) "
                  f"courses {e['median_courses']:.2f}, above {e['n_above']} below {e['n_below']} outside {e['n_outside_fixed_range']}")
    print("sweep verdict:", res["sweep"]["verdict_informed"])


if __name__ == "__main__":
    main()
