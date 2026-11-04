#!/usr/bin/env python3
"""
Script 203: Track C Stage 2, test H: which patient-level quantity is worth measuring?
====================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. The design below was written and committed BEFORE any
result of this script existed (it was written while the regime maps of the mixed-cap configurations were still running).

Question. Script 194 and script 202 show that, in the controlled worlds, a controller that knows every patient's parameters
(the "informed" arms: the best of 56 open-loop schedules for each patient) reaches far more restricted time than any rule that
reads scans, at the same mean number of courses. Which single patient-level quantity carries that value? A quantity that
carries most of it is the one worth measuring; a quantity that carries little is not.

Method (open loop, so no measurement noise and no filter). For every controlled world of a configuration (script 193) draw two
sets of simulated patients: a selection set and an evaluation set. For each of the 56 candidate schedules of
src/trackc_regime.py compute the time to progression of every patient (RANO-like rule on the true burden). The score of
schedule s for patient i is  min(ttp, 730 d) - lambda * courses  with lambda = 0, 5 and 10 days per course (an ASSUMED grid;
a tie goes to fewer courses). Then, for each lambda:
  V_fixed    the schedule with the best mean score on the selection patients, evaluated on the evaluation patients;
  V_perfect  for every evaluation patient the best of the 56 schedules (the patient's parameters are known exactly);
  V_X        the patients are split into groups by quantity X (K = 4 quantile groups for one quantity, 3 x 3 groups for a
             pair of quantities; the responder class of the mixture family has 2 groups); in each group the schedule with the
             best mean score on the selection patients of that group is applied to the evaluation patients of that group.
  share_X = (V_X - V_fixed) / (V_perfect - V_fixed)  (undefined when the denominator is below 1e-6 days).
Quantities: r (growth rate), m (kill strength relative to growth), f0 (initial resistant fraction), eps (residual kill of
resistant cells), B0 (initial burden), and for the mixture family the responder class. B0 is read from the baseline scan, so
it costs nothing to measure; m shows in the response to the first courses; r shows in the regrowth between courses; f0 and
eps are not visible on a scan. These links are descriptive and are not tested here.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev): patients seed [selection seed, world index] and
[evaluation seed, world index] per configuration, see H_SEEDS below. No noise key (open loop).
Written to output/trackc_stage2_H_value_of_information<suffix>/ (manifest.json first, then value_by_world.csv, results.json).

Run:  python src/203_trackc_value_of_information.py --config single_mix --workers 4
      python src/203_trackc_value_of_information.py --config single_mix --smoke
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import hashlib
import itertools
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

SCRIPT = "src/203_trackc_value_of_information.py"
OUT_NAME = "trackc_stage2_H_value_of_information"
H_SEEDS = {"single": (20261850, 20261851), "single_rep": (20261852, 20261853), "mixture": (20261854, 20261855),
           "single_mix": (20261856, 20261857), "mixture_mix": (20261858, 20261859)}      # (selection, evaluation)
N_SEL = 1000
N_EVAL = 1000
LAMBDAS = (0.0, 5.0, 10.0)                  # days of restricted time one course is worth (ASSUMED grid)
TAU = 730.0
ENDPOINT = "rano"
K_SINGLE = 4
K_PAIR = 3
QUANTITIES = ("r", "m", "f0", "eps", "B0")
MIN_DENOMINATOR = 1e-6
WEAKEST = [
    "the 'perfect' value is the best of 56 open-loop schedules per patient, a lower bound of what perfect control could reach; every share is relative to that family, not to the true optimum",
    "each quantity is known exactly and without delay; a real measurement has error, and r, f0 and eps are not visible on a scan, so a share is an upper bound on what a measurement could deliver",
    "the groups are chosen from 1,000 selection patients, so a group policy carries selection noise (a pair has 9 groups of about 110 patients); the shares are conservative for pairs",
    "the patient-level quantities are drawn independently inside each world; in real patients they correlate, and the value of one quantity depends on that structure",
    "open loop only: this is the value of knowing a quantity before treatment, not the value of any feedback rule; scan noise and filters are not involved",
    "one endpoint (the RANO-like rule on the true burden), one horizon (730 d), an assumed cost grid, and controlled worlds only (weak kill by construction, test C and D)",
    "the controlled worlds are a sample of an assumed family; the cross-world summaries are descriptive and uncorrected for the number of quantities",
]


def _smoke(args):
    global N_SEL, N_EVAL
    if args.smoke:
        N_SEL, N_EVAL = 60, 60


def read_controlled(d_dir):
    with open(d_dir / "controlled_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


# ---------------------------------------------------------------- pure helpers (tested with synthetic data)
def quantile_labels(x, k, reference=None):
    """Group index 0..k-1 of each value of x by the quantiles of `reference` (default x). Ties go to the higher group.
    The cut points are shared by the selection and the evaluation patients, so both use the same groups."""
    x = np.asarray(x, float)
    ref = x if reference is None else np.asarray(reference, float)
    edges = np.quantile(ref, np.linspace(0, 1, k + 1)[1:-1])
    return np.searchsorted(edges, x, side="right")


def group_value(score_sel, score_ev, lab_sel, lab_ev):
    """Mean evaluation score when each group uses the schedule that is best on the selection patients of that group.
    score_*: (n_schedules, n_patients); lab_*: group index per patient. A group with no selection patient uses the schedule
    that is best over all selection patients."""
    n_ev = score_ev.shape[1]
    overall = int(np.argmax(score_sel.mean(axis=1)))
    total = 0.0
    for g in np.unique(lab_ev):
        me = lab_ev == g
        ms = lab_sel == g
        best = int(np.argmax(score_sel[:, ms].mean(axis=1))) if ms.any() else overall
        total += float(score_ev[best, me].sum())
    return total / n_ev


def perfect_value(score_ev):
    """Mean of the best schedule per evaluation patient."""
    return float(score_ev.max(axis=0).mean())


def share(v_x, v_fixed, v_perfect):
    den = v_perfect - v_fixed
    return float((v_x - v_fixed) / den) if den > MIN_DENOMINATOR else float("nan")


def score_matrix(ttp, courses, lam):
    """Score of every schedule for every patient: restricted time minus lambda per course (tiny tie-break to fewer courses)."""
    return np.minimum(ttp, TAU) - lam * courses - 1e-6 * courses


def quantity_groups(sel, ev, w):
    """Group labels of the selection and evaluation patients for every single quantity and every pair of quantities.
    Returns {name: (lab_sel, lab_ev)}; names are 'r', 'm', ..., 'r+m', ... and, for the mixture family, 'class'."""
    out = {}
    pooled = {q: np.r_[sel[q], ev[q]] for q in QUANTITIES}
    for q in QUANTITIES:
        out[q] = (quantile_labels(sel[q], K_SINGLE, pooled[q]), quantile_labels(ev[q], K_SINGLE, pooled[q]))
    for a, b in itertools.combinations(QUANTITIES, 2):
        la = (quantile_labels(sel[a], K_PAIR, pooled[a]), quantile_labels(ev[a], K_PAIR, pooled[a]))
        lb = (quantile_labels(sel[b], K_PAIR, pooled[b]), quantile_labels(ev[b], K_PAIR, pooled[b]))
        out[f"{a}+{b}"] = (la[0] * K_PAIR + lb[0], la[1] * K_PAIR + lb[1])
    if "q" in w:
        out["class"] = (cal.responder_mask(sel, w).astype(int), cal.responder_mask(ev, w).astype(int))
    return out


def evaluate_world(w, n_sel, n_ev, seeds):
    """All quantities, all lambdas for one world. Pure function of the world and the seeds; returns long-format rows."""
    world = int(w["world"])
    sel = cal.draw_patients(w, n_sel, [seeds[0], world])
    ev = cal.draw_patients(w, n_ev, [seeds[1], world])
    sched, _ = rg.candidate_schedules()
    ttp_sel, courses_sel = rg.evaluate_schedules(sel, sched, ENDPOINT)
    ttp_ev, courses_ev = rg.evaluate_schedules(ev, sched, ENDPOINT)
    groups = quantity_groups(sel, ev, w)
    rows = []
    for lam in LAMBDAS:
        s_sel, s_ev = score_matrix(ttp_sel, courses_sel, lam), score_matrix(ttp_ev, courses_ev, lam)
        zeros = (np.zeros(n_sel, int), np.zeros(n_ev, int))
        v_fixed = group_value(s_sel, s_ev, *zeros)
        v_perfect = perfect_value(s_ev)
        for name, (lab_sel, lab_ev) in groups.items():
            v_x = group_value(s_sel, s_ev, lab_sel, lab_ev)
            rows.append({"world": world, "lam": lam, "quantity": name, "v_fixed": v_fixed, "v_perfect": v_perfect,
                         "v_group": v_x, "gain_days": v_x - v_fixed, "share": share(v_x, v_fixed, v_perfect)})
    return rows


def _job(args):
    return evaluate_world(*args)


def analyse(rows):
    """Summary over worlds for every lambda and quantity. Pure function, no I/O."""
    out = {"n_worlds": len({r["world"] for r in rows}), "lambdas": sorted({r["lam"] for r in rows}), "by_lambda": {}}
    for lam in out["lambdas"]:
        sel = [r for r in rows if r["lam"] == lam]
        worlds = sorted({r["world"] for r in sel})
        evpi = [next(r["v_perfect"] - r["v_fixed"] for r in sel if r["world"] == w) for w in worlds]
        block = {"perfect_minus_fixed_days": _q(evpi), "quantities": {}}
        for q in sorted({r["quantity"] for r in sel}):
            qs = sorted([r for r in sel if r["quantity"] == q], key=lambda r: r["world"])
            sh = np.array([r["share"] for r in qs], float)
            gain = np.array([r["gain_days"] for r in qs], float)
            ok = np.isfinite(sh)
            block["quantities"][q] = {
                "n_worlds_defined": int(ok.sum()),
                "gain_days_q25_q50_q75": _q(gain), "share_q25_q50_q75": _q(sh[ok]) if ok.any() else None,
                "n_worlds_share_ge_half": int((sh[ok] >= 0.5).sum()) if ok.any() else 0}
        singles = [q for q in block["quantities"] if "+" not in q and q != "class"]
        best = {}
        for w in worlds:
            cand = [(next(r["share"] for r in sel if r["world"] == w and r["quantity"] == q), q) for q in singles]
            cand = [c for c in cand if np.isfinite(c[0])]
            if cand:
                best[max(cand)[1]] = best.get(max(cand)[1], 0) + 1
        block["best_single_quantity_counts"] = dict(sorted(best.items(), key=lambda kv: -kv[1]))
        out["by_lambda"][f"{lam:g}"] = block
    return out


def _q(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    return [float(v) for v in np.percentile(x, [25, 50, 75])] if x.size else None


def write_csv(rows, out_dir):
    cols = list(rows[0].keys())
    with open(out_dir / "value_by_world.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            w.writerow([f"{r[c]:.6g}" if isinstance(r[c], float) else r[c] for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="single", help="run configuration of src/trackc_calibration.py")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    if cfg["name"] not in H_SEEDS:
        sys.exit(f"no seeds for config {cfg['name']}")
    seeds = H_SEEDS[cfg["name"]]
    OUT_DIR, D_DIR = PROJECT_ROOT / "output" / (OUT_NAME + cfg["suffix"] + ("_smoke" if args.smoke else "")), cal.stage_dir("D", cfg, args.smoke)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (D_DIR / "results.json").exists():
        sys.exit(f"script 193 has not finished: {D_DIR / 'results.json'} is missing")
    worlds = read_controlled(D_DIR)
    if not worlds:
        sys.exit("script 193 found no controlled world")
    sched, names = rg.candidate_schedules()
    protocol = {"script": "203_trackc_value_of_information", "stage": "2H", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_sel": N_SEL, "n_eval": N_EVAL, "selection_seed": seeds[0], "evaluation_seed": seeds[1],
                "lambdas": list(LAMBDAS), "tau": TAU, "endpoint": ENDPOINT, "k_single": K_SINGLE, "k_pair": K_PAIR,
                "quantities": list(QUANTITIES), "candidate_schedules": names, "scan_every": rg.SCAN_EVERY,
                "controlled_worlds": [w["world"] for w in worlds],
                "input": str((D_DIR / "controlled_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((D_DIR / "results.json").read_bytes()).hexdigest(),
                "rano_factor": cal.tc.RANO_FACTOR, "rano_floor": cal.tc.RANO_FLOOR}
    manifest = man.build_manifest("203_trackc_value_of_information", protocol, purpose="stage2_dev",
                                  scenario=seeds, noise=None, script=SCRIPT,
                                  require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py", "src/trackc_regime.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(worlds)} controlled worlds", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = list(ex.map(_job, [(w, N_SEL, N_EVAL, seeds) for w in worlds]))
    rows = [r for part in parts for r in part]
    write_csv(rows, OUT_DIR)
    res = analyse(rows)
    res.update({"script": "203_trackc_value_of_information",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; open-loop value of knowing a patient-level quantity",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"], "weakest_points": WEAKEST,
                "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for lam, block in res["by_lambda"].items():
        print(f"== lambda {lam} d per course: perfect minus fixed (25th, 50th, 75th percentile over worlds) {block['perfect_minus_fixed_days']}")
        for q, e in sorted(block["quantities"].items(), key=lambda kv: -(kv[1]["share_q25_q50_q75"] or [0, -9, 0])[1]):
            if e["share_q25_q50_q75"]:
                print(f"   {q:8s} share median {e['share_q25_q50_q75'][1]:5.2f}  gain median {e['gain_days_q25_q50_q75'][1]:6.1f} d  "
                      f"share>=0.5 in {e['n_worlds_share_ge_half']}/{e['n_worlds_defined']} worlds")
        print("   best single quantity per world:", block["best_single_quantity_counts"])


if __name__ == "__main__":
    main()
