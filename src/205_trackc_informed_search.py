"""Script 205, test G (added 2026-11-04): how far below the best open-loop schedule is the "informed" arm, and how far below
the best fixed dosing is the fixed frontier? SIMULATION ONLY; every range and rule is an assumption; development seeds.

Why. The informed arm of scripts 194, 202, 203 and 204 takes, for every patient, the best of the 56 open-loop schedules of
src/trackc_regime.py. Without measurement noise and with known parameters the dynamics are deterministic, so no closed-loop
controller can beat the best open-loop schedule of the same patient. The informed arm is therefore a LOWER bound of the best any
controller can do; at the same time the fixed frontier (about a dozen or 56 fixed schedules) is a lower bound of the best fixed
dosing. The gap between the two can therefore err in either direction, and neither an upper nor a lower bound holds for it.
This script searches beyond the 56 schedules on both sides with the same method and the same accounting.

Method (open loop, no measurement noise). For every controlled world of one configuration (script 193) draw selection patients
and evaluation patients (separate seeds). Schedules live on the first NC = 27 cycles (days 0 to 756; nothing later can change a
restricted time of 730 days); every schedule of the family of 56 is cut to these cycles and duplicates are removed, so courses
are counted within these cycles only (this differs from scripts 194 to 204 for the schedules of the family that have courses
after cycle 26; nothing else changes). Courses are limited to 12. Endpoints: "baseline" (progression when the true burden
reaches 1.2 x baseline, detected continuously) and "rano" (RANO-like rule at scans every 84 days). Value of a patient under a
schedule: min(time to progression, 730 d) - (lam + 1e-6) x courses, so that ties go to fewer courses.
G1, informed optimum. Per patient: start from each of the N_START best schedules of the family and climb by steepest ascent over
the moves "toggle one cycle" (27 moves) and "swap two neighbouring cycles" (26 moves) until no move improves or ROUNDS_G1 rounds
are used; keep the best result. The cost lam is 0 and LAM_COST days of restricted time per course (an ASSUMED value).
G2, fixed frontier. For mu in MUS (days per course, assumed) climb the shared schedule that maximises the mean restricted time
minus (mu + 1e-6) x mean courses on the selection patients, starting from the best schedule of the family, and keep every
schedule on the climbing path. Evaluate these schedules and the 56 of the family on the evaluation patients; the searched
frontier is the upper concave envelope of all their (mean courses, mean restricted time) points; the family frontier uses the
56 alone. Gaps are taken at the courses of the informed arm: gap = restricted mean time of the informed arm minus the frontier.
Reported per world, endpoint and lam: the mean gain of the search over the best of the family, the share of patients improved,
the share still improving at the round limit, the courses and restricted times, and the gaps above the family frontier (before
and after the search) and above the searched frontier.

Readout fixed before any result (median over worlds, lam = 0):
  G1 label per endpoint:   "best of 56 is near the optimum"   if the median over worlds of the mean gain is below 2 d;
                           "best of 56 understates the informed optimum"   if it is 5 d or more;   "intermediate" otherwise;
  room label per endpoint: "room above fixed dosing stays under 10 d"   if the median over worlds of the gap of the searched
                           informed arm above the searched frontier is below 10 d;   "room above fixed dosing is 10 d or more"
                           otherwise.
These thresholds are arbitrary and are stated only so that the labels cannot be chosen after seeing the numbers.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev): patients seed [selection seed, world index] and
[evaluation seed, world index], see G_SEEDS. No noise key (open loop).
Written to output/trackc_stage2_G_informed_search<suffix>/ (manifest.json first, then search_by_world.csv,
frontier_by_world.csv, results.json).

Run:  python src/205_trackc_informed_search.py --config mixture --workers 4
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
import trackc_policies as pol  # noqa: E402
import trackc_regime as rg  # noqa: E402

SCRIPT = "src/205_trackc_informed_search.py"
OUT_NAME = "trackc_stage2_G_informed_search"
G_SEEDS = {"single": (20261870, 20261871), "single_rep": (20261872, 20261873), "mixture": (20261874, 20261875),
           "single_mix": (20261876, 20261877), "mixture_mix": (20261878, 20261879)}      # (selection, evaluation)
N_SEL = 100
N_EVAL = 150
TAU = 730.0
NC = 27                              # cycles 0..26 start at days 0..728; later cycles cannot change a restricted time of 730 d
EPS = 1e-6                           # tie-break: fewer courses
TOL = 1e-9
N_START = 2
ROUNDS_G1 = 8
ROUNDS_G2 = 10
LAM_COST = 10.0                      # ASSUMED cost of one course in days of restricted time (as in script 204)
LAMS = (0.0, LAM_COST)
MUS = (0.0, 5.0, 20.0)               # ASSUMED prices of a course used to trace the fixed frontier
ENDPOINTS = ("baseline", "rano")
NEAR_GAIN = 2.0                      # label thresholds (arbitrary)
FAR_GAIN = 5.0
ROOM_LIMIT = 10.0
WEAKEST = [
    "the search is a local steepest-ascent search with toggle and neighbour-swap moves from the two best schedules of the family; it gives a lower bound of the best schedule of a patient, not the optimum, and other starts could reach higher values",
    "open loop and noise-free: with known parameters and deterministic dynamics no closed-loop controller beats the best open-loop schedule, but a rule that has to read noisy scans cannot know the parameters, so these values are not realisable",
    "the informed arm is chosen on the evaluation patients themselves (it knows them); the fixed schedules of the search are found on separate selection patients and evaluated on the evaluation patients, but the envelope takes the best points over the evaluation patients, so the frontier carries a small in-sample optimism, as in part M2 of script 204",
    "courses are counted within the first 27 cycles, which differs from scripts 194 to 204 for the schedules of the family that have courses after cycle 26",
    "the costs lam and mu are assumed values; the labels use arbitrary thresholds; only the numbers carry information",
    "controlled worlds only (weak kill by construction, tests C and D); one seed; 150 evaluation patients and 100 selection patients per world",
    "moves are limited to the first 27 cycles and to 12 courses; the closed-loop-equals-open-loop argument holds for this noise-free model only",
]


def _smoke(args):
    global N_SEL, N_EVAL, ROUNDS_G1, ROUNDS_G2
    if args.smoke:
        N_SEL, N_EVAL, ROUNDS_G1, ROUNDS_G2 = 16, 24, 2, 2


def read_controlled(d_dir):
    with open(d_dir / "controlled_worlds.csv", encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


# ---------------------------------------------------------------- pure helpers (tested with synthetic data)
def truncate(schedules, nc=NC):
    """Schedules cut to the first nc cycles, duplicates removed (rows keep their first-seen order)."""
    s = np.array(schedules, bool, copy=True)
    s[:, nc:] = False
    _, first = np.unique(s, axis=0, return_index=True)
    return s[np.sort(first)]


def neighbor_masks(masks, nc=NC, budget=None):
    """masks: (n, N_CYCLES) bool -> (M, n, N_CYCLES) bool, M = nc + (nc - 1): toggle each of the first nc cycles, then swap each pair of
    neighbouring cycles. A move that would give more than `budget` courses (default 12) returns the mask unchanged."""
    budget = rg.tc.BUDGET if budget is None else budget
    masks = np.asarray(masks, bool)
    n, width = masks.shape
    m = nc + (nc - 1)
    out = np.repeat(masks[None], m, axis=0)
    for c in range(nc):
        out[c, :, c] = ~masks[:, c]
    for c in range(nc - 1):
        out[nc + c, :, c] = masks[:, c + 1]
        out[nc + c, :, c + 1] = masks[:, c]
    bad = out.sum(axis=2) > budget
    if bad.any():
        base = np.broadcast_to(masks[None], out.shape)
        out[bad] = base[bad]
    return out


def climb_patients(cur, cur_val, value_fn, rounds):
    """Steepest ascent per patient. cur: (n, width) bool, cur_val: (n,). value_fn(nb, ids) -> (M, len(ids)) values of the neighbour masks
    nb (M, len(ids), width) of the patients `ids`. A patient leaves the search when no neighbour improves its value by more than TOL.
    Returns (masks, values, rounds used, converged flags); values never decrease."""
    cur = np.array(cur, bool, copy=True)
    cur_val = np.array(cur_val, float, copy=True)
    n = len(cur_val)
    active = np.ones(n, bool)
    used = np.zeros(n, int)
    for _ in range(rounds):
        ids = np.flatnonzero(active)
        if ids.size == 0:
            break
        nb = neighbor_masks(cur[ids])
        val = value_fn(nb, ids)
        k = np.arange(ids.size)
        j = val.argmax(axis=0)
        best = val[j, k]
        better = best > cur_val[ids] + TOL
        sel = ids[better]
        cur[sel] = nb[j[better], k[better]]
        cur_val[sel] = best[better]
        used[sel] += 1
        active[ids[~better]] = False
    return cur, cur_val, used, ~active


def climb_fixed(mask0, objective_fn, rounds):
    """Steepest ascent of one shared schedule. objective_fn(masks (M, width)) -> (M,) objective. Returns the path of
    (mask, objective) from the start to the last improvement."""
    cur = np.array(mask0, bool, copy=True)
    cur_val = float(objective_fn(cur[None])[0])
    path = [(cur.copy(), cur_val)]
    for _ in range(rounds):
        nb = neighbor_masks(cur[None])[:, 0, :]
        val = np.asarray(objective_fn(nb), float)
        j = int(val.argmax())
        if val[j] <= cur_val + TOL:
            break
        cur, cur_val = nb[j].copy(), float(val[j])
        path.append((cur.copy(), cur_val))
    return path


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


def label_g1(median_gain):
    if median_gain < NEAR_GAIN:
        return "best of 56 is near the optimum"
    if median_gain >= FAR_GAIN:
        return "best of 56 understates the informed optimum"
    return "intermediate"


def label_room(median_gap):
    return "room above fixed dosing stays under 10 d" if median_gap < ROOM_LIMIT else "room above fixed dosing is 10 d or more"


# ---------------------------------------------------------------- simulation
def subset(scen, ids):
    return {k: np.asarray(v)[ids] for k, v in scen.items()}


def eval_per_patient(scen, nb, endpoint, scan_every=rg.SCAN_EVERY):
    """nb: (M, k, width) schedules, one per patient and move, scen: k patients. Returns ttp and courses, both (M, k)."""
    m, k = nb.shape[0], nb.shape[1]
    tiled = {key: np.tile(v, m) for key, v in scen.items()}
    tiled["scan_every"] = np.full(m * k, int(scan_every))
    out = rg.tc.simulate(tiled, pol.Schedule(nb.reshape(m * k, -1)), [0, 0], endpoint=endpoint)
    return out["ttp"].reshape(m, k), out["courses"].reshape(m, k)


def search_informed(ev, s56, ttp56, c56, endpoint, lam, n_start, rounds):
    """Per-patient search (G1). Returns the final masks, their values, rounds used and converged flags, and the values of the best of 56."""
    n = ttp56.shape[1]
    cols = np.arange(n)
    val56 = np.minimum(ttp56, TAU) - (lam + EPS) * c56
    order = np.argsort(-val56, axis=0, kind="stable")

    def value_fn(nb, ids):
        ttp, c = eval_per_patient(subset(ev, ids), nb, endpoint)
        return np.minimum(ttp, TAU) - (lam + EPS) * c

    best_mask = np.zeros((n, s56.shape[1]), bool)
    best_val = np.full(n, -np.inf)
    best_used = np.zeros(n, int)
    best_conv = np.zeros(n, bool)
    for s in range(min(n_start, len(s56))):
        idx = order[s]
        mask, val, used, conv = climb_patients(s56[idx], val56[idx, cols], value_fn, rounds)
        better = val > best_val + TOL          # the first start always counts (best_val is -inf); a later start must be strictly better
        best_mask[better], best_val[better] = mask[better], val[better]
        best_used[better], best_conv[better] = used[better], conv[better]
    return best_mask, best_val, best_used, best_conv, val56[order[0], cols]


def fixed_frontier_points(sel, s56, ttp_sel, c_sel, endpoint, mus, rounds):
    """Masks found by climbing the shared schedule on the selection patients for every price mu (G2), with the path kept."""
    found = []
    for mu in mus:
        start = s56[int((np.minimum(ttp_sel, TAU).mean(axis=1) - (mu + EPS) * c_sel.mean(axis=1)).argmax())]

        def objective(masks, mu=mu):
            ttp, c = rg.evaluate_schedules(sel, masks, endpoint)
            return np.minimum(ttp, TAU).mean(axis=1) - (mu + EPS) * c.mean(axis=1)

        found += [m for m, _ in climb_fixed(start, objective, rounds)]
    return truncate(np.array(found))


def evaluate_world(args):
    w, n_sel, n_ev, seeds, rounds1, rounds2 = args
    world = int(w["world"])
    sel = cal.draw_patients(w, n_sel, [seeds[0], world])
    ev = cal.draw_patients(w, n_ev, [seeds[1], world])
    sched, _ = rg.candidate_schedules()
    s56 = truncate(sched)
    cols = np.arange(n_ev)
    search_rows, front_rows = [], []
    for endpoint in ENDPOINTS:
        ttp56, c56 = rg.evaluate_schedules(ev, s56, endpoint)
        capped56 = np.minimum(ttp56, TAU)
        hull56 = upper_envelope(c56.mean(axis=1), capped56.mean(axis=1))
        ttp_sel, c_sel = rg.evaluate_schedules(sel, s56, endpoint)
        fixed = fixed_frontier_points(sel, s56, ttp_sel, c_sel, endpoint, MUS, rounds2)
        ttp_f, c_f = rg.evaluate_schedules(ev, fixed, endpoint)
        px = np.concatenate([c56.mean(axis=1), c_f.mean(axis=1)])
        py = np.concatenate([capped56.mean(axis=1), np.minimum(ttp_f, TAU).mean(axis=1)])
        hull_ls = upper_envelope(px, py)
        for lam in LAMS:
            mask, val, used, conv, val56 = search_informed(ev, s56, ttp56, c56, endpoint, lam, N_START, rounds1)
            ttp_i, c_i = eval_per_patient(ev, mask[None], endpoint)
            ttp_i, c_i = ttp_i[0], c_i[0]
            order = np.argsort(-(capped56 - (lam + EPS) * c56), axis=0, kind="stable")[0]
            x56, y56 = float(c56[order, cols].mean()), float(capped56[order, cols].mean())
            xls, yls = float(c_i.mean()), float(np.minimum(ttp_i, TAU).mean())
            gain = val - val56
            search_rows.append({"world": world, "endpoint": endpoint, "lam": lam, "n": n_ev, "mean_value_family": float(val56.mean()),
                                "mean_value_search": float(val.mean()), "mean_gain": float(gain.mean()), "frac_improved": float((gain > TOL).mean()),
                                "frac_not_converged": float((~conv).mean()), "courses_family": x56, "courses_search": xls,
                                "rmst_family": y56, "rmst_search": yls})
            f56_at56, f56_atls = frontier_at(hull56, x56), frontier_at(hull56, xls)
            fls_at56, fls_atls = frontier_at(hull_ls, x56), frontier_at(hull_ls, xls)
            front_rows.append({"world": world, "endpoint": endpoint, "lam": lam, "courses_family": x56, "rmst_family": y56,
                               "courses_search": xls, "rmst_search": yls,
                               "gap_family_vs_family_frontier": y56 - f56_at56, "gap_search_vs_family_frontier": yls - f56_atls,
                               "gap_family_vs_searched_frontier": y56 - fls_at56, "gap_search_vs_searched_frontier": yls - fls_atls,
                               "frontier_gain_at_family_courses": fls_at56 - f56_at56, "frontier_gain_at_search_courses": fls_atls - f56_atls,
                               "n_fixed_searched": int(len(fixed))})
    return {"search": search_rows, "front": front_rows}


# ---------------------------------------------------------------- summary
def _q(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    return [float(v) for v in np.percentile(x, [25, 50, 75])] if x.size else None


def analyse(search_rows, front_rows):
    """Median over worlds per endpoint and lam. Pure function, no I/O."""
    out = {"n_worlds": len({r["world"] for r in search_rows}), "by_endpoint": {}}
    for ep in ENDPOINTS:
        block = {}
        for lam in LAMS:
            s = [r for r in search_rows if r["endpoint"] == ep and r["lam"] == lam]
            f = [r for r in front_rows if r["endpoint"] == ep and r["lam"] == lam]
            e = {"mean_gain_q25_q50_q75": _q([r["mean_gain"] for r in s]), "frac_improved_q25_q50_q75": _q([r["frac_improved"] for r in s]),
                 "frac_not_converged_q25_q50_q75": _q([r["frac_not_converged"] for r in s]),
                 "courses_family_q25_q50_q75": _q([r["courses_family"] for r in s]), "courses_search_q25_q50_q75": _q([r["courses_search"] for r in s])}
            for key in ("gap_family_vs_family_frontier", "gap_search_vs_family_frontier", "gap_family_vs_searched_frontier",
                        "gap_search_vs_searched_frontier", "frontier_gain_at_family_courses", "frontier_gain_at_search_courses"):
                vals = np.array([r[key] for r in f], float)
                ok = np.isfinite(vals)
                e[key] = {"q25_q50_q75": _q(vals), "n_defined": int(ok.sum()), "n_above": int((vals[ok] > 0.05).sum()), "n_below": int((vals[ok] < -0.05).sum())}
            block[f"lam{lam:g}"] = e
        g1 = block["lam0"]["mean_gain_q25_q50_q75"]
        room = block["lam0"]["gap_search_vs_searched_frontier"]["q25_q50_q75"]
        block["label_g1"] = label_g1(g1[1]) if g1 else None
        block["label_room"] = label_room(room[1]) if room else None
        out["by_endpoint"][ep] = block
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
    ap.add_argument("--config", default="mixture", help="run configuration of src/trackc_calibration.py")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    if cfg["name"] not in G_SEEDS:
        sys.exit(f"no seeds for config {cfg['name']}")
    seeds = G_SEEDS[cfg["name"]]
    OUT_DIR, D_DIR = PROJECT_ROOT / "output" / (OUT_NAME + cfg["suffix"] + ("_smoke" if args.smoke else "")), cal.stage_dir("D", cfg, args.smoke)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    if not (D_DIR / "results.json").exists():
        sys.exit(f"script 193 has not finished: {D_DIR / 'results.json'} is missing")
    worlds = read_controlled(D_DIR)
    if not worlds:
        sys.exit("script 193 found no controlled world")
    if args.smoke:
        worlds = worlds[:2]
    sched, names = rg.candidate_schedules()
    s56 = truncate(sched)
    protocol = {"script": "205_trackc_informed_search", "stage": "2G", "config": cfg["name"], "family": cfg["family"], "smoke": bool(args.smoke),
                "n_sel": N_SEL, "n_eval": N_EVAL, "selection_seed": seeds[0], "evaluation_seed": seeds[1], "tau": TAU, "nc": NC,
                "n_start": N_START, "rounds_g1": ROUNDS_G1, "rounds_g2": ROUNDS_G2, "lams": list(LAMS), "mus": list(MUS),
                "endpoints": list(ENDPOINTS), "n_family_schedules": len(sched), "n_family_after_truncation": int(len(s56)),
                "n_family_with_late_courses": int(np.asarray(sched)[:, NC:].any(axis=1).sum()), "moves": "toggle one cycle; swap neighbouring cycles",
                "scan_every": rg.SCAN_EVERY, "near_gain": NEAR_GAIN, "far_gain": FAR_GAIN, "room_limit": ROOM_LIMIT,
                "controlled_worlds": [w["world"] for w in worlds], "input": str((D_DIR / "controlled_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((D_DIR / "results.json").read_bytes()).hexdigest()}
    manifest = man.build_manifest("205_trackc_informed_search", protocol, purpose="stage2_dev", scenario=seeds, noise=None, script=SCRIPT,
                                  require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py", "src/trackc_regime.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(worlds)} controlled worlds", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        parts = list(ex.map(evaluate_world, [(w, N_SEL, N_EVAL, seeds, ROUNDS_G1, ROUNDS_G2) for w in worlds]))
    search_rows = [r for part in parts for r in part["search"]]
    front_rows = [r for part in parts for r in part["front"]]
    write_csv(search_rows, OUT_DIR, "search_by_world.csv")
    write_csv(front_rows, OUT_DIR, "frontier_by_world.csv")
    res = analyse(search_rows, front_rows)
    res.update({"script": "205_trackc_informed_search",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; how far below the best open-loop schedule the informed arm and the fixed frontier lie",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"], "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for ep, b in res["by_endpoint"].items():
        for lam in LAMS:
            e = b[f"lam{lam:g}"]
            g = e["mean_gain_q25_q50_q75"]
            print(f"{ep:9s} lam {lam:4.1f}: search gain {g[1]:+6.2f} d (q25 {g[0]:+.2f}, q75 {g[2]:+.2f}); improved {e['frac_improved_q25_q50_q75'][1]:.2f}; "
                  f"not converged {e['frac_not_converged_q25_q50_q75'][1]:.2f}; courses {e['courses_family_q25_q50_q75'][1]:.2f} -> {e['courses_search_q25_q50_q75'][1]:.2f}; "
                  f"gap family/family {e['gap_family_vs_family_frontier']['q25_q50_q75'][1]:+.1f}, search/family {e['gap_search_vs_family_frontier']['q25_q50_q75'][1]:+.1f}, "
                  f"search/searched {e['gap_search_vs_searched_frontier']['q25_q50_q75'][1]:+.1f} d")
        print(f"{ep}: {b['label_g1']}; {b['label_room']}")


if __name__ == "__main__":
    main()
