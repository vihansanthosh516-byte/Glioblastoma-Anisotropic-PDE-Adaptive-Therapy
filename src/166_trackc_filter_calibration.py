#!/usr/bin/env python3
"""
Script 166: Track C particle-filter calibration check (audit C1/C2, 2026-10-26).
==============================================================================
SIMULATION ONLY. Does the filter's posterior mean what it says, when the truth is drawn
(a) from the filter's own prior and (b) from outside it?

Separation of draws (all seeds new; none used before for tuning or scoring):
  development  seed 20261111 (n = 150 per suite/schedule): free to look at, may inform a
               FUTURE settings change, which would then need a new sealed seed
  sealed       seed 20261112 (n = 300): scored once with the FROZEN settings below. The
               script has no options that change the filter, and refuses to overwrite an
               existing sealed result.
Frozen filter settings = the ones used by scripts 150-158: 500 particles, Liu-West h = 0.90,
prior RANGES, Gaussian log-burden likelihood with the true sigma, resampling at ESS < P/2.

Truth suites: in_prior (RANGES); out-of-prior weak_drug (m 3-6), fast_growth (r 0.05-0.10),
high_resistance (f0 0.1-0.3). Action schedules are fixed (periodic q = 1 and q = 4), so the
check does not depend on any policy.

Per live cycle (before progression) it records, for r, m, f0, c, eps, B0, the resistant
fraction and the burden: 90% weighted equal-tailed interval coverage and width; ESS after
the new reading (before resampling); and a one-step predictive check of each new reading
(PIT under the particle mixture; 90% predictive coverage and KS distance from uniform).

Run:  python src/166_trackc_filter_calibration.py            (development + sealed)
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.stats import norm

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_sim as tc  # noqa: E402

OUT_DIR = PROJECT_ROOT / "output" / "trackc_filter_calibration"
FROZEN = {"n_particles": 500, "h": 0.90, "ranges": "trackc_sim.RANGES", "resample_ess_fraction": 0.5}
DRAWS = {"development": (20261111, 150), "sealed": (20261112, 300)}
SUITES = {"in_prior": tc.RANGES,
          "weak_drug_m3to6": dict(tc.RANGES, m=(3.0, 6.0, "log")),
          "fast_growth_r05to10": dict(tc.RANGES, r=(0.05, 0.10, "log")),
          "high_resistance_f0_10to30": dict(tc.RANGES, f0=(0.1, 0.3, "log"))}
SCHEDULES = {"q1_frontloaded": 1, "q4_paced": 4}
TARGETS = list(tc.PARAMS) + ["fr", "burden"]


def wquant(x, w, q):
    """Weighted quantiles along axis 1 (rows independent)."""
    order = np.argsort(x, 1)
    xs = np.take_along_axis(x, order, 1)
    cw = np.cumsum(np.take_along_axis(w, order, 1), 1)
    rows = np.arange(x.shape[0])
    return [xs[rows, np.minimum(np.argmax(cw >= qq, 1), x.shape[1] - 1)] for qq in q]


def run_block(args):
    """One (suite, schedule, draw) block: true dynamics + filter, per-cycle records."""
    suite, sched, draw, seed, n = args
    k = list(SUITES).index(suite) * 10 + list(SCHEDULES).index(sched)
    scen = tc.draw_scenarios(n, seed + 1000 * k, ranges=SUITES[suite])
    rng = np.random.default_rng(seed + 1000 * k + 1)
    q = SCHEDULES[sched]
    S, R = (1 - scen["f0"]) * scen["B0"], scen["f0"] * scen["B0"]
    level = tc.PROG_FACTOR * scen["B0"]
    ttp = np.full(n, np.nan)
    noise = rng.standard_normal((n, tc.N_CYCLES)) * scen["sigma"][:, None]
    pf = tc.ParticleFilter(scen["sigma"], FROZEN["n_particles"], rng, h=FROZEN["h"])
    rec = {t: {"hit": [], "width": []} for t in TARGETS}
    rec.update(ess=[], pit=[], cycle=[])
    used = np.zeros(n, int)
    for j in range(tc.N_CYCLES):
        live = np.isnan(ttp)
        if not live.any():
            break
        y = np.log(S + R) + noise[:, j]
        mu_pred = np.log(np.maximum(pf.S + pf.R, 1e-300))           # particles propagated to cycle j
        w_prior = pf.w
        if j > 0:
            pit = (w_prior * norm.cdf((y[:, None] - mu_pred) / pf.sigma)).sum(1)
            rec["pit"].append(pit[live])
        lw = pf.logw - 0.5 * ((y[:, None] - mu_pred) / pf.sigma) ** 2
        ww = np.exp(lw - lw.max(1, keepdims=True))
        ww /= ww.sum(1, keepdims=True)
        rec["ess"].append((1.0 / (ww ** 2).sum(1))[live])
        pf.update(j, y)
        w = pf.w
        truth = {name: scen[name] for name in tc.PARAMS}
        truth["fr"], truth["burden"] = R / (S + R), S + R
        est = {name: pf.theta[name] for name in tc.PARAMS}
        est["fr"], est["burden"] = pf.fr(), pf.S + pf.R
        for t in TARGETS:
            lo, hi = wquant(est[t], w, (0.05, 0.95))
            rec[t]["hit"].append(((truth[t] >= lo) & (truth[t] <= hi))[live])
            rec[t]["width"].append((hi - lo)[live])
        rec["cycle"].append(np.full(int(live.sum()), j))
        treat = np.full(n, j % q == 0) & (used < tc.BUDGET) & live
        used += treat
        S, R = tc.run_cycle(S, R, scen, treat, j * tc.CYCLE_DAYS, prog_level=level, ttp=ttp)
        pf.propagate(j, treat)
    cyc = np.concatenate(rec["cycle"])
    out = {"n_scenarios": n, "n_live_cycles": int(cyc.size), "targets": {}}
    for t in TARGETS:
        hit, width = np.concatenate(rec[t]["hit"]), np.concatenate(rec[t]["width"])
        lo_p, hi_p, kind = SUITES["in_prior"].get(t, (None, None, None))
        prior_w = None
        if kind == "log":
            prior_w = float(np.exp(np.log(lo_p) + 0.95 * np.log(hi_p / lo_p)) - np.exp(np.log(lo_p) + 0.05 * np.log(hi_p / lo_p)))
        elif kind == "lin":
            prior_w = 0.9 * (hi_p - lo_p)
        bins = {"cycles_0_5": cyc < 6, "cycles_6_17": (cyc >= 6) & (cyc < 18), "cycles_18_plus": cyc >= 18}
        out["targets"][t] = {"coverage90": float(hit.mean()), "median_width": float(np.median(width)),
                             "median_width_over_prior90_width": (float(np.median(width)) / prior_w) if prior_w else None,
                             "coverage90_by_cycle": {b: (float(hit[m].mean()) if m.any() else None) for b, m in bins.items()}}
    ess = np.concatenate(rec["ess"])
    pit = np.concatenate(rec["pit"]) if rec["pit"] else np.array([])
    ks = float(np.max(np.abs(np.sort(pit) - (np.arange(1, pit.size + 1) - 0.5) / pit.size))) if pit.size else None
    out["ess_after_reading"] = {"median": float(np.median(ess)), "p10": float(np.quantile(ess, 0.1)),
                                "frac_below_P_over_10": float(np.mean(ess < FROZEN["n_particles"] / 10))}
    out["predictive_check"] = {"n": int(pit.size),
                               "coverage90": float(np.mean((pit >= 0.05) & (pit <= 0.95))) if pit.size else None,
                               "ks_from_uniform": ks,
                               "frac_pit_below_0.05": float(np.mean(pit < 0.05)) if pit.size else None,
                               "frac_pit_above_0.95": float(np.mean(pit > 0.95)) if pit.size else None}
    return (suite, sched, draw), out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--draws", nargs="+", default=["development", "sealed"], choices=list(DRAWS))
    ap.add_argument("--n-scale", type=float, default=1.0, help="smoke runs only; < 1 shrinks n")
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    sealed_path = out_dir / "sealed_results.json"
    if "sealed" in args.draws and sealed_path.exists():
        sys.exit(f"refusing to overwrite {sealed_path}: the sealed draw is scored once")
    if "sealed" in args.draws and args.n_scale != 1.0 and Path(args.out) == OUT_DIR:
        sys.exit("the sealed draw is only scored at full size")
    t0 = time.time()
    jobs = [(s, sc, d, DRAWS[d][0], max(5, int(DRAWS[d][1] * args.n_scale)))
            for d in args.draws for s in SUITES for sc in SCHEDULES]
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        res = dict(ex.map(run_block, jobs))
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        commit = "unknown"
    for d in args.draws:
        block = {f"{s}/{sc}": res[(s, sc, d)] for s in SUITES for sc in SCHEDULES}
        result = {"script": "166_trackc_filter_calibration", "draw": d,
                  "label": ("SIMULATION filter calibration; frozen settings; out-of-prior suites are failure "
                            "tests by design" + ("; SEALED, scored once" if d == "sealed" else "; development draw")),
                  "git_commit_at_run": commit, "frozen_filter_settings": FROZEN,
                  "seed": DRAWS[d][0], "n_scale": args.n_scale,
                  "suites": {k: {kk: list(vv) for kk, vv in v.items()} for k, v in SUITES.items()},
                  "results": block, "runtime_s": time.time() - t0}
        (out_dir / f"{d}_results.json").write_text(json.dumps(result, indent=2, default=float))
        for key, b in block.items():
            print(d, key, "fr cov", round(b["targets"]["fr"]["coverage90"], 3),
                  "m cov", round(b["targets"]["m"]["coverage90"], 3),
                  "pred cov", b["predictive_check"]["coverage90"] and round(b["predictive_check"]["coverage90"], 3),
                  "ESS med", round(b["ess_after_reading"]["median"], 1))


if __name__ == "__main__":
    main()
