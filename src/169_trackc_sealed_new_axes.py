#!/usr/bin/env python3
"""
Script 169: Track C sealed robustness test on new misspecification axes.
=======================================================================
SIMULATION ONLY. Plan: paper/research_plan_trackC_draft.md. Containment arms and the non-filter
comparators are read from output/trackc_sealed_containment/frozen_protocol.json (no retuning).
The filter keeps its own model and prior; only the truth changes. None of the axes below was
used in scripts 150-168.

PROTOCOL (fixed in code and committed before tuning or scoring; SHA-256 written to the outputs)
  Stage "tune": the ablation (filtered single threshold, wide prior, no risk rule) gets the
    same 12-config budget as each script-167 family, on the script-167 tuning draw
    (seed 20261121, n 300, noise 1_121_001/2), objective mean capped TTP. Winner frozen.
  Stage "seal": 13 suites (nominal + 6 new axes x 2 severities), n 600 each, plus a sampled
    mixture suite, n 1,500. Scenario seed 20261300 + k; noise seeds [1_300_000 + 10k + rep,
    chunk] for rep 1-2, chunks of 150 scenarios (same chunks and seeds for every arm).
    Axes (truth only):
      switching    S -> R under drug, mu in [1e-4, 1e-3] / [1e-3, 1e-2] per day (log)
      drift        drug strength m falls 2% / 5% per treated course
      missed       reading lost with prob 0.2 / 0.4 (last reading repeated)
      delayed      reading shows burden 7 / 14 days before the decision
      endpoint     true progression at 1.15 / 1.25 x B0 (controller assumes 1.2)
      start burden B0 in [0.5, 0.6] / [0.6, 0.7] (filter prior [0.05, 0.5])
    Mixture (per patient, independent): noise sd x {1, 2, 3}; reading bias {0, +10%, +20%};
      m log-uniform [3, 20]; r log-uniform [0.001, 0.10]; competition (1, 1) or (2, 0.5);
      with prob 0.5 one axis above (uniform) at one severity (uniform).
  Comparator: primary_nonfilter from script 167 (raw single threshold 0.6).
  H1: wide containment minus comparator, mean TTP, 95% lower bound > -7 days in all 13 suites.
  H2: mixture suite gain >= 14 days and 95% lower bound > 0.
  Two-stage paired bootstrap (scenarios, then replicates), 4,000 draws, seed 169.
  Everything else (front-loaded, raw on/off, ablation, Zhang, CVaR) is descriptive. Narrow
  containment is not re-scored (runtime; its weak-drug failure is shown in script 168).
  Scored once: refuses to overwrite the sealed result.

Run:  python src/169_trackc_sealed_new_axes.py tune --workers 7
      python src/169_trackc_sealed_new_axes.py seal --workers 7
"""
from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import hashlib
import json
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import trackc_sim as tc  # noqa: E402

FROZEN = PROJECT_ROOT / "output" / "trackc_sealed_containment" / "frozen_protocol.json"
OUT_DIR = PROJECT_ROOT / "output" / "trackc_sealed_new_axes"
N, N_MIX, CHUNK, N_REP, N_BOOT = 600, 1500, 150, 2, 4000
MARGIN, MIN_GAIN = -7.0, 14.0
AXES = {  # name: (moderate, severe)
    "switching": ((1e-4, 1e-3), (1e-3, 1e-2)),
    "drift": (0.02, 0.05),
    "missed": (0.2, 0.4),
    "delayed": (7.0, 14.0),
    "endpoint": (1.15, 1.25),
    "start_burden": ((0.5, 0.6), (0.6, 0.7)),
}
SUITES = ["nominal"] + [f"{a}_{s}" for a in AXES for s in ("mod", "sev")]
MIXTURE = "mixture"
ABLATION_GRID = [("filtw", a) for a in (0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.18)]
TUNE = {"seed": 20261121, "n": 300, "noise": [1_121_001, 1_121_002]}
PROTOCOL = {"axes": AXES, "suites": SUITES + [MIXTURE], "n": N, "n_mix": N_MIX, "chunk": CHUNK,
            "n_rep": N_REP, "n_boot": N_BOOT, "boot_seed": 169, "margin": MARGIN, "min_gain": MIN_GAIN,
            "ablation_grid": ABLATION_GRID, "tune": TUNE, "event_mode": "refine", "dt": 0.25,
            "seed_rule": "scenario 20261300 + k; noise [1_300_000 + 10k + rep, chunk]",
            "mixture": "sd x {1,2,3}; bias {0,.1,.2}; m logU[3,20]; r logU[.001,.1]; comp (1,1)/(2,.5); "
                       "p .5 one axis (uniform) at one severity (uniform)",
            "H1": "wide containment - primary_nonfilter, mean TTP, CI95 lower > -7 in all 13 suites",
            "H2": "mixture: diff >= 14 and CI95 lower > 0"}
PROTOCOL_HASH = hashlib.sha256(json.dumps(PROTOCOL, sort_keys=True, default=str).encode()).hexdigest()


def _apply_axis(s, idx, axis, sev, rng):
    """Set one axis at severity sev (0 mod, 1 sev) for scenarios idx (in place)."""
    v = AXES[axis][sev]
    m = len(idx)
    if axis == "switching":
        s["mu"][idx] = np.exp(rng.uniform(np.log(v[0]), np.log(v[1]), m))
    elif axis == "drift":
        s["drift"][idx] = v
    elif axis == "missed":
        s["miss_p"][idx] = v
    elif axis == "delayed":
        s["delay"][idx] = v
    elif axis == "endpoint":
        s["prog_factor"][idx] = v
    elif axis == "start_burden":
        s["B0"][idx] = rng.uniform(v[0], v[1], m)


def _blank(s, n):
    s.update(mu=np.zeros(n), drift=np.zeros(n), miss_p=np.zeros(n), delay=np.zeros(n),
             prog_factor=np.full(n, tc.PROG_FACTOR))
    return s


def build(k, name, n):
    rng = np.random.default_rng([20261300 + k, 2])
    if name != MIXTURE:
        s = _blank(tc.draw_scenarios(n, 20261300 + k), n)
        if name != "nominal":
            axis, sev = name.rsplit("_", 1)
            _apply_axis(s, np.arange(n), axis, ("mod", "sev").index(sev), rng)
        return s
    comp = rng.integers(0, 2, n)
    s = _blank(tc.draw_scenarios(n, 20261300 + k, ranges=dict(tc.RANGES, m=(3.0, 20.0, "log"),
                                                              r=(0.001, 0.10, "log")),
                                 a_sr=np.where(comp, 2.0, 1.0), a_rs=np.where(comp, 0.5, 1.0)), n)
    s["sigma_true"] = s["sigma"] * rng.choice([1.0, 2.0, 3.0], n)
    s["bias"] = np.log1p(rng.choice([0.0, 0.1, 0.2], n))
    has = rng.random(n) < 0.5
    axis = rng.integers(0, len(AXES), n)
    sev = rng.integers(0, 2, n)
    for a_i, a in enumerate(AXES):
        for v in (0, 1):
            _apply_axis(s, np.flatnonzero(has & (axis == a_i) & (sev == v)), a, v, rng)
    return s


def _job(args):
    """One arm on one chunk of one suite: (spec, k, name, n, c) -> TTP (N_REP, chunk size)."""
    spec, k, name, n, c, n_rep = args
    s = build(k, name, n)
    idx = np.arange(c * CHUNK, min(n, (c + 1) * CHUNK))
    sub = tc.subset(s, idx)
    pol = tc.make_policy(tuple(spec))
    reps = np.stack([tc.simulate(sub, pol, [1_300_000 + 10 * k + r, c])["ttp"] for r in range(1, n_rep + 1)])
    return (name, tuple(spec), c), reps


def _untreated(args):
    k, name, n = args
    return name, float(tc.simulate(build(k, name, n), tc.Periodic(10 ** 6, 10 ** 6), 0)["capped"].mean())


def _commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
    except Exception:
        return "unknown"


def tune(out_dir, workers, n_scale):
    fp = out_dir / "frozen_ablation.json"
    if fp.exists():
        sys.exit(f"refusing to overwrite {fp}")
    scen = tc.draw_scenarios(max(7, int(TUNE["n"] * n_scale)), TUNE["seed"])
    with ProcessPoolExecutor(max_workers=workers) as ex:
        reps = dict(ex.map(tc.score_spec_reps, [(s, scen, TUNE["noise"]) for s in ABLATION_GRID]))
    mean = {s: float(reps[s].mean()) for s in ABLATION_GRID}
    win = max(ABLATION_GRID, key=mean.get)
    fp.write_text(json.dumps({"script": "169_trackc_sealed_new_axes", "stage": "tune",
                              "protocol_sha256": PROTOCOL_HASH, "git_commit_at_run": _commit(),
                              "n_scale": n_scale, "tuning_mean_ttp": {str(list(s)): mean[s] for s in ABLATION_GRID},
                              "ablation_winner": list(win)}, indent=2))
    print("ablation winner", win, round(mean[win], 1))


def seal(out_dir, workers, n_scale):
    sp = out_dir / "sealed_results.json"
    if sp.exists():
        sys.exit(f"refusing to overwrite {sp}: the sealed test is scored once")
    abl = json.loads((out_dir / "frozen_ablation.json").read_text())
    if abl["protocol_sha256"] != PROTOCOL_HASH or abl["n_scale"] != n_scale:
        sys.exit("protocol or n_scale differs from the tuning stage; refusing to score")
    frozen = json.loads(FROZEN.read_text())
    arms = {"wide_containment": frozen["family_winners"]["containment_wide"],
            "comparator_raw_single": frozen["primary_nonfilter"],
            "raw_onoff": frozen["family_winners"]["raw_onoff"],
            "front_loaded": frozen["family_winners"]["periodic"],
            "ablation_filtered_threshold": abl["ablation_winner"],
            "reference_zhang": ["raw", 1.0, 0.5]}
    sizes = {name: max(7, int((N_MIX if name == MIXTURE else N) * n_scale)) for name in SUITES + [MIXTURE]}
    kk = {name: k for k, name in enumerate(SUITES + [MIXTURE])}
    t0 = time.time()
    jobs = [(spec, kk[nm], nm, sizes[nm], c, N_REP) for nm in sizes
            for c in range(-(-sizes[nm] // CHUNK)) for spec in dict.fromkeys(tuple(v) for v in arms.values())]
    jobs.sort(key=lambda j: {"risk": 0, "filtw": 1}.get(j[0][0], 2))      # slow filter jobs first
    with ProcessPoolExecutor(max_workers=workers) as ex:
        parts = dict(ex.map(_job, jobs))
        untreated = dict(ex.map(_untreated, [(kk[nm], nm, sizes[nm]) for nm in sizes]))
    reps = {(nm, tuple(spec)): np.concatenate([parts[(nm, tuple(spec), c)] for c in range(-(-sizes[nm] // CHUNK))], 1)
            for nm in sizes for spec in arms.values()}
    stats = {"mean": lambda x: float(np.mean(x)), "cvar10_scenario_mean": tc.cvar_scenario_mean,
             "cvar10_joint": tc.cvar_joint}
    out = {}
    for nm in sizes:
        base = reps[(nm, tuple(arms["comparator_raw_single"]))]
        row = {"n": sizes[nm], "untreated_capped_frac": untreated[nm]}
        for a, spec in arms.items():
            x = reps[(nm, tuple(spec))]
            r = {k: f(x) for k, f in stats.items()}
            if a != "comparator_raw_single":
                d, ci = tc.hier_bootstrap_diff(x, base, stats["mean"], n_boot=N_BOOT, seed=169)
                r["mean_diff_vs_comparator"] = {"diff": d, "ci95": ci}
            row[a] = r
        w = reps[(nm, tuple(arms["wide_containment"]))]
        for other in ("front_loaded", "ablation_filtered_threshold"):
            d, ci = tc.hier_bootstrap_diff(w, reps[(nm, tuple(arms[other]))], stats["mean"], n_boot=N_BOOT, seed=169)
            row[f"wide_minus_{other}"] = {"diff": d, "ci95": ci}
        out[nm] = row
    h1_rows = {nm: out[nm]["wide_containment"]["mean_diff_vs_comparator"] for nm in SUITES}
    mix = out[MIXTURE]["wide_containment"]["mean_diff_vs_comparator"]
    verdict = {"H1_pass": all(r["ci95"][0] > MARGIN for r in h1_rows.values()),
               "H1_failing_suites": [nm for nm, r in h1_rows.items() if r["ci95"][0] <= MARGIN],
               "H1_worst_lower_bound": min((r["ci95"][0], nm) for nm, r in h1_rows.items()),
               "H2_pass": mix["diff"] >= MIN_GAIN and mix["ci95"][0] > 0, "H2_mixture": mix}
    result = {"script": "169_trackc_sealed_new_axes",
              "label": ("SIMULATION; frozen script-167 arms; new misspecification axes + sampled mixture; "
                        "sealed, scored once. Axes and mixture weights are author choices"),
              "protocol_sha256": PROTOCOL_HASH, "protocol": PROTOCOL, "arms": arms,
              "frozen_167_sha256": frozen["protocol_sha256"], "git_commit_at_run": _commit(),
              "n_scale": n_scale, "verdict": verdict, "suites": out, "runtime_s": time.time() - t0}
    sp.write_text(json.dumps(result, indent=2, default=float))
    for nm in sizes:
        r = out[nm]
        print(f"{nm:20s} wide {r['wide_containment']['mean_diff_vs_comparator']['diff']:+7.1f} "
              f"{np.round(r['wide_containment']['mean_diff_vs_comparator']['ci95'], 1)} "
              f"ablation {r['ablation_filtered_threshold']['mean_diff_vs_comparator']['diff']:+7.1f} "
              f"front {r['front_loaded']['mean_diff_vs_comparator']['diff']:+7.1f} untreated_capped {r['untreated_capped_frac']:.3f}")
    print(json.dumps(verdict, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=("tune", "seal"))
    ap.add_argument("--workers", type=int, default=7)
    ap.add_argument("--n-scale", type=float, default=1.0)
    ap.add_argument("--out", type=str, default=str(OUT_DIR))
    args = ap.parse_args()
    out_dir = Path(args.out)
    if args.n_scale != 1.0 and out_dir == OUT_DIR:
        sys.exit("reduced-size runs must write elsewhere (--out)")
    out_dir.mkdir(parents=True, exist_ok=True)
    (tune if args.stage == "tune" else seal)(out_dir, args.workers, args.n_scale)


if __name__ == "__main__":
    main()
