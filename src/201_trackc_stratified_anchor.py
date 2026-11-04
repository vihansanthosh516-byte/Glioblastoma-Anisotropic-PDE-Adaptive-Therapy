#!/usr/bin/env python3
"""
Script 201: Track C Stage 2, test K: does a mixture world also reproduce the MGMT-stratified PFS of RTOG 0525?
===========================================================================================================
SIMULATION ONLY. DEVELOPMENT seeds only (never a test). Exploratory. Fixed in paper/trackc_anchored_spec_2026-11-04.md
section 7 BEFORE any mixture-family result was seen. Plan: paper/trackc_restart_plan_2026-11-03.md.

Question. The mixture family (config "mixture", src/trackc_calibration.py) has two classes of patients: responders keep
the earlier kill range, non-responders have a weak one. RTOG 0525 (Gilbert 2013, PMC3816958; both arms pooled, PFS from
random assignment) reports a median PFS of 8.7 months (95% CI 6.6 to 11.2) for tumours with a methylated MGMT promoter and
5.7 months (5.1 to 6.1) for unmethylated tumours, with 245 methylated and 517 unmethylated among the 762 patients with a
result. If the responders of a world are mapped to the methylated tumours, the same simulated trial should give a
responder median inside 6.6 to 11.2 months, a non-responder median inside 5.1 to 6.1 months, and the world's responder
fraction q should lie inside the Wilson 95% interval of 245/762. The mapping is an ASSUMPTION; no data say that MGMT
status is the response class of this model.

Method. For each kept world of script 192 (config mixture), draw 3,000 NEW patients, run the standard arm (6 courses) and
the dose-dense arm (drug effect x2.0, the assumed ratio of script 193) on the same patients, pool the two arms (equal
sizes, as in the trial), and compute the interpolated Kaplan-Meier median PFS of the responders, of the non-responders
and of everybody. Flags: unmeth_ok, meth_ok, q_ok; all_ok when all three hold. The check selects nothing and tunes
nothing: it describes which kept worlds are "strictly anchored". No dosing policy is compared here.

Seeds (registry: src/trackc_seed_registry.py, purpose stage2_dev), config "mixture":
    patients   seed [20261827, world index]
    noise      key [1_890_400, world index]   (the arms ignore readings; the key is only passed to the simulator)
Written to output/trackc_stage2_K_stratified_anchor_mixture/ (manifest.json first, then strata_by_world.csv, results.json).

Run:  python src/201_trackc_stratified_anchor.py --config mixture --workers 4
      python src/201_trackc_stratified_anchor.py --config mixture --smoke
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
import trackc_survival as sv  # noqa: E402

SCRIPT = "src/201_trackc_stratified_anchor.py"
N_PAT = 3000
WEAKEST = [
    "the mapping responder = methylated is an assumption; the model's responders are a kill-strength class and nothing in the published data ties MGMT status to it",
    "the published stratified medians pool the standard and the dose-dense arm; the dose-dense arm here is the standard arm with the drug effect multiplied by an assumed 2.0, which also scales switching",
    "one published stratified pair of medians and one fraction; MGMT status was available for 91% of the samples, so the published fraction describes the tested patients only",
    "the model has no early deaths or drop-outs and no patients who progressed during chemoradiation; the trial population was selected",
    "3,000 patients per world; the stratum medians of a world near an interval edge can flip",
    "the check was fixed before the results were seen, but the mixture family itself (responders at the earlier kill range) was chosen after the single-family result; it is one robustness probe, not a model of MGMT",
]


def _smoke(args):
    global N_PAT
    if args.smoke:
        N_PAT = 400


def read_rows(path):
    with open(path, encoding="utf-8") as fh:
        return [{k: (int(v) if k == "world" else float(v)) for k, v in row.items()} for row in csv.DictReader(fh)]


def _in(x, lo_hi):
    return bool(np.isfinite(x) and lo_hi[0] <= x <= lo_hi[1])


def evaluate_world(w, n_pat, cfg):
    """One kept mixture world. Pure function; returns a flat row."""
    scen = cal.draw_patients(w, n_pat, [cfg["k_seed"], int(w["world"])])
    seed = [cfg["noise_k"], int(w["world"])]
    resp = cal.responder_mask(scen, w)
    caps = cal.std_caps(cfg["std_courses"], n_pat)
    t0, e0, _ = cal.run_arm(scen, courses=caps, seed=seed)
    t1, e1, _ = cal.run_arm(scen, courses=caps, dose=cal.DOSE_DENSE_MULT, seed=seed)
    t, e, r = np.r_[t0, t1], np.r_[e0, e1], np.r_[resp, resp]
    med = lambda m: float(sv.km_median_interp(t[m], e[m])) if m.sum() > 1 else float("nan")  # noqa: E731
    row = {"world": int(w["world"]), "q": float(w["q"]), "responder_share_sample": float(resp.mean()),
           "median_responders": med(r), "median_non_responders": med(~r), "median_pooled": med(np.ones_like(r)),
           "median_standard_arm_all": float(sv.km_median_interp(t0, e0)),
           "event_fraction_responders": float(np.mean(np.r_[e0, e1][r])) if r.any() else float("nan"),
           "event_fraction_non_responders": float(np.mean(np.r_[e0, e1][~r])) if (~r).any() else float("nan")}
    lo, hi = cal.MGMT_METHYLATED_FRACTION_CI
    row["unmeth_ok"] = _in(row["median_non_responders"], cal.MGMT_UNMETHYLATED_PFS_DAYS)
    row["meth_ok"] = _in(row["median_responders"], cal.MGMT_METHYLATED_PFS_DAYS)
    row["q_ok"] = bool(lo <= row["q"] <= hi)
    row["all_ok"] = bool(row["unmeth_ok"] and row["meth_ok"] and row["q_ok"])
    return row


def _job(args):
    return evaluate_world(*args)


def _q(x):
    x = np.asarray([v for v in x if np.isfinite(v)], float)
    return [float(v) for v in np.percentile(x, [5, 25, 50, 75, 95])] if x.size else None


def analyse(rows, controlled_ids=()):
    """Counts and quantiles over the kept worlds. Pure function, no I/O."""
    ctl = set(int(i) for i in controlled_ids)
    n = len(rows)
    res = {"n_kept_worlds": n, "n_controlled_worlds": len(ctl),
           "targets_days": {"unmethylated": list(cal.MGMT_UNMETHYLATED_PFS_DAYS), "methylated": list(cal.MGMT_METHYLATED_PFS_DAYS),
                            "methylated_fraction_wilson_95": list(cal.MGMT_METHYLATED_FRACTION_CI),
                            "published_counts_methylated_with_result": list(cal.MGMT_COUNTS)}}
    if n == 0:
        return res

    def block(sel):
        m = len(sel)
        out = {"n_worlds": m}
        if m == 0:
            return out
        for key in ("unmeth_ok", "meth_ok", "q_ok", "all_ok"):
            out[f"n_{key}"] = int(sum(bool(r[key]) for r in sel))
        out["n_unmeth_and_meth_ok"] = int(sum(bool(r["unmeth_ok"] and r["meth_ok"]) for r in sel))
        out["median_responders_q05_q25_q50_q75_q95"] = _q([r["median_responders"] for r in sel])
        out["median_non_responders_q05_q25_q50_q75_q95"] = _q([r["median_non_responders"] for r in sel])
        out["q_q05_q25_q50_q75_q95"] = _q([r["q"] for r in sel])
        return out

    res["kept"] = block(rows)
    res["controlled"] = block([r for r in rows if r["world"] in ctl])
    res["strictly_anchored_world_ids"] = sorted(r["world"] for r in rows if r["all_ok"])
    res["strictly_anchored_and_controlled_world_ids"] = sorted(r["world"] for r in rows if r["all_ok"] and r["world"] in ctl)
    return res


def write_csv(rows, out_dir, controlled_ids):
    ctl = set(int(i) for i in controlled_ids)
    cols = list(rows[0].keys()) + ["controlled"]
    with open(out_dir / "strata_by_world.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            vals = dict(r, controlled=r["world"] in ctl)
            w.writerow([f"{vals[c]:.6g}" if isinstance(vals[c], float) else int(vals[c]) if isinstance(vals[c], (bool, np.bool_)) else vals[c] for c in cols])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--config", default="mixture", help="run configuration with a mixture family (mixture, mixture_mix)")
    args = ap.parse_args()
    _smoke(args)
    cfg = cal.get_config(args.config)
    if cfg["family"] != "mixture" or cfg.get("k_seed") is None:
        sys.exit(f"config {cfg['name']!r} has no responder class, so there is no stratified check")
    OUT_DIR = cal.stage_dir("K", cfg, args.smoke)
    C_DIR, D_DIR = cal.stage_dir("C", cfg, args.smoke), cal.stage_dir("D", cfg, args.smoke)
    if (OUT_DIR / "results.json").exists():
        sys.exit(f"refusing to overwrite {OUT_DIR / 'results.json'}")
    for need in (C_DIR / "results.json", C_DIR / "kept_worlds.csv", D_DIR / "results.json"):
        if not need.exists():
            sys.exit(f"missing input {need}")
    kept = read_rows(C_DIR / "kept_worlds.csv")
    controlled_path = D_DIR / "controlled_worlds.csv"
    controlled_ids = [r["world"] for r in read_rows(controlled_path)] if controlled_path.exists() else []
    if not kept:
        sys.exit("script 192 kept no world")
    protocol = {"script": "201_trackc_stratified_anchor", "stage": "2K", "config": cfg["name"], "family": cfg["family"],
                "smoke": bool(args.smoke), "n_patients_per_world": N_PAT, "scenario_seed": cfg["k_seed"],
                "noise_base": cfg["noise_k"], "dose_dense_multiplier": cal.DOSE_DENSE_MULT,
                **cal.exposure_protocol(cfg),
                "targets_days": {"unmethylated": list(cal.MGMT_UNMETHYLATED_PFS_DAYS), "methylated": list(cal.MGMT_METHYLATED_PFS_DAYS),
                                 "methylated_fraction_wilson_95": list(cal.MGMT_METHYLATED_FRACTION_CI)},
                "kept_worlds": [k["world"] for k in kept], "controlled_worlds": controlled_ids,
                "input": str((C_DIR / "kept_worlds.csv").relative_to(PROJECT_ROOT)),
                "input_results_sha256": hashlib.sha256((C_DIR / "results.json").read_bytes()).hexdigest()}
    manifest = man.build_manifest("201_trackc_stratified_anchor", protocol, purpose="stage2_dev",
                                  scenario=(cfg["k_seed"], cfg["k_seed"]), noise=(cfg["noise_k"], cfg["noise_k"]),
                                  script=SCRIPT, require_clean=not args.smoke, check_seeds=not args.smoke,
                                  paths=man.CORE_CODE + ("src/trackc_survival.py", "src/trackc_calibration.py"))
    man.write_manifest(OUT_DIR, manifest)                  # before any simulation
    print(f"manifest written; protocol sha256 {manifest['protocol_sha256'][:12]}; git clean {manifest['git']['clean']}; "
          f"{len(kept)} kept worlds, {len(controlled_ids)} controlled; config {cfg['name']}", flush=True)

    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        rows = list(ex.map(_job, [(w, N_PAT, cfg) for w in kept]))
    write_csv(rows, OUT_DIR, controlled_ids)
    res = analyse(rows, controlled_ids)
    res.update({"script": "201_trackc_stratified_anchor",
                "label": "SIMULATION; DEVELOPMENT seeds only; exploratory; stratified anchor check, nothing selected or tuned",
                "protocol": protocol, "protocol_sha256": manifest["protocol_sha256"],
                "weakest_points": WEAKEST, "runtime_s": time.time() - t0})
    (OUT_DIR / "results.json").write_text(json.dumps(res, indent=2, default=float) + "\n", encoding="utf-8")

    print(f"runtime {res['runtime_s']:.0f} s; wrote {OUT_DIR}")
    for name in ("kept", "controlled"):
        b = res[name]
        if b.get("n_worlds"):
            print(f"{name}: {b['n_worlds']} worlds; unmethylated ok {b['n_unmeth_ok']}, methylated ok {b['n_meth_ok']}, "
                  f"q ok {b['n_q_ok']}, all three {b['n_all_ok']}")


if __name__ == "__main__":
    main()
