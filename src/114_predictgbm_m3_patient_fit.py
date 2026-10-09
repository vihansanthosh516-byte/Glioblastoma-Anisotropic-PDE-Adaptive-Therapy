#!/usr/bin/env python3
"""Script 114: H-2 model M3, per-patient fit from the pre-op scan only (Amendment 7 A7.2, thresholds A7.2a).

Development (TUM) only; test recurrences stay locked (predictgbm_io). Same tissue model, atlas, grid and exact solver as
script 112. Per patient: seed u = 1 in the 2 mm voxel at the centre of mass of the pre-op core; for each arm
(iso, aniso r1, aniso r10) and lambda grow until the volume u >= 0.6 (T1Gd) reaches the core volume; fit score =
Dice(u >= 0.6, core) + Dice(u >= 0.35, core + edema). The fit never reads the recurrence. Per arm family the best
lambda (and r) is chosen per patient (ties -> smaller lambda). Then growth continues until u >= 0.1 holds the
standard-plan volume (as script 112) and the plan is top-k at that volume. Coverage is stored for every cell so the
selection can be audited; only the fit-selected cell is the M3 result.
Output: output/predictgbm/cache_114/*.json, output/predictgbm/m3_dev.json
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm"
CACHE = OUT / "cache_114"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
from solver_monotone import TensorFKMonotone  # noqa: E402

_spec = importlib.util.spec_from_file_location("s112", PROJECT_ROOT / "src" / "112_predictgbm_models_dev.py")
s112 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(s112)

T_T1, T_FLAIR = 0.6, 0.35   # A7.2a
LAMBDAS, RHO, ARMS, T_MAX, U_FRONT = s112.LAMBDAS, s112.RHO, s112.ARMS, s112.T_MAX, s112.U_FRONT
BOX_MM = 60   # seed is one voxel, so the front travels further than in script 112


def dice(a, b):
    s = a.sum() + b.sum()
    return 1.0 if s == 0 else float(2 * (a & b).sum() / s)


def run_case(pid):
    c = pg.load_case(pid, with_recurrence=True)   # recurrence used for scoring only, after the fit
    core1 = pg.core(c["seg"])
    std = pg.standard_plan(c["seg"], c["brain"])
    re, ra = pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])
    wm, gm, brain2 = s112.to2(c["wm"]), s112.to2(c["gm"]), s112.to2(c["brain"]) >= 0.5
    core2 = s112.to2(core1) >= 0.5
    vis2 = s112.to2(np.isin(c["seg"], (1, 2, 3))) >= 0.5
    if not core2.any():
        return {"pid": pid, "skip": "empty core at 2 mm"}
    lo = np.maximum(np.argwhere(core2).min(0) - BOX_MM // 2, 0)
    hi = np.minimum(np.argwhere(core2).max(0) + BOX_MM // 2 + 1, core2.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    scale = (wm + gm / 10.0)[box]
    scale[core2[box]] = np.maximum(scale[core2[box]], 1.0)
    dom = (brain2[box] & (scale > 0.05)) | core2[box]
    cb, vb = core2[box], vis2[box]
    pts = np.argwhere(cb)
    seed = tuple(pts[np.argmin(((pts - pts.mean(0)) ** 2).sum(1))])   # core voxel nearest the centre of mass
    ncore, target2 = int(cb.sum()), std.sum() / 8.0
    rec = {"pid": pid, "plan_voxels": int(std.sum()), "n_rec_enh": int(re.sum()), "n_rec_all": int(ra.sum()),
           "standard_enh": pg.coverage(re, std), "standard_all": pg.coverage(ra, std), "runs": {}}
    for arm, r in ARMS:
        base = TensorFKMonotone(s112.tensor(r, scale, box), dom, h=2.0)
        cp0, cm0, lam0 = base.cp, base.cm, base.lam_max
        for lam in LAMBDAS:
            D = lam * lam * RHO
            base.cp = [(e, k * D) for e, k in cp0]
            base.cm = [(e, k * D) for e, k in cm0]
            base.lam_max = lam0 * D
            u = np.zeros(cb.shape, np.float32)
            u[seed] = 1.0
            t = 0.0
            while t < T_MAX and (u >= T_T1).sum() < ncore:
                u = base.run(u, [RHO], 2.0)[0][0]
                t += 2.0
            fit = dice(u >= T_T1, cb) + dice(u >= T_FLAIR, vb)
            t_fit = t
            while t < T_MAX and (u >= U_FRONT).sum() < target2:
                u = base.run(u, [RHO], 2.0)[0][0]
                t += 2.0
            full = np.zeros(core2.shape, np.float32)
            full[box] = u
            plan = pg.model_plan(s112.up1(full), c["seg"], c["brain"], k=int(std.sum()))
            rec["runs"][f"{arm}|lam={lam:g}"] = {"fit": fit, "t_fit": t_fit, "days": t,
                                                 "enh": pg.coverage(re, plan), "all": pg.coverage(ra, plan)}
    return rec


def pick(runs, fam):
    cands = [k for k in runs if k.startswith(fam)]
    return max(cands, key=lambda k: (runs[k]["fit"], -float(k.split("lam=")[1])))


def stage_run(shard, n):
    CACHE.mkdir(parents=True, exist_ok=True)
    ids = pg.dev_ids()[shard::n]
    t0 = time.time()
    for k, pid in enumerate(ids):
        f = CACHE / f"{pid}.json"
        if f.exists():
            continue
        t1 = time.time()
        r = run_case(pid)
        r["seconds"] = round(time.time() - t1, 1)
        f.with_suffix(".tmp").write_text(json.dumps(r))
        f.with_suffix(".tmp").replace(f)
        print(f"[{k + 1}/{len(ids)}] {pid} {r.get('skip', '')} {r['seconds']}s total {time.time() - t0:.0f}s", flush=True)


def stage_select():
    from src import patient_stats as ps
    recs = [json.loads(f.read_text()) for f in sorted(CACHE.glob("*.json"))]
    recs = [r for r in recs if "skip" not in r and r["n_rec_enh"] > 0]
    res = {"script": "114_predictgbm_m3_patient_fit", "split": "development (TUM)", "n_patients": len(recs),
           "thresholds": {"t1gd": T_T1, "flair": T_FLAIR, "source": "arXiv 2311.16536 sec 2.2 (A7.2a)"},
           "standard_enh_mean": float(np.mean([r["standard_enh"] for r in recs])), "families": {}}
    for fam in ("M3i_iso", "M3a_aniso"):
        prefix = "M1" if fam == "M3i_iso" else "M2"
        sel = [r["runs"][pick(r["runs"], prefix)] for r in recs]
        res["families"][fam] = {
            "enh_mean": float(np.mean([s["enh"] for s in sel])), "all_mean": float(np.mean([s["all"] for s in sel])),
            "fit_mean": float(np.mean([s["fit"] for s in sel])),
            "delta_vs_standard_enh": ps.summarize_delta(np.array([s["enh"] - r["standard_enh"] for s, r in zip(sel, recs)]))}
    a = [r["runs"][pick(r["runs"], "M2")]["enh"] for r in recs]
    i = [r["runs"][pick(r["runs"], "M1")]["enh"] for r in recs]
    res["M3a_minus_M3i_enh"] = ps.summarize_delta(np.array(a) - np.array(i))
    res["note"] = ("Development data. The fit uses the pre-op scan only, so the per-patient choice does not see the "
                   "recurrence; but the grid and thresholds were set by the plan, not tuned here. Test set locked.")
    (OUT / "m3_dev.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(res, indent=1, default=float)[:2000])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["run", "select"], required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    a = ap.parse_args()
    stage_run(a.shard, a.n_shards) if a.stage == "run" else stage_select()
