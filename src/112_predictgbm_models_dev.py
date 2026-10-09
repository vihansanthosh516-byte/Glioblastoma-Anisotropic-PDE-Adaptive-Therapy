#!/usr/bin/env python3
"""Script 112: H-2 growth models M1 / M2 on PREDICT-GBM development patients (Amendment 6 A6.6, Amendment 7 A7.1).

Development (TUM) only; test recurrences stay locked (predictgbm_io). Exact positivity-preserving solver (Selling).
Models (A6.6): M1 isotropic tissue model, D(x) = D (wm + gm / 10) I; M2 the same scale times the UCSF-PDGM DTI atlas
tensor normalised to unit mean diffusivity and sharpened by r in {1, 10} (run_improved_aniso.arm_tensor "aniso").
Space: PREDICT-GBM arrays are the SRI24 240 x 240 x 155 grid with an identity affine; MU / atlas arrays are the same
grid flipped in x and y (MU affine diag(-1, -1, 1); brain-mask correlation 0.969 flipped xy vs 0.638 as is; left-right
decided by the headers, brains are nearly symmetric). The atlas is therefore flipped in x and y, and the tensor
components xz and yz change sign (xy keeps its sign after two flips).
Grid: 2 mm (block mean of the 1 mm maps), box = pre-op core +/- 40 mm. Seed u0 = pre-op core (labels 1, 3).
lambda = sqrt(D / rho) in {1, 2, 4, 8, 16} mm with rho = 0.1 / day (only lambda shapes the front once the volume is
fixed). Integrate in 2-day chunks until the region u >= 0.1 holds at least the standard-plan volume (or 3,000 days).
The density u (nearest-upsampled to 1 mm) is the model score; plan = predictgbm_io.model_plan (top-k at the
standard-plan volume inside the brain, as the benchmark). Coverage of enhancing and of all recurrence is stored.
Population choice (A6.6) per arm: the lambda (and r for M2) with the best mean development coverage; done in
--stage select. Sharded: --stage run --shard i --n-shards N; per-patient cache.
Output: output/predictgbm/cache_112/*.json, output/predictgbm/models_dev.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm"
CACHE = OUT / "cache_112"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402
import run_improved_aniso as ria  # noqa: E402
from solver_monotone import TensorFKMonotone  # noqa: E402

LAMBDAS = [1.0, 2.0, 4.0, 8.0, 16.0]
RHO = 0.1
ARMS = [("M1_iso", None), ("M2_aniso_r1", 1.0), ("M2_aniso_r10", 10.0)]
BOX_MM = 40
U_FRONT = 0.1
T_MAX = 3000.0


def to2(a):
    return np.asarray(a, np.float32)[:240, :240, :154].reshape(120, 2, 120, 2, 77, 2).mean((1, 3, 5))


def up1(a):
    u = np.repeat(np.repeat(np.repeat(a, 2, 0), 2, 1), 2, 2)
    return np.concatenate([u, u[:, :, -1:]], 2)   # 154 -> 155 slices


_ATLAS = None


def atlas_pg():
    """Atlas eigen-decomposition flipped into PREDICT-GBM array order (x and y flipped, xz / yz sign change)."""
    global _ATLAS
    if _ATLAS is None:
        A = ria.load_atlas()
        v = A["v"][::-1, ::-1].copy()
        v[..., 0, :] *= -1   # x component of each eigenvector
        v[..., 1, :] *= -1   # y component: Dxy = vx vy keeps sign, Dxz and Dyz flip
        _ATLAS = {"w": A["w"][::-1, ::-1].copy(), "v": v, "tissue": A["tissue"][::-1, ::-1].copy()}
    return _ATLAS


def tensor(arm_r, scale, box):
    if arm_r is None:
        one = np.ones(scale.shape, np.float32)
        T = np.stack([one, one, one, 0 * one, 0 * one, 0 * one], -1)
    else:
        A = atlas_pg()
        T = ria.arm_tensor(A, "aniso", arm_r, box).astype(np.float32)
        md = T[..., :3].mean(-1, keepdims=True)
        T = np.where(md > 1e-6, T / np.maximum(md, 1e-6), np.array([1, 1, 1, 0, 0, 0], np.float32))
    return T * scale[..., None]


def run_case(pid):
    c = pg.load_case(pid, with_recurrence=True)   # dev only: lock refuses test ids
    core1 = pg.core(c["seg"])
    std = pg.standard_plan(c["seg"], c["brain"])
    re, ra = pg.rec_enhancing(c["rec"]), pg.rec_all(c["rec"])
    wm, gm, brain2 = to2(c["wm"]), to2(c["gm"]), to2(c["brain"]) >= 0.5
    core2 = to2(core1) >= 0.5
    if not core2.any():
        return {"pid": pid, "skip": "empty core at 2 mm"}
    lo = np.maximum(np.argwhere(core2).min(0) - BOX_MM // 2, 0)
    hi = np.minimum(np.argwhere(core2).max(0) + BOX_MM // 2 + 1, core2.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    scale = (wm + gm / 10.0)[box]
    scale[core2[box]] = np.maximum(scale[core2[box]], 1.0)
    dom = (brain2[box] & (scale > 0.05)) | core2[box]
    target2 = std.sum() / 8.0
    rec = {"pid": pid, "plan_voxels": int(std.sum()), "n_rec_enh": int(re.sum()), "n_rec_all": int(ra.sum()),
           "standard_enh": pg.coverage(re, std), "standard_all": pg.coverage(ra, std), "runs": {}}
    for arm, r in ARMS:
        base = TensorFKMonotone(tensor(r, scale, box), dom, h=2.0)   # split once; D only rescales the coefficients
        cp0, cm0, lam0 = base.cp, base.cm, base.lam_max
        for lam in LAMBDAS:
            D = lam * lam * RHO
            base.cp = [(e, c * D) for e, c in cp0]
            base.cm = [(e, c * D) for e, c in cm0]
            base.lam_max = lam0 * D
            u = core2[box].astype(np.float32)
            t, clamp = 0.0, 0.0
            while t < T_MAX and (u >= U_FRONT).sum() < target2:
                u = base.run(u, [RHO], 2.0)[0][0]
                clamp += float(base.clamp_added[0])
                t += 2.0
            full = np.zeros(core2.shape, np.float32)
            full[box] = u
            plan = pg.model_plan(up1(full), c["seg"], c["brain"], k=int(std.sum()))
            rec["runs"][f"{arm}|lam={lam:g}"] = {"enh": pg.coverage(re, plan), "all": pg.coverage(ra, plan),
                                                 "days": t, "clamp_added": clamp}
    return rec


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
    keys = sorted(recs[0]["runs"])
    mean = {k: float(np.mean([r["runs"][k]["enh"] for r in recs])) for k in keys}
    std_mean = float(np.mean([r["standard_enh"] for r in recs]))
    best = {}
    for arm, _ in ARMS:
        fam = "M2" if arm.startswith("M2") else "M1"
        cands = [k for k in keys if k.startswith(fam)]
        best[fam] = max(cands, key=lambda k: (mean[k], -float(k.split("lam=")[1])))
    res = {"script": "112_predictgbm_models_dev", "split": "development (TUM)", "n_patients": len(recs),
           "standard_enh_mean": std_mean, "mean_enh_by_cell": mean, "selected": best,
           "delta_vs_standard_enh_selected": {fam: ps.summarize_delta(np.array([r["runs"][k]["enh"] - r["standard_enh"] for r in recs]))
                                              for fam, k in best.items()},
           "note": "Development data; selection uses these same patients, so these deltas are optimistic. Test set locked."}
    (OUT / "models_dev.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("n_patients", "standard_enh_mean", "selected")}, indent=1))
    for fam, d in res["delta_vs_standard_enh_selected"].items():
        print(fam, f"{d['mean']:+.4f}", d["mean_ci95"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["run", "select"], required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    a = ap.parse_args()
    stage_run(a.shard, a.n_shards) if a.stage == "run" else stage_select()
