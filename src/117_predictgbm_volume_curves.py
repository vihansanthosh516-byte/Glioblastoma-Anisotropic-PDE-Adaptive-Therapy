#!/usr/bin/env python3
"""Script 117: GRAND_PLAN 13 X1 (barrier-aware geodesic plan) and X2 (coverage-volume curves), PREDICT-GBM development.

Declared in GRAND_PLAN section 13 before this script ran. Development (TUM) only; test recurrences stay locked.
Score maps (higher = treat first), each turned into a plan by predictgbm_io.topk_plan inside the brain at
k = fraction x standard-plan volume, fractions {0.5, ..., 1.0}:
  euclid   : minus Euclidean distance from the core (fraction 1.0 = the standard 15 mm plan, up to ties)
  geodesic : minus geodesic distance from the core inside the brain, paths blocked where CSF probability > 0.5
             (ESTRO-EANO 2023 barriers: ventricles, falx / tentorium spaces, sulci); 26-neighbour Dijkstra at 1 mm in a
             box of core +/- 40 mm; unreachable voxels rank last
  gliodil, sbtc (LOTI), pinngbm : released prediction maps, resampled and clipped as the benchmark
  unet_raw : released U-Net prediction, ranked by its raw value (no clip; the clip destroys its ranking, ledger 2am)
Endpoints: mean enhancing coverage per (map, fraction); X2 = smallest fraction whose mean coverage >= euclid at 1.0;
X1 = geodesic - euclid at fraction 1.0 (equal volume). Patient bootstrap from patient_stats.
Output: output/predictgbm/volume_curves.json, volume_curves_pairs.csv
"""
from __future__ import annotations

import json
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "predictgbm"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402

FRACS = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
BOX = 40
CSF_BLOCK = 0.5
MAPS = ["euclid", "geodesic", "gliodil", "sbtc", "pinngbm", "unet_raw"]


def geodesic(core, passable):
    """Shortest-path distance (mm) from the core through passable voxels; 26-neighbour graph in a box."""
    lo = np.maximum(np.argwhere(core).min(0) - BOX, 0)
    hi = np.minimum(np.argwhere(core).max(0) + BOX + 1, core.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    P, C = passable[box] | core[box], core[box]
    idx = -np.ones(P.shape, np.int64)
    idx[P] = np.arange(P.sum())
    rows, cols, w = [], [], []
    for dx in (-1, 0, 1):
        for dy in (-1, 0, 1):
            for dz in (-1, 0, 1):
                if (dx, dy, dz) <= (0, 0, 0):
                    continue
                a = idx[max(dx, 0):P.shape[0] + min(dx, 0), max(dy, 0):P.shape[1] + min(dy, 0), max(dz, 0):P.shape[2] + min(dz, 0)]
                b = idx[max(-dx, 0):P.shape[0] + min(-dx, 0), max(-dy, 0):P.shape[1] + min(-dy, 0), max(-dz, 0):P.shape[2] + min(-dz, 0)]
                m = (a >= 0) & (b >= 0)
                rows.append(a[m]); cols.append(b[m]); w.append(np.full(m.sum(), np.sqrt(dx * dx + dy * dy + dz * dz)))
    r, c, ww = np.concatenate(rows), np.concatenate(cols), np.concatenate(w)
    n = int(P.sum())
    src = n   # one super-source joined to every core voxel at cost 0 (tiny epsilon)
    cid = idx[C]
    r = np.concatenate([r, np.full(cid.size, src)]); c = np.concatenate([c, cid]); ww = np.concatenate([ww, np.full(cid.size, 1e-9)])
    G = coo_matrix((ww, (r, c)), shape=(n + 1, n + 1)).tocsr()
    d = dijkstra(G, directed=False, indices=src)[:n]
    out = np.full(core.shape, np.inf, np.float32)
    sub = np.full(P.shape, np.inf, np.float32)
    sub[P] = d
    out[box] = sub
    return out


def score_maps(pid, c):
    core = pg.core(c["seg"])
    eu = distance_transform_edt(~core)
    passable = c["brain"] & (c["csf"] <= CSF_BLOCK)
    gd = geodesic(core, passable)
    gd = np.where(np.isfinite(gd), gd, 1e6 + eu)   # unreachable: last, ordered by Euclidean distance
    maps = {"euclid": -eu, "geodesic": -gd}
    for m in ("gliodil", "sbtc", "pinngbm"):
        if (pg.DATA / pid / f"{m}_pred.nii.gz").exists():
            maps[m] = np.clip(pg.resample_to(pg._load(pid, f"{m}_pred.nii.gz").astype(np.float32), c["seg"].shape), 0, 1)
    if (pg.DATA / pid / "unet_pred.nii.gz").exists():
        maps["unet_raw"] = pg.resample_to(pg._load(pid, "unet_pred.nii.gz").astype(np.float32), c["seg"].shape)
    return maps


def rank_plan(score, k, brain):
    """Top-k by score inside the brain; scores may be negative (distance maps), so no fade fallback here."""
    cand = np.flatnonzero(brain)
    s = score.flat[cand]
    out = np.zeros(score.shape, bool)
    out.flat[cand[np.argsort(-s, kind="stable")[:k]]] = True
    return out


def main():
    from src import patient_stats as ps
    rows = []
    ids = pg.dev_ids()
    for i, pid in enumerate(ids):
        c = pg.load_case(pid, with_recurrence=True)
        std = pg.standard_plan(c["seg"], c["brain"])
        re = pg.rec_enhancing(c["rec"])
        if not pg.core(c["seg"]).any() or re.sum() == 0:
            continue
        maps = score_maps(pid, c)
        r = {"pid": pid, "plan_voxels": int(std.sum()), "standard_enh": pg.coverage(re, std)}
        for name, s in maps.items():
            for f in FRACS:
                k = int(round(f * std.sum()))
                plan = rank_plan(s, k, c["brain"]) if name in ("euclid", "geodesic", "unet_raw") else pg.topk_plan(s, k, c["brain"])
                r[f"{name}@{f:g}"] = pg.coverage(re, plan)
        rows.append(r)
        if (i + 1) % 10 == 0:
            print(f"{i + 1}/{len(ids)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "volume_curves_pairs.csv", index=False)
    base = df["euclid@1"]
    res = {"script": "117_predictgbm_volume_curves", "split": "development (TUM)", "declared": "GRAND_PLAN 13 X1, X2",
           "n_patients": int(len(df)), "check_euclid1_vs_standard_max_abs": float((df["euclid@1"] - df["standard_enh"]).abs().max()),
           "mean": {}, "x2_smallest_fraction_matching_standard": {}, "delta_vs_euclid1": {}}
    for name in MAPS:
        cols = [f"{name}@{f:g}" for f in FRACS if f"{name}@{f:g}" in df]
        if not cols:
            continue
        sub = df.dropna(subset=cols)
        res["mean"][name] = {col.split("@")[1]: float(sub[col].mean()) for col in cols}
        ok = [float(col.split("@")[1]) for col in cols if sub[col].mean() >= sub["euclid@1"].mean()]
        res["x2_smallest_fraction_matching_standard"][name] = {"fraction": min(ok) if ok else None, "n": int(len(sub))}
        res["delta_vs_euclid1"][name] = {col.split("@")[1]: ps.summarize_delta((sub[col] - sub["euclid@1"]).to_numpy())
                                         for col in cols}
    res["x1_geodesic_minus_euclid_equal_volume"] = res["delta_vs_euclid1"]["geodesic"]["1"]
    res["weakest_points"] = ["Development patients; published models may have seen TUM.",
                             "CSF > 0.5 is a proxy for ESTRO barriers (no falx / tentorium segmentation).",
                             "Smallest-fraction endpoint uses mean coverage on a 0.1 grid; uncertainty in delta_vs_euclid1."]
    (OUT / "volume_curves.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("n_patients", "check_euclid1_vs_standard_max_abs", "mean",
                                          "x2_smallest_fraction_matching_standard")}, indent=1, default=float))


if __name__ == "__main__":
    main()
