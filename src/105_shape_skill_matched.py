#!/usr/bin/env python3
"""Script 105: volume-matched shape test (Amendment 6 A6.4; GRAND_PLAN #2). Exploratory.

Question: with volume error removed, does the PDE know more than the geometric rule about WHERE the core will be?
Both forecasts are cut to the TRUE next-scan volume n = |T| (an oracle on volume, so neither is a usable forecast):
  geo_matched  script 98's reshape: grow/shrink S along the signed distance inside brain | S, top-n voxels
               (identical to script 98's oracle_volume_matched).
  pde_matched  the fold's selected PDE cell (registered primary arm aniso r=1, selected_cells.json, chosen on other
               folds' first pairs) run from S over the interval; top-n voxels by PDE density u inside the same domain;
               ties (e.g. u = 0 outside the PDE box) broken by the same signed distance, so the rules differ only where
               the PDE density orders voxels differently from distance to the boundary.
Contrast (per patient mean over primary pairs): Dice(pde_matched, T) - Dice(geo_matched, T); GBM-only first, then all.
Sanity checks reported: thresholded PDE Dice vs script 100 cache; geo_matched Dice vs script 98 oracle rows.

Stages: forecast (sharded: --shard i --n-shards N; cache per pair) | analyze.
Output: output/shape_skill/cache/*.json, output/shape_skill_matched.json (+ manifest)
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
MAN = PROJECT_ROOT / "data" / "manifests"
CACHE = OUT / "shape_skill" / "cache"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import run_improved_aniso as ria  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

ARM_FAMILY = "aniso|r=1"   # registered primary arm (A1.5 / A3.1)


def _load(name, fname):
    spec = importlib.util.spec_from_file_location(name, PROJECT_ROOT / "src" / fname)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


m100 = _load("m100", "100_pde_manifest.py")
m98 = _load("m98", "98_baseline_ladder.py")


def all_primary_pairs():
    return m100.get_pairs(2.0, stage="forecast") + m100.get_pairs(2.0, stage="selected")


def pde_matched(u, sd, n):
    cand = np.flatnonzero(np.isfinite(sd).ravel())
    n = min(int(n), cand.size)
    out = np.zeros(sd.shape, bool)
    if n <= 0:
        return out
    order = np.lexsort((sd.ravel()[cand], -u.ravel()[cand]))   # by u descending, then signed distance ascending
    out.flat[cand[order[:n]]] = True
    return out


def forecast_one(pair, sel):
    pid = pair["patient_id"]
    S, T, brain = m100.load_inputs(pid, pair["tp_in"], pair["tp_out"], 2.0)
    rec = {k: pair[k] for k in ("patient_id", "tp_in", "tp_out", "dt_days", "fold", "has_prior_scan")}
    rec.update({"v_in": int(S.sum()), "v_out": int(T.sum())})
    if not S.any():
        rec["skip"] = "empty input core"
        return rec
    A = ria._ATLAS
    cell = sel[str(int(pair["fold"]))][ARM_FAMILY]
    kv = dict(x.split("=") for x in cell.split("|")[1:])
    arm = cell.split("|")[0]
    lo = np.maximum(np.argwhere(S).min(0) - ria.MARGIN_VOX, 0)
    hi = np.minimum(np.argwhere(S).max(0) + ria.MARGIN_VOX + 1, S.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    dom = (brain[box] & A["tissue"][box]) | S[box]
    T0 = ria.arm_tensor(A, arm, float(kv["r"]), box)
    T0[S[box] & ~A["tissue"][box]] = np.array([1, 1, 1, 0, 0, 0], np.float32)
    solver = ria.TensorFK(T0 * float(kv["d"]), dom, h=2.0)
    uu, _ = solver.run(S[box].astype(np.float32), [float(kv["rho"])], pair["dt_days"],
                       kill_schedule=ria.build_schedule(pair.get("treatment")), t_start=pair["day_in"])
    u = np.zeros(S.shape, np.float32)
    u[box] = uu[0]
    sd = m98.signed_domain_distance(S, brain)
    n = int(T.sum())
    rec.update({"cell": cell,
                "dice_persistence": float(ria.dice(S, T)),
                "dice_pde_threshold": float(ria.dice(u > ria.U_VISIBLE, T)),
                "dice_pde_matched": float(ria.dice(pde_matched(u, sd, n), T)),
                "dice_geo_matched": float(ria.dice(m98.reshape_to_volume(S, sd, n), T)),
                "n_voxels_u_pos": int((u > 0).sum())})
    return rec


def stage_forecast(shard, n_shards):
    CACHE.mkdir(parents=True, exist_ok=True)
    ria._forecast_init()
    sel = json.loads((OUT / "pde_manifest" / "selected_cells.json").read_text())
    pairs = all_primary_pairs()[shard::n_shards]
    print(f"shape-skill shard {shard}/{n_shards}: {len(pairs)} pairs", flush=True)
    t0 = time.time()
    for k, p in enumerate(pairs):
        f = CACHE / f"{p['patient_id']}_{p['tp_in']}_{p['tp_out']}.json"
        if f.exists():
            continue
        t1 = time.time()
        rec = forecast_one(p, sel)
        rec["seconds"] = round(time.time() - t1, 1)
        tmp = f.with_suffix(".tmp")
        tmp.write_text(json.dumps(rec, default=float))
        tmp.replace(f)
        print(f"[{k + 1}/{len(pairs)}] {p['patient_id']} {p['tp_in']}->{p['tp_out']} {rec.get('skip', '')} "
              f"{rec['seconds']}s total {time.time() - t0:.0f}s", flush=True)


def analyze():
    recs = [json.loads(f.read_text()) for f in sorted(CACHE.glob("*.json"))]
    df = pd.DataFrame([r for r in recs if "skip" not in r])
    pats = pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)
    gbm = set(pats.loc[pats["is_gbm"].astype(bool), "patient_id"]) if "is_gbm" in pats else None
    bl = pd.read_csv(OUT / "baseline_ladder_pairs.csv")
    if gbm is None:
        gbm = set(bl.loc[bl["is_gbm"].astype(bool), "patient_id"])
    res = {"script": "105_shape_skill_matched", "rule": "Amendment 6 A6.4 (exploratory)", "arm": ARM_FAMILY,
           "n_pairs": int(len(df)), "n_patients": int(df["patient_id"].nunique()), "n_skipped": int(len(recs) - len(df)),
           "n_boot": ps.DEFAULT_N_BOOT, "contrasts": {}, "sanity": {}}
    pop = {"gbm_only": df[df["patient_id"].isin(gbm)], "all_eligible": df}
    for name, d in pop.items():
        pm = d.groupby("patient_id")[["dice_pde_matched", "dice_geo_matched", "dice_persistence", "dice_pde_threshold"]].mean()
        res["contrasts"][name] = {
            "pde_matched_minus_geo_matched": ps.summarize_delta((pm["dice_pde_matched"] - pm["dice_geo_matched"]).to_numpy()),
            "pde_matched_minus_persistence": ps.summarize_delta((pm["dice_pde_matched"] - pm["dice_persistence"]).to_numpy()),
            "geo_matched_minus_persistence": ps.summarize_delta((pm["dice_geo_matched"] - pm["dice_persistence"]).to_numpy()),
            "mean_dice": {c: float(pm[c].mean()) for c in pm.columns}}
        first = d[~d["has_prior_scan"].astype(bool)].groupby("patient_id")[["dice_pde_matched", "dice_geo_matched"]].mean()
        res["contrasts"][name]["first_pairs_pde_matched_minus_geo_matched"] = ps.summarize_delta(
            (first["dice_pde_matched"] - first["dice_geo_matched"]).to_numpy())
    # sanity 1: thresholded PDE Dice vs script 100 caches (first pairs: grid cell; later pairs: dice_sel)
    diffs = []
    for r in df.itertuples(index=False):
        key = f"{r.patient_id}_{r.tp_in}_{r.tp_out}.json"
        c2 = OUT / "pde_manifest" / ("cache_sel" if r.has_prior_scan else "cache_h2") / key
        if c2.exists():
            c = json.loads(c2.read_text())
            ref = c.get("dice_sel", {}).get(ARM_FAMILY) if r.has_prior_scan else c.get("dice", {}).get("grid", {}).get(r.cell)
            if ref is not None:
                diffs.append(abs(r.dice_pde_threshold - ref))
    res["sanity"]["pde_threshold_vs_script100"] = {"n": len(diffs), "max_abs_diff": float(max(diffs)) if diffs else None}
    # sanity 2: geo_matched vs script 98 oracle rows
    orc = bl[bl["method"] == "oracle_volume_matched"].set_index(["patient_id", "tp_in", "tp_out"])["dice"]
    j = df.set_index(["patient_id", "tp_in", "tp_out"])["dice_geo_matched"]
    jj = pd.concat([j, orc], axis=1, keys=["s105", "s98"]).dropna()
    res["sanity"]["geo_matched_vs_script98_oracle"] = {"n": int(len(jj)), "max_abs_diff": float((jj["s105"] - jj["s98"]).abs().max())}
    res["weakest_points"] = [
        "Both forecasts use the true next volume, so this isolates shape skill; neither is a usable forecast.",
        "One population cell per fold; patient-specific fitting is not tested.",
        "The PDE density is computed only inside a box around the input core; outside it u = 0 and the geometric tie-break decides.",
        "2 mm grid; exploratory, not in the registered hypothesis family."]
    (OUT / "shape_skill_matched.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("shape_skill_105", OUT / "shape_skill_matched.manifest.json", script="src/105_shape_skill_matched.py",
                       seed=ps.DEFAULT_SEED, config={"arm": ARM_FAMILY, "n_boot": ps.DEFAULT_N_BOOT},
                       inputs=[MAN / "forecast_pairs_mu.csv", MAN / "split_mu.csv", OUT / "pde_manifest" / "selected_cells.json",
                               OUT / "baseline_ladder_pairs.csv", ria.ATLAS_NPZ],
                       dataset="MU-Glioma-Post", patient_split="data/manifests/split_mu.csv",
                       primary_endpoint="Dice(pde_matched) - Dice(geo_matched), per patient, GBM-only")
    for name in pop:
        c = res["contrasts"][name]["pde_matched_minus_geo_matched"]
        print(name, f"n={c['n']} mean {c['mean']:+.4f} CI {c['mean_ci95'][0]:+.4f}..{c['mean_ci95'][1]:+.4f} med {c['median']:+.4f}")
    print(json.dumps(res["sanity"]))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["forecast", "analyze"], required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    a = ap.parse_args()
    if a.stage == "forecast":
        stage_forecast(a.shard, a.n_shards)
    else:
        analyze()
