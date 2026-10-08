#!/usr/bin/env python3
"""Script 100: PDE forecasts on the manifest pairs and split, extended grid, clamp log (Phase 2.3; plan A1.1, A1.5).

Differences from script 81 (which this replaces for the plan):
  * pairs   : data/manifests/forecast_pairs_mu.csv, in_primary only (14-365 d, non-empty core input), rolling-origin
              (every consecutive pair), not only scan 1 -> 2
  * folds   : data/manifests/split_mu.csv (patient-level, hash based)
  * grid    : EXTENDED, to test the grid-edge problem (W3): rho {0, .01, .03, .1, .2}, d {.003, .01, .03, .1, .3},
              sharpening r {1, 10}
  * logging : mass added by the lower clamp / removed by the upper clamp for every run (W21)
  * target  : core (labels 1, 3), 2 mm, as script 81
Arms: aniso (r = 1, 10), iso_same, iso_homog, plus aniso_global_z (homogeneous diagonal tensor diag(1, 1, r),
r = 10, normalised to unit trace/3, a direction that is NOT patient- or tract-informed; negative control for H3/H4).
Mode --h 1: pilot at 1 mm (inputs at native 1 mm, atlas repeated 2x by nearest neighbour), reduced grid,
stratified first pairs (A4.4). Answers whether the h = 2 mm optimum moves (W22).

Stages: forecast | analyze | pilot_analyze.  Output dir: output/pde_manifest/{cache_h2, cache_h1}, results.json
Selection (analyze): per arm and fold, the (r, d, rho) cell with the best mean Dice over TRAINING patients (patient
means of the pair Dice); applied to the test fold. All reported Dice are out-of-fold.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output" / "pde_manifest"
MAN = PROJECT_ROOT / "data" / "manifests"
sys.path.insert(0, str(PROJECT_ROOT))

import run_improved_aniso as ria  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

RHOS = [0.0, 0.01, 0.03, 0.1, 0.2]
DS = [0.003, 0.01, 0.03, 0.1, 0.3]
SHARPEN = [1.0, 10.0]
GLOBAL_R = 10.0
CORE = (1, 3)
PILOT_RHOS = [0.03, 0.1, 0.2]
PILOT_DS = [0.01, 0.03, 0.1]
PILOT_ARMS = [("aniso", 1.0), ("iso_homog", 1.0)]
PILOT_PER_STRATUM = 3   # 3 x 3 strata -> at most 27 first pairs, all folds (Amendment 4.4)


def arms_full():
    return ([("aniso", r) for r in SHARPEN] + [("iso_same", 1.0), ("iso_homog", 1.0), ("aniso_global_z", GLOBAL_R)])


def global_z_tensor(shape, r):
    w = np.array([1.0, 1.0, r], np.float32)
    w = w * (3.0 / w.sum())                   # trace 3 -> mean diffusivity 1, as the atlas arms
    return np.tile(np.array([w[0], w[1], w[2], 0, 0, 0], np.float32), tuple(shape) + (1,))


def mask_native(pid, tp):
    p = ria.MU_DIR / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_tumorMask.nii.gz"
    return np.isin(np.asarray(nib.load(str(p)).dataobj), CORE)[:240, :240, :154]


def brain_native(pid, tp):
    p = ria.MU_DIR / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_brain_t1n.nii.gz"
    return (np.asarray(nib.load(str(p)).dataobj) > 0)[:240, :240, :154]


def up2(a):
    for ax in range(3):
        a = np.repeat(a, 2, axis=ax)
    return a


def load_inputs(pid, tp_in, tp_out, h):
    # 2 mm: block-average of the 1 mm binary mask, >= 0.5 (identical to script 81)
    if h == 2.0:
        raw1 = np.isin(np.asarray(nib.load(str(ria.mu_paths(pid, tp_in)[0])).dataobj), CORE)
        raw2 = np.isin(np.asarray(nib.load(str(ria.mu_paths(pid, tp_out)[0])).dataobj), CORE)
        t1 = ria.to_2mm(raw1) >= 0.5
        t2 = ria.to_2mm(raw2) >= 0.5
        brain = ria.to_2mm(np.asarray(nib.load(str(ria.mu_paths(pid, tp_in)[1])).dataobj) > 0) >= 0.5
        return t1, t2, brain
    return mask_native(pid, tp_in), mask_native(pid, tp_out), brain_native(pid, tp_in)


def forecast_pair(pair, atlas, h, grid):
    pid = pair["patient_id"]
    tum1, tum2, brain = load_inputs(pid, pair["tp_in"], pair["tp_out"], h)
    rec = {k: pair[k] for k in ("patient_id", "tp_in", "tp_out", "dt_days", "day_in")}
    rec.update({"h_mm": h, "v1_voxels": int(tum1.sum()), "v2_voxels": int(tum2.sum())})
    if tum1.sum() == 0:
        rec["skip"] = "empty target mask at input scan (cannot seed)"
        return rec
    A = atlas
    gshape = tum1.shape
    if h == 1.0:
        A = {"w": up2(atlas["w"]), "v": up2(atlas["v"]), "md": up2(atlas["md"]), "tissue": up2(atlas["tissue"]),
             "md_ref": atlas["md_ref"]}
    margin = ria.MARGIN_VOX * (2 if h == 1.0 else 1)
    lo = np.maximum(np.argwhere(tum1).min(0) - margin, 0)
    hi = np.minimum(np.argwhere(tum1).max(0) + margin + 1, gshape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    dom = (brain[box] & A["tissue"][box]) | tum1[box]
    u0 = tum1[box].astype(np.float32)
    schedule = ria.build_schedule(pair.get("treatment"))
    rec["dice"] = {"no_change": ria.dice(tum1, tum2), "grid": {}}
    rec["clamp"] = {}
    steps = 0

    def tensor_for(arm, r):
        if arm == "aniso_global_z":
            return global_z_tensor(tum1[box].shape, r)
        T0 = ria.arm_tensor(A, arm, r, box)
        T0[tum1[box] & ~A["tissue"][box]] = np.array([1, 1, 1, 0, 0, 0], np.float32)
        return T0

    def run(arm, r, dval, rhos):
        solver = ria.TensorFK(tensor_for(arm, r) * dval, dom, h=h)
        u, n = solver.run(u0, rhos, pair["dt_days"], kill_schedule=schedule, t_start=pair["day_in"])
        outs = []
        for i, rho in enumerate(rhos):
            full = np.zeros(gshape, np.float32)
            full[box] = u[i]
            outs.append((rho, ria.dice(full > ria.U_VISIBLE, tum2),
                         [float(solver.clamp_added[i]), float(solver.clamp_removed[i]), float(u0.sum())]))
        return outs, n

    if grid.get("selected") is not None:
        # later (rolling-origin) pairs: only the cell chosen for this patient's fold on first-pair training data
        sel = grid["selected"][str(int(pair["fold"]))]
        rec["dice_sel"], rec["clamp_sel"] = {}, {}
        for fam, cell in sel.items():
            kv = dict(x.split("=") for x in cell.split("|")[1:])
            arm = cell.split("|")[0]
            outs, n = run(arm, float(kv["r"]), float(kv["d"]), [float(kv["rho"])])
            rec["dice_sel"][fam] = outs[0][1]
            rec["clamp_sel"][fam] = outs[0][2]
            steps += n
    else:
        for arm, r in grid["arms"]:
            for dval in grid["ds"]:
                outs, n = run(arm, r, dval, grid["rhos"])
                steps += n
                for rho, dv, cl in outs:
                    key = f"{arm}|r={r:g}|d={dval:g}|rho={rho:g}"
                    rec["dice"]["grid"][key] = dv
                    rec["clamp"][key] = cl
    rec["solver_steps"] = steps
    return rec


_GRID = None
_H = 2.0
_CACHE = None


def _init(h, grid, cache):
    global _GRID, _H, _CACHE
    ria._forecast_init()
    _GRID, _H, _CACHE = grid, h, cache


def worker(pair):
    key = f"{pair['patient_id']}_{pair['tp_in']}_{pair['tp_out']}"
    f = _CACHE / f"{key}.json"
    if f.exists():
        rec = json.loads(f.read_text())
        if rec.get("atlas_sha") == ria._ATLAS_SHA and rec.get("dt_days") == pair["dt_days"]:
            return rec
    t0 = time.time()
    rec = forecast_pair(pair, ria._ATLAS, _H, _GRID)
    rec["seconds"] = round(time.time() - t0, 1)
    rec["atlas_sha"] = ria._ATLAS_SHA
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.replace(f)
    return rec


def get_pairs(h, limit=None, stage="forecast"):
    """forecast / pilot: first primary pair of every patient (no earlier scan). selected: all later primary pairs."""
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")
    pairs = pairs[pairs["in_primary"]].copy()
    cohort = {p["patient_id"]: p.get("treatment_schedule")
              for p in json.loads((PROJECT_ROOT / "output" / "mu_glioma_cohort.json").read_text())}
    split = pd.read_csv(MAN / "split_mu.csv").set_index("patient_id")["fold"]
    pairs["fold"] = pairs["patient_id"].map(split)
    pairs = pairs[pairs["has_prior_scan"]] if stage == "selected" else pairs[~pairs["has_prior_scan"]]
    if h == 1.0:   # pilot (Amendment 4.4): stratified subset of first pairs, all folds. Strata = input-core-volume tertile x
        # interval tertile (both known before the forecast); up to PILOT_PER_STRATUM per stratum, ordered by sha256(patient_id).
        import hashlib
        pairs = pairs.copy()
        pairs["vq"] = pd.qcut(pairs["core_vox_in"], 3, labels=False, duplicates="drop")
        pairs["dq"] = pd.qcut(pairs["dt_days"], 3, labels=False, duplicates="drop")
        pairs["hk"] = pairs["patient_id"].map(lambda s: hashlib.sha256(("pilot|" + s).encode()).hexdigest())
        pairs = pairs.sort_values("hk").groupby(["vq", "dq"], group_keys=False).head(PILOT_PER_STRATUM)
    pairs = pairs.sort_values(["patient_id", "day_in"])
    if limit:
        pairs = pairs.head(limit)
    out = pairs.to_dict("records")
    for p in out:
        p["treatment"] = cohort.get(p["patient_id"])
    return out


def stage_forecast(h, limit, shard, n_shards, stage):
    """Sequential loop over this shard's pairs. Parallelism = several processes started with --shard i --n-shards N
    (multiprocessing.spawn is refused on this Windows setup: PermissionError WinError 5)."""
    selected = None
    if stage == "selected":
        cache = OUT / "cache_sel"
        selected = json.loads((OUT / "selected_cells.json").read_text())
    else:
        cache = OUT / f"cache_h{int(h)}"
    cache.mkdir(parents=True, exist_ok=True)
    arms = arms_full() if h == 2.0 else PILOT_ARMS
    ds = DS if h == 2.0 else PILOT_DS
    rhos = RHOS if h == 2.0 else PILOT_RHOS
    grid = {"arms": arms, "ds": ds, "rhos": rhos, "selected": selected}
    pairs = get_pairs(h, limit, stage)[shard::n_shards]
    what = "selected cells only" if selected else f"{len(arms)} arms x {len(ds)} d x {len(rhos)} rho"
    print(f"{stage} shard {shard}/{n_shards}: {len(pairs)} pairs, h={h}, {what}", flush=True)
    _init(h, grid, cache)
    t0 = time.time()
    for k, pair in enumerate(pairs):
        rec = worker(pair)
        print(f"[{k + 1}/{len(pairs)}] {rec['patient_id']} {rec['tp_in']}->{rec['tp_out']} {rec.get('skip', '')} "
              f"{rec.get('seconds')}s total {time.time() - t0:.0f}s", flush=True)


def load_records(cache_name):
    recs = []
    for f in sorted((OUT / cache_name).glob("*.json")):
        r = json.loads(f.read_text())
        if "skip" not in r:
            recs.append(r)
    return recs


def family_of(cell):
    arm = cell.split("|")[0]
    return f"aniso|{cell.split('|')[1]}" if arm == "aniso" else arm


def select_cells():
    """Per fold and arm family: the (r, d, rho) cell with the best mean Dice over TRAINING patients, using first pairs.
    Patient means of the pair Dice; ties go to the first cell in sorted order."""
    recs = load_records("cache_h2")
    split = pd.read_csv(MAN / "split_mu.csv").set_index("patient_id")["fold"]
    all_cells = sorted(next(iter(recs))["dice"]["grid"].keys())
    df = pd.DataFrame([{"patient_id": r["patient_id"], **r["dice"]["grid"]} for r in recs])
    pm = df.groupby("patient_id")[all_cells].mean()
    pm["fold"] = pm.index.map(split)
    fam_cells = {}
    for c in all_cells:
        fam_cells.setdefault(family_of(c), []).append(c)
    sel = {str(f): {fam: pm[pm["fold"] != f][cl].mean().idxmax() for fam, cl in fam_cells.items()} for f in range(5)}
    (OUT / "selected_cells.json").write_text(json.dumps(sel, indent=1))
    print(json.dumps(sel, indent=1))
    return sel


def analyze():
    recs = load_records("cache_h2")
    recs_sel = load_records("cache_sel") if (OUT / "cache_sel").exists() else []
    split = pd.read_csv(MAN / "split_mu.csv").set_index("patient_id")["fold"]
    pats = pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)
    gbm = pats[pats["dataset"] == "MU"].set_index("patient_id")["is_gbm"]
    sel = json.loads((OUT / "selected_cells.json").read_text()) if (OUT / "selected_cells.json").exists() else select_cells()
    fams = sorted(sel["0"])
    rows = []
    for r in recs:
        f = int(split[r["patient_id"]])
        rows.append({"patient_id": r["patient_id"], "tp_in": r["tp_in"], "tp_out": r["tp_out"], "fold": f, "dt_days": r["dt_days"],
                     "has_prior_scan": False, "no_change": r["dice"]["no_change"],
                     **{fam: r["dice"]["grid"][sel[str(f)][fam]] for fam in fams}})
    for r in recs_sel:
        f = int(split[r["patient_id"]])
        rows.append({"patient_id": r["patient_id"], "tp_in": r["tp_in"], "tp_out": r["tp_out"], "fold": f, "dt_days": r["dt_days"],
                     "has_prior_scan": True, "no_change": r["dice"]["no_change"], **{fam: r["dice_sel"][fam] for fam in fams}})
    df2 = pd.DataFrame(rows)
    lad = pd.read_csv(PROJECT_ROOT / "output" / "baseline_ladder_pairs.csv")
    g = lad[lad["method"] == "geometric_train_rate"][["patient_id", "tp_in", "tp_out", "dice"]].rename(columns={"dice": "geometric"})
    df2 = df2.merge(g, on=["patient_id", "tp_in", "tp_out"], how="left")
    df2["is_gbm"] = df2["patient_id"].map(gbm)
    df2.to_csv(OUT / "oof_pairs.csv", index=False)

    edge = {fam: [] for fam in fams}
    for f in range(5):
        for fam in fams:
            kv = dict(x.split("=") for x in sel[str(f)][fam].split("|")[1:])
            rho, d = float(kv["rho"]), float(kv["d"])
            edge[fam].append({"cell": sel[str(f)][fam], "rho_edge": rho in (RHOS[0], RHOS[-1]), "rho_top": rho == RHOS[-1],
                              "d_bottom": d == DS[0], "d_top": d == DS[-1]})

    def delta(sub, a, b):
        pa, pb = ps.patient_means(sub, a), ps.patient_means(sub, b)
        j = pd.concat([pa, pb], axis=1, keys=["a", "b"]).dropna()
        return (j["a"] - j["b"]).to_numpy()

    res = {"script": "100_pde_manifest", "n_boot": ps.DEFAULT_N_BOOT, "n_pairs": int(len(df2)), "n_first_pairs": int((~df2["has_prior_scan"]).sum()),
           "n_later_pairs": int(df2["has_prior_scan"].sum()), "n_patients": int(df2["patient_id"].nunique()),
           "grid": {"rhos": RHOS, "ds": DS, "sharpen": SHARPEN, "global_r": GLOBAL_R},
           "selection": "per fold and arm, best mean Dice over training patients' FIRST pairs; later pairs run at that cell only",
           "selected_cells_per_fold": sel, "grid_edge": edge, "contrasts": {}}
    res["grid_edge_summary"] = {fam: {"folds_rho_at_top": int(sum(x["rho_top"] for x in e)),
                                      "folds_rho_at_edge": int(sum(x["rho_edge"] for x in e)),
                                      "folds_d_at_bottom": int(sum(x["d_bottom"] for x in e)),
                                      "folds_d_at_top": int(sum(x["d_top"] for x in e))} for fam, e in edge.items()}
    strata = {"all_pairs_primary": df2, "gbm_only": df2[df2["is_gbm"] == True],  # noqa: E712
              "first_pair_only": df2[~df2["has_prior_scan"]]}
    contrasts = ([(f"{fam}_vs_persistence", fam, "no_change") for fam in fams]
                 + [(f"{fam}_vs_geometric_train_rate", fam, "geometric") for fam in fams]
                 + [("aniso|r=1_vs_iso_same", "aniso|r=1", "iso_same"), ("aniso|r=10_vs_iso_same", "aniso|r=10", "iso_same"),
                    ("aniso|r=1_vs_iso_homog", "aniso|r=1", "iso_homog"),
                    ("aniso|r=10_vs_aniso_global_z", "aniso|r=10", "aniso_global_z"),
                    ("iso_homog_vs_iso_same", "iso_homog", "iso_same")])
    for sname, sub in strata.items():
        blk = {"n_patients": int(sub["patient_id"].nunique())}
        for name, a, b in contrasts:
            d = delta(sub, a, b)
            if len(d) >= 5:
                blk[name] = ps.summarize_delta(d)
        res["contrasts"][sname] = blk
    cl = []
    for r in recs:
        f = int(split[r["patient_id"]])
        for fam in fams:
            a, rm, m0 = r["clamp"][sel[str(f)][fam]]
            cl.append({"fam": fam, "added_frac": a / m0, "removed_frac": rm / m0})
    cl = pd.DataFrame(cl)
    res["clamp_selected_cells"] = {fam: {"median_added_frac": float(x["added_frac"].median()),
                                         "p95_added_frac": float(x["added_frac"].quantile(.95)),
                                         "max_added_frac": float(x["added_frac"].max())} for fam, x in cl.groupby("fam")}
    allc = [v[0] / v[2] for r in recs for v in r["clamp"].values()]
    res["clamp_all_cells"] = {"median_added_frac": float(np.median(allc)), "p95": float(np.percentile(allc, 95)),
                              "max": float(np.max(allc))}
    (OUT / "results.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("pde_manifest_100", OUT / "results.manifest.json", script="src/100_pde_manifest.py", seed=ps.DEFAULT_SEED,
                       config={"grid": res["grid"], "families": fams},
                       inputs=[MAN / "forecast_pairs_mu.csv", MAN / "split_mu.csv", ria.ATLAS_NPZ],
                       dataset="MU-Glioma-Post core", patient_split="data/manifests/split_mu.csv",
                       primary_endpoint="PDE Dice delta vs persistence (per patient)")
    print(json.dumps(res["grid_edge_summary"], indent=1))
    for sname, b in res["contrasts"].items():
        for k in ("aniso|r=1_vs_persistence", "aniso|r=1_vs_geometric_train_rate", "iso_homog_vs_persistence",
                  "iso_homog_vs_geometric_train_rate"):
            if k in b:
                d = b[k]
                print(sname, k, f"n={b['n_patients']} mean {d['mean']:+.3f} CI {d['mean_ci95'][0]:+.3f}..{d['mean_ci95'][1]:+.3f} "
                                f"med {d['median']:+.3f}")


def pilot_analyze():
    r1 = load_records("cache_h1")
    r2 = {f"{r['patient_id']}_{r['tp_in']}_{r['tp_out']}": r for r in load_records("cache_h2")}
    rows = []
    for r in r1:
        k = f"{r['patient_id']}_{r['tp_in']}_{r['tp_out']}"
        if k not in r2:
            continue
        for c, v1 in r["dice"]["grid"].items():
            v2 = r2[k]["dice"]["grid"].get(c)
            if v2 is None:
                continue
            rows.append({"pair": k, "cell": c, "h1": v1, "h2": v2, "clamp_h1": r["clamp"][c][0] / r["clamp"][c][2]})
    df = pd.DataFrame(rows)
    m = df.groupby("cell")[["h1", "h2"]].mean()
    m["diff_h2_minus_h1"] = m["h2"] - m["h1"]
    out = {"n_pairs": int(df["pair"].nunique()), "mean_dice_by_cell": m.round(4).to_dict(orient="index"),
           "best_cell_h1": m["h1"].idxmax(), "best_cell_h2": m["h2"].idxmax(),
           "best_dice_h1": float(m["h1"].max()), "best_dice_h2": float(m["h2"].max()),
           "max_abs_cell_diff": float(m["diff_h2_minus_h1"].abs().max()),
           "mean_abs_cell_diff": float(m["diff_h2_minus_h1"].abs().mean()),
           "cell_rank_spearman": float(m["h1"].rank().corr(m["h2"].rank())),
           "no_change_mean": float(np.mean([r["dice"]["no_change"] for r in r1]))}
    # the cell most folds chose at 2 mm, as a fixed reference: PDE minus persistence at each resolution, per pair
    ref = {}
    for cell in ("aniso|r=1|d=0.01|rho=0.1", "iso_homog|r=1|d=0.01|rho=0.1"):
        d1, d2, nc1, nc2 = [], [], [], []
        for r in r1:
            k = f"{r['patient_id']}_{r['tp_in']}_{r['tp_out']}"
            if k in r2 and cell in r["dice"]["grid"] and cell in r2[k]["dice"]["grid"]:
                d1.append(r["dice"]["grid"][cell] - r["dice"]["no_change"])
                d2.append(r2[k]["dice"]["grid"][cell] - r2[k]["dice"]["no_change"])
                nc1.append(r["dice"]["no_change"]); nc2.append(r2[k]["dice"]["no_change"])
        if d1:
            ref[cell] = {"n": len(d1), "delta_vs_persistence_h1_mean": float(np.mean(d1)), "delta_vs_persistence_h2_mean": float(np.mean(d2)),
                         "per_pair_diff_h1_minus_h2_mean": float(np.mean(np.array(d1) - np.array(d2))),
                         "per_pair_diff_abs_max": float(np.max(np.abs(np.array(d1) - np.array(d2)))),
                         "persistence_dice_h1_mean": float(np.mean(nc1)), "persistence_dice_h2_mean": float(np.mean(nc2))}
    out["reference_cell_resolution_effect"] = ref
    (OUT / "pilot_h1_vs_h2.json").write_text(json.dumps(out, indent=1, default=float))
    print(json.dumps({k: v for k, v in out.items() if k != "mean_dice_by_cell"}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["forecast", "select", "selected", "analyze", "pilot_analyze"], required=True)
    ap.add_argument("--h", type=float, default=2.0)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage in ("forecast", "selected"):
        stage_forecast(a.h, a.limit, a.shard, a.n_shards, a.stage)
    elif a.stage == "select":
        select_cells()
    elif a.stage == "analyze":
        analyze()
    else:
        pilot_analyze()
