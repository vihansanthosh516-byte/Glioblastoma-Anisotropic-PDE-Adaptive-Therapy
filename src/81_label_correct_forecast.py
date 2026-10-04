#!/usr/bin/env python3
"""
Script 81: Forecast scored on the correct tumour target (Track B negative 1, measurement validity).
=================================================================================================
run_improved_aniso.py scores "whole tumour = mask > 0". MU-Glioma-Post masks use BraTS-GLI
post-treatment labels: 1 NETC (non-enhancing/necrotic core), 2 SNFH (FLAIR signal / edema),
3 ET (enhancing tumour), 4 RC (resection cavity). mask > 0 therefore scores the surgical
cavity (which collapses after surgery) and edema (which changes with steroids/radiation)
as "tumour". A tumour-cell growth model should be scored on the cellular tumour.
Probe (volumes only, 152 pairs, scan index 0 -> 1): grew / shrank-or-same
   mask > 0      47% / 53%  median log ratio -0.04
   labels 1,3    61% / 39%  median log ratio +0.50   (core; 130 pairs with core > 0 at both scans)
   label 4       28% / 72%  median log ratio -0.51   (cavity)
That probe was read BEFORE this design was written; the target definition below is justified
by the label definitions, not by choosing the best-looking row, and every definition is reported.

PRE-SPECIFIED DESIGN (fixed before the first run; not changed after seeing data)
-------------------------------------------------------------------------------
Pipeline  identical to run_improved_aniso.forecast_patient: scan index 0 -> 1 pairs (152),
          same DTI atlas, same parameter grid (rho, d, r), same treatment kill term, same
          5-fold patient-level CV selection of (rho, d, r), same u > 0.5 visibility, all Dice
          out-of-fold. ONLY the tumour masks (initial condition u0 and score target) change.
Target    --target core (PRIMARY): labels {1, 3}  (cellular tumour core)
          --target wt   (SECONDARY): labels {1, 2, 3} (whole tumour without the cavity)
          Patients with an empty target mask at scan 1 cannot be seeded and are listed as
          skipped. An empty target at scan 2 scores Dice 0 for every arm (kept, intention to
          forecast).
Arms      anisotropic, iso_same, iso_homog, no_change (as in run_improved_aniso).
Primary endpoints, target core (Holm-corrected over the two):
          P1  out-of-fold Dice, anisotropic vs no_change
          P2  out-of-fold Dice, anisotropic vs iso_same
Also reported (descriptive; not in the Holm family): iso_homog vs no_change, grown vs
          shrank-or-same subgroups, the subset with core > 0 at both scans, and the committed
          mask > 0 result next to it for the same patients.
Stages    forecast (resumable per-patient cache), analyze.

Output: output/forecast_labels_<target>/{forecast_grid.json, results.json}
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import nibabel as nib
import numpy as np
from scipy import ndimage

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
import run_improved_aniso as ria  # noqa: E402

LABELS = {"core": (1, 3), "wt": (1, 2, 3)}
TARGET = "core"


def paths(target: str):
    out = PROJECT_ROOT / "output" / f"forecast_labels_{target}"
    return out, out / "cache", out / "forecast_grid.json", out / "results.json"


def label_mask(path, labels):
    return ria.to_2mm(np.isin(np.asarray(nib.load(str(path)).dataobj), labels)) >= 0.5


def forecast_patient(pair, atlas, labels):
    """run_improved_aniso.forecast_patient with the tumour masks taken from `labels`."""
    m1p, t1p = ria.mu_paths(pair["patient_id"], pair["tp1"])
    m2p, _ = ria.mu_paths(pair["patient_id"], pair["tp2"])
    tum1, tum2 = label_mask(m1p, labels), label_mask(m2p, labels)
    brain = ria.to_2mm(np.asarray(nib.load(str(t1p)).dataobj) > 0) >= 0.5
    rec = {**pair, "v1_voxels": int(tum1.sum()), "v2_voxels": int(tum2.sum()), "grew": bool(tum2.sum() > tum1.sum())}
    if tum1.sum() == 0:
        rec["skip"] = "empty target mask at scan 1 (cannot seed)"
        return rec
    rec["dice"] = {"no_change": ria.dice(tum1, tum2), "grid": {}}
    lo = np.maximum(np.argwhere(tum1).min(0) - ria.MARGIN_VOX, 0)
    hi = np.minimum(np.argwhere(tum1).max(0) + ria.MARGIN_VOX + 1, ria.GRID2)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    dom = (brain[box] & atlas["tissue"][box]) | tum1[box]
    u0 = tum1[box].astype(np.float32)
    schedule = ria.build_schedule(pair.get("treatment"))
    steps = 0
    for arm, r in [("aniso", r_) for r_ in ria.SHARPEN] + [("iso_same", 1.0), ("iso_homog", 1.0)]:
        T0 = ria.arm_tensor(atlas, arm, r, box)
        T0[tum1[box] & ~atlas["tissue"][box]] = np.array([1, 1, 1, 0, 0, 0], np.float32)
        for dval in ria.DS:
            solver = ria.TensorFK(T0 * dval, dom)
            u, n = solver.run(u0, ria.RHOS, pair["dt_days"], kill_schedule=schedule, t_start=pair["t1_day"])
            steps += n
            for i, rho in enumerate(ria.RHOS):
                full = np.zeros(ria.GRID2, np.float32)
                full[box] = u[i]
                rec["dice"]["grid"][f"{arm}|r={r:g}|d={dval:g}|rho={rho:g}"] = ria.dice(full > ria.U_VISIBLE, tum2)
    rec["solver_steps"] = steps
    return rec


def worker(args):
    pair, target = args
    _, cache_dir, _, _ = paths(target)
    cache = cache_dir / f"forecast_{pair['patient_id']}.json"
    if cache.exists():
        rec = json.loads(cache.read_text())
        if rec.get("atlas_sha") == ria._ATLAS_SHA and rec.get("dt_days") == pair["dt_days"]:
            return rec
    t0 = time.time()
    rec = forecast_patient(pair, ria._ATLAS, LABELS[target])
    rec["seconds"] = round(time.time() - t0, 1)
    rec["atlas_sha"] = ria._ATLAS_SHA
    tmp = cache.with_suffix(".tmp")
    tmp.write_text(json.dumps(rec))
    tmp.replace(cache)
    return rec


def stage_forecast(target: str, limit: int | None) -> None:
    import multiprocessing as mp
    out, cache_dir, grid_json, _ = paths(target)
    cache_dir.mkdir(parents=True, exist_ok=True)
    pairs, excluded = ria.forecast_pairs()
    if limit:
        pairs = pairs[:limit]
    recs = {}
    with mp.get_context("spawn").Pool(ria.FORECAST_WORKERS, initializer=ria._forecast_init) as pool:
        for k, rec in enumerate(pool.imap_unordered(worker, [(p, target) for p in pairs])):
            recs[rec["patient_id"]] = rec
            print(f"[{k + 1}/{len(pairs)}] {rec['patient_id']} {rec.get('skip', '')} "
                  f"no_change={rec.get('dice', {}).get('no_change')} {rec['seconds']}s", flush=True)
    grid_json.write_text(json.dumps({"target": target, "labels": LABELS[target], "rhos": ria.RHOS,
                                     "ds_mm2_per_day": ria.DS, "sharpen": ria.SHARPEN, "u_visible": ria.U_VISIBLE,
                                     "atlas_sha": ria._atlas_sha(), "excluded": excluded,
                                     "patients": [recs[p["patient_id"]] for p in pairs]}, indent=1))


def analyze(target: str) -> None:
    _, _, grid_json, results_json = paths(target)
    g = json.loads(grid_json.read_text())
    recs = [r for r in g["patients"] if "skip" not in r]
    skipped = [(r["patient_id"], r["skip"]) for r in g["patients"] if "skip" in r]
    grew = np.array([r["grew"] for r in recs])
    both = np.array([r["v1_voxels"] > 0 and r["v2_voxels"] > 0 for r in recs])
    oof, chosen = ria.cv_select(recs, "dice")
    oof["no_change"] = np.array([r["dice"]["no_change"] for r in recs])
    arms = ["anisotropic", "iso_same", "iso_homog", "no_change"]

    def summ(x):
        return {"mean": float(np.mean(x)), "mean_ci95": ria.boot_ci(x), "median": float(np.median(x)), "n": int(len(x))}

    def block(mask):
        return {a: summ(oof[a][mask]) for a in arms} if mask.sum() > 0 else None

    def tests(mask):
        if mask.sum() <= 5:
            return None
        return {"anisotropic_vs_no_change": ria.paired(oof["anisotropic"][mask], oof["no_change"][mask]),
                "anisotropic_vs_iso_same": ria.paired(oof["anisotropic"][mask], oof["iso_same"][mask]),
                "iso_homog_vs_no_change": ria.paired(oof["iso_homog"][mask], oof["no_change"][mask])}

    allm = np.ones(len(recs), bool)
    prim = tests(allm)
    adj = ria.holm({k: prim[k]["wilcoxon_p"] for k in ("anisotropic_vs_no_change", "anisotropic_vs_iso_same")})
    for k, v in adj.items():
        prim[k]["wilcoxon_p_holm"] = v
    res = {"script": "81_label_correct_forecast", "target": target, "labels": LABELS[target],
           "n_scored": len(recs), "n_skipped": len(skipped), "skipped": skipped, "n_grew": int(grew.sum()),
           "n_shrank_or_same": int((~grew).sum()), "n_target_nonempty_both_scans": int(both.sum()),
           "interval_days_median": float(np.median([r["dt_days"] for r in recs])),
           "dice_out_of_fold_all": block(allm), "primary_tests_all": prim,
           "subset_target_nonempty_both_scans": {"dice": block(both), "tests": tests(both)},
           "descriptive_grew": {"dice": block(grew), "tests": tests(grew)},
           "descriptive_shrank_or_same": {"dice": block(~grew), "tests": tests(~grew)},
           "chosen_parameters_per_fold": chosen}
    results_json.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "chosen_parameters_per_fold"}, indent=2))
    print(f"[saved] {results_json}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("stage", choices=["forecast", "analyze"])
    ap.add_argument("--target", choices=list(LABELS), default="core")
    ap.add_argument("--limit", type=int, default=None)
    a = ap.parse_args()
    stage_forecast(a.target, a.limit) if a.stage == "forecast" else analyze(a.target)


if __name__ == "__main__":
    main()
