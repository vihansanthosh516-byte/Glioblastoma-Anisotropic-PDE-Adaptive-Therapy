#!/usr/bin/env python3
"""Script 119: GRAND_PLAN 14 Y1, per-patient treatment response fitted on the previous pair (MU, exploratory).

Declared in GRAND_PLAN section 14 before this script ran. Uses the later (rolling-origin) primary pairs of script 100
(pairs with a prior scan). For each pair, the previous pair of the same patient (its output scan = this pair's input
scan) is used to fit ONE kill multiplier k that scales the assumed alpha (TMZ) and beta (radiation) of
run_improved_aniso; grid k in {0, 0.25, 0.5, 1, 2, 4}; pick the best Dice on the previous pair (ties -> k closest to 1,
then smaller). The current pair is then forecast with k* and with the fixed k = 1. Model: exact solver, isotropic
tissue arm iso_homog with the (d, rho) cell selected for the patient's fold in script 100 (direction does not help,
ledger 2ak). Compared with persistence (no change). Subgroup: shrink / stable / growth from the previous pair
(cutoffs +/-10%, analysis plan subgroup rule). Development only; exploratory, labelled.
Output: output/mu_fitted_response/results.json, pairs.csv, cache/*.json
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
OUT = PROJECT_ROOT / "output" / "mu_fitted_response"
CACHE = OUT / "cache"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import run_improved_aniso as ria  # noqa: E402
from solver_monotone import TensorFKMonotone  # noqa: E402

_s = importlib.util.spec_from_file_location("s100", PROJECT_ROOT / "src" / "100_pde_manifest.py")
s100 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(s100)
ria.TensorFK = TensorFKMonotone

KS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
ARM = "iso_homog"
SEL = PROJECT_ROOT / "output" / "pde_manifest_monotone" / "selected_cells.json"


def previous_pair(pair, allp):
    m = allp[(allp["patient_id"] == pair["patient_id"]) & (allp["tp_out"] == pair["tp_in"])]
    return None if m.empty else m.sort_values("day_in").iloc[-1].to_dict()


def forecast(pair, cell, ks):
    """Dice at the output scan for each k (kill multiplier)."""
    pid = pair["patient_id"]
    tum1, tum2, brain = s100.load_inputs(pid, pair["tp_in"], pair["tp_out"], 2.0)
    if tum1.sum() == 0:
        return None, None, None
    A = ria._ATLAS
    lo = np.maximum(np.argwhere(tum1).min(0) - ria.MARGIN_VOX, 0)
    hi = np.minimum(np.argwhere(tum1).max(0) + ria.MARGIN_VOX + 1, tum1.shape)
    box = tuple(slice(int(a), int(b)) for a, b in zip(lo, hi))
    dom = (brain[box] & A["tissue"][box]) | tum1[box]
    kv = dict(x.split("=") for x in cell.split("|")[1:])
    T0 = ria.arm_tensor(A, ARM, float(kv["r"]), box)
    T0[tum1[box] & ~A["tissue"][box]] = np.array([1, 1, 1, 0, 0, 0], np.float32)
    solver = ria.TensorFK(T0 * float(kv["d"]), dom, h=2.0)
    sched = ria.build_schedule(pair.get("treatment"))
    out = {}
    for k in ks:
        u, _ = solver.run(tum1[box].astype(np.float32), [float(kv["rho"])], pair["dt_days"], kill_schedule=sched,
                          t_start=pair["day_in"], alpha=ria.ALPHA_TMZ * k, beta=ria.BETA_RT * k)
        full = np.zeros(tum1.shape, np.float32)
        full[box] = np.asarray(u)[0]
        out[k] = ria.dice(full > ria.U_VISIBLE, tum2)
    return out, ria.dice(tum1, tum2), (int(tum1.sum()), int(tum2.sum()))


def run(shard, n):
    ria._forecast_init()
    CACHE.mkdir(parents=True, exist_ok=True)
    sel = json.loads(SEL.read_text())
    allp = pd.read_csv(s100.MAN / "forecast_pairs_mu.csv")
    cohort = {p["patient_id"]: p.get("treatment_schedule")
              for p in json.loads((PROJECT_ROOT / "output" / "mu_glioma_cohort.json").read_text())}
    pairs = s100.get_pairs(2.0, stage="selected")[shard::n]
    t0 = time.time()
    for i, p in enumerate(pairs):
        f = CACHE / f"{p['patient_id']}_{p['tp_in']}_{p['tp_out']}.json"
        if f.exists():
            continue
        cell = sel[str(int(p["fold"]))][ARM]
        prev = previous_pair(p, allp)
        rec = {k: p[k] for k in ("patient_id", "tp_in", "tp_out", "dt_days", "fold")}
        if prev is None:
            rec["skip"] = "no previous pair"
        else:
            prev["treatment"] = cohort.get(p["patient_id"])
            fit, _, (pv1, pv2) = forecast(prev, cell, KS)
            if fit is None:
                rec["skip"] = "empty mask on previous pair"
            else:
                kstar = max(KS, key=lambda k: (round(fit[k], 6), -abs(np.log((k + 1e-3) / 1.0)), -k))
                cur, nc, (v1, v2) = forecast(p, cell, sorted({kstar, 1.0}))
                if cur is None:
                    rec["skip"] = "empty mask on input"
                else:
                    c = (pv2 - pv1) / max(pv1, 1)
                    rec.update({"fit_dice": fit, "k_star": kstar, "dice_fitted": cur[kstar], "dice_fixed": cur[1.0],
                                "dice_persistence": nc, "prev_change": c,
                                "subgroup": "shrink" if c < -0.10 else ("growth" if c > 0.10 else "stable")})
        f.write_text(json.dumps(rec, default=float))
        print(f"[{i + 1}/{len(pairs)}] {f.stem} {rec.get('skip', rec.get('k_star'))} {time.time() - t0:.0f}s", flush=True)


def analyze():
    from src import patient_stats as ps
    recs = [json.loads(f.read_text()) for f in sorted(CACHE.glob("*.json"))]
    df = pd.DataFrame([r for r in recs if "skip" not in r])
    df.to_csv(OUT / "pairs.csv", index=False)

    def per_patient(d, a, b):   # one value per patient (mean over that patient's pairs), as the registered analysis
        return d.assign(x=d[a] - d[b]).groupby("patient_id")["x"].mean().to_numpy()

    res = {"script": "119_mu_fitted_treatment_response", "declared": "GRAND_PLAN 14 Y1", "label": "exploratory",
           "n_pairs": int(len(df)), "n_patients": int(df["patient_id"].nunique()),
           "n_skipped": int(sum("skip" in r for r in recs)), "k_star_counts": df["k_star"].value_counts().to_dict(),
           "all": {}, "by_subgroup": {}}
    for a, b in (("dice_fitted", "dice_persistence"), ("dice_fixed", "dice_persistence"), ("dice_fitted", "dice_fixed")):
        res["all"][f"{a}_minus_{b}"] = ps.summarize_delta(per_patient(df, a, b))
    for g, d in df.groupby("subgroup"):
        res["by_subgroup"][g] = {"n_pairs": int(len(d)), "n_patients": int(d["patient_id"].nunique())}
        for a, b in (("dice_fitted", "dice_persistence"), ("dice_fitted", "dice_fixed")):
            if d["patient_id"].nunique() >= 5:
                res["by_subgroup"][g][f"{a}_minus_{b}"] = ps.summarize_delta(per_patient(d, a, b))
    res["weakest_points"] = ["Exploratory; one scalar k per patient from one previous pair.",
                             "Treatment schedules are recorded yes/no with assumed timing; alpha and beta are assumed.",
                             "Later pairs only (needs a previous pair), so fewer patients than the primary analysis."]
    (OUT / "results.json").write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: res[k] for k in ("n_pairs", "n_patients", "k_star_counts")}, indent=1, default=str))
    for k, v in res["all"].items():
        print(k, round(v["mean"], 4), [round(x, 4) for x in v["mean_ci95"]])
    for g, v in res["by_subgroup"].items():
        for k, d in v.items():
            if isinstance(d, dict):
                print(g, k, round(d["mean"], 4), [round(x, 4) for x in d["mean_ci95"]], "n", d["n"])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["run", "analyze"], required=True)
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--n-shards", type=int, default=1)
    a = ap.parse_args()
    run(a.shard, a.n_shards) if a.stage == "run" else analyze()
