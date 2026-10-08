#!/usr/bin/env python3
"""Script 103: how much of the forecast effect is within segmentation / registration noise? (Phase 3.3; plan A1.x, review round 2 item 7)

Declared before the run. The noise model is an ASSUMPTION: a 1-voxel (2 mm) boundary error is plausible for a glioma core
mask from two scans, but it is not measured on these data.
Perturbations of a 2 mm core mask (all seeded by sha256(patient|tp|draw)):
  erode1   binary erosion by 1 voxel (6-connectivity)         deterministic
  dilate1  binary dilation by 1 voxel                          deterministic
  jitter   every boundary voxel (inside or outside) flips with p = 0.25
  shift    translation by 1 voxel (2 mm) along one random axis and sign (registration error)
Experiment A (ceiling): Dice(S, perturb(S)) per pair, patient mean. This is the Dice that a perfect 'no change' truth would
  score if only measurement noise separated the two scans.
Experiment B (stability of the effect): K = 5 draws. In each draw the INPUT mask and the TARGET mask are perturbed
  independently by shift then jitter; persistence and the geometric rule (volume V' * exp(g dt), g from baseline_ladder.json,
  computed from the training folds of the unperturbed data) are re-scored against the perturbed target. Report the
  patient-level mean delta (geometric - persistence) per draw, its mean and SD over draws, for GBM-only primary pairs and for
  all primary pairs. Question: is the unperturbed delta (+0.0149 GBM-only, baseline_ladder.json) larger than what noise of this
  size does to it? The PDE is not re-run here (cost); its delta is within 0.002 of the geometric rule (results.json), so the
  conclusion carries over only as an approximation. This is stated as a limit.
Input : output/baseline_ladder.json, data/manifests/*.csv, MU masks   Output: output/segmentation_perturbation.json
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path as _Path

import numpy as np
import pandas as pd
from scipy import ndimage

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
MAN = PROJECT_ROOT / "data" / "manifests"
sys.path.insert(0, str(PROJECT_ROOT))

from src import forecast_metrics as fm  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

_spec = importlib.util.spec_from_file_location("b98", PROJECT_ROOT / "src" / "98_baseline_ladder.py")
b98 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(b98)

K_DRAWS = 5
P_FLIP = 0.25
STRUCT = ndimage.generate_binary_structure(3, 1)


def rng_for(*parts):
    h = hashlib.sha256("|".join(str(p) for p in parts).encode()).hexdigest()
    return np.random.default_rng(int(h[:16], 16))


def erode1(m):
    return ndimage.binary_erosion(m, STRUCT)


def dilate1(m):
    return ndimage.binary_dilation(m, STRUCT)


def jitter(m, rng):
    ring = dilate1(m) & ~erode1(m)
    flip = ring & (rng.random(m.shape) < P_FLIP)
    return m ^ flip


def shift1(m, rng):
    ax, sg = int(rng.integers(3)), int(rng.choice([-1, 1]))
    out = np.zeros_like(m)
    src = [slice(None)] * 3
    dst = [slice(None)] * 3
    if sg > 0:
        src[ax], dst[ax] = slice(0, -1), slice(1, None)
    else:
        src[ax], dst[ax] = slice(1, None), slice(0, -1)
    out[tuple(dst)] = m[tuple(src)]
    return out


def dice(a, b):
    return fm.dice(a, b)


def main():
    t0 = time.time()
    lad = json.loads((OUT / "baseline_ladder.json").read_text())
    g_train = {int(k): v for k, v in lad["g_train_per_day"].items()}
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")
    pats = pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)
    pats = pats[pats["dataset"] == "MU"].set_index("patient_id")
    split = pd.read_csv(MAN / "split_mu.csv").set_index("patient_id")["fold"]
    pairs["fold"] = pairs["patient_id"].map(split)
    pairs["is_gbm"] = pairs["patient_id"].map(pats["is_gbm"]).astype(bool)
    pairs = pairs[pairs["in_primary"]].reset_index(drop=True)
    rows_a, rows_b = [], []
    for i, r in enumerate(pairs.itertuples(index=False)):
        pid, f = r.patient_id, int(r.fold)
        S = b98.core_mask(pid, r.tp_in)
        T = b98.core_mask(pid, r.tp_out)
        if not S.any():
            continue
        brain = b98.brain_mask(pid, r.tp_in)
        a = {"patient_id": pid, "tp_in": r.tp_in, "tp_out": r.tp_out, "is_gbm": r.is_gbm,
             "erode1": dice(S, erode1(S)), "dilate1": dice(S, dilate1(S)),
             "jitter": dice(S, jitter(S, rng_for(pid, r.tp_in, "A", "j"))),
             "shift": dice(S, shift1(S, rng_for(pid, r.tp_in, "A", "s")))}
        rows_a.append(a)
        for k in range(K_DRAWS):
            S2 = jitter(shift1(S, rng_for(pid, r.tp_in, k, "S1")), rng_for(pid, r.tp_in, k, "S2"))
            T2 = jitter(shift1(T, rng_for(pid, r.tp_out, k, "T1")), rng_for(pid, r.tp_out, k, "T2"))
            if not S2.any():
                S2 = S
            sd = b98.signed_domain_distance(S2, brain)
            F = b98.reshape_to_volume(S2, sd, float(S2.sum()) * np.exp(g_train[f] * r.dt_days))
            dp, dg = dice(S2, T2), dice(F, T2)
            rows_b.append({"patient_id": pid, "tp_in": r.tp_in, "tp_out": r.tp_out, "is_gbm": r.is_gbm, "draw": k,
                           "persistence": dp, "geometric": dg, "delta": dg - dp})
        if i % 40 == 0:
            print(f"[{time.time() - t0:.0f}s] pair {i}/{len(pairs)}", flush=True)
    A, B = pd.DataFrame(rows_a), pd.DataFrame(rows_b)
    B.to_csv(OUT / "segmentation_perturbation_pairs.csv", index=False)
    res = {"script": "103_segmentation_perturbation", "assumption": "1-voxel (2 mm) boundary and registration error, not measured",
           "n_pairs": int(len(A)), "p_flip": P_FLIP, "k_draws": K_DRAWS, "A_ceiling_patient_mean_dice": {}, "B_delta_by_draw": {}}
    for pop, sub in (("gbm_only", A[A["is_gbm"]]), ("all", A)):
        res["A_ceiling_patient_mean_dice"][pop] = {c: float(ps.patient_means(sub, c).mean()) for c in ("erode1", "dilate1", "jitter", "shift")}
    for pop, sub in (("gbm_only", B[B["is_gbm"]]), ("all", B)):
        per_draw = []
        for k in range(K_DRAWS):
            d = sub[sub["draw"] == k]
            pm = ps.patient_means(d, "delta")
            s = ps.summarize_delta(pm.to_numpy(), n_boot=2000)
            per_draw.append({"draw": k, "mean_delta": s["mean"], "ci95": s["mean_ci95"], "median": s["median"],
                             "persistence_mean_dice": float(ps.patient_means(d, "persistence").mean())})
        md = np.array([x["mean_delta"] for x in per_draw])
        res["B_delta_by_draw"][pop] = {"draws": per_draw, "mean_over_draws": float(md.mean()), "sd_over_draws": float(md.std(ddof=1)),
                                       "n_draws_ci_excludes_zero": int(sum(x["ci95"][0] > 0 for x in per_draw))}
    (OUT / "segmentation_perturbation.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("segmentation_perturbation_103", OUT / "segmentation_perturbation.manifest.json",
                       script="src/103_segmentation_perturbation.py", seed=0, config={"k_draws": K_DRAWS, "p_flip": P_FLIP},
                       inputs=[OUT / "baseline_ladder.json", MAN / "forecast_pairs_mu.csv", MAN / "split_mu.csv"],
                       dataset="MU-Glioma-Post core", patient_split="data/manifests/split_mu.csv", primary_endpoint="n/a (robustness)")
    print(json.dumps({k: v for k, v in res.items() if k in ("A_ceiling_patient_mean_dice",)}, indent=1))
    for pop in res["B_delta_by_draw"]:
        b = res["B_delta_by_draw"][pop]
        print(pop, "mean delta over draws", round(b["mean_over_draws"], 4), "SD", round(b["sd_over_draws"], 4), "draws CI>0:", b["n_draws_ci_excludes_zero"])


if __name__ == "__main__":
    main()
