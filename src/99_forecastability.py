#!/usr/bin/env python3
"""Script 99: forecastability ceiling and signal-to-noise of observed change (Phase 2.5; masterplan §181-182).

Question before any model: how much of the next scan can be predicted at all?
  * Noise floor (assumption, stated): the change in core volume caused by moving the boundary one 2 mm voxel
    (dilation / erosion of the input mask). It is a proxy for segmentation uncertainty, not a measured one.
  * SNR_i = |ln((V_out + 1) / (V_in + 1))| / noise_floor_i, with noise_floor_i = ln((V_dil + 1) / (V_ero + 1)) / 2.
  * Forecastability ceiling: the Dice that persistence would get if the next mask were the input mask moved by
    the noise floor only (one-voxel erosion vs dilation vs the input), i.e. the best a "no-change plus noise" null can do.
  * Persistence Dice vs the oracle volume-matched Dice (script 98) shows how much a perfect volume forecast could add.

Input : data/manifests/forecast_pairs_mu.csv, output/baseline_ladder_pairs.csv, MU masks
Output: output/forecastability.json, output/forecastability_pairs.csv
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy import ndimage

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
MAN = PROJECT_ROOT / "data" / "manifests"
sys.path.insert(0, str(PROJECT_ROOT))

import run_improved_aniso as ria  # noqa: E402
from src import forecast_metrics as fm  # noqa: E402
from src import patient_stats as ps  # noqa: E402
from src.run_manifest import write_run_manifest  # noqa: E402

MU_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"


def core_mask(pid, tp):
    p = MU_DIR / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_tumorMask.nii.gz"
    return ria.to_2mm(np.isin(np.asarray(nib.load(str(p)).dataobj), (1, 3))) >= 0.5


def main():
    t0 = time.time()
    pairs = pd.read_csv(MAN / "forecast_pairs_mu.csv")
    pairs = pairs[pairs["in_primary"]]
    rows = []
    for i, r in enumerate(pairs.itertuples(index=False)):
        S, T = core_mask(r.patient_id, r.tp_in), core_mask(r.patient_id, r.tp_out)
        if not S.any():
            continue
        dil, ero = ndimage.binary_dilation(S), ndimage.binary_erosion(S)
        vS, vT, vD, vE = int(S.sum()), int(T.sum()), int(dil.sum()), int(ero.sum())
        floor = np.log((vD + 1) / (vE + 1)) / 2.0
        rows.append({"patient_id": r.patient_id, "tp_in": r.tp_in, "tp_out": r.tp_out, "dt_days": r.dt_days,
                     "v_in": vS, "v_out": vT, "noise_floor_ln": float(floor),
                     "abs_ln_change": float(abs(np.log((vT + 1) / (vS + 1)))),
                     "snr": float(abs(np.log((vT + 1) / (vS + 1))) / floor),
                     "dice_in_vs_dilated": fm.dice(S, dil), "dice_in_vs_eroded": fm.dice(S, ero),
                     "dice_persistence": fm.dice(S, T)})
        if i % 80 == 0:
            print(f"[{time.time() - t0:.0f}s] {i}/{len(pairs)}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "forecastability_pairs.csv", index=False)
    df["null_ceiling_dice"] = (df["dice_in_vs_dilated"] + df["dice_in_vs_eroded"]) / 2.0

    lad = pd.read_csv(OUT / "baseline_ladder_pairs.csv")
    orc = lad[lad["method"] == "oracle_volume_matched"][["patient_id", "tp_in", "tp_out", "dice"]].rename(columns={"dice": "dice_oracle_volume"})
    df = df.merge(orc, on=["patient_id", "tp_in", "tp_out"], how="left")

    pm = lambda col: ps.patient_means(df, col)  # noqa: E731
    res = {"script": "99_forecastability", "n_pairs": int(len(df)), "n_patients": int(df["patient_id"].nunique()),
           "assumption": "noise floor = half the ln-volume spread between one-voxel (2 mm) erosion and dilation of the input mask",
           "noise_floor_ln_median": float(df["noise_floor_ln"].median()),
           "abs_ln_change_median": float(df["abs_ln_change"].median()),
           "share_pairs_snr_lt_1": float((df["snr"] < 1).mean()),
           "share_pairs_snr_lt_2": float((df["snr"] < 2).mean()),
           "snr_quartiles": [float(x) for x in df["snr"].quantile([0.25, 0.5, 0.75])],
           "patient_mean_dice": {"persistence": float(pm("dice_persistence").mean()),
                                 "null_ceiling_noise_only": float(pm("null_ceiling_dice").mean()),
                                 "oracle_volume_matched": float(pm("dice_oracle_volume").mean())},
           "by_snr_bin": {}}
    df["snr_bin"] = pd.cut(df["snr"], [0, 1, 2, 4, 1e9], labels=["<1", "1-2", "2-4", ">4"])
    for b, g in df.groupby("snr_bin", observed=True):
        res["by_snr_bin"][str(b)] = {"n_pairs": int(len(g)), "n_patients": int(g["patient_id"].nunique()),
                                     "persistence_dice_mean": float(g["dice_persistence"].mean()),
                                     "oracle_volume_dice_mean": float(g["dice_oracle_volume"].mean())}
    (OUT / "forecastability.json").write_text(json.dumps(res, indent=1, default=float))
    write_run_manifest("forecastability_99", OUT / "forecastability.manifest.json", script="src/99_forecastability.py",
                       seed=ps.DEFAULT_SEED, config={"noise": "1-voxel erosion/dilation"},
                       inputs=[MAN / "forecast_pairs_mu.csv", OUT / "baseline_ladder_pairs.csv"], dataset="MU-Glioma-Post",
                       patient_split="n/a", primary_endpoint="n/a (descriptive)")
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
