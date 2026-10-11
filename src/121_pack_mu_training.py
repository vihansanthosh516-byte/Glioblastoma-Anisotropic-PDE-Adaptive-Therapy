#!/usr/bin/env python3
"""Script 121: pack MU-Glioma-Post pairs as extra U-Net TRAINING samples (GRAND_PLAN 15 Z2, declared before the run).

Per MU primary pair (data/manifests/forecast_pairs_mu.csv, in_primary): input-scan core (labels 1, 3), edema (2), brain,
T1c and FLAIR (0.5-99.5 percentile inside the brain -> 0-1), target all recurrence (labels 1-3 at the later scan, inside
the input brain) and enhancing (label 3); all at 2 mm (block mean of the 1 mm maps, 240 x 240 x 154 crop), flipped in
x and y into PREDICT-GBM array order. k2 = 15 mm standard-plan voxel count of the input core / 8.
Tissue proxy for MU: mean of the TUM development wm / gm / csf maps at 2 mm (dev only; no test patient is read).
Private Kaggle dataset vihansanthosh/mu-train-2mm-private. Output: kaggle_run/pack_mu/mu_train_2mm.npz + manifest
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path as _Path

import nibabel as nib
import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
PACK = PROJECT_ROOT / "kaggle_run" / "pack_mu"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402
import run_improved_aniso as ria  # noqa: E402


def to2(a):
    return np.asarray(a, np.float32)[:240, :240, :154].reshape(120, 2, 120, 2, 77, 2).mean((1, 3, 5))


def flip(a):
    return a[::-1, ::-1].copy()


def load(pid, tp, kind):
    return np.asarray(nib.load(str(ria.MU_DIR / pid / f"Timepoint_{tp}" / f"{pid}_Timepoint_{tp}_{kind}.nii.gz")).dataobj)


def norm(img, brain):
    v = img[brain]
    lo, hi = np.percentile(v, [0.5, 99.5]) if v.size else (0, 1)
    return np.clip((img - lo) / max(hi - lo, 1e-6), 0, 1) * brain


def main():
    PACK.mkdir(parents=True, exist_ok=True)
    arr = {}
    tis = {t: np.zeros((120, 120, 77), np.float64) for t in ("wm", "gm", "csf")}
    dev = pg.dev_ids()
    for pid in dev:
        for t in tis:
            tis[t] += to2(pg._load(pid, f"{t}_pbmap.nii.gz"))
    for t in tis:
        arr[f"atlas|{t}"] = (tis[t] / len(dev)).astype(np.float16)
    pairs = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "forecast_pairs_mu.csv")
    pairs = pairs[pairs["in_primary"]]
    keys = []
    for _, r in pairs.iterrows():
        seg_in = np.rint(load(r.patient_id, r.tp_in, "tumorMask"))[:240, :240, :154]
        seg_out = np.rint(load(r.patient_id, r.tp_out, "tumorMask"))[:240, :240, :154]
        brain = (load(r.patient_id, r.tp_in, "brain_t1n") > 0)[:240, :240, :154]
        core = np.isin(seg_in, (1, 3))
        if not core.any():
            continue
        std = (distance_transform_edt(~core) <= 15) & brain
        key = f"{r.patient_id}|{r.tp_in}|{r.tp_out}"
        ch = {"core": core, "edema": seg_in == 2, "brain": brain,
              "t1c": norm(load(r.patient_id, r.tp_in, "brain_t1c")[:240, :240, :154].astype(np.float32), brain),
              "flair": norm(load(r.patient_id, r.tp_in, "brain_t2f")[:240, :240, :154].astype(np.float32), brain),
              "y": np.isin(seg_out, (1, 2, 3)) & brain, "ye": (seg_out == 3) & brain}
        for n, a in ch.items():
            arr[f"{key}|{n}"] = flip(to2(a)).astype(np.float16)
        arr[f"{key}|k2"] = np.array(int(round(std.sum() / 8.0)))
        keys.append(key)
    out = PACK / "mu_train_2mm.npz"
    np.savez_compressed(out, keys=np.array(keys), **arr)
    (PACK / "dataset-metadata.json").write_text(json.dumps({"title": "mu-train-2mm-private", "id": "vihansanthosh/mu-train-2mm-private",
                                                            "licenses": [{"name": "other"}]}, indent=1))
    man = {"script": "121_pack_mu_training", "n_pairs": len(keys), "n_patients": len({k.split("|")[0] for k in keys}),
           "tissue_proxy_from_dev_patients": len(dev), "test_patients_read": 0, "pack_mb": round(out.stat().st_size / 1e6, 1),
           "pack_sha256": hashlib.sha256(out.read_bytes()).hexdigest(), "visibility": "private"}
    (PROJECT_ROOT / "kaggle_run" / "pack_mu_manifest.json").write_text(json.dumps(man, indent=1))
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
