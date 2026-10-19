#!/usr/bin/env python3
"""Script 131: H-3 test inputs for the final model (C1 = script 115 arm hybrid_mri_rank_mu), from script 124 outputs.

Per UCSF patient (script 124 npz, PREDICT-GBM array order, 1 mm): the same 10 channels (tissue = dev-mean atlas, A11.4) as s115.build(hybrid=True,
mri=True), built from pre-op data only: core (labels 1, 3), edema (2), wm, gm, csf, brain, distance / 50 (cap 3),
M1 density (lambda 2 mm), T1c, FLAIR. T1c/FLAIR are scaled as the PREDICT-GBM *_bet_normalized maps (BrainLes
percentile normalizer 0.1 / 99.9 of the skull-stripped volume, clipped to [0, 1]; script 125) and then block-averaged
to 2 mm as script 116. No ALPTDG label is read (script 124 never opened one).
Output: kaggle_run/pack_ucsf/ucsf_h3_inputs_2mm.npz (+ dataset-metadata.json, private), manifest
output/predictgbm/h3/ucsf_inputs_2mm_manifest.json
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path as _Path

import numpy as np
from scipy.ndimage import distance_transform_edt

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
PACK = PROJECT_ROOT / "kaggle_run" / "pack_ucsf"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402

_s = importlib.util.spec_from_file_location("s115", PROJECT_ROOT / "src" / "115_predictgbm_unet_dev.py")
s115 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(s115)
s112 = s115.s112


def bl_norm(a, brain):
    a = a * brain
    lo, hi = np.percentile(a, [0.1, 99.9])
    return np.clip((a - lo) / max(hi - lo, 1e-6), 0, 1) * brain


def build(z):
    seg, brain = np.rint(z["seg"]).astype(np.uint8), z["brain"] > 0
    c = {"seg": seg, "brain": brain, "wm": z["wm"].astype(np.float32), "gm": z["gm"].astype(np.float32)}
    core1 = pg.core(seg)
    std = pg.standard_plan(seg, brain)
    core2 = s112.to2(core1) >= 0.5
    dist = distance_transform_edt(~core1)[:240, :240, :154]
    dist2 = dist.reshape(120, 2, 120, 2, 77, 2).min((1, 3, 5)) / 50.0
    wm, gm, csf = ATLAS["wm"], ATLAS["gm"], ATLAS["csf"]   # A11.4: dev-mean tissue atlas, as the MU training pairs
    brain2 = s112.to2(brain)
    x = [core2, s112.to2(seg == 2), wm, gm, csf, brain2, np.minimum(dist2, 3.0)]
    x.append(s115.m1_density2(wm, gm, brain2 >= 0.5, core2, std.sum() / 8.0) if core2.any() else np.zeros(s115.SHAPE, np.float32))
    x += [s112.to2(bl_norm(z["t1c"].astype(np.float32), brain)), s112.to2(bl_norm(z["flair"].astype(np.float32), brain))]
    return np.stack(x).astype(np.float16)


ATLAS = {}


def main():
    mu = np.load(os.environ.get("GBM_MU_PACK", str(PROJECT_ROOT / "kaggle_run" / "pack_mu" / "mu_train_2mm.npz")))
    ATLAS.update({t: mu[f"atlas|{t}"].astype(np.float32) for t in ("wm", "gm", "csf")})
    PACK.mkdir(parents=True, exist_ok=True)
    arr = {}
    for f in sorted((H3 / "inputs").glob("*.npz")):
        arr[f"ucsf_{f.stem}"] = build(np.load(f))
        print(f.stem, arr[f"ucsf_{f.stem}"].shape, flush=True)
    out = PACK / "ucsf_h3_inputs_2mm.npz"
    np.savez_compressed(out, **arr)
    (PACK / "dataset-metadata.json").write_text(json.dumps({"title": "ucsf-h3-inputs-private", "id": "vihansanthosh/ucsf-h3-inputs-private",
                                                            "licenses": [{"name": "other"}]}, indent=1))
    man = {"script": "131_pack_ucsf_h3_model_inputs", "n_patients": len(arr), "channels": 10, "alptdg_labels_read": 0,
           "pack_mb": round(out.stat().st_size / 1e6, 1), "pack_sha256": hashlib.sha256(out.read_bytes()).hexdigest(),
           "visibility": "private"}
    (PROJECT_ROOT / "output" / "predictgbm" / "h3").mkdir(parents=True, exist_ok=True)
    (PROJECT_ROOT / "output" / "predictgbm" / "h3" / "ucsf_inputs_2mm_manifest.json").write_text(json.dumps(man, indent=1))
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
