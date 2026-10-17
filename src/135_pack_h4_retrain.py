#!/usr/bin/env python3
"""Script 135: H-4 step 2 retrain pack (Amendment 13, A13.2).

Assembles the training samples of the three post-hoc sets for the pooled retrain of arm hybrid_mri_rank_mu:
  LUMIERE (61) + RHUH (40): sample = s115.build(pid, hybrid=True, mri=True) with the MRI read from the
    PREDICT-GBM carrier files t1c_bet_normalized.nii.gz / flair_bet_normalized.nii.gz, block-averaged to 2 mm
    exactly as script 116; recurrence time1 comes from the carrier recurrence_preop.nii.gz through the pg
    loader (allowed: configs/h2_frozen.yaml + tag h2-frozen exist).
  UCSF (59): sample = script 131 build of the script 124 pre-op arrays (same 10 channels, dev-mean tissue
    atlas, A11.4); targets = ALPTDG time1 seg (labels 1-3 all / 3 enhancing inside the pre-op brain), read
    here for the first time as TRAINING targets (the set is post-hoc dev since A12.1).
No time2 target is packed (endpoints are time1). TUM and MU are rebuilt at train time and are not packed here.
Output: kaggle_run/pack_h4/h4_retrain_pack.npz (+ dataset-metadata.json, private), manifest
output/predictgbm/h4/retrain_pack_manifest.json
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import os
import sys
import time
import zipfile
from pathlib import Path as _Path

import nibabel as nib
import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
PACK_DIR = PROJECT_ROOT / "kaggle_run" / "pack_h4"
MANIFEST = PROJECT_ROOT / "output" / "predictgbm" / "h4" / "retrain_pack_manifest.json"
H3 = _Path(os.environ.get("UCSF_H3_DIR", r"C:\Users\vihan\data\ucsf_h3"))
ZIP = _Path(os.environ.get("UCSF_ALPTDG_ZIP", r"C:\Users\vihan\Downloads\UCSF_POSTOP_GLIOMA_DATASET_FINAL_v1.0.zip"))
MU_PACK = _Path(os.environ.get("GBM_MU_PACK", str(PROJECT_ROOT / "kaggle_run" / "pack_mu" / "mu_train_2mm.npz")))
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import predictgbm_io as pg  # noqa: E402

_s = importlib.util.spec_from_file_location("s115", PROJECT_ROOT / "src" / "115_predictgbm_unet_dev.py")
s115 = importlib.util.module_from_spec(_s)
_s.loader.exec_module(s115)
s112 = s115.s112

_s131 = importlib.util.spec_from_file_location("s131", PROJECT_ROOT / "src" / "131_pack_ucsf_h3_model_inputs.py")
s131 = importlib.util.module_from_spec(_s131)
_s131.loader.exec_module(s131)

MODS = ["t1c_bet_normalized.nii.gz", "flair_bet_normalized.nii.gz"]


def to2(a):
    return np.asarray(a, np.float32)[:240, :240, :154].reshape(120, 2, 120, 2, 77, 2).mean((1, 3, 5))


def alptdg_label(z, sid, tp="time1"):
    a = np.asarray(nib.Nifti1Image.from_bytes(gzip.decompress(z.read(f"{sid}/{sid}_{tp}_seg.nii.gz"))).dataobj)
    return np.ascontiguousarray(np.rint(a).astype(np.uint8)[::-1, ::-1])


def build_carrier_mri(pids):
    """Extended 2 mm MRI pack (t1c, flair) for the LUMIERE/RHUH patients, same normalisation as script 116."""
    arr = {f"{pid}|{m}": to2(nib.load(str(pg.DATA / pid / m)).dataobj).astype(np.float16) for pid in pids for m in MODS}
    out = PACK_DIR / "pg_mri_ext_2mm.npz"
    np.savez_compressed(out, **arr)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="smoke test: max patients per set (0 = all)")
    args = ap.parse_args()
    PACK_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    lum = sorted(p["patient_id"] for p in json.loads((pg.LISTS / "lumiere.json").read_text())["patients"])
    rhu = sorted(p["patient_id"] for p in json.loads((pg.LISTS / "rhuh.json").read_text())["patients"])
    ucsf = sorted(f.stem for f in (H3 / "inputs").glob("*.npz"))
    if args.limit:
        lum, rhu, ucsf = lum[: args.limit], rhu[: args.limit], ucsf[: args.limit]

    ext = build_carrier_mri(lum + rhu)
    os.environ["GBM_PG_MRI"] = str(ext)
    print(f"carrier MRI pack {ext.name}: {2 * (len(lum) + len(rhu))} arrays", flush=True)

    arr, stats = {}, {}
    for setname, ids in (("lumiere", lum), ("rhuh", rhu)):
        n_empty_y = n_ye = 0
        for pid in ids:
            d = s115.build(pid, hybrid=True, mri=True)
            x = d["x"]
            assert x.shape == (10, *s115.SHAPE), (pid, x.shape)
            arr[f"{pid}|x"] = x
            arr[f"{pid}|y"] = np.asarray(d["y"], np.float32)
            arr[f"{pid}|ye"] = np.asarray(d["ye"], np.float32)
            arr[f"{pid}|k2"] = np.array(d["k2"], np.int64)
            arr[f"{pid}|brain2"] = np.asarray(d["brain2"], np.float32)
            n_empty_y += int(d["y"].sum() == 0)
            n_ye += int(d["ye"].sum() > 0)
            print(f"  {setname} {pid} y={int(d['y'].sum())} ye={int(d['ye'].sum())} k2={d['k2']}", flush=True)
        stats[setname] = {"n": len(ids), "n_empty_recurrence": n_empty_y, "n_with_enhancing": n_ye}

    mu = np.load(MU_PACK)
    s131.ATLAS.update({t: mu[f"atlas|{t}"].astype(np.float32) for t in ("wm", "gm", "csf")})
    z = zipfile.ZipFile(ZIP)
    names = set(z.namelist())
    n_empty_y = n_ye = 0
    for sid in ucsf:
        assert f"{sid}/{sid}_time1_seg.nii.gz" in names, f"no ALPTDG time1 label for {sid}"
        inp = np.load(H3 / "inputs" / f"{sid}.npz")
        seg, brain = np.rint(inp["seg"]).astype(np.uint8), inp["brain"] > 0
        std = pg.standard_plan(seg, brain)
        x = s131.build(inp)
        lab = alptdg_label(z, sid)
        y = (s112.to2((np.isin(lab, (1, 2, 3))) & brain) >= 0.5).astype(np.float32)
        ye = (s112.to2((lab == 3) & brain) >= 0.5).astype(np.float32)
        assert x.shape == (10, *s115.SHAPE), (sid, x.shape)
        arr[f"{sid}|x"] = x
        arr[f"{sid}|y"] = y
        arr[f"{sid}|ye"] = ye
        arr[f"{sid}|k2"] = np.array(int(round(std.sum() / 8.0)), np.int64)
        arr[f"{sid}|brain2"] = (s112.to2(brain) >= 0.5).astype(np.float32)
        n_empty_y += int(y.sum() == 0)
        n_ye += int(ye.sum() > 0)
        print(f"  ucsf {sid} y={int(y.sum())} ye={int(ye.sum())}", flush=True)
    stats["ucsf"] = {"n": len(ucsf), "n_empty_recurrence": n_empty_y, "n_with_enhancing": n_ye}

    out = PACK_DIR / "h4_retrain_pack.npz"
    np.savez_compressed(out, **arr)
    (PACK_DIR / "dataset-metadata.json").write_text(json.dumps(
        {"title": "h4-retrain-pack-private", "id": "vihansanthosh/h4-retrain-pack-private",
         "licenses": [{"name": "other"}]}, indent=1))
    man = {"script": "135_pack_h4_retrain", "amendment": "A13.2", "n_patients": sum(s["n"] for s in stats.values()),
           "stats": stats, "channels": 10, "threads": "time1 only, no time2",
           "pack_file": str(out.relative_to(PROJECT_ROOT)), "pack_mb": round(out.stat().st_size / 1e6, 1),
           "pack_sha256": hashlib.sha256(out.read_bytes()).hexdigest(), "minutes": round((time.time() - t0) / 60, 1),
           "visibility": "private"}
    MANIFEST.write_text(json.dumps(man, indent=1))
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
