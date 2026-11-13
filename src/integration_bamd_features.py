"""BAMD step 2a: fixed MRI features from BraTS 2021 images and expert masks.

Protocol: paper/integration_protocol_bamd_2026-11-12.md (59 features). Reads images and
masks only; no survival or expression data. Writes output/integration_bamd/features.csv.

  python -I src/integration_bamd_features.py [--workers 6]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BRATS = PROJECT_ROOT / "data" / "brats"
MAP_PATH = PROJECT_ROOT / "data" / "external" / "brats2021_map" / "BraTS2021_MappingToTCIA.xlsx"
OUT = PROJECT_ROOT / "output" / "integration_bamd"
COLLECTIONS = ("TCGA-GBM", "UPENN-GBM")
REGIONS = {"ET": (4,), "NCR": (1,), "ED": (2,)}
IMAGES = ("t1", "t1ce", "t2", "flair")
VOXEL_CM3 = 0.001  # BraTS 2021 is resampled to 1 mm isotropic


def feature_names() -> list[str]:
    names = [f"logvol_{r}" for r in ("ET", "NCR", "ED", "WT")]
    names += ["ratio_ncr_core", "ratio_et_wt", "ratio_ed_wt"]
    names += ["sv_wt", "sph_wt", "sv_tc", "sph_tc"]
    names += [f"{s}_{r}_{im}" for r in REGIONS for im in IMAGES for s in ("mean", "sd", "p10", "p90")]
    return names


def surface_voxels(mask: np.ndarray) -> int:
    """Count exposed voxel faces (6-neighbour), a simple surface-area proxy in mm^2."""
    m = mask.astype(np.int8)
    faces = 0
    for ax in range(3):
        d = np.diff(np.pad(m, [(1, 1) if a == ax else (0, 0) for a in range(3)]), axis=ax)
        faces += int(np.abs(d).sum())
    return faces


def shape(mask: np.ndarray) -> tuple[float, float]:
    v = float(mask.sum())
    if v == 0:
        return 0.0, 0.0
    s = float(surface_voxels(mask))
    sph = (np.pi ** (1 / 3)) * (6 * v) ** (2 / 3) / s
    return s / v, sph


def case_features(case_id: str, root: Path = BRATS) -> dict:
    import nibabel as nib

    d = root / case_id
    seg = np.asarray(nib.load(str(d / f"{case_id}_seg.nii.gz")).dataobj).astype(np.int16)
    imgs = {im: np.asarray(nib.load(str(d / f"{case_id}_{im}.nii.gz")).dataobj, dtype=np.float32) for im in IMAGES}
    feats = {"brats_id": case_id}
    masks = {r: np.isin(seg, labs) for r, labs in REGIONS.items()}
    wt = seg > 0
    tc = masks["ET"] | masks["NCR"]
    vols = {r: float(m.sum()) * VOXEL_CM3 for r, m in masks.items()}
    vols["WT"] = float(wt.sum()) * VOXEL_CM3
    for r in ("ET", "NCR", "ED", "WT"):
        feats[f"logvol_{r}"] = float(np.log1p(vols[r]))
    core = vols["NCR"] + vols["ET"]
    feats["ratio_ncr_core"] = vols["NCR"] / core if core > 0 else 0.0
    feats["ratio_et_wt"] = vols["ET"] / vols["WT"] if vols["WT"] > 0 else 0.0
    feats["ratio_ed_wt"] = vols["ED"] / vols["WT"] if vols["WT"] > 0 else 0.0
    feats["sv_wt"], feats["sph_wt"] = shape(wt)
    feats["sv_tc"], feats["sph_tc"] = shape(tc)
    for im, x in imgs.items():
        brain = x != 0
        mu, sd = float(x[brain].mean()), float(x[brain].std())
        z = (x - mu) / sd if sd > 0 else x * 0
        for r, m in masks.items():
            vals = z[m]
            if vals.size == 0:
                for s in ("mean", "sd", "p10", "p90"):
                    feats[f"{s}_{r}_{im}"] = np.nan  # imputed with training median later
                continue
            feats[f"mean_{r}_{im}"] = float(vals.mean())
            feats[f"sd_{r}_{im}"] = float(vals.std())
            feats[f"p10_{r}_{im}"] = float(np.percentile(vals, 10))
            feats[f"p90_{r}_{im}"] = float(np.percentile(vals, 90))
    return feats


def case_table() -> pd.DataFrame:
    m = pd.read_excel(MAP_PATH)
    m = m.rename(columns={m.columns[0]: "collection", "PatientID on TCIA Radiology Portal": "patient_id",
                          "BraTS2021 ID": "brats_id"})
    m = m[m["collection"].isin(COLLECTIONS)][["collection", "patient_id", "brats_id"]].copy()
    m["patient_id"] = m["patient_id"].astype(str)
    m["on_disk"] = [(BRATS / b / f"{b}_seg.nii.gz").exists() for b in m["brats_id"]]
    return m[m["on_disk"]].drop(columns="on_disk").reset_index(drop=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=None, help="testing only")
    args = ap.parse_args(argv)
    cases = case_table()
    if args.limit:
        cases = cases.head(args.limit)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(case_features, cases["brats_id"], chunksize=4))
    feats = pd.DataFrame(rows)
    assert set(feats.columns[1:]) == set(feature_names())
    feats = feats[["brats_id", *feature_names()]]
    out = cases.merge(feats, on="brats_id", validate="one_to_one")
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "features.csv"
    out.to_csv(path, index=False)
    meta = {"n_cases": len(out), "by_collection": out["collection"].value_counts().to_dict(),
            "n_features": len(feature_names()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / "features_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
