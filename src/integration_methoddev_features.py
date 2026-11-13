"""Method development: richer MRI features (location, texture, habitats, peritumoral ring).

Protocol: paper/integration_protocol_methoddev_2026-11-12.md. Images and masks only.
--cohort selects which collection to extract. Confirmation cohorts (CPTAC-GBM, IvyGAP)
must not be extracted until the lock note exists (checked here).

  python -I src/integration_methoddev_features.py --cohort TCGA-GBM [--workers 6]
"""
from __future__ import annotations

import argparse
import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import ndimage

PROJECT_ROOT = Path(__file__).resolve().parent.parent
BRATS = PROJECT_ROOT / "data" / "brats"
MAP_PATH = PROJECT_ROOT / "data" / "external" / "brats2021_map" / "BraTS2021_MappingToTCIA.xlsx"
OUT = PROJECT_ROOT / "output" / "integration_methoddev"
LOCK_NOTE = PROJECT_ROOT / "paper" / "integration_methoddev_lock_2026-11-12.md"
CONFIRM = ("CPTAC-GBM", "IvyGAP")
IMAGES = ("t1", "t1ce", "t2", "flair")
REGIONS = {"ET": (4,), "NCR": (1,), "ED": (2,)}
LEVELS = 32


def normalize(x, brain, tumor):
    """Scale by normal-appearing brain (brain minus tumor): median and IQR."""
    ref = x[brain & ~tumor]
    med = np.median(ref)
    iqr = np.subtract(*np.percentile(ref, [75, 25]))
    return (x - med) / (iqr if iqr > 0 else 1.0)


def glcm_features(vals_img, mask, lo=-3.0, hi=6.0):
    """3D GLCM over 13 directions, distance 1, symmetric, averaged. Returns 5 features."""
    if mask.sum() < 27:
        return dict.fromkeys(("contrast", "homogeneity", "energy", "correlation", "entropy"), np.nan)
    q = np.clip(((vals_img - lo) / (hi - lo) * LEVELS).astype(int), 0, LEVELS - 1)
    dirs = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1),
            (0, 1, 1), (0, 1, -1), (1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1)]
    zz, yy, xx = np.nonzero(mask)
    P = np.zeros((LEVELS, LEVELS))
    sh = mask.shape
    for d in dirs:
        z2, y2, x2 = zz + d[0], yy + d[1], xx + d[2]
        ok = (z2 >= 0) & (z2 < sh[0]) & (y2 >= 0) & (y2 < sh[1]) & (x2 >= 0) & (x2 < sh[2])
        ok[ok] = mask[z2[ok], y2[ok], x2[ok]]
        a, b = q[zz[ok], yy[ok], xx[ok]], q[z2[ok], y2[ok], x2[ok]]
        np.add.at(P, (a, b), 1)
        np.add.at(P, (b, a), 1)
    P /= P.sum()
    i, j = np.indices(P.shape)
    mu_i, mu_j = (i * P).sum(), (j * P).sum()
    sd_i, sd_j = np.sqrt(((i - mu_i) ** 2 * P).sum()), np.sqrt(((j - mu_j) ** 2 * P).sum())
    nz = P[P > 0]
    return {"contrast": float(((i - j) ** 2 * P).sum()), "homogeneity": float((P / (1 + np.abs(i - j))).sum()),
            "energy": float((P ** 2).sum()),
            "correlation": float(((i - mu_i) * (j - mu_j) * P).sum() / (sd_i * sd_j)) if sd_i * sd_j > 0 else 0.0,
            "entropy": float(-(nz * np.log2(nz)).sum())}


def case_features(case_id: str, root: Path = BRATS) -> dict:
    import nibabel as nib
    d = root / case_id
    seg = np.asarray(nib.load(str(d / f"{case_id}_seg.nii.gz")).dataobj).astype(np.int16)
    raw = {im: np.asarray(nib.load(str(d / f"{case_id}_{im}.nii.gz")).dataobj, dtype=np.float32) for im in IMAGES}
    brain = raw["t1"] != 0
    wt = seg > 0
    tc = np.isin(seg, (1, 4))
    img = {im: normalize(x, brain, wt) for im, x in raw.items()}
    f = {"brats_id": case_id}
    # location (BraTS images are co-registered to the SRI24 atlas, 240 x 240 x 155, 1 mm)
    bc = np.array(ndimage.center_of_mass(brain))
    for name, m in (("wt", wt), ("tc", tc)):
        if m.sum() == 0:
            for k in ("x", "y", "z", "absx", "dist_center"):
                f[f"loc_{name}_{k}"] = np.nan
            continue
        c = np.array(ndimage.center_of_mass(m)) - bc
        f[f"loc_{name}_x"], f[f"loc_{name}_y"], f[f"loc_{name}_z"] = map(float, c)
        f[f"loc_{name}_absx"] = float(abs(c[0]))
        f[f"loc_{name}_dist_center"] = float(np.linalg.norm(c))
    # depth: mean distance of tumor core voxels to the brain surface
    edt = ndimage.distance_transform_edt(brain)
    f["loc_tc_depth"] = float(edt[tc].mean()) if tc.any() else np.nan
    f["loc_wt_min_depth"] = float(edt[wt].min()) if wt.any() else np.nan
    # texture: T1ce in ET, FLAIR in ED, T2 in NCR, T1ce in NCR
    for im, r in (("t1ce", "ET"), ("flair", "ED"), ("t2", "NCR"), ("t1ce", "NCR"), ("flair", "ET")):
        for k, v in glcm_features(img[im], np.isin(seg, REGIONS[r])).items():
            f[f"glcm_{k}_{r}_{im}"] = v
    # habitats: T1ce high/low x FLAIR high/low inside whole tumor (fixed thresholds on normalized scale)
    if wt.any():
        hi_c, hi_f = img["t1ce"][wt] > 2.0, img["flair"][wt] > 2.0
        for a in (0, 1):
            for b in (0, 1):
                f[f"hab_c{a}_f{b}"] = float(((hi_c == a) & (hi_f == b)).mean())
    else:
        for a in (0, 1):
            for b in (0, 1):
                f[f"hab_c{a}_f{b}"] = np.nan
    # peritumoral ring: 1-5 mm outside the whole tumor, inside brain
    ring = ndimage.binary_dilation(wt, iterations=5) & ~wt & brain
    for im in IMAGES:
        v = img[im][ring]
        f[f"ring_mean_{im}"] = float(v.mean()) if v.size else np.nan
        f[f"ring_sd_{im}"] = float(v.std()) if v.size else np.nan
    # normalized intensity summaries per region (NAWM-referenced, unlike the Gate 1 z-score)
    for r, labs in REGIONS.items():
        m = np.isin(seg, labs)
        for im in IMAGES:
            v = img[im][m]
            f[f"nmean_{r}_{im}"] = float(v.mean()) if v.size else np.nan
            f[f"nsd_{r}_{im}"] = float(v.std()) if v.size else np.nan
    return f


def cases(cohort: str) -> pd.DataFrame:
    m = pd.read_excel(MAP_PATH)
    m = m.rename(columns={m.columns[0]: "collection", "PatientID on TCIA Radiology Portal": "patient_id",
                          "BraTS2021 ID": "brats_id"})
    m = m[m["collection"] == cohort][["collection", "patient_id", "brats_id"]]
    m["patient_id"] = m["patient_id"].astype(str)
    return m[[(BRATS / b / f"{b}_seg.nii.gz").exists() for b in m["brats_id"]]].reset_index(drop=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cohort", required=True)
    ap.add_argument("--workers", type=int, default=6)
    args = ap.parse_args(argv)
    if args.cohort in CONFIRM and not LOCK_NOTE.exists():
        raise SystemExit("confirmation cohort refused: lock note does not exist yet")
    c = cases(args.cohort)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(case_features, c["brats_id"], chunksize=2))
    out = c.merge(pd.DataFrame(rows), on="brats_id", validate="one_to_one")
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / f"features_{args.cohort}.csv"
    out.to_csv(path, index=False)
    meta = {"cohort": args.cohort, "n": len(out), "n_features": out.shape[1] - 3,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (OUT / f"features_{args.cohort}_meta.json").write_text(json.dumps(meta, indent=2))
    print(meta)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
