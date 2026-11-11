"""Frozen .tfm direction check (paper/trackB_tfm_direction_check_protocol_2026-11-10.md). Brain masks only.

Reads, per scan: native ct1_brain_mask, atlas brain_mask (both non-tumor masks), ct1.tfm text. No CT1 intensities,
no tumor labels. Usage: python -m src.fk_lumiere_tfm_direction_check
"""
from __future__ import annotations

import gzip
import json
import zipfile

import nibabel as nib
import numpy as np
from scipy.ndimage import affine_transform

from src.fk_lumiere_transform_audit import INVENTORY, OUT_DIR, TFM, ZIP_PATH, _visit_dirs, parse_itk_tfm, tfm_to_ras_4x4

NATIVE_MASK = "DeepBraTumIA-segmentation/native/skull_strip/ct1_brain_mask.nii.gz"
ATLAS_MASK = "DeepBraTumIA-segmentation/atlas/skull_strip/brain_mask.nii.gz"
ALLOWED = (NATIVE_MASK, ATLAS_MASK, TFM)
MIN_MEDIAN_DICE, MIN_WINS, MIN_MEDIAN_GAP = 0.80, 80, 0.10


def resample_mask(native_mask, native_affine, atlas_shape, atlas_affine, atlas_to_native_ras):
    m = np.linalg.inv(native_affine) @ atlas_to_native_ras @ atlas_affine
    return affine_transform(np.asarray(native_mask, np.uint8), m[:3, :3], offset=m[:3, 3], output_shape=tuple(atlas_shape),
                            order=0, mode="constant", cval=0, prefilter=False).astype(bool)


def dice(a, b):
    s = a.sum() + b.sum()
    return 1.0 if s == 0 else float(2 * (a & b).sum() / s)


def compare(native_mask, native_affine, atlas_mask, atlas_affine, T):
    out = {}
    for name, M in (("F_T_atlas_to_native", T), ("I_inverse_T", np.linalg.inv(T))):
        out[name] = dice(resample_mask(native_mask, native_affine, atlas_mask.shape, atlas_affine, M), atlas_mask)
    return out


def decide(rows, n_total):
    f = np.array([r["F_T_atlas_to_native"] for r in rows])
    i = np.array([r["I_inverse_T"] for r in rows])
    res = {"median_dice": {"F": float(np.median(f)), "I": float(np.median(i))},
           "wins": {"F": int((f > i).sum()), "I": int((i > f).sum()), "ties": int((f == i).sum())},
           "median_gap": {"F_minus_I": float(np.median(f - i)), "I_minus_F": float(np.median(i - f))},
           "n_scans": len(rows)}
    resolved = None
    for x, a, b in (("F", f, i), ("I", i, f)):
        if (len(rows) == n_total and np.median(a) >= MIN_MEDIAN_DICE and (a > b).sum() >= MIN_WINS
                and np.median(a - b) >= MIN_MEDIAN_GAP):
            resolved = x
    res["decision"] = resolved or "unresolved"
    return res


def _load(zf, name):
    with zf.open(name) as fh:
        img = nib.Nifti1Image.from_bytes(gzip.decompress(fh.read()))
    return np.asarray(img.dataobj) > 0, np.asarray(img.affine, float)


def run():
    inv = json.loads(INVENTORY.read_text())
    rows, opened = [], []
    with zipfile.ZipFile(ZIP_PATH) as zf:
        names = zf.namelist()
        for pid, rec in inv["patients"].items():
            dirs = _visit_dirs(names, pid)
            for k, v in enumerate(rec["visits"]):
                d = dirs[(v["week"], v["order"])]
                nm, na = _load(zf, d + NATIVE_MASK)
                am, aa = _load(zf, d + ATLAS_MASK)
                tf = parse_itk_tfm(zf.read(d + TFM).decode())
                opened += [d + NATIVE_MASK, d + ATLAS_MASK, d + TFM]
                rows.append({"patient": pid, "visit": k + 1, **compare(nm, na, am, aa, tfm_to_ras_4x4(tf["leaves"]))})
    assert all(any(o.endswith(a) for a in ALLOWED) for o in opened), "opened a file outside the allowed set"
    return {"protocol": "paper/trackB_tfm_direction_check_protocol_2026-11-10.md",
            "read": "native ct1_brain_mask, atlas brain_mask, ct1.tfm only; no intensities, no tumor labels",
            "n_files_opened": len(opened), "rule": {"min_median_dice": MIN_MEDIAN_DICE, "min_wins": MIN_WINS,
                                                    "min_median_gap": MIN_MEDIAN_GAP},
            **decide(rows, 84), "rows": rows}


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res = run()
    (OUT_DIR / "tfm_direction_check.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "rows"}, indent=1))
