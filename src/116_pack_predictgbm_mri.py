#!/usr/bin/env python3
"""Script 116: pack pre-op MRI (T1c, FLAIR; already brain-extracted and scaled 0-1 by PREDICT-GBM) at 2 mm for the
development (TUM) patients, for the script 115 MRI arms on Kaggle (GRAND_PLAN 12 B2, user option 1, 2026-10-10).
Only pre-op images. t1c_warped_longitudinal (the follow-up scan) is never read: it would leak the recurrence.
Test patients are not packed. Private Kaggle dataset vihansanthosh/predictgbm-mri2mm-private.
Output: kaggle_run/pack_mri/pg_mri2mm.npz, kaggle_run/pack_mri_manifest.json
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path as _Path

import nibabel as nib
import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
PACK = PROJECT_ROOT / "kaggle_run" / "pack_mri"
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402

MODS = ["t1c_bet_normalized.nii.gz", "flair_bet_normalized.nii.gz"]


def to2(a):
    return np.asarray(a, np.float32)[:240, :240, :154].reshape(120, 2, 120, 2, 77, 2).mean((1, 3, 5))


def main():
    PACK.mkdir(parents=True, exist_ok=True)
    ids, test = pg.dev_ids(), pg.test_ids()
    assert not set(ids) & test
    arr = {f"{pid}|{m}": to2(nib.load(str(pg.DATA / pid / m)).dataobj).astype(np.float16) for pid in ids for m in MODS}
    out = PACK / "pg_mri2mm.npz"
    np.savez_compressed(out, **arr)
    (PACK / "dataset-metadata.json").write_text(json.dumps({"title": "predictgbm-mri2mm-private",
                                                            "id": "vihansanthosh/predictgbm-mri2mm-private",
                                                            "licenses": [{"name": "other"}]}, indent=1))
    man = {"script": "116_pack_predictgbm_mri", "n_patients": len(ids), "modalities": MODS, "test_patients_packed": 0,
           "followup_scan_packed": False, "pack_mb": round(out.stat().st_size / 1e6, 1),
           "pack_sha256": hashlib.sha256(out.read_bytes()).hexdigest(), "visibility": "private"}
    (PROJECT_ROOT / "kaggle_run" / "pack_mri_manifest.json").write_text(json.dumps(man, indent=1))
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
