#!/usr/bin/env python3
"""Script 113: private Kaggle pack of PREDICT-GBM inputs for the GPU runs of script 112 (MIT-licensed public data).

All 243 patients: tumor_seg (int8), t1c_bet_mask (bits), wm/gm/csf probability maps (float16, 1 mm).
Recurrence masks ONLY for development (TUM) patients; test recurrences are not packed, so the A6.6 lock also holds on
Kaggle. Test ids are written to test_ids.json for the loader (GBM_PG_TEST_IDS). Round-trip checked against the NIfTI.
Output: kaggle_run/pack_pg/ (git-ignored), kaggle_run/pack_pg_manifest.json
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
PACK = PROJECT_ROOT / "kaggle_run" / "pack_pg"
sys.path.insert(0, str(PROJECT_ROOT / "src"))
import predictgbm_io as pg  # noqa: E402

CODE = ["run_improved_aniso.py", "src/__init__.py", "src/radiation_model.py", "src/treatment_aware_pde.py", "src/tmz_pk.py",
        "src/solver_monotone.py", "src/patient_stats.py", "src/run_manifest.py", "src/predictgbm_io.py",
        "src/112_predictgbm_models_dev.py", "output/forecast_validation/dti_atlas_2mm.npz"]


def main():
    if PACK.exists():
        shutil.rmtree(PACK)
    (PACK / "repo").mkdir(parents=True)
    test = pg.test_ids()
    arr, n_rec = {}, 0
    for pid in pg.all_ids():
        arr[f"{pid}|tumor_seg.nii.gz"] = pg.seg_labels(pg._load(pid, "tumor_seg.nii.gz")).astype(np.int8)
        arr[f"{pid}|t1c_bet_mask.nii.gz|bits"] = np.packbits((pg._load(pid, "t1c_bet_mask.nii.gz") > 0).ravel())
        for t in ("wm", "gm", "csf"):
            arr[f"{pid}|{t}_pbmap.nii.gz"] = pg._load(pid, f"{t}_pbmap.nii.gz").astype(np.float16)
        if pid not in test:
            arr[f"{pid}|recurrence_preop.nii.gz"] = pg.seg_labels(pg._load(pid, "recurrence_preop.nii.gz")).astype(np.int8)
            n_rec += 1
    out = PACK / "pg_inputs.npz"
    np.savez_compressed(out, **arr)
    (PACK / "test_ids.json").write_text(json.dumps(sorted(test)))
    z = np.load(out)
    leaked = [k for k in z.files if k.endswith("recurrence_preop.nii.gz") and k.split("|")[0] in test]
    assert not leaked, leaked
    pid = pg.dev_ids()[0]
    assert np.array_equal(z[f"{pid}|tumor_seg.nii.gz"], pg.seg_labels(pg._load(pid, "tumor_seg.nii.gz")).astype(np.int8))
    for rel in CODE:
        (PACK / "repo" / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(PROJECT_ROOT / rel, PACK / "repo" / rel)
    (PACK / "dataset-metadata.json").write_text(json.dumps({"title": "predictgbm-inputs-private",
                                                            "id": "vihansanthosh/predictgbm-inputs-private",
                                                            "licenses": [{"name": "other"}]}, indent=1))
    man = {"script": "113_pack_predictgbm_for_kaggle", "n_patients": len(pg.all_ids()), "n_test": len(test),
           "n_recurrence_packed_dev_only": n_rec, "test_recurrences_packed": 0, "pack_mb": round(out.stat().st_size / 1e6, 1),
           "pack_sha256": hashlib.sha256(out.read_bytes()).hexdigest(), "visibility": "private"}
    (PROJECT_ROOT / "kaggle_run" / "pack_pg_manifest.json").write_text(json.dumps(man, indent=1))
    print(json.dumps(man, indent=1))


if __name__ == "__main__":
    main()
