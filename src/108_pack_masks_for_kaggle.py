#!/usr/bin/env python3
"""Script 108: build the private Kaggle dataset for the GPU forecast runs (GRAND_PLAN L3, L2a).

What is uploaded (author approved a PRIVATE Kaggle dataset of derived masks on 2026-10-08):
  mu_masks_1mm.npz   binary masks only, no MRI intensities: core (labels 1+3) and brain (> 0) for every scan in a
                     primary forecast pair, 1 mm, cropped to 240 x 240 x 154, bit-packed (np.packbits)
  repo/              the code script 100 needs (solver, monotone solver, stats), the pair / split manifests, the
                     cohort treatment schedules and the 2 mm DTI atlas, in the repo layout so PROJECT_ROOT works
Check: every packed mask is unpacked and compared with the NIfTI path (script 100 mask_native / brain_native).
Output: kaggle_run/pack/ (git-ignored) with dataset-metadata.json; kaggle_run/pack_manifest.json (tracked)
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import sys
from pathlib import Path as _Path

import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
PACK = PROJECT_ROOT / "kaggle_run" / "pack"
sys.path.insert(0, str(PROJECT_ROOT))

spec = importlib.util.spec_from_file_location("m100", PROJECT_ROOT / "src" / "100_pde_manifest.py")
m100 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m100)

CODE = ["run_improved_aniso.py", "src/__init__.py", "src/radiation_model.py", "src/treatment_aware_pde.py",
        "src/tmz_pk.py", "src/tmz_pk.py", "src/100_pde_manifest.py", "src/solver_monotone.py", "src/patient_stats.py", "src/run_manifest.py"]
DATA = ["data/manifests/forecast_pairs_mu.csv", "data/manifests/split_mu.csv", "output/mu_glioma_cohort.json",
        "output/forecast_validation/dti_atlas_2mm.npz"]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    assert not os.environ.get("GBM_MASK_PACK"), "unset GBM_MASK_PACK: masks must be read from the NIfTI files"
    if PACK.exists():
        shutil.rmtree(PACK)
    (PACK / "repo").mkdir(parents=True)
    pairs = pd.read_csv(PROJECT_ROOT / "data" / "manifests" / "forecast_pairs_mu.csv")
    pairs = pairs[pairs["in_primary"]]
    scans = sorted({(r.patient_id, int(t)) for r in pairs.itertuples() for t in (r.tp_in, r.tp_out)})
    arrays = {}
    for pid, tp in scans:
        arrays[f"core|{pid}|{tp}"] = np.packbits(m100.mask_native(pid, tp).ravel())
        arrays[f"brain|{pid}|{tp}"] = np.packbits(m100.brain_native(pid, tp).ravel())
    out = PACK / "mu_masks_1mm.npz"
    np.savez_compressed(out, **arrays)
    # round-trip check through the loader used on Kaggle
    os.environ["GBM_MASK_PACK"] = str(out)
    m100._PACK = None
    bad = 0
    for pid, tp in scans:
        a = m100._packed("core", pid, tp)
        b = m100._packed("brain", pid, tp)
        os.environ.pop("GBM_MASK_PACK")
        bad += int(not np.array_equal(a, m100.mask_native(pid, tp))) + int(not np.array_equal(b, m100.brain_native(pid, tp)))
        os.environ["GBM_MASK_PACK"] = str(out)
    os.environ.pop("GBM_MASK_PACK")
    assert bad == 0, f"{bad} masks differ after packing"
    files = {}
    for rel in CODE + DATA:
        dst = PACK / "repo" / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(PROJECT_ROOT / rel, dst)
        files[rel] = sha(dst)
    (PACK / "dataset-metadata.json").write_text(json.dumps(
        {"title": "gbm-forecast-masks-private", "id": "vihansanthosh/gbm-forecast-masks-private",
         "licenses": [{"name": "other"}]}, indent=1))
    man = {"script": "108_pack_masks_for_kaggle", "n_scans": len(scans), "n_patients": len({p for p, _ in scans}),
           "n_arrays": len(arrays), "pack_mb": round(out.stat().st_size / 1e6, 1), "pack_sha256": sha(out),
           "roundtrip_mismatches": bad, "contents": "binary core (labels 1+3) and brain masks only; no MRI intensities",
           "files_sha256": files, "visibility": "private"}
    (PROJECT_ROOT / "kaggle_run" / "pack_manifest.json").write_text(json.dumps(man, indent=1))
    print(json.dumps({k: v for k, v in man.items() if k != "files_sha256"}, indent=1))


if __name__ == "__main__":
    main()
