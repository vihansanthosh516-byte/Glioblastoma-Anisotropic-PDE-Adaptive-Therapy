"""MRI-only rigid registration of LUMIERE visits 1 and 3 onto visit 2 (amendment G3; frozen section 5b method/QC).

Reads per patient only: atlas ct1_skull_strip for visits 1, 2, 3 and the visit-2 atlas brain_mask.
Never reads seg_mask, native segmentations or HD-GLIO files. Saves parameters and QC only (no warped images).
Usage: python -m src.fk_lumiere_registration_run [--workers N]
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed

import nibabel as nib
import numpy as np

from src.fk_longitudinal_prep import register_rigid_mri
from src.fk_lumiere_transform_audit import ATLAS_CT1, INVENTORY, PROJECT_ROOT, ZIP_PATH, _visit_dirs

OUT_DIR = PROJECT_ROOT / "output" / "track_b_lumiere_registration_20261110"
ATLAS_MASK = "DeepBraTumIA-segmentation/atlas/skull_strip/brain_mask.nii.gz"
ALLOWED = (ATLAS_CT1, ATLAS_MASK)
EXPECTED_SHAPE = (182, 218, 182)


def _load(zf, name):
    raw = zf.open(name).read()
    img = nib.Nifti1Image.from_bytes(gzip.decompress(raw))
    return np.asarray(img.dataobj, dtype=float), tuple(float(z) for z in img.header.get_zooms()[:3]), \
        hashlib.sha256(raw).hexdigest()


def register_patient(pid: str, visits: list) -> dict:
    with zipfile.ZipFile(ZIP_PATH) as zf:
        dirs = _visit_dirs(zf.namelist(), pid)
        paths = [dirs[(v["week"], v["order"])] for v in visits]
        names = [p + ATLAS_CT1 for p in paths] + [paths[1] + ATLAS_MASK]
        assert all(n.endswith(ALLOWED) for n in names)
        vols, hashes, sanity = [], {}, []
        for n in names:
            arr, zooms, h = _load(zf, n)
            hashes[n] = h
            if arr.shape != EXPECTED_SHAPE:
                sanity.append(f"{n}: shape {arr.shape}")
            if not np.isfinite(arr).all():
                sanity.append(f"{n}: nonfinite values")
            if not np.allclose(zooms, 1.0):
                sanity.append(f"{n}: spacing {zooms}")
            vols.append(arr)
    v1, v2, v3, mask = vols[0], vols[1], vols[2], vols[3] > 0
    if not mask.any():
        sanity.append("visit-2 brain mask empty")
    out = {"patient": pid, "weeks": [v["week"] for v in visits], "input_sha256": hashes, "sanity_failures": sanity}
    if sanity:
        out["passed"] = False
        return out
    for key, moving in (("visit1", v1), ("visit3", v3)):
        r = register_rigid_mri(v2, moving, fixed_support=mask, spacing_mm=(1.0, 1.0, 1.0))
        out[key] = {"params": [float(x) for x in r["params"]], "ncc_before": float(r["ncc_before"]),
                    "ncc_after": float(r["ncc_after"]), "qc": r["qc"]}
    out["passed"] = bool(out["visit1"]["qc"]["passed"] and out["visit3"]["qc"]["passed"])
    return out


def main(workers: int):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    cases = OUT_DIR / "cases"
    cases.mkdir(exist_ok=True)
    inv = json.loads(INVENTORY.read_text())
    todo = {pid: r["visits"] for pid, r in inv["patients"].items() if not (cases / f"{pid}.json").exists()}
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(register_patient, pid, v): pid for pid, v in todo.items()}
        for f in as_completed(futs):
            res = f.result()
            (cases / f"{futs[f]}.json").write_text(json.dumps(res, indent=1))
            print(futs[f], res["passed"], flush=True)
    allres = [json.loads((cases / f"{pid}.json").read_text()) for pid in inv["patients"]]
    summary = {"method": "paper/trackB_lumiere_geometry_amendment_2026-11-10.md G3; "
                         "paper/track_b_no_contact_dataset_plan_2026-10-31.md section 5b",
               "read": "atlas ct1_skull_strip visits 1-3 and visit-2 atlas brain_mask only; no tumor labels",
               "registration_code_sha256": hashlib.sha256((PROJECT_ROOT / "src" / "fk_longitudinal_prep.py")
                                                          .read_bytes()).hexdigest(),
               "n_patients": len(allres), "n_passed": sum(r["passed"] for r in allres),
               "failed": {r["patient"]: (r["sanity_failures"] or
                                         {k: r[k]["qc"]["reasons"] for k in ("visit1", "visit3")
                                          if not r[k]["qc"]["passed"]}) for r in allres if not r["passed"]},
               "cases": allres}
    (OUT_DIR / "registration_qc.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps({k: v for k, v in summary.items() if k != "cases"}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    main(ap.parse_args().workers)
