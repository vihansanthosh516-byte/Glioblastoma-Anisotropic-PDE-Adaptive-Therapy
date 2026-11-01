"""Header-only LUMIERE reader for Track B (exploratory). Never reads image or mask voxel arrays.

Each NIfTI-1 file is opened inside the release zip and only its first 348 decompressed bytes (the header)
are read. Usage: python -m src.fk_lumiere_headers  (writes output/track_b_lumiere_inventory_20261031/header_check.json)
"""
from __future__ import annotations

import gzip
import json
import re
import struct
import zipfile
from pathlib import Path as _Path

from src.fk_longitudinal_prep import Visit, primary_triplets

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
ZIP_PATH = PROJECT_ROOT / "data" / "external" / "lumiere" / "Imaging-v202211.zip"
OUT_DIR = PROJECT_ROOT / "output" / "track_b_lumiere_inventory_20261031"

SEQUENCES = ("CT1.nii.gz", "T1.nii.gz", "T2.nii.gz", "FLAIR.nii.gz")
ATLAS = "DeepBraTumIA-segmentation/atlas/"
LABEL_MAP = ATLAS + "segmentation/seg_mask.nii.gz"
BRAIN_MASK = ATLAS + "skull_strip/brain_mask.nii.gz"
ATLAS_IMAGES = tuple(ATLAS + f"skull_strip/{s}_skull_strip.nii.gz" for s in ("ct1", "t1", "t2", "flair"))
HEADER_FILES = (LABEL_MAP, BRAIN_MASK) + ATLAS_IMAGES
_STUDY = re.compile(r"(Patient-\d+)/week-(\d+)(?:-(\d+))?/(.*)$")
NIFTI1_HEADER_BYTES = 348


def parse_nifti1_header(raw: bytes) -> dict:
    """Parse the fields Track B needs from a NIfTI-1 header (no voxel data)."""
    if len(raw) < NIFTI1_HEADER_BYTES:
        raise ValueError("short NIfTI header")
    endian = "<" if struct.unpack("<i", raw[:4])[0] == 348 else ">"
    if struct.unpack(endian + "i", raw[:4])[0] != 348:
        raise ValueError("not a NIfTI-1 header")
    dim = struct.unpack(endian + "8h", raw[40:56])
    datatype = struct.unpack(endian + "h", raw[70:72])[0]
    pixdim = struct.unpack(endian + "8f", raw[76:108])
    qform, sform = struct.unpack(endian + "2h", raw[252:256])
    magic = raw[344:348].rstrip(b"\x00").decode("ascii", "replace")
    n = dim[0]
    return {"shape": list(dim[1:1 + n]), "spacing_mm": [round(float(p), 4) for p in pixdim[1:1 + min(n, 3)]],
            "datatype": int(datatype), "qform_code": int(qform), "sform_code": int(sform), "magic": magic}


def read_header(zf: zipfile.ZipFile, name: str) -> dict:
    with zf.open(name) as fh, gzip.GzipFile(fileobj=fh) as gz:
        return parse_nifti1_header(gz.read(NIFTI1_HEADER_BYTES))


def list_visits(zf: zipfile.ZipFile) -> tuple[list[Visit], dict]:
    files: dict[tuple[str, int, int], set[str]] = {}
    prefix: dict[tuple[str, int, int], str] = {}
    for n in zf.namelist():
        m = _STUDY.search(n)
        if m:
            key = (m[1], int(m[2]), int(m[3] or 0))
            files.setdefault(key, set()).add(m[4])
            prefix[key] = n[: len(n) - len(m[4])]
    visits = [Visit(p, w, o, all(s in f for s in SEQUENCES), LABEL_MAP in f, BRAIN_MASK in f)
              for (p, w, o), f in files.items()]
    return visits, prefix


def header_check(zip_path=ZIP_PATH, *, min_gap_weeks=11, max_gap_weeks=15) -> dict:
    with zipfile.ZipFile(zip_path) as zf:
        visits, prefix = list_visits(zf)
        tri = primary_triplets(visits, min_gap_weeks=min_gap_weeks, max_gap_weeks=max_gap_weeks)
        patients, problems, grids = {}, [], {}
        for pid, vs in tri["eligible"].items():
            rows = []
            for v in vs:
                base = prefix[(v.patient_id, v.week, v.order)]
                heads = {}
                for f in HEADER_FILES:
                    try:
                        heads[f] = read_header(zf, base + f)
                    except (KeyError, ValueError, OSError) as exc:
                        problems.append(f"{pid} week {v.week}: {f}: {exc}")
                shapes = {tuple(h["shape"][:3]) for h in heads.values()}
                spacings = {tuple(h["spacing_mm"]) for h in heads.values()}
                if len(shapes) > 1 or len(spacings) > 1:
                    problems.append(f"{pid} week {v.week}: atlas files disagree on grid")
                for s, sp in zip(shapes, spacings):
                    grids[str((s, sp))] = grids.get(str((s, sp)), 0) + 1
                rows.append({"week": v.week, "order": v.order, "n_headers_read": len(heads),
                             "grid": [list(s) for s in shapes], "spacing_mm": [list(s) for s in spacings],
                             "label_datatype": heads.get(LABEL_MAP, {}).get("datatype")})
            patients[pid] = {"visits": rows, "gap_weeks": vs[2].week - vs[1].week}
    return {"source": str(_Path(zip_path).name), "read": "NIfTI-1 headers only (first 348 decompressed bytes); "
            "no image or mask voxel arrays read",
            "files_per_visit": list(HEADER_FILES), "n_eligible": len(patients),
            "n_excluded": len(tri["excluded"]), "grid_counts": grids, "problems": problems, "patients": patients}


def main():
    out = OUT_DIR / "header_check.json"
    if out.exists():
        raise SystemExit(f"{out} exists; refusing to overwrite")
    rec = header_check()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1))
    print({k: rec[k] for k in ("n_eligible", "n_excluded", "grid_counts", "problems")})


if __name__ == "__main__":
    main()
