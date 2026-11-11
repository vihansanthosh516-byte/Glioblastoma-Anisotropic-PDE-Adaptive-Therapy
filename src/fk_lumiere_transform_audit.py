"""Metadata-only audit of the LUMIERE coverage transform chain (Track B). No voxel arrays, intensities or labels.

Reads, for the 28 planned triplets in output/track_b_lumiere_inventory_20261031/header_check.json:
native ``CT1.nii.gz`` header, atlas ``ct1_skull_strip.nii.gz`` header (348 bytes each) and
``DeepBraTumIA-segmentation/native/transformation/ct1.tfm`` (ITK text). Nothing else is opened.
Usage: python -m src.fk_lumiere_transform_audit
"""
from __future__ import annotations

import gzip
import json
import re
import struct
import zipfile
from pathlib import Path as _Path

import numpy as np

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
ZIP_PATH = PROJECT_ROOT / "data" / "external" / "lumiere" / "Imaging-v202211.zip"
INVENTORY = PROJECT_ROOT / "output" / "track_b_lumiere_inventory_20261031" / "header_check.json"
OUT_DIR = PROJECT_ROOT / "output" / "track_b_lumiere_transform_audit_20261110"
NATIVE_CT1 = "CT1.nii.gz"
ATLAS_CT1 = "DeepBraTumIA-segmentation/atlas/skull_strip/ct1_skull_strip.nii.gz"
TFM = "DeepBraTumIA-segmentation/native/transformation/ct1.tfm"
ALLOWED = (NATIVE_CT1, ATLAS_CT1, TFM)
LPS = np.diag([-1.0, -1.0, 1.0, 1.0])            # RAS <-> LPS (self-inverse)
AFFINE_ITK_TYPES = ("MatrixOffsetTransformBase", "AffineTransform", "Euler3DTransform",
                    "VersorRigid3DTransform", "Similarity3DTransform", "ScaleVersor3DTransform")


# ---------------------------------------------------------------- NIfTI-1 header (348 bytes, no voxels)
def parse_header_geometry(raw: bytes) -> dict:
    if len(raw) < 348:
        raise ValueError("short NIfTI header")
    e = "<" if struct.unpack("<i", raw[:4])[0] == 348 else ">"
    if struct.unpack(e + "i", raw[:4])[0] != 348:
        raise ValueError("not a NIfTI-1 header")
    dim = struct.unpack(e + "8h", raw[40:56])
    pix = struct.unpack(e + "8f", raw[76:108])
    units = raw[123]
    qcode, scode = struct.unpack(e + "2h", raw[252:256])
    b, c, d, qx, qy, qz = struct.unpack(e + "6f", raw[256:280])
    srow = np.array(struct.unpack(e + "12f", raw[280:328]), float).reshape(3, 4)
    a2 = 1.0 - (b * b + c * c + d * d)
    a = float(np.sqrt(a2)) if a2 > 0 else 0.0
    rot = np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                    [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                    [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-c*c-b*b]])
    qfac = -1.0 if pix[0] < 0 else 1.0
    q = np.eye(4)
    q[:3, :3] = rot @ np.diag([pix[1], pix[2], pix[3] * qfac])
    q[:3, 3] = (qx, qy, qz)
    s = np.eye(4)
    s[:3] = srow
    chosen = "sform" if scode > 0 else ("qform" if qcode > 0 else "none")
    affine = s if chosen == "sform" else q if chosen == "qform" else None
    return {"shape": [int(x) for x in dim[1:1 + dim[0]]][:3], "spacing_mm": [float(x) for x in pix[1:4]],
            "spatial_units_code": int(units & 0x07), "qform_code": int(qcode), "sform_code": int(scode),
            "qform": q.tolist(), "sform": s.tolist(), "affine_used": chosen,
            "affine": None if affine is None else affine.tolist(),
            "qform_sform_agree": bool(qcode > 0 and scode > 0 and np.allclose(q, s, atol=1e-3))}


# ---------------------------------------------------------------- ITK transform text
def parse_itk_tfm(text: str) -> dict:
    """Parse an ITK Insight Transform File. Returns each component; never guesses missing parts."""
    if not text.startswith("#Insight Transform File"):
        raise ValueError("not an ITK transform file")
    blocks, cur = [], None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("Transform:"):
            cur = {"type": line.split(":", 1)[1].strip()}
            blocks.append(cur)
        elif cur is not None and line.startswith("Parameters:"):
            cur["parameters"] = [float(x) for x in line.split(":", 1)[1].split()]
        elif cur is not None and line.startswith("FixedParameters:"):
            cur["fixed_parameters"] = [float(x) for x in line.split(":", 1)[1].split()]
    leaves = [b for b in blocks if not b["type"].startswith("CompositeTransform")]
    return {"blocks": blocks, "leaves": leaves}


def itk_leaf_to_4x4(leaf: dict) -> np.ndarray:
    """4x4 for an ITK MatrixOffset-family leaf in LPS mm: y = A (x - c) + c + t. Raises for anything else."""
    kind = leaf["type"].split("_")[0]
    if kind not in AFFINE_ITK_TYPES or kind != "MatrixOffsetTransformBase" and kind != "AffineTransform":
        raise ValueError(f"unsupported or non-matrix ITK transform: {leaf['type']}")
    p, c = leaf.get("parameters"), leaf.get("fixed_parameters")
    if p is None or len(p) != 12 or c is None or len(c) != 3:
        raise ValueError("matrix transform needs 12 parameters and a 3-value centre")
    A = np.array(p[:9]).reshape(3, 3)
    t, c = np.array(p[9:]), np.array(c)
    m = np.eye(4)
    m[:3, :3] = A
    m[:3, 3] = t + c - A @ c
    return m


def tfm_to_ras_4x4(leaves) -> np.ndarray:
    """Compose ITK leaves (applied in ITK order: last listed is applied first) and convert to RAS."""
    m = np.eye(4)
    for leaf in leaves:                   # CompositeTransform applies the last-added transform first
        m = m @ itk_leaf_to_4x4(leaf)
    return LPS @ m @ LPS


# ---------------------------------------------------------------- composition (voxel index maps)
def downsample_to_1mm(factor: int = 2, offset: float = 0.5) -> np.ndarray:
    """Scoring voxel i (block of `factor`) -> 1 mm voxel centre factor*i + offset. Offset is a protocol choice."""
    m = np.diag([float(factor)] * 3 + [1.0])
    m[:3, 3] = offset
    return m


def rigid_voxel_map(params, *, shape, spacing_mm) -> np.ndarray:
    """4x4 output-to-input voxel map identical to fk_longitudinal_prep.apply_rigid."""
    from src.fk_longitudinal_prep import _rotation
    sp = np.asarray(spacing_mm, float)
    rot = _rotation(np.asarray(params[:3], float))
    centre = (np.asarray(shape, float) - 1) / 2.0
    mat = np.diag(1 / sp) @ rot @ np.diag(sp)
    m = np.eye(4)
    m[:3, :3] = mat
    m[:3, 3] = centre + np.asarray(params[3:], float) / sp - mat @ centre
    return m


def chain_output_to_input(*, atlas_affine, native_affine, atlas_to_native_ras, rigid=None,
                          scoring_to_atlas=None) -> np.ndarray:
    """Scoring voxel -> native voxel. Every factor is an exact 4x4; the direction of ``atlas_to_native_ras``
    must already be verified (it maps atlas physical points to native physical points)."""
    m = np.eye(4) if scoring_to_atlas is None else np.asarray(scoring_to_atlas, float)
    if rigid is not None:
        m = np.asarray(rigid, float) @ m          # fixed (visit-2 atlas) voxel -> moving (visit-3 atlas) voxel
    m = np.asarray(atlas_affine, float) @ m       # atlas voxel -> atlas RAS mm
    m = np.asarray(atlas_to_native_ras, float) @ m
    return np.linalg.inv(np.asarray(native_affine, float)) @ m


# ---------------------------------------------------------------- audit
def _visit_dirs(names, pid):
    pat = re.compile(rf"Imaging/{pid}/(week-(\d+)(?:-(\d+))?)/")
    out = {}
    for n in names:
        m = pat.search(n)
        if m:
            out[(int(m[2]), int(m[3] or 0))] = f"Imaging/{pid}/{m[1]}/"
    return out


def run_audit(zip_path=ZIP_PATH, inventory=INVENTORY) -> dict:
    inv = json.loads(_Path(inventory).read_text())
    rows, opened = {}, []
    with zipfile.ZipFile(zip_path) as zf:
        names = zf.namelist()
        for pid, rec in inv["patients"].items():
            dirs = _visit_dirs(names, pid)
            vrows = []
            for i, v in enumerate(rec["visits"]):
                d = dirs.get((v["week"], v["order"]))
                row = {"visit": i + 1, "week": v["week"], "order": v["order"], "dir_found": d is not None}
                if d is not None:
                    for key, rel in (("native", NATIVE_CT1), ("atlas", ATLAS_CT1)):
                        name = d + rel
                        if name in names:
                            with zf.open(name) as fh, gzip.GzipFile(fileobj=fh) as gz:
                                row[key] = parse_header_geometry(gz.read(348))
                            opened.append(name)
                        else:
                            row[key] = None
                    name = d + TFM
                    if name in names:
                        text = zf.read(name).decode("utf-8", "replace")
                        opened.append(name)
                        tf = parse_itk_tfm(text)
                        row["tfm"] = {"blocks": [{k: b.get(k) for k in ("type", "parameters", "fixed_parameters")}
                                                 for b in tf["blocks"]],
                                      "leaf_types": [b["type"] for b in tf["leaves"]]}
                        try:
                            A = np.array(tf["leaves"][0]["parameters"][:9]).reshape(3, 3) if len(tf["leaves"]) == 1 else None
                            row["tfm"]["single_matrix_leaf"] = len(tf["leaves"]) == 1
                            row["tfm"]["representable_4x4"] = bool(tf["leaves"]) and all(
                                b["type"].split("_")[0] in ("MatrixOffsetTransformBase", "AffineTransform")
                                for b in tf["leaves"])
                            if A is not None:
                                sv = np.linalg.svd(A, compute_uv=False)
                                row["tfm"]["singular_values"] = sv.tolist()
                                row["tfm"]["det"] = float(np.linalg.det(A))
                                row["tfm"]["rigid_within_1e-3"] = bool(np.allclose(sv, 1, atol=1e-3))
                        except Exception as exc:          # noqa: BLE001 - recorded, not hidden
                            row["tfm"]["error"] = repr(exc)
                    else:
                        row["tfm"] = None
                vrows.append(row)
            rows[pid] = vrows
    assert all(any(n.endswith(a) for a in ALLOWED) for n in opened), "opened a file outside the allowed set"
    return summarise(rows, opened)


def summarise(rows, opened) -> dict:
    flat = [v for vs in rows.values() for v in vs]
    leaf = sorted({t for v in flat if v.get("tfm") for t in v["tfm"]["leaf_types"]})
    return {
        "read": "NIfTI-1 headers (348 bytes) of native CT1 and atlas CT1, and ct1.tfm text only; "
                "no voxel arrays, intensities or labels",
        "n_patients": len(rows), "n_visits": len(flat), "n_files_opened": len(opened),
        "missing": {k: sum(1 for v in flat if not v.get(k)) for k in ("native", "atlas", "tfm")},
        "native_shapes": sorted({str(v["native"]["shape"]) for v in flat if v.get("native")}),
        "native_spacings": sorted({str([round(s, 3) for s in v["native"]["spacing_mm"]]) for v in flat if v.get("native")}),
        "native_affine_used": sorted({v["native"]["affine_used"] for v in flat if v.get("native")}),
        "native_qform_sform_disagree": sum(1 for v in flat if v.get("native") and v["native"]["sform_code"] > 0
                                           and v["native"]["qform_code"] > 0 and not v["native"]["qform_sform_agree"]),
        "atlas_shapes": sorted({str(v["atlas"]["shape"]) for v in flat if v.get("atlas")}),
        "atlas_affines_distinct": len({json.dumps(np.round(v["atlas"]["affine"], 4).tolist())
                                       for v in flat if v.get("atlas") and v["atlas"]["affine"]}),
        "tfm_leaf_types": leaf,
        "tfm_all_representable_4x4": all(v["tfm"].get("representable_4x4") for v in flat if v.get("tfm")),
        "tfm_rigid_count": sum(1 for v in flat if v.get("tfm") and v["tfm"].get("rigid_within_1e-3")),
        "tfm_singular_value_range": [min(min(v["tfm"]["singular_values"]) for v in flat if v.get("tfm", {}) and "singular_values" in v["tfm"]),
                                     max(max(v["tfm"]["singular_values"]) for v in flat if v.get("tfm", {}) and "singular_values" in v["tfm"])],
        "patients": rows,
    }


if __name__ == "__main__":
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    res = run_audit()
    (OUT_DIR / "transform_audit.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "patients"}, indent=1))
