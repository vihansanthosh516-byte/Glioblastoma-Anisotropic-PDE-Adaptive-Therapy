"""Generate outcome-blind, MRI-only 5 mm slice review montages for CFB-GBM.

This utility reads only t0/t1/t2 T1-Gd paths supplied by the candidate
crosswalk. It does not construct or access GTV paths, score predictions, or
assign a clinical pass/fail status.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import nibabel as nib
import numpy as np
from PIL import Image, ImageDraw, ImageFont


VISITS = ("t0", "t1", "t2")
VISIT_COLORS = {
    "t0": (55, 210, 85),
    "t1": (255, 92, 75),
    "t2": (85, 160, 255),
}
PLANES = (("axial", 2), ("coronal", 1), ("sagittal", 0))
PANEL_SIZE = 112
CELL_PADDING = 8
GRID_COLUMNS = 4
TARGET_SLICE_MM = 5.0
HEADER_ATOL = 1e-4


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_crosswalk(path: Path) -> tuple[list[dict[str, str]], str]:
    digest = sha256_file(path)
    with path.open("r", newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    required = {"patient_id", "gap_days", *(f"{visit}_t1gd_path" for visit in VISITS)}
    if not rows or not required.issubset(rows[0]):
        missing = sorted(required - set(rows[0] if rows else []))
        raise ValueError(f"Crosswalk is empty or missing required columns: {missing}")
    candidates = [row for row in rows if 76 <= float(row["gap_days"]) <= 104]
    ids = [str(int(float(row["patient_id"]))) for row in candidates]
    if len(candidates) != 38 or len(set(ids)) != len(ids):
        raise ValueError(f"Expected 38 unique 76–104 day candidates; found {len(ids)}")
    for row in candidates:
        row["patient_id"] = str(int(float(row["patient_id"])))
    return candidates, digest


def checked_mri_path(row: dict[str, str], visit: str) -> Path:
    path = Path(row[f"{visit}_t1gd_path"]).expanduser()
    name = path.name.casefold()
    if "gtv" in name or not name.endswith(f"_{visit}_t1gd.nii.gz"):
        raise ValueError(f"Not an allowed {visit} T1-Gd path: {path.name}")
    if not path.is_file():
        raise FileNotFoundError(path)
    return path.resolve()


def canonical_header(path: Path) -> tuple[nib.Nifti1Image, dict[str, Any]]:
    image = nib.as_closest_canonical(nib.load(str(path), mmap=True))
    header = {
        "shape": list(image.shape),
        "affine": np.asarray(image.affine, dtype=float).tolist(),
        "spacing_mm": [float(x) for x in image.header.get_zooms()[:3]],
        "axis_codes": list(nib.aff2axcodes(image.affine)),
        "spatial_units": image.header.get_xyzt_units()[0],
    }
    return image, header


def read_patient_images(
    row: dict[str, str],
) -> tuple[dict[str, np.ndarray], dict[str, dict[str, Any]], dict[str, Path]]:
    arrays: dict[str, np.ndarray] = {}
    headers: dict[str, dict[str, Any]] = {}
    paths: dict[str, Path] = {}
    for visit in VISITS:
        path = checked_mri_path(row, visit)
        image, header = canonical_header(path)
        if header["spatial_units"] != "mm":
            raise ValueError(f"Spatial units must be declared mm for patient {row['patient_id']}/{visit}")
        if "unknown" in header["axis_codes"]:
            raise ValueError(f"Unknown axis orientation for patient {row['patient_id']}/{visit}")
        array = np.asarray(image.dataobj, dtype=np.float32)
        if not np.isfinite(array).all():
            raise ValueError(f"Non-finite T1-Gd values for patient {row['patient_id']}/{visit}")
        arrays[visit] = array
        headers[visit] = header
        paths[visit] = path

    reference = headers["t0"]
    for visit in ("t1", "t2"):
        current = headers[visit]
        same_shape = current["shape"] == reference["shape"]
        same_affine = np.allclose(
            np.asarray(current["affine"]),
            np.asarray(reference["affine"]),
            rtol=0.0,
            atol=HEADER_ATOL,
        )
        if not same_shape or not same_affine:
            raise ValueError(
                f"No resampling is performed; {visit} grid differs from t0 for "
                f"patient {row['patient_id']}"
            )
    return arrays, headers, paths


def normalize_shared(arrays: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    # One deterministic shared p1–p99 window per patient keeps visit brightness
    # comparable while avoiding storage of a second copy of each volume.
    sampled = np.concatenate([arrays[v].reshape(-1)[::64] for v in VISITS])
    low, high = (float(x) for x in np.percentile(sampled, [1, 99]))
    if not math.isfinite(low + high) or high <= low:
        raise ValueError("Shared T1-Gd display window is degenerate")
    normalized: dict[str, np.ndarray] = {}
    for visit in VISITS:
        normalized[visit] = np.clip((arrays[visit] - low) / (high - low), 0, 1)
    return normalized, {"p1": low, "p99": high, "sample_stride_voxels": 64}


def slice_indices(size: int, spacing_mm: float) -> tuple[list[int], float]:
    step = max(1, int(round(TARGET_SLICE_MM / spacing_mm)))
    indices = list(range(0, size, step))
    if indices[-1] != size - 1:
        indices.append(size - 1)
    return indices, step * spacing_mm


def display_slice(volume: np.ndarray, axis: int, index: int, spacing_mm: list[float]) -> Image.Image:
    if axis == 2:
        plane = volume[:, :, index].T
        row_spacing, column_spacing = spacing_mm[1], spacing_mm[0]
    elif axis == 1:
        plane = volume[:, index, :].T
        plane = np.flipud(plane)  # superior up in canonical RAS display
        row_spacing, column_spacing = spacing_mm[2], spacing_mm[0]
    else:
        plane = volume[index, :, :].T
        plane = np.flipud(plane)  # superior up in canonical RAS display
        row_spacing, column_spacing = spacing_mm[2], spacing_mm[1]
    pixels = np.uint8(np.clip(plane * 255.0, 0, 255))
    image = Image.fromarray(pixels, mode="L")
    aspect = (plane.shape[1] * column_spacing) / (plane.shape[0] * row_spacing)
    if aspect >= 1:
        display_size = (PANEL_SIZE, max(1, round(PANEL_SIZE / aspect)))
    else:
        display_size = (max(1, round(PANEL_SIZE * aspect)), PANEL_SIZE)
    image = image.resize(display_size, Image.Resampling.BILINEAR)
    panel = Image.new("RGB", (PANEL_SIZE, PANEL_SIZE), (12, 12, 12))
    panel.paste(image.convert("RGB"), ((PANEL_SIZE - image.width) // 2, (PANEL_SIZE - image.height) // 2))
    return panel


def make_plane_montage(
    patient_id: str,
    plane_name: str,
    axis: int,
    arrays: dict[str, np.ndarray],
    headers: dict[str, dict[str, Any]],
    normalization: dict[str, float],
    output_path: Path,
) -> dict[str, Any]:
    shape = arrays["t0"].shape
    spacing = headers["t0"]["spacing_mm"][axis]
    indices, actual_step_mm = slice_indices(shape[axis], spacing)
    tile_w = 3 * PANEL_SIZE + 2 * CELL_PADDING + 8
    tile_h = PANEL_SIZE + 34
    columns = GRID_COLUMNS
    rows = math.ceil(len(indices) / columns)
    width = columns * tile_w + CELL_PADDING
    height = rows * tile_h + 38
    canvas = Image.new("RGB", (width, height), (28, 30, 34))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    axis_code = headers["t0"]["axis_codes"][axis]
    draw.text(
        (CELL_PADDING, 8),
        f"Patient {patient_id} | {plane_name.title()} | t0/t1/t2 T1-Gd | "
        f"~{actual_step_mm:g} mm spacing | shared p1–p99 display window",
        fill=(245, 245, 245),
        font=font,
    )
    for tile_index, slice_index in enumerate(indices):
        x0 = CELL_PADDING + (tile_index % columns) * tile_w
        y0 = 32 + (tile_index // columns) * tile_h
        offset = slice_index * spacing
        draw.text(
            (x0, y0),
            f"{axis_code}-axis index {slice_index} | offset {offset:.1f} mm from voxel origin",
            fill=(225, 225, 225),
            font=font,
        )
        for visit_index, visit in enumerate(VISITS):
            px = x0 + visit_index * (PANEL_SIZE + CELL_PADDING)
            draw.text((px + 3, y0 + 13), visit, fill=VISIT_COLORS[visit], font=font)
            panel = display_slice(arrays[visit], axis, slice_index, headers[visit]["spacing_mm"])
            canvas.paste(panel, (px, y0 + 27))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite existing montage: {output_path}")
    canvas.save(output_path, format="PNG", optimize=True)
    return {
        "file": output_path.name,
        "sha256": sha256_file(output_path),
        "plane": plane_name,
        "axis": axis,
        "axis_code": axis_code,
        "source_spacing_mm": spacing,
        "target_spacing_mm": TARGET_SLICE_MM,
        "nearest_native_grid_interval_mm": actual_step_mm,
        "slice_indices": indices,
        "slice_count": len(indices),
        "image_size_px": list(canvas.size),
        "shared_display_window_p1_p99": normalization,
    }


def process_patient(row: dict[str, str], output_dir: Path) -> dict[str, Any]:
    started = time.perf_counter()
    patient_id = row["patient_id"]
    arrays, headers, paths = read_patient_images(row)
    normalized, display_window = normalize_shared(arrays)
    patient_dir = output_dir / f"patient_{patient_id}"
    montages = [
        make_plane_montage(
            patient_id,
            plane_name,
            axis,
            normalized,
            headers,
            display_window,
            patient_dir / f"patient_{patient_id}_{plane_name}_5mm.png",
        )
        for plane_name, axis in PLANES
    ]
    inputs = [
        {
            "visit": visit,
            "filename": paths[visit].name,
            "sha256": sha256_file(paths[visit]),
            "bytes": paths[visit].stat().st_size,
            "shape_canonical": headers[visit]["shape"],
            "affine_canonical": headers[visit]["affine"],
            "spacing_mm_canonical": headers[visit]["spacing_mm"],
            "axis_codes_canonical": headers[visit]["axis_codes"],
            "spatial_units": headers[visit]["spatial_units"],
        }
        for visit in VISITS
    ]
    result = {
        "patient_id": patient_id,
        "gap_days_from_crosswalk": float(row["gap_days"]),
        "inputs": inputs,
        "shared_display_window": display_window,
        "montages": montages,
        "elapsed_seconds": time.perf_counter() - started,
    }
    del arrays, normalized
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", required=True, type=Path, help="Current CFB v4 candidate crosswalk CSV")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for per-patient montage PNGs and manifest")
    parser.add_argument("--patient-ids", nargs="+", help="Optional subset from the 38 fixed-window candidates")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    crosswalk = args.crosswalk.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    candidates, crosswalk_sha = read_crosswalk(crosswalk)
    selected_ids = {str(int(x)) for x in args.patient_ids} if args.patient_ids else None
    if selected_ids is not None:
        candidate_ids = {row["patient_id"] for row in candidates}
        missing = sorted(selected_ids - candidate_ids)
        if missing:
            raise ValueError(f"Patient IDs are not in the fixed 38-candidate window: {missing}")
        candidates = [row for row in candidates if row["patient_id"] in selected_ids]
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    results = []
    errors = []
    for row in candidates:
        try:
            result = process_patient(row, output_dir)
            results.append(result)
            print(f"rendered patient {row['patient_id']} ({len(results)}/{len(candidates)})", flush=True)
        except Exception as exc:
            errors.append({"patient_id": row["patient_id"], "error_type": type(exc).__name__, "message": str(exc)})
            print(f"failed patient {row['patient_id']}: {type(exc).__name__}", flush=True)
    manifest = {
        "study": "Outcome-blind MRI-only T1-Gd slice display; no target scoring",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": Path(__file__).name,
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "crosswalk": {"filename": crosswalk.name, "sha256": crosswalk_sha},
        "candidate_window_days": [76, 104],
        "selected_patient_ids": [row["patient_id"] for row in candidates],
        "input_roles_read": ["t0_t1gd", "t1_t1gd", "t2_t1gd"],
        "forbidden_inputs_opened": {"gtv": False, "brain_mask": False},
        "geometry": "Canonical axis order for display; no resampling. t1/t2 shapes and affines must match t0 (atol 1e-4 mm; rtol 0).",
        "display_normalization": "One shared per-patient p1–p99 linear window across deterministic every-64th flattened voxel samples from all three visits; clipped to [0,1]. Display-only; no metric or score is computed.",
        "slice_design": "Use the nearest achievable native-grid interval to the 5 mm target along each canonical axis, including both volume ends; offsets are voxel-index times spacing from the voxel origin and the final interval may be shorter. Axial, coronal, sagittal montages show t0, t1, t2 panels side by side, with panel aspect ratios scaled from physical voxel spacing.",
        "display_orientation": "Canonical RAS axes; superior is up in coronal and sagittal panels. No radiological left/right convention is claimed.",
        "clinical_or_registration_pass": None,
        "thresholds": None,
        "results": results,
        "errors": errors,
        "summary": {
            "patients_requested": len(candidates),
            "patients_completed": len(results),
            "montages_written": sum(len(result["montages"]) for result in results),
            "elapsed_seconds": sum(result["elapsed_seconds"] for result in results),
        },
    }
    manifest_path = output_dir / "manifest.json"
    with manifest_path.open("x", encoding="utf-8") as stream:
        stream.write(json.dumps(manifest, indent=2))
    print(json.dumps(manifest["summary"], indent=2))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
