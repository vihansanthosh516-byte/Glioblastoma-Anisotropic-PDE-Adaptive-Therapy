"""Render paginated, MRI-only 5 mm T1-Gd review panels for CFB-GBM.

Each slice page contains at most eight slices. Every slice shows t0, t1, t2
side by side at 320 px per visit. This is a sampled visual aid, not registration
certification; no target or brain-mask image is opened or scored.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import cfb_mri_5mm_review as review


VISIT_PANEL_SIZE = 320
CELL_PADDING = 10
MAX_SLICES_PER_PAGE = 8
PAGE_COLUMNS = 2


def hash_file(path: Path) -> str:
    return review.sha256_file(path)


def shared_display_window(arrays: dict[str, np.ndarray]) -> dict[str, float]:
    sampled = np.concatenate([arrays[visit].reshape(-1)[::64] for visit in review.VISITS])
    low, high = (float(value) for value in np.percentile(sampled, [1, 99]))
    if not math.isfinite(low + high) or high <= low:
        raise ValueError("Shared T1-Gd display window is degenerate")
    return {"p1": low, "p99": high, "sample_stride_voxels": 64}


def display_raw_slice(
    volume: np.ndarray,
    axis: int,
    index: int,
    spacing_mm: list[float],
    display_window: dict[str, float],
) -> Image.Image:
    if axis == 2:
        plane = volume[:, :, index].T
        row_spacing, column_spacing = spacing_mm[1], spacing_mm[0]
    elif axis == 1:
        plane = np.flipud(volume[:, index, :].T)
        row_spacing, column_spacing = spacing_mm[2], spacing_mm[0]
    else:
        plane = np.flipud(volume[index, :, :].T)
        row_spacing, column_spacing = spacing_mm[2], spacing_mm[1]
    low, high = display_window["p1"], display_window["p99"]
    pixels = np.uint8(np.clip((plane - low) / (high - low), 0, 1) * 255.0)
    image = Image.fromarray(pixels, mode="L")
    aspect = (plane.shape[1] * column_spacing) / (plane.shape[0] * row_spacing)
    if aspect >= 1:
        display_size = (VISIT_PANEL_SIZE, max(1, round(VISIT_PANEL_SIZE / aspect)))
    else:
        display_size = (max(1, round(VISIT_PANEL_SIZE * aspect)), VISIT_PANEL_SIZE)
    image = image.resize(display_size, Image.Resampling.BILINEAR)
    panel = Image.new("RGB", (VISIT_PANEL_SIZE, VISIT_PANEL_SIZE), (12, 12, 12))
    panel.paste(image.convert("RGB"), ((VISIT_PANEL_SIZE - image.width) // 2, (VISIT_PANEL_SIZE - image.height) // 2))
    return panel


def make_page(
    patient_id: str,
    plane_name: str,
    axis: int,
    slice_indices: list[int],
    page_number: int,
    actual_step_mm: float,
    arrays: dict[str, Any],
    headers: dict[str, dict[str, Any]],
    display_window: dict[str, float],
    output_path: Path,
) -> dict[str, Any]:
    tile_width = 3 * VISIT_PANEL_SIZE + 4 * CELL_PADDING
    tile_height = VISIT_PANEL_SIZE + 44
    rows = math.ceil(len(slice_indices) / PAGE_COLUMNS)
    width = PAGE_COLUMNS * tile_width + CELL_PADDING
    height = 42 + rows * tile_height + CELL_PADDING
    canvas = Image.new("RGB", (width, height), (28, 30, 34))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.load_default()
    axis_code = headers["t0"]["axis_codes"][axis]
    draw.text(
        (CELL_PADDING, 8),
        f"Patient {patient_id} | {plane_name.title()} | T1-Gd | page {page_number:03d} | ~{actual_step_mm:g} mm slices",
        fill=(245, 245, 245),
        font=font,
    )
    for local_index, slice_index in enumerate(slice_indices):
        col = local_index % PAGE_COLUMNS
        row = local_index // PAGE_COLUMNS
        x0 = CELL_PADDING + col * tile_width
        y0 = 34 + row * tile_height
        offset_mm = slice_index * headers["t0"]["spacing_mm"][axis]
        draw.text(
            (x0, y0),
            f"{axis_code}-axis index {slice_index} | {offset_mm:.1f} mm from voxel origin",
            fill=(225, 225, 225),
            font=font,
        )
        for visit_index, visit in enumerate(review.VISITS):
            panel_x = x0 + visit_index * (VISIT_PANEL_SIZE + CELL_PADDING)
            draw.text(
                (panel_x + 3, y0 + 13), visit,
                fill=review.VISIT_COLORS[visit], font=font,
            )
            panel = display_raw_slice(
                arrays[visit], axis, slice_index, headers[visit]["spacing_mm"],
                display_window,
            )
            canvas.paste(panel, (panel_x, y0 + 27))
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite page: {output_path}")
    canvas.save(output_path, format="PNG", optimize=True)
    return {
        "file": output_path.name,
        "sha256": hash_file(output_path),
        "bytes": output_path.stat().st_size,
        "plane": plane_name,
        "axis": axis,
        "axis_code": axis_code,
        "page_number": page_number,
        "page_slice_count": len(slice_indices),
        "slice_indices": slice_indices,
        "source_spacing_mm": headers["t0"]["spacing_mm"][axis],
        "target_slice_interval_mm": review.TARGET_SLICE_MM,
        "nearest_native_grid_interval_mm": actual_step_mm,
        "visit_panel_size_px": VISIT_PANEL_SIZE,
        "image_size_px": [width, height],
    }


def render_patient(
    row: dict[str, str],
    staging_dir: Path,
    panel_limit: int,
) -> dict[str, Any]:
    patient_id = row["patient_id"]
    paths = {visit: review.checked_mri_path(row, visit) for visit in review.VISITS}
    before_hashes = {visit: hash_file(path) for visit, path in paths.items()}
    arrays, headers, read_paths = review.read_patient_images(row)
    if paths != read_paths:
        raise RuntimeError("MRI source paths changed during preflight")
    display_window = shared_display_window(arrays)
    after_hashes = {visit: hash_file(path) for visit, path in paths.items()}
    if before_hashes != after_hashes:
        raise RuntimeError("MRI source file changed during image load")

    old_panel_size = review.PANEL_SIZE
    review.PANEL_SIZE = VISIT_PANEL_SIZE
    page_records = []
    try:
        for plane_name, axis in review.PLANES:
            indices, actual_step_mm = review.slice_indices(
                arrays["t0"].shape[axis], headers["t0"]["spacing_mm"][axis]
            )
            for start in range(0, len(indices), panel_limit):
                page_indices = indices[start : start + panel_limit]
                page_number = start // panel_limit + 1
                page_path = staging_dir / (
                    f"patient_{patient_id}_{plane_name}_5mm_page_{page_number:03d}.png"
                )
                page_records.append(
                    make_page(
                        patient_id, plane_name, axis, page_indices, page_number,
                        actual_step_mm, arrays, headers, display_window, page_path,
                    )
                )
    finally:
        review.PANEL_SIZE = old_panel_size

    result = {
        "patient_id": patient_id,
        "gap_days_from_crosswalk": float(row["gap_days"]),
        "inputs": [
            {
                "visit": visit,
                "filename": paths[visit].name,
                "sha256": before_hashes[visit],
                "bytes": paths[visit].stat().st_size,
                "shape_canonical": headers[visit]["shape"],
                "affine_canonical": headers[visit]["affine"],
                "spacing_mm_canonical": headers[visit]["spacing_mm"],
                "axis_codes_canonical": headers[visit]["axis_codes"],
                "spatial_units": headers[visit]["spatial_units"],
            }
            for visit in review.VISITS
        ],
        "shared_display_window": display_window,
        "pages": page_records,
    }
    del arrays
    return result


def run(args: argparse.Namespace) -> dict[str, Any]:
    crosswalk = args.crosswalk.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output_dir}")
    candidates, crosswalk_hash = review.read_crosswalk(crosswalk)
    selected = {str(int(patient_id)) for patient_id in args.patient_ids}
    by_id = {row["patient_id"]: row for row in candidates}
    missing = sorted(selected - set(by_id))
    if missing:
        raise ValueError(f"Patient IDs not in fixed candidate window: {missing}")
    rows = [by_id[patient_id] for patient_id in sorted(selected, key=int)]

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=False, exist_ok=False)
    results = [render_patient(row, output_dir, args.max_slices_per_page) for row in rows]
    script_path = Path(__file__).resolve()
    manifest = {
            "study": "MRI-only paginated T1-Gd visual review aid; no target scoring",
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "script": script_path.name,
            "script_sha256": hash_file(script_path),
            "base_reader_sha256": hash_file(Path(review.__file__).resolve()),
            "crosswalk": {"filename": crosswalk.name, "sha256": crosswalk_hash},
            "crosswalk_handling": "The crosswalk CSV is parsed and contains GTV-path metadata; GTV paths are not dereferenced or opened. Only selected t0/t1/t2 T1-Gd NIfTI volumes are opened and used.",
            "selected_patient_ids": [row["patient_id"] for row in rows],
            "input_roles_read": ["t0_t1gd", "t1_t1gd", "t2_t1gd"],
            "forbidden_inputs_opened": {"gtv": False, "brain_mask": False},
            "geometry": "Canonical display axes; no resampling. T1-Gd visit shapes and affines must match t0 within 1e-4 mm.",
            "display_normalization": "Per-patient shared p1-p99 window from deterministic every-64th flattened voxel samples across all three visits; clipped to [0,1]. Display only.",
            "slice_design": "Nearest native-grid interval to 5 mm per axis; includes both volume ends, so final interval can be shorter. Pages contain at most the configured slice count; each slice shows t0, t1, t2 at 320 px per visit.",
            "display_orientation": "Canonical RAS axes; superior is up in coronal and sagittal panels. No radiological left/right convention is claimed.",
            "clinical_or_registration_pass": None,
            "thresholds": None,
            "metrics": None,
            "results": results,
            "summary": {
                "patients_completed": len(results),
                "pages_written": sum(len(result["pages"]) for result in results),
                "slices_rendered": sum(page["page_slice_count"] for result in results for page in result["pages"]),
            },
    }
    manifest_path = output_dir / "manifest.json"
    with manifest_path.open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2)
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--patient-ids", nargs="+", required=True)
    parser.add_argument("--max-slices-per-page", type=int, default=8)
    args = parser.parse_args()
    if not 1 <= args.max_slices_per_page <= MAX_SLICES_PER_PAGE:
        raise ValueError("max-slices-per-page must be between 1 and 8")
    result = run(args)
    print(json.dumps(result["summary"], indent=2))


if __name__ == "__main__":
    main()
