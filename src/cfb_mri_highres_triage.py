"""Create larger MRI-only review panels and descriptive registration triage.

Reads only t0/t1/t2 T1-Gd images from the fixed CFB candidate crosswalk.
Metrics describe image support and intensity similarity; they do not assign
registration status or authorize target-mask access.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import cfb_mri_5mm_review as review


def descriptive_pair_metrics(reference: np.ndarray, moving: np.ndarray) -> dict[str, Any]:
    """Return unthresholded, MRI-only pair metrics on a shared native grid."""
    if reference.shape != moving.shape:
        raise ValueError("MRI shapes differ; no resampling is performed")
    ref_ok = np.isfinite(reference)
    mov_ok = np.isfinite(moving)
    if not ref_ok.all() or not mov_ok.all():
        raise ValueError("MRI contains non-finite values")
    ref_support = reference != 0
    mov_support = moving != 0
    intersection = int(np.count_nonzero(ref_support & mov_support))
    support_denominator = int(ref_support.sum()) + int(mov_support.sum())
    dice = 2.0 * intersection / support_denominator if support_denominator else None
    sample_count = 0
    if intersection < 2:
        correlation = None
        mutual_information = None
    else:
        common = (ref_support & mov_support).reshape(-1)
        sample_positions = np.flatnonzero(common)[::64]
        sample_count = int(sample_positions.size)
        if sample_count < 2:
            correlation = None
            mutual_information = None
        else:
            x = reference.reshape(-1)[sample_positions].astype(np.float64, copy=False)
            y = moving.reshape(-1)[sample_positions].astype(np.float64, copy=False)
            x = np.clip((x - np.percentile(x, 1)) / max(float(np.percentile(x, 99) - np.percentile(x, 1)), np.finfo(float).eps), 0, 1)
            y = np.clip((y - np.percentile(y, 1)) / max(float(np.percentile(y, 99) - np.percentile(y, 1)), np.finfo(float).eps), 0, 1)
            sx, sy = float(np.std(x)), float(np.std(y))
            correlation = float(np.corrcoef(x, y)[0, 1]) if sx > 0 and sy > 0 else None
            histogram, _, _ = np.histogram2d(x, y, bins=64, range=((0, 1), (0, 1)))
            joint = histogram / max(float(histogram.sum()), 1.0)
            px, py = joint.sum(axis=1), joint.sum(axis=0)
            expected = px[:, None] * py[None, :]
            nonzero = joint > 0
            mutual_information = float(np.sum(joint[nonzero] * np.log(joint[nonzero] / expected[nonzero])))
    return {
        "reference_nonzero_voxels": int(ref_support.sum()),
        "moving_nonzero_voxels": int(mov_support.sum()),
        "nonzero_support_dice": float(dice) if dice is not None else None,
        "shared_nonzero_voxels_sampled": sample_count,
        "intensity_sampling_rule": "flat-index order within shared nonzero voxels, then every 64th; deterministic, not spatially uniform",
        "shared_nonzero_intensity_pearson": correlation,
        "shared_nonzero_intensity_mutual_information_nats": mutual_information,
        "interpretation": "descriptive only; no acceptance threshold or registration pass",
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    crosswalk = args.crosswalk.expanduser().resolve()
    output = args.output_dir.expanduser().resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite output directory: {output}")
    candidates, crosswalk_hash = review.read_crosswalk(crosswalk)
    selected = {str(int(patient_id)) for patient_id in args.patient_ids}
    by_id = {row["patient_id"]: row for row in candidates}
    missing = sorted(selected - set(by_id))
    if missing:
        raise ValueError(f"Patient IDs not in fixed candidate window: {missing}")
    rows = [by_id[patient_id] for patient_id in sorted(selected, key=int)]
    output.mkdir(parents=True, exist_ok=False)
    old_panel_size = review.PANEL_SIZE
    review.PANEL_SIZE = args.panel_size
    results = []
    try:
        for row in rows:
            arrays, headers, paths = review.read_patient_images(row)
            metrics = {
                "t0_to_t1": descriptive_pair_metrics(arrays["t0"], arrays["t1"]),
                "t0_to_t2": descriptive_pair_metrics(arrays["t0"], arrays["t2"]),
            }
            normalized, display_window = review.normalize_shared(arrays)
            patient_dir = output / f"patient_{row['patient_id']}"
            montages = [
                review.make_plane_montage(
                    row["patient_id"], plane_name, axis, normalized, headers,
                    display_window,
                    patient_dir / f"patient_{row['patient_id']}_{plane_name}_5mm_{args.panel_size}px.png",
                )
                for plane_name, axis in review.PLANES
            ]
            results.append({
                "patient_id": row["patient_id"],
                "gap_days_from_crosswalk": float(row["gap_days"]),
                "inputs": [
                    {"visit": visit, "filename": paths[visit].name,
                     "sha256": review.sha256_file(paths[visit]),
                     "shape_canonical": headers[visit]["shape"],
                     "spacing_mm_canonical": headers[visit]["spacing_mm"],
                     "spatial_units": headers[visit]["spatial_units"]}
                    for visit in review.VISITS
                ],
                "display_panel_size_px": args.panel_size,
                "descriptive_registration_triage": metrics,
                "montages": montages,
            })
            del arrays, normalized
    finally:
        review.PANEL_SIZE = old_panel_size
    result = {
        "study": "MRI-only high-resolution display and descriptive triage; no target scoring",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": Path(__file__).name,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_renderer_sha256": review.sha256_file(Path(review.__file__).resolve()),
        "crosswalk": {"filename": crosswalk.name, "sha256": crosswalk_hash},
        "selected_patient_ids": [row["patient_id"] for row in rows],
        "crosswalk_handling": "The crosswalk CSV is parsed and contains GTV-path metadata; GTV paths are not dereferenced or opened. Only the three T1-Gd NIfTI volumes per selected patient are opened and used.",
        "input_roles_read": ["t0_t1gd", "t1_t1gd", "t2_t1gd"],
        "forbidden_inputs_opened": {"gtv": False, "brain_mask": False},
        "metric_limits": [
            "nonzero-support Dice measures only nonzero image-intensity support; it is not anatomy, a brain mask, or evidence of valid negative labels",
            "intensity metrics can change with anatomy, treatment, acquisition, and normalization",
            "intensity metrics sample every 64th voxel in flattened shared-support order; this deterministic sample is not spatially uniform",
            "both-empty support makes Dice undefined; fewer than two sampled points makes correlation and mutual information undefined",
            "metrics have no pass threshold and do not certify full-volume or tumor-boundary alignment",
        ],
        "results": results,
    }
    with (output / "manifest.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--crosswalk", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--patient-ids", nargs="+", required=True)
    parser.add_argument("--panel-size", type=int, default=160)
    args = parser.parse_args()
    if args.panel_size < 112 or args.panel_size > 1024:
        raise ValueError("panel-size must be between 112 and 1024 pixels")
    result = run(args)
    print(json.dumps({"patients": len(result["results"]), "montages": sum(len(x["montages"]) for x in result["results"]), "panel_size_px": args.panel_size}, indent=2))


if __name__ == "__main__":
    main()
