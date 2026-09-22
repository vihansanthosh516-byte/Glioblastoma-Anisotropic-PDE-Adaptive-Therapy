#!/usr/bin/env python3
"""
Load MU-Glioma-Post cohort and compute tumor volumes over time.

Reads:  data/tcia/MU-Glioma-Post/PatientID_XXXX/Timepoint_N/*tumorMask.nii.gz
Writes: output/mu_glioma_volumes.csv
        output/mu_glioma_summary.json
"""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import nibabel as nib
import pandas as pd

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
DATA_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_patient_timepoints(patient_dir: Path) -> dict:
    """Load all timepoints for one patient, return volumes + voxel counts."""
    result = {"patient_id": patient_dir.name, "timepoints": {}}

    timepoint_dirs = sorted([t for t in patient_dir.iterdir() if t.is_dir()])
    for tp_dir in timepoint_dirs:
        tp_name = tp_dir.name  # e.g., Timepoint_1
        mask_files = list(tp_dir.glob("*tumorMask*.nii.gz"))
        if not mask_files:
            continue

        try:
            img = nib.load(mask_files[0])
            data = img.get_fdata()
            header = img.header

            # Get voxel dimensions in mm
            zooms = header.get_zooms()[:3]
            voxel_volume_mm3 = float(np.prod(zooms))

            # Count tumor voxels (any positive label)
            n_voxels = int((data > 0).sum())
            volume_mm3 = n_voxels * voxel_volume_mm3

            result["timepoints"][tp_name] = {
                "n_voxels": n_voxels,
                "voxel_volume_mm3": voxel_volume_mm3,
                "volume_mm3": volume_mm3,
                "shape": list(data.shape),
            }
        except Exception as e:
            result["timepoints"][tp_name] = {"error": str(e)}

    return result


def main():
    print("=" * 60)
    print("MU-GLIOMA-POST COHORT LOADING")
    print("=" * 60)
    print()

    patient_dirs = sorted([p for p in DATA_DIR.iterdir()
                           if p.is_dir() and p.name.startswith("PatientID_")])
    print(f"[LOAD] Found {len(patient_dirs)} patient directories")

    rows = []
    full_data = {}

    for p in patient_dirs:
        print(f"  Loading {p.name}...", end=" ")
        result = load_patient_timepoints(p)
        full_data[p.name] = result

        tp_names = sorted(result["timepoints"].keys())
        print(f"{len(tp_names)} timepoints: {tp_names}")

        # Build a row with volumes
        row = {"patient_id": p.name}
        for tp_name, tp_data in result["timepoints"].items():
            if "volume_mm3" in tp_data:
                row[tp_name] = tp_data["volume_mm3"]
                row[f"{tp_name}_voxels"] = tp_data["n_voxels"]
        rows.append(row)

    # Write CSV
    df = pd.DataFrame(rows)
    csv_path = OUTPUT_DIR / "mu_glioma_volumes.csv"
    df.to_csv(csv_path, index=False)
    print(f"\n[SAVE] {csv_path}")
    print(f"  Columns: {list(df.columns)}")
    print(f"  Rows: {len(df)}")

    # Summary stats
    tp_cols = [c for c in df.columns if c.startswith("Timepoint_") and not c.endswith("_voxels")]
    print(f"\n[STATS] Volume statistics per timepoint (mm³):")
    for tp in tp_cols:
        vals = df[tp].dropna()
        print(f"  {tp}: n={len(vals)}, "
              f"median={vals.median():.0f}, "
              f"mean={vals.mean():.0f}, "
              f"range=[{vals.min():.0f}, {vals.max():.0f}]")

    # Trajectory analysis: how many patients have decreasing volumes?
    if len(tp_cols) >= 2:
        decreasing = 0
        increasing = 0
        for _, row in df.iterrows():
            vals = [row[tp] for tp in tp_cols if pd.notna(row[tp])]
            if len(vals) >= 2:
                if vals[-1] < vals[0]:
                    decreasing += 1
                else:
                    increasing += 1
        print(f"\n[TRAJECTORY] Patients with shrinking tumor: {decreasing}")
        print(f"[TRAJECTORY] Patients with growing tumor: {increasing}")

    # Save summary
    summary = {
        "n_patients": len(df),
        "timepoints": tp_cols,
        "volume_stats": {
            tp: {
                "n": int(df[tp].notna().sum()),
                "median_mm3": float(df[tp].median()) if df[tp].notna().any() else None,
                "mean_mm3": float(df[tp].mean()) if df[tp].notna().any() else None,
                "min_mm3": float(df[tp].min()) if df[tp].notna().any() else None,
                "max_mm3": float(df[tp].max()) if df[tp].notna().any() else None,
            }
            for tp in tp_cols
        },
    }
    summary_path = OUTPUT_DIR / "mu_glioma_summary.json"
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n[SAVE] {summary_path}")

    print()
    print("[SUCCESS] MU-Glioma data loaded")


if __name__ == "__main__":
    main()