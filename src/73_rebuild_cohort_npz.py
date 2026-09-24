#!/usr/bin/env python3
"""
Rebuild output/spatial_recurrence_profiles.npz using REAL tumor masks as
initial conditions from MU-Glioma-Post Timepoint_1 imaging.

Unlike the previous version (Gaussian init for all patients), this version:
  - Loads each patient's real tumor mask from data/tcia/MU-Glioma-Post/
  - Resamples 3D mask to a 2D 100x100 grid (mid-axial slice)
  - Uses the real mask as the initial condition for the PDE

Result: per-patient geometry varies (D_f, P/A, spatial extent).

Reads:  output/mu_glioma_params_real.csv
        data/tcia/MU-Glioma-Post/PatientID_XXXX/Timepoint_*/..._tumorMask.nii.gz
Writes: output/spatial_recurrence_profiles.npz (BACKED UP first)
"""

from __future__ import annotations

import shutil
from pathlib import Path
import numpy as np
import pandas as pd
import nibabel as nib
from scipy.ndimage import zoom, gaussian_filter

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
DATA_DIR = PROJECT_ROOT / "data" / "tcia" / "MU-Glioma-Post"

GRID_SIZE = 100
TARGET_GENES = ["S100A6", "S100A11", "S100A8", "CCL3L1"]
ZONE_REGIONS = {
    "Leading Edge": [0, 33],
    "Cellular Tumor": [33, 66],
    "Infiltrating Tumor": [66, 100],
}


def find_first_mask(patient_id: str) -> Path | None:
    """Return path to the earliest available tumor mask for a patient."""
    pdir = DATA_DIR / patient_id
    if not pdir.exists():
        return None

    # Find all timepoints, take the lowest-numbered one
    tps = sorted([p for p in pdir.iterdir() if p.is_dir() and p.name.startswith("Timepoint_")],
                 key=lambda p: int(p.name.split("_")[-1]))
    for tp in tps:
        masks = list(tp.glob("*_tumorMask.nii.gz"))
        if masks:
            return masks[0]
    return None


def extract_mid_slice(mask_path: Path, target_size: int = GRID_SIZE) -> np.ndarray | None:
    """
    Load 3D NIfTI mask, take the mid-axial slice, resample to target_size x target_size,
    return as a binary float array.
    """
    try:
        img = nib.load(str(mask_path))
        data = img.get_fdata()
    except Exception as e:
        print(f"    ERROR loading {mask_path}: {e}")
        return None

    if data.ndim != 3:
        return None

    # Take mid-axial slice
    mid = data.shape[2] // 2
    slice_2d = data[:, :, mid]

    # Resample to target grid using NEAREST NEIGHBOR (order=0)
    # order=1 (bilinear) would smear the tiny 2.6% tumor region into the whole grid
    zoom_factors = (target_size / slice_2d.shape[0], target_size / slice_2d.shape[1])
    resampled = zoom(slice_2d, zoom_factors, order=0)  # order=0 = nearest neighbor

    # Binarize — mask has values 0/1/2/3, so > 0 catches all tumor labels
    resampled = (resampled > 0).astype(np.float64)

    # Guard against empty masks
    if resampled.sum() == 0:
        # Try other axial slices
        for mid in range(data.shape[2]):
            if data[:, :, mid].sum() > 0:
                slice_2d = data[:, :, mid]
                resampled = zoom(slice_2d,
                                (target_size / slice_2d.shape[0],
                                 target_size / slice_2d.shape[1]),
                                order=0)  # nearest neighbor
                resampled = (resampled > 0).astype(np.float64)
                break

    if resampled.sum() == 0:
        return None

    return resampled


def build_rho_field_from_mask(mask: np.ndarray, rho_per_day: float) -> np.ndarray:
    """
    Build rho field using the real mask as the spatial pattern.
    Tumor interior has high rho (from real growth rate), background has ~0.

    Uses abs(rho_per_day) so shrinking tumors (negative rho) still
    produce a positive metabolic field — the sign only affects growth
    direction, not the spatial distribution of activity.
    """
    from scipy.ndimage import distance_transform_edt

    # Distance from tumor boundary
    dist = distance_transform_edt(mask > 0.5)
    if dist.max() > 0:
        dist_norm = dist / dist.max()
    else:
        dist_norm = np.zeros_like(dist)

    # Use absolute value: tumor is metabolically active regardless of sign
    rho_abs = abs(rho_per_day)
    if rho_abs < 1e-6:
        rho_abs = 1e-6  # guard against zero

    # rho_field: high inside tumor (peak at center), ~0 in background
    mask_binary = (mask > 0.5).astype(float)
    rho_field = rho_abs * (0.3 + 0.7 * dist_norm) * mask_binary

    # Add tiny background value so the field is positive everywhere
    # but 1000x smaller than tumor interior
    rho_field = rho_field + (rho_abs * 1e-3) * (1 - mask_binary)

    # Soften edges
    rho_field = gaussian_filter(rho_field, sigma=2.0)

    return rho_field.astype(np.float64)


def build_D_field_from_mask(mask: np.ndarray, rho_per_day: float) -> np.ndarray:
    """
    Build D field using real mask. High D inside tumor (infiltrative),
    lower outside.
    """
    D_base = abs(rho_per_day) / 10.0
    if D_base == 0:
        D_base = 1e-5

    # Tumor interior: full diffusion; background: 10% of that
    D_field = np.where(mask > 0.5, D_base, D_base * 0.1)
    D_field = gaussian_filter(D_field, sigma=1.5)
    return D_field.astype(np.float64)


def main():
    print("=" * 60)
    print("REBUILDING NPZ WITH REAL TUMOR MASKS AS INITIAL CONDITIONS")
    print("=" * 60)
    print()

    npz_path = OUTPUT_DIR / "spatial_recurrence_profiles.npz"

    # Backup
    if npz_path.exists():
        backup = npz_path.with_suffix(".npz.bak_gaussian")
        shutil.copy(npz_path, backup)
        print(f"[BACKUP] {backup}")
    print()

    # Load real patients
    params = pd.read_csv(OUTPUT_DIR / "mu_glioma_params_real.csv")
    params = params[params["r_squared"] >= 0.5].copy()
    print(f"[LOAD] {len(params)} real patients with R^2 >= 0.5")

    # Build arrays
    patient_ids = []
    rho_fields = []
    D_fields = []
    density_maps = []
    risk_maps = []
    skipped = []

    for i, (_, row) in enumerate(params.iterrows()):
        pid = row["patient_id"]
        rho = float(row["rho_per_day"])

        mask_path = find_first_mask(pid)
        if mask_path is None:
            skipped.append((pid, "no mask"))
            continue

        mask = extract_mid_slice(mask_path, target_size=GRID_SIZE)
        if mask is None:
            skipped.append((pid, "empty mask"))
            continue

        rho_field = build_rho_field_from_mask(mask, rho)
        D_field = build_D_field_from_mask(mask, rho)
        density_map = mask / (mask.sum() + 1e-10) * 1000  # normalized density
        risk_map = rho_field / (rho_field.max() + 1e-10)

        patient_ids.append(pid)
        rho_fields.append(rho_field)
        D_fields.append(D_field)
        density_maps.append(density_map)
        risk_maps.append(risk_map)

        if (i + 1) % 25 == 0:
            print(f"  Processed {i + 1}/{len(params)} patients")

    n_kept = len(patient_ids)
    print()
    print(f"[BUILD] Kept {n_kept}/{len(params)} patients with usable masks")
    if skipped:
        print(f"[SKIP] {len(skipped)} patients skipped (first 5): {skipped[:5]}")

    # Stack
    rho_fields = np.stack(rho_fields)
    D_fields = np.stack(D_fields)
    density_maps = np.stack(density_maps)
    risk_maps = np.stack(risk_maps)
    patient_ids = np.array(patient_ids)

    # Save
    print()
    print(f"[SAVE] Writing {npz_path}")
    np.savez_compressed(
        npz_path,
        patient_ids=patient_ids,
        risk_maps=risk_maps,
        density_maps=density_maps,
        D_fields=D_fields,
        rho_fields=rho_fields,
        grid_size=np.array(GRID_SIZE),
        zone_regions=np.array([list(v) for v in ZONE_REGIONS.values()]),
        target_genes=np.array(TARGET_GENES),
        gene_weights=np.array([1.0, 1.0, 1.0, 1.0]),
    )

    # Verify
    check = np.load(npz_path, allow_pickle=True)
    print(f"[VERIFY] Keys: {list(check.keys())}")
    print(f"[VERIFY] n patients: {len(check['patient_ids'])}")
    print(f"[VERIFY] patient_ids[:5]: {list(check['patient_ids'][:5])}")

    # Per-patient geometry stats
    print()
    print("[GEOMETRY] Per-patient mask extent (first 5):")
    for i in range(min(5, n_kept)):
        mask_sum = (check["density_maps"][i] > 0).sum()
        nonzero_frac = mask_sum / (GRID_SIZE * GRID_SIZE)
        print(f"  {check['patient_ids'][i]}: tumor fraction = {nonzero_frac:.3f}")

    # Global stats
    all_fracs = [(check["density_maps"][i] > 0).sum() / (GRID_SIZE * GRID_SIZE)
                 for i in range(n_kept)]
    print()
    print(f"[STATS] Tumor fraction across {n_kept} patients:")
    print(f"  Mean:   {np.mean(all_fracs):.3f}")
    print(f"  Median: {np.median(all_fracs):.3f}")
    print(f"  Range:  [{np.min(all_fracs):.3f}, {np.max(all_fracs):.3f}]")

    print()
    print(f"[SUCCESS] Rebuilt npz with {n_kept} real patients using real tumor masks")


if __name__ == "__main__":
    main()