#!/usr/bin/env python3
"""Spatial Validation Metrics: Dice Similarity Coefficient & Hausdorff Distance.

Adds rigorous 3D/2D spatial validation beyond 1D volume scalars using the
medical imaging standard metrics:

  - Dice Similarity Coefficient (DSC): overlap agreement in [0, 1]
  - Hausdorff Distance (HD) / HD95: boundary distance in mm
  - Mean Surface Distance (MSD): average boundary-to-boundary distance
  - Normalized Volume Difference: relative volume error

These metrics quantify spatial agreement between a SIMULATED tumor mask and
an OBSERVED (ground-truth / reference) tumor mask.

Grid spacing DX=1.0 mm is assumed throughout, so all distance outputs are in mm.

Usage:
    python src/53_spatial_metrics.py --validate
"""
from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.ndimage import binary_erosion, distance_transform_edt
from scipy.spatial.distance import directed_hausdorff

warnings.filterwarnings("ignore")

# Grid spacing (mm) — must match src/42_anisotropic_pde.py / src/45_validation_synthesis.py
DX_MM = 1.0
DX = DX_MM

# Default threshold for converting continuous density to a binary mask
DEFAULT_THRESHOLD = 0.1

# Clinical thresholds
DSC_CLINICAL_THRESHOLD = 0.7   # DSC >= 0.7 indicates good spatial agreement
HD_CLINICAL_THRESHOLD_MM = 5.0  # HD <= 5mm indicates good boundary alignment


def dice_similarity_coefficient(
    mask_a: np.ndarray,
    mask_b: np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
) -> float:
    """Compute the Dice Similarity Coefficient (DSC) between two tumor masks.

        DSC = 2|A intersection B| / (|A| + |B|)

    Returns a value in [0, 1], where 1 = perfect spatial overlap.
    A value >= 0.7 indicates good spatial agreement (clinical threshold).

    Args:
        mask_a: First continuous or binary tumor mask.
        mask_b: Second continuous or binary tumor mask.
        threshold: Intensity threshold to binarize continuous masks.

    Returns:
        DSC in [0, 1]; 1.0 if both masks are empty (perfect agreement).
    """
    binary_a = (mask_a > threshold).astype(bool)
    binary_b = (mask_b > threshold).astype(bool)

    intersection = np.logical_and(binary_a, binary_b).sum()
    union = binary_a.sum() + binary_b.sum()

    if union == 0:
        return 1.0  # Both empty = perfect agreement
    return float(2.0 * intersection / union)


def extract_boundary(
    mask: np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
) -> np.ndarray:
    """Extract surface voxels from a binary mask using morphological erosion.

    The boundary is the set of voxels inside the mask that are adjacent to a
    background voxel (i.e., the surface/shell).

    Args:
        mask: Continuous or binary tumor mask.
        threshold: Intensity threshold to binarize continuous masks.

    Returns:
        Array of shape (N, ndim) with the coordinates of boundary voxels.
        Returns an empty array (shape (0, ndim)) if the mask is empty.
    """
    binary = mask > threshold
    if not binary.any():
        return np.empty((0, binary.ndim), dtype=int)

    eroded = binary_erosion(binary)
    boundary = binary & ~eroded
    return np.column_stack(np.where(boundary))


def hausdorff_distance(
    mask_a: np.ndarray,
    mask_b: np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    percentile: Optional[float] = None,
) -> float:
    """Compute the (symmetric) Hausdorff Distance between two tumor masks.

        HD = max{ sup inf ||a-b||, sup inf ||b-a|| }

    If `percentile` is specified (e.g., 95), returns the robust HD95 instead of
    the maximum distance.

    Units: mm (based on grid spacing DX_MM=1.0 mm).
    A value <= 5mm indicates good boundary alignment (clinical threshold).

    Args:
        mask_a: First tumor mask.
        mask_b: Second tumor mask.
        threshold: Intensity threshold to binarize continuous masks.
        percentile: If given (0-100), compute the percentile-based robust HD
            (e.g., 95 for HD95) instead of the maximum distance.

    Returns:
        Hausdorff distance in mm (0.0 if either mask is empty).
    """
    boundary_a = extract_boundary(mask_a, threshold)
    boundary_b = extract_boundary(mask_b, threshold)

    if len(boundary_a) == 0 or len(boundary_b) == 0:
        return 0.0

    def _distances(src: np.ndarray, dst: np.ndarray) -> np.ndarray:
        """All nearest-neighbor distances from src to dst (vectorized)."""
        # Use distance transform for efficiency: build a grid marking dst voxels,
        # compute the EDT from any non-dst voxel, then sample at src.
        shape = mask_b.shape
        dst_grid = np.zeros(shape, dtype=bool)
        # Clamp indices to valid range
        idxs = np.clip(dst, 0, np.array(shape) - 1)
        dst_grid[tuple(idxs.T)] = True
        # edt computes distance from True voxels; sample at src locations
        ed = distance_transform_edt(~dst_grid, sampling=DX_MM)
        src_idx = np.clip(src, 0, np.array(shape) - 1)
        return ed[tuple(src_idx.T)]

    d_ab = _distances(boundary_a, boundary_b)  # a -> b nearest
    d_ba = _distances(boundary_b, boundary_a)  # b -> a nearest

    if percentile is not None:
        # Robust: use requested percentile instead of the maximum
        all_d = np.concatenate([d_ab, d_ba])
        return float(np.percentile(all_d, percentile))

    return float(max(d_ab.max(), d_ba.max()))


def mean_surface_distance(
    mask_a: np.ndarray,
    mask_b: np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
) -> float:
    """Compute the Mean Surface Distance (MSD) between two tumor masks.

    MSD = mean over boundaries of (mean distance from a-boundary to b-boundary
          and mean distance from b-boundary to a-boundary).

    Units: mm.

    Args:
        mask_a: First tumor mask.
        mask_b: Second tumor mask.
        threshold: Intensity threshold to binarize continuous masks.

    Returns:
        Mean surface distance in mm (0.0 if either mask is empty).
    """
    boundary_a = extract_boundary(mask_a, threshold)
    boundary_b = extract_boundary(mask_b, threshold)

    if len(boundary_a) == 0 or len(boundary_b) == 0:
        return 0.0

    shape = mask_a.shape

    def _mean_dist(src: np.ndarray, dst: np.ndarray) -> float:
        dst_grid = np.zeros(shape, dtype=bool)
        idxs = np.clip(dst, 0, np.array(shape) - 1)
        dst_grid[tuple(idxs.T)] = True
        ed = distance_transform_edt(~dst_grid, sampling=DX_MM)
        src_idx = np.clip(src, 0, np.array(shape) - 1)
        return float(np.mean(ed[tuple(src_idx.T)]))

    msd_ab = _mean_dist(boundary_a, boundary_b)
    msd_ba = _mean_dist(boundary_b, boundary_a)
    return float(0.5 * (msd_ab + msd_ba))


def normalized_volume_difference(
    volume_a: float,
    volume_b: float,
) -> float:
    """Compute the normalized volume difference: |V_a - V_b| / max(V_a, V_b).

    Returns a value in [0, 1]; 0 = identical volumes.
    Returns 0.0 if both volumes are zero (perfect agreement).
    """
    denom = max(abs(volume_a), abs(volume_b))
    if denom == 0:
        return 0.0
    return float(abs(volume_a - volume_b) / denom)


def compute_spatial_metrics(
    simulated: np.ndarray,
    observed: np.ndarray,
    threshold: float = DEFAULT_THRESHOLD,
    compute_hd95: bool = True,
) -> Dict[str, Any]:
    """Compute the full spatial validation metric suite.

    Args:
        simulated: Simulated/predicted tumor mask (continuous or binary).
        observed: Observed/reference tumor mask (continuous or binary).
        threshold: Intensity threshold to binarize continuous masks.
        compute_hd95: If True, also compute the robust 95th-percentile HD.

    Returns:
        Dictionary with:
            - dice_coefficient: DSC in [0, 1]
            - hausdorff_distance_mm: max HD in mm
            - hausdorff_distance_95_mm: HD95 in mm (if compute_hd95)
            - mean_surface_distance_mm: MSD in mm
            - normalized_volume_difference: in [0, 1]
            - simulated_volume_mm3: voxel count * DX^3
            - observed_volume_mm3: voxel count * DX^3
            - meets_clinical_threshold: bool (DSC >= 0.7 AND HD <= 5mm)
    """
    dsc = dice_similarity_coefficient(simulated, observed, threshold)
    hd = hausdorff_distance(simulated, observed, threshold)
    msd = mean_surface_distance(simulated, observed, threshold)

    sim_vol = float(np.sum(simulated > threshold) * (DX_MM ** simulated.ndim))
    obs_vol = float(np.sum(observed > threshold) * (DX_MM ** observed.ndim))
    nvd = normalized_volume_difference(sim_vol, obs_vol)

    result = {
        "dice_coefficient": float(dsc),
        "hausdorff_distance_mm": float(hd),
        "mean_surface_distance_mm": float(msd),
        "normalized_volume_difference": float(nvd),
        "simulated_volume_mm3": float(sim_vol),
        "observed_volume_mm3": float(obs_vol),
        "meets_clinical_threshold": bool(
            dsc >= DSC_CLINICAL_THRESHOLD and hd <= HD_CLINICAL_THRESHOLD_MM
        ),
    }
    if compute_hd95:
        hd95 = hausdorff_distance(simulated, observed, threshold, percentile=95)
        result["hausdorff_distance_95_mm"] = float(hd95)
    return result


def run_self_validation() -> Dict[str, Any]:
    """Self-contained validation of spatial metrics on synthetic shapes.

    Verifies the metrics behave correctly on known geometries:
      - Identical masks -> DSC=1, HD=0
      - Disjoint masks  -> DSC=0, HD>0
      - Half-shifted masks -> DSC in (0, 1), HD = shift distance
    """
    print(f"\n{'='*70}")
    print("SPATIAL METRICS — SELF VALIDATION")
    print(f"{'='*70}\n")

    results: Dict[str, Any] = {}

    # Test 1: Identical masks
    shape = (50, 50)
    mask = np.zeros(shape)
    yy, xx = np.mgrid[0:50, 0:50]
    mask[(xx - 25) ** 2 + (yy - 25) ** 2 <= 10 ** 2] = 1.0
    m = compute_spatial_metrics(mask, mask)
    print(f"Identical masks:    DSC={m['dice_coefficient']:.3f}, "
          f"HD={m['hausdorff_distance_mm']:.2f}mm, MSD={m['mean_surface_distance_mm']:.2f}mm")
    results["identical"] = m
    assert abs(m["dice_coefficient"] - 1.0) < 1e-9, "Identical masks must give DSC=1"
    assert m["hausdorff_distance_mm"] < 1e-9, "Identical masks must give HD=0"

    # Test 2: Disjoint masks
    mask2 = np.zeros(shape)
    mask2[(xx - 2) ** 2 + (yy - 25) ** 2 <= 8 ** 2] = 1.0  # far left, no overlap
    m2 = compute_spatial_metrics(mask, mask2)
    print(f"Disjoint masks:     DSC={m2['dice_coefficient']:.3f}, "
          f"HD={m2['hausdorff_distance_mm']:.2f}mm, MSD={m2['mean_surface_distance_mm']:.2f}mm")
    results["disjoint"] = m2
    assert m2["dice_coefficient"] == 0.0, "Disjoint masks must give DSC=0"
    assert m2["hausdorff_distance_mm"] > 0, "Disjoint masks must give HD>0"

    # Test 3: Small boundary shift (HD equals shift distance)
    mask3 = np.zeros(shape)
    mask3[(xx - 27) ** 2 + (yy - 25) ** 2 <= 10 ** 2] = 1.0  # shifted 2px right
    m3 = compute_spatial_metrics(mask, mask3)
    print(f"Shifted +2mm:        DSC={m3['dice_coefficient']:.3f}, "
          f"HD={m3['hausdorff_distance_mm']:.2f}mm, MSD={m3['mean_surface_distance_mm']:.2f}mm")
    results["shifted_2mm"] = m3
    assert m3["dice_coefficient"] > 0.7, "Small shift should keep DSC >= 0.7"
    assert 1.0 <= m3["hausdorff_distance_mm"] <= 3.0, (
        "HD for a 2mm shift should be ~2mm"
    )

    # Test 4: Empty observed, non-empty simulated
    empty = np.zeros(shape)
    m4 = compute_spatial_metrics(mask, empty)
    print(f"Observed empty:     DSC={m4['dice_coefficient']:.3f}, HD={m4['hausdorff_distance_mm']:.2f}mm")
    results["empty_observed"] = m4
    assert m4["dice_coefficient"] == 0.0, "Empty observed must give DSC=0"
    assert m4["hausdorff_distance_mm"] == 0.0, "Empty observed must give HD=0"

    # Test 5: Both empty
    m5 = compute_spatial_metrics(empty, empty)
    print(f"Both empty:          DSC={m5['dice_coefficient']:.3f}, HD={m5['hausdorff_distance_mm']:.2f}mm")
    results["both_empty"] = m5
    assert m5["dice_coefficient"] == 1.0, "Both empty must give DSC=1"

    print(f"\n{'='*70}")
    print("ALL SELF-VALIDATION ASSERTIONS PASSED")
    print(f"{'='*70}\n")
    return results


def main():
    parser = argparse.ArgumentParser(
        description="Spatial validation metrics (Dice + Hausdorff)"
    )
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Run self-validation on synthetic shapes",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output JSON file for validation results",
    )
    args = parser.parse_args()

    if args.validate:
        results = run_self_validation()
        if args.output:
            out = Path(args.output)
            out.parent.mkdir(parents=True, exist_ok=True)
            with open(out, "w") as f:
                json.dump(results, f, indent=2)
            print(f"Results saved to: {out}")
        return 0

    parser.print_help()
    print("\nExamples:")
    print("  python src/53_spatial_metrics.py --validate")
    return 0


if __name__ == "__main__":
    sys.exit(main())
