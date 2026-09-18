#!/usr/bin/env python3
"""
BraTS / NIfTI / DICOM Loader for GBM Digital Twin
==================================================
Loads real 3D brain tumor imaging data and prepares simulation-ready arrays.
Supports BraTS, TCGA-GBM, and standard neuroimaging formats.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import nibabel as nib
import numpy as np
from scipy.ndimage import zoom

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
BRATS_LABELS = {
    "background": 0,
    "necrotic_core": 1,      # NCR/NET - label 1
    "edema": 2,              # ED - label 2
    "enhancing_core": 3,     # ET - label 3
}

# BraTS 2021/2023 file naming convention
BRATS_MODALITIES = {
    "t1": "_t1.nii.gz",
    "t1c": "_t1ce.nii.gz",   # post-contrast T1
    "t2": "_t2.nii.gz",
    "flair": "_flair.nii.gz",
    "seg": "_seg.nii.gz",
}

DEFAULT_TARGET_SHAPE = (64, 64, 64)
DEFAULT_VOXEL_SIZE_MM = 1.0


# --------------------------------------------------------------------------- #
# Core Loading Functions
# --------------------------------------------------------------------------- #
def load_nifti(path: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Load NIfTI file and return data, affine, and metadata.
    
    Returns:
        data: numpy array (H, W, D) - raw voxel data
        affine: 4x4 affine transformation matrix (voxel -> world mm)
        meta: dict with shape, zooms, dtype
    """
    path = Path(path)
    img = nib.load(str(path))
    data = img.get_fdata(dtype=np.float32)
    affine = img.affine.copy()
    
    meta = {
        "shape": data.shape,
        "zooms": img.header.get_zooms()[:3],  # (dx, dy, dz) in mm
        "dtype": str(data.dtype),
        "affine": affine,
    }
    return data, affine, meta


def load_brats_patient(
    patient_dir: Union[str, Path],
    target_shape: Tuple[int, int, int] = DEFAULT_TARGET_SHAPE,
) -> Dict:
    """
    Load all BraTS modalities for a single patient.
    
    Args:
        patient_dir: Directory containing BraTS files (e.g., Brats2021_00001/)
        target_shape: Target (D, H, W) for resampling
    
    Returns:
        Dict with keys: t1, t1c, t2, flair, seg, affine, meta
    """
    patient_dir = Path(patient_dir)
    
    # Find patient ID from directory name
    patient_id = patient_dir.name
    
    results = {"patient_id": patient_id}
    
    # Load each modality
    for modality, suffix in BRATS_MODALITIES.items():
        # Try multiple naming patterns
        patterns = [
            f"{patient_id}{suffix}",
            f"*{suffix}",
        ]
        
        file_path = None
        for pattern in patterns:
            matches = list(patient_dir.glob(pattern))
            if matches:
                file_path = matches[0]
                break
        
        if file_path and file_path.exists():
            data, affine, meta = load_nifti(file_path)
            results[modality] = data
            results[f"{modality}_affine"] = affine
            results[f"{modality}_meta"] = meta
        else:
            results[modality] = None
            print(f"  [Warning] {modality} not found for {patient_id}")
    
    # Use segmentation affine as reference
    if results.get("seg_affine") is not None:
        results["reference_affine"] = results["seg_affine"]
    elif results.get("t1c_affine") is not None:
        results["reference_affine"] = results["t1c_affine"]
    else:
        results["reference_affine"] = np.eye(4)
    
    return results


def extract_tumor_masks(seg_data: np.ndarray) -> Dict[str, np.ndarray]:
    """
    Extract tumor subregions from BraTS segmentation.
    
    Args:
        seg_data: Segmentation array with labels 0, 1, 2, 3
    
    Returns:
        Dict with binary masks for each region
    """
    return {
        "whole_tumor": (seg_data > 0).astype(np.float32),           # 1+2+3
        "tumor_core": np.isin(seg_data, [1, 3]).astype(np.float32), # 1+3 (NCR+ET)
        "enhancing_tumor": (seg_data == 3).astype(np.float32),       # 3 (ET)
        "edema": (seg_data == 2).astype(np.float32),                 # 2 (ED)
        "necrotic_core": (seg_data == 1).astype(np.float32),         # 1 (NCR)
    }


def resample_to_grid(
    volume: np.ndarray,
    source_affine: np.ndarray,
    target_shape: Tuple[int, int, int] = DEFAULT_TARGET_SHAPE,
    target_voxel_size: float = DEFAULT_VOXEL_SIZE_MM,
    order: int = 1,  # 1=linear, 0=nearest (for masks)
    mode: str = "constant",
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Resample volume to target grid with isotropic voxels.
    
    Args:
        volume: Input volume (H, W, D) or (D, H, W)
        source_affine: 4x4 affine from source (voxel -> world mm)
        target_shape: Target (D, H, W)
        target_voxel_size: Target isotropic voxel size in mm
        order: Interpolation order
        mode: Padding mode
    
    Returns:
        resampled: (D, H, W) on target grid
        target_affine: 4x4 affine for target grid
    """
    # Ensure (D, H, W) ordering for 3D solver
    if volume.shape[-1] < volume.shape[0]:  # Likely (H, W, D)
        volume = np.transpose(volume, (2, 0, 1))
    
    D_src, H_src, W_src = volume.shape
    D_tgt, H_tgt, W_tgt = target_shape
    
    # Get source voxel sizes from affine
    src_zooms = np.abs([source_affine[0,0], source_affine[1,1], source_affine[2,2]])
    
    # Compute zoom factors
    zoom_factors = (
        src_zooms[2] / target_voxel_size * D_src / D_tgt,  # D (slice)
        src_zooms[1] / target_voxel_size * H_src / H_tgt,  # H
        src_zooms[0] / target_voxel_size * W_src / W_tgt,  # W
    )
    
    # Resample
    resampled = zoom(volume, zoom_factors, order=order, mode=mode, prefilter=False)
    
    # Crop or pad to exact target shape
    resampled = _crop_or_pad(resampled, target_shape)
    
    # Build target affine (isotropic, centered)
    target_affine = np.array([
        [target_voxel_size, 0, 0, -target_voxel_size * W_tgt / 2],
        [0, target_voxel_size, 0, -target_voxel_size * H_tgt / 2],
        [0, 0, target_voxel_size, -target_voxel_size * D_tgt / 2],
        [0, 0, 0, 1]
    ], dtype=np.float32)
    
    return resampled.astype(np.float32), target_affine


def _crop_or_pad(arr: np.ndarray, target_shape: Tuple[int, int, int]) -> np.ndarray:
    """Center-crop or pad array to target shape."""
    result = np.zeros(target_shape, dtype=arr.dtype)
    
    slices_src = []
    slices_tgt = []
    for i, (src_len, tgt_len) in enumerate(zip(arr.shape, target_shape)):
        if src_len <= tgt_len:
            start_tgt = (tgt_len - src_len) // 2
            slices_src.append(slice(0, src_len))
            slices_tgt.append(slice(start_tgt, start_tgt + src_len))
        else:
            start_src = (src_len - tgt_len) // 2
            slices_src.append(slice(start_src, start_src + tgt_len))
            slices_tgt.append(slice(0, tgt_len))
    
    result[tuple(slices_tgt)] = arr[tuple(slices_src)]
    return result


def extract_dti_tensor_field(
    dwi_data: np.ndarray,
    dwi_affine: np.ndarray,
    bvals: np.ndarray,
    bvecs: np.ndarray,
    mask: Optional[np.ndarray] = None,
    target_shape: Tuple[int, int, int] = DEFAULT_TARGET_SHAPE,
) -> np.ndarray:
    """
    Extract diffusion tensor field from DWI data using DTIFIT-style fitting.
    
    Simplified implementation - for production use FSL dtifit or DIPY.
    
    Args:
        dwi_data: 4D DWI data (X, Y, Z, N_gradients)
        dwi_affine: Affine for DWI
        bvals: b-values array (N_gradients,)
        bvecs: b-vectors array (3, N_gradients) or (N_gradients, 3)
        mask: Optional brain mask (X, Y, Z)
        target_shape: Target grid shape (D, H, W)
    
    Returns:
        tensor_field: (3, 3, D, H, W) symmetric positive-definite tensors
    """
    # This is a placeholder - real DTI fitting requires FSL/DIPY
    # For now, return a synthetic tensor field matching target shape
    D, H, W = target_shape
    tensor_field = np.zeros((3, 3, D, H, W), dtype=np.float32)
    
    # Isotropic baseline
    tensor_field[0, 0] = 0.0013  # D_xx
    tensor_field[1, 1] = 0.0013  # D_yy
    tensor_field[2, 2] = 0.0013  # D_zz
    
    return tensor_field


def prepare_simulation_inputs(
    patient_data: Dict,
    target_shape: Tuple[int, int, int] = DEFAULT_TARGET_SHAPE,
    target_voxel_size: float = DEFAULT_VOXEL_SIZE_MM,
) -> Dict:
    """
    Prepare all simulation inputs from loaded patient data.
    
    Args:
        patient_data: Output from load_brats_patient()
        target_shape: Target grid shape (D, H, W)
        target_voxel_size: Target isotropic voxel size in mm
    
    Returns:
        Dict with simulation-ready arrays:
            u0: Initial tumor density (D, H, W) float32 [0,1]
            tensor_field: Diffusion tensor (3, 3, D, H, W) float32
            masks: Dict of tumor subregion masks
            affine: Target affine
    """
    # Get reference affine
    ref_affine = patient_data.get("reference_affine", np.eye(4))
    
    # 1. Initial tumor density from segmentation
    if patient_data["seg"] is not None:
        masks = extract_tumor_masks(patient_data["seg"])
        # Use whole tumor as initial condition
        u0 = resample_to_grid(
            masks["whole_tumor"],
            patient_data["seg_affine"],
            target_shape=target_shape,
            target_voxel_size=target_voxel_size,
            order=0,  # nearest for mask
        )[0]
    elif patient_data["t1c"] is not None:
        # Fallback: threshold T1ce for enhancing tumor
        t1c = patient_data["t1c"]
        t1c_norm = (t1c - t1c.min()) / (t1c.max() - t1c.min() + 1e-8)
        u0 = resample_to_grid(
            (t1c_norm > 0.5).astype(np.float32),
            patient_data["t1c_affine"],
            target_shape=target_shape,
            target_voxel_size=target_voxel_size,
            order=0,
        )[0]
    else:
        # Synthetic spherical seed
        D, H, W = target_shape
        z, y, x = np.mgrid[0:D, 0:H, 0:W]
        center = (D//2, H//2, W//2)
        dist = np.sqrt((x-center[2])**2 + (y-center[1])**2 + (z-center[0])**2)
        u0 = (dist < 3).astype(np.float32)
    
    # 2. Tensor field - placeholder (real DTI needs FSL/DIPY)
    # For now, create synthetic anisotropic tensor field
    tensor_field = _create_synthetic_tensor_field(target_shape)
    
    # 3. Resample all masks
    resampled_masks = {}
    if patient_data["seg"] is not None:
        masks = extract_tumor_masks(patient_data["seg"])
        for name, mask in masks.items():
            resampled_masks[name] = resample_to_grid(
                mask, patient_data["seg_affine"],
                target_shape=target_shape,
                target_voxel_size=target_voxel_size,
                order=0,
            )[0]
    
    # Build target affine
    target_affine = np.array([
        [target_voxel_size, 0, 0, -target_voxel_size * target_shape[2] / 2],
        [0, target_voxel_size, 0, -target_voxel_size * target_shape[1] / 2],
        [0, 0, target_voxel_size, -target_voxel_size * target_shape[0] / 2],
        [0, 0, 0, 1]
    ], dtype=np.float32)
    
    return {
        "u0": u0,
        "tensor_field": tensor_field,
        "masks": resampled_masks,
        "affine": target_affine,
        "patient_id": patient_data.get("patient_id", "UNKNOWN"),
    }


def _create_synthetic_tensor_field(
    target_shape: Tuple[int, int, int],
    tract_orientation: str = "corpus_callosum",
) -> np.ndarray:
    """
    Create synthetic 3D tensor field with white matter tract anisotropy.
    
    This is a placeholder for real DTI-derived tensors.
    Real implementation would use FSL dtifit or DIPY.
    """
    D, H, W = target_shape
    tensor_field = np.zeros((3, 3, D, H, W), dtype=np.float32)
    
    # Base isotropic diffusivity (gray matter)
    D_iso = 0.0013  # mm²/day
    
    # Tract parameters
    if tract_orientation == "corpus_callosum":
        # Corpus callosum: left-right (x-axis) at mid-sagittal
        tract_center_z = D // 2
        tract_center_y = H // 2
        tract_width = H // 4
        direction = np.array([1.0, 0.0, 0.0])  # x-axis
    elif tract_orientation == "cingulum":
        # Cingulum: anterior-posterior (y-axis) near midline
        tract_center_z = D // 2
        tract_center_y = H // 3
        tract_width = H // 5
        direction = np.array([0.0, 1.0, 0.0])  # y-axis
    else:
        # Diagonal
        tract_center_z = D // 2
        tract_center_y = H // 2
        tract_width = H // 4
        direction = np.array([1.0, 1.0, 0.0]) / np.sqrt(2)
    
    D_parallel = 0.013   # mm²/day along tract
    D_perp = 0.0013      # mm²/day perpendicular
    
    # Build tensor field
    z, y, x = np.mgrid[0:D, 0:H, 0:W]
    center_z, center_y = tract_center_z, tract_center_y
    
    # Distance to tract centerline
    pos = np.stack([x - W/2, y - center_y, z - center_z], axis=-1)
    proj = np.sum(pos * direction, axis=-1, keepdims=True) * direction
    dist_perp = np.sqrt(np.sum((pos - proj)**2, axis=-1))
    
    in_tract = dist_perp < tract_width
    
    # Tensor components
    n = direction
    delta_D = D_parallel - D_perp
    
    tensor_field[0, 0] = D_perp + delta_D * n[0]**2
    tensor_field[1, 1] = D_perp + delta_D * n[1]**2
    tensor_field[2, 2] = D_perp + delta_D * n[2]**2
    tensor_field[0, 1] = tensor_field[1, 0] = delta_D * n[0] * n[1]
    tensor_field[0, 2] = tensor_field[2, 0] = delta_D * n[0] * n[2]
    tensor_field[1, 2] = tensor_field[2, 1] = delta_D * n[1] * n[2]
    
    # Outside tract: isotropic
    tensor_field[0, 0][~in_tract] = D_iso
    tensor_field[1, 1][~in_tract] = D_iso
    tensor_field[2, 2][~in_tract] = D_iso
    tensor_field[0, 1][~in_tract] = 0
    tensor_field[0, 2][~in_tract] = 0
    tensor_field[1, 2][~in_tract] = 0
    
    return tensor_field


# --------------------------------------------------------------------------- #
# CLI Entry Point
# --------------------------------------------------------------------------- #
def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Load BraTS patient and prepare simulation inputs"
    )
    parser.add_argument("patient_dir", type=str, help="Path to BraTS patient directory")
    parser.add_argument("--target-shape", type=int, nargs=3, default=DEFAULT_TARGET_SHAPE,
                        help="Target grid shape (D H W)")
    parser.add_argument("--voxel-size", type=float, default=DEFAULT_VOXEL_SIZE_MM,
                        help="Target isotropic voxel size (mm)")
    parser.add_argument("--output", type=str, default=None,
                        help="Output NPZ file for simulation inputs")
    parser.add_argument("--plot", action="store_true",
                        help="Save visualization plots")
    
    args = parser.parse_args()
    
    print(f"Loading patient from {args.patient_dir}...")
    patient_data = load_brats_patient(args.patient_dir)
    
    print(f"Preparing simulation inputs...")
    sim_inputs = prepare_simulation_inputs(
        patient_data,
        target_shape=tuple(args.target_shape),
        target_voxel_size=args.voxel_size,
    )
    
    if args.output:
        output_path = Path(args.output)
        np.savez_compressed(output_path, **sim_inputs)
        print(f"Saved simulation inputs to {output_path}")
    
    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        
        # Plot mid-slices
        D, H, W = sim_inputs["u0"].shape
        mid = D // 2
        
        fig, axes = plt.subplots(2, 3, figsize=(12, 8))
        
        axes[0, 0].imshow(sim_inputs["u0"][mid], cmap='hot', origin='lower')
        axes[0, 0].set_title('Initial Tumor (u0)')
        axes[0, 0].axis('off')
        
        axes[0, 1].imshow(sim_inputs["tensor_field"][0, 0, mid], cmap='viridis', origin='lower')
        axes[0, 1].set_title('D_xx (mid-slice)')
        axes[0, 1].axis('off')
        
        axes[0, 2].imshow(sim_inputs["tensor_field"][1, 1, mid], cmap='viridis', origin='lower')
        axes[0, 2].set_title('D_yy (mid-slice)')
        axes[0, 2].axis('off')
        
        for name, mask in sim_inputs["masks"].items():
            axes[1, 0].imshow(mask[mid], cmap='Blues', alpha=0.5, origin='lower')
        axes[1, 0].set_title('Tumor Masks')
        axes[1, 0].axis('off')
        
        axes[1, 1].axis('off')
        axes[1, 2].axis('off')
        
        plt.tight_layout()
        plot_path = Path(args.output or "sim_inputs").with_suffix(".png")
        plt.savefig(plot_path, dpi=200, bbox_inches="tight")
        print(f"Saved plot to {plot_path}")
        plt.close()


if __name__ == "__main__":
    main()