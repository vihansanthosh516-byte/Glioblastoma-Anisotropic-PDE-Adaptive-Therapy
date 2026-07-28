#!/usr/bin/env python3
"""
Real 3D Medical Imaging Loader for GBM Digital Twin
====================================================
Loads NIfTI/DICOM brain scans and converts to simulation-ready arrays.
Supports BraTS, TCGA-GBM, and standard neuroimaging formats.
"""
from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import nibabel as nib
import numpy as np
from scipy import ndimage
from scipy.ndimage import zoom

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Constants & Type Aliases
# --------------------------------------------------------------------------- #
Array3D = np.ndarray  # (D, H, W) convention for 3D solver
TensorField = np.ndarray  # (3, 3, D, H, W) SPD tensor field

# BraTS label mapping
BRATS_LABELS = {
    "background": 0,
    "necrotic_core": 1,      # NCR/NET
    "edema": 2,              # ED
    "enhancing_core": 3,     # ET
}

# Target simulation grid
DEFAULT_TARGET_SHAPE = (64, 64, 64)  # (D, H, W) for Track C 3D solver
DEFAULT_VOXEL_SIZE_MM = 1.0  # mm isotropic


# --------------------------------------------------------------------------- #
# Core Loading Functions
# --------------------------------------------------------------------------- #
def load_nifti(path: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray, Dict]:
    """
    Load NIfTI file and return data array, affine, and metadata.
    
    Returns:
        data: numpy array (H, W, D) or (D, H, W) - raw voxel data
        affine: 4x4 affine transformation matrix (voxel -> world mm)
        meta: dict with shape, zooms, dtype, etc.
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


def load_nifti_seg(path: Union[str, Path], labels: List[int] = None) -> np.ndarray:
    """
    Load BraTS-style segmentation mask.
    
    Args:
        path: Path to segmentation NIfTI
        labels: Labels to include in mask (default: all tumor [1, 2, 3])
    
    Returns:
        tumor_mask: (D, H, W) float32 binary mask [0, 1] 
        (1 = any tumor: necrotic core + edema + enhancing)
    """
    if labels is None:
        labels = [1, 2, 3]  # All tumor subregions
    
    data, _, _ = load_nifti(path)
    # BraTS: combine specified labels into single binary mask
    tumor_mask = np.isin(data, labels).astype(np.float32)
    return _reorder_to_dhw(tumor_mask)


def load_nifti_subregion(
    path: Union[str, Path], 
    label: int
) -> np.ndarray:
    """
    Load specific subregion from BraTS segmentation.
    
    Args:
        path: Path to segmentation NIfTI
        label: Single label value (1=necrotic, 2=edema, 3=enhancing)
    
    Returns:
        mask: (D, H, W) float32 binary mask
    """
    data, _, _ = load_nifti(path)
    mask = (data == label).astype(np.float32)
    return _reorder_to_dhw(mask)


def load_nifti_t1c(path: Union[str, Path], threshold: float = 0.5) -> np.ndarray:
    """
    Load T1ce (post-contrast T1) and extract enhancing core mask.
    
    Args:
        path: Path to T1ce NIfTI
        threshold: Intensity threshold for enhancing region (normalized 0-1)
    
    Returns:
        enhancing_mask: (D, H, W) float32 binary mask
    """
    data, _, _ = load_nifti(path)
    # Normalize to 0-1
    data_norm = (data - data.min()) / (data.max() - data.min() + 1e-8)
    enhancing_mask = (data_norm > threshold).astype(np.float32)
    return _reorder_to_dhw(enhancing_mask)


def load_nifti_flair(path: Union[str, Path], threshold: float = 0.5) -> np.ndarray:
    """
    Load FLAIR and extract edema mask.
    
    Args:
        path: Path to FLAIR NIfTI
        threshold: Intensity threshold for edema region (normalized 0-1)
    
    Returns:
        edema_mask: (D, H, W) float32 binary mask
    """
    data, _, _ = load_nifti(path)
    data_norm = (data - data.min()) / (data.max() - data.min() + 1e-8)
    edema_mask = (data_norm > threshold).astype(np.float32)
    return _reorder_to_dhw(edema_mask)


def load_nifti_t1(path: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load T1-weighted (pre-contrast) NIfTI.
    
    Returns:
        volume: (D, H, W) float32
        affine: 4x4 affine
    """
    data, affine, _ = load_nifti(path)
    return _reorder_to_dhw(data), affine


def load_dicom_series(dicom_dir: Union[str, Path]) -> Tuple[np.ndarray, np.ndarray]:
    """
    Load DICOM series from directory and reconstruct 3D volume.
    
    Returns:
        volume: (D, H, W) float32
        affine: 4x4 voxel-to-world transformation
    """
    try:
        import pydicom
    except ImportError:
        raise ImportError("pydicom required for DICOM loading: pip install pydicom")
    
    dicom_dir = Path(dicom_dir)
    files = list(dicom_dir.glob("*.dcm"))
    if not files:
        files = list(dicom_dir.rglob("*.dcm"))
    
    slices = []
    for f in files:
        ds = pydicom.dcmread(str(f))
        if hasattr(ds, 'InstanceNumber'):
            slices.append(ds)
    
    # Sort by slice location
    slices.sort(key=lambda s: float(getattr(s, 'SliceLocation', s.InstanceNumber)))
    
    # Stack pixel arrays
    volume = np.stack([s.pixel_array.astype(np.float32) for s in slices], axis=0)
    
    # Build affine from DICOM tags (simplified)
    if len(slices) > 1:
        dz = abs(float(slices[1].SliceLocation) - float(slices[0].SliceLocation))
    else:
        dz = float(getattr(slices[0], 'SliceThickness', 1.0))
    
    dx = float(getattr(slices[0], 'PixelSpacing', [1.0, 1.0])[0])
    dy = float(getattr(slices[0], 'PixelSpacing', [1.0, 1.0])[1])
    
    affine = np.array([
        [dx, 0, 0, 0],
        [0, dy, 0, 0],
        [0, 0, dz, 0],
        [0, 0, 0, 1]
    ], dtype=np.float32)
    
    return _reorder_to_dhw(volume), affine


# --------------------------------------------------------------------------- #
# Reordering Utilities (DHW <-> HWD)
# --------------------------------------------------------------------------- #
def _reorder_to_dhw(arr: np.ndarray) -> np.ndarray:
    """Convert (H, W, D) -> (D, H, W) for 3D solver convention."""
    if arr.ndim != 3:
        return arr
    # Nibabel typically returns (H, W, D) - transpose to (D, H, W)
    return np.transpose(arr, (2, 0, 1))


def _reorder_to_hwd(arr: np.ndarray) -> np.ndarray:
    """Convert (D, H, W) -> (H, W, D) for nibabel saving."""
    if arr.ndim != 3:
        return arr
    return np.transpose(arr, (1, 2, 0))


# --------------------------------------------------------------------------- #
# Resampling & Alignment
# --------------------------------------------------------------------------- #
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
        volume: Input volume (D, H, W) or (H, W, D)
        source_affine: 4x4 affine from source (voxel -> world mm)
        target_shape: (D, H, W) target grid size
        target_voxel_size: Target isotropic voxel size in mm
        order: Interpolation order (1=linear, 0=nearest)
        mode: Padding mode
    
    Returns:
        resampled: (D, H, W) on target grid
        target_affine: 4x4 affine for target grid
    """
    # Ensure (D, H, W) ordering
    if volume.shape[-1] < volume.shape[0]:  # Likely (H, W, D)
        volume = _reorder_to_dhw(volume)
    
    D_src, H_src, W_src = volume.shape
    D_tgt, H_tgt, W_tgt = target_shape
    
    # Get source voxel sizes from affine
    src_zooms = np.abs([source_affine[0,0], source_affine[1,1], source_affine[2,2]])
    
    # Compute zoom factors
    zoom_factors = (
        src_zooms[2] / target_voxel_size * D_src / D_tgt,  # D (slice) dimension
        src_zooms[1] / target_voxel_size * H_src / H_tgt,  # H dimension
        src_zooms[0] / target_voxel_size * W_src / W_tgt,  # W dimension
    )
    
    # Resample
    resampled = zoom(volume, zoom_factors, order=order, mode=mode, prefilter=False)
    
    # Crop or pad to exact target shape
    resampled = _crop_or_pad(resampled, target_shape)
    
    # Build target affine (isotropic, centered at origin)
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
            # Pad
            start_tgt = (tgt_len - src_len) // 2
            slices_src.append(slice(0, src_len))
            slices_tgt.append(slice(start_tgt, start_tgt + src_len))
        else:
            # Crop
            start_src = (src_len - tgt_len) // 2
            slices_src.append(slice(start_src, start_src + tgt_len))
            slices_tgt.append(slice(0, tgt_len))
    
    result[tuple(slices_tgt)] = arr[tuple(slices_src)]
    return result


# --------------------------------------------------------------------------- #
# DTI Tensor Field Extraction (Simplified)
# --------------------------------------------------------------------------- #
def extract_dti_tensors_simple(
    dwi_path: Union[str, Path],
    bval_path: Union[str, Path],
    bvec_path: Union[str, Path],
    mask_path: Optional[Union[str, Path]] = None,
    target_shape: Tuple[int, int, int] = DEFAULT_TARGET_SHAPE,
    d_parallel: float = 0.013,
    d_perpendicular: float = 0.0013,
) -> TensorField:
    """
    Extract diffusion tensor field from DWI data.
    
    This is a SIMPLIFIED placeholder. For production use:
    - FSL dtifit (dtifit)
    - DIPY (dipy.reconst.dti.TensorModel)
    - MRTRIX3 (dwi2tensor)
    
    This function creates a synthetic tensor field aligned to the principal
    diffusion direction estimated from the b=0 and high b-value images.
    
    Args:
        dwi_path: 4D DWI NIfTI (X, Y, Z, N_gradients)
        bval_path: b-values text file
        bvec_path: b-vectors text file (3 x N_gradients or N_gradients x 3)
        mask_path: Optional brain mask NIfTI
        target_shape: Target grid (D, H, W)
        d_parallel: Diffusion along principal direction (mm^2/day)
        d_perpendicular: Diffusion perpendicular (mm^2/day)
    
    Returns:
        tensor_field: (3, 3, D, H, W) symmetric positive-definite tensors
    """
    # Load DWI
    dwi_data, dwi_affine, _ = load_nifti(dwi_path)
    dwi_data = _reorder_to_dhw(dwi_data)  # (D, H, W, N_grad)
    
    # Load bvals/bvecs
    bvals = np.loadtxt(bval_path)
    bvecs = np.loadtxt(bvec_path)
    if bvecs.shape[0] == 3:
        bvecs = bvecs.T  # Ensure (N_grad, 3)
    
    D, H, W, N_grad = dwi_data.shape
    
    # Find b=0 indices
    b0_indices = np.where(bvals < 50)[0]
    high_b_indices = np.where(bvals > 500)[0]
    
    if len(b0_indices) == 0 or len(high_b_indices) == 0:
        raise ValueError("Need both b=0 and high b-value volumes")
    
    # Average b=0 images
    S0 = np.mean(dwi_data[..., b0_indices], axis=-1)  # (D, H, W)
    
    # Estimate principal direction per voxel using simple ADC
    # ADC = -log(S/S0) / b  (for single high b-value)
    if len(high_b_indices) > 0:
        # Use first high b-value for direction estimation
        b_high = bvals[high_b_indices[0]]
        S_high = dwi_data[..., high_b_indices[0]]
        adc = -np.log(np.maximum(S_high / (S0 + 1e-8), 1e-8)) / b_high
        adc = np.clip(adc, 0, 0.01)
    else:
        adc = np.ones((D, H, W)) * 0.001
    
    # For simplicity, create synthetic tract-like principal directions
    # In practice, you'd fit full tensor model and get eigenvectors
    # Here we create a corpus callosum-like tract structure
    tensor_field = _create_synthetic_tract_tensor_field(
        shape=(D, H, W),
        d_parallel=d_parallel,
        d_perpendicular=d_perpendicular,
    )
    
    # Resample to target grid if needed
    if (D, H, W) != target_shape:
        # Resample each tensor component
        resampled = np.zeros((3, 3) + target_shape, dtype=np.float32)
        for i in range(3):
            for j in range(3):
                resampled[i, j], _ = resample_to_grid(
                    tensor_field[i, j], dwi_affine, target_shape, order=1
                )
        tensor_field = resampled
    
    return tensor_field


def _create_synthetic_tract_tensor_field(
    shape: Tuple[int, int, int],
    d_parallel: float = 0.013,
    d_perpendicular: float = 0.0013,
) -> TensorField:
    """
    Create synthetic corpus callosum-like tract tensor field.
    This is a placeholder for real DTI fitting.
    """
    D, H, W = shape
    tensor_field = np.zeros((3, 3, D, H, W), dtype=np.float32)
    
    # Create coordinate grids
    z, y, x = np.mgrid[0:D, 0:H, 0:W]
    center = np.array([D/2, H/2, W/2])
    
    # Corpus callosum: horizontal tract through mid-sagittal
    # Elliptical cross-section
    dist_cc = np.sqrt((y - center[1])**2 + (x - center[0])**2)
    cc_mask = dist_cc < min(H, W) / 4
    
    # Principal direction: left-right (x-axis) in corpus callosum
    n_cc = np.array([1.0, 0.0, 0.0])  # Left-right
    
    # Cingulum: curved tract
    # Approximate as curved tract in superior region
    cing_mask = (z > center[0] * 0.7) & (dist_cc < min(H, W) / 3)
    # Direction follows anterior-posterior curve
    n_cing = np.array([0.0, 1.0, 0.0])  # Anterior-posterior
    
    # Initialize with isotropic gray matter
    for i in range(3):
        for j in range(3):
            if i == j:
                tensor_field[i, j] = d_perpendicular
            else:
                tensor_field[i, j] = 0.0
    
    # Corpus callosum: anisotropic along x
    if np.any(cc_mask):
        delta = d_parallel - d_perpendicular
        tensor_field[0, 0][cc_mask] = d_perpendicular + delta * n_cc[0]**2
        tensor_field[1, 1][cc_mask] = d_perpendicular + delta * n_cc[1]**2
        tensor_field[2, 2][cc_mask] = d_perpendicular + delta * n_cc[2]**2
        tensor_field[0, 1][cc_mask] = tensor_field[1, 0][cc_mask] = delta * n_cc[0] * n_cc[1]
        tensor_field[0, 2][cc_mask] = tensor_field[2, 0][cc_mask] = delta * n_cc[0] * n_cc[2]
        tensor_field[1, 2][cc_mask] = tensor_field[2, 1][cc_mask] = delta * n_cc[1] * n_cc[2]
    
    # Cingulum: anisotropic along y
    if np.any(cing_mask):
        delta = d_parallel - d_perpendicular
        tensor_field[0, 0][cing_mask] = d_perpendicular + delta * n_cing[0]**2
        tensor_field[1, 1][cing_mask] = d_perpendicular + delta * n_cing[1]**2
        tensor_field[2, 2][cing_mask] = d_perpendicular + delta * n_cing[2]**2
        tensor_field[0, 1][cing_mask] = tensor_field[1, 0][cing_mask] = delta * n_cing[0] * n_cing[1]
        tensor_field[0, 2][cing_mask] = tensor_field[2, 0][cing_mask] = delta * n_cing[0] * n_cing[2]
        tensor_field[1, 2][cing_mask] = tensor_field[2, 1][cing_mask] = delta * n_cing[1] * n_cing[2]
    
    # Verify SPD
    _verify_spd_tensor_field(tensor_field)
    
    return tensor_field


def _verify_spd_tensor_field(tensor_field: TensorField) -> bool:
    """Verify all tensors are symmetric positive-definite."""
    D, H, W = tensor_field.shape[2:]
    for i in range(D):
        for j in range(H):
            for k in range(W):
                t = np.array([
                    [tensor_field[0,0,i,j,k], tensor_field[0,1,i,j,k], tensor_field[0,2,i,j,k]],
                    [tensor_field[1,0,i,j,k], tensor_field[1,1,i,j,k], tensor_field[1,2,i,j,k]],
                    [tensor_field[2,0,i,j,k], tensor_field[2,1,i,j,k], tensor_field[2,2,i,j,k]],
                ])
                # Check symmetry
                if not np.allclose(t, t.T, atol=1e-6):
                    return False
                # Check positive definite
                eigs = np.linalg.eigvalsh(t)
                if eigs.min() <= 1e-12:
                    return False
    return True


# --------------------------------------------------------------------------- #
# Patient Data Organization
# --------------------------------------------------------------------------- #
def load_brats_patient(
    patient_dir: Union[str, Path],
    target_shape: Tuple[int, int, int] = DEFAULT_TARGET_SHAPE,
) -> Dict[str, np.ndarray]:
    """
    Load all modalities for a single BraTS patient.
    
    Expected directory structure:
    patient_dir/
        BRATS_XXXX_t1.nii.gz
        BRATS_XXXX_t1ce.nii.gz
        BRATS_XXXX_t2.nii.gz
        BRATS_XXXX_flair.nii.gz
        BRATS_XXXX_seg.nii.gz
    
    Returns:
        dict with keys: 't1', 't1ce', 't2', 'flair', 'seg', 
                        'tumor_mask', 'enhancing_mask', 'edema_mask', 'necrotic_mask',
                        'affine'
    """
    patient_dir = Path(patient_dir)
    
    # Find files
    files = {}
    for modality in ['t1', 't1ce', 't2', 'flair', 'seg']:
        matches = list(patient_dir.glob(f"*{modality}*.nii*"))
        if matches:
            files[modality] = matches[0]
    
    result = {}
    
    # Load segmentation first (reference affine)
    if 'seg' in files:
        seg_data, seg_affine, _ = load_nifti(files['seg'])
        result['seg'] = _reorder_to_dhw(seg_data)
        result['affine'] = seg_affine
        
        # Extract subregions
        result['tumor_mask'] = load_nifti_seg(files['seg'])
        result['necrotic_mask'] = load_nifti_subregion(files['seg'], 1)
        result['edema_mask'] = load_nifti_subregion(files['seg'], 2)
        result['enhancing_mask'] = load_nifti_subregion(files['seg'], 3)
    
    # Load modalities
    for modality in ['t1', 't1ce', 't2', 'flair']:
        if modality in files:
            data, affine, _ = load_nifti(files[modality])
            result[modality] = _reorder_to_dhw(data)
            if 'affine' not in result:
                result['affine'] = affine
    
    # Resample all to target grid
    if 'affine' in result:
        for key in ['t1', 't1ce', 't2', 'flair', 'seg', 'tumor_mask', 
                    'necrotic_mask', 'edema_mask', 'enhancing_mask']:
            if key in result:
                resampled, _ = resample_to_grid(
                    result[key], result['affine'], target_shape, order=0 if 'mask' in key else 1
                )
                result[key] = resampled
    
    return result


def discover_brats_patients(data_root: Union[str, Path]) -> List[Path]:
    """
    Discover BraTS patient directories in data root.
    
    Returns:
        List of patient directory paths
    """
    data_root = Path(data_root)
    patients = []
    
    # Look for directories containing segmentation files
    for seg_file in data_root.rglob("*seg*.nii*"):
        patient_dir = seg_file.parent
        if patient_dir not in patients:
            patients.append(patient_dir)
    
    return sorted(patients)


# --------------------------------------------------------------------------- #
# CLI / Testing
# --------------------------------------------------------------------------- #
def main():
    """Test with a sample file if available."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python 50_dicom_loader.py <path_to_nifti>")
        print("       python 50_dicom_loader.py --brats <brats_data_root>")
        return
    
    if sys.argv[1] == "--brats":
        if len(sys.argv) < 3:
            print("Provide BraTS data root path")
            return
        patients = discover_brats_patients(sys.argv[2])
        print(f"Found {len(patients)} BraTS patients")
        for p in patients[:5]:
            print(f"  {p}")
            try:
                data = load_brats_patient(p)
                print(f"    Keys: {list(data.keys())}")
                for k, v in data.items():
                    if isinstance(v, np.ndarray):
                        print(f"    {k}: shape={v.shape}, dtype={v.dtype}, range=[{v.min():.3f}, {v.max():.3f}]")
            except Exception as e:
                print(f"    Error: {e}")
    else:
        path = sys.argv[1]
        print(f"Loading: {path}")
        try:
            data, affine, meta = load_nifti(path)
            print(f"Shape: {data.shape}, Affine:\n{affine}")
            print(f"Meta: {meta}")
        except Exception as e:
            print(f"Error: {e}")


if __name__ == "__main__":
    main()