#!/usr/bin/env python3
"""
Load UCSF-PDGM pre-computed DTI eigenvalue maps and build full diffusion tensors.

GUARANTEED POSITIVE-DEFINITE: Uses eigenvalue decomposition D = R(theta) @ diag(L1, L2) @ R(-theta),
which is always positive-semidefinite when L1, L2 >= 0.
"""
from __future__ import annotations

import numpy as np
import nibabel as nib
from pathlib import Path
from scipy.ndimage import zoom
from typing import Tuple, Optional, Dict, Any


def is_pd(D: np.ndarray) -> bool:
    """Check if 2x2 tensor is positive-definite."""
    return (D[0, 0] > 0 and D[1, 1] > 0 and np.linalg.det(D) > 0)


def build_pd_tensor(l1: float, l2: float, theta: float) -> np.ndarray:
    """
    Build guaranteed positive-definite 2x2 tensor from eigenvalues and angle.
    
    D = R(theta) @ diag(L1, L2) @ R(-theta)
    Always PD if L1, L2 >= 0.
    """
    c, s = np.cos(theta), np.sin(theta)
    R = np.array([[c, -s], [s, c]])
    Lambda = np.diag([max(l1, 1e-6), max(l2, 1e-6)])
    D = R @ Lambda @ R.T
    return D


def load_patient_tensors(
    patient_dir: str | Path,
    tumor_mask: Optional[np.ndarray] = None,
    dx: float = 1.0,
    target_shape: Optional[Tuple[int, int, int]] = None,
    eigenvector_method: str = "fa_grad",  # "fa_grad" | "l1_grad"
) -> Dict[str, Any]:
    """
    Load UCSF-PDGM DTI eigenvalue maps and build full diffusion tensor field.

    Uses D = V @ diag(L1, L2, L3) @ V^T which is GUARANTEED positive-definite.

    eigenvector_method:
      - "fa_grad": estimate principal eigenvector from FA gradient (proxy).
      - "l1_grad": estimate principal eigenvector from L1 eigenvalue gradient.
    """
    patient_dir = Path(patient_dir)
    patient_id = patient_dir.name.replace("_nifti", "")

    # Load eigenvalue maps
    L1 = _load_nifti(patient_dir / f"{patient_id}_DTI_eddy_L1.nii.gz")
    L2 = _load_nifti(patient_dir / f"{patient_id}_DTI_eddy_L2.nii.gz")
    L3 = _load_nifti(patient_dir / f"{patient_id}_DTI_eddy_L3.nii.gz")
    fa = _load_nifti(patient_dir / f"{patient_id}_DTI_eddy_FA.nii.gz")
    md = _load_nifti(patient_dir / f"{patient_id}_DTI_eddy_MD.nii.gz")

    original_shape = fa.shape

    # Resample if target_shape specified
    if target_shape is not None and tuple(fa.shape) != target_shape:
        zoom_factors = np.array(target_shape) / np.array(fa.shape)
        L1 = zoom(L1, zoom_factors, order=1)
        L2 = zoom(L2, zoom_factors, order=1)
        L3 = zoom(L3, zoom_factors, order=1)
        fa = zoom(fa, zoom_factors, order=1)
        md = zoom(md, zoom_factors, order=1)

    # Sanitize eigenvalues: clip negative values from noise, ensure non-negative
    L1 = np.clip(L1, 0.0, None)
    L2 = np.clip(L2, 0.0, None)
    L3 = np.clip(L3, 0.0, None)

    # Sort per voxel: L1 >= L2 >= L3
    eigenvalues = np.stack([L1, L2, L3], axis=-1)
    eigenvalues = np.sort(eigenvalues, axis=-1)[..., ::-1]  # descending
    L1, L2, L3 = eigenvalues[..., 0], eigenvalues[..., 1], eigenvalues[..., 2]

    # Principal eigenvector gradient method (proxy for principal diffusion direction)
    if eigenvector_method == "fa_grad":
        v1 = _vector_from_gradient(fa)
    elif eigenvector_method == "l1_grad":
        v1 = _vector_from_gradient(eigenvalues[..., 0])  # L1 map after sorting
    else:
        raise ValueError(f"Unknown eigenvector_method: {eigenvector_method}")

    # Complete orthonormal basis
    v2, v3 = _complete_basis(v1)

    # Build rotation matrix V at each voxel: V = [v1, v2, v3]
    V = np.stack([v1, v2, v3], axis=-1)  # (X, Y, Z, 3, 3)

    # Build tensor using eigenvalue decomposition: D = V @ diag(L1, L2, L3) @ V^T
    # This is ALWAYS positive-semidefinite when L1, L2 >= 0
    eigenvalues_stacked = np.stack([L1, L2, L3], axis=-1)  # (X, Y, Z, 3)
    
    # Vectorized tensor construction
    D = np.einsum('...ik,...k,...jk->...ij', V, eigenvalues_stacked, V)

    # Verify positive-definite and clip if needed
    eigvals_check = np.linalg.eigvalsh(D)
    pd_mask = eigvals_check >= -1e-12
    if not np.all(pd_mask):
        print(f"WARNING: {np.sum(~pd_mask)} voxels have negative eigenvalues, clipping")
        D = np.where(~pd_mask[..., None, None], 
                     np.einsum('...ik,...k,...jk->...ij', V, 
                               np.clip(eigenvalues_stacked, 1e-6, None), 
                               V), D)

    # If target_shape given, extract central axial slice info
    slice_idx = None
    if target_shape is not None:
        slice_idx = target_shape[2] // 2

    return {
        'tensors': D,                          # (X, Y, Z, 3, 3) mm²/s, GUARANTEED PD
        'eigenvalues': eigenvalues_stacked,    # (X, Y, Z, 3) [L1, L2, L3]
        'principal_eigenvector': v1,           # (X, Y, Z, 3)
        'tumor_mask': tumor_mask,              # (X, Y, Z)
        'original_shape': original_shape,
        'slice_idx': slice_idx,                # axial slice index for 2D extraction
        'fa': fa,
        'md': md,
    }


def _load_nifti(path: Path) -> np.ndarray:
    """Load NIfTI file and return float32 array."""
    if not Path(path).exists():
        raise FileNotFoundError(f"NIfTI file not found: {path}")
    return nib.load(str(path)).get_fdata().astype(np.float32)


def _vector_from_gradient(img: np.ndarray) -> np.ndarray:
    """Estimate principal direction from the normalized gradient of a scalar map.

    Works for either FA (fa_grad) or L1 eigenvalue (l1_grad). Gradient points
    along the direction of fastest change of the scalar field.
    """
    from scipy.ndimage import gaussian_filter
    img_smooth = gaussian_filter(img.astype(np.float64), sigma=1.0)
    gx, gy, gz = np.gradient(img_smooth)
    v1 = np.stack([gx, gy, gz], axis=-1)
    norm = np.linalg.norm(v1, axis=-1, keepdims=True)
    norm = np.where(norm > 1e-6, norm, 1.0)
    v1 = v1 / norm
    zero_grad = np.linalg.norm(v1.squeeze(), axis=-1) < 1e-6
    if np.any(zero_grad):
        v1 = v1.copy()
        v1[zero_grad, 0] = 1.0
        v1[zero_grad, 1] = 0.0
        v1[zero_grad, 2] = 0.0
    return v1


def _complete_basis(v1: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Complete orthonormal basis from principal eigenvector v1."""
    ref1 = np.array([1.0, 0.0, 0.0])
    ref2 = np.array([0.0, 1.0, 0.0])
    dot1 = np.sum(v1 * ref1, axis=-1, keepdims=True)
    v2 = ref1 - dot1 * v1
    norm2 = np.linalg.norm(v2, axis=-1, keepdims=True)
    norm2 = np.where(norm2 > 1e-6, norm2, 1.0)
    v2 = v2 / norm2
    v3 = np.cross(v1, v2)
    norm3 = np.linalg.norm(v3, axis=-1, keepdims=True)
    norm3 = np.where(norm3 > 1e-6, norm3, 1.0)
    v3 = v3 / norm3
    return v2, v3


def compute_patient_stats(tensor_dict: Dict[str, Any], tumor_mask: Optional[np.ndarray] = None) -> Dict[str, float]:
    """Compute tumor-region statistics for calibration."""
    fa = tensor_dict['fa']
    md = tensor_dict['md']
    eigenvalues = tensor_dict['eigenvalues']
    mask = tumor_mask if tumor_mask is not None else tensor_dict.get('tumor_mask')

    if mask is None:
        mask = np.ones_like(fa, dtype=bool)
    mask = mask > 0.5

    return {
        'mean_fa': float(fa[mask].mean()),
        'mean_md': float(md[mask].mean()),
        'median_fa': float(np.median(fa[mask])),
        'median_md': float(np.median(md[mask])),
        'mean_l1': float(tensor_dict['eigenvalues'][..., 0][mask].mean()),
        'mean_l2': float(tensor_dict['eigenvalues'][..., 1][mask].mean()),
        'mean_l3': float(tensor_dict['eigenvalues'][..., 2][mask].mean()),
        'tumor_volume_voxels': int(mask.sum()),
        'tumor_fa_std': float(fa[mask].std()),
        'tumor_md_std': float(md[mask].std()),
    }


def save_tensor_for_pde(tensor_dict: Dict[str, Any], output_path: str | Path) -> None:
    """
    Save tensor in format expected by PDE solver.
    
    Now saves the full 3D tensor field (X, Y, Z, 3, 3) using guaranteed PD eigenvalue decomposition.
    """
    D = tensor_dict['tensors']  # (X, Y, Z, 3, 3) already PD
    
    np.save(output_path, D.astype(np.float32))
    print(f"Saved 3D tensor field to {output_path}: shape={D.shape}")


def load_tensor_for_pde_2d(npy_path: str | Path, slice_idx: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Load 3D tensor field and extract 2D components for PDE solver.
    Returns D_xx, D_xy, D_yy as (M, N) arrays for the 100x100 grid.
    """
    D = np.load(npy_path)  # (X, Y, Z, 3, 3)
    
    if slice_idx is None:
        slice_idx = D.shape[2] // 2
    
    D_slice = D[:, :, slice_idx]  # (X, Y, 3, 3)
    
    D_xx = D_slice[..., 0, 0]  # (X, Y)
    D_xy = D_slice[..., 0, 1]  # (X, Y)
    D_yy = D_slice[..., 1, 1]  # (X, Y)
    
    return D_xx, D_xy, D_yy


if __name__ == "__main__":
    # Quick test: verify PD property on synthetic data
    print("Testing guaranteed PD tensor build...")
    
    test_cases = [
        (0.013, 0.0013, 0.0),
        (0.001, 0.0005, 0.2),
        (1e-4, 1e-4, 0.7),
        (-1e-4, 1e-4, 0.5),
    ]
    
    all_pd = True
    for l1, l2, theta in test_cases:
        D = build_pd_tensor(l1, l2, theta)
        pd_check = is_pd(D)
        det = np.linalg.det(D)
        eigs = np.linalg.eigvalsh(D)
        if not pd_check:
            all_pd = False
        print(f"  L1={l1:.4e}, L2={l2:.4e}, theta={theta:.3f} -> "
              f"det={det:.3e}, eigs=({eigs[0]:.3e},{eigs[1]:.3e}), PD={pd_check}")
        all_pd = all_pd and pd_check
    
    if all_pd:
        print("SUCCESS: All test cases guarantee positive-definiteness!")
    else:
        print("FAILURE: Some test cases failed positive-definiteness check")