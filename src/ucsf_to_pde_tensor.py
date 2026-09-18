#!/usr/bin/env python3
"""
Convert UCSF-PDGM 3D DTI data to 2D tensor fields for the PDE solver (42_anisotropic_pde.py).

The solver expects D_xx, D_xy, D_yy built via:
    D = R(theta) diag(lam1, lam2) R(-theta)
    D_xx = lam1*cos^2(theta) + lam2*sin^2(theta)
    D_yy = lam1*sin^2(theta) + lam2*cos^2(theta)
    D_xy = (lam1-lam2)*sin(theta)*cos(theta)
"""
import numpy as np
import nibabel as nib
import sys
from pathlib import Path
from scipy.ndimage import gaussian_filter, zoom

# Load the PDE solver module by executing its source
import importlib.util
spec = importlib.util.spec_from_file_location("aniso_mod", "src/42_anisotropic_pde.py")
aniso_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aniso_mod)
AnisotropicFKSolver = aniso_mod.AnisotropicFKSolver
TensorFieldBuilder = aniso_mod.TensorFieldBuilder


def ucsf_to_pde_tensor_fields(
    patient_dir: str | Path,
    slice_axis: int = 2,  # 0=sagittal, 1=coronal, 2=axial (default)
    slice_idx: Optional[int] = None,
    target_grid: int = 100,  # PDE solver grid size (GRID_SIZE)
) -> Dict[str, np.ndarray]:
    """
    Convert UCSF-PDGM 3D DTI data to 2D tensor fields for the PDE solver.

    Returns dict with keys: D_xx, D_xy, D_yy (each shape: (grid_size, grid_size))
    """
    patient_dir = Path(patient_dir)
    patient_id = patient_dir.name.replace("_nifti", "")

    # Load eigenvalue maps
    L1 = nib.load(str(patient_dir / f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata().astype(np.float32)
    L2 = nib.load(str(patient_dir / f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata().astype(np.float32)
    L3 = nib.load(str(patient_dir / f"{patient_id}_DTI_eddy_L3.nii.gz")).get_fdata().astype(np.float32)
    fa = nib.load(str(patient_dir / f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata().astype(np.float32)
    md = nib.load(str(patient_dir / f"{patient_id}_DTI_eddy_MD.nii.gz")).get_fdata().astype(np.float32)

    # Extract central slice along specified axis
    if slice_idx is None:
        slice_idx = L1.shape[slice_axis] // 2

    if slice_axis == 0:  # sagittal
        L1_slice = L1[slice_idx, :, :]
        L2_slice = L2[slice_idx, :, :]
        L3_slice = L3[slice_idx, :, :]
        fa_slice = fa[slice_idx, :, :]
    elif slice_axis == 1:  # coronal
        L1_slice = L1[:, slice_idx, :]
        L2_slice = L2[:, slice_idx, :]
        L3_slice = L3[:, slice_idx, :]
        fa_slice = fa[:, slice_idx, :]
    else:  # axial (default)
        L1_slice = L1[:, :, slice_idx]
        L2_slice = L2[:, :, slice_idx]
        L3_slice = L3[:, :, slice_idx]
        fa_slice = fa[:, :, slice_idx]

    # Compute principal eigenvector direction from FA gradient
    # Gaussian smooth first to reduce noise
    fa_smooth = gaussian_filter(fa_slice, sigma=1.0)
    gx, gy = np.gradient(fa_smooth)

    # Principal eigenvector direction angle theta = atan2(gy, gx)
    theta = np.arctan2(gy, gx)  # shape: (M, N)

    # Eigenvalues: use L1 (max) and L2 (min) from the slice
    lam1 = L1_slice
    lam2 = L2_slice

    # Create tract mask based on FA threshold (typical GBM FA > 0.2)
    fa_threshold = 0.2
    tract_mask = fa_slice > fa_threshold

    # Set eigenvalues: inside tract use target values, outside use base
    d_parallel = 0.013  # D_white from solver: 0.013 mm²/day
    d_perpendicular = 0.0013  # D_perpendicular ≈ D_white/10
    d_base = 0.0013  # base isotropic diffusivity

    lam1_field = np.where(tract_mask, d_parallel, d_base)
    lam2_field = np.where(tract_mask, d_perpendicular, d_base)

    # Compute tensor components using same formulas as TensorFieldBuilder
    c = np.cos(theta)
    s = np.sin(theta)

    D_xx = lam1_field * c * c + lam2_field * s * s
    D_yy = lam1_field * s * s + lam2_field * c * c
    D_xy = (lam1_field - lam2_field) * s * c
    D_yx = D_xy  # symmetry by construction

    # Resize to target grid (100x100 = GRID_SIZE) using area interpolation
    target = (target_grid, target_grid)

    D_xx = zoom(D_xx, (target[0] / D_xx.shape[0], target[1] / D_xx.shape[1]), order=3)
    D_yy = zoom(D_yy, (target[0] / D_yy.shape[0], target[1] / D_yy.shape[1]), order=3)
    D_xy = zoom(D_xy, (target[0] / D_xy.shape[0], target[1] / D_xy.shape[1]), order=3)

    return {
        'D_xx': D_xx.astype(np.float32),
        'D_xy': D_xy.astype(np.float32),
        'D_yy': D_yy.astype(np.float32),
        'theta': theta,
        'lam1': lam1_field,
        'lam2': lam2_field,
        'fa_mean': float(fa_slice.mean()),
        'patient_id': patient_id,
        'slice_axis': slice_axis,
        'slice_idx': slice_idx,
    }


def run_pde_with_ucsf_tensor(patient_dir: str | Path, output_dir: str | Path) -> Dict[str, Any]:
    """
    Full pipeline: load UCSF data -> build tensor fields -> run PDE simulation.
    """
    patient_dir = Path(patient_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Build tensor fields from UCSF data
    print(f"Building tensor fields for {patient_dir.name}...")
    tensor_fields = ucsf_to_pde_tensor_fields(patient_dir)

    # Step 2: Create TensorFieldBuilder and override with UCSF tensors
    builder = TensorFieldBuilder()  # default 100x100 grid with tract

    # Replace builder's tensor fields with UCSF-derived ones
    builder.D_xx = tensor_fields['D_xx']
    builder.D_yy = tensor_fields['D_yy']
    builder.D_xy = tensor_fields['D_xy']
    builder.D_yx = tensor_fields['D_xy']  # symmetry by construction
    builder.lambda_1 = tensor_fields['lam1']
    builder.lambda_2 = tensor_fields['lam2']

    # Validate the constructed tensor
    val_metrics = builder.validate_tensor()

    # Step 3: Run PDE simulation
    print("Running PDE simulation...")
    solver = AnisotropicFKSolver(
        D_xx=builder.D_xx,
        D_xy=builder.D_xy,
        D_yy=builder.D_yy,
        dt=0.1,  # days
        rho=0.025,  # typical GBM growth rate
        carrying_capacity=1.0,
    )

    # Run 90 steps (90 days) using the step method
    print("Running 90-day simulation (step-by-step)...")
    for step in range(90):
        solver.step()
    
    # Save results
    save_path = output_dir / f"{patient_dir.name}_ucsf_pde_results.npz"
    np.savez(save_path,
             D_xx=builder.D_xx,
             D_yy=builder.D_yy,
             D_xy=builder.D_xy,
             lambda_1=builder.lambda_1,
             lambda_2=builder.lambda_2,
             final_volume=float(solver.u.sum()),
             cfl_ok=builder.cfl_ok,
             validation_metrics=val_metrics,
    )

    print(f"Saved results to {save_path}")
    print(f"Validation: sym_err={val_metrics['symmetry_max_error']:.3e}, "
          f"min_eig={val_metrics['min_eigenvalue']:.3e}")

    return {
        'tensor_fields': tensor_fields,
        'validation': val_metrics,
        'pde_results': solver.u,
        'save_path': save_path,
    }


if __name__ == "__main__":
    # Test on patient 0004
    patient_dir = Path(r"C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5\UCSF-PDGM-0004_nifti")

    print("=" * 60)
    print("UCSF-PDGM to PDE Solver Integration Test")
    print("=" * 60)

    result = run_pde_with_ucsf_tensor(patient_dir, Path("output"))

    print("\n=== Integration Test Complete ===")
    print(f"Patient: {result['tensor_fields']['patient_id']}")
    print(f"FA mean: {result['tensor_fields']['fa_mean']:.4f}")
    print(f"Tensor validation symmetry error: {result['validation']['symmetry_max_error']:.2e}")
    print(f"Positive-definite: {'PASS' if result['validation']['positive_definite_pass'] else 'FAIL'}")
    print(f"Final volume: {result['pde_results'].sum():.4f}")
    print(f"CFL OK: {result['validation']['cfl_ok']}")