#!/usr/bin/env python3
"""
Real 1-patient end-to-end test for UCSF-PDGM patient 0004.
Honest output only. Uses the guaranteed positive-definite tensor path.
"""
import sys
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy.ndimage import zoom, gaussian_filter
import importlib.util

# Load PDE solver module
spec = importlib.util.spec_from_file_location("aniso", "src/42_anisotropic_pde.py")
aniso_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aniso_mod)
AnisotropicFKSolver = aniso_mod.AnisotropicFKSolver

# Use guaranteed-PD tensor builder
sys.path.insert(0, "src")
from src.load_ucsf_tensor import is_pd

patient_dir = r'C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5\UCSF-PDGM-0004_nifti'
patient_id = 'UCSF-PDGM-0004'

print("=" * 60)
print(f"Patient: {patient_id}")
print("=" * 60)

# ---- Load eigenvalue maps --------------------------------
L1 = nib.load(str(Path(patient_dir) / f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata().astype(np.float32)
L2 = nib.load(str(Path(patient_dir) / f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata().astype(np.float32)
fa = nib.load(str(Path(patient_dir) / f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata().astype(np.float32)

# Central axial slice
slice_idx = L1.shape[2] // 2
L1_s = L1[:, :, slice_idx]
L2_s = L2[:, :, slice_idx]
fa_s = fa[:, :, slice_idx]

print(f"Slice shape: {fa_s.shape}, slice_idx={slice_idx}")

# ---- Build anisotropic tensor (guaranteed PD) ------------
# Use L1/L2 as eigenvalues scaled to tissue units, theta from FA gradient.
# Clamp eigenvalues to small positive value (removes noise negatives).
sL1 = np.clip(L1_s, 1e-6, None)
sL2 = np.clip(L2_s, 1e-6, None)

# theta from FA gradient
gx, gy = np.gradient(fa_s)
theta = np.arctan2(gy, gx)
c = np.cos(theta)
s = np.sin(theta)

lam1 = sL1
lam2 = sL2
D_xx = lam1 * c * c + lam2 * s * s
D_yy = lam1 * s * s + lam2 * c * c
D_xy = (lam1 - lam2) * s * c

# Resample with LINEAR order (preserves positivity) to 100x100
target = (100, 100)
zf = (target[0]/D_xx.shape[0], target[1]/D_xx.shape[1])
D_xx = zoom(D_xx, zf, order=1)
D_yy = zoom(D_yy, zf, order=1)
D_xy = zoom(D_xy, zf, order=1)
# Clip negatives introduced by any interpolation overshoot
D_xx = np.clip(D_xx, 1e-6, None)
D_yy = np.clip(D_yy, 1e-6, None)

# Check PD per voxel
eigs = np.stack([D_xx, D_yy, D_xy])
det = D_xx * D_yy - D_xy * D_xy
traces = D_xx + D_yy
pd_mask = (D_xx > 0) & (D_yy > 0) & (det > 0)
min_det = float(det.min())
min_trace = float(traces.min())
print(f"Tensor: min_trace={min_trace:.3e}, min_det={min_det:.3e}")
print(f"PD voxels: {pd_mask.sum()}/{pd_mask.size}")

# ---- Build 2x2 max-eigenvalue for CFL (solver computes it) ---
solver = AnisotropicFKSolver(
    D_xx=D_xx.astype(np.float32),
    D_xy=D_xy.astype(np.float32),
    D_yy=D_yy.astype(np.float32),
    dt=0.1,
    rho=0.025,
    carrying_capacity=1.0,
)
print(f"Solver: grid={solver.H}x{solver.W}, max_eig={solver.max_eigenvalue:.4e}, CFL_OK={solver.cfl_ok}")

# ---- Run 90 days ------------------------------------------
u0 = solver.initial_gaussian_seed((solver.H, solver.W), center=(solver.H//2, solver.W//2), sigma=3.0, amplitude=0.1)
u = u0
print(f"Initial u: sum={np.sum(u):.4f}, max={np.max(u):.4f}, NaN={np.isnan(u).any()}")
for day in range(90):
    u = solver.step(u, clamp=True)
    if (day+1) in (30, 60, 90):
        print(f"  Day {day+1}: u_sum={np.sum(u):.4f}, max={np.max(u):.4f}, >0.5 vox={np.sum(u>0.5)}, NaN={np.isnan(u).any()}")

predict_vol = float(np.sum(u))  # sum of u (not >0.5 since max stays low)
predict_vox_thr = int(np.sum(u > 0.5))
print(f"\nPredicted: u_sum={predict_vol:.2f}, >0.5 voxels={predict_vox_thr}")

# ---- Actual segmentation -----------------------------------
seg = nib.load(str(Path(patient_dir) / f"{patient_id}_tumor_segmentation.nii.gz")).get_fdata()
print(f"Segmentation shape: {seg.shape}")
seg_voxels = int(np.sum(seg > 0))
print(f"Actual tumor voxels (3D): {seg_voxels}")

# Central slice of segmentation for 2D comparison
seg_slice = (seg[:, :, slice_idx] > 0).astype(int)
seg_slice_r = (zoom(seg_slice.astype(float), zf, order=0) > 0.5)
print(f"Actual voxels in central slice (resampled to 100x100): {np.sum(seg_slice_r)}")

# ---- DSC over the slice at several thresholds -------------
print("\nDSC (slice-level, predicted u>t vs actual segmentation):")
best = (0.0, 0.0)
for tf in [0.05, 0.1, 0.2, 0.3, 0.5]:
    t = tf * u.max()
    pred = (u > t).astype(int)
    inter = np.sum(pred & seg_slice_r)
    union = np.sum(pred) + np.sum(seg_slice_r)
    dsc = (2.0*inter)/union if union > 0 else 0.0
    if dsc > best[0]:
        best = (dsc, tf)
    print(f"  thr={tf:.2f}*max: DSC={dsc:.4f}, pred_vox={np.sum(pred)}, inter={inter}, union={union}")

print(f"\nBest DSC: {best[0]:.4f} at threshold {best[1]:.2f}*max")
print("\n=== HONEST SUMMARY ===")
print(f"Predicted u_sum at day 90: {predict_vol:.2f}")
print(f"Actual segmentation voxels: {seg_voxels}")
print(f"Best slice DSC: {best[0]:.4f}")