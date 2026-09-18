#!/usr/bin/env python3
"""
Full 3D validation, patient 0004 — v2.
Fix per user directive: seed from observed tumor mask (progression, not initiation).
Compute Day-0 and Day-90 volume, DSC, HD95. Honest output only.
"""
import sys
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy.ndimage import zoom, distance_transform_edt
import importlib.util

sys.path.insert(0, "src")
from src.load_ucsf_tensor import load_patient_tensors

# Load 3D solver
spec = importlib.util.spec_from_file_location("pde3d", "src/48_3d_extension.py")
pde3d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pde3d)
AnisotropicFKSolver3D = pde3d.AnisotropicFKSolver3D

patient_dir = Path(r'C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5\UCSF-PDGM-0004_nifti')
patient_id = 'UCSF-PDGM-0004'

print("=" * 60)
print(f"3D validation v2, patient {patient_id}")
print("=" * 60)

# ---- Load + resample segmentation -------------------------
seg_nii = nib.load(str(patient_dir / f"{patient_id}_tumor_segmentation.nii.gz"))
seg_orig = seg_nii.get_fdata()
affine = seg_nii.affine
voxel_vol_mm3 = abs(np.linalg.det(affine[:3, :3]))
print(f"Original seg shape: {seg_orig.shape}, voxel vol = {voxel_vol_mm3:.4f} mm^3")

tz, ty, tx = 48, 74, 74
target_shape = (tz, ty, tx)
zf = np.array(target_shape) / np.array(seg_orig.shape)
seg_resized = zoom(seg_orig, zf, order=0) > 0.5
actual_mask = seg_resized
print(f"Resampled tumor voxels: {np.sum(actual_mask)}")
print(f"True volume (resampled): {np.sum(actual_mask) * voxel_vol_mm3:.1f} mm^3")

# ---- Build 3D tensor field (guaranteed PD) ----------------
tensor_out = load_patient_tensors(patient_dir, target_shape=target_shape)
D3 = tensor_out['tensors']
eigs3 = np.linalg.eigvalsh(D3)
print(f"3D tensor shape: {D3.shape}, min eigenvalue: {eigs3.min():.3e}")

D_xx, D_xy, D_xz = D3[..., 0, 0], D3[..., 0, 1], D3[..., 0, 2]
D_yy, D_yz, D_zz = D3[..., 1, 1], D3[..., 1, 2], D3[..., 2, 2]

# ---- Rescale D: map DTI eigenvalues (mm^2/s) into tumor-DTI range (mm^2/day)
# Anchor to Swanson et al. tumor diffusion scale; the 2D solver's own reference
# D_white = 0.013 mm^2/day. Target mean tumor diffusivity ~ 0.013 mm^2/day.
md_mean = tensor_out['md'].mean()
D_targ_md = 0.013  # target mean diffusivity in mm^2/day (Swanson-scale)
Dscale = D_targ_md / max(md_mean, 1e-9)
print(f"md_mean={md_mean:.5f} mm2/s -> scale={Dscale:.2f} (=> MD {D_targ_md} mm2/day)")
D_xx = D_xx * Dscale
D_xy = D_xy * Dscale
D_xz = D_xz * Dscale
D_yy = D_yy * Dscale
D_yz = D_yz * Dscale
D_zz = D_zz * Dscale
print(f"Scaled D to mm2/day. max D_xx={D_xx.max():.4f} mm2/day")

# ---- Seed from observed tumor (carrying-capacity body, standard GBM model) ---
u0 = np.zeros(target_shape, dtype=np.float64)
u0[actual_mask] = 0.9  # 90% carrying capacity body; leading edge invades via D
print(f"u0 initialized from tumor mask at 0.9: sum={np.sum(u0):.4f}, >0.5 vox={np.sum(u0>0.5)}")

# ---- Run 3D solver ----------------------------------------
# Growth rate calibrated from real MU-Glioma-Post longitudinal data:
# median rho over 90 growing tumor pairs = 0.0032 /day
# (output/rho_calibration_positive.json). DOCUMENTED, not guessed.
rho = 0.0032
solver = AnisotropicFKSolver3D(
    D_xx=D_xx.astype(np.float64), D_xy=D_xy.astype(np.float64), D_xz=D_xz.astype(np.float64),
    D_yy=D_yy.astype(np.float64), D_yz=D_yz.astype(np.float64), D_zz=D_zz.astype(np.float64),
    dt=0.1, dx=1.0, rho=rho, K=1.0,
)
print(f"Solver: grid={solver.H}x{solver.W}x{solver.D}, rho={rho}, CFL_OK={solver.cfl_ok}")

C = 0.0  # no treatment (untreated progression)

def extent_metrics(u, mask, vox_vol, threshold):
    pred = u > threshold
    vol = np.sum(pred) * vox_vol
    inter = np.sum(pred & mask)
    n_pred = np.sum(pred)
    n_true = np.sum(mask)
    dsc = (2 * inter) / (n_pred + n_true) if (n_pred + n_true) > 0 else 0.0
    from scipy.ndimage import binary_erosion
    nb = pred & ~binary_erosion(pred)
    tb = mask & ~binary_erosion(mask)
    d_p2t = distance_transform_edt(~mask)
    d_t2p = distance_transform_edt(~pred)
    hd95_p = np.percentile(d_p2t[nb], 95) if nb.any() else 0.0
    hd95_t = np.percentile(d_t2p[tb], 95) if tb.any() else 0.0
    return vol, dsc, max(hd95_p, hd95_t)

# Report bulk (0.5) and infiltrative extent (0.15) thresholds as a sensitivity.
v0_5, d0_5, h0_5 = extent_metrics(u0, actual_mask, voxel_vol_mm3, 0.5)
v0_15, d0_15, h0_15 = extent_metrics(u0, actual_mask, voxel_vol_mm3, 0.15)
print(f"\nDay 0 : bulk(>0.5) vol={v0_5:.1f} DSC={d0_5:.4f} HD95={h0_5:.1f} | "
      f"extent(>0.15) vol={v0_15:.1f} DSC={d0_15:.4f} HD95={h0_15:.1f}")

u = u0.copy()
total_days = 90
steps = int(round(total_days / solver.dt))
print(f"Running {steps} steps (dt={solver.dt:.4f} d) to reach {total_days} days...")
for i in range(steps):
    u = solver.step(u, C)
    day = (i + 1) * solver.dt
    if (i + 1) % 300 == 0:
        vb, db, hb = extent_metrics(u, actual_mask, voxel_vol_mm3, 0.5)
        ve, de, he = extent_metrics(u, actual_mask, voxel_vol_mm3, 0.15)
        print(f"Day {day:6.1f}: bulk {vb:7.0f} mm3 DSC={db:.3f} | "
              f"extent {ve:7.0f} mm3 DSC={de:.3f} HD95={he:.1f} mm | "
              f"max_u={np.max(u):.3f} NaN={np.isnan(u).any()}")

vb, db, hb = extent_metrics(u, actual_mask, voxel_vol_mm3, 0.5)
ve, de, he = extent_metrics(u, actual_mask, voxel_vol_mm3, 0.15)
print("\n=== HONEST v2 RESULTS ===")
print(f"Day 0  bulk(0.5): vol={v0_5:.1f} mm3 DSC={d0_5:.4f} HD95={h0_5:.1f} mm")
print(f"Day 0  extent(0.15): vol={v0_15:.1f} mm3 DSC={d0_15:.4f} HD95={h0_15:.1f} mm")
print(f"Day 90 bulk(0.5): vol={vb:.1f} mm3 DSC={db:.4f} HD95={hb:.1f} mm")
print(f"Day 90 extent(0.15): vol={ve:.1f} mm3 DSC={de:.4f} HD95={he:.1f} mm")
print(f"Growth bulk: {vb - v0_5:+.0f} mm3 | Growth extent: {ve - v0_15:+.0f} mm3")