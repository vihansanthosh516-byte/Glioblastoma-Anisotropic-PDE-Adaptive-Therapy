#!/usr/bin/env python3
"""
Method 3: Full DTI tensor reconstruction from raw 4D DWI + rotated bvecs + bvals.
Fits the true diffusion tensor (Stejskal-Tanner log-linear LSQ) per voxel, gets
the real principal eigenvector, then tests anisotropic vs isotropic prediction
on patient 0004. Honest numbers only.
"""
import sys
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy.ndimage import zoom, distance_transform_edt, binary_erosion
import importlib.util

sys.path.insert(0, "src")
spec = importlib.util.spec_from_file_location("pde3d", "src/48_3d_extension.py")
pde3d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pde3d)
AnisotropicFKSolver3D = pde3d.AnisotropicFKSolver3D

UCSF_ROOT = Path(r"C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5")
BVAL = Path(r"C:\Users\vihan\Downloads\UCSF-bval-unzip\UCSF-PDGM_DTI.bval")
RHO = 0.0032
D_BASE = 0.013
EXTENT = 0.15
DAYS = 90
TARGET = (48, 74, 74)


def fit_tensor_field(patient_dir, pid):
    """Fit a 3x3 symmetric diffusion tensor per voxel via log-linear LSQ."""
    dwi = nib.load(str(patient_dir / f"{pid}_DTI_eddy_noreg.nii.gz")).get_fdata()
    bvec = np.array([l.split() for l in open(
        str(patient_dir / f"{pid}_DTI_eddy.eddy_rotated_bvecs")).read().splitlines()], dtype=float)
    bval = np.array([float(x) for x in open(str(BVAL)).read().split()])

    ny, nx, nz, nt = dwi.shape
    b0_idx = np.where(bval == 0)[0]
    dwi_idx = np.where(bval > 0)[0]
    S0 = dwi[..., b0_idx[0]]  # reference b0 volume

    # Stejskal-Tanner: ln(S/S0) = -b * g D g^T  (linear in D's 6 components)
    # Du = -(1/b) ln(S/S0) ; for tensor with dx,xx terms.
    bv = bval[dwi_idx]
    g = bvec[:, dwi_idx]  # (3, ndwi)
    # Design matrix rows: [gx^2, gy^2, gz^2, 2gxgy, 2gxgz, 2gygz]
    A = np.stack([
        g[0]**2, g[1]**2, g[2]**2,
        2*g[0]*g[1], 2*g[0]*g[2], 2*g[1]*g[2]
    ], axis=1)  # (ndwi, 6)
    AtA = A.T @ A
    AtA_inv = np.linalg.inv(AtA)

    # Fit at each voxel (vectorized over spatial dims)
    X = dwi[..., dwi_idx].astype(np.float64)  # (ny,nx,nz,ndwi)
    S0 = np.maximum(S0, 1e-12)
    X = np.maximum(X, 1e-12)
    y = -(1.0 / bv[None, None, None, :]) * np.log(X / S0[..., None])  # (ny,nx,nz,ndwi)
    # Dvec[..., k] = sum_b inv[k,b] * sum_d A.T[b,d]*y[...,d]
    Aty = np.tensordot(A.T, y, axes=([1], [3]))  # (6, ny,nx,nz)
    Dvec = np.tensordot(AtA_inv, Aty, axes=([1], [0]))  # (6, ny,nx,nz)
    Dvec = np.moveaxis(Dvec, 0, -1)  # (ny,nx,nz,6)
    # Dvec = [Dxx, Dyy, Dzz, Dxy, Dxz, Dyz]
    Dxx, Dyy, Dzz = Dvec[..., 0], Dvec[..., 1], Dvec[..., 2]
    Dxy, Dxz, Dyz = Dvec[..., 3], Dvec[..., 4], Dvec[..., 5]
    # Assemble (ny,nx,nz,3,3)
    shape = (ny, nx, nz)
    D = np.zeros(shape + (3, 3))
    D[..., 0, 0], D[..., 1, 1], D[..., 2, 2] = Dxx, Dyy, Dzz
    D[..., 0, 1] = D[..., 1, 0] = Dxy
    D[..., 0, 2] = D[..., 2, 0] = Dxz
    D[..., 1, 2] = D[..., 2, 1] = Dyz
    return D, dwi.shape[:3]


def resample_tensor(D, src_shape, target):
    # Resample each of 6 components from native DWI grid to TARGET
    zf = np.array(target) / np.array(src_shape)
    out = np.zeros(target + (3, 3))
    for i in range(3):
        for j in range(3):
            out[..., i, j] = zoom(D[..., i, j], zf, order=1)
    return out


def extent_metrics(u, mask, vox_vol, threshold=EXTENT):
    pred = u > threshold
    vol = float(np.sum(pred)) * vox_vol
    dsc = (2*np.sum(pred & mask))/(np.sum(pred)+np.sum(mask)) if np.sum(pred)+np.sum(mask) > 0 else 0.0
    nb = pred & ~binary_erosion(pred); tb = mask & ~binary_erosion(mask)
    h95p = np.percentile(distance_transform_edt(~mask)[nb], 95) if nb.any() else 0.0
    h95t = np.percentile(distance_transform_edt(~pred)[tb], 95) if tb.any() else 0.0
    return vol, dsc, float(max(h95p, h95t))


pid = "UCSF-PDGM-0004"
patient_dir = UCSF_ROOT / f"{pid}_nifti"
seg_nii = nib.load(str(patient_dir / f"{pid}_tumor_segmentation.nii.gz"))
seg_orig = seg_nii.get_fdata()
vox_vol = abs(np.linalg.det(seg_nii.affine[:3, :3]))
zf = np.array(TARGET) / np.array(seg_orig.shape)
mask = zoom(seg_orig, zf, order=0) > 0.5

print("Fitting true diffusion tensor (Method 3)...")
Df, src_shape = fit_tensor_field(patient_dir, pid)
print(f"Fitted tensor field from 4D DWI; native grid {src_shape}")
DT = resample_tensor(Df, src_shape, TARGET)
eig = np.linalg.eigvalsh(DT)
print(f"Tensor min eigenvalue (mm2/s): {eig.min():.2e}  PD={eig.min() > -1e-12}")

# Rescale so mean tumor diffusivity -> D_BASE mm2/day
# Use nonzero-brain mask from segmentation upsampled? Use all voxels for scale.
md_t = (DT[...,0,0]+DT[...,1,1]+DT[...,2,2])/3.0
ok = md_t > 1e-9
md_mean = md_t[ok].mean() if ok.any() else 1e-6
scale = D_BASE / md_mean
DT = DT * scale
print(f"md_mean={md_mean:.3e} mm2/s -> scale={scale:.1f} (=> MD {D_BASE} mm2/day)")

def run(Dt, aniso):
    if aniso:
        D_xx, D_yy, D_zz = Dt[...,0,0], Dt[...,1,1], Dt[...,2,2]
        D_xy, D_xz, D_yz = Dt[...,0,1], Dt[...,0,2], Dt[...,1,2]
    else:
        Dm = (Dt[...,0,0]+Dt[...,1,1]+Dt[...,2,2])/3.0
        D_xx=D_yy=D_zz=Dm; D_xy=D_xz=D_yz=np.zeros_like(Dm)
    solver = AnisotropicFKSolver3D(
        D_xx=D_xx.astype(np.float64), D_xy=D_xy.astype(np.float64), D_xz=D_xz.astype(np.float64),
        D_yy=D_yy.astype(np.float64), D_yz=D_yz.astype(np.float64), D_zz=D_zz.astype(np.float64),
        dt=0.1, dx=1.0, rho=RHO, K=1.0)
    u0 = np.zeros(TARGET, np.float64); u0[mask] = 0.9
    u = u0.copy()
    for _ in range(int(round(DAYS/solver.dt))):
        u = solver.step(u, 0.0)
    vol, dsc, hd = extent_metrics(u, mask, vox_vol)
    return dsc, hd, vol, float(u.max())

print(f"\nPatient {pid}: TRUE TENSOR (Method 3) aniso vs iso\n")
dsc_a, hd_a, vol_a, mxa = run(DT, True)
dsc_i, hd_i, vol_i, mxi = run(DT, False)
print(f"  true_tensor aniso : DSC={dsc_a:.4f}  HD95={hd_a:.2f} mm  vol90={vol_a:.0f}  max_u={mxa:.3f}")
print(f"  iso               : DSC={dsc_i:.4f}  HD95={hd_i:.2f} mm  vol90={vol_i:.0f}  max_u={mxi:.3f}")

# Summary vs prior methods
print("\n=== COMPARISON (patient 0004) ===")
print(f"  true_tensor aniso: DSC={dsc_a:.4f}")
print(f"  l1_grad    aniso: DSC=0.6736 (prior)")
print(f"  fa_grad    aniso: DSC=0.6713 (prior)")
print(f"  isotropic       : DSC={dsc_i:.4f}")