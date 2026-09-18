#!/usr/bin/env python3
"""
Single-patient A/B test: L1-gradient anisotropic tensor vs isotropic baseline.
Tests whether L1-gradient eigenvector construction captures anisotropy better
than the FA-gradient proxy. Patient 0004. Real numbers only.
"""
import sys
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy.ndimage import zoom, distance_transform_edt, binary_erosion
import importlib.util

sys.path.insert(0, "src")
from src.load_ucsf_tensor import load_patient_tensors

spec = importlib.util.spec_from_file_location("pde3d", "src/48_3d_extension.py")
pde3d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pde3d)
AnisotropicFKSolver3D = pde3d.AnisotropicFKSolver3D

UCSF_ROOT = Path(r"C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5")
RHO = 0.0032
D_BASE = 0.013
EXTENT = 0.15
DAYS = 90
TARGET = (48, 74, 74)


def extent_metrics(u, mask, vox_vol, threshold=EXTENT):
    pred = u > threshold
    vol = float(np.sum(pred)) * vox_vol
    dsc = (2 * np.sum(pred & mask)) / (np.sum(pred) + np.sum(mask)) if np.sum(pred) + np.sum(mask) > 0 else 0.0
    nb = pred & ~binary_erosion(pred)
    tb = mask & ~binary_erosion(mask)
    h95p = np.percentile(distance_transform_edt(~mask)[nb], 95) if nb.any() else 0.0
    h95t = np.percentile(distance_transform_edt(~pred)[tb], 95) if tb.any() else 0.0
    return vol, dsc, float(max(h95p, h95t))


def run(pid, method):
    patient_dir = UCSF_ROOT / f"{pid}_nifti"
    seg_nii = nib.load(str(patient_dir / f"{pid}_tumor_segmentation.nii.gz"))
    seg_orig = seg_nii.get_fdata()
    vox_vol = abs(np.linalg.det(seg_nii.affine[:3, :3]))
    zf = np.array(TARGET) / np.array(seg_orig.shape)
    mask = zoom(seg_orig, zf, order=0) > 0.5

    # Load tensor once with the given eigenvector method (l1_grad or fa_grad).
    # iso is derived from the same anisotropic tensor (mean diffusivity) so the
    # comparison isolates anisotropy.
    evm = method if method != "iso" else "l1_grad"
    out = load_patient_tensors(patient_dir, target_shape=TARGET, eigenvector_method=evm)
    D3 = out['tensors'] * (D_BASE / max(out['md'].mean(), 1e-9))

    if method == "iso":
        Dm = (D3[..., 0, 0] + D3[..., 1, 1] + D3[..., 2, 2]) / 3.0
        D_xx = D_yy = D_zz = Dm
        D_xy = D_xz = D_yz = np.zeros_like(Dm)
    else:
        D_xx, D_yy, D_zz = D3[..., 0, 0], D3[..., 1, 1], D3[..., 2, 2]
        D_xy, D_xz, D_yz = D3[..., 0, 1], D3[..., 0, 2], D3[..., 1, 2]

    solver = AnisotropicFKSolver3D(
        D_xx=D_xx.astype(np.float64), D_xy=D_xy.astype(np.float64), D_xz=D_xz.astype(np.float64),
        D_yy=D_yy.astype(np.float64), D_yz=D_yz.astype(np.float64), D_zz=D_zz.astype(np.float64),
        dt=0.1, dx=1.0, rho=RHO, K=1.0)
    u0 = np.zeros(TARGET, np.float64); u0[mask] = 0.9
    u = u0.copy()
    for _ in range(int(round(DAYS / solver.dt))):
        u = solver.step(u, 0.0)
    vol, dsc, hd = extent_metrics(u, mask, vox_vol)
    return dsc, hd, vol, float(u.max())


pid = "UCSF-PDGM-0004"
print(f"Patient {pid}: L1-gradient aniso vs isotropic (rho={RHO}, D->{D_BASE} mm2/day, extent>{EXTENT})\n")
for method in ["l1_grad", "iso"]:
    dsc, hd, vol, mx = run(pid, method)
    print(f"  {method:8s}: DSC={dsc:.4f}  HD95={hd:.2f} mm  vol90={vol:.0f} mm3  max_u={mx:.3f}")

# Also include prior fa_grad result for reference
dsc_fa, hd_fa, vol_fa, _ = run(pid, "fa_grad")
print(f"  fa_grad  : DSC={dsc_fa:.4f}  HD95={hd_fa:.2f} mm  vol90={vol_fa:.0f} mm3  (prior proxy)")