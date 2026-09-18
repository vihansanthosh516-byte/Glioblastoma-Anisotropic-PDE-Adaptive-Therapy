#!/usr/bin/env python3
"""
Option C: 5 grade-2/3 diffuse-glioma patients, anisotropic (true tensor,
Method 3) vs isotropic validation. Same rho/D across patients; only the
patient-specific DTI tensor varies.

Diffuse astrocytomas (grade 2-3) are defined by white-matter tract
infiltration, the regime where anisotropy should matter most.
"""
import sys
import numpy as np
import nibabel as nib
import pandas as pd
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
META = pd.read_csv(r"C:\Users\vihan\Downloads\UCSF-PDGM-metadata_v5.csv")

RHO = 0.0032
D_BASE = 0.013
EXTENT = 0.15
DAYS = 90
TARGET = (48, 74, 74)

# Grade 2-3 IDH-mutant diffuse astrocytomas (available in download, small-to-moderate tumors)
# From metadata these are diffuse/fibrillary low-grade tumors -> tract-infiltration regime
PATIENTS = ["UCSF-PDGM-0231", "UCSF-PDGM-0232", "UCSF-PDGM-0234",
            "UCSF-PDGM-0237", "UCSF-PDGM-0241"]


def fit_tensor_field(patient_dir, pid):
    """Fit true diffusion tensor (method-3) via Stejskal-Tanner log-linear LSQ."""
    dwi = nib.load(str(patient_dir / f"{pid}_DTI_eddy_noreg.nii.gz")).get_fdata()
    bvec = np.array([l.split() for l in open(
        str(patient_dir / f"{pid}_DTI_eddy.eddy_rotated_bvecs")).read().splitlines()], dtype=float)
    bval = np.array([float(x) for x in open(str(BVAL)).read().split()])
    ny, nx, nz, nt = dwi.shape
    b0_idx = np.where(bval == 0)[0][0]
    dwi_idx = np.where(bval > 0)[0]
    S0 = np.maximum(dwi[..., b0_idx], 1e-12)
    bv = bval[dwi_idx]
    g = bvec[:, dwi_idx]
    A = np.stack([g[0]**2, g[1]**2, g[2]**2, 2*g[0]*g[1], 2*g[0]*g[2], 2*g[1]*g[2]], axis=1)
    AtA_inv = np.linalg.inv(A.T @ A)
    X = np.maximum(dwi[..., dwi_idx].astype(np.float64), 1e-12)
    y = -(1.0 / bv[None, None, None, :]) * np.log(X / S0[..., None])

    shape = (ny, nx, nz)
    D = np.zeros(shape + (3, 3))
    # Solve per z-slice to bound memory; PD-enforce each voxel individually.
    for zi in range(nz):
        yz = y[..., zi, :]                       # (ny, nx, ndwi)
        Aty = np.tensordot(A.T, yz, axes=([1], [2]))          # (6, ny, nx)
        Dvec = np.tensordot(AtA_inv, Aty, axes=([1], [0]))    # (6, ny, nx)
        Dvec = np.moveaxis(Dvec, 0, -1)          # (ny, nx, 6)
        Dxx, Dyy, Dzz = Dvec[..., 0], Dvec[..., 1], Dvec[..., 2]
        Dxy, Dxz, Dyz = Dvec[..., 3], Dvec[..., 4], Dvec[..., 5]
        D[..., zi, 0, 0] = Dxx; D[..., zi, 1, 1] = Dyy; D[..., zi, 2, 2] = Dzz
        D[..., zi, 0, 1] = D[..., zi, 1, 0] = Dxy
        D[..., zi, 0, 2] = D[..., zi, 2, 0] = Dxz
        D[..., zi, 1, 2] = D[..., zi, 2, 1] = Dyz
        # PD enforcement per voxel on this slice
        flat = D[..., zi, :, :].reshape(-1, 3, 3)
        eigvals, eigvecs = np.linalg.eigh(flat)          # eigh -> (w, v)
        eigvals = np.maximum(eigvals, 1e-6)
        Dp = np.einsum('nij,nj,nkj->nik', eigvecs, eigvals, eigvecs)  # V diag(e) V^T
        D[..., zi, :, :] = Dp.reshape(ny, nx, 3, 3)
    return D, shape


def resample_tensor(D, src, target):
    zf = np.array(target) / np.array(src)
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


def run(Dt, mask, vox_vol, aniso):
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
    return dsc, hd, vol


rows = []
for pid in PATIENTS:
    patient_dir = UCSF_ROOT / f"{pid}_nifti"
    if not patient_dir.is_dir():
        print(f"SKIP {pid}: not found"); continue
    num = int(pid.split("-")[-1])
    meta_id = f"UCSF-PDGM-{num:03d}"
    mrow = META[META["ID"] == meta_id]
    meta_row = mrow.iloc[0] if not mrow.empty else pd.Series({"WHO CNS Grade": "NA", "Final pathologic diagnosis (WHO 2021)": "NA"})

    seg_nii = nib.load(str(patient_dir / f"{pid}_tumor_segmentation.nii.gz"))
    seg_orig = seg_nii.get_fdata()
    vox_vol = abs(np.linalg.det(seg_nii.affine[:3, :3]))
    zf = np.array(TARGET) / np.array(seg_orig.shape)
    mask = zoom(seg_orig, zf, order=0) > 0.5

    Df, src = fit_tensor_field(patient_dir, pid)
    Dt = resample_tensor(Df, src, TARGET)
    eig = np.linalg.eigvalsh(Dt)
    pd_ok = bool(eig.min() > -1e-12)
    md_t = (Dt[...,0,0]+Dt[...,1,1]+Dt[...,2,2])/3.0
    ok = md_t > 1e-9
    md_mean = md_t[ok].mean() if ok.any() else 1e-6
    Dt = Dt * (D_BASE / md_mean)

    dsc_a, hd_a, vol_a = run(Dt, mask, vox_vol, True)
    dsc_i, hd_i, vol_i = run(Dt, mask, vox_vol, False)
    rows.append({
        "patient": pid, "grade": meta_row["WHO CNS Grade"],
        "diag": meta_row["Final pathologic diagnosis (WHO 2021)"],
        "dsc_aniso": dsc_a, "dsc_iso": dsc_i,
        "hd95_aniso": hd_a, "hd95_iso": hd_i,
        "vol_d90_aniso": vol_a, "vol_d90_iso": vol_i,
        "pd_ok": pd_ok, "tumor_vox": int(mask.sum()),
    })
    print(f"  {pid} done: aniso DSC={dsc_a:.3f} iso DSC={dsc_i:.3f} (PD={pd_ok}, vox={mask.sum()})")

df = pd.DataFrame(rows)
print("\n" + "=" * 100)
print("OPTION C: GRADE 2-3 DIFFUSE GLIOMAS — aniso (true tensor) vs iso")
print("=" * 100)
print(df.to_string(index=False))

dfn = df.dropna(subset=["dsc_aniso", "dsc_iso"])
if len(dfn) >= 5:
    from scipy.stats import wilcoxon
    try:
        stat, p = wilcoxon(dfn["dsc_aniso"], dfn["dsc_iso"])
        print(f"\nAniso vs Iso Wilcoxon signed-rank (n={len(dfn)}): W={stat:.1f}, p={p:.4f}")
    except ValueError as e:
        print(f"\nWilcoxon n/a: {e}")
    print(f"Mean DSC aniso: {dfn['dsc_aniso'].mean():.4f} | Mean DSC iso: {dfn['dsc_iso'].mean():.4f}")
    delta = dfn["dsc_aniso"] - dfn["dsc_iso"]
    print(f"Mean difference (aniso-iso): {delta.mean():+.4f}, aniso wins {int((delta>0).sum())}/{len(delta)} patients")

df.to_csv("output/optionC_grade23_validation.csv", index=False)
print("\nSaved output/optionC_grade23_validation.csv")