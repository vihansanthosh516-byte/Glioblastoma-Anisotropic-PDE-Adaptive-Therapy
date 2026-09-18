#!/usr/bin/env python3
"""
5-patient anisotropic vs isotropic validation.
Same rho and D for all patients; only the patient-specific DTI tensor varies.

rho    = 0.0032 /day  (calibrated from MU-Glioma-Post, 90 growing pairs)
D_base = 0.013 mm2/day (Swanson-anchored mean tumor diffusivity)

For each patient, run anisotropic (patient DTI tensor) and isotropic
(D = D_mean * I) forward PDE for 90 days. Report extent metrics at the
documented infiltration threshold, plus MGMT/EOR/OS from clinical metadata.

NOTE: patient UCSF-PDGM-0006 is NOT present in the downloaded data; the
requested slot missing-0006 is filled by UCSF-PDGM-0008 (available).
"""
import sys
import numpy as np
import nibabel as nib
import pandas as pd
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
META = pd.read_csv(r"C:\Users\vihan\Downloads\UCSF-PDGM-metadata_v5.csv")

RHO = 0.0032
D_BASE = 0.013     # target mean diffusivity mm2/day
EXTENT = 0.15      # documented infiltration-extent threshold
DAYS = 90
TARGET = (48, 74, 74)

# Requested set; 0006 absent -> substituted with 0008 (flagged below).
PATIENTS = ["UCSF-PDGM-0004", "UCSF-PDGM-0005", "UCSF-PDGM-0007",
            "UCSF-PDGM-0008", "UCSF-PDGM-0010"]
SUBSTITUTION = {"UCSF-PDGM-0008": "0006 not in download -> 0008 used"}


def extent_metrics(u, mask, vox_vol, threshold=EXTENT):
    pred = u > threshold
    vol = float(np.sum(pred)) * vox_vol
    inter = np.sum(pred & mask)
    dsc = (2 * inter) / (np.sum(pred) + np.sum(mask)) if (np.sum(pred) + np.sum(mask)) > 0 else 0.0
    nb = pred & ~binary_erosion(pred)
    tb = mask & ~binary_erosion(mask)
    d_p2t = distance_transform_edt(~mask)
    d_t2p = distance_transform_edt(~pred)
    h95p = np.percentile(d_p2t[nb], 95) if nb.any() else 0.0
    h95t = np.percentile(d_t2p[tb], 95) if tb.any() else 0.0
    return vol, dsc, float(max(h95p, h95t))


def run_patient(patient_dir, meta_row, aniso):
    pid = patient_dir.name.replace("_nifti", "")
    # Segmentation
    seg_nii = nib.load(str(patient_dir / f"{pid}_tumor_segmentation.nii.gz"))
    seg_orig = seg_nii.get_fdata()
    vox_vol = abs(np.linalg.det(seg_nii.affine[:3, :3]))
    zf = np.array(TARGET) / np.array(seg_orig.shape)
    mask = zoom(seg_orig, zf, order=0) > 0.5

    # Tensor
    out = load_patient_tensors(patient_dir, target_shape=TARGET)
    D3 = out['tensors']  # (48,74,74,3,3) mm2/s
    md_mean = out['md'].mean()
    scale = D_BASE / max(md_mean, 1e-9)  # rescale to D_BASE mm2/day
    D3 = D3 * scale

    if aniso:
        D_xx, D_yy, D_zz = D3[..., 0, 0], D3[..., 1, 1], D3[..., 2, 2]
        D_xy, D_xz, D_yz = D3[..., 0, 1], D3[..., 0, 2], D3[..., 1, 2]
    else:
        D_mean = (D3[..., 0, 0] + D3[..., 1, 1] + D3[..., 2, 2]) / 3.0
        D_xx = D_yy = D_zz = D_mean
        D_xy = D_xz = D_yz = np.zeros_like(D_mean)

    solver = AnisotropicFKSolver3D(
        D_xx=D_xx.astype(np.float64), D_xy=D_xy.astype(np.float64), D_xz=D_xz.astype(np.float64),
        D_yy=D_yy.astype(np.float64), D_yz=D_yz.astype(np.float64), D_zz=D_zz.astype(np.float64),
        dt=0.1, dx=1.0, rho=RHO, K=1.0)
    u0 = np.zeros(TARGET, np.float64)
    u0[mask] = 0.9
    v0, _, _ = extent_metrics(u0, mask, vox_vol)

    u = u0.copy()
    C = 0.0
    steps = int(round(DAYS / solver.dt))
    for _ in range(steps):
        u = solver.step(u, C)

    v90, dsc90, hd9590 = extent_metrics(u, mask, vox_vol)
    return {
        "day0_vol": v0, "day90_vol": v90, "dsc": dsc90, "hd95": hd9590,
        "mgmt": meta_row["MGMT status"], "eor": meta_row["EOR"],
        "os": meta_row["OS"],
    }


def meta_for(pid):
    # Map 4-digit folder id -> 3-digit metadata id (metadata uses 3-digit zero-padded)
    num = int(pid.replace("UCSF-PDGM-", ""))
    return f"UCSF-PDGM-{num:03d}"


rows = []
for pid in PATIENTS:
    patient_dir = UCSF_ROOT / f"{pid}_nifti"
    if not patient_dir.is_dir():
        print(f"SKIP {pid}: not found")
        continue
    meta_id = meta_for(pid)
    mrow = META[META["ID"] == meta_id]
    meta_row = mrow.iloc[0] if not mrow.empty else pd.Series({"MGMT status": "NA", "EOR": "NA", "OS": None})

    an = run_patient(patient_dir, meta_row, aniso=True)
    iso = run_patient(patient_dir, meta_row, aniso=False)
    rows.append({
        "patient": pid, "subst": SUBSTITUTION.get(pid, ""),
        "dsc_aniso": an["dsc"], "dsc_iso": iso["dsc"],
        "hd95_aniso": an["hd95"], "hd95_iso": iso["hd95"],
        "vol_d0": an["day0_vol"], "vol_d90_aniso": an["day90_vol"], "vol_d90_iso": iso["day90_vol"],
        "mgmt": an["mgmt"], "eor": an["eor"], "os": an["os"],
    })

df = pd.DataFrame(rows)
print("\n" + "=" * 100)
print("5-PATIENT VALIDATION  (rho=0.0032, D scaled to 0.013 mm2/day, extent>0.15)")
print("NOTE: 0006 absent in download; substituted with 0008")
print("=" * 100)
print(df.to_string(index=False))

dfn = df[df["dsc_aniso"].notna() & df["dsc_iso"].notna()]
if len(dfn) >= 5:
    from scipy.stats import wilcoxon
    try:
        stat, p = wilcoxon(dfn["dsc_aniso"], dfn["dsc_iso"])
        print(f"\nAniso vs Iso (Wilcoxon signed-rank, n={len(dfn)}): W={stat:.1f}, p={p:.4f}")
    except ValueError as e:
        print(f"\nWilcoxon not applicable: {e}")
    print(f"Mean DSC aniso: {dfn['dsc_aniso'].mean():.3f}  (median {dfn['dsc_aniso'].median():.3f})")
    print(f"Mean DSC iso   : {dfn['dsc_iso'].mean():.3f}   (median {dfn['dsc_iso'].median():.3f})")
    if dfn["mgmt"].nunique() >= 2:
        for g, grp in dfn.groupby("mgmt"):
            print(f"\nMGMT={g} (n={len(grp)}): aniso DSC {grp['dsc_aniso'].mean():.3f} | iso {grp['dsc_iso'].mean():.3f}")

df.to_csv("output/five_patient_validation.csv", index=False)
print("\nSaved output/five_patient_validation.csv")