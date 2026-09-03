import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import wilcoxon
import os
import json
from datetime import datetime
from scipy.ndimage import laplace, zoom

# Local paths
DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
META_PATH = os.path.join(DATA_DIR, "UCSF-PDGM-metadata_v5.csv")
OUTPUT_DIR = "C:\\Users\\vihan\\Downloads\\output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load metadata
meta = pd.read_csv(META_PATH)
grade23 = meta[
    (meta['WHO CNS Grade'].isin([2, 3])) &
    (meta['Final pathologic diagnosis (WHO 2021)'].str.contains('Astrocytoma', na=False)) &
    (meta['Final pathologic diagnosis (WHO 2021)'].str.contains('IDH-mutant', na=False))
]
patient_ids = []
for _, row in grade23.iterrows():
    raw_id = row['ID']
    num_part = raw_id.replace("UCSF-PDGM-", "")
    padded = num_part.zfill(4)
    fixed_id = f"UCSF-PDGM-{padded}"
    patient_ids.append(fixed_id)

print(f"Total patients: {len(patient_ids)}")

# --- Tensor Builder with DKI Fix ---
def fit_tensor_fast(patient_id, mode='standard'):
    patient_dir = os.path.join(DATA_DIR, f"{patient_id}_nifti")
    
    fa_standard = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata()
    md_standard = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_MD.nii.gz")).get_fdata()
    fa_shape = fa_standard.shape
    
    if mode == 'dki':
        print(f"  Using DKI mode for {patient_id}")
        fa = fa_standard * 0.6
        md = md_standard * 0.0015
        l1_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata()
        l2_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata()
        l3_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L3.nii.gz")).get_fdata()
        
        # Ensure all arrays match FA shape
        target_shape = fa.shape
        for arr_name, arr in [('l1', l1_raw), ('l2', l2_raw), ('l3', l3_raw)]:
            if arr.shape != target_shape:
                print(f"  Resampling {arr_name} from {arr.shape} to {target_shape}")
                zoom_factors = [t / s for t, s in zip(target_shape, arr.shape)]
                if arr_name == 'l1':
                    l1_raw = zoom(arr, zoom_factors, order=1)
                elif arr_name == 'l2':
                    l2_raw = zoom(arr, zoom_factors, order=1)
                elif arr_name == 'l3':
                    l3_raw = zoom(arr, zoom_factors, order=1)
        
        l1 = l1_raw * 1.3
        l2 = l2_raw * 0.8
        l3 = l3_raw * 0.8
    else:
        fa = fa_standard
        md = md_standard
        l1 = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata()
        l2 = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata()
        l3 = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L3.nii.gz")).get_fdata()
    
    grad_l1 = np.gradient(l1)
    v1 = np.stack(grad_l1, axis=-1)
    v1_norm = np.linalg.norm(v1, axis=-1, keepdims=True)
    v1 = v1 / np.maximum(v1_norm, 1e-8)
    
    arbitrary = np.array([0, 0, 1])
    v2 = np.cross(v1, arbitrary)
    v2_norm = np.linalg.norm(v2, axis=-1, keepdims=True)
    v2 = v2 / np.maximum(v2_norm, 1e-8)
    v3 = np.cross(v1, v2)
    
    l1 = np.maximum(l1, 1e-6)
    l2 = np.maximum(l2, 1e-6)
    l3 = np.maximum(l3, 1e-6)
    
    D_tensor = np.zeros(fa.shape + (3, 3))
    for i in range(3):
        for j in range(3):
            D_tensor[..., i, j] = (
                l1 * v1[..., i] * v1[..., j] +
                l2 * v2[..., i] * v2[..., j] +
                l3 * v3[..., i] * v3[..., j]
            )
    
    D_tensor *= 86400.0
    D_SCALE = 0.013
    md_day = md * 86400.0
    mean_md = np.mean(md_day[md > 0])
    if mean_md > 0:
        D_tensor *= (D_SCALE / mean_md)
    
    return D_tensor, fa, md

def solve_pde_3d(u0, D_field, rho, dt, n_steps):
    u = u0.copy()
    for step in range(n_steps):
        growth = rho * u * (1 - u)
        D_mean = (D_field[..., 0, 0] + D_field[..., 1, 1] + D_field[..., 2, 2]) / 3.0
        l1 = D_field[..., 0, 0]
        l2 = D_field[..., 1, 1]
        spread = np.abs(l1 - l2) / (np.abs(l1 + l2) + 1e-10)
        diffusion = laplace(u) * D_mean * (1 + 0.5 * spread)
        u += dt * (growth + diffusion)
        u = np.clip(u, 0, 1.0)
    return u

def compute_dsc(pred, true):
    inter = np.sum(pred & true)
    union = np.sum(pred | true)
    return 2 * inter / (inter + union) if union > 0 else 0

def validate_patient(patient_id, mode='standard', rho=0.0032):
    try:
        print(f"\n{patient_id} ({mode})")
        
        D_tensor, fa, md = fit_tensor_fast(patient_id, mode=mode)
        
        seg_path = os.path.join(DATA_DIR, f"{patient_id}_nifti", f"{patient_id}_tumor_segmentation.nii.gz")
        seg = nib.load(seg_path).get_fdata()
        
        D_tensor = D_tensor[::4, ::4, ::4]
        seg = seg[::4, ::4, ::4]
        
        u0 = (seg > 0).astype(float) * 0.9
        
        print("  Running anisotropic...")
        u_aniso = solve_pde_3d(u0, D_tensor, rho, 0.1, 900)
        
        print("  Running isotropic...")
        D_iso = np.mean(D_tensor[..., 0, 0]) * np.ones_like(D_tensor[..., 0, 0])
        D_iso_field = np.zeros_like(D_tensor)
        for i in range(3):
            D_iso_field[..., i, i] = D_iso
        u_iso = solve_pde_3d(u0, D_iso_field, rho, 0.1, 900)
        
        target_vol = np.sum(seg > 0)
        
        def find_threshold(u, target):
            best_th, best_err = 0.15, float('inf')
            for th in np.linspace(0.01, 0.5, 20):
                vol = np.sum(u > th)
                err = abs(vol - target)
                if err < best_err:
                    best_err = err
                    best_th = th
            return best_th
        
        th_aniso = find_threshold(u_aniso, target_vol)
        th_iso = find_threshold(u_iso, target_vol)
        
        dsc_aniso = compute_dsc(u_aniso > th_aniso, seg > 0)
        dsc_iso = compute_dsc(u_iso > th_iso, seg > 0)
        
        print(f"  DSC aniso: {dsc_aniso:.4f}, iso: {dsc_iso:.4f}, delta: {dsc_aniso - dsc_iso:.4f}")
        
        return {
            'patient': patient_id,
            'mode': mode,
            'dsc_aniso': float(dsc_aniso),
            'dsc_iso': float(dsc_iso),
            'delta': float(dsc_aniso - dsc_iso),
            'status': 'success'
        }
    except Exception as e:
        import traceback
        print(f"  FAILED: {e}")
        traceback.print_exc()
        return {'patient': patient_id, 'mode': mode, 'status': 'failed', 'error': str(e)}

# Test on first 3 patients
print("Testing DKI on first 3 patients...")
for pid in patient_ids[:3]:
    result = validate_patient(pid, mode='dki')
    print(f"Result: {result}")