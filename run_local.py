import numpy as np
import pandas as pd
import nibabel as nib
from scipy.stats import wilcoxon
import os
import json
import time
from datetime import datetime
from scipy.ndimage import laplace

# Local paths
DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
META_PATH = os.path.join(DATA_DIR, "UCSF-PDGM-metadata_v5.csv")
OUTPUT_DIR = "C:\\Users\\vihan\\Downloads\\output"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("Data loaded.")

# Load metadata
meta = pd.read_csv(META_PATH)

# Filter for grade 2-3 IDH-mutant astrocytomas
grade23 = meta[
    (meta['WHO CNS Grade'].isin([2, 3])) &
    (meta['Final pathologic diagnosis (WHO 2021)'].str.contains('Astrocytoma', na=False)) &
    (meta['Final pathologic diagnosis (WHO 2021)'].str.contains('IDH-mutant', na=False))
]

# Fix patient IDs by adding leading zeros
patient_ids = []
for _, row in grade23.iterrows():
    raw_id = row['ID']
    num_part = raw_id.replace("UCSF-PDGM-", "")
    padded = num_part.zfill(4)
    fixed_id = f"UCSF-PDGM-{padded}"
    patient_ids.append(fixed_id)

print(f"Found {len(patient_ids)} grade 2-3 patients")
print(f"First 5 (fixed): {patient_ids[:5]}")

def fit_tensor_fast(patient_id):
    patient_dir = os.path.join(DATA_DIR, f"{patient_id}_nifti")
    
    fa = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata()
    md = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_MD.nii.gz")).get_fdata()
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
        D_trace = (D_field[..., 0, 0] + D_field[..., 1, 1] + D_field[..., 2, 2]) / 3.0
        diffusion = laplace(u) * D_trace
        u += dt * (growth + diffusion)
        u = np.clip(u, 0, 1.0)
    return u

def compute_dsc(pred, true):
    inter = np.sum(pred & true)
    union = np.sum(pred | true)
    return 2 * inter / (inter + union) if union > 0 else 0

def validate_patient(patient_id, rho=0.0032):
    try:
        print(f"\n{patient_id}")
        
        D_tensor, fa, md = fit_tensor_fast(patient_id)
        
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
            'dsc_aniso': float(dsc_aniso),
            'dsc_iso': float(dsc_iso),
            'delta': float(dsc_aniso - dsc_iso),
            'status': 'success'
        }
    except Exception as e:
        print(f"  FAILED: {e}")
        return {'patient': patient_id, 'status': 'failed', 'error': str(e)}

def run_batch(patient_ids):
    results = []
    checkpoint_path = os.path.join(OUTPUT_DIR, "checkpoint.json")
    
    if os.path.exists(checkpoint_path):
        with open(checkpoint_path) as f:
            checkpoint = json.load(f)
        done = set(checkpoint['done'])
        results = checkpoint['results']
        print(f"Resuming: {len(done)}/{len(patient_ids)} done")
    else:
        done = set()
        print(f"Starting fresh: {len(patient_ids)} patients")
    
    remaining = [p for p in patient_ids if p not in done]
    
    for i, pid in enumerate(remaining):
        print(f"\n[{i+1}/{len(remaining)}] {pid}")
        result = validate_patient(pid)
        results.append(result)
        
        if result['status'] == 'success':
            done.add(pid)
        
        with open(checkpoint_path, 'w') as f:
            json.dump({'done': list(done), 'results': results}, f)
        
        df = pd.DataFrame([r for r in results if r['status'] == 'success'])
        if len(df) > 0:
            df.to_csv(os.path.join(OUTPUT_DIR, "results.csv"), index=False)
            if len(df) > 1:
                try:
                    w, p = wilcoxon(df['dsc_aniso'], df['dsc_iso'])
                    print(f"  N={len(df)}, aniso={df['dsc_aniso'].mean():.4f}, iso={df['dsc_iso'].mean():.4f}, p={p:.4f}")
                except:
                    pass
    
    return results

# Run batch
results = run_batch(patient_ids)

# Final stats
df = pd.DataFrame([r for r in results if r['status'] == 'success'])

if len(df) > 0:
    w, p = wilcoxon(df['dsc_aniso'], df['dsc_iso'])
    print(f"\n{'='*60}")
    print(f"FINAL RESULTS")
    print(f"{'='*60}")
    print(f"N:           {len(df)}")
    print(f"Aniso mean:  {df['dsc_aniso'].mean():.4f}")
    print(f"Iso mean:    {df['dsc_iso'].mean():.4f}")
    print(f"Delta mean:  {df['delta'].mean():.4f}")
    print(f"Aniso wins:  {sum(df['delta'] > 0)}/{len(df)}")
    print(f"Wilcoxon:    W={w:.1f}, p={p:.4f}")
    print(f"{'='*60}")
    
    df.to_csv(os.path.join(OUTPUT_DIR, f"final_results_{datetime.now().strftime('%Y%m%d_%H%M')}.csv"), index=False)

print(f"\n✅ Results saved to: {OUTPUT_DIR}")