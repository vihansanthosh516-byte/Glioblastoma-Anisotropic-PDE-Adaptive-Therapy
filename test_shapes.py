import numpy as np
import nibabel as nib
import os
from scipy.ndimage import zoom

DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
patient_id = "UCSF-PDGM-0231"
patient_dir = os.path.join(DATA_DIR, f"{patient_id}_nifti")

# Load all arrays
fa_standard = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata()
md_standard = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_MD.nii.gz")).get_fdata()
l1_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata()
l2_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata()
l3_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L3.nii.gz")).get_fdata()

print(f"FA shape: {fa_standard.shape}")
print(f"MD shape: {md_standard.shape}")
print(f"L1 shape: {l1_raw.shape}")
print(f"L2 shape: {l2_raw.shape}")
print(f"L3 shape: {l3_raw.shape}")

# Test resampling
target_shape = fa_standard.shape
for name, arr in [('L1', l1_raw), ('L2', l2_raw), ('L3', l3_raw)]:
    if arr.shape != target_shape:
        print(f"Resampling {name} from {arr.shape} to {target_shape}")
        zoom_factors = [t / s for t, s in zip(target_shape, arr.shape)]
        resampled = zoom(arr, zoom_factors, order=1)
        print(f"  Resampled shape: {resampled.shape}")
    else:
        print(f"{name} already matches target shape")

print("Shape verification complete!")