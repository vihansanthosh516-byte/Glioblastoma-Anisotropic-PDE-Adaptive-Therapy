import nibabel as nib
import numpy as np
import os

DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
patient_id = "UCSF-PDGM-0231"
patient_dir = os.path.join(DATA_DIR, f"{patient_id}_nifti")

# Load patient eigenvalues
l1_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata()
l2_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata()
l3_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L3.nii.gz")).get_fdata()

# Check if they vary
print(f"l1_raw min: {l1_raw.min():.6f}, max: {l1_raw.max():.6f}, mean: {l1_raw.mean():.6f}")
print(f"l2_raw min: {l2_raw.min():.6f}, max: {l2_raw.max():.6f}, mean: {l2_raw.mean():.6f}")
print(f"l3_raw min: {l3_raw.min():.6f}, max: {l3_raw.max():.6f}, mean: {l3_raw.mean():.6f}")
print(f"l1_raw std: {l1_raw.std():.6f}")

# Check atlas scaled values
l1_atlas = l1_raw * 1.5
l2_atlas = l2_raw * 0.5
l3_atlas = l3_raw * 0.5

print(f"\nAtlas scaled - l1 mean: {l1_atlas.mean():.6f}, std: {l1_atlas.std():.6f}")
print(f"Atlas scaled - l2 mean: {l2_atlas.mean():.6f}, std: {l2_atlas.std():.6f}")
print(f"Atlas scaled - l3 mean: {l3_atlas.mean():.6f}, std: {l3_atlas.std():.6f}")

# Check after clipping
l1_clip = np.maximum(l1_atlas, 1e-6)
l2_clip = np.maximum(l2_atlas, 1e-6)
l3_clip = np.maximum(l3_atlas, 1e-6)

print(f"\nAfter clip - l1 mean: {l1_clip.mean():.6f}, std: {l1_clip.std():.6f}")
print(f"After clip - l2 mean: {l2_clip.mean():.6f}, std: {l2_clip.std():.6f}")
print(f"After clip - l3 mean: {l3_clip.mean():.6f}, std: {l3_clip.std():.6f}")