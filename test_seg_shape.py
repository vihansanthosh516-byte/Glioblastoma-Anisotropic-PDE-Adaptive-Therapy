import numpy as np
import nibabel as nib
import os

DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
patient_id = "UCSF-PDGM-0231"
patient_dir = os.path.join(DATA_DIR, f"{patient_id}_nifti")

# Check segmentation shape
seg_path = os.path.join(patient_dir, f"{patient_id}_tumor_segmentation.nii.gz")
seg = nib.load(seg_path).get_fdata()
print(f"Segmentation shape: {seg.shape}")
print(f"Segmentation unique values: {np.unique(seg)}")