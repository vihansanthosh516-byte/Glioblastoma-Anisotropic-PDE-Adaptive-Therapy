import nibabel as nib
import os
import numpy as np

DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
patient_ids = ["UCSF-PDGM-0231", "UCSF-PDGM-0232", "UCSF-PDGM-0233", "UCSF-PDGM-0234", "UCSF-PDGM-0235"]

for pid in patient_ids:
    patient_dir = os.path.join(DATA_DIR, f"{pid}_nifti")
    try:
        fa = nib.load(os.path.join(patient_dir, f"{pid}_DTI_eddy_FA.nii.gz")).get_fdata()
        seg = nib.load(os.path.join(patient_dir, f"{pid}_tumor_segmentation.nii.gz")).get_fdata()
        print(f"{pid}: FA={fa.shape}, Seg={seg.shape}, Match={fa.shape == seg.shape}")
    except Exception as e:
        print(f"{pid}: ERROR - {e}")