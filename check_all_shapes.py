import nibabel as nib
import os
import numpy as np

DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
META_PATH = os.path.join(DATA_DIR, "UCSF-PDGM-metadata_v5.csv")
import pandas as pd
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

# Check shapes for all
shapes = {}
for pid in patient_ids:
    patient_dir = os.path.join(DATA_DIR, f"{pid}_nifti")
    try:
        fa = nib.load(os.path.join(patient_dir, f"{pid}_DTI_eddy_FA.nii.gz")).get_fdata()
        seg = nib.load(os.path.join(patient_dir, f"{pid}_tumor_segmentation.nii.gz")).get_fdata()
        l1 = nib.load(os.path.join(patient_dir, f"{pid}_DTI_eddy_L1.nii.gz")).get_fdata()
        key = (fa.shape, seg.shape, l1.shape)
        if key not in shapes:
            shapes[key] = []
        shapes[key].append(pid)
    except Exception as e:
        print(f"{pid}: ERROR - {e}")

print("\nShape groups:")
for shape, pids in shapes.items():
    print(f"  {shape}: {len(pids)} patients - {pids[:5]}{'...' if len(pids) > 5 else ''}")