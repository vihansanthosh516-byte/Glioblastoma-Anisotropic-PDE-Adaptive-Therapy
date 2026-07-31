import os
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt

# Path to the real patient segmentation file
nii_path = "data/brats/BraTS2021_00000/BraTS2021_00000_seg.nii.gz"

if not os.path.exists(nii_path):
    import glob
    found = glob.glob("data/brats/*/*seg.nii.gz")
    if found:
        nii_path = found[0]

# Load 3D NIfTI volume
nii = nib.load(nii_path)
tumor_volume = nii.get_fdata()

print(f"Loaded Real Patient Volume Shape: {tumor_volume.shape}")

# Extract binary mask for whole tumor region
voxels = tumor_volume > 0

# Downsample grid by factor of 2 for fast 3D rendering
ds = 2
voxels_ds = voxels[::ds, ::ds, ::ds]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Render 3D voxels
ax.voxels(voxels_ds, facecolors="crimson", edgecolor="k", alpha=0.5)

ax.set_title("Real Patient 3D Tumor Segmentation (BraTS2021_00000)")
ax.set_xlabel("X (voxels)")
ax.set_ylabel("Y (voxels)")
ax.set_zlabel("Z (voxels)")

plt.tight_layout()
plt.show()