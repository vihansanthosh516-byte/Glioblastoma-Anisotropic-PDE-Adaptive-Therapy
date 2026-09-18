import numpy as np
import nibabel as nib
import os
from scipy.ndimage import laplace, zoom

DATA_DIR = "C:\\Users\\vihan\\Downloads\\ucsf-pdgm-grade23"
patient_id = "UCSF-PDGM-0231"
patient_dir = os.path.join(DATA_DIR, f"{patient_id}_nifti")

# Load standard DTI first to get the correct shape
fa_standard = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata()
md_standard = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_MD.nii.gz")).get_fdata()
fa_shape = fa_standard.shape

print(f"FA shape: {fa_shape}")

# DKI mode
fa = fa_standard * 0.6
md = md_standard * 0.0015
l1_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata()
l2_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata()
l3_raw = nib.load(os.path.join(patient_dir, f"{patient_id}_DTI_eddy_L3.nii.gz")).get_fdata()

print(f"L1 raw shape: {l1_raw.shape}")
print(f"L2 raw shape: {l2_raw.shape}")
print(f"L3 raw shape: {l3_raw.shape}")

# Apply kurtosis correction
l1 = l1_raw * 1.3
l2 = l2_raw * 0.8
l3 = l3_raw * 0.8

print(f"L1 shape after scaling: {l1.shape}")

# Build tensor
grad_l1 = np.gradient(l1)
print(f"grad_l1[0] shape: {grad_l1[0].shape}")

v1 = np.stack(grad_l1, axis=-1)
print(f"v1 shape: {v1.shape}")

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
print(f"D_tensor shape: {D_tensor.shape}")

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
print(f"mean_md: {mean_md}")
if mean_md > 0:
    D_tensor *= (D_SCALE / mean_md)

print(f"D_tensor final shape: {D_tensor.shape}")

# Downsample
D_tensor_ds = D_tensor[::4, ::4, ::4]
print(f"D_tensor downsampled: {D_tensor_ds.shape}")

# Load and downsample seg
seg = nib.load(os.path.join(patient_dir, f"{patient_id}_tumor_segmentation.nii.gz")).get_fdata()
seg_ds = seg[::4, ::4, ::4]
print(f"Seg downsampled: {seg_ds.shape}")

u0 = (seg_ds > 0).astype(float) * 0.9
print(f"u0 shape: {u0.shape}")

# Test PDE step
print("Testing PDE step...")
u = u0.copy()
rho = 0.0032
dt = 0.1

for step in range(5):  # Just 5 steps
    growth = rho * u * (1 - u)
    D_mean = (D_tensor_ds[..., 0, 0] + D_tensor_ds[..., 1, 1] + D_tensor_ds[..., 2, 2]) / 3.0
    l1_ds = D_tensor_ds[..., 0, 0]
    l2_ds = D_tensor_ds[..., 1, 1]
    spread = np.abs(l1_ds - l2_ds) / (np.abs(l1_ds + l2_ds) + 1e-10)
    diffusion = laplace(u) * D_mean * (1 + 0.5 * spread)
    u += dt * (growth + diffusion)
    u = np.clip(u, 0, 1.0)
    if step == 0:
        print(f"  step 0: u shape={u.shape}, D_mean shape={D_mean.shape}, spread shape={spread.shape}")

print("PDE test passed!")