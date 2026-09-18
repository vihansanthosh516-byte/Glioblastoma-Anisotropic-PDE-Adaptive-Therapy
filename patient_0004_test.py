import sys
sys.path.insert(0, 'src')
import numpy as np
import nibabel as nib
from pathlib import Path
from scipy.ndimage import zoom
from scipy.ndimage import gaussian_filter
import importlib.util

# Load PDE solver module
spec = importlib.util.spec_from_file_location("aniso", "src/42_anisotropic_pde.py")
aniso_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aniso_mod)
AnisotropicFKSolver = aniso_mod.AnisotropicFKSolver

# Patient directory
patient_dir = r'C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5\UCSF-PDGM-0004_nifti'
patient_id = 'UCSF-PDGM-0004'

# Step 1: Load eigenvalue maps
L1 = nib.load(str(Path(patient_dir) / f"{patient_id}_DTI_eddy_L1.nii.gz")).get_fdata().astype(np.float32)
L2 = nib.load(str(Path(patient_dir) / f"{patient_id}_DTI_eddy_L2.nii.gz")).get_fdata().astype(np.float32)
fa = nib.load(str(Path(patient_dir) / f"{patient_id}_DTI_eddy_FA.nii.gz")).get_fdata().astype(np.float32)

# Extract central axial slice
slice_idx = L1.shape[2] // 2
L1_s = L1[:, :, slice_idx]
L2_s = L2[:, :, slice_idx]
fa_s = fa[:, :, slice_idx]

# Compute principal eigenvector from FA gradient
fa_s_smooth = gaussian_filter(fa_s, sigma=1.0)
gx, gy, gz = np.gradient(fa_s_smooth)
theta = np.arctan2(gy, gx)

# Eigenvalue fields
d_parallel = 0.013
d_perp = 0.0013
d_base = 0.0013
fa_threshold = 0.2
tract_mask = fa_s > fa_threshold
lam1_field = np.where(tract_mask, d_parallel, d_base)
lam2_field = np.where(tract_mask, d_perp, d_base)

# Build tensor components
c = np.cos(theta)
s = np.sin(theta)
D_xx = lam1_field * c * c + lam2_field * s * s
D_yy = lam1_field * s * s + lam2_field * c * c
D_xy = (lam1_field - lam2_field) * s * c

# Resize to 100x100
target = (100, 100)
D_xx = zoom(D_xx, (target[0]/D_xx.shape[0], target[1]/D_xx.shape[1]), order=3)
D_yy = zoom(D_yy, (target[0]/D_yy.shape[0], target[1]/D_yy.shape[1]), order=3)
D_xy = zoom(D_xy, (target[0]/D_xy.shape[0], target[1]/D_xy.shape[1]), order=3)

# Step 2: Run PDE solver
builder = AnisotropicFKSolver(
    D_xx=D_xx,
    D_xy=D_xy,
    D_yy=D_yy,
    dt=0.1,
    rho=0.025,
    carrying_capacity=1.0,
)

print(f"Solver init: grid={builder.H}x{builder.W}")

# Run 90 days
for step in range(90):
    builder.step()

predicted_volume_mm3 = float(builder.u.sum()) * (builder.dx ** 3)
print(f"\n1. Predicted volume at day 90: {predicted_volume_mm3:.2f} mm³")

# Step 2: Actual volume from segmentation
tumor_mask_path = Path(patient_dir) / f"{patient_id}_tumor_segmentation.nii.gz"
tumor_mask_nifti = nib.load(str(tumor_mask_path))
tumor_mask_data = tumor_mask_nifti.get_fdata()

# Voxel volume from affine
affine = tumor_mask_nifti.affine
voxel_volume = abs(np.linalg.det(affine)) ** (1/3)
actual_volume_mm3 = float(tumor_mask_data.sum()) * voxel_volume
print(f"\n2. Actual volume from segmentation: {actual_volume_mm3:.2f} mm³")

# Step 3: DSC between predicted and actual shape
u_field = builder.u
actual_binary = (tumor_mask_data > 0).astype(int)

# Resize actual to match PDE grid (100x100) if needed
if tumor_mask_data.shape != (100, 100):
    actual_binary_resized = zoom(tumor_mask_data.astype(float), (100/tumor_mask_data.shape[0], 100/tumor_mask_data.shape[1])) > 0.5
else:
    actual_binary_resized = (tumor_mask_data > 0).astype(int)

# Try several thresholds on PDE solution
print("\n3. DSC between predicted and actual shape:")
for thresh_factor in [0.1, 0.2, 0.3, 0.4, 0.5]:
    threshold = thresh_factor * builder.u.max()
    pred_binary = (builder.u > threshold).astype(int)
    
    # Make sure both are 100x100
    if pred_binary.shape != (100, 100):
        pred_binary = zoom(pred_binary.astype(float), (100/pred_binary.shape[0], 100/pred_binary.shape[1])) > 0.5
    
    intersection = np.sum(pred_binary * actual_binary_resized)
    union = np.sum(pred_binary) + np.sum(actual_binary_resized)
    dsc = (2.0 * intersection) / union if union > 0 else 0
    print(f"   Threshold {thresh_factor:.1f}*max: DSC={dsc:.4f}")

print("\n" + "="*60)
print("END-TO-END TEST COMPLETE")
PYEOF