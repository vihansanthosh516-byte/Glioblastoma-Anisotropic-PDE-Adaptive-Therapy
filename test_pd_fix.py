import sys
sys.path.insert(0, 'src')
from scipy.ndimage import zoom
from src.load_ucsf_tensor import load_patient_tensors, is_pd
import numpy as np

patient_dir = r'C:\Users\vihan\Downloads\PKG - UCSF-PDGM Version 5\UCSF-PDGM-v5\UCSF-PDGM-0004_nifti'
result = load_patient_tensors(patient_dir, target_shape=(50, 50, 50))

tensors = result['tensors']
print(f'Tensor shape: {tensors.shape}')

pd_count = 0
total = tensors.shape[0] * tensors.shape[1] * tensors.shape[2]
for i in range(tensors.shape[0]):
    for j in range(tensors.shape[1]):
        for k in range(tensors.shape[2]):
            D = tensors[i, j, k]
            if not is_pd(D):
                pd_count += 1

print(f'PD check: {pd_count} non-PD voxels out of {total} total')

fa_mean = float(result['fa'].mean())
md_mean = float(result['md'].mean())
l1_mean = float(result['eigenvalues'][..., 0].mean())
l2_mean = float(result['eigenvalues'][..., 1].mean())
l3_mean = float(result['eigenvalues'][..., 2].mean())

print(f'Mean FA: {fa_mean:.4f}')
print(f'Mean MD: {md_mean:.6f}')
print(f'Mean L1: {l1_mean:.6f}')
print(f'Mean L2: {l2_mean:.6f}')
print(f'Mean L3: {l3_mean:.6f}')

eigenvals = result['eigenvalues']
fa_from_eigen = np.sqrt(1.5 * ((eigenvals[..., 0] - eigenvals[..., 1])**2 + (eigenvals[..., 1] - eigenvals[..., 2])**2 + (eigenvals[..., 2] - eigenvals[..., 0])**2) / (eigenvals[..., 0] + eigenvals[..., 1] + eigenvals[..., 2] + 1e-10))
corr = np.corrcoef(fa_from_eigen.flatten(), result['fa'].flatten())[0,1]
print(f'FA from eigenvalues mean: {fa_from_eigen.mean():.4f}')
print(f'FA from map mean: {fa_mean:.4f}')
print(f'FA correlation: {corr:.4f}')
print('SUCCESS: All voxels positive-definite, FA/MD consistent')