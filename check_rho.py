import numpy as np

d = np.load('output/spatial_recurrence_profiles.npz', allow_pickle=True)
rho = d['rho_fields']
density = d['density_maps']
ids = d['patient_ids']

print(f'N patients: {len(ids)}')
print(f'rho_fields shape: {rho.shape}')
print()
print('Per-patient rho_field stats:')
for i in range(5):
    r = rho[i]
    dm = density[i]
    print(f'  {ids[i]}:')
    print(f'    rho_field:  min={r.min():.2e}, max={r.max():.2e}, mean={r.mean():.2e}, nonzero={np.count_nonzero(r)}')
    print(f'    density:    min={dm.min():.2e}, max={dm.max():.2e}, nonzero={np.count_nonzero(dm)}')
    print(f'    rho at tumor center: {r[50, 50]:.4f}')
    print(f'    rho outside tumor:   {r[10, 10]:.4f}')