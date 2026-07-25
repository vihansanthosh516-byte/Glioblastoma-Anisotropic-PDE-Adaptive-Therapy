import numpy as np
import matplotlib.pyplot as plt

# Load 3D anisotropic concentration array
c = np.load("output/phase3_aniso_concentration.npy")

print(f"Loaded 3D concentration shape: {c.shape}")

# Define a threshold for tumor volume (e.g., top 20% concentration)
threshold = c.max() * 0.2
voxels = c > threshold

# Downsample grid slightly for fast 3D interactive rendering
ds = 2  # Downsample factor
voxels_ds = voxels[::ds, ::ds, ::ds]

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection="3d")

# Render 3D voxels
ax.voxels(voxels_ds, facecolors="#d62728", edgecolor="k", alpha=0.6)

ax.set_title("Phase 3: Interactive 3D Tumor Concentration Volume")
ax.set_xlabel("X (mm)")
ax.set_ylabel("Y (mm)")
ax.set_zlabel("Z (mm)")

plt.tight_layout()
plt.show()

