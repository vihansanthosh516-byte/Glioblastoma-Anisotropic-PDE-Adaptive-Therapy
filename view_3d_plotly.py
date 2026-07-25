import numpy as np
import plotly.graph_objects as go

# 1. Load 3D concentration
c = np.load("output/phase3_aniso_concentration.npy")

# 2. Set up spatial grid coordinates (0.78125 mm voxel size)
dx = 0.78125
x, y, z = np.mgrid[0:c.shape[0], 0:c.shape[1], 0:c.shape[2]] * dx

# 3. Create interactive 3D Isosurface
fig = go.Figure(data=go.Isosurface(
    x=x.flatten(),
    y=y.flatten(),
    z=z.flatten(),
    value=c.flatten(),
    isomin=c.max() * 0.15,
    isomax=c.max(),
    surface_count=4,
    colorscale="YlOrRd",
    caps=dict(x_show=False, y_show=False, z_show=False)
))

fig.update_layout(
    title="Phase 3: Interactive 3D Tumor Concentration Isosurface",
    scene=dict(
        xaxis_title="X (mm)",
        yaxis_title="Y (mm)",
        zaxis_title="Z (mm)"
    )
)

fig.write_html("output/3d_tumor_interactive.html")
print("Saved interactive 3D model to output/3d_tumor_interactive.html")

