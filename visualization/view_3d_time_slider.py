#!/usr/bin/env python3
"""
Interactive 4D Tumor Evolution Visualizer
=========================================
Loads sequential 3D tumor snapshots and creates a Plotly HTML dashboard
with a time slider to scrub through tumor evolution day-by-day.
"""

import os
import glob
import numpy as np
import plotly.graph_objects as go


def load_time_series(time_series_dir: str = "output/time_series"):
    """Load all time-step .npy files from the time series directory."""
    files = sorted(glob.glob(os.path.join(time_series_dir, "tumor_3d_day_*.npy")))
    if not files:
        raise FileNotFoundError(f"No tumor_3d_day_*.npy files found in {time_series_dir}")
    
    grids = []
    day_labels = []
    for f_path in files:
        grid = np.load(f_path)
        grids.append(grid)
        day_num = os.path.basename(f_path).replace("tumor_3d_day_", "").replace(".npy", "")
        day_labels.append(int(day_num))
    
    return grids, day_labels, files


def create_interactive_dashboard(grids, day_labels, output_path: str = "output/true_3d_time_series_dashboard.html"):
    """Create Plotly figure with 3D isosurface rendering and time slider."""
    
    if not grids:
        print("[ERROR] No grids to visualize")
        return
    
    initial_grid = grids[0]
    nx, ny, nz = initial_grid.shape
    
    # Create spatial coordinate grids
    X, Y, Z = np.mgrid[:nx, :ny, :nz]
    
    # Downsample for lighter Plotly rendering (factor of 2)
    ds = 2
    X = X[::ds, ::ds, ::ds]
    Y = Y[::ds, ::ds, ::ds]
    Z = Z[::ds, ::ds, ::ds]
    
    # Downsample all grids
    grids_ds = [g[::ds, ::ds, ::ds] for g in grids]
    grids = grids_ds
    nx, ny, nz = grids[0].shape
    
    x_flat = X.flatten()
    y_flat = Y.flatten()
    z_flat = Z.flatten()
    
    # Determine value range for consistent isosurfaces across frames
    all_nonzero = np.concatenate([g[g > 0].flatten() for g in grids if np.any(g > 0)])
    if len(all_nonzero) > 0:
        vmin = float(np.percentile(all_nonzero, 10))
        vmax = float(np.percentile(all_nonzero, 95))
    else:
        vmin, vmax = 0.05, 0.8
    
    # Use fewer isosurface levels for lighter HTML
    iso_levels = np.linspace(vmin, vmax, 4)
    
    print(f"[INFO] Grid shape: {initial_grid.shape}")
    print(f"[INFO] Value range: [{vmin:.4f}, {vmax:.4f}]")
    print(f"[INFO] Time points: {day_labels}")
    print(f"[INFO] Isosurface levels: {len(iso_levels)}")
    
    # Create initial frame with multiple isosurfaces
    fig = go.Figure()
    
    for i, level in enumerate(iso_levels):
        fig.add_trace(go.Isosurface(
            x=x_flat,
            y=y_flat,
            z=z_flat,
            value=grids[0].flatten(),
            isomin=level,
            isomax=level,
            surface_count=1,
            colorscale='Hot',
            opacity=0.3,
            caps=dict(x_show=False, y_show=False, z_show=False),
            showscale=bool(i == 0),
            colorbar=dict(title="Tumor Density", thickness=20, len=0.75) if i == 0 else None,
            visible=True
        ))
    
    # Build frames for time slider
    frames = []
    for grid, day in zip(grids, day_labels):
        frame_data = []
        for i, level in enumerate(iso_levels):
            frame_data.append(go.Isosurface(
                x=x_flat,
                y=y_flat,
                z=z_flat,
                value=grid.flatten(),
                isomin=level,
                isomax=level,
                surface_count=1,
                colorscale='Hot',
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showscale=bool(i == 0)
            ))
        frames.append(go.Frame(data=frame_data, name=f"Day {day}"))
    
    fig.frames = frames
    
    # Layout with time slider and play/pause buttons
    fig.update_layout(
        title={
            'text': f"4D Tumor Evolution: Interactive 3D Time-Series (Days {day_labels[0]}–{day_labels[-1]})",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 16}
        },
        scene=dict(
            xaxis_title="X (voxels)",
            yaxis_title="Y (voxels)",
            zaxis_title="Z (voxels)",
            aspectmode='data',
            camera=dict(eye=dict(x=1.5, y=1.5, z=1.2))
        ),
        updatemenus=[{
            "type": "buttons",
            "direction": "left",
            "pad": {"r": 10, "t": 87},
            "showactive": False,
            "x": 0.1,
            "xanchor": "right",
            "y": 0,
            "yanchor": "top",
            "buttons": [
                {
                    "args": [None, {"frame": {"duration": 500, "redraw": True}, "fromcurrent": True, "transition": {"duration": 200, "easing": "cubic-in-out"}, "mode": "immediate"}],
                    "label": "▶ Play",
                    "method": "animate"
                },
                {
                    "args": [[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate", "transition": {"duration": 0}}],
                    "label": "⏸ Pause",
                    "method": "animate"
                }
            ]
        }],
        sliders=[{
            "active": 0,
            "yanchor": "top",
            "xanchor": "left",
            "currentvalue": {
                "font": {"size": 16},
                "prefix": "Simulation Day: ",
                "visible": True,
                "xanchor": "right"
            },
            "transition": {"duration": 200, "easing": "cubic-in-out"},
            "pad": {"b": 10, "t": 50},
            "len": 0.9,
            "x": 0.1,
            "y": 0,
            "steps": [
                {
                    "args": [
                        [f"Day {day}"],
                        {"frame": {"duration": 0, "redraw": True}, "mode": "immediate", "transition": {"duration": 0}}
                    ],
                    "label": f"Day {day}",
                    "method": "animate"
                }
                for day in day_labels
            ]
        }],
        width=1000,
        height=700,
        margin=dict(l=0, r=0, t=60, b=0)
    )
    
    # Save HTML
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_html(output_path, include_plotlyjs='cdn')
    print(f"[SUCCESS] Interactive 4D dashboard saved to: {output_path}")
    
    return output_path


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Generate interactive 4D tumor evolution dashboard from saved time-series snapshots"
    )
    parser.add_argument("--input-dir", type=str, default="output/time_series",
                        help="Directory containing tumor_3d_day_*.npy files")
    parser.add_argument("--output", type=str, default="output/true_3d_time_series_dashboard.html",
                        help="Output HTML file path")
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("  4D INTERACTIVE TUMOR EVOLUTION VISUALIZER")
    print("=" * 60)
    
    try:
        grids, day_labels, files = load_time_series(args.input_dir)
        print(f"[LOAD] Loaded {len(grids)} time steps from {args.input_dir}")
        print(f"[LOAD] Days: {day_labels}")
        
        create_interactive_dashboard(grids, day_labels, args.output)
        
    except FileNotFoundError as e:
        print(f"[ERROR] {e}")
        print("[HINT] Run the simulation first: python src/timed_drug_infusion.py --days 120")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())