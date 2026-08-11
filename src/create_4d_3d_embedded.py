#!/usr/bin/env python3
"""
4D 3D Tumor Viewer - Embedded Data Version (No CORS Issues)
============================================================
Embeds all 3D snapshots directly in HTML as base64/JSON.
Shows tumor at day 0 with proper visualization.
"""

import os
import json
import numpy as np
import base64
from pathlib import Path


BASE_OUTPUT = Path("output/batch_digital_twins")


def load_patient_4d_data(patient_id):
    """Load patient data and run 3D simulation if needed."""
    patient_dir = BASE_OUTPUT / patient_id
    
    sim_file = patient_dir / "clinical_simulation_data.npz"
    schedule_file = patient_dir / "treatment_schedule.json"
    inverse_file = patient_dir / "inverse_estimation.json"
    
    if not sim_file.exists() or not schedule_file.exists():
        return None
    
    try:
        sim_data = np.load(sim_file, allow_pickle=True)
        volume_history = sim_data.get('volume_history')
        drug_conc = sim_data.get('drug_concentration')
        rt_dose = sim_data.get('radiation_dose')
        final_state = sim_data.get('final_state')
        
        with open(schedule_file) as f:
            schedule = json.load(f)
        events = schedule.get('events', [])
        
        rho, D = 0.015, 0.01
        if inverse_file.exists():
            with open(inverse_file) as f:
                inv = json.load(f)
            rho = float(inv.get('rho', 0.015))
            D = float(inv.get('D', 0.01))
        
        vol_hist = volume_history.tolist() if volume_history is not None else []
        drug_list = drug_conc.tolist() if drug_conc is not None else []
        rt_list = rt_dose.tolist() if rt_dose is not None else []
        
        vol_hist = [float(v) for v in vol_hist]
        drug_list = [float(v) for v in drug_list]
        rt_list = [float(v) for v in rt_list]
        
        init_vol = float(vol_hist[0]) if vol_hist else 0.0
        final_vol = float(vol_hist[-1]) if vol_hist else 0.0
        
        return {
            'patient_id': patient_id,
            'volume_history': vol_hist,
            'drug_concentration': drug_list,
            'radiation_dose': rt_list,
            'events': events,
            'rho': float(rho),
            'D': float(D),
            'initial_volume': init_vol,
            'final_volume': final_vol,
        }
    except Exception as e:
        print(f"Error loading {patient_id}: {e}")
        return None


def encode_npy_base64(arr):
    """Encode numpy array as base64 string."""
    return base64.b64encode(arr.astype(np.float32).tobytes()).decode('ascii')


def main():
    import sys
    import base64
    
    target_patient = sys.argv[1] if len(sys.argv) > 1 else "PatientID_0045"
    
    print(f"Building embedded 3D viewer for {target_patient}...")
    
    patient_data = load_patient_4d_data(target_patient)
    if not patient_data:
        print("Failed to load patient data")
        return 1
    
    # Run quick 3D simulation to generate snapshots
    print("Running 3D simulation...")
    
    # Load patient mask
    patient_dir = Path(f"data/tcia/MU-Glioma-Post/{target_patient}")
    mask_file = patient_dir / "Timepoint_1" / f"{target_patient}_Timepoint_1_tumorMask.nii.gz"
    
    if not mask_file.exists():
        mask_file = Path(f"data/brats/{target_patient}/{target_patient}_seg.nii.gz")
        if not mask_file.exists():
            print(f"  [WARNING] No mask found for {target_patient}")
            return 1
    
    import nibabel as nib
    from scipy.ndimage import gaussian_filter
    
    nii = nib.load(str(mask_file))
    u = (nii.get_fdata() > 0).astype(np.float32)
    u = u[::2, ::2, ::2]  # 2x downsample
    print(f"  [LOAD] Mask shape: {u.shape}, initial voxels: {np.sum(u > 0.05)}")
    
    # Parameters
    rho = patient_data['rho']
    D = patient_data['D']
    alpha_drug = 0.08
    alpha_rt = 0.05
    beta_rt = 0.03
    dt = 1.0
    total_days = 120
    save_interval = 5
    
    # Build treatment schedule
    drug_conc = np.zeros(total_days + 1)
    rt_dose = np.zeros(total_days + 1)
    surgery_flags = np.zeros(total_days + 1, dtype=bool)
    
    for ev in patient_data['events']:
        start = int(ev['start_day'])
        end = int(ev['end_day']) if ev.get('end_day') else start
        ev_type = ev['type']
        
        if ev_type == 'surgery':
            if start <= total_days:
                surgery_flags[start] = True
        elif ev_type == 'chemotherapy':
            drug_name = ev.get('drug', '').lower()
            if 'temozolomide' in drug_name or 'tmz' in drug_name:
                for day in range(start, min(end + 1, total_days + 1)):
                    cycle_day = (day - start) % 28
                    if cycle_day < 5 or (day <= end and day <= 60):
                        drug_conc[day] = 1.0
            elif 'avastin' in drug_name or 'bevacizumab' in drug_name:
                for day in range(start, min(end + 1, total_days + 1), 14):
                    drug_conc[day] = 1.0
            else:
                drug_conc[start:min(end + 1, total_days + 1)] = 0.5
        elif ev_type == 'radiation':
            for day in range(start, min(end + 1, total_days + 1)):
                if day % 7 < 5:
                    rt_dose[day] = 2.0
    
    # Run simulation and collect snapshots
    os.makedirs("output/time_series_3d_embedded", exist_ok=True)
    
    # Clear old
    for f in Path("output/time_series_3d_embedded").glob("*.npy"):
        f.unlink()
    
    np.save("output/time_series_3d_embedded/tumor_3d_day_000.npy", u)
    
    snapshots = []
    day_labels = []
    snapshots.append(u.copy())
    day_labels.append(0)
    
    for day in range(1, total_days + 1):
        if day <= total_days and surgery_flags[day]:
            print(f"  [EVENT] Day {day}: Surgical resection (90% debulking)")
            u = u * 0.1
        
        C_t = drug_conc[day] if day < len(drug_conc) else 0
        rt = rt_dose[day] if day < len(rt_dose) else 0
        
        diffusion = gaussian_filter(u, sigma=0.5) - u
        growth = rho * u * (1.0 - u)
        drug_kill = 0.08 * C_t * u
        rt_kill = (0.05 * rt + 0.03 * rt**2) * u if rt > 0 else 0
        
        u = u + dt * (diffusion + growth - drug_kill - rt_kill)
        u = np.clip(u, 0.0, 1.0)
        u[u < 0.05] = 0.0
        
        if day % save_interval == 0:
            snapshots.append(u.copy())
            day_labels.append(day)
    
    print(f"  [SIM] Collected {len(snapshots)} snapshots")
    
    # Encode snapshots as base64 for embedding
    print("Encoding snapshots...")
    encoded_snapshots = []
    shapes = []
    for snap in snapshots:
        encoded = encode_npy_base64(snap)
        encoded_snapshots.append(encoded)
        shapes.append(snap.shape)
    
    # Get grid info
    first_snap = snapshots[0]
    nx, ny, nz = first_snap.shape
    ds = 2
    nx_ds, ny_ds, nz_ds = nx // ds, ny // ds, nz // ds
    
    # Create coordinate grids (downsampled)
    X, Y, Z = np.mgrid[:nx_ds, :ny_ds, :nz_ds]
    x_flat = X.flatten().tolist()
    y_flat = Y.flatten().tolist()
    z_flat = Z.flatten().tolist()
    
    # Determine value range from ALL snapshots
    all_vals = np.concatenate([s[::ds, ::ds, ::ds].flatten() for s in snapshots])
    all_nonzero = all_vals[all_vals > 0]
    if len(all_nonzero) > 0:
        vmin = float(np.percentile(all_nonzero, 5))
        vmax = float(np.percentile(all_nonzero, 95))
    else:
        vmin, vmax = 0.05, 1.0
    
    iso_levels = np.linspace(vmin, vmax, 4).tolist()
    
    print(f"  Value range: [{vmin:.4f}, {vmax:.4f}]")
    print(f"  Iso levels: {iso_levels}")
    print(f"  Grid: {nx_ds}x{ny_ds}x{nz_ds}")
    
    # Build HTML with embedded data
    html = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>4D 3D Tumor Viewer - {target_patient}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; text-align: center; }}
        .controls {{ margin: 20px 0; padding: 15px; background: #f8f9fa; border-radius: 8px; }}
        button {{ padding: 10px 20px; margin: 5px; background: #3498db; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 14px; }}
        button:hover {{ background: #2980b9; }}
        .slider-container {{ margin: 20px 0; }}
        #plot {{ width: 100%; height: 750px; }}
        .info {{ margin: 15px 0; padding: 15px; background: #e8f4fd; border-radius: 6px; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>4D 3D Tumor Viewer - {target_patient}</h1>
        <div class="info">
            <strong>Patient:</strong> {target_patient} | 
            <strong>rho:</strong> {patient_data['rho']:.6f} /day | 
            <strong>D:</strong> {patient_data['D']:.6f} mm²/day | 
            <strong>Initial Vol:</strong> {patient_data['initial_volume']:.0f} voxels | 
            <strong>Final Vol:</strong> {patient_data['final_volume']:.0f} voxels
        </div>
        
        <div class="controls">
            <button onclick="playAnimation()">Play</button>
            <button onclick="pauseAnimation()">Pause</button>
            <button onclick="prevFrame()">Prev</button>
            <button onclick="nextFrame()">Next</button>
            <label style="margin-left:20px;">Day: <span id="currentDay">0</span></label>
        </div>
        
        <div class="slider-container">
            <input type="range" id="timeSlider" min="0" max="{len(snapshots)-1}" value="0" step="1" style="width:100%;" oninput="goToFrame(this.value)">
            <div style="display:flex; justify-content:space-between; margin-top:5px; font-size:12px; color:#666;">
                {' '.join([f'<span>Day {d}</span>' for d in day_labels])}
            </div>
        </div>
        
        <div id="plot" style="width:100%; height:750px;"></div>
    </div>

    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script>
        // EMBEDDED DATA - No fetch needed!
        const encodedSnapshots = {json.dumps(encoded_snapshots)};
        const dayLabels = {json.dumps(day_labels)};
        const isoLevels = {json.dumps(iso_levels)};
        const xFlat = {json.dumps(x_flat)};
        const yFlat = {json.dumps(y_flat)};
        const zFlat = {json.dumps(z_flat)};
        const nx = {nx_ds}, ny = {ny_ds}, nz = {nz_ds};
        const vmin = {vmin}, vmax = {vmax};
        
        let currentFrame = 0;
        let animationTimer = null;
        const snapshots = [];
        
        // Decode base64 snapshots on load
        function decodeSnapshots() {{
            for (const encoded of encodedSnapshots) {{
                const binary = atob(encoded);
                const bytes = new Uint8Array(binary.length);
                for (let i = 0; i < binary.length; i++) {{
                    bytes[i] = binary.charCodeAt(i);
                }}
                const float32 = new Float32Array(bytes.buffer);
                snapshots.push(float32);
            }}
            console.log('Decoded', snapshots.length, 'snapshots');
        }}
        
        function createIsosurfaceTraces(snapshotData) {{
            const traces = [];
            for (let i = 0; i < isoLevels.length; i++) {{
                traces.push({{
                    type: 'isosurface',
                    x: xFlat, y: yFlat, z: zFlat,
                    value: Array.from(snapshotData),
                    isomin: isoLevels[i],
                    isomax: isoLevels[i],
                    surfaceCount: 1,
                    colorscale: 'Hot',
                    opacity: 0.3,
                    caps: {{x_show: false, y_show: false, z_show: false}},
                    showscale: i === 0,
                    colorbar: i === 0 ? {{title: 'Tumor Density', thickness: 20, len: 0.75}} : undefined
                }});
            }}
            return traces;
        }}
        
        function updatePlot(frameIdx) {{
            if (frameIdx >= snapshots.length) return;
            currentFrame = frameIdx;
            document.getElementById('currentDay').textContent = dayLabels[frameIdx];
            document.getElementById('timeSlider').value = frameIdx;
            
            const traces = createIsosurfaceTraces(snapshots[frameIdx]);
            Plotly.react('plot', traces, {{
                scene: {{
                    xaxis: {{title: 'X', range: [0, nx-1]}},
                    yaxis: {{title: 'Y', range: [0, ny-1]}},
                    zaxis: {{title: 'Z', range: [0, nz-1]}},
                    aspectmode: 'data',
                    camera: {{eye: {{x: 1.5, y: 1.5, z: 1.2}}}}
                }},
                title: {{text: '4D Tumor Evolution - {target_patient} (Day ' + dayLabels[frameIdx] + ')', font: {{size: 16}}}},
                margin: {{l: 0, r: 0, t: 50, b: 0}},
                height: 750
            }});
        }}
        
        function playAnimation() {{
            if (animationTimer) return;
            animationTimer = setInterval(() => {{
                const next = (currentFrame + 1) % dayLabels.length;
                updatePlot(next);
            }}, 500);
        }}
        
        function pauseAnimation() {{
            clearInterval(animationTimer);
            animationTimer = null;
        }}
        
        function prevFrame() {{
            pauseAnimation();
            updatePlot((currentFrame - 1 + dayLabels.length) % dayLabels.length);
        }}
        
        function nextFrame() {{
            pauseAnimation();
            updatePlot((currentFrame + 1) % dayLabels.length);
        }}
        
        function goToFrame(idx) {{
            pauseAnimation();
            updatePlot(parseInt(idx));
        }}
        
        // Initialize
        decodeSnapshots();
        if (snapshots.length > 0) {{
            updatePlot(0);
            document.getElementById('timeSlider').max = dayLabels.length - 1;
        }}
    </script>
</body>
</html>'''
    
    output_path = BASE_OUTPUT / f"4d_3d_embedded_{target_patient}.html"
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"[SUCCESS] Embedded 4D 3D viewer saved to: {output_path}")
    print(f"[INFO] Open directly in browser (no server needed!)")
    
    import webbrowser
    webbrowser.open(f"file://{output_path.absolute()}")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())