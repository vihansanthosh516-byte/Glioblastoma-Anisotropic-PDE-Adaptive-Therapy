#!/usr/bin/env python3
"""
4D 3D Tumor Viewer - Interactive 3D Isosurfaces with Time Slider + Patient Dropdown
===================================================================================
Shows actual 3D volumetric tumor evolution (isosurfaces) with time slider,
not just 2D charts. Patient dropdown to switch between patients.
"""

import os
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
from scipy.ndimage import gaussian_filter


BASE_OUTPUT = Path("output/batch_digital_twins")


def load_patient_data(patient_id):
    """Load all data for a specific patient."""
    patient_dir = BASE_OUTPUT / patient_id
    
    # Load clinical simulation data
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


def run_3d_simulation(patient_data, patient_id, total_days=120, save_interval=5):
    """Run 3D simulation and save time series snapshots."""
    print(f"[3D SIM] Running 3D simulation for {patient_id}...")
    
    # Load patient mask
    patient_dir = Path(f"data/tcia/MU-Glioma-Post/{patient_id}")
    mask_file = patient_dir / "Timepoint_1" / f"{patient_id}_Timepoint_1_tumorMask.nii.gz"
    
    if not mask_file.exists():
        # Try BraTS format
        mask_file = Path(f"data/brats/{patient_id}/{patient_id}_seg.nii.gz")
        if not mask_file.exists():
            print(f"  [WARNING] No mask found for {patient_id}")
            return None
    
    import nibabel as nib
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
    n_steps = total_days
    
    # Build treatment schedule
    drug_conc = np.zeros(n_steps + 1)
    rt_dose = np.zeros(n_steps + 1)
    surgery_flags = np.zeros(n_steps + 1, dtype=bool)
    
    for ev in patient_data['events']:
        start = int(ev['start_day'])
        end = int(ev['end_day']) if ev.get('end_day') else start
        ev_type = ev['type']
        
        if ev_type == 'surgery':
            if start <= n_steps:
                surgery_flags[start] = True
        elif ev_type == 'chemotherapy':
            drug_name = ev.get('drug', '').lower()
            if 'temozolomide' in drug_name or 'tmz' in drug_name:
                # Daily during concurrent RT, then 5/28 cycles
                for day in range(start, min(end + 1, n_steps + 1)):
                    cycle_day = (day - start) % 28
                    if cycle_day < 5 or day <= end and day <= 60:  # Daily during RT
                        drug_conc[day] = 1.0
            elif 'avastin' in drug_name or 'bevacizumab' in drug_name:
                for day in range(start, min(end + 1, n_steps + 1), 14):
                    drug_conc[day] = 1.0
            else:
                drug_conc[start:min(end + 1, n_steps + 1)] = 0.5
        
        elif ev_type == 'radiation':
            for day in range(start, min(end + 1, n_steps + 1)):
                if day % 7 < 5:  # Weekdays only
                    rt_dose[day] = 2.0
    
    # Run 3D simulation with snapshot saving
    os.makedirs("output/time_series_3d", exist_ok=True)
    
    # Clear old snapshots
    for f in Path("output/time_series_3d").glob("tumor_3d_day_*.npy"):
        f.unlink()
    
    # Save initial state
    np.save("output/time_series_3d/tumor_3d_day_000.npy", u)
    
    volume_history = []
    
    for day in range(1, n_steps + 1):
        # Surgery
        if day <= n_steps and surgery_flags[day]:
            print(f"  [EVENT] Day {day}: Surgical resection (90% debulking)")
            u = u * 0.1
        
        # Treatment
        C_t = drug_conc[day] if day < len(drug_conc) else 0
        rt = rt_dose[day] if day < len(rt_dose) else 0
        
        # PDE step (isotropic diffusion for speed)
        diffusion = gaussian_filter(u, sigma=0.5) - u
        growth = rho * u * (1.0 - u)
        drug_kill = 0.08 * C_t * u
        rt_kill = (alpha_rt * rt + beta_rt * rt**2) * u if rt > 0 else 0
        
        u = u + dt * (diffusion + growth - drug_kill - rt_kill)
        u = np.clip(u, 0.0, 1.0)
        u[u < 0.05] = 0.0
        
        vol = np.sum(u > 0.05)
        volume_history.append(float(vol))
        
        # Save snapshot
        if day % save_interval == 0:
            np.save(f"output/time_series_3d/tumor_3d_day_{day:03d}.npy", u)
    
    # Save final
    np.save(f"output/time_series_3d/tumor_3d_day_{n_steps:03d}.npy", u)
    
    print(f"  [SIM] Final volume: {volume_history[-1]:.0f} voxels")
    print(f"  [SIM] Saved 3D snapshots every {save_interval} days")
    
    return volume_history


def create_4d_viewer_html(patient_data, patient_ids, patient_id):
    """Create HTML with 3D isosurface viewer + time slider for one patient."""
    pid = patient_id
    p = patient_data
    
    # Check if 3D snapshots exist
    snapshot_files = sorted(Path("output/time_series_3d").glob("tumor_3d_day_*.npy"))
    if not snapshot_files:
        print(f"  [WARNING] No 3D snapshots found for {pid}")
        return None
    
    # Load first snapshot to get grid info
    first_snap = np.load(snapshot_files[0])
    nx, ny, nz = first_snap.shape
    
    # Downsample for Plotly
    ds = 2
    nx_ds, ny_ds, nz_ds = nx // ds, ny // ds, nz // ds
    
    # Create coordinate grids
    X, Y, Z = np.mgrid[:nx_ds, :ny_ds, :nz_ds]
    x_flat = X.flatten()
    y_flat = Y.flatten()
    z_flat = Z.flatten()
    
    # Determine value range
    all_vals = []
    for f in snapshot_files:
        snap = np.load(f)
        if ds > 1:
            snap = snap[::ds, ::ds, ::ds]
        all_vals.append(snap.flatten())
    all_vals = np.concatenate(all_vals)
    all_nonzero = all_vals[all_vals > 0]
    if len(all_nonzero) > 0:
        vmin = float(np.percentile(all_nonzero, 10))
        vmax = float(np.percentile(all_nonzero, 95))
    else:
        vmin, vmax = 0.05, 0.8
    
    iso_levels = np.linspace(vmin, vmax, 4)
    
    # Build frames
    frames = []
    day_labels = []
    
    for f_path in snapshot_files:
        snap = np.load(f_path)
        if ds > 1:
            snap = snap[::ds, ::ds, ::ds]
        
        day_num = int(f_path.stem.split('_')[-1])
        day_labels.append(day_num)
        
        frame_data = []
        for i, level in enumerate(iso_levels):
            frame_data.append(go.Isosurface(
                x=x_flat, y=y_flat, z=z_flat,
                value=snap.flatten(),
                isomin=level, isomax=level,
                surface_count=1,
                colorscale='Hot',
                opacity=0.3,
                caps=dict(x_show=False, y_show=False, z_show=False),
                showscale=bool(i == 0)
            ))
        
        frames.append(go.Frame(data=frame_data, name=f"Day {day_num}"))
    
    # Initial data
    initial_data = []
    first_snap_ds = first_snap[::ds, ::ds, ::ds] if ds > 1 else first_snap
    for i, level in enumerate(iso_levels):
        initial_data.append(go.Isosurface(
            x=x_flat, y=y_flat, z=z_flat,
            value=first_snap_ds.flatten(),
            isomin=level, isomax=level,
            surface_count=1,
            colorscale='Hot',
            opacity=0.3,
            caps=dict(x_show=False, y_show=False, z_show=False),
            showscale=bool(i == 0),
            colorbar=dict(title="Tumor Density", thickness=20, len=0.75) if i == 0 else None
        ))
    
    # Build dropdown options
    patient_options = ''.join([f'<option value="{pid}">{pid}</option>' for pid in patient_ids])
    
    html = f'''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>4D 3D Tumor Viewer - {pid}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 0; padding: 20px; background: #f5f5f5; }}
        .container {{ max-width: 1200px; margin: 0 auto; background: white; border-radius: 8px; padding: 20px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        h1 {{ color: #2c3e50; text-align: center; }}
        .patient-selector {{ background: #f8f9fa; padding: 15px; border-radius: 8px; margin-bottom: 20px; border: 1px solid #e9ecef; }}
        select {{ width: 100%; max-width: 300px; padding: 10px; font-size: 16px; border: 2px solid #ddd; border-radius: 6px; }}
        .info-panel {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 15px; margin: 15px 0; padding: 15px; background: #e8f4fd; border-radius: 6px; }}
        .info-box {{ text-align: center; }}
        .info-label {{ font-size: 12px; color: #7f8c8d; text-transform: uppercase; }}
        .info-value {{ font-size: 18px; font-weight: 600; color: #2c3e50; }}
        #plot {{ width: 100%; height: 700px; }}
        .controls {{ margin: 15px 0; padding: 15px; background: #f8f9fa; border-radius: 8px; }}
        button {{ padding: 10px 20px; margin: 5px; background: #3498db; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 14px; }}
        button:hover {{ background: #2980b9; }}
        .slider-container {{ margin: 20px 0; }}
    </style>
</head>
<body>
    <div class="container">
        <h1>4D 3D Tumor Viewer</h1>
        
        <div class="patient-selector">
            <label for="patientDropdown" style="font-weight:600; margin-right:10px;">Select Patient:</label>
            <select id="patientDropdown" onchange="switchPatient(this.value)">
                {''.join([f'<option value="{pid}">{pid}</option>' for pid in patient_ids])}
            </select>
        </div>
        
        <div id="infoPanel" class="info-panel" style="display:none;">
            <div class="info-box"><div class="info-label">Patient ID</div><div class="info-value" id="infoPatientId">-</div></div>
            <div class="info-box"><div class="info-label">rho</div><div class="info-value" id="infoRho">-</div></div>
            <div class="info-box"><div class="info-label">D</div><div class="info-value" id="infoD">-</div></div>
            <div class="info-box"><div class="info-label">Initial Vol</div><div class="info-value" id="infoVol0">-</div></div>
            <div class="info-box"><div class="info-label">Final Vol</div><div class="info-value" id="infoVolF">-</div></div>
            <div class="info-box"><div class="info-label">Change</div><div class="info-value" id="infoVolChange">-</div></div>
        </div>
        
        <div class="controls">
            <button onclick="playAnimation()">Play</button>
            <button onclick="pauseAnimation()">Pause</button>
            <button onclick="prevFrame()">Previous</button>
            <button onclick="nextFrame()">Next</button>
            <span style="margin-left:20px;">Day: <span id="currentDay">0</span></span>
        </div>
        
        <div id="plot" style="width:100%; height:700px;"></div>
        
        <div class="slider-container">
            <input type="range" id="timeSlider" min="0" max="{len(day_labels)-1}" value="0" step="1" style="width:100%;" oninput="goToFrame(this.value)">
            <div style="display:flex; justify-content:space-between; margin-top:5px;">
                {''.join([f'<span>Day {d}</span>' for d in day_labels[::max(1,len(day_labels)//10)]])}
            </div>
        </div>
        
        <div id="plot" style="width:100%; height:700px;"></div>
    </div>

    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script>
        const patientData = {json.dumps({pid: p for pid, p in patient_data.items()})};
        const patientIds = {json.dumps(patient_ids)};
        const currentPatient = "{pid}";
        
        // Load snapshots
        const snapshotFiles = {json.dumps([str(f) for f in snapshot_files])};
        const dayLabels = {json.dumps(day_labels)};
        const isoLevels = {json.dumps(iso_levels.tolist())};
        const vmin = {vmin};
        const vmax = {vmax};
        
        let currentFrame = 0;
        let animationTimer = null;
        let framesLoaded = false;
        
        // Pre-load all snapshots
        const snapshots = [];
        async function loadAllSnapshots() {{
            for (const f of snapshotFiles) {{
                try {{
                    const response = await fetch(f);
                    const buffer = await response.arrayBuffer();
                    const data = new Float32Array(buffer);
                    snapshots.push(data);
                }} catch (e) {{
                    console.error('Error loading', f, e);
                }}
            }}
            framesLoaded = true;
            console.log('Loaded', snapshots.length, 'snapshots');
        }}
        
        // Coordinate grids
        const nx = {nx_ds}, ny = {ny_ds}, nz = {nz_ds};
        const X = [], Y = [], Z = [];
        for (let i = 0; i < nx; i++) {{
            for (let j = 0; j < ny; j++) {{
                for (let k = 0; k < nz; k++) {{
                    X.push(i); Y.push(j); Z.push(k);
                }}
            }}
        }}
        const xFlat = X, yFlat = Y, zFlat = Z;
        
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
                title: {{text: '4D Tumor Evolution - {pid} (Day ' + dayLabels[frameIdx] + ')', font: {{size: 16}}}},
                margin: {{l: 0, r: 0, t: 50, b: 0}},
                height: 700
            }});
        }}
        
        function switchPatient(newPid) {{
            if (newPid === currentPatient) return;
            currentPatient = newPid;
            window.location.href = window.location.pathname + '?patient=' + newPid;
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
        async function init() {{
            await loadAllSnapshots();
            if (snapshots.length > 0) {{
                updatePlot(0);
                document.getElementById('timeSlider').max = dayLabels.length - 1;
            }}
        }}
        
        init();
    </script>
</body>
</html>'''
    
    return html


def main():
    print("=" * 60)
    print("  4D 3D TUMOR VIEWER WITH PATIENT DROPDOWN")
    print("=" * 60)
    
    # Get all patient IDs
    patient_ids = sorted([d.name for d in BASE_OUTPUT.glob("PatientID_*") if d.is_dir()])
    print(f"[FOUND] {len(patient_ids)} patients")
    
    # Default patient
    import sys
    target_patient = sys.argv[1] if len(sys.argv) > 1 else patient_ids[0]
    
    if target_patient not in patient_ids:
        print(f"Patient {target_patient} not found. Available: {patient_ids[:5]}...")
        target_patient = patient_ids[0]
    
    print(f"[TARGET] Building 4D viewer for {target_patient}")
    
    # Load patient data
    patient_data = {}
    for pid in patient_ids:
        data = load_patient_data(pid)
        if data:
            patient_data[pid] = data
    
    if target_patient not in patient_data:
        print(f"[ERROR] No data for {target_patient}")
        return 1
    
    # Run 3D simulation for target patient
    print(f"[SIM] Running 3D simulation for {target_patient}...")
    run_3d_simulation(patient_data[target_patient], target_patient)
    
    # Create HTML
    print(f"[BUILD] Creating 4D viewer HTML...")
    html = create_4d_viewer_html(patient_data, patient_ids, target_patient)
    
    if html is None:
        print("[ERROR] Failed to create viewer")
        return 1
    
    output_path = BASE_OUTPUT / f"4d_3d_viewer_{target_patient}.html"
    with open(output_path, 'w') as f:
        f.write(html)
    
    print(f"[SUCCESS] 4D 3D viewer saved to: {output_path}")
    
    import webbrowser
    webbrowser.open(f"file://{output_path.absolute()}")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())