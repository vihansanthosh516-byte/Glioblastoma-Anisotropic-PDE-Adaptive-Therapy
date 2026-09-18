#!/usr/bin/env python3
"""
Integrated Clinical-Driven Digital Twin Pipeline
=================================================
Automated pipeline: Clinical spreadsheet + Imaging -> Inverse estimation -> 
Treatment schedule from clinical events -> 4D Simulation -> Interactive Dashboard

This pipeline reads the MU-Glioma-Post clinical spreadsheet to extract actual
treatment events (surgery, chemo, radiation, immunotherapy) with timing, and
uses them to drive the reaction-diffusion PDE simulation.
"""

import os
import sys
import json
import argparse
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Dict, Tuple, Optional

try:
    from src.treatment_aware_pde import TreatmentSchedule
except ImportError:
    import sys, pathlib
    sys.path.insert(0, str(pathlib.Path(__file__).parent.parent))
    from src.treatment_aware_pde import TreatmentSchedule


CLINICAL_EXCEL = "data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx"
SHEET_NAME = "MU Glioma Post"


def load_clinical_data(patient_id: str) -> Optional[Dict]:
    """Load and parse clinical data for a specific patient."""
    try:
        df = pd.read_excel(CLINICAL_EXCEL, sheet_name=SHEET_NAME)
        patient_row = df[df['Patient_ID'] == patient_id]
        
        if len(patient_row) == 0:
            print(f"[WARNING] No clinical record found for {patient_id}")
            return None
        
        row = patient_row.iloc[0]
        return row.to_dict()
    except Exception as e:
        print(f"[ERROR] Failed to load clinical data: {e}")
        return None


def extract_treatment_events(clinical_row: Dict) -> List[Dict]:
    """Extract treatment events from clinical data with start/end days."""
    events = []
    
    # Helper to safely get numeric values
    def get_val(key):
        val = clinical_row.get(key)
        if pd.isna(val):
            return None
        return float(val)
    
    def get_str(key):
        val = clinical_row.get(key)
        if pd.isna(val):
            return None
        return str(val).strip()
    
    # 1. Surgery (first procedure)
    surgery_day = get_val('Number of days from Diagnosis to First surgery or procedure ')
    if surgery_day is not None:
        events.append({
            'type': 'surgery',
            'start_day': max(0, surgery_day),  # Surgery can be before diagnosis day 0
            'end_day': max(0, surgery_day),
            'description': 'Surgical resection',
            'effect': 'volume_reduction_90pct'  # 90% debulking
        })
    
    # 2. Initial Chemo Therapy
    chemo_start = get_val(' Number of days from Diagnosis to Initial Chemo Therapy Start date')
    chemo_end = get_val(' Number of days from Diagnosis to Initial Chemo Therapy end date')
    chemo_name = get_str('Name of Initial Chemo Therapy')
    
    if chemo_start is not None and chemo_end is not None:
        events.append({
            'type': 'chemotherapy',
            'start_day': chemo_start,
            'end_day': chemo_end,
            'description': f'Chemo: {chemo_name}',
            'effect': 'drug_kill',
            'drug': chemo_name,
            'schedule': 'daily'  # TMZ is typically daily during radiation then 5/28 cycles
        })
    
    # 3. Radiation Therapy
    rt_start = get_val('Number of days from Diagnosis to Radiation Therapy Start date')
    rt_end = get_val('Number of days from Diagnosis to Radiation Therapy end date')
    rt_type = get_str('Radiation Therapy')
    
    if rt_start is not None and rt_end is not None:
        events.append({
            'type': 'radiation',
            'start_day': rt_start,
            'end_day': rt_end,
            'description': f'Radiation: {rt_type}',
            'effect': 'radiation_kill',
            'dose_per_fraction': 2.0,  # Gy
            'fractions_per_week': 5
        })
    
    # 4. Additional Therapy (adjuvant chemo cycles)
    add_start = get_val(' Number of Days from Diagnosis to Starting Additional Therapy ')
    add_end = get_val(' Number of Days from Diagnosis to Complete Additional Therapy ')
    add_name = get_str('Additional Therapy')
    add_cycles = clinical_row.get('Number of Cycles of Additional Therapy')
    add_cycle_len = get_val('Cycle length of Additional Therapy (q days)')
    
    if add_start is not None and add_end is not None:
        events.append({
            'type': 'chemotherapy',
            'start_day': add_start,
            'end_day': add_end,
            'description': f'Adjuvant: {add_name}',
            'effect': 'drug_kill',
            'drug': add_name,
            'cycles': add_cycles,
            'cycle_length': add_cycle_len
        })
    
    # 5. 2nd Additional Therapy
    add2_start = get_val('Number of Days from Diagnosis to Starting 2nd_Additional Therapy ')
    add2_end = get_val('Number of Days from Dagnosis to Complete 2nd_Additional Therapy ')
    add2_name = get_str('2nd_Additional Therapy')
    add2_cycle_len = get_val('Cycle length of 2nd_Additional Therapy (q days)')
    
    if add2_start is not None and add2_end is not None:
        events.append({
            'type': 'chemotherapy',
            'start_day': add2_start,
            'end_day': add2_end,
            'description': f'2nd line: {add2_name}',
            'effect': 'drug_kill',
            'drug': add2_name,
            'cycle_length': add2_cycle_len
        })
    
    # 6. Immunotherapy
    immuno_start = get_val('Number of Days from Diagnosis to Start Immunotherapy ')
    immuno_end = get_val('Number of Days from Diagnosis to Complete Immunotherapy ')
    immuno_name = get_str('Immuno therapy')
    immuno_cycle_len = get_val('Cycle length of Immunotherapy (q days)')
    
    if immuno_start is not None and immuno_end is not None:
        events.append({
            'type': 'immunotherapy',
            'start_day': immuno_start,
            'end_day': immuno_end,
            'description': f'Immuno: {immuno_name}',
            'effect': 'drug_kill',
            'drug': immuno_name,
            'cycle_length': immuno_cycle_len
        })
    
    # 7. Brachytherapy
    brachy_day = get_val('Number of Days from Diagnosis to the day of Insertion of Brachytherapy ')
    brachy_name = get_str('Brachy therapy')
    
    if brachy_day is not None:
        events.append({
            'type': 'brachytherapy',
            'start_day': brachy_day,
            'end_day': brachy_day,
            'description': f'Brachy: {brachy_name}',
            'effect': 'local_radiation'
        })
    
    # 8. Other therapies (LITT, Optune, etc.)
    other_start = get_val('Number of Days from Diagnosis to Start Other Additional Therapy ')
    other_end = get_val('Number of Days from Diagnosis to Complete Other Additional Therapy ')
    other_name = get_str('Other Types of Therapy (LITT, more chemo, proton therapy)')
    
    if other_start is not None and other_end is not None:
        events.append({
            'type': 'other',
            'start_day': other_start,
            'end_day': other_end,
            'description': f'Other: {other_name}',
            'effect': 'custom'
        })
    
    # Sort by start day
    events.sort(key=lambda x: x['start_day'] if x['start_day'] is not None else 9999)
    
    return events


def build_treatment_schedule(events: List[Dict], total_days: int, dt: float = 1.0) -> Dict:
    """Build treatment schedule arrays for the PDE simulation."""
    n_steps = int(total_days / dt) + 1
    
    # Drug concentration over time (for chemo/immuno)
    drug_conc = np.zeros(n_steps)
    # Radiation dose over time
    rt_dose = np.zeros(n_steps)
    # Surgery flags
    surgery_applied = np.zeros(n_steps, dtype=bool)
    
    for event in events:
        start_idx = int(max(0, event['start_day']) / dt)
        end_idx = int(min(total_days, event['end_day']) / dt) if event['end_day'] else start_idx
        
        if event['type'] == 'surgery':
            if start_idx < n_steps:
                surgery_applied[start_idx] = True
        
        elif event['type'] in ['chemotherapy', 'immunotherapy']:
            drug_name = event.get('drug', '').lower()
            cycle_len = event.get('cycle_length', 28)
            
            if 'temozolomide' in drug_name or 'tmz' in drug_name:
                if event.get('schedule') == 'daily':
                    drug_conc[start_idx:end_idx+1] = 1.0
                else:
                    for day in range(int(event['start_day']), int(event['end_day']) + 1):
                        cycle_day = (day - int(event['start_day'])) % int(cycle_len)
                        if cycle_day < 5:
                            idx = int(day / dt)
                            if idx < n_steps:
                                drug_conc[idx] = 1.0
            
            elif 'lomustine' in drug_name or 'ccnu' in drug_name:
                for day in range(int(event['start_day']), int(event['end_day']) + 1, int(cycle_len)):
                    idx = int(day / dt)
                    if idx < n_steps:
                        drug_conc[idx] = 1.0
            
            elif 'avastin' in drug_name or 'bevacizumab' in drug_name:
                for day in range(int(event['start_day']), int(event['end_day']) + 1, 14):
                    idx = int(day / dt)
                    if idx < n_steps:
                        drug_conc[idx] = 1.0
            
            else:
                drug_conc[start_idx:end_idx+1] = 0.5
        
        elif event['type'] == 'radiation':
            for day in range(int(event['start_day']), int(event['end_day']) + 1):
                if day % 7 < 5:
                    idx = int(day / dt)
                    if idx < n_steps:
                        rt_dose[idx] = 2.0
        
        elif event['type'] == 'brachytherapy':
            if start_idx < n_steps:
                rt_dose[start_idx] = 50.0
    
    # Convert numpy arrays to lists for JSON serialization
    return {
        'drug_concentration': drug_conc.tolist(),
        'radiation_dose': rt_dose.tolist(),
        'surgery_flags': surgery_applied.tolist(),
        'events': events
    }


def find_tumor_masks(patient_dir: str) -> Tuple[List[str], List[int]]:
    """Find tumor mask NIfTI files in patient directory."""
    patient_path = Path(patient_dir)
    
    # MU-Glioma-Post format
    timepoint_dirs = sorted([d for d in patient_path.iterdir() 
                            if d.is_dir() and d.name.startswith("Timepoint_")])
    masks = []
    tp_nums = []
    
    for tp_dir in timepoint_dirs:
        mask_files = list(tp_dir.glob("*_tumorMask.nii.gz"))
        if mask_files:
            masks.append(str(mask_files[0]))
            tp_nums.append(int(tp_dir.name.split("_")[1]))
    
    if len(masks) >= 1:
        return masks, tp_nums
    
    # BraTS format
    brats_masks = list(patient_path.glob("*_seg.nii.gz"))
    if brats_masks:
        return [str(brats_masks[0])], [0]
    
    return [], []


def run_inverse_estimation(mask_t0: str, mask_t1: str, delta_t: float, output_json: str, treatment_schedule: Optional[TreatmentSchedule] = None) -> Tuple[Optional[float], Optional[float]]:
    """Run inverse parameter estimation with optional treatment schedule."""
    try:
        import nibabel as nib
    except ImportError:
        return None, None
    
    try:
        # Load volumes from NIfTI
        img0 = nib.load(mask_t0)
        data0 = img0.get_fdata(dtype=np.float32)
        voxel_vol0 = abs(np.linalg.det(img0.affine[:3, :3]))
        V0 = float(np.sum(data0 > 0) * voxel_vol0)
        
        img1 = nib.load(mask_t1)
        data1 = img1.get_fdata(dtype=np.float32)
        voxel_vol1 = abs(np.linalg.det(img1.affine[:3, :3]))
        V1 = float(np.sum(data1 > 0) * voxel_vol1)
        
        # Call estimation API directly with treatment schedule
        import sys, pathlib, importlib.util
        spec = importlib.util.spec_from_file_location(
            "inv_est", pathlib.Path(__file__).with_name("51_inverse_parameter_estimation.py")
        )
        inv_est = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(inv_est)
        
        result = inv_est.estimate_patient_parameters(
            t0_volume=V0,
            t1_volume=V1,
            delta_t_days=delta_t,
            treatment_schedule=treatment_schedule,
        )
        
        # Save result
        import json
        result_json = {k: v for k, v in result.items() if k != "bootstrap_samples"}
        Path(output_json).parent.mkdir(parents=True, exist_ok=True)
        with open(output_json, 'w') as f:
            json.dump(result_json, f, indent=2)
        
        return result['rho'], result['D']
    except Exception as e:
        print(f"[ERROR] Inverse estimation failed: {e}")
        return None, None


def run_clinical_4d_simulation(
    mask_path: str,
    rho: float,
    D: float,
    treatment_schedule: Dict,
    total_days: int,
    output_dir: str,
    patient_id: str
) -> bool:
    """Run 4D simulation with clinical treatment schedule."""
    
    # Extract arrays from treatment schedule (already lists from build_treatment_schedule)
    drug_conc_list = treatment_schedule['drug_concentration']
    rt_dose_list = treatment_schedule['radiation_dose']
    surgery_flags_list = treatment_schedule['surgery_flags']
    
    # Save arrays to a temporary JSON file for the simulation script
    arrays_file = Path(output_dir) / "treatment_arrays.json"
    with open(arrays_file, 'w') as f:
        json.dump({
            'drug_concentration': drug_conc_list,
            'radiation_dose': rt_dose_list,
            'surgery_flags': surgery_flags_list
        }, f)
    
    # Create a modified simulation script that loads from the JSON file
    sim_script = f"""
import numpy as np
import nibabel as nib
from scipy.ndimage import gaussian_filter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import os
import json

# Load treatment schedule arrays
with open(r'{arrays_file}', 'r') as f:
    schedule = json.load(f)
drug_conc = np.array(schedule['drug_concentration'])
rt_dose = np.array(schedule['radiation_dose'])
surgery_flags = np.array(schedule['surgery_flags'])

# Load patient mask
nii = nib.load(r'{mask_path}')
u = (nii.get_fdata() > 0).astype(np.float32)
u = u[::2, ::2, ::2]  # 2x downsample
print(f"[LOAD] Mask shape: {{u.shape}}, volume: {{np.sum(u > 0.05)}} voxels")

# Parameters
rho = {rho}
D = {D}
alpha_drug = 0.08
alpha_rt = 0.05
beta_rt = 0.03
k_elim = 0.1
dt = 1.0
total_days = {total_days}
n_steps = len(drug_conc)

# Save snapshots
os.makedirs("output/time_series", exist_ok=True)
np.save("output/time_series/tumor_3d_day_000.npy", u)

volume_history = []
surgery_done = False

for day in range(1, total_days + 1):
    idx = day - 1
    
    # Surgery (apply once at start_day)
    if idx < len(surgery_flags) and surgery_flags[idx] and not surgery_done:
        print(f"[EVENT] Day {{day}}: Surgical resection (90% debulking)")
        u = u * 0.1
        surgery_done = True
    
    # Drug PK: concentration from schedule with decay
    if idx < len(drug_conc):
        C_t = drug_conc[idx]
    else:
        C_t = 0
    
    # Radiation dose
    rt = rt_dose[idx] if idx < len(rt_dose) else 0
    
    # PDE step
    # Diffusion
    diffusion = gaussian_filter(u, sigma=0.5) - u
    
    # Growth
    growth = rho * u * (1.0 - u)
    
    # Drug kill
    drug_kill = alpha_drug * C_t * u
    
    # Radiation kill (LQ model)
    if rt > 0:
        rt_kill = (alpha_rt * rt + beta_rt * rt**2) * u
    else:
        rt_kill = 0
    
    # Update
    u = u + dt * (diffusion + growth - drug_kill - rt_kill)
    u = np.clip(u, 0.0, 1.0)
    u[u < 0.05] = 0.0
    
    # Track volume
    vol = np.sum(u > 0.05)
    volume_history.append(vol)
    
    # Save snapshots every 5 days
    if day % 5 == 0:
        np.save(f"output/time_series/tumor_3d_day_{{day:03d}}.npy", u)

# Save final results
print(f"[SIM] Final volume: {{vol}} voxels")
print(f"[SIM] Total drug exposure: {{np.sum(drug_conc):.1f}} day-units")
print(f"[SIM] Total radiation: {{np.sum(rt_dose):.1f}} Gy")

# Save volume history
np.savez_compressed(r"{output_dir}/clinical_simulation_data.npz",
    volume_history=np.array(volume_history),
    drug_concentration=drug_conc,
    radiation_dose=rt_dose,
    final_state=u,
    rho=rho, D=D, total_days=total_days
)

# Plot
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes[0,0].plot(range(1, total_days+1), volume_history, 'r-', linewidth=2)
axes[0,0].set_xlabel('Day'); axes[0,0].set_ylabel('Volume (voxels)')
axes[0,0].set_title('Tumor Volume Over Time'); axes[0,0].grid(True, alpha=0.3)

axes[0,1].plot(range(total_days), drug_conc[:total_days], 'b-', label='Drug Conc.')
axes[0,1].plot(range(total_days), rt_dose[:total_days]*0.1, 'r-', label='RT Dose (x0.1)')
axes[0,1].set_xlabel('Day'); axes[0,1].set_ylabel('Intensity')
axes[0,1].set_title('Treatment Schedule'); axes[0,1].legend(); axes[0,1].grid(True, alpha=0.3)

# 3D snapshots
snap_days = [0, 30, 60, 90, 120]
for i, d in enumerate(snap_days):
    if d <= total_days:
        try:
            snap = np.load(f"output/time_series/tumor_3d_day_{{d:03d}}.npy")
            ax = axes[1, i//3] if i < 4 else axes[1, 1]
            vox = snap > 0.1
            if np.any(vox):
                ax.voxels(vox[::2, ::2, ::2], facecolors='#d62728', edgecolor='k', alpha=0.3)
            ax.set_title(f'Day {{d}}')
        except:
            pass

plt.tight_layout()
plt.savefig(r"{output_dir}/clinical_4d_simulation.png", dpi=300, bbox_inches="tight")
plt.close()

print("[SUCCESS] Clinical 4D simulation complete")
"""
    
    script_path = Path(output_dir) / "run_clinical_sim.py"
    with open(script_path, 'w') as f:
        f.write(sim_script)
    
    cmd = [sys.executable, str(script_path)]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=".", timeout=300)
    
    if result.returncode != 0:
        print(f"[ERROR] Clinical simulation failed:")
        print(result.stderr)
        return False
    
    print(result.stdout)
    return True


def generate_dashboard(time_series_dir: str, output_html: str) -> bool:
    """Generate interactive 4D dashboard."""
    cmd = [
        sys.executable, "visualization/view_3d_time_slider.py",
        "--input-dir", time_series_dir,
        "--output", output_html
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
    
    if result.returncode != 0:
        print(f"[ERROR] Dashboard generation failed: {result.stderr}")
        return False
    
    print(result.stdout)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Clinical-Driven Digital Twin Pipeline"
    )
    parser.add_argument("--patient-dir", type=str, required=True,
                        help="Path to patient directory (MU-Glioma-Post or BraTS)")
    parser.add_argument("--days", type=int, default=120,
                        help="Total simulation days")
    parser.add_argument("--output-dir", type=str, default="output/clinical_digital_twin",
                        help="Output directory")
    parser.add_argument("--open-browser", action="store_true",
                        help="Open dashboard in browser")
    
    args = parser.parse_args()
    
    patient_dir = Path(args.patient_dir)
    patient_id = patient_dir.name
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("  CLINICAL-DRIVEN DIGITAL TWIN PIPELINE")
    print("=" * 70)
    print(f"Patient: {patient_id}")
    print(f"Input: {patient_dir}")
    print(f"Output: {output_dir}")
    print(f"Simulation days: {args.days}")
    print("=" * 70)
    
    # Step 1: Load clinical data
    print("\n[STEP 1] Loading clinical metadata...")
    clinical_row = load_clinical_data(patient_id)
    
    if clinical_row is None:
        print("[WARNING] No clinical data found, using generic protocol")
        events = []
    else:
        print(f"[SUCCESS] Loaded clinical record for {patient_id}")
        # Extract treatment events
        print("\n[STEP 2] Extracting treatment events from clinical data...")
        events = extract_treatment_events(clinical_row)
        
        print(f"Found {len(events)} treatment events:")
        for ev in events:
            print(f"  Day {ev['start_day']:.0f}-{ev['end_day']:.0f}: {ev['description']} ({ev['type']})")
    
    # Build treatment schedule
    print("\n[STEP 3] Building treatment schedule...")
    treatment_schedule = build_treatment_schedule(events, args.days)
    tmz_days = tuple(float(i) for i, value in enumerate(treatment_schedule['drug_concentration']) if value > 0)
    pde_schedule = TreatmentSchedule(tmz_bolus_days=tmz_days)
    
    # Save schedule for reference
    with open(output_dir / "treatment_schedule.json", 'w') as f:
        json.dump(treatment_schedule, f, indent=2, default=str)
    
    # Step 4: Find tumor masks
    print("\n[STEP 4] Finding tumor segmentation masks...")
    masks, tp_nums = find_tumor_masks(str(patient_dir))
    
    if len(masks) == 0:
        print("[ERROR] No tumor masks found")
        return 1
    
    print(f"Found {len(masks)} timepoint(s)")
    
    # Step 5: Inverse estimation (if 2+ timepoints)
    rho, D = 0.015, 0.01  # Population defaults
    
    if len(masks) >= 2:
        print("\n[STEP 5] Running inverse estimation...")
        
        # Get clinical timing for delta_t
        tp1_day = clinical_row.get('Number of Days from Diagnosis to 1st MRI (Timepoint_1) ', 0) if clinical_row else 0
        tp2_day = clinical_row.get('Number of Days from Diagnosis to 2nd MRI (Timepoint_2) ', 56) if clinical_row else 56
        
        if pd.isna(tp1_day) or pd.isna(tp2_day):
            delta_t = 56
        else:
            delta_t = float(tp2_day) - float(tp1_day)
        
        est_output = output_dir / "inverse_estimation.json"
        est_rho, est_D = run_inverse_estimation(masks[0], masks[1], delta_t, str(est_output), pde_schedule)
        
        if est_rho is not None:
            rho, D = est_rho, est_D
            print(f"[SUCCESS] Personalized: rho={rho:.6f}/day, D={D:.6f} mm2/day")
        else:
            print(f"[FALLBACK] Using defaults: rho={rho}, D={D}")
    else:
        print(f"\n[STEP 5] Single timepoint - using population defaults: rho={rho}, D={D}")
    
    # Step 6: Run clinical 4D simulation
    print(f"\n[STEP 6] Running clinical 4D simulation ({args.days} days)...")
    success = run_clinical_4d_simulation(
        mask_path=masks[0],
        rho=rho,
        D=D,
        treatment_schedule=treatment_schedule,
        total_days=args.days,
        output_dir=str(output_dir),
        patient_id=patient_id
    )
    
    if not success:
        return 1
    
    # Step 7: Generate dashboard
    print("\n[STEP 7] Generating interactive 4D dashboard...")
    time_series_dir = "output/time_series"
    dashboard_path = output_dir / "clinical_4d_dashboard.html"
    
    success = generate_dashboard(time_series_dir, str(dashboard_path))
    
    if not success:
        return 1
    
    # Summary
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print("=" * 70)
    print(f"Patient: {patient_id}")
    print(f"Clinical events used: {len(events)}")
    print(f"Personalized parameters: rho={rho:.6f}/day, D={D:.6f} mm2/day")
    print(f"Treatment schedule: {output_dir}/treatment_schedule.json")
    print(f"4D simulation snapshots: {time_series_dir}/")
    print(f"Interactive dashboard: {dashboard_path}")
    print("=" * 70)
    
    if args.open_browser:
        import webbrowser
        webbrowser.open(f"file://{dashboard_path.absolute()}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())