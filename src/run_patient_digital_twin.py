#!/usr/bin/env python3
"""
End-to-End Patient Digital Twin Pipeline
=========================================
Automated pipeline: Clinical imaging -> Inverse estimation -> 4D Simulation -> Interactive Dashboard

Usage:
    python src/run_patient_digital_twin.py --patient-dir data/tcia/MU-Glioma-Post/PatientID_0041 --days 120
    python src/run_patient_digital_twin.py --patient-dir data/brats/BraTS2021_00000 --days 120 --infusion-days 30 60 90
"""

import os
import sys
import json
import argparse
import subprocess
import numpy as np
import nibabel as nib
from pathlib import Path


def find_tumor_masks(patient_dir: str):
    """Find tumor mask NIfTI files in patient directory (MU-Glioma-Post or BraTS format)."""
    patient_path = Path(patient_dir)
    
    # MU-Glioma-Post format: Timepoint_N/PatientID_Timepoint_N_tumorMask.nii.gz
    timepoint_dirs = sorted([d for d in patient_path.iterdir() if d.is_dir() and d.name.startswith("Timepoint_")])
    masks = []
    days = []
    
    for tp_dir in timepoint_dirs:
        mask_files = list(tp_dir.glob("*_tumorMask.nii.gz"))
        if mask_files:
            masks.append(str(mask_files[0]))
            # Extract timepoint number
            tp_num = int(tp_dir.name.split("_")[1])
            days.append(tp_num)
    
    if len(masks) >= 2:
        return masks, days
    
    # BraTS format: single segmentation file
    brats_masks = list(patient_path.glob("*_seg.nii.gz"))
    if brats_masks:
        return [str(brats_masks[0])], [0]
    
    return [], []


def load_clinical_timing(clinical_excel: str, patient_id: str):
    """Load timepoint days from clinical Excel file."""
    try:
        import pandas as pd
        df = pd.read_excel(clinical_excel, sheet_name='MU Glioma Post')
        p = df[df['Patient_ID'] == patient_id]
        if len(p) == 0:
            return None
        
        cols = {
            'tp1': 'Number of Days from Diagnosis to 1st MRI (Timepoint_1) ',
            'tp2': 'Number of Days from Diagnosis to 2nd MRI (Timepoint_2) ',
            'tp3': 'Number of Days from Diagnosis to 3rd MRI (Timepoint_3) ',
            'tp4': 'Number of Days from Diagnosis to 4th MRI (Timepoint_4) ',
            'tp5': 'Number of Days from Diagnosis to 5th MRI (Timepoint_5) ',
            'tp6': 'Number of Days from Diagnosis to 6th MRI (Timepoint_6) ',
        }
        days = []
        for key, col in cols.items():
            val = p[col].values[0]
            if not pd.isna(val):
                days.append(float(val))
        return days
    except Exception as e:
        print(f"[WARNING] Could not load clinical timing: {e}")
        return None


def run_inverse_estimation(mask_t0: str, mask_t1: str, delta_t: float, output_json: str):
    """Run inverse parameter estimation and return rho, D."""
    cmd = [
        sys.executable, "src/51_inverse_parameter_estimation.py",
        "--nifti-t0", mask_t0,
        "--nifti-t1", mask_t1,
        "--delta-t", str(delta_t),
        "--output", output_json
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
    
    if result.returncode != 0:
        print(f"[ERROR] Inverse estimation failed: {result.stderr}")
        return None, None
    
    # Parse output for rho, D
    try:
        with open(output_json, 'r') as f:
            data = json.load(f)
        rho = data.get('rho', 0.015)
        D = data.get('D', 0.01)
        return rho, D
    except Exception as e:
        print(f"[WARNING] Could not parse inverse estimation output: {e}")
        return 0.015, 0.01


def run_4d_simulation(patient_id: str, mask_path: str, rho: float, D: float, days: int, infusion_days: list, output_dir: str):
    """Run 4D simulation with personalized parameters."""
    # Use the timed_drug_infusion.py with custom parameters
    cmd = [
        sys.executable, "src/timed_drug_infusion.py",
        "--patient", patient_id,
        "--mask-path", mask_path,
        "--days", str(days),
        "--rho", str(rho),
        "--alpha", "0.08",
        "--k_elim", "0.1",
        "--infusion-days", *[str(d) for d in infusion_days],
        "--output", f"{output_dir}/patient_4d_simulation.png"
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=".")
    
    if result.returncode != 0:
        print(f"[ERROR] 4D simulation failed: {result.stderr}")
        return False
    
    print(result.stdout)
    return True


def generate_dashboard(time_series_dir: str, output_html: str):
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
        description="End-to-End Patient Digital Twin Pipeline"
    )
    parser.add_argument("--patient-dir", type=str, required=True,
                        help="Path to patient directory (MU-Glioma-Post or BraTS format)")
    parser.add_argument("--days", type=int, default=120,
                        help="Total simulation days")
    parser.add_argument("--infusion-days", type=int, nargs="*", default=[30, 60, 90],
                        help="Days when drug infusions occur")
    parser.add_argument("--clinical-excel", type=str, default=None,
                        help="Path to clinical data Excel (for MU-Glioma-Post)")
    parser.add_argument("--output-dir", type=str, default="output/patient_digital_twin",
                        help="Output directory for results")
    parser.add_argument("--open-browser", action="store_true",
                        help="Open dashboard in browser after generation")
    
    args = parser.parse_args()
    
    patient_dir = Path(args.patient_dir)
    patient_id = patient_dir.name
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 70)
    print("  END-TO-END PATIENT DIGITAL TWIN PIPELINE")
    print("=" * 70)
    print(f"Patient: {patient_id}")
    print(f"Input: {patient_dir}")
    print(f"Output: {output_dir}")
    print(f"Simulation days: {args.days}")
    print(f"Infusion days: {args.infusion_days}")
    print("=" * 70)
    
    # Step 1: Find tumor masks
    print("\n[STEP 1] Finding tumor segmentation masks...")
    masks, tp_indices = find_tumor_masks(str(patient_dir))
    
    if len(masks) == 0:
        print("[ERROR] No tumor masks found")
        return 1
    
    print(f"Found {len(masks)} timepoint(s): {masks}")
    
    # Step 2: Inverse estimation (if we have 2+ timepoints)
    rho, D = 0.015, 0.01  # Population defaults
    
    if len(masks) >= 2:
        # Step 2: Get clinical timing
        clinical_days = None
        if args.clinical_excel and os.path.exists(args.clinical_excel):
            clinical_days = load_clinical_timing(args.clinical_excel, patient_id)
            if clinical_days:
                print(f"Clinical timepoints (days from diagnosis): {clinical_days}")
        
        # Use first two timepoints for inverse estimation
        mask_t0 = masks[0]
        mask_t1 = masks[1]
        
        if clinical_days and len(clinical_days) >= 2:
            delta_t = clinical_days[1] - clinical_days[0]
        else:
            delta_t = 56  # default fallback
            print(f"[WARNING] Using default delta_t = {delta_t} days")
        
        print(f"\n[STEP 2] Running inverse estimation (TP1 -> TP2, dt={delta_t} days)...")
        est_output = output_dir / "inverse_estimation.json"
        est_rho, est_D = run_inverse_estimation(mask_t0, mask_t1, delta_t, str(est_output))
        
        if est_rho is not None:
            rho, D = est_rho, est_D
            print(f"[SUCCESS] Personalized parameters: rho={rho:.6f} /day, D={D:.6f} mm2/day")
        else:
            print(f"[FALLBACK] Using population defaults: rho={rho}, D={D}")
    else:
        print(f"\n[STEP 2] Single timepoint detected (BraTS format). Using population defaults: rho={rho}, D={D}")
        # For BraTS, use first mask for simulation
        masks = [masks[0]]
    
    # Step 3: Run 4D simulation
    print(f"\n[STEP 3] Running 4D simulation ({args.days} days)...")
    success = run_4d_simulation(
        patient_id=patient_id,
        mask_path=masks[0],
        rho=rho,
        D=D,
        days=args.days,
        infusion_days=args.infusion_days,
        output_dir=str(output_dir)
    )
    
    if not success:
        return 1
    
    # Step 4: Generate interactive dashboard
    print(f"\n[STEP 4] Generating interactive 4D dashboard...")
    time_series_dir = "output/time_series"
    dashboard_path = output_dir / "patient_4d_dashboard.html"
    
    success = generate_dashboard(time_series_dir, str(dashboard_path))
    
    if not success:
        return 1
    
    # Summary
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE")
    print("=" * 70)
    if 'est_output' in locals():
        print(f"Inverse estimation: {est_output}")
    else:
        print("Inverse estimation: SKIPPED (single timepoint - used population defaults)")
    print(f"4D simulation snapshots: {time_series_dir}/")
    print(f"Interactive dashboard: {dashboard_path}")
    print("=" * 70)
    
    if args.open_browser:
        import webbrowser
        webbrowser.open(f"file://{dashboard_path.absolute()}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())