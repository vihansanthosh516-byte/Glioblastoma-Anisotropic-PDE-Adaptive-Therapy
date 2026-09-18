#!/usr/bin/env python3
"""
Batch Digital Twin Pipeline - Process All 203+ Patients
=======================================================
Runs the clinical-driven digital twin pipeline for all patients in MU-Glioma-Post dataset.
Generates personalized 4D dashboards for each patient.
"""

import os
import sys
import json
import subprocess
import time
import pandas as pd
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading


CLINICAL_EXCEL = "data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx"
SHEET_NAME = "MU Glioma Post"
BASE_OUTPUT = Path("output/batch_digital_twins")
BASE_OUTPUT.mkdir(parents=True, exist_ok=True)

# Thread-safe logging
log_lock = threading.Lock()
results = []
failures = []


def log(msg: str):
    """Thread-safe logging with timestamp."""
    with log_lock:
        timestamp = datetime.now().strftime("%H:%M:%S")
        # Replace emojis with ASCII for Windows console compatibility
        msg = msg.replace('\u2705', '[OK]').replace('\u274c', '[FAIL]').replace('\u23f1', '[TIMEOUT]').replace('\U0001f4a5', '[EXC]').replace('\u26a0', '[WARN]')
        print(f"[{timestamp}] {msg}")


def get_all_patients() -> list:
    """Get list of all patient directories."""
    patients = []
    base = Path("data/tcia/MU-Glioma-Post")
    for d in base.iterdir():
        if d.is_dir() and d.name.startswith("PatientID_"):
            # Check if has at least one tumor mask
            has_mask = False
            for tp in d.iterdir():
                if tp.is_dir() and tp.name.startswith("Timepoint_"):
                    if list(tp.glob("*_tumorMask.nii.gz")):
                        has_mask = True
                        break
            if has_mask:
                patients.append(d.name)
    return sorted(patients)


def run_single_patient(patient_id: str, days: int = 120) -> dict:
    """Run digital twin pipeline for a single patient."""
    patient_dir = f"data/tcia/MU-Glioma-Post/{patient_id}"
    output_dir = BASE_OUTPUT / patient_id
    output_dir.mkdir(parents=True, exist_ok=True)
    
    log(f"[{patient_id}] Starting pipeline...")
    start_time = time.time()
    
    try:
        # Run the pipeline
        cmd = [
            sys.executable, "src/run_digital_twin_pipeline.py",
            "--patient-dir", patient_dir,
            "--days", str(days),
            "--output-dir", str(output_dir)
        ]
        
        result = subprocess.run(
            cmd, 
            capture_output=True, 
            text=True, 
            cwd=".",
            timeout=600  # 10 min per patient
        )
        
        elapsed = time.time() - start_time
        
        if result.returncode == 0:
            # Check outputs
            dashboard = output_dir / "clinical_4d_dashboard.html"
            schedule = output_dir / "treatment_schedule.json"
            inverse = output_dir / "inverse_estimation.json"
            
            # Extract key metrics
            rho, D = None, None
            if inverse.exists():
                try:
                    with open(inverse) as f:
                        data = json.load(f)
                    rho = data.get('rho')
                    D = data.get('D')
                except:
                    pass
            
            log(f"[{patient_id}] ✅ SUCCESS ({elapsed:.1f}s) - rho={rho}, D={D}")
            
            return {
                "patient_id": patient_id,
                "status": "success",
                "elapsed_seconds": round(elapsed, 1),
                "rho": rho,
                "D": D,
                "dashboard": str(dashboard) if dashboard.exists() else None,
                "schedule": str(schedule) if schedule.exists() else None,
                "inverse_estimation": str(inverse) if inverse.exists() else None,
                "stdout_tail": result.stdout[-500:] if result.stdout else "",
            }
        else:
            log(f"[{patient_id}] ❌ FAILED ({elapsed:.1f}s)")
            log(f"  stderr: {result.stderr[-500:]}")
            
            return {
                "patient_id": patient_id,
                "status": "failed",
                "elapsed_seconds": round(elapsed, 1),
                "error": result.stderr[-1000:] if result.stderr else "Unknown error",
            }
            
    except subprocess.TimeoutExpired:
        log(f"[{patient_id}] ⏱️ TIMEOUT (>10 min)")
        return {
            "patient_id": patient_id,
            "status": "timeout",
            "error": "Pipeline exceeded 10 minute timeout"
        }
    except Exception as e:
        log(f"[{patient_id}] 💥 EXCEPTION: {e}")
        return {
            "patient_id": patient_id,
            "status": "exception",
            "error": str(e)
        }


def save_progress(result: dict):
    """Save result to JSONL log file."""
    log_file = BASE_OUTPUT / "batch_results.jsonl"
    with log_lock:
        with open(log_file, 'a') as f:
            f.write(json.dumps(result) + '\n')


def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Batch Digital Twin Pipeline for All Patients")
    parser.add_argument("--days", type=int, default=120, help="Simulation days per patient")
    parser.add_argument("--max-workers", type=int, default=4, help="Parallel workers (default 4)")
    parser.add_argument("--start-from", type=str, default=None, help="Start from specific patient ID")
    parser.add_argument("--limit", type=int, default=None, help="Limit number of patients")
    parser.add_argument("--skip-existing", action="store_true", help="Skip patients with existing dashboards")
    
    args = parser.parse_args()
    
    # Get all patients
    all_patients = get_all_patients()
    log(f"Found {len(all_patients)} patients with tumor masks")
    
    # Filter if starting from specific patient
    if args.start_from:
        try:
            idx = all_patients.index(args.start_from)
            all_patients = all_patients[idx:]
            log(f"Starting from {args.start_from} (index {idx})")
        except ValueError:
            log(f"Patient {args.start_from} not found")
            return 1
    
    # Skip existing
    if args.skip_existing:
        filtered = []
        for p in all_patients:
            dash = BASE_OUTPUT / p / "clinical_4d_dashboard.html"
            if not dash.exists():
                filtered.append(p)
            else:
                log(f"[{p}] ⏭️ SKIP (dashboard exists)")
        all_patients = filtered
        log(f"After skip-existing: {len(all_patients)} patients remaining")
    
    # Limit
    if args.limit:
        all_patients = all_patients[:args.limit]
        log(f"Limited to first {args.limit} patients")
    
    if not all_patients:
        log("No patients to process")
        return 0
    
    log(f"Processing {len(all_patients)} patients with {args.max_workers} workers")
    log(f"Output directory: {BASE_OUTPUT}")
    log("=" * 60)
    
    # Process patients
    successful = 0
    failed = 0
    
    if args.max_workers == 1:
        # Sequential
        for patient_id in all_patients:
            result = run_single_patient(patient_id, args.days)
            save_progress(result)
            if result["status"] == "success":
                successful += 1
            else:
                failed += 1
    else:
        # Parallel
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            future_to_patient = {
                executor.submit(run_single_patient, pid, args.days): pid 
                for pid in all_patients
            }
            
            for future in as_completed(future_to_patient):
                patient_id = future_to_patient[future]
                try:
                    result = future.result()
                    save_progress(result)
                    if result["status"] == "success":
                        successful += 1
                    else:
                        failed += 1
                except Exception as e:
                    log(f"[{patient_id}] 💥 EXCEPTION in future: {e}")
                    failed += 1
    
    # Summary
    log("=" * 60)
    log(f"BATCH COMPLETE: {successful} success, {failed} failed out of {len(all_patients)}")
    log(f"Results saved to: {BASE_OUTPUT}/batch_results.jsonl")
    
    # Generate summary CSV
    log_file = BASE_OUTPUT / "batch_results.jsonl"
    if log_file.exists():
        df = pd.read_json(log_file, lines=True)
        summary_csv = BASE_OUTPUT / "batch_summary.csv"
        df.to_csv(summary_csv, index=False)
        log(f"Summary CSV: {summary_csv}")
        
        # Print stats
        if len(df) > 0:
            success_df = df[df['status'] == 'success']
            if len(success_df) > 0:
                log(f"Average time per patient: {success_df['elapsed_seconds'].mean():.1f}s")
                if 'rho' in success_df.columns:
                    valid_rho = success_df['rho'].dropna()
                    if len(valid_rho) > 0:
                        log(f"Rho range: {valid_rho.min():.6f} - {valid_rho.max():.6f} /day")
                    valid_D = success_df['D'].dropna()
                    if len(valid_D) > 0:
                        log(f"D range: {valid_D.min():.6f} - {valid_D.max():.6f} mm2/day")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())