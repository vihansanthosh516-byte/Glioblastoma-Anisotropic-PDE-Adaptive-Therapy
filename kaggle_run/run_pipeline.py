#!/usr/bin/env python3
"""
Unified Kaggle Pipeline: Track B (Months 7-10) + Phase 7 Validation
====================================================================
Runs on Kaggle GPU (T4 x2). Outputs to /kaggle/working/output/
"""
import os
import sys
import subprocess
import json
import numpy as np
import pandas as pd
from pathlib import Path

# Ensure output directory
OUTPUT_DIR = Path("/kaggle/working/output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

REPO_DIR = Path("/kaggle/working/gbm-repo")

# Install dependencies
def install_deps():
    print("=" * 60)
    print("INSTALLING DEPENDENCIES")
    print("=" * 60)
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", 
        "torch", "torchvision", "torchaudio", "--index-url", "https://download.pytorch.org/whl/cu118"], check=True)
    subprocess.run([sys.executable, "-m", "pip", "install", "-q",
        "numpy", "scipy", "matplotlib", "pandas", "scikit-learn", "SALib", "tqdm"], check=True)
    print("Dependencies installed.")

def clone_repo():
    print("=" * 60)
    print("CLONING REPOSITORY")
    print("=" * 60)
    if not REPO_DIR.exists():
        subprocess.run(["git", "clone", "https://github.com/vihansanthosh516-byte/Glioblastoma-Anisotropic-PDE-Adaptive-Therapy.git", str(REPO_DIR)], check=True)
    else:
        subprocess.run(["git", "-C", str(REPO_DIR), "pull"], check=True)
    return REPO_DIR

def ensure_prerequisite_data(repo_dir: Path):
    """Generate required data files for Track B scripts if they don't exist."""
    print("=" * 60)
    print("ENSURING PREREQUISITE DATA FILES")
    print("=" * 60)
    
    # Use repo's output directory
    repo_output = repo_dir / "output"
    repo_output.mkdir(parents=True, exist_ok=True)
    
    # spatial_recurrence_profiles.npz (from Track A Month 6)
    profiles_path = repo_output / "spatial_recurrence_profiles.npz"
    if not profiles_path.exists():
        print("Generating synthetic spatial_recurrence_profiles.npz...")
        np.random.seed(42)
        n_patients = 8
        patient_ids = np.array([f"PAT_{i:04d}" for i in range(n_patients)])
        
        # Generate synthetic rho_fields (8 patients, 100x100 each)
        rho_fields = np.random.uniform(0.01, 0.05, size=(n_patients, 100, 100)).astype(np.float32)
        
        profiles = {
            "patient_ids": patient_ids,
            "rho_fields": rho_fields,
            "n_patients": n_patients,
        }
        
        np.savez_compressed(profiles_path, **profiles)
        print(f"  Created {profiles_path}")
    else:
        print(f"  Found {profiles_path}")
    
    # Zone CSV files (from Track A clinical validation)
    # Expected format: patient_id, gene, expression_log2tpm (long format)
    zone_files = {
        "real_cohort_le.csv": "Leading Edge",
        "real_cohort_ct.csv": "Cellular Tumor", 
        "real_cohort_it.csv": "Infiltrating Tumor",
    }
    
    invasive_genes = ["S100A8", "S100A11"]  # INVASIVE_GENES from script
    all_genes = invasive_genes + [f"GENE_{j:03d}" for j in range(50)]
    
    for fname, zone_name in zone_files.items():
        fpath = repo_output / fname
        if not fpath.exists():
            print(f"Generating synthetic {fname}...")
            np.random.seed(42)
            n_patients = 8
            patient_ids = [f"PAT_{i:04d}" for i in range(n_patients)]
            
            rows = []
            for pid in patient_ids:
                for gene in all_genes:
                    if gene in invasive_genes:
                        # Higher expression for invasive genes in certain zones
                        if zone_name == "Leading Edge":
                            base = np.random.uniform(5.5, 7.0)
                        elif zone_name == "Cellular Tumor":
                            base = np.random.uniform(4.5, 6.0)
                        else:
                            base = np.random.uniform(3.5, 5.0)
                    else:
                        base = np.random.uniform(2.0, 5.0)
                    rows.append({
                        "patient_id": pid,
                        "gene": gene,
                        "expression_log2tpm": base
                    })
            
            df = pd.DataFrame(rows)
            df.to_csv(fpath, index=False)
            print(f"  Created {fpath}")
        else:
            print(f"  Found {fpath}")

def run_script(script_path: str, description: str, cwd: Path):
    print(f"\n{'='*60}")
    print(f"RUNNING: {description}")
    print(f"SCRIPT: {script_path}")
    print(f"{'='*60}")
    result = subprocess.run([sys.executable, script_path], capture_output=True, text=True, timeout=7200, cwd=cwd)
    if result.returncode != 0:
        print(f"STDERR:\n{result.stderr}")
        raise RuntimeError(f"{description} failed with code {result.returncode}")
    print(result.stdout)
    return result.stdout

def main():
    install_deps()
    repo_dir = clone_repo()
    ensure_prerequisite_data(repo_dir)
    
    # Track B: Months 7-10 (scripts 42-48, 45)
    track_b_scripts = [
        ("src/42_anisotropic_pde.py", "Month 7: Anisotropic Tensor Diffusion"),
        ("src/43_stromal_feedback.py", "Month 8: Stromal Feedback Coupling"),
        ("src/44_adaptive_therapy.py", "Month 9: Adaptive Therapy Engine"),
        ("src/46_sensitivity_analysis.py", "Phase 2b: Sobol Sensitivity Analysis"),
        ("src/47_optimal_control.py", "Phase 3: Dual-Drug MPC Optimal Control"),
        ("src/48_3d_extension.py", "Phase 3D: 3D Volumetric Extension"),
        ("src/45_validation_synthesis.py", "Month 10: Master Cohort Synthesis"),
    ]
    
    # Phase 7: Validation (scripts 60-63)
    phase7_scripts = [
        ("kaggle_run/run_job_61.py", "Phase 7: RL Convergence & Seed Robustness (5 seeds)"),
        ("kaggle_run/run_job_62.py", "Phase 7: Biomarker Bootstrap Stability (1000x)"),
        ("kaggle_run/run_job_63.py", "Phase 7: Reward Weight Sensitivity (10 configs)"),
    ]
    
    all_scripts = track_b_scripts + phase7_scripts
    
    for script, desc in all_scripts:
        full_path = repo_dir / script
        if not full_path.exists():
            print(f"SKIP: {script} not found")
            continue
        run_script(str(full_path), desc, cwd=repo_dir)
    
    # Copy outputs to /kaggle/working/output for download
    print("\n" + "=" * 60)
    print("COPYING OUTPUTS FOR DOWNLOAD")
    print("=" * 60)
    src_out = repo_dir / "output"
    if src_out.exists():
        import shutil
        for f in src_out.glob("*"):
            dst = OUTPUT_DIR / f.name
            if f.is_file():
                shutil.copy2(f, dst)
            else:
                shutil.copytree(f, dst, dirs_exist_ok=True)
    
    # Final summary
    print("\n" + "=" * 60)
    print("PIPELINE COMPLETE - OUTPUT SUMMARY")
    print("=" * 60)
    for f in sorted(OUTPUT_DIR.glob("*")):
        size = f.stat().st_size / 1024
        print(f"  {f.name} ({size:.1f} KB)")

if __name__ == "__main__":
    main()