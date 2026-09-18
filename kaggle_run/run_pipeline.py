#!/usr/bin/env python3
"""
Unified Kaggle Pipeline: Phase 13 (PPO Training) + Phase 15 (Virtual Trial)
====================================================================
Runs on Kaggle GPU (T4 x2). Outputs to /kaggle/working/output/
"""
import os
import sys
import subprocess
import json
import numpy as np
import pandas as pd
import shutil
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
        "numpy", "scipy", "matplotlib", "pandas", "scikit-learn", "SALib", "tqdm",
        "stable-baselines3", "gymnasium", "pyserial", "scipy"], check=True)
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
    """Generate required data files for Phase 13/15 scripts if they don't exist."""
    print("=" * 60)
    print("ENSURING PREREQUISITE DATA FILES")
    print("=" * 60)
    
    repo_output = repo_dir / "output"
    repo_output.mkdir(parents=True, exist_ok=True)
    
    # spatial_recurrence_profiles.npz
    profiles_path = repo_output / "spatial_recurrence_profiles.npz"
    if not profiles_path.exists():
        print("Generating synthetic spatial_recurrence_profiles.npz...")
        np.random.seed(42)
        n_patients = 8
        patient_ids = np.array([f"PAT_{i:04d}" for i in range(n_patients)])
        rho_fields = np.random.uniform(0.01, 0.05, size=(n_patients, 100, 100)).astype(np.float32)
        profiles = {"patient_ids": patient_ids, "rho_fields": rho_fields, "n_patients": n_patients}
        np.savez_compressed(profiles_path, **profiles)
        print(f"  Created {profiles_path}")
    else:
        print(f"  Found {profiles_path}")
    
    # Zone CSV files
    zone_files = {
        "real_cohort_le.csv": "Leading Edge",
        "real_cohort_ct.csv": "Cellular Tumor", 
        "real_cohort_it.csv": "Infiltrating Tumor",
    }
    invasive_genes = ["S100A8", "S100A11"]
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
                    if gene in ["S100A8", "S100A11"]:
                        if zone_name == "Leading Edge":
                            base = np.random.uniform(5.5, 7.0)
                        elif zone_name == "Cellular Tumor":
                            base = np.random.uniform(4.5, 6.0)
                        else:
                            base = np.random.uniform(3.5, 5.0)
                    else:
                        base = np.random.uniform(2.0, 5.0)
                    rows.append({"patient_id": pid, "gene": gene, "expression_log2tpm": base})
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
    
    # Phase 13: Circadian-Aware RL Training
    phase13_scripts = [
        ("src/rl/train_chronotherapy.py", "Phase 13: Circadian-Aware PPO Training"),
    ]
    
    # Phase 15: Virtual Clinical Trial
    phase15_scripts = [
        ("src/phase15_virtual_trial.py", "Phase 15: 1000-Patient Virtual Clinical Trial"),
    ]
    
    all_phases = [
        ("Phase 13 (PPO)", phase13_scripts),
        ("Phase 15 (Trial)", phase15_scripts),
    ]
    
    for phase_name, scripts in all_phases:
        print(f"\n{'='*60}")
        print(f"PHASE: {phase_name}")
        print(f"{'='*60}")
        for script, desc in scripts:
            full_path = REPO_DIR / script
            if not full_path.exists():
                print(f"SKIP: {script} not found")
                continue
            run_script(str(full_path), desc, cwd=REPO_DIR)
    
    # Copy outputs to /kaggle/working/output for download
    print("\n" + "=" * 60)
    print("COPYING OUTPUTS FOR DOWNLOAD")
    print("=" * 60)
    src_out = REPO_DIR / "output"
    if src_out.exists():
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