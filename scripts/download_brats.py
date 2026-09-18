#!/usr/bin/env python3
"""
BraTS Dataset Downloader for GBM Digital Twin
==============================================
Downloads BraTS 2021/2023 training data from public sources.

Usage:
    python scripts/download_brats.py --year 2021 --output-dir data/brats
    python scripts/download_brats.py --year 2023 --output-dir data/brats --n-patients 50

Note: BraTS data requires registration at https://www.med.upenn.edu/sbia/brats2021/
For automated download, use the BraTS 2023 Kaggle dataset:
    kaggle datasets download -d awsaf49/brats2023-dataset
"""
import argparse
import os
import sys
import subprocess
from pathlib import Path


def download_brats_kaggle(output_dir: Path, year: int = 2023, unzip: bool = True):
    """
    Download BraTS dataset from Kaggle.
    
    Requires: kaggle CLI configured with API token
    pip install kaggle
    kaggle config set -n path -v /path/to/kaggle.json
    """
    print(f"Downloading BraTS {year} from Kaggle...")
    
    if year == 2021:
        dataset = "awsaf49/brats2021-dataset"
    elif year == 2023:
        dataset = "awsaf49/brats2023-dataset"
    else:
        raise ValueError(f"Unsupported year: {year}")
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Download
    cmd = ["kaggle", "datasets", "download", "-d", dataset, "-p", str(output_dir)]
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
        print("\nMake sure you have:")
        print("1. Installed kaggle: pip install kaggle")
        print("2. Configured API token: kaggle config set -n path -v ~/.kaggle/kaggle.json")
        print("3. Accepted competition rules on Kaggle website")
        return False
    
    # Find and unzip
    zip_files = list(output_dir.glob("*.zip"))
    if zip_files and unzip:
        import zipfile
        for zip_file in zip_files:
            print(f"Extracting {zip_file}...")
            with zipfile.ZipFile(zip_file, 'r') as z:
                z.extractall(output_dir)
            zip_file.unlink()
    
    print(f"Downloaded and extracted to {output_dir}")
    return True


def download_tcga_gdc(output_dir: Path, n_patients: int = 50):
    """
    Download TCGA-GBM data via GDC API.
    
    Note: This is a simplified version. For full download, use GDC Data Transfer Tool.
    """
    print("TCGA-GBM download via GDC API...")
    print("Note: This requires the GDC Data Transfer Tool for bulk download.")
    print("For programmatic access, see: https://docs.gdc.cancer.gov/API/Users_Guide/")
    
    # Example query for TCGA-GBM
    # Would need to query files with:
    # - cases.project.project_id = TCGA-GBM
    # - files.data_category = Transcriptome Profiling / DNA Methylation / etc.
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save example query
    query = {
        "filters": {
            "op": "and",
            "content": [
                {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-GBM"}},
                {"op": "=", "content": {"field": "files.data_type", "value": "Gene Expression Quantification"}}
            ]
        },
        "format": "JSON",
        "fields": "file_id,file_name,cases.submitter_id",
        "size": n_patients
    }
    
    import json
    query_path = output_dir / "gdc_query.json"
    with open(query_path, 'w') as f:
        json.dump(query, f, indent=2)
    
    print(f"Saved example GDC query to {query_path}")
    print("To download, use: gdc-client download -m gdc_query.json")
    
    return True


def organize_brats_data(data_root: Path):
    """
    Organize downloaded BraTS data into patient directories.
    
    BraTS typically comes as:
    Brats2021_00001_t1.nii.gz
    Brats2021_00001_t1ce.nii.gz
    Brats2021_00001_t2.nii.gz
    Brats2021_00001_flair.nii.gz
    Brats2021_00001_seg.nii.gz
    
    This organizes into:
    data/brats/Brats2021_00001/
        Brats2021_00001_t1.nii.gz
        ...
    """
    data_root = Path(data_root)
    
    # Find all NIfTI files
    nifti_files = list(data_root.rglob("*.nii.gz"))
    
    if not nifti_files:
        print("No NIfTI files found")
        return
    
    print(f"Found {len(nifti_files)} NIfTI files")
    
    # Group by patient ID
    patients = {}
    for f in nifti_files:
        # Extract patient ID (e.g., Brats2021_00001)
        name = f.name
        if '_' in name:
            # Typical format: Brats2021_00001_t1.nii.gz
            parts = name.split('_')
            if len(parts) >= 2:
                patient_id = f"{parts[0]}_{parts[1]}"
            else:
                patient_id = parts[0]
        else:
            patient_id = f.stem
        
        if patient_id not in patients:
            patients[patient_id] = []
        patients[patient_id].append(f)
    
    print(f"Found {len(patients)} patients")
    
    # Move to organized structure
    for patient_id, files in patients.items():
        patient_dir = data_root / patient_id
        patient_dir.mkdir(exist_ok=True)
        
        for f in files:
            dest = patient_dir / f.name
            if not dest.exists():
                f.rename(dest)
    
    print(f"Organized into {len(patients)} patient directories")
    return True


def main():
    parser = argparse.ArgumentParser(
        description="Download BraTS/TCGA data for GBM Digital Twin"
    )
    parser.add_argument(
        "--source", type=str, default="kaggle",
        choices=["kaggle", "tcga"],
        help="Data source (kaggle for BraTS, tcga for TCGA-GBM)"
    )
    parser.add_argument(
        "--year", type=int, default=2023,
        choices=[2021, 2023],
        help="BraTS challenge year (for Kaggle source)"
    )
    parser.add_argument(
        "--output-dir", type=str, default="data/brats",
        help="Output directory for downloaded data"
    )
    parser.add_argument(
        "--n-patients", type=int, default=50,
        help="Number of patients to download (TCGA only)"
    )
    parser.add_argument(
        "--organize", action="store_true",
        help="Organize existing BraTS files into patient directories"
    )
    
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    
    if args.organize:
        print(f"Organizing BraTS data in {args.output_dir}...")
        organize_brats_data(Path(args.output_dir))
        return
    
    if args.source == "kaggle":
        download_brats_kaggle(output_dir, args.year)
    elif args.source == "tcga":
        download_tcga_gdc(output_dir, args.n_patients)


if __name__ == "__main__":
    main()