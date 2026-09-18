#!/usr/bin/env python3
"""Single entry point for the real-patient digital-twin workflow.

Examples:
  python run_digital_twin.py --mode index
  python run_digital_twin.py --mode patient --patient-id PatientID_0003 --days 90
  python run_digital_twin.py --mode cohort --limit 5
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List

from src.mu_glioma_loader import (
    DEFAULT_CLINICAL_EXCEL,
    DEFAULT_DATA_ROOT,
    build_cohort_metadata,
    list_patient_ids,
    load_patient_record,
)


def _patient_command(args: argparse.Namespace, patient_id: str) -> List[str]:
    output = Path(args.output_dir) / patient_id
    return [
        sys.executable,
        "src/run_digital_twin_pipeline.py",
        "--patient-dir", str(Path(args.data_root) / patient_id),
        "--days", str(args.days),
        "--output-dir", str(output),
    ]


def run_patient(args: argparse.Namespace) -> int:
    record = load_patient_record(args.patient_id, args.data_root, args.clinical_excel)
    if not record.timepoints:
        print(f"No tumor masks found for {args.patient_id}", file=sys.stderr)
        return 1
    print(f"{args.patient_id}: {record.n_timepoints} timepoint(s)")
    completed = subprocess.run(_patient_command(args, args.patient_id), check=False)
    return completed.returncode


def run_index(args: argparse.Namespace) -> int:
    metadata = build_cohort_metadata(args.data_root, args.clinical_excel, include_volumes=not args.no_volumes)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(metadata, indent=2, default=list), encoding="utf-8")
    print(f"Indexed {len(metadata)} patients; saved to {output}")
    return 0


def run_cohort(args: argparse.Namespace) -> int:
    patient_ids = list_patient_ids(args.data_root)
    if args.limit is not None:
        patient_ids = patient_ids[:args.limit]
    results: List[Dict[str, Any]] = []
    for patient_id in patient_ids:
        record = load_patient_record(patient_id, args.data_root, args.clinical_excel)
        item: Dict[str, Any] = {
            "patient_id": patient_id,
            "n_timepoints": record.n_timepoints,
            "has_longitudinal_pair": record.has_longitudinal_pair,
            "status": "indexed",
        }
        if args.run:
            completed = subprocess.run(_patient_command(args, patient_id), check=False)
            item["status"] = "success" if completed.returncode == 0 else "failed"
            item["return_code"] = completed.returncode
        results.append(item)
        print(f"{patient_id}: {item['status']}")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"Processed {len(results)} patients; saved to {output}")
    return 0 if all(item["status"] == "success" or not args.run for item in results) else 1


def main() -> int:
    parser = argparse.ArgumentParser(description="Real-patient GBM digital-twin entry point")
    parser.add_argument("--mode", choices=("index", "patient", "cohort"), required=True)
    parser.add_argument("--patient-id")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--clinical-excel", type=Path, default=DEFAULT_CLINICAL_EXCEL)
    parser.add_argument("--days", type=int, default=120)
    parser.add_argument("--output-dir", type=Path, default=Path("output/digital_twins"))
    parser.add_argument("--output", type=Path, default=Path("output/mu_glioma_cohort.json"))
    parser.add_argument("--limit", type=int)
    parser.add_argument("--run", action="store_true", help="Run the existing simulation for cohort mode")
    parser.add_argument("--no-volumes", action="store_true")
    args = parser.parse_args()
    if args.mode == "patient":
        if not args.patient_id:
            parser.error("--patient-id is required in patient mode")
        return run_patient(args)
    if args.mode == "index":
        return run_index(args)
    return run_cohort(args)


if __name__ == "__main__":
    raise SystemExit(main())
