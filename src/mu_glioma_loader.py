#!/usr/bin/env python3
"""Utilities for reading the MU-Glioma-Post longitudinal cohort."""
from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

DEFAULT_DATA_ROOT = Path("data/tcia/MU-Glioma-Post")
DEFAULT_CLINICAL_EXCEL = Path("data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx")


@dataclass(frozen=True)
class TimepointRecord:
    number: int
    directory: str
    mask_path: str
    volume_mm3: Optional[float] = None
    day_from_diagnosis: Optional[float] = None


@dataclass(frozen=True)
class PatientRecord:
    patient_id: str
    patient_dir: str
    timepoints: Tuple[TimepointRecord, ...]
    n_timepoints: int
    has_longitudinal_pair: bool
    treatment_schedule: Dict[str, Any] | None = None


def _volume_mm3(path: Path) -> Optional[float]:
    try:
        import nibabel as nib
        image = nib.load(str(path))
        voxel_volume = abs(float(np.linalg.det(image.affine[:3, :3])))
        return float(np.count_nonzero(image.get_fdata() > 0) * voxel_volume)
    except (ImportError, OSError, ValueError):
        return None


def clinical_treatment_schedule(path: Path, patient_id: str) -> Dict[str, Any]:
    """Extract normalized treatment timing fields from the clinical workbook."""
    if not path.exists():
        return {}
    try:
        import pandas as pd
        frame = pd.read_excel(path, sheet_name="MU Glioma Post")
        rows = frame[frame["Patient_ID"].astype(str) == patient_id]
        if rows.empty:
            return {}
        row = rows.iloc[0]
        def number(*names: str) -> Optional[float]:
            for name in names:
                if name in row.index and not pd.isna(row[name]):
                    return float(row[name])
            return None
        return {
            "surgery_day": number("Number of days from Diagnosis to First surgery or procedure "),
            "radiation_start_day": number("Number of days from Diagnosis to Radiation Therapy Start date"),
            "radiation_end_day": number("Number of days from Diagnosis to Radiation Therapy end date"),
            "tmz_start_day": number(" Number of days from Diagnosis to Initial Chemo Therapy Start date"),
            "tmz_end_day": number(" Number of days from Diagnosis to Initial Chemo Therapy end date"),
        }
    except (ImportError, OSError, ValueError, KeyError):
        return {}


def _clinical_timing(path: Path, patient_id: str) -> Dict[int, float]:
    if not path.exists():
        return {}
    try:
        import pandas as pd
        frame = pd.read_excel(path, sheet_name="MU Glioma Post")
        rows = frame[frame["Patient_ID"].astype(str) == patient_id]
        if rows.empty:
            return {}
        row = rows.iloc[0]
        timing: Dict[int, float] = {}
        for number in range(1, 7):
            prefix = f"Number of Days from Diagnosis to {number}"
            ordinal = {1: "1st", 2: "2nd", 3: "3rd", 4: "4th", 5: "5th", 6: "6th"}[number]
            candidates = [
                f"Number of Days from Diagnosis to {ordinal} MRI (Timepoint_{number}) ",
                f"Number of Days from Diagnosis to {ordinal} MRI (Timepoint_{number})",
            ]
            for column in candidates:
                if column in row.index and not pd.isna(row[column]):
                    timing[number] = float(row[column])
                    break
        return timing
    except (ImportError, OSError, ValueError, KeyError):
        return {}


def find_timepoints(patient_dir: Path) -> List[Tuple[int, Path, Path]]:
    records: List[Tuple[int, Path, Path]] = []
    for directory in patient_dir.glob("Timepoint_*"):
        if not directory.is_dir():
            continue
        try:
            number = int(directory.name.split("_")[-1])
        except ValueError:
            continue
        masks = sorted(directory.glob("*_tumorMask.nii.gz"))
        if masks:
            records.append((number, directory, masks[0]))
    return sorted(records, key=lambda item: item[0])


def load_patient_record(
    patient_id: str,
    data_root: Path = DEFAULT_DATA_ROOT,
    clinical_excel: Path = DEFAULT_CLINICAL_EXCEL,
    include_volumes: bool = True,
) -> PatientRecord:
    patient_dir = data_root / patient_id
    if not patient_dir.is_dir():
        raise FileNotFoundError(f"Patient directory not found: {patient_dir}")
    timing = _clinical_timing(clinical_excel, patient_id)
    timepoints = tuple(
        TimepointRecord(
            number=number,
            directory=str(directory),
            mask_path=str(mask),
            volume_mm3=_volume_mm3(mask) if include_volumes else None,
            day_from_diagnosis=timing.get(number),
        )
        for number, directory, mask in find_timepoints(patient_dir)
    )
    return PatientRecord(
        patient_id=patient_id,
        patient_dir=str(patient_dir),
        timepoints=timepoints,
        n_timepoints=len(timepoints),
        has_longitudinal_pair=len(timepoints) >= 2,
        treatment_schedule=clinical_treatment_schedule(clinical_excel, patient_id),
    )


def list_patient_ids(data_root: Path = DEFAULT_DATA_ROOT) -> List[str]:
    if not data_root.is_dir():
        return []
    return sorted(
        path.name for path in data_root.glob("PatientID_*")
        if path.is_dir() and find_timepoints(path)
    )


def build_cohort_metadata(
    data_root: Path = DEFAULT_DATA_ROOT,
    clinical_excel: Path = DEFAULT_CLINICAL_EXCEL,
    patient_ids: Optional[Sequence[str]] = None,
    include_volumes: bool = True,
) -> List[Dict[str, Any]]:
    ids = list(patient_ids) if patient_ids is not None else list_patient_ids(data_root)
    return [asdict(load_patient_record(pid, data_root, clinical_excel, include_volumes)) for pid in ids]


def _json_default(value: Any) -> Any:
    if isinstance(value, tuple):
        return list(value)
    raise TypeError(f"Not JSON serializable: {type(value)!r}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Index MU-Glioma-Post longitudinal tumor masks")
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--clinical-excel", type=Path, default=DEFAULT_CLINICAL_EXCEL)
    parser.add_argument("--output", type=Path, default=Path("output/mu_glioma_cohort.json"))
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--no-volumes", action="store_true")
    args = parser.parse_args()
    ids = list_patient_ids(args.data_root)
    if args.limit is not None:
        ids = ids[: args.limit]
    metadata = build_cohort_metadata(args.data_root, args.clinical_excel, ids, not args.no_volumes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(metadata, indent=2, default=_json_default), encoding="utf-8")
    longitudinal = sum(item["has_longitudinal_pair"] for item in metadata)
    print(f"Indexed {len(metadata)} patients; {longitudinal} have at least two timepoints")
    print(f"Saved metadata to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
