"""Audit GBM-Space sectioning metadata without reading expression matrices."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import zipfile
from pathlib import Path
from typing import Any

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output" / "trackA_v2_data_check"
OUTPUT_PATH = OUTPUT_DIR / "s14_metadata_check.json"
SOURCE_URL = (
    "https://www.biorxiv.org/content/biorxiv/early/2025/05/14/"
    "2025.05.13.653495/DC3/embed/media-3.xlsx?download=true"
)
PRECHECK_SIZE_BYTES = 136_557
CORE_TAGS = {"CT", "HBV", "MVP", "NE", "PNZ"}
EDGE_TAGS = {"IT", "LE", "INFILTRATED BRAIN", "INFILTRATED CORTEX"}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tokens_from_features(value: Any) -> set[str]:
    if pd.isna(value):
        return set()
    tokens = re.split(r"[;,.]", str(value))
    cleaned = set()
    for token in tokens:
        token = re.sub(r"\s*\([^)]*\)", "", token).strip().upper()
        if token:
            cleaned.add(token)
    return cleaned


def label_class(features: set[str]) -> str:
    has_core = bool(features & CORE_TAGS)
    has_edge = bool(features & EDGE_TAGS) or any(
        token.startswith("INFILTRATED ") for token in features
    )
    if has_core and has_edge:
        return "mixed"
    if has_core:
        return "core_only"
    if has_edge:
        return "edge_only"
    return "unresolved"


def audit(workbook_path: Path) -> dict[str, Any]:
    workbook_path = workbook_path.expanduser().resolve()
    if not workbook_path.is_file():
        raise FileNotFoundError(workbook_path)

    sheets = pd.ExcelFile(workbook_path).sheet_names
    required_sheets = {"Index", "Table S1", "Table S14"}
    if not required_sheets.issubset(sheets):
        raise ValueError(f"Missing expected metadata sheets: {sorted(required_sheets - set(sheets))}")

    # Read only cohort and serial-section metadata. Do not open expression-result tables.
    table_s1 = pd.read_excel(workbook_path, sheet_name="Table S1")
    table_s14 = pd.read_excel(workbook_path, sheet_name="Table S14")
    table_s14.columns = [str(column).strip() for column in table_s14.columns]
    required_columns = {"tumor_id", "site_id", "method", "Histopathological features"}
    if not required_columns.issubset(table_s14.columns):
        raise ValueError(f"Missing expected Table S14 columns: {sorted(required_columns - set(table_s14.columns))}")

    method = table_s14["method"].fillna("").astype(str)
    table_s14 = table_s14.assign(
        _method=method,
        _is_nuclei=method.str.contains("single nuclei", case=False),
        _features=table_s14["Histopathological features"].map(tokens_from_features),
    )
    valid_ids = table_s14["tumor_id"].notna() & table_s14["site_id"].notna()
    table_s14 = table_s14.loc[valid_ids].copy()

    sites: list[dict[str, Any]] = []
    for (donor, site), group in table_s14.groupby(["tumor_id", "site_id"], dropna=False):
        if not bool(group["_is_nuclei"].any()):
            continue
        all_features: set[str] = set()
        for feature_set in group["_features"]:
            all_features.update(feature_set)
        sites.append(
            {
                "donor_id": str(donor),
                "site_id": str(site),
                "histopathology_tags": sorted(all_features),
                "conservative_label": label_class(all_features),
            }
        )

    site_frame = pd.DataFrame(sites)
    donors_with_sites = sorted(site_frame["donor_id"].unique().tolist()) if not site_frame.empty else []
    donors_with_annotated_sites = sorted(
        site_frame.loc[site_frame["histopathology_tags"].map(bool), "donor_id"].unique().tolist()
    ) if not site_frame.empty else []
    paired_donors = []
    if not site_frame.empty:
        for donor, group in site_frame.groupby("donor_id"):
            labels = set(group["conservative_label"])
            if {"core_only", "edge_only"}.issubset(labels):
                paired_donors.append(str(donor))

    donor_column = next(
        (column for column in table_s1.columns if str(column).strip().casefold() == "patient id"),
        None,
    )
    patient_count = int(table_s1[donor_column].dropna().astype(str).nunique()) if donor_column else None
    class_counts = (
        {str(key): int(value) for key, value in site_frame["conservative_label"].value_counts().items()}
        if not site_frame.empty
        else {}
    )

    return {
        "input_file": workbook_path.name,
        "source_url": SOURCE_URL,
        "measured_size_bytes": workbook_path.stat().st_size,
        "earlier_precheck_size_bytes": PRECHECK_SIZE_BYTES,
        "size_matches_earlier_precheck": workbook_path.stat().st_size == PRECHECK_SIZE_BYTES,
        "sha256": sha256_file(workbook_path),
        "valid_xlsx_zip": zipfile.is_zipfile(workbook_path),
        "workbook_sheets": sheets,
        "sheets_read": ["Index", "Table S1", "Table S14"],
        "expression_sheets_read": False,
        "expression_values_read": False,
        "expression_scored": False,
        "patient_count_table_s1": patient_count,
        "donors_with_single_nuclei_sites_table_s14": len(donors_with_sites),
        "single_nuclei_site_count": len(sites),
        "donors_with_histopathology_features_at_single_nuclei_sites": len(donors_with_annotated_sites),
        "site_label_counts_under_conservative_screen": class_counts,
        "conservative_core_edge_paired_donor_count": len(paired_donors),
        "conservative_core_edge_paired_donors": paired_donors,
        "core_tags": sorted(CORE_TAGS),
        "edge_tags": sorted(EDGE_TAGS),
        "label_rule": (
            "Group serial-section histopathology tags by donor and site. Classify a site as core-only or edge-only "
            "only when all observed tags support one side; classify sites with both types as mixed."
        ),
        "eligibility": "FAIL: no validated nucleus-level core/edge crosswalk; conservative site screen is below 10 donors.",
        "weakest_point": (
            "Histopathology annotations come from serial Visium/H&E sections, while single-nucleus material is pooled "
            "from thick sections. A donor-site label does not establish each nucleus's anatomical region. "
            "The conservative tag rule is an audit screen, not an author-validated binary label."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", required=True, type=Path)
    args = parser.parse_args()
    result = audit(args.workbook)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
