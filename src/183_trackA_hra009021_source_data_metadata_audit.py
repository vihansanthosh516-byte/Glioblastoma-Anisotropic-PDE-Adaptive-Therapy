"""Inspect source workbook text metadata only; never read numeric source values."""

from __future__ import annotations

import hashlib
import json
import re
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parent.parent
WORKBOOK_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "41467_2024_55243_MOESM4_ESM.xlsx"
METADATA_AUDIT_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "hra009021_metadata_audit.json"
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "hra009021_source_data_audit.json"
EXPECTED_BYTES = 50_969_169


def inspect_string(value: str, identifiers: set[str], region_labels: set[str], identifier_hits: set[str]) -> None:
    value = value.strip()
    if value in identifiers:
        identifier_hits.add(value)
    if re.fullmatch(r"(?:HBI|LBI|Tumor)\d{2}", value):
        region_labels.add(value)


def main() -> None:
    metadata = json.loads(METADATA_AUDIT_PATH.read_text(encoding="utf-8"))
    identifiers = {
        str(row[field]).strip()
        for row in metadata["sample_crosswalk"]
        for field in (
            "individual_accession",
            "individual_identifier",
            "sample_accession",
            "sample_name",
            "experiment_accession",
            "run_accession",
        )
        if row.get(field)
    }
    identifier_hits: set[str] = set()
    region_labels: set[str] = set()

    if WORKBOOK_PATH.stat().st_size != EXPECTED_BYTES:
        raise ValueError(f"Unexpected workbook size: {WORKBOOK_PATH.stat().st_size}")

    with zipfile.ZipFile(WORKBOOK_PATH) as archive:
        workbook_xml = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet_names = [
            element.attrib["name"]
            for element in workbook_xml.iter()
            if element.tag.endswith("sheet") and "name" in element.attrib
        ]

        # Inspect shared text only. Numeric <v> values are never read or emitted.
        if "xl/sharedStrings.xml" in archive.namelist():
            with archive.open("xl/sharedStrings.xml") as strings_file:
                for _, item in ET.iterparse(strings_file, events=("end",)):
                    if item.tag.endswith("si"):
                        inspect_string("".join(item.itertext()), identifiers, region_labels, identifier_hits)
                        item.clear()

        # Also check inline string cells, without reading numeric cell values.
        for member in archive.namelist():
            if not (member.startswith("xl/worksheets/sheet") and member.endswith(".xml")):
                continue
            with archive.open(member) as sheet_file:
                for _, element in ET.iterparse(sheet_file, events=("end",)):
                    if element.tag.endswith("c") and element.attrib.get("t") == "inlineStr":
                        inspect_string("".join(element.itertext()), identifiers, region_labels, identifier_hits)
                        element.clear()

    sample_names = {
        row["sample_name"]
        for row in metadata["sample_crosswalk"]
        if row.get("sample_name")
    }
    region_counts = {
        region: sum(label.startswith(region) for label in region_labels)
        for region in ("HBI", "LBI", "Tumor")
    }
    result = {
        "dataset": "HRA009021",
        "source_workbook": WORKBOOK_PATH.name,
        "source_url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-55243-5/MediaObjects/41467_2024_55243_MOESM4_ESM.xlsx",
        "source_file_bytes": WORKBOOK_PATH.stat().st_size,
        "source_file_sha256": hashlib.sha256(WORKBOOK_PATH.read_bytes()).hexdigest(),
        "sheet_names": sheet_names,
        "metadata_only_inspection": True,
        "numeric_source_values_read": False,
        "expression_scored": False,
        "gsa_sample_names_tested": len(sample_names),
        "gsa_sample_names_found": sorted(sample_names.intersection(identifier_hits)),
        "gsa_accession_or_run_identifiers_found": sorted(identifier_hits - sample_names),
        "region_label_strings_found": sorted(region_labels),
        "region_label_counts": region_counts,
        "direct_gsa_sample_to_region_map_found": bool(sample_names.intersection(identifier_hits)),
        "public_record_search": {
            "direct_crosswalk_found": False,
            "sources": [
                "https://ngdc.cncb.ac.cn/gsa-human/browse/experimentDetail/HRX1592058",
                "https://ngdc.cncb.ac.cn/gsa-human/browse/sampleDetail/HRS1407319",
                "https://pmc.ncbi.nlm.nih.gov/articles/PMC11685534/",
            ],
        },
        "weakest_point": "Source Data workbook and public records show donor-region labels, but no GSA sample name is linked to HBI, LBI, or Tumor. Do not infer GBMxx01/02/03 suffix order.",
    }
    OUTPUT_PATH.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
