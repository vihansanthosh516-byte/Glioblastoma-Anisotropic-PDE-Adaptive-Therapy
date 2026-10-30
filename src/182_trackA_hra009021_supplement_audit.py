"""Audit HRA009021 paper supplement labels against GSA sample metadata."""

from __future__ import annotations

import json
import hashlib
import re
from collections import Counter
from pathlib import Path

from pypdf import PdfReader


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PDF_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "41467_2024_55243_MOESM1_ESM.pdf"
METADATA_AUDIT_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "hra009021_metadata_audit.json"
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "hra009021_supplement_audit.json"


def main() -> None:
    pdf = PdfReader(str(PDF_PATH))
    page_text = [page.extract_text() or "" for page in pdf.pages]
    metadata = json.loads(METADATA_AUDIT_PATH.read_text(encoding="utf-8"))

    # Supplementary Table 4 begins on PDF page 16 and continues on page 17.
    table4_text = "\n".join(page_text[15:17])
    labels = [
        {"region": region, "donor_code": donor}
        for region, donor in re.findall(r"^\s*(HBI|LBI|Tumor)(\d{2})\b", table4_text, flags=re.MULTILINE)
    ]
    region_counts = dict(sorted(Counter(row["region"] for row in labels).items()))
    supplement_donor_codes = sorted({row["donor_code"] for row in labels})

    gsa_identifiers = sorted(
        {
            row["individual_identifier"]
            for row in metadata["sample_crosswalk"]
            if row["individual_identifier"]
        }
    )
    gsa_donor_codes = sorted(
        identifier.removeprefix("GBM")
        for identifier in gsa_identifiers
        if identifier.startswith("GBM") and identifier.removeprefix("GBM").isdigit()
    )
    gsa_sample_names = sorted({row["sample_name"] for row in metadata["sample_crosswalk"]})
    pdf_text = "\n".join(page_text)
    sample_names_in_pdf = [name for name in gsa_sample_names if name in pdf_text]

    result = {
        "dataset": "HRA009021",
        "article_doi": "10.1038/s41467-024-55243-5",
        "source_file": PDF_PATH.name,
        "source_url": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-55243-5/MediaObjects/41467_2024_55243_MOESM1_ESM.pdf",
        "source_file_bytes": PDF_PATH.stat().st_size,
        "source_file_sha256": hashlib.sha256(PDF_PATH.read_bytes()).hexdigest(),
        "pdf_pages": len(page_text),
        "supplementary_table_4_pdf_pages": [16, 17],
        "region_label_source": "Supplement labels HBI/LBI/Tumor samples; article describes region selection from preoperative MRI perfusion and image-guided sampling, with histology checks.",
        "supplementary_table_4_region_labels": labels,
        "region_label_counts": region_counts,
        "supplement_donor_codes": supplement_donor_codes,
        "gsa_donor_identifiers": gsa_identifiers,
        "gsa_donor_codes": gsa_donor_codes,
        "donor_codes_match_exactly": supplement_donor_codes == gsa_donor_codes,
        "gsa_sample_names_found_in_pdf": sample_names_in_pdf,
        "gsa_sample_names_tested": len(gsa_sample_names),
        "direct_gsa_sample_to_region_map_found": bool(sample_names_in_pdf),
        "can_assign_gsa_sample_suffix_order": False,
        "weakest_point": "The supplement labels HBI, LBI, and Tumor samples by donor code, but does not map each GBMxx01/02/03 sample name to a region. Do not guess suffix order.",
        "analysis_status": "metadata audit only; no expression values read or scored",
    }
    OUTPUT_PATH.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
