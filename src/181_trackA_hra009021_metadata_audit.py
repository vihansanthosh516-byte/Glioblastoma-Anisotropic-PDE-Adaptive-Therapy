"""Audit user-provided HRA009021 metadata workbook without reading expression data."""

from __future__ import annotations

import argparse
import hashlib
import json
import warnings
from pathlib import Path
from urllib.parse import urlparse

import openpyxl


PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackA_v2_data_check" / "hra009021_metadata_audit.json"


def read_records(workbook: openpyxl.Workbook, sheet_name: str) -> tuple[list[str], list[dict[str, object]]]:
    sheet = workbook[sheet_name]
    values = list(sheet.iter_rows(values_only=True))
    if not values:
        return [], []
    raw_headers = [str(value).strip() if value is not None else "" for value in values[0]]
    counts: dict[str, int] = {}
    headers = []
    for header in raw_headers:
        if not header:
            headers.append("")
            continue
        counts[header] = counts.get(header, 0) + 1
        headers.append(header if counts[header] == 1 else f"{header} [{counts[header]}]")
    records = []
    for row in values[1:]:
        if not any(value is not None and str(value).strip() for value in row):
            continue
        records.append({header: value for header, value in zip(headers, row) if header})
    return headers, records


def text(value: object) -> str:
    return "" if value is None else str(value).strip()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("workbook", type=Path, help="Local HRA009021.xlsx metadata workbook")
    args = parser.parse_args()
    source = args.workbook.expanduser().resolve(strict=True)

    warnings.filterwarnings("ignore", message="Workbook contains no default style")
    workbook = openpyxl.load_workbook(source, read_only=False, data_only=True)
    sheets = {
        name: read_records(workbook, name)
        for name in ("Individual", "Sample", "Experiment", "Run")
    }
    all_headers = {sheet: result[0] for sheet, result in sheets.items()}
    all_records = {sheet: result[1] for sheet, result in sheets.items()}

    individuals = all_records["Individual"]
    samples = all_records["Sample"]
    experiments = all_records["Experiment"]
    runs = all_records["Run"]
    individual_by_accession = {
        text(record.get("Accession")): record
        for record in individuals
        if text(record.get("Accession"))
    }
    sample_by_accession = {
        text(record.get("Accession")): record
        for record in samples
        if text(record.get("Accession"))
    }
    experiment_by_accession = {
        text(record.get("Accession")): record
        for record in experiments
        if text(record.get("Accession"))
    }
    runs_by_experiment: dict[str, list[dict[str, object]]] = {}
    for record in runs:
        runs_by_experiment.setdefault(text(record.get("Experiment accession")), []).append(record)

    sample_crosswalk = []
    for sample in samples:
        sample_accession = text(sample.get("Accession"))
        individual_accession = text(sample.get("Individual Accession"))
        linked_individual = individual_by_accession.get(individual_accession, {})
        linked_experiments = [
            experiment
            for experiment in experiments
            if text(experiment.get("BioSample accession")) == sample_accession
        ]
        for experiment in linked_experiments:
            experiment_accession = text(experiment.get("Accession"))
            linked_runs = runs_by_experiment.get(experiment_accession, [])
            for run in linked_runs:
                file1_name = text(run.get("Read filename 1"))
                file2_name = text(run.get("Read filename 2"))
                file1_url = text(run.get("DownLoad Read file1"))
                file2_url = text(run.get("DownLoad Read file2"))
                file1_url_name = urlparse(file1_url).path.rsplit("/", 1)[-1]
                file2_url_name = urlparse(file2_url).path.rsplit("/", 1)[-1]
                sample_crosswalk.append(
                    {
                        "individual_accession": individual_accession,
                        "individual_identifier": text(linked_individual.get("Individual Identifier")),
                        "sample_accession": sample_accession,
                        "sample_name": text(sample.get("Sample Name")),
                        "experiment_accession": experiment_accession,
                        "run_accession": text(run.get("Accession")),
                        "read1_filename": file1_name,
                        "read1_url_filename": file1_url_name,
                        "read2_filename": file2_name,
                        "read2_url_filename": file2_url_name,
                    }
                )

    region_terms = ("region", "hbi", "lbi", "core", "edge", "margin", "peritum")
    region_columns = {
        sheet: [header for header in headers if any(term in header.lower() for term in region_terms)]
        for sheet, headers in all_headers.items()
    }
    region_value_hits = sorted(
        {
            term
            for records in all_records.values()
            for record in records
            for value in record.values()
            for term in region_terms
            if term in text(value).lower()
        }
    )
    read2_url_mismatches = [
        row for row in sample_crosswalk if row["read2_filename"] != row["read2_url_filename"]
    ]
    donors_per_sample = {
        text(record.get("Individual Accession"))
        for record in samples
        if text(record.get("Individual Accession"))
    }
    samples_by_donor: dict[str, int] = {}
    for sample in samples:
        donor = text(sample.get("Individual Accession"))
        if donor:
            samples_by_donor[donor] = samples_by_donor.get(donor, 0) + 1

    result = {
        "dataset": "HRA009021",
        "source_file_name": source.name,
        "source_file_bytes": source.stat().st_size,
        "source_file_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "audit_scope": "metadata workbook only; no FASTQ, count matrix, or expression values read",
        "sheet_counts": {
            "individuals": len(individuals),
            "samples": len(samples),
            "experiments": len(experiments),
            "runs": len(runs),
            "crosswalk_rows": len(sample_crosswalk),
        },
        "unique_individual_accessions": len(individual_by_accession),
        "unique_sample_accessions": len(sample_by_accession),
        "unique_experiment_accessions": len(experiment_by_accession),
        "unique_run_accessions": len({text(record.get("Accession")) for record in runs}),
        "samples_per_donor": dict(sorted(samples_by_donor.items())),
        "sample_donor_run_links_complete": (
            len(sample_crosswalk) == len(runs)
            and all(text(row["individual_accession"]) in individual_by_accession for row in sample_crosswalk)
        ),
        "region_label_columns": region_columns,
        "region_label_terms_found_in_values": region_value_hits,
        "region_labels_present": bool(region_value_hits or any(region_columns.values())),
        "read2_filename_url_mismatches": {
            "count": len(read2_url_mismatches),
            "examples": read2_url_mismatches[:3],
            "weakest_point": "The workbook's second-read URL field may be wrong; compare each URL to the public archive listing before using it.",
        },
        "sample_crosswalk": sample_crosswalk,
        "weakest_point": "The workbook links donor, sample, experiment, and run records, but does not identify HBI, LBI, or tumor-core regions.",
    }

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(OUTPUT_PATH)


if __name__ == "__main__":
    main()
