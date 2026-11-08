from __future__ import annotations

import argparse
import json
import platform
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import numpy as np

from src import trackA_tmz_response_score as base


PROJECT_ROOT = base.PROJECT_ROOT
AMENDMENT_PATH = PROJECT_ROOT / "paper" / "trackA_tmz_response_119of121_amendment_2026-11-08.md"
AMENDMENT_FREEZE_PATH = PROJECT_ROOT / "paper" / "trackA_tmz_response_119of121_SHA256SUMS"
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackA_tmz_response_proxy_119of121" / "result.json"
EXPECTED_UNMAPPED_SYMBOLS = ("AF131215.2", "PPIAP22")
EXPECTED_MAPPED_SIZE = 119
EXPECTED_MAPPED_POSITIVE = 47
EXPECTED_MAPPED_NEGATIVE = 72
EXPECTED_LENIENT_LOADER_DONORS = ("GS11", "GS12")
AMENDMENT_SEAL_FILES = (
    "paper/trackA_tmz_response_121_gene_panel_2026-11-07.tsv",
    "paper/trackA_tmz_response_score_protocol_2026-11-07.md",
    "paper/trackA_tmz_response_score_SHA256SUMS",
    "paper/trackA_tmz_response_119of121_amendment_2026-11-08.md",
    "src/trackA_tmz_response_score.py",
    "src/trackA_tmz_response_score_119of121.py",
    "tests/test_trackA_tmz_response_score.py",
    "tests/test_trackA_tmz_response_119of121.py",
    "src/189_trackA_yu2020_lodo.py",
    "paper/trackA_yu2020/PROTOCOL.md",
    "paper/trackA_yu2020/SHA256SUMS",
)


def map_panel_to_matrix(panel: list[dict[str, str]], gene_ids) -> tuple[list[dict[str, str]], tuple[str, ...]]:
    """Keep only exact source symbols authorized by the 119-of-121 amendment."""
    symbols = [str(row["gene_symbol"]).strip() for row in panel]
    if len(symbols) != base.EXPECTED_PANEL_SIZE or len(symbols) != len(set(symbols)):
        raise base.ProxyStop("Amended mapping requires the unchanged 121-gene source panel")
    ids = [str(gene) for gene in gene_ids]
    if len(ids) != len(set(ids)):
        raise base.ProxyStop("Count matrix contains duplicate gene identifiers")
    present = set(ids)
    missing = tuple(sorted(set(symbols) - present))
    if missing != EXPECTED_UNMAPPED_SYMBOLS:
        raise base.ProxyStop(
            "Exact missing symbols differ from the frozen amendment: " + ", ".join(missing)
        )
    mapped = [row for row in panel if row["gene_symbol"] in present]
    positive = sum(float(row["rho_viability_100uM"]) > 0 for row in mapped)
    negative = sum(float(row["rho_viability_100uM"]) < 0 for row in mapped)
    if (len(mapped), positive, negative) != (
        EXPECTED_MAPPED_SIZE, EXPECTED_MAPPED_POSITIVE, EXPECTED_MAPPED_NEGATIVE
    ):
        raise base.ProxyStop(
            "Mapped panel does not match the frozen 119-gene composition: "
            f"{len(mapped)} genes, {positive} positive, {negative} negative"
        )
    return mapped, missing


def verify_amendment_freeze(
    root: Path = PROJECT_ROOT, freeze_path: Path | None = None
) -> dict[str, str]:
    root = Path(root)
    freeze_path = Path(AMENDMENT_FREEZE_PATH if freeze_path is None else freeze_path)
    if not freeze_path.is_file():
        raise base.ProxyStop("Amendment freeze manifest is missing")
    entries: dict[str, str] = {}
    for line in freeze_path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        parts = text.split(maxsplit=1)
        if len(parts) != 2:
            raise base.ProxyStop("Malformed amendment freeze manifest")
        digest, relative_path = parts
        relative_path = relative_path.strip().lstrip("*")
        if relative_path in entries:
            raise base.ProxyStop(f"Duplicate amendment freeze entry: {relative_path}")
        entries[relative_path] = digest.lower()
    expected = set(AMENDMENT_SEAL_FILES)
    if set(entries) != expected:
        missing = sorted(expected - set(entries))
        extra = sorted(set(entries) - expected)
        raise base.ProxyStop(f"Amendment freeze manifest mismatch; missing={missing}; extra={extra}")
    for relative_path, expected_digest in entries.items():
        path = root / relative_path
        if not path.is_file():
            raise base.ProxyStop(f"Amendment-frozen file is missing: {relative_path}")
        if base.sha256_file(path) != expected_digest:
            raise base.ProxyStop(f"Frozen file hash mismatch: {relative_path}")
    return entries


def run_proxy_analysis(output_path: Path = OUTPUT_PATH) -> dict[str, Any]:
    output_path = Path(output_path)
    if output_path.exists():
        raise base.ProxyStop(f"Refusing to overwrite existing result: {output_path}")

    baseline_hashes = base.verify_freeze()
    amendment_hashes = verify_amendment_freeze()
    panel = base.read_panel()
    dataset, donors, flow, input_hashes = base._load_frozen_a5_dataset()
    if tuple(donors) != base.PRIMARY_DONORS:
        raise base.ProxyStop(f"A5 primary donor list differs from the amendment: {donors}")
    loaded_donors = np.asarray(dataset.donors, dtype=str)
    loaded_labels = np.asarray(dataset.labels)
    loaded_counts = np.asarray(dataset.counts)
    if (
        loaded_counts.ndim != 2
        or loaded_donors.ndim != 1
        or loaded_labels.ndim != 1
        or loaded_counts.shape[0] != len(loaded_donors)
        or len(loaded_labels) != len(loaded_donors)
    ):
        raise base.ProxyStop("A5 loader returned arrays with inconsistent cell rows")
    loaded_donor_set = set(loaded_donors.tolist())
    allowed_loader_donors = set(donors) | set(EXPECTED_LENIENT_LOADER_DONORS)
    unexpected_donors = sorted(loaded_donor_set - allowed_loader_donors)
    if unexpected_donors:
        raise base.ProxyStop(f"Unexpected donor IDs from A5 loader: {unexpected_donors}")
    missing_primary_donors = sorted(set(donors) - loaded_donor_set)
    if missing_primary_donors:
        raise base.ProxyStop(f"A5 loader is missing primary donors: {missing_primary_donors}")
    primary_keep = np.isin(loaded_donors, donors)
    primary_donors = loaded_donors[primary_keep]
    if set(primary_donors.tolist()) != set(donors):
        raise base.ProxyStop("Could not isolate the exact six primary donors before scoring")
    mapped_panel, missing = map_panel_to_matrix(panel, dataset.gene_ids)
    result = base.compute_proxy_scores(
        loaded_counts[primary_keep],
        primary_donors,
        loaded_labels[primary_keep],
        dataset.gene_ids,
        mapped_panel,
    )
    result["record_type"] = "exploratory_trackA_tmz_response_direction_score_119_of_121"
    result["panel"] = {
        "n_source_genes": base.EXPECTED_PANEL_SIZE,
        "n_scored_genes": len(mapped_panel),
        "coverage": "119/121 exact matrix symbols",
        "positive_resistance_associated": EXPECTED_MAPPED_POSITIVE,
        "negative_sensitivity_associated": EXPECTED_MAPPED_NEGATIVE,
        "unmapped_source_symbols": list(missing),
        "substitutions": [],
        "directions": "original Ntafoulis in-vitro viability correlation directions",
    }
    result["interpretation"] = (
        "mixed-cell regional expression proxy using 119 of 121 source genes with exact matrix symbols; "
        "the two unmapped genes were omitted without substitution; not a validated predictor or a direct "
        "drug-response measurement; no donor is labeled sensitive or resistant; does not calibrate Track C "
        "or link Track A to Duke"
    )

    package_versions: dict[str, str] = {}
    for package in ("numpy", "scipy", "pandas", "scikit-learn", "torch"):
        try:
            package_versions[package] = version(package)
        except PackageNotFoundError:
            package_versions[package] = "not-installed"
    result["provenance"] = {
        "git_commit": base._git_commit(),
        "baseline_freeze_sha256": base.sha256_file(base.FREEZE_PATH),
        "baseline_frozen_file_sha256": baseline_hashes,
        "amendment_freeze_sha256": base.sha256_file(AMENDMENT_FREEZE_PATH),
        "amendment_frozen_file_sha256": amendment_hashes,
        "amendment_sha256": base.sha256_file(AMENDMENT_PATH),
        "original_panel_sha256": base.sha256_file(base.PANEL_PATH),
        "original_protocol_sha256": base.sha256_file(base.PROTOCOL_PATH),
        "a5_input_hashes": input_hashes,
        "a5_cell_flow": {
            "loader_scope": "The frozen A5 loader may return the six strict primary donors plus GS11/GS12; scoring uses only the six strict primary donors.",
            "loaded_donors": sorted(loaded_donor_set),
            "scored_donors": list(donors),
            "scored_cells": int(primary_keep.sum()),
            "loader_flow": {key: value for key, value in flow.items() if key != "input_hashes"},
        },
        "source_supplement_pmcid": "PMC10575865",
        "source_supplement_sha256": base.SOURCE_SUPPLEMENT_SHA256,
        "source_supplement_md5": base.SOURCE_SUPPLEMENT_MD5,
        "python": sys.version,
        "platform": platform.platform(),
        "package_versions": package_versions,
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if temporary_path.exists():
        raise base.ProxyStop(f"Temporary result path already exists: {temporary_path}")
    temporary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary_path.replace(output_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Exploratory Track A TMZ response-direction score, 119 of 121 genes")
    parser.add_argument("--run", action="store_true", required=True, help="Run once under the frozen amendment")
    args = parser.parse_args()
    if not args.run:
        return 2
    try:
        result = run_proxy_analysis()
    except base.ProxyStop as exc:
        print(f"STOP: {exc}", file=sys.stderr, flush=True)
        return 2
    print(json.dumps({"status": result["status"], "n_donors": result["cohort"]["n_donors"], "result_path": str(OUTPUT_PATH)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
