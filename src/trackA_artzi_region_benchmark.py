"""Locked patient-held-out benchmark for Artzi primary core vs infiltrative edge.

The --preflight path reads HDF5 structure and identifier arrays only. It never
reads an expression/count dataset. Counts are accessed only by --score after
the plan, manifest, and this file have been committed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
import platform
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PLAN_PATH = PROJECT_ROOT / "paper" / "trackA_artzi_preanalysis_2026-11-16.md"
MANIFEST_PATH = PROJECT_ROOT / "paper" / "trackA_artzi_benchmark_manifest.json"
SIGNATURE_DIR = PROJECT_ROOT / "paper" / "trackA_artzi_signatures"
SCRIPT_PATH = Path(__file__).resolve()
ANALYSIS_CODE_PATHS = [SCRIPT_PATH, PROJECT_ROOT / "src" / "trackA_margin_evaluation.py"]
MINIMUM_PATIENTS = 10
MAX_COMPARATORS = 41
ALPHA = 0.05
RANK_GENE_CAP = 30
RANK_PAIR_COUNT = 5
# Source release metadata recorded in paper/trackA_artzi_zenodo_metadata_2026-11-16.md.
EXPECTED_MATRIX_SIZE = 2_994_664_853
EXPECTED_MATRIX_MD5 = "24b06f248245cd6788da20ed0ea7a1fb"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def project_path(path: Path) -> Path:
    """Resolve relative inputs and outputs from the repository root."""
    return path if path.is_absolute() else PROJECT_ROOT / path


def combined_code_hash() -> str:
    digest = hashlib.sha256()
    for path in sorted(ANALYSIS_CODE_PATHS):
        digest.update(str(path.relative_to(PROJECT_ROOT)).replace("\\", "/").encode("utf-8"))
        digest.update(b"\0")
        digest.update(bytes.fromhex(sha256_file(path)))
    return digest.hexdigest().upper()


def _signature_genes() -> dict[str, set[str]]:
    c17 = set()
    with (SIGNATURE_DIR / "c17_bastola_edge_core_markers.tsv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            c17.add(row["gene_symbol"].strip())
    c18 = set()
    with (SIGNATURE_DIR / "c18_5ala_positive_vs_core_leading_edge.tsv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            c18.add(row["gene_symbol"].strip())
    c26_le, c26_core = set(), set()
    with (SIGNATURE_DIR / "c26_saadawy_le_vs_cellular_core.tsv").open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["Signature"].strip().casefold() == "leading edge":
                c26_le.add(row["Symbol"].strip())
            elif row["Signature"].strip().casefold() == "cellular core":
                c26_core.add(row["Symbol"].strip())
    return {"C17": c17, "C18": c18, "C26_LE": c26_le, "C26_CORE": c26_core}


ARM_SPECS: dict[str, dict[str, Any]] = {
    "C01": {"name": "TSP", "eligibility": "portable", "source": "Geman et al.; full scoped inventory and frozen adapter in paper/trackA_prior_art_full_2026-11-16.md and paper/trackA_artzi_preanalysis_2026-11-16.md"},
    "C02": {"name": "k-TSP", "eligibility": "portable", "source": "Tan et al.; full scoped inventory and frozen adapter in paper/trackA_prior_art_full_2026-11-16.md and paper/trackA_artzi_preanalysis_2026-11-16.md"},
    "C03": {"name": "AUCTSP", "eligibility": "portable", "source": "AUCTSP; full scoped inventory and frozen adapter in paper/trackA_prior_art_full_2026-11-16.md and paper/trackA_artzi_preanalysis_2026-11-16.md"},
    "C04": {"name": "Top-scoring triplet", "eligibility": "portable", "source": "Top-scoring triplet; full scoped inventory and frozen adapter in paper/trackA_prior_art_full_2026-11-16.md and paper/trackA_artzi_preanalysis_2026-11-16.md"},
    "C05": {"name": "Rank-pair tree ensemble", "eligibility": "portable", "source": "ranktree-style gene-pair classifier; full scoped inventory and frozen adapter in project Track A plan"},
    "C06": {"name": "GBM rank-pair/Robencla", "eligibility": "pre-score structurally ineligible", "source": "Cross-dataset cluster-replication method; one study available for this benchmark; see prior-art inventory"},
    "C07": {"name": "Penalized pair-order GLMM", "eligibility": "pre-score structurally ineligible", "source": "Between-study random-effects method; one study available for this benchmark; see prior-art inventory"},
    "C08": {"name": "Paired conditional logistic regression", "eligibility": "portable", "source": "Matched-method review and conditional-logistic prior art; see prior-art inventory"},
    "C09": {"name": "Stability-selected conditional logistic regression", "eligibility": "portable", "source": "Stable-CLR and patient-resampled stability-selection prior art; see prior-art inventory"},
    "C10": {"name": "Paired centroid", "eligibility": "portable", "source": "Standard paired-difference classifier; see prior-art inventory"},
    "C11": {"name": "Ridge logistic regression", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C12": {"name": "Linear SVM", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C13": {"name": "Random forest", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C14": {"name": "Paired differential-expression score", "eligibility": "portable", "source": "Paired differential-expression prior art; see prior-art inventory"},
    "C15": {"name": "Ivy GAP regional admixture", "eligibility": "pre-score structurally ineligible", "source": "Ivy GAP four-compartment admixture has no fixed signed P-CT vs P-IT scalar rule; see prior-art inventory"},
    "C16": {"name": "Darmanis infiltration programs", "eligibility": "pre-score structurally ineligible", "source": "Published cell-state/infiltration programs have no fixed pathology-directed predictor; see prior-art inventory"},
    "C17": {"name": "Bastola edge/core markers", "eligibility": "portable if frozen signature coverage passes", "source": "Bastola et al.; frozen source genes in paper/trackA_artzi_signatures/c17_bastola_edge_core_markers.tsv"},
    "C18": {"name": "5-ALA-positive vs core leading edge", "eligibility": "portable if frozen signature coverage passes", "source": "5-ALA regional study; frozen source genes in paper/trackA_artzi_signatures/c18_5ala_positive_vs_core_leading_edge.tsv"},
    "C19": {"name": "Greenwald outer-layer programs", "eligibility": "pre-score structurally ineligible", "source": "Spatial layer/state association has no fixed signed patient-level predictor; see prior-art inventory"},
    "C20": {"name": "Harwood regional modules", "eligibility": "pre-score structurally ineligible", "source": "Regional module/neighborhood analysis has no frozen scalar predictor; see prior-art inventory"},
    "C21": {"name": "Pai infiltration programs", "eligibility": "pre-score structurally ineligible", "source": "Invasion biology study supplies no fixed pathology-directed patient-level score; see prior-art inventory"},
    "C22": {"name": "GBM matrisome regional biology", "eligibility": "pre-score structurally ineligible", "source": "Regional matrix biology has no frozen signed classifier; see prior-art inventory"},
    "C23": {"name": "Petry iron/invasion program", "eligibility": "pre-score structurally ineligible", "source": "Region/functional analysis supplies no frozen expression-only classifier; see prior-art inventory"},
    "C24": {"name": "Sonpatki spatial niche programs", "eligibility": "pre-score structurally ineligible", "source": "Spatial niche method requires unavailable spatial inputs and has no pseudobulk score; see prior-art inventory"},
    "C25": {"name": "Always-tie reference", "eligibility": "portable", "source": "Frozen no-information reference in project Track A plan"},
    "C26": {"name": "Saadawy leading-edge vs cellular-core invasion index", "eligibility": "portable if both frozen signature coverages pass", "source": "Saadawy posted supplement; exact source table in paper/trackA_artzi_signatures/c26_saadawy_2026_supplementary_table_s1.csv"},
    "C27": {"name": "Kersch Ivy-structure L1 logistic classifier", "eligibility": "portable", "source": "Kersch et al.; binary adaptation of source Ivy structure classifier; see prior-art inventory"},
    "C28": {"name": "Sum-zero gene-ratio logistic regression", "eligibility": "portable", "source": "Constrained gene-ratio regression; see prior-art inventory"},
    "C29": {"name": "Gaussian Naive Bayes", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C30": {"name": "Gradient-boosted stumps", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C31": {"name": "BMC/MBMC monotone gene-pair ensemble", "eligibility": "portable adapter", "source": "Fourquet et al. BMC/MBMC; source revision 9b81000bbd6fb48c05b54b5a857f15577dda6a2f; binary monotone logistic adapter, not an exact upstream reproduction"},
    "C32": {"name": "Decision tree", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C33": {"name": "Gaussian process classifier", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C34": {"name": "RBF SVM", "eligibility": "portable", "source": "Standard expression classifier; see prior-art inventory"},
    "C35": {"name": "Pair t-score", "eligibility": "portable", "source": "Pair t-score gene-pair classifier; see 2025 gene-pair review in prior-art inventory"},
    "C36": {"name": "Gene expression ratios (GER)", "eligibility": "portable", "source": "Gene expression ratio method; see 2025 gene-pair review in prior-art inventory"},
    "C37": {"name": "Pairwise SVM ensemble (PSVM)", "eligibility": "portable", "source": "Pairwise SVM ensemble; see 2025 gene-pair review in prior-art inventory"},
    "C38": {"name": "Top-scoring pair of groups (TSPG)", "eligibility": "portable", "source": "TSPG; see 2025 gene-pair review in prior-art inventory"},
    "C39": {"name": "k-TSP plus SVM", "eligibility": "portable", "source": "k-TSP plus machine learning; see 2025 gene-pair review in prior-art inventory"},
    "C40": {"name": "REO plus machine learning", "eligibility": "portable", "source": "Rank-based gene-order plus machine learning; see 2025 gene-pair review in prior-art inventory"},
    "C41": {"name": "TSP plus machine learning", "eligibility": "portable", "source": "TSP plus machine learning; see 2025 gene-pair review in prior-art inventory"},
}


def build_manifest() -> dict[str, Any]:
    required = ["C%02d" % i for i in range(1, MAX_COMPARATORS + 1)]
    if list(ARM_SPECS) != required:
        raise RuntimeError("Comparator specification must contain ordered C01-C41")
    signature_files = [
        "c17_bastola_edge_core_markers.tsv",
        "c18_5ala_positive_vs_core_leading_edge.tsv",
        "c26_saadawy_2026_supplementary_table_s1.csv",
        "c26_saadawy_le_vs_cellular_core.tsv",
    ]
    plan_rows = {}
    for line in PLAN_PATH.read_text(encoding="utf-8").splitlines():
        match = re.match(r"^\|\s*(C\d{2})\s*\|", line)
        if match:
            plan_rows[match.group(1)] = line.strip()
    if set(plan_rows) != set(required):
        raise RuntimeError("Plan must contain exactly the frozen C01-C41 comparator rows")
    signatures = _signature_genes()
    arms = []
    for arm_id, values in ARM_SPECS.items():
        arm = {"id": arm_id, **values}
        arm["frozen_plan_row"] = plan_rows[arm_id]
        if arm_id == "C17":
            with (SIGNATURE_DIR / "c17_bastola_edge_core_markers.tsv").open(encoding="utf-8-sig", newline="") as handle:
                arm["gene_weights"] = {row["gene_symbol"].strip(): int(row["direction"]) for row in csv.DictReader(handle, delimiter="\t")}
        elif arm_id == "C18":
            arm["genes"] = sorted(signatures["C18"])
        elif arm_id == "C26":
            arm["leading_edge_genes"] = sorted(signatures["C26_LE"])
            arm["cellular_core_genes"] = sorted(signatures["C26_CORE"])
        arms.append(arm)
    return {
        "title": "Track A Artzi patient-held-out region benchmark",
        "analysis_label": "exploratory",
        "target": "Primary P-IT vs primary P-CT whole-tissue regional ecosystem; patients are outer-held-out as pairs.",
        "plan_path": str(PLAN_PATH.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "plan_sha256": sha256_file(PLAN_PATH),
        "analysis_code_path": str(SCRIPT_PATH.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "analysis_code_sha256": combined_code_hash(),
        "analysis_code_files": {str(path.relative_to(PROJECT_ROOT)).replace("\\", "/"): sha256_file(path) for path in sorted(ANALYSIS_CODE_PATHS)},
        "minimum_complete_patient_pairs": MINIMUM_PATIENTS,
        "maximum_comparator_arms": MAX_COMPARATORS,
        "alpha": ALPHA,
        "primary_score": "Mean patient-level paired concordance; P-IT higher=1, P-CT higher=0, tie=0.5.",
        "pass_rule": "Strictly exceed every eligible comparator and reject every one-sided exact paired sign test after Holm correction at alpha 0.05.",
        "candidate_method": {
            "name": "Maximin patient-paired rank-pair vote (MP-RPV)",
            "gene_cap": RANK_GENE_CAP,
            "pair_cap": RANK_PAIR_COUNT,
            "pair_orientation": "Maximize the minimum training-patient signed P-IT minus P-CT pair-order contrast; mean contrast then symbol order break ties.",
            "zero_pair_fallback": "Tie scores if no pair has a strictly positive worst-patient contrast.",
            "novelty_status": "Narrow candidate; novelty unresolved and not a claim of global uniqueness.",
        },
        "data": {
            "zenodo_record": "https://zenodo.org/records/22791462",
            "metadata_filename": "Final_meta.csv",
            "matrix_filename": "large_sparse_matrices2.h5",
            "matrix_size_bytes": EXPECTED_MATRIX_SIZE,
            "matrix_release_md5": EXPECTED_MATRIX_MD5,
            "matrix_values_read_by_manifest_generation": False,
        },
        "features": {
            "rank_based_feature_cap": RANK_GENE_CAP,
            "feature_selection": "Training-fold absolute mean paired P-IT minus P-CT log1p(CPM); symbol-order tie break.",
            "nonrank_transform": "log1p(raw_count / specimen_total_count * 10000)",
            "gene_mapping": "Exact released gene-symbol strings only; duplicate or ambiguous assay symbols stop scoring.",
            "signature_files": {name: sha256_file(SIGNATURE_DIR / name) for name in signature_files},
        },
        "outer_split": "Leave one complete patient pair out; both specimens stay in the same fold.",
        "arms": arms,
        "software": {
            "python": platform.python_version(),
            "numpy": _module_version("numpy"),
            "pandas": _module_version("pandas"),
            "scipy": _module_version("scipy"),
            "scikit_learn": _module_version("sklearn"),
            "h5py": _module_version("h5py"),
        },
        "source_revisions": {
            "MBMC": "9b81000bbd6fb48c05b54b5a857f15577dda6a2f",
        },
    }


def _module_version(name: str) -> str:
    try:
        module = __import__(name)
    except ImportError:
        return "not-installed"
    return str(getattr(module, "__version__", "unknown"))


def _verify_frozen_files() -> dict[str, Any]:
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    actual = {
        "plan_sha256": sha256_file(PLAN_PATH),
        "analysis_code_sha256": combined_code_hash(),
    }
    for name, value in actual.items():
        if manifest.get(name) != value:
            raise RuntimeError(f"Frozen file hash mismatch for {name}")
    for name, expected in manifest["features"]["signature_files"].items():
        if sha256_file(SIGNATURE_DIR / name) != expected:
            raise RuntimeError(f"Frozen signature hash mismatch for {name}")
    for name, expected in manifest["analysis_code_files"].items():
        if sha256_file(PROJECT_ROOT / name) != expected:
            raise RuntimeError(f"Frozen analysis-code file hash mismatch for {name}")
    if [arm["id"] for arm in manifest["arms"]] != list(ARM_SPECS):
        raise RuntimeError("Manifest arm order differs from executable arm order")
    return manifest


def _text_dataset(dataset: Any) -> bool:
    import h5py

    return h5py.check_string_dtype(dataset.dtype) is not None or dataset.dtype.kind in {"S", "U", "O"}


def _string_values(dataset: Any) -> list[str]:
    values = dataset[...]
    if getattr(values, "ndim", 0) != 1:
        raise ValueError(f"Identifier dataset must be one-dimensional: {dataset.name}")
    out = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="strict")
        out.append(str(value).strip())
    if any(not value for value in out):
        raise ValueError(f"Blank identifier in {dataset.name}")
    return out


def _identifier_digest(values: list[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest().upper()


def _identifier_kind(path: str) -> str | None:
    clean = path.casefold().rstrip("/")
    leaf = clean.split("/")[-1]
    parts = set(clean.replace("-", "_").split("/"))
    gene_leaves = {"gene", "genes", "gene_name", "gene_names", "gene_symbol", "gene_symbols", "feature", "features", "feature_name", "feature_names", "var_names", "rownames"}
    cell_leaves = {"cell", "cells", "cell_id", "cell_ids", "cellid", "barcodes", "barcode", "obs_names", "colnames"}
    if leaf in gene_leaves or leaf == "_index" and "var" in parts:
        return "gene"
    if leaf in cell_leaves or leaf == "_index" and "obs" in parts:
        return "cell"
    return None


def _read_identifier_candidates(handle: Any) -> dict[str, list[str]]:
    import h5py

    found: dict[str, list[str]] = {"gene": [], "cell": []}

    def visitor(_name: str, obj: Any) -> None:
        if not isinstance(obj, h5py.Dataset) or not _text_dataset(obj):
            return
        kind = _identifier_kind(obj.name)
        if kind:
            found[kind].append(obj.name)

    handle.visititems(visitor)
    return found


def _shape_from_sparse_group(group: Any) -> tuple[int, int] | None:
    if not all(key in group for key in ("data", "indices", "indptr")):
        if not all(key in group for key in ("dims", "i", "p", "x", "rownames", "colnames")):
            return None
        dims = group["dims"][...]
        if getattr(dims, "ndim", 0) != 1 or len(dims) != 2:
            return None
        return tuple(int(value) for value in dims)
    shape = group.attrs.get("shape")
    if shape is None and "shape" in group and hasattr(group["shape"], "__getitem__"):
        shape = group["shape"][...]
    if shape is None:
        return None
    return tuple(int(value) for value in shape)


def _hdf_inventory(handle: Any) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    import h5py

    entries: list[dict[str, Any]] = []
    matrices: list[dict[str, Any]] = []

    def visitor(_name: str, obj: Any) -> None:
        if isinstance(obj, h5py.Group):
            shape = _shape_from_sparse_group(obj)
            if shape is not None:
                kind = "r_csc_block" if all(key in obj for key in ("dims", "i", "p", "x", "rownames", "colnames")) else "sparse_group"
                matrices.append({"path": obj.name, "kind": kind, "shape": list(shape), "keys": sorted(obj.keys())})
            return
        entry = {"path": obj.name, "shape": list(obj.shape), "dtype": str(obj.dtype)}
        entries.append(entry)
        if obj.ndim == 2 and obj.dtype.kind in "biuf":
            matrices.append({"path": obj.name, "kind": "dense_dataset", "shape": list(obj.shape), "dtype": str(obj.dtype)})

    handle.visititems(visitor)
    return entries, matrices


def _metadata_columns(header: list[str]) -> tuple[str, str, str]:
    lookup = {name.strip().casefold().replace(" ", "_"): name for name in header}
    patient = next((lookup[key] for key in ("patient", "patient_id", "patientid") if key in lookup), None)
    setting = next((lookup[key] for key in ("setting", "region", "sample_setting") if key in lookup), None)
    cell_candidates = [lookup[key] for key in ("cell_id", "cellid", "cell", "barcode", "barcodes", "cell_identifier") if key in lookup]
    cell_candidates = list(dict.fromkeys(cell_candidates))
    if not patient or not setting or len(cell_candidates) != 1:
        raise ValueError(f"Need unique Patient, Setting, and cell-ID fields; got header fields {header}")
    return patient, setting, cell_candidates[0]


def _setting(value: str) -> str | None:
    value = value.strip().casefold().replace("_", "-").replace(" ", "")
    aliases = {"p-ct": "P-CT", "pct": "P-CT", "p-it": "P-IT", "pit": "P-IT"}
    return aliases.get(value)


def _read_metadata_for_pairs(path: Path, fields: tuple[str, str, str]) -> tuple[set[str], dict[str, tuple[str, str]]]:
    patient_field, setting_field, cell_field = fields
    core_patients: set[str] = set()
    edge_patients: set[str] = set()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            patient = row.get(patient_field, "").strip()
            region = _setting(row.get(setting_field, ""))
            if patient and region == "P-CT":
                core_patients.add(patient)
            elif patient and region == "P-IT":
                edge_patients.add(patient)
    pairs = core_patients & edge_patients
    if not pairs:
        raise ValueError("Metadata contains no shared primary P-CT/P-IT patient codes")
    cells: dict[str, tuple[str, str]] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            patient = row.get(patient_field, "").strip()
            region = _setting(row.get(setting_field, ""))
            if patient not in pairs or region not in {"P-CT", "P-IT"}:
                continue
            cell_id = row.get(cell_field, "").strip()
            if not cell_id:
                raise ValueError(f"Blank cell identifier for eligible patient {patient} ({region})")
            if cell_id in cells:
                raise ValueError(f"Metadata cell identifier is duplicated: {cell_id}")
            cells[cell_id] = (patient, region)
    return pairs, cells


def _preflight(matrix_path: Path, metadata_path: Path) -> dict[str, Any]:
    import h5py

    manifest = _verify_frozen_files()
    if matrix_path.stat().st_size != EXPECTED_MATRIX_SIZE:
        raise ValueError("Matrix byte size differs from the Zenodo record")
    with metadata_path.open("r", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle))
    fields = _metadata_columns(header)
    patient_pairs, eligible_cells = _read_metadata_for_pairs(metadata_path, fields)
    eligible_cell_ids = set(eligible_cells)

    with h5py.File(matrix_path, "r") as handle:
        datasets, matrices = _hdf_inventory(handle)
        id_paths = _read_identifier_candidates(handle)
        gene_paths, cell_paths = set(id_paths["gene"]), set(id_paths["cell"])
        matrix_blocks = sorted(
            (item for item in matrices if item["kind"] == "r_csc_block"),
            key=lambda item: item["path"],
        )
        if not matrix_blocks:
            raise ValueError("No supported R CSC matrix blocks with dims/i/p/x/rownames/colnames were found")

        blocks = []
        genes_reference: list[str] | None = None
        gene_digest: str | None = None
        matrix_cell_count = 0
        matched_eligible: set[str] = set()
        duplicate_eligible: set[str] = set()
        for matrix in matrix_blocks:
            group_path = matrix["path"]
            gene_path, cell_path = f"{group_path}/rownames", f"{group_path}/colnames"
            if gene_path not in gene_paths or cell_path not in cell_paths:
                raise ValueError(f"Sparse matrix block lacks text rownames/colnames: {group_path}")
            genes = _string_values(handle[gene_path])
            cells = _string_values(handle[cell_path])
            shape = tuple(matrix["shape"])
            if shape != (len(genes), len(cells)):
                raise ValueError(f"R CSC dimensions do not match identifiers in {group_path}")
            if len(set(genes)) != len(genes):
                raise ValueError(f"Duplicate gene identifiers in {gene_path}")
            block_gene_digest = _identifier_digest(genes)
            if genes_reference is None:
                genes_reference = genes
                gene_digest = block_gene_digest
            elif block_gene_digest != gene_digest:
                raise ValueError("R CSC blocks do not share one ordered gene identifier list")

            for cell_id in cells:
                if cell_id in eligible_cell_ids:
                    if cell_id in matched_eligible:
                        duplicate_eligible.add(cell_id)
                    else:
                        matched_eligible.add(cell_id)
            matrix_cell_count += len(cells)
            blocks.append({
                "path": group_path,
                "kind": matrix["kind"],
                "shape": list(shape),
                "orientation": "genes_by_cells",
                "gene_id_path": gene_path,
                "cell_id_path": cell_path,
                "gene_count": len(genes),
                "cell_count": len(cells),
                "gene_ids_sha256": block_gene_digest,
                "cell_ids_sha256": _identifier_digest(cells),
            })
    if genes_reference is None or gene_digest is None:
        raise ValueError("No gene identifiers were available in supported matrix blocks")
    genes = genes_reference
    missing_count = len(eligible_cell_ids - matched_eligible)

    donor_regions: dict[str, Counter[str]] = {patient: Counter() for patient in patient_pairs}
    for cell_id, (patient, region) in eligible_cells.items():
        donor_regions[patient][region] += 1
    incomplete = {patient: dict(counts) for patient, counts in donor_regions.items() if not counts["P-CT"] or not counts["P-IT"]}
    signatures = _signature_genes()
    coverage = {}
    gene_set = set(genes)
    for key, source in signatures.items():
        coverage[key] = {"available": len(source & gene_set), "required": len(source), "fraction": len(source & gene_set) / max(1, len(source))}
    matrix_layout = {
        "kind": "r_csc_blocks",
        "blocks": blocks,
        "block_count": len(blocks),
        "shape": [len(genes), matrix_cell_count],
        "orientation": "genes_by_cells",
        "gene_ids_sha256": gene_digest,
    }
    join_ok = not incomplete and not missing_count and not duplicate_eligible
    count_gate = len(patient_pairs) >= MINIMUM_PATIENTS
    return {
        "analysis_label": "exploratory",
        "mode": "metadata_only",
        "expression_values_opened": False,
        "scores_calculated": False,
        "plan_sha256": manifest["plan_sha256"],
        "analysis_code_sha256": manifest["analysis_code_sha256"],
        "manifest_sha256": sha256_file(MANIFEST_PATH),
        "matrix_file": {"name": matrix_path.name, "size_bytes": matrix_path.stat().st_size, "expected_release_md5": EXPECTED_MATRIX_MD5},
        "metadata_file": {"name": metadata_path.name, "size_bytes": metadata_path.stat().st_size, "sha256": sha256_file(metadata_path)},
        "hdf5_inventory": datasets,
        "matrix_layout": matrix_layout,
        "gene_id_path": blocks[0]["gene_id_path"],
        "gene_count": len(genes),
        "matrix_cell_count": matrix_cell_count,
        "metadata_columns": {"patient": fields[0], "setting": fields[1], "cell_id": fields[2]},
        "eligible_patient_codes": sorted(patient_pairs),
        "eligible_patient_count": len(patient_pairs),
        "eligible_cell_count": len(eligible_cells),
        "eligible_cells_per_patient_region": {patient: dict(counts) for patient, counts in sorted(donor_regions.items())},
        "exact_cell_join": join_ok,
        "matched_eligible_cell_count": len(matched_eligible),
        "missing_metadata_cell_count": missing_count,
        "duplicate_eligible_matrix_cell_count": len(duplicate_eligible),
        "incomplete_patient_region_counts": incomplete,
        "signature_coverage": coverage,
        "minimum_patient_gate_passed": count_gate,
        "score_gate_passed": bool(join_ok and count_gate),
        "status": "metadata_preflight_pass" if join_ok and count_gate else "metadata_gate_stop",
        "next_step": "If score_gate_passed is true, commit this metadata-only audit before running --score; otherwise stop and report the gate failure.",
    }


def _read_counts_by_patient_region(matrix_path: Path, metadata_path: Path, preflight: dict[str, Any]):
    """Read selected raw-count cells from locked R CSC blocks into patient pseudobulks."""
    import h5py
    import numpy as np
    from scipy.sparse import csc_matrix

    layout = preflight["matrix_layout"]
    if layout.get("kind") != "r_csc_blocks" or not layout.get("blocks"):
        raise RuntimeError("Metadata-only preflight did not lock supported R CSC blocks")
    with h5py.File(matrix_path, "r") as handle:
        genes = _string_values(handle[preflight["gene_id_path"]])
        if len(genes) != preflight["gene_count"] or _identifier_digest(genes) != layout["gene_ids_sha256"]:
            raise RuntimeError("HDF5 gene identifiers changed after metadata-only preflight")
        fields = tuple(preflight["metadata_columns"][key] for key in ("patient", "setting", "cell_id"))
        patients, cell_rows = _read_metadata_for_pairs(metadata_path, fields)
        ordered_patients = sorted(patients)
        if ordered_patients != preflight["eligible_patient_codes"]:
            raise RuntimeError("Eligible patient codes changed after metadata-only preflight")
        patient_index = {patient: index for index, patient in enumerate(ordered_patients)}
        group_count, gene_count = len(ordered_patients) * 2, len(genes)
        totals = np.zeros((group_count, gene_count), dtype=np.float64)
        matched_cells: set[str] = set()
        for block_spec in layout["blocks"]:
            target = handle[block_spec["path"]]
            if not isinstance(target, h5py.Group) or not all(key in target for key in ("dims", "i", "p", "x", "rownames", "colnames")):
                raise RuntimeError(f"R CSC block changed after metadata-only preflight: {block_spec['path']}")
            cells = _string_values(handle[block_spec["cell_id_path"]])
            shape = tuple(int(value) for value in target["dims"][...])
            if shape != tuple(block_spec["shape"]) or len(cells) != block_spec["cell_count"]:
                raise RuntimeError(f"R CSC dimensions or cell identifiers changed: {block_spec['path']}")
            if _identifier_digest(cells) != block_spec["cell_ids_sha256"]:
                raise RuntimeError(f"R CSC cell identifiers changed: {block_spec['path']}")

            group_codes = np.full(len(cells), -1, dtype=np.int32)
            for cell_index, cell_id in enumerate(cells):
                pair = cell_rows.get(cell_id)
                if pair is None:
                    continue
                if cell_id in matched_cells:
                    raise RuntimeError(f"Eligible metadata cell occurs more than once in matrix blocks: {cell_id}")
                matched_cells.add(cell_id)
                patient, region = pair
                region_offset = 0 if region == "P-CT" else 1
                group_codes[cell_index] = patient_index[patient] * 2 + region_offset

            if shape != (gene_count, len(cells)):
                raise RuntimeError(f"R CSC block is not genes-by-cells: {block_spec['path']}")
            values_ds, indices_ds, indptr_ds = target["x"], target["i"], target["p"]
            if len(indptr_ds) != len(cells) + 1 or len(values_ds) != len(indices_ds):
                raise RuntimeError(f"R CSC pointer/value lengths are inconsistent: {block_spec['path']}")
            if int(indptr_ds[-1]) != len(values_ds):
                raise RuntimeError(f"R CSC final pointer differs from sparse payload length: {block_spec['path']}")
            for start in range(0, len(cells), 20000):
                stop = min(len(cells), start + 20000)
                chunk_codes = group_codes[start:stop]
                if not np.any(chunk_codes >= 0):
                    continue
                ptr = np.asarray(indptr_ds[start : stop + 1], dtype=np.int64)
                if len(ptr) != stop - start + 1 or (np.diff(ptr) < 0).any():
                    raise RuntimeError(f"R CSC pointers are invalid: {block_spec['path']}")
                value_start, value_stop = int(ptr[0]), int(ptr[-1])
                if value_start < 0 or value_stop > len(values_ds):
                    raise RuntimeError(f"R CSC pointers exceed sparse payload length: {block_spec['path']}")
                vals = np.asarray(values_ds[value_start:value_stop], dtype=np.float64)
                indices = np.asarray(indices_ds[value_start:value_stop], dtype=np.int64)
                local_ptr = ptr - value_start
                if len(vals) != len(indices) or (indices < 0).any() or (indices >= gene_count).any():
                    raise RuntimeError(f"R CSC indices do not match the gene dimension: {block_spec['path']}")
                if not np.isfinite(vals).all() or (vals < 0).any() or not np.equal(vals, np.floor(vals)).all():
                    raise ValueError("Raw matrix contains non-finite, negative, or non-integer values")
                sparse_block = csc_matrix((vals, indices, local_ptr), shape=(gene_count, stop - start))
                for code in np.unique(chunk_codes[chunk_codes >= 0]):
                    totals[code] += np.asarray(sparse_block[:, chunk_codes == code].sum(axis=1)).ravel()
        if matched_cells != set(cell_rows):
            raise RuntimeError("Some eligible metadata cells were absent from matrix blocks after preflight")

    names = [(patient, region) for patient in ordered_patients for region in ("P-CT", "P-IT")]
    if not np.isfinite(totals).all() or (totals < 0).any():
        raise ValueError("Aggregated counts are non-finite or negative")
    empty = [names[index] for index, row in enumerate(totals) if row.sum() <= 0]
    if empty:
        raise ValueError(f"Empty regional pseudobulk after count aggregation: {empty}")
    return ordered_patients, names, genes, totals


def _log_cpm(counts):
    import numpy as np

    totals = counts.sum(axis=1)
    if (totals <= 0).any() or not np.isfinite(totals).all():
        raise ValueError("Every pseudobulk must have a positive finite library total")
    values = np.log1p(counts / totals[:, None] * 10000.0)
    if not np.isfinite(values).all():
        raise ValueError("log1p(CPM) transform produced non-finite values")
    return values


def _stack(core, edge):
    import numpy as np

    return np.vstack((core, edge)), np.r_[np.zeros(len(core), dtype=int), np.ones(len(edge), dtype=int)]


def _top_genes(core, edge, genes, count):
    import numpy as np

    delta = (edge - core).mean(axis=0)
    order = np.lexsort((np.asarray(genes, dtype=str), -np.abs(delta)))
    return order[: min(count, len(order))]


def _patient_concordance(core_scores, edge_scores):
    import numpy as np

    return (edge_scores > core_scores).astype(float) + 0.5 * (edge_scores == core_scores)


def _score_estimator(model, train_x, train_y, test_x):
    import numpy as np

    model.fit(train_x, train_y)
    if hasattr(model, "decision_function"):
        result = np.asarray(model.decision_function(test_x), dtype=float).reshape(-1)
    elif hasattr(model, "predict_proba"):
        result = np.asarray(model.predict_proba(test_x)[:, 1], dtype=float)
    else:
        result = np.asarray(model.predict(test_x), dtype=float).reshape(-1)
    if not np.isfinite(result).all():
        raise ValueError("Classifier produced non-finite prediction scores")
    return result


def _pair_arrays(indices, genes):
    import numpy as np

    indices = np.asarray(indices, dtype=int)
    local_genes = [genes[index] for index in indices]
    sort_order = np.argsort(np.asarray(local_genes, dtype=str), kind="stable")
    left, right = np.triu_indices(len(indices), k=1)
    return indices[sort_order[left]], indices[sort_order[right]]


def _pair_votes(values, first, second):
    import numpy as np

    return np.sign(values[:, first] - values[:, second])


def _tsp_statistics(core, edge, genes, pool):
    import numpy as np

    first, second = _pair_arrays(pool, genes)
    x, y = _stack(core, edge)
    votes = _pair_votes(x, first, second)
    class_difference = votes[y == 1].mean(axis=0) - votes[y == 0].mean(axis=0)
    return first, second, class_difference


def _rank_disjoint_pairs(first, second, values, genes, maximum):
    order = sorted(range(len(values)), key=lambda index: (-float(values[index]), genes[first[index]], genes[second[index]]))
    selected, used = [], set()
    for index in order:
        a, b = genes[first[index]], genes[second[index]]
        if a in used or b in used:
            continue
        selected.append(index)
        used.update((a, b))
        if len(selected) == maximum:
            break
    return selected


def _fit_conditional_beta(core, edge, *, l1=False, penalty=1.0):
    import numpy as np
    from scipy.optimize import minimize
    from scipy.special import expit

    difference = edge - core
    difference = difference[:, np.any(difference != 0, axis=0)]
    if difference.shape[1] == 0:
        return np.zeros(edge.shape[1]), np.zeros(edge.shape[1], dtype=bool)
    if l1:
        def objective(beta):
            margins = difference @ beta
            loss = np.logaddexp(0.0, -margins).sum() + penalty * np.abs(beta).sum()
            gradient = -(difference.T @ expit(-margins))
            # L-BFGS uses a deterministic subgradient at zero.
            gradient += penalty * np.sign(beta)
            return float(loss), gradient
    else:
        def objective(beta):
            margins = difference @ beta
            loss = np.logaddexp(0.0, -margins).sum() + 0.5 * penalty * np.dot(beta, beta)
            gradient = -(difference.T @ expit(-margins)) + penalty * beta
            return float(loss), gradient
    result = minimize(objective, np.zeros(difference.shape[1]), jac=True, method="L-BFGS-B", options={"maxiter": 2000, "ftol": 1e-12})
    if not result.success and not np.isfinite(result.fun):
        raise RuntimeError(f"Conditional logistic optimizer failed: {result.message}")
    keep = np.any((edge - core) != 0, axis=0)
    beta = np.zeros(edge.shape[1])
    beta[keep] = result.x
    return beta, np.abs(beta) > 1e-8


def _candidate_scores(core, edge, test_core, test_edge, genes):
    import numpy as np

    pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
    first, second = _pair_arrays(pool, genes)
    if len(first) == 0:
        return 0.0, 0.0, {"selected_pairs": []}
    delta = _pair_votes(edge, first, second) - _pair_votes(core, first, second)
    options = []
    for index in range(len(first)):
        positive = float(np.min(delta[:, index]))
        negative = float(np.min(-delta[:, index]))
        if positive > negative:
            orientation, worst = 1.0, positive
        elif negative > positive:
            orientation, worst = -1.0, negative
        else:
            orientation, worst = (1.0 if float(delta[:, index].mean()) >= 0 else -1.0), positive
        mean_oriented = float(np.mean(orientation * delta[:, index]))
        options.append((worst, mean_oriented, genes[first[index]], genes[second[index]], index, orientation))
    options.sort(key=lambda item: (-item[0], -item[1], item[2], item[3]))
    selected = [item for item in options if item[0] > 0][:RANK_PAIR_COUNT]
    if not selected:
        return 0.0, 0.0, {"selected_pairs": []}
    indices = [item[4] for item in selected]
    directions = np.asarray([item[5] for item in selected])
    core_vote = _pair_votes(test_core[None, :], first[indices], second[indices])[0]
    edge_vote = _pair_votes(test_edge[None, :], first[indices], second[indices])[0]
    info = {"selected_pairs": [{"gene_1": item[2], "gene_2": item[3], "orientation": item[5], "worst_training_contrast": item[0], "mean_training_contrast": item[1]} for item in selected]}
    return float(np.mean(core_vote * directions)), float(np.mean(edge_vote * directions)), info


def _percentile_scores(row, genes):
    from scipy.stats import rankdata

    return rankdata(row, method="average") / len(genes)


def _signature_weights(genes, arm_id):
    import numpy as np

    weights = np.zeros(len(genes), dtype=float)
    gene_index = {gene: index for index, gene in enumerate(genes)}
    signatures = _signature_genes()
    if arm_id == "C17":
        with (SIGNATURE_DIR / "c17_bastola_edge_core_markers.tsv").open(encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                gene = row["gene_symbol"].strip()
                if gene in gene_index:
                    weights[gene_index[gene]] = float(row["direction"])
    elif arm_id == "C18":
        for gene in signatures["C18"]:
            if gene in gene_index:
                weights[gene_index[gene]] = 1.0
    elif arm_id == "C26":
        for gene in signatures["C26_LE"]:
            if gene in gene_index:
                weights[gene_index[gene]] = 1.0
        for gene in signatures["C26_CORE"]:
            if gene in gene_index:
                weights[gene_index[gene]] = -1.0
    return weights


def _signature_score(row, genes, arm_id):
    import numpy as np

    ranks = _percentile_scores(row, genes)
    weights = _signature_weights(genes, arm_id)
    if arm_id == "C26":
        positive = weights > 0
        negative = weights < 0
        return float(ranks[positive].mean() - ranks[negative].mean())
    positive = weights > 0
    negative = weights < 0
    if positive.any() and negative.any():
        return float(ranks[positive].mean() - ranks[negative].mean())
    return float(ranks[positive].mean())


def _scores_from_pair_selection(core, edge, test_core, test_edge, genes, pool, maximum, criterion="tsp"):
    import numpy as np

    first, second, difference = _tsp_statistics(core, edge, genes, pool)
    values = np.abs(difference)
    if criterion == "auc":
        from sklearn.metrics import roc_auc_score

        x, y = _stack(core, edge)
        votes = _pair_votes(x, first, second)
        aucs = np.asarray([roc_auc_score(y, votes[:, index]) for index in range(votes.shape[1])])
        values = np.maximum(aucs, 1 - aucs)
        directions = np.where(aucs >= 0.5, 1.0, -1.0)
    elif criterion == "t":
        x, y = _stack(core, edge)
        votes = _pair_votes(x, first, second)
        directions = np.where(difference >= 0, 1.0, -1.0)
        values = np.asarray([_welch_abs_t(votes[y == 1, j], votes[y == 0, j]) for j in range(votes.shape[1])])
    else:
        directions = np.where(difference >= 0, 1.0, -1.0)
    selected = _rank_disjoint_pairs(first, second, values, genes, maximum)
    if not selected:
        return 0.0, 0.0, []
    core_votes = _pair_votes(test_core[None, :], first[selected], second[selected])[0] * directions[selected]
    edge_votes = _pair_votes(test_edge[None, :], first[selected], second[selected])[0] * directions[selected]
    details = [(genes[first[index]], genes[second[index]], float(directions[index]), float(values[index])) for index in selected]
    return float(core_votes.mean()), float(edge_votes.mean()), details


def _welch_abs_t(left, right):
    import numpy as np

    left, right = np.asarray(left, dtype=float), np.asarray(right, dtype=float)
    difference = float(left.mean() - right.mean())
    variance = float(left.var(ddof=1) / len(left) + right.var(ddof=1) / len(right)) if len(left) > 1 and len(right) > 1 else 0.0
    if variance <= 0:
        return math.inf if difference != 0 else 0.0
    return abs(difference) / math.sqrt(variance)


def _linear_model_score(arm_id, core, edge, test_core, test_edge, genes, donor_ids):
    import numpy as np
    from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
    from sklearn.gaussian_process import GaussianProcessClassifier
    from sklearn.gaussian_process.kernels import RBF
    from sklearn.linear_model import LogisticRegression
    from sklearn.naive_bayes import GaussianNB
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC
    from sklearn.tree import DecisionTreeClassifier

    train_x, train_y = _stack(core, edge)
    test_x = np.vstack((test_core, test_edge))
    info: dict[str, Any] = {}

    if arm_id in {"C01", "C02", "C03", "C36"}:
        pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
        maximum = 1 if arm_id in {"C01", "C03"} else RANK_PAIR_COUNT
        criterion = "auc" if arm_id == "C03" else "t" if arm_id == "C36" else "tsp"
        cs, es, selected = _scores_from_pair_selection(core, edge, test_core, test_edge, genes, pool, maximum, criterion)
        info["selected_pairs"] = selected
        return cs, es, info

    if arm_id == "C04":
        return _triplet_model_score(core, edge, test_core, test_edge, genes, donor_ids)

    if arm_id == "C05":
        pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
        first, second = _pair_arrays(pool, genes)
        train_features = _pair_votes(train_x, first, second)
        test_features = _pair_votes(test_x, first, second)
        model = RandomForestClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced", random_state=0, n_jobs=1)
        scores = _score_estimator(model, train_features, train_y, test_features)
        info["pair_feature_count"] = int(train_features.shape[1])
        return float(scores[0]), float(scores[1]), info

    if arm_id in {"C08", "C09"}:
        pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
        core_pool, edge_pool = core[:, pool], edge[:, pool]
        if arm_id == "C08":
            beta, _ = _fit_conditional_beta(core_pool, edge_pool, l1=False, penalty=1.0)
            info["selected_gene_count"] = int(np.count_nonzero(beta))
            return float(test_core[pool] @ beta), float(test_edge[pool] @ beta), info
        else:
            rng = np.random.default_rng(0)
            selected_count = np.zeros(len(pool), dtype=int)
            for _ in range(100):
                sample = rng.integers(0, len(core), size=len(core))
                _, selected = _fit_conditional_beta(core_pool[sample], edge_pool[sample], l1=True, penalty=1.0)
                selected_count += selected
            keep = selected_count >= math.ceil(0.75 * 100)
            info["selected_gene_count"] = int(keep.sum())
            if not keep.any():
                return 0.0, 0.0, info
            beta, _ = _fit_conditional_beta(core_pool[:, keep], edge_pool[:, keep], l1=False, penalty=1.0)
            return float(test_core[pool[keep]] @ beta), float(test_edge[pool[keep]] @ beta), info

    if arm_id == "C10":
        weights = (edge - core).mean(axis=0)
        return float(test_core @ weights), float(test_edge @ weights), info

    if arm_id in {"C11", "C12", "C13", "C27"}:
        if arm_id == "C11":
            model = make_pipeline(StandardScaler(), LogisticRegression(C=1.0, class_weight="balanced", solver="lbfgs", max_iter=2000, random_state=0))
        elif arm_id == "C12":
            model = make_pipeline(StandardScaler(), SVC(kernel="linear", C=1.0, class_weight="balanced"))
        elif arm_id == "C13":
            model = RandomForestClassifier(n_estimators=500, max_features="sqrt", class_weight="balanced", random_state=0, n_jobs=1)
        else:
            model = make_pipeline(StandardScaler(), LogisticRegression(penalty="l1", C=1.0, solver="liblinear", class_weight="balanced", max_iter=5000, random_state=0))
        scores = _score_estimator(model, train_x, train_y, test_x)
        return float(scores[0]), float(scores[1]), info

    if arm_id == "C14":
        delta = (edge - core).mean(axis=0)
        order = np.lexsort((np.asarray(genes, dtype=str), -np.abs(delta)))[: min(50, len(genes))]
        weights = np.sign(delta[order])
        info["selected_gene_count"] = len(order)
        return float(test_core[order] @ weights), float(test_edge[order] @ weights), info

    if arm_id in {"C17", "C18", "C26"}:
        cs = _signature_score(test_core, genes, arm_id)
        es = _signature_score(test_edge, genes, arm_id)
        return cs, es, info

    if arm_id == "C25":
        return 0.0, 0.0, info

    pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
    if arm_id == "C28":
        train_x = train_x[:, pool] - train_x[:, pool].mean(axis=1, keepdims=True)
        test_x = test_x[:, pool] - test_x[:, pool].mean(axis=1, keepdims=True)
        model = LogisticRegression(C=1.0, solver="lbfgs", max_iter=2000, random_state=0)
    elif arm_id == "C29":
        train_x, test_x = train_x[:, pool], test_x[:, pool]
        model = GaussianNB(var_smoothing=1e-9)
    elif arm_id == "C30":
        train_x, test_x = train_x[:, pool], test_x[:, pool]
        model = GradientBoostingClassifier(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=0)
    elif arm_id == "C32":
        train_x, test_x = train_x[:, pool], test_x[:, pool]
        model = DecisionTreeClassifier(max_depth=3, random_state=0)
    elif arm_id == "C33":
        train_x, test_x = train_x[:, pool], test_x[:, pool]
        model = make_pipeline(StandardScaler(), GaussianProcessClassifier(kernel=RBF(length_scale=1.0, length_scale_bounds="fixed"), optimizer=None, random_state=0))
    elif arm_id == "C34":
        train_x, test_x = train_x[:, pool], test_x[:, pool]
        model = make_pipeline(StandardScaler(), SVC(kernel="rbf", C=1.0, gamma="scale", class_weight="balanced"))
    elif arm_id == "C35":
        return _pair_tscore_model(core, edge, test_core, test_edge, genes, pool)
    elif arm_id == "C36":
        raise AssertionError("C36 handled above")
    elif arm_id == "C37":
        return _pairwise_svm_model(core, edge, test_core, test_edge, genes, pool, donor_ids)
    elif arm_id == "C38":
        return _tspg_model(core, edge, test_core, test_edge, genes, donor_ids)
    elif arm_id == "C39":
        return _ktsp_svm_model(core, edge, test_core, test_edge, genes, donor_ids)
    elif arm_id == "C40":
        return _reo_svm_model(core, edge, test_core, test_edge, genes, donor_ids)
    elif arm_id == "C41":
        return _tsp_ml_model(core, edge, test_core, test_edge, genes, donor_ids)
    elif arm_id == "C31":
        return _monotone_pair_model(core, edge, test_core, test_edge, genes, pool, donor_ids)
    else:
        raise KeyError(f"No executable adapter for {arm_id}")
    scores = _score_estimator(model, train_x, train_y, test_x)
    return float(scores[0]), float(scores[1]), info


def _triplet_patterns(values, triplets):
    import numpy as np

    permutation_ids = {permutation: index for index, permutation in enumerate(itertools.permutations(range(3)))}
    patterns = np.full((len(values), len(triplets)), -1, dtype=np.int8)
    for column, triple in enumerate(triplets):
        block = values[:, triple]
        orders = np.argsort(block, axis=1, kind="stable")
        tied = (block[:, 0] == block[:, 1]) | (block[:, 0] == block[:, 2]) | (block[:, 1] == block[:, 2])
        for row_index, order in enumerate(orders):
            if not tied[row_index]:
                patterns[row_index, column] = permutation_ids[tuple(int(value) for value in order)]
    return patterns


def _triplet_fit(core, edge, test_core, test_edge, genes, pool, number):
    import numpy as np

    x, y = _stack(core, edge)
    pool = sorted(pool, key=lambda index: genes[index])
    triples = list(itertools.combinations(pool, 3))
    if not triples:
        return 0.0, 0.0, []
    patterns = _triplet_patterns(x, triples)
    category_delta = np.zeros((6, len(triples)))
    for category in range(6):
        chosen = patterns == category
        category_delta[category] = chosen[y == 1].mean(axis=0) - chosen[y == 0].mean(axis=0)
    best_category = np.argmax(np.abs(category_delta), axis=0)
    strengths = np.abs(category_delta[best_category, np.arange(len(triples))])
    directions = np.where(category_delta[best_category, np.arange(len(triples))] >= 0, 1.0, -1.0)
    ranking = sorted(range(len(triples)), key=lambda i: (-float(strengths[i]), tuple(genes[g] for g in triples[i])))
    selected = ranking[: min(number, len(ranking))]
    test_patterns = _triplet_patterns(np.vstack((test_core, test_edge)), [triples[index] for index in selected])
    votes = np.zeros_like(test_patterns, dtype=float)
    for column, index in enumerate(selected):
        valid = test_patterns[:, column] >= 0
        votes[valid, column] = np.where(test_patterns[valid, column] == best_category[index], directions[index], -directions[index])
    score = votes.mean(axis=1) if votes.shape[1] else np.zeros(2)
    details = [{"genes": [genes[g] for g in triples[index]], "category": int(best_category[index]), "direction": float(directions[index]), "strength": float(strengths[index])} for index in selected]
    return float(score[0]), float(score[1]), details


def _triplet_model_score(core, edge, test_core, test_edge, genes, donor_ids):
    import numpy as np

    candidates = [2, 3, 4]
    inner = {number: [] for number in candidates}
    for heldout in range(len(core)):
        train = np.asarray([index for index in range(len(core)) if index != heldout], dtype=int)
        if len(train) < 2:
            continue
        pool = _top_genes(core[train], edge[train], genes, RANK_GENE_CAP)
        for number in candidates:
            cs, es, _ = _triplet_fit(core[train], edge[train], core[heldout], edge[heldout], genes, pool, number)
            inner[number].append(float(es > cs) + 0.5 * float(es == cs))
    chosen = min(candidates, key=lambda number: (-float(np.mean(inner[number])) if inner[number] else math.inf, number))
    pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
    cs, es, details = _triplet_fit(core, edge, test_core, test_edge, genes, pool, chosen)
    return cs, es, {"chosen_triplet_count": chosen, "inner_mean_concordance": {str(key): float(np.mean(value)) if value else None for key, value in inner.items()}, "selected_triplets": details}


def _fit_monotone_logistic(x, y, signs, penalty=1.0):
    import numpy as np
    from scipy.optimize import minimize
    from scipy.special import expit
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler().fit(x)
    scaled = scaler.transform(x)
    labels = np.asarray(y, dtype=float)
    signed = labels * 2.0 - 1.0

    def objective(parameters):
        weights, intercept = parameters[:-1], parameters[-1]
        margin = signed * (scaled @ weights + intercept)
        loss = np.logaddexp(0.0, -margin).sum() + 0.5 * penalty * np.dot(weights, weights)
        residual = -signed * expit(-margin)
        gradient = np.r_[scaled.T @ residual + penalty * weights, residual.sum()]
        return float(loss), gradient

    bounds = [(0.0, None) if sign > 0 else (None, 0.0) for sign in signs] + [(None, None)]
    result = minimize(objective, np.zeros(3), jac=True, bounds=bounds, method="L-BFGS-B", options={"maxiter": 2000, "ftol": 1e-12})
    if not np.isfinite(result.fun):
        raise RuntimeError("Monotone logistic optimization produced a non-finite loss")
    return scaler, result.x[:-1], float(result.x[-1])


def _monotone_probability(model, x):
    import numpy as np
    from scipy.special import expit

    scaler, weights, intercept = model
    return expit(scaler.transform(x) @ weights + intercept)


def _monotone_pair_model(core, edge, test_core, test_edge, genes, pool, donor_ids):
    import numpy as np

    first, second = _pair_arrays(pool, genes)
    directions = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
    patient_count = len(core)
    pair_errors = np.zeros((len(first), len(directions)), dtype=float)
    for heldout in range(patient_count):
        train_ids = [index for index in range(patient_count) if index != heldout]
        train_x, train_y = _stack(core[train_ids], edge[train_ids])
        test_x, test_y = _stack(core[heldout : heldout + 1], edge[heldout : heldout + 1])
        for pair_index, (g1, g2) in enumerate(zip(first, second)):
            pair_train = train_x[:, [g1, g2]]
            pair_test = test_x[:, [g1, g2]]
            for direction_index, signs in enumerate(directions):
                model = _fit_monotone_logistic(pair_train, train_y, signs)
                probability = _monotone_probability(model, pair_test)
                pair_errors[pair_index, direction_index] += float(np.abs((probability >= 0.5).astype(float) - test_y).mean()) / patient_count
    best_direction = np.argmin(pair_errors, axis=1)
    best_error = pair_errors[np.arange(len(first)), best_direction]
    ordered = sorted(range(len(first)), key=lambda i: (float(best_error[i]), genes[first[i]], genes[second[i]]))
    selected, used = [], set()
    for index in ordered:
        a, b = genes[first[index]], genes[second[index]]
        if a in used or b in used:
            continue
        selected.append(index)
        used.update((a, b))
        if len(selected) == 5:
            break
    if not selected:
        return 0.0, 0.0, {"selected_pairs": []}
    test_x = np.vstack((test_core, test_edge))
    scores = np.zeros((2, len(selected)))
    for selected_index, pair_index in enumerate(selected):
        g1, g2 = first[pair_index], second[pair_index]
        train_x, train_y = _stack(core, edge)
        model = _fit_monotone_logistic(train_x[:, [g1, g2]], train_y, directions[best_direction[pair_index]])
        scores[:, selected_index] = _monotone_probability(model, test_x[:, [g1, g2]])
    details = [{"genes": [genes[first[index]], genes[second[index]]], "directions": directions[int(best_direction[index])], "inner_mae": float(best_error[index])} for index in selected]
    return float(scores[0].mean()), float(scores[1].mean()), {"selected_pairs": details}


def _pair_tscore_model(core, edge, test_core, test_edge, genes, pool):
    import numpy as np

    x, y = _stack(core, edge)
    first, second = _pair_arrays(pool, genes)
    ranked, axes, means, scales, projection_centers, projection_scales = [], [], [], [], [], []
    for index, (g1, g2) in enumerate(zip(first, second)):
        pair = x[:, [g1, g2]]
        center = pair.mean(axis=0)
        scale = pair.std(axis=0, ddof=1)
        scale[scale == 0] = 1.0
        z = (pair - center) / scale
        axis = z[y == 1].mean(axis=0) - z[y == 0].mean(axis=0)
        norm = np.linalg.norm(axis)
        if norm == 0:
            strength = 0.0
            axis = np.zeros(2)
        else:
            axis = axis / norm
            projection = z @ axis
            strength = _welch_abs_t(projection[y == 1], projection[y == 0])
        projection_center = float((z @ axis).mean())
        projection_scale = float((z @ axis).std(ddof=1))
        if projection_scale == 0:
            projection_scale = 1.0
        ranked.append(strength)
        axes.append(axis)
        means.append(center)
        scales.append(scale)
        projection_centers.append(projection_center)
        projection_scales.append(projection_scale)
    selected = _rank_disjoint_pairs(first, second, np.asarray(ranked), genes, 5)
    if not selected:
        return 0.0, 0.0, {"selected_pairs": []}
    test = np.vstack((test_core, test_edge))
    projections = []
    for index in selected:
        pair_x = (test[:, [first[index], second[index]]] - means[index]) / scales[index]
        projections.append((pair_x @ axes[index] - projection_centers[index]) / projection_scales[index])
    scores = np.column_stack(projections).mean(axis=1)
    details = [{"genes": [genes[first[index]], genes[second[index]]], "t_score": float(ranked[index])} for index in selected]
    return float(scores[0]), float(scores[1]), {"selected_pairs": details}


def _pair_features(values, pairs):
    import numpy as np

    if not pairs:
        return np.zeros((len(values), 0))
    first = np.asarray([pair[0] for pair in pairs], dtype=int)
    second = np.asarray([pair[1] for pair in pairs], dtype=int)
    directions = np.asarray([pair[2] for pair in pairs], dtype=float)
    return np.sign(values[:, first] - values[:, second]) * directions


def _pair_svm_score(core, edge, test_core, test_edge, pairs, *, class_weight=None):
    import numpy as np
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC

    if not pairs:
        return 0.0, 0.0
    train_x, train_y = _stack(core, edge)
    test_x = np.vstack((test_core, test_edge))
    model = make_pipeline(StandardScaler(), SVC(kernel="linear", C=1.0, class_weight=class_weight))
    scores = _score_estimator(model, _pair_features(train_x, pairs), train_y, _pair_features(test_x, pairs))
    return float(scores[0]), float(scores[1])


def _pairwise_svm_model(core, edge, test_core, test_edge, genes, pool, donor_ids):
    import numpy as np
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC

    rng = np.random.default_rng(0)
    correct: dict[tuple[int, int], int] = defaultdict(int)
    tested: dict[tuple[int, int], int] = defaultdict(int)
    patient_count = len(core)
    for _ in range(50):
        train_size = max(2, int(round(0.8 * patient_count)))
        sample = np.sort(rng.choice(patient_count, size=train_size, replace=False))
        validation = np.asarray([index for index in range(patient_count) if index not in set(sample)], dtype=int)
        if not len(validation):
            continue
        local_pool = _top_genes(core[sample], edge[sample], genes, RANK_GENE_CAP)
        first, second = _pair_arrays(local_pool, genes)
        train_x, train_y = _stack(core[sample], edge[sample])
        validation_x, validation_y = _stack(core[validation], edge[validation])
        for g1, g2 in zip(first, second):
            key = (int(g1), int(g2))
            model = make_pipeline(StandardScaler(), SVC(kernel="linear", C=1.0))
            model.fit(train_x[:, [g1, g2]], train_y)
            prediction = model.predict(validation_x[:, [g1, g2]])
            correct[key] += int(np.count_nonzero(prediction == validation_y))
            tested[key] += len(validation_y)
    selected = [key for key in tested if tested[key] and correct[key] / tested[key] >= 0.90]
    selected.sort(key=lambda pair: (-correct[pair] / tested[pair], genes[pair[0]], genes[pair[1]]))
    full_train, _ = _stack(core, edge)
    full_test = np.vstack((test_core, test_edge))
    scores = []
    train_y = np.r_[np.zeros(len(core), dtype=int), np.ones(len(edge), dtype=int)]
    for g1, g2 in selected:
        model = make_pipeline(StandardScaler(), SVC(kernel="linear", C=1.0))
        model.fit(full_train[:, [g1, g2]], train_y)
        scores.append(model.decision_function(full_test[:, [g1, g2]]))
    if not scores:
        return 0.0, 0.0, {"selected_pair_count": 0}
    mean_scores = np.column_stack(scores).mean(axis=1)
    details = [{"genes": [genes[g1], genes[g2]], "validation_accuracy": correct[(g1, g2)] / tested[(g1, g2)], "validation_cases": tested[(g1, g2)]} for g1, g2 in selected]
    return float(mean_scores[0]), float(mean_scores[1]), {"selected_pair_count": len(selected), "selected_pairs": details}


def _tspg_model(core, edge, test_core, test_edge, genes, donor_ids):
    import numpy as np

    rng = np.random.default_rng(0)
    side_one, side_two, appearances = Counter(), Counter(), Counter()
    patient_count = len(core)
    for _ in range(1000):
        sample_size = max(2, int(math.floor(0.9 * patient_count)))
        sample = np.sort(rng.choice(patient_count, size=sample_size, replace=False))
        local_core, local_edge = core[sample], edge[sample]
        pool = _top_genes(local_core, local_edge, genes, 40)
        first, second, difference = _tsp_statistics(local_core, local_edge, genes, pool)
        selected = _rank_disjoint_pairs(first, second, np.abs(difference), genes, 20)
        for index in selected:
            g1, g2 = int(first[index]), int(second[index])
            direction = 1 if difference[index] >= 0 else -1
            left, right = (g1, g2) if direction > 0 else (g2, g1)
            side_one[genes[left]] += 1
            side_two[genes[right]] += 1
            appearances[genes[left]] += 1
            appearances[genes[right]] += 1
    group_one, group_two = [], []
    for gene, count in appearances.items():
        left_fraction = side_one[gene] / count
        right_fraction = side_two[gene] / count
        if left_fraction >= 0.80:
            group_one.append(gene)
        elif right_fraction >= 0.80:
            group_two.append(gene)
    if not group_one or not group_two:
        return 0.0, 0.0, {"group_one": sorted(group_one), "group_two": sorted(group_two)}
    gene_index = {gene: index for index, gene in enumerate(genes)}
    test = np.vstack((test_core, test_edge))
    ranks = np.vstack([_percentile_scores(row, genes) for row in test])
    scores = ranks[:, [gene_index[gene] for gene in group_one]].mean(axis=1) - ranks[:, [gene_index[gene] for gene in group_two]].mean(axis=1)
    return float(scores[0]), float(scores[1]), {"group_one": sorted(group_one), "group_two": sorted(group_two)}


def _select_tsp_pairs(core, edge, genes, maximum=5):
    import numpy as np

    pool = _top_genes(core, edge, genes, RANK_GENE_CAP)
    first, second, difference = _tsp_statistics(core, edge, genes, pool)
    indices = _rank_disjoint_pairs(first, second, np.abs(difference), genes, maximum)
    return [(int(first[index]), int(second[index]), 1.0 if difference[index] >= 0 else -1.0) for index in indices]


def _inner_pair_svm_concordance(core, edge, genes, selector, feature_counts, *, class_weight=None):
    import numpy as np

    scores = {count: [] for count in feature_counts}
    for heldout in range(len(core)):
        train = np.asarray([index for index in range(len(core)) if index != heldout], dtype=int)
        pairs = selector(core[train], edge[train])
        if not pairs:
            for count in feature_counts:
                scores[count].append(0.5)
            continue
        x, y = _stack(core[train], edge[train])
        test_x = np.vstack((core[heldout], edge[heldout]))
        for count in feature_counts:
            selected = pairs[: min(count, len(pairs))]
            if not selected:
                scores[count].append(0.5)
                continue
            model = _pair_svc_model()
            model.fit(_pair_features(x, selected), y)
            prediction = model.decision_function(_pair_features(test_x, selected))
            scores[count].append(float(prediction[1] > prediction[0]) + 0.5 * float(prediction[1] == prediction[0]))
    return {count: float(np.mean(values)) for count, values in scores.items()}


def _pair_svc_model():
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    from sklearn.svm import SVC

    return make_pipeline(StandardScaler(), SVC(kernel="linear", C=1.0))


def _ktsp_svm_model(core, edge, test_core, test_edge, genes, donor_ids):
    import numpy as np

    candidate_counts = [1, 3, 5]
    inner = _inner_pair_svm_concordance(core, edge, genes, lambda c, e: _select_tsp_pairs(c, e, genes, 5), candidate_counts)
    chosen = min(candidate_counts, key=lambda count: (-inner[count], count))
    pairs = _select_tsp_pairs(core, edge, genes, 5)[:chosen]
    cs, es = _pair_svm_score(core, edge, test_core, test_edge, pairs)
    return cs, es, {"selected_pair_count": chosen, "inner_mean_concordance": inner, "selected_pairs": [(genes[a], genes[b], direction) for a, b, direction in pairs]}


def _reo_pairs(core, edge, genes):
    import numpy as np

    pool = sorted(_top_genes(core, edge, genes, RANK_GENE_CAP), key=lambda index: genes[index])
    first, second = _pair_arrays(pool, genes)
    x, y = _stack(core, edge)
    votes = _pair_votes(x, first, second)
    selected = []
    for index in range(votes.shape[1]):
        core_positive = np.mean(votes[y == 0, index] > 0)
        core_negative = np.mean(votes[y == 0, index] < 0)
        edge_positive = np.mean(votes[y == 1, index] > 0)
        edge_negative = np.mean(votes[y == 1, index] < 0)
        if core_positive >= 0.95 and edge_negative >= 0.95:
            selected.append((int(first[index]), int(second[index]), 1.0))
        elif core_negative >= 0.95 and edge_positive >= 0.95:
            selected.append((int(first[index]), int(second[index]), -1.0))
    return selected


def _mrmr_order(core, edge, genes, pairs):
    import numpy as np
    from sklearn.feature_selection import mutual_info_classif
    from sklearn.metrics import mutual_info_score

    if not pairs:
        return []
    x, y = _stack(core, edge)
    features = _pair_features(x, pairs)
    relevance = mutual_info_classif(features, y, discrete_features=True, random_state=0)
    remaining = set(range(len(pairs)))
    chosen = []
    while remaining:
        ranked = []
        for index in remaining:
            redundancy = np.mean([mutual_info_score(features[:, index], features[:, previous]) for previous in chosen]) if chosen else 0.0
            g1, g2, _ = pairs[index]
            ranked.append((-(float(relevance[index]) - redundancy), genes[g1], genes[g2], index))
        selected = min(ranked)[-1]
        chosen.append(selected)
        remaining.remove(selected)
    return [pairs[index] for index in chosen]


def _reo_svm_model(core, edge, test_core, test_edge, genes, donor_ids):
    import numpy as np

    pairs = _mrmr_order(core, edge, genes, _reo_pairs(core, edge, genes))
    if not pairs:
        return 0.0, 0.0, {"selected_pair_count": 0, "inner_mean_concordance": {}}
    candidate_counts = list(range(1, min(30, len(pairs)) + 1))
    inner = _inner_pair_svm_concordance(core, edge, genes, lambda c, e: _mrmr_order(c, e, genes, _reo_pairs(c, e, genes)), candidate_counts)
    chosen = min(candidate_counts, key=lambda count: (-inner[count], count))
    selected = pairs[:chosen]
    cs, es = _pair_svm_score(core, edge, test_core, test_edge, selected)
    return cs, es, {"selected_pair_count": chosen, "inner_mean_concordance": inner, "selected_pairs": [(genes[a], genes[b], d) for a, b, d in selected]}


def _tsp_chi2_pairs(core, edge, genes):
    import numpy as np
    from sklearn.feature_selection import chi2

    pool = sorted(_top_genes(core, edge, genes, RANK_GENE_CAP), key=lambda index: genes[index])
    first, second, difference = _tsp_statistics(core, edge, genes, pool)
    x, y = _stack(core, edge)
    votes = _pair_votes(x, first, second)
    chi_values, _ = chi2((votes + 1.0) / 2.0, y)
    tsp_order = np.argsort(-np.abs(difference), kind="stable")
    chi_order = np.argsort(-np.nan_to_num(chi_values, nan=-np.inf), kind="stable")
    tsp_ranks = np.empty(len(first), dtype=float)
    chi_ranks = np.empty(len(first), dtype=float)
    tsp_ranks[tsp_order] = np.arange(len(first))
    chi_ranks[chi_order] = np.arange(len(first))
    combined = tsp_ranks + chi_ranks
    pair_order = sorted(range(len(first)), key=lambda i: (combined[i], genes[first[i]], genes[second[i]]))
    direction = np.where(difference >= 0, 1.0, -1.0)
    selected, used = [], set()
    for index in pair_order:
        g1, g2 = int(first[index]), int(second[index])
        if genes[g1] in used or genes[g2] in used:
            continue
        selected.append((g1, g2, float(direction[index])))
        used.update((genes[g1], genes[g2]))
    return selected


def _tsp_ml_model(core, edge, test_core, test_edge, genes, donor_ids):
    import numpy as np

    candidate_counts = list(range(1, min(20, len(_tsp_chi2_pairs(core, edge, genes))) + 1))
    if not candidate_counts:
        return 0.0, 0.0, {"selected_pair_count": 0, "inner_mean_concordance": {}}
    inner = _inner_pair_svm_concordance(core, edge, genes, lambda c, e: _tsp_chi2_pairs(c, e, genes), candidate_counts)
    chosen = min(candidate_counts, key=lambda count: (-inner[count], count))
    pairs = _tsp_chi2_pairs(core, edge, genes)[:chosen]
    cs, es = _pair_svm_score(core, edge, test_core, test_edge, pairs)
    return cs, es, {"selected_pair_count": chosen, "inner_mean_concordance": inner, "selected_pairs": [(genes[a], genes[b], d) for a, b, d in pairs]}


def _md5_file(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _eligible_comparators(preflight: dict[str, Any]) -> tuple[list[str], dict[str, dict[str, Any]]]:
    coverage = preflight["signature_coverage"]
    eligible, decisions = [], {}
    for arm_id, spec in ARM_SPECS.items():
        reason = spec["eligibility"]
        if reason == "portable" or reason == "portable adapter":
            eligible.append(arm_id)
            decisions[arm_id] = {"eligible": True, "basis": reason}
        elif arm_id == "C17":
            fraction = coverage["C17"]["fraction"]
            passed = fraction >= 0.70
            decisions[arm_id] = {"eligible": passed, "coverage": fraction, "basis": "70% exact-symbol coverage gate"}
            if passed:
                eligible.append(arm_id)
        elif arm_id == "C18":
            fraction = coverage["C18"]["fraction"]
            passed = fraction >= 0.70
            decisions[arm_id] = {"eligible": passed, "coverage": fraction, "basis": "70% exact-symbol coverage gate"}
            if passed:
                eligible.append(arm_id)
        elif arm_id == "C26":
            le_fraction = coverage["C26_LE"]["fraction"]
            core_fraction = coverage["C26_CORE"]["fraction"]
            passed = le_fraction >= 0.70 and core_fraction >= 0.70
            decisions[arm_id] = {"eligible": passed, "coverage_le": le_fraction, "coverage_core": core_fraction, "basis": "70% exact-symbol coverage in both sets"}
            if passed:
                eligible.append(arm_id)
        else:
            decisions[arm_id] = {"eligible": False, "basis": reason}
    return eligible, decisions


def _score_locked_benchmark(matrix_path: Path, metadata_path: Path, preflight_path: Path, output_dir: Path) -> dict[str, Any]:
    import numpy as np
    from trackA_margin_evaluation import evaluate

    manifest = _verify_frozen_files()
    preflight = json.loads(preflight_path.read_text(encoding="utf-8"))
    if preflight.get("mode") != "metadata_only" or preflight.get("expression_values_opened") is not False or preflight.get("scores_calculated") is not False:
        raise RuntimeError("Input is not a clean metadata-only preflight")
    if preflight.get("score_gate_passed") is not True:
        raise RuntimeError("Metadata-only join and patient-count gates did not pass")
    if preflight.get("plan_sha256") != manifest["plan_sha256"] or preflight.get("analysis_code_sha256") != manifest["analysis_code_sha256"]:
        raise RuntimeError("Preflight plan/code hashes differ from committed manifest")
    if preflight.get("manifest_sha256") != sha256_file(MANIFEST_PATH):
        raise RuntimeError("Metadata-only preflight manifest hash differs from committed manifest")
    if preflight.get("matrix_file", {}).get("size_bytes") != EXPECTED_MATRIX_SIZE:
        raise RuntimeError("Downloaded matrix file size differs from source record")
    if preflight.get("metadata_file", {}).get("sha256") != sha256_file(metadata_path):
        raise RuntimeError("Metadata CSV changed after metadata-only preflight")
    if _md5_file(matrix_path) != EXPECTED_MATRIX_MD5:
        raise RuntimeError("Downloaded matrix MD5 differs from Zenodo source record")
    matrix_md5 = EXPECTED_MATRIX_MD5

    # Rescan only HDF5 structure and identifier arrays. No count dataset is touched here.
    fresh_preflight = _preflight(matrix_path, metadata_path)
    for key in ("matrix_layout", "gene_id_path", "gene_count", "matrix_cell_count", "eligible_patient_codes", "eligible_cell_count", "signature_coverage"):
        if fresh_preflight.get(key) != preflight.get(key):
            raise RuntimeError(f"Metadata-only preflight changed after lock: {key}")

    patients, sample_names, genes, counts = _read_counts_by_patient_region(matrix_path, metadata_path, preflight)
    if len(patients) < MINIMUM_PATIENTS:
        raise RuntimeError("Patient minimum missed after count aggregation; no score was calculated")
    expression = _log_cpm(counts)
    core, edge = expression[0::2], expression[1::2]
    if len(core) != len(patients) or len(edge) != len(patients):
        raise RuntimeError("Pseudobulk patient-region rows do not match patient list")

    comparators, decisions = _eligible_comparators(preflight)
    if len(comparators) > MAX_COMPARATORS:
        raise RuntimeError("Eligible comparator count exceeds frozen upper bound")
    arm_scores: dict[str, dict[str, tuple[float, float]]] = {"MP-RPV": {}}
    for arm_id in comparators:
        arm_scores[arm_id] = {}
    details: dict[str, dict[str, Any]] = defaultdict(dict)
    method_failures: dict[str, str] = {}

    for heldout, patient in enumerate(patients):
        train = np.asarray([index for index in range(len(patients)) if index != heldout], dtype=int)
        test_core, test_edge = core[heldout], edge[heldout]
        candidate_core, candidate_edge, candidate_info = _candidate_scores(core[train], edge[train], test_core, test_edge, genes)
        arm_scores["MP-RPV"][patient] = (candidate_core, candidate_edge)
        details["MP-RPV"][patient] = candidate_info
        for arm_id in comparators:
            if arm_id in method_failures:
                continue
            try:
                score_core, score_edge, info = _linear_model_score(arm_id, core[train], edge[train], test_core, test_edge, genes, [patients[index] for index in train])
                arm_scores[arm_id][patient] = (score_core, score_edge)
                details[arm_id][patient] = info
            except Exception as error:
                method_failures[arm_id] = f"{type(error).__name__}: {error}"

    missing = {arm: sorted(set(patients) - set(scores)) for arm, scores in arm_scores.items() if set(scores) != set(patients)}
    benchmark_complete = not method_failures and not missing
    evaluation = None
    if benchmark_complete:
        evaluation = evaluate(arm_scores, "MP-RPV", comparators, minimum_donors=MINIMUM_PATIENTS, alpha=ALPHA)

    concordance = {}
    score_rows = []
    for arm, by_patient in arm_scores.items():
        per_patient = {}
        for patient, (core_score, edge_score) in by_patient.items():
            value = 1.0 if edge_score > core_score else 0.0 if edge_score < core_score else 0.5
            per_patient[patient] = value
            score_rows.append({"patient": patient, "arm": arm, "P_CT_score": core_score, "P_IT_score": edge_score, "concordance": value})
        concordance[arm] = {"per_patient": per_patient, "mean": float(np.mean(list(per_patient.values()))) if per_patient else None}

    result = {
        "analysis_label": "exploratory",
        "question": manifest["target"],
        "plan_sha256": manifest["plan_sha256"],
        "analysis_code_sha256": manifest["analysis_code_sha256"],
        "matrix_md5": matrix_md5,
        "metadata_sha256": sha256_file(metadata_path),
        "patient_count": len(patients),
        "patients": patients,
        "gene_count": len(genes),
        "comparator_eligibility": decisions,
        "eligible_comparator_count": len(comparators),
        "eligible_comparators": comparators,
        "benchmark_complete": benchmark_complete,
        "method_failures": method_failures,
        "missing_patient_predictions": missing,
        "concordance": concordance,
        "evaluation": evaluation,
        "passed": bool(evaluation and evaluation["passed"]),
        "caveat": "Exploratory within-study held-out patient benchmark. Artzi was selected after prior Track A negatives. Whole-tissue ecosystem endpoint, not tumor-cell-intrinsic biology or drug response.",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_json(output_dir / "result.json", result)
    _write_json(output_dir / "model_details.json", details)
    if score_rows:
        with (output_dir / "patient_scores.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["patient", "arm", "P_CT_score", "P_IT_score", "concordance"], delimiter="\t")
            writer.writeheader()
            writer.writerows(score_rows)
    return result


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--write-manifest", action="store_true")
    mode.add_argument("--preflight", action="store_true")
    mode.add_argument("--score", action="store_true")
    parser.add_argument("--matrix", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--preflight-json", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.write_manifest:
        _write_json(MANIFEST_PATH, build_manifest())
        return
    if args.preflight:
        if not args.matrix or not args.metadata or not args.output:
            parser.error("--preflight requires --matrix, --metadata, and --output")
        try:
            matrix_path = project_path(args.matrix)
            metadata_path = project_path(args.metadata)
            result = _preflight(matrix_path, metadata_path)
        except Exception as error:
            result = {
                "analysis_label": "exploratory",
                "mode": "metadata_only",
                "expression_values_opened": False,
                "scores_calculated": False,
                "plan_sha256": sha256_file(PLAN_PATH),
                "analysis_code_sha256": combined_code_hash(),
                "manifest_sha256": sha256_file(MANIFEST_PATH) if MANIFEST_PATH.exists() else None,
                "status": "metadata_preflight_failed",
                "matrix_file": {"name": args.matrix.name, "size_bytes": matrix_path.stat().st_size if matrix_path.exists() else None},
                "metadata_file": {"name": args.metadata.name, "size_bytes": metadata_path.stat().st_size if metadata_path.exists() else None},
                "score_gate_passed": False,
                "failure": f"{type(error).__name__}: {error}",
            }
        _write_json(project_path(args.output), result)
        print(json.dumps({key: result[key] for key in ("mode", "expression_values_opened", "eligible_patient_count", "score_gate_passed", "signature_coverage") if key in result}, indent=2))
        return
    if args.score:
        if not args.matrix or not args.metadata or not args.preflight_json or not args.output:
            parser.error("--score requires --matrix, --metadata, --preflight-json, and --output")
        matrix_path = project_path(args.matrix)
        metadata_path = project_path(args.metadata)
        preflight_path = project_path(args.preflight_json)
        output_path = project_path(args.output)
        if not matrix_path.exists() or not metadata_path.exists() or not preflight_path.exists():
            raise FileNotFoundError("Matrix, metadata, and committed metadata-only preflight must exist")
        result = _score_locked_benchmark(matrix_path, metadata_path, preflight_path, output_path)
        print(json.dumps({"benchmark_complete": result["benchmark_complete"], "patient_count": result["patient_count"], "eligible_comparator_count": result["eligible_comparator_count"], "passed": result["passed"]}, indent=2))
        return


if __name__ == "__main__":
    main()
