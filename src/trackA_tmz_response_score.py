from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
from importlib.metadata import PackageNotFoundError, version
import json
import platform
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from scipy.stats import rankdata


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PANEL_PATH = PROJECT_ROOT / "paper" / "trackA_tmz_response_121_gene_panel_2026-11-07.tsv"
PROTOCOL_PATH = PROJECT_ROOT / "paper" / "trackA_tmz_response_score_protocol_2026-11-07.md"
FREEZE_PATH = PROJECT_ROOT / "paper" / "trackA_tmz_response_score_SHA256SUMS"
A5_RUNNER_PATH = PROJECT_ROOT / "src" / "189_trackA_yu2020_lodo.py"
OUTPUT_PATH = PROJECT_ROOT / "output" / "trackA_tmz_response_proxy" / "result.json"
PRIMARY_DONORS = ("GS1", "GS2", "GS3", "GS5", "GS6", "GS13")
N_CONTROL_PANELS = 1_000
CONTROL_SEED = 20261107
MINIMUM_CELLS_PER_REGION = 20
EXPECTED_PANEL_SIZE = 121
EXPECTED_POSITIVE = 48
EXPECTED_NEGATIVE = 73
SOURCE_SUPPLEMENT_SHA256 = "a21740f50dfc096550b2576f97bff90b70eb6ade62eeba8d25acd3f2b092c7ce"
SOURCE_SUPPLEMENT_MD5 = "d4c65e85ad84c441856590f355a591e3"


class ProxyStop(RuntimeError):
    """Raised when a frozen proxy-analysis gate fails."""


def _validated_count_matrix(counts: np.ndarray) -> np.ndarray:
    matrix = np.asarray(counts)
    if matrix.ndim != 2 or not np.issubdtype(matrix.dtype, np.number) or np.issubdtype(matrix.dtype, np.bool_):
        raise ProxyStop("Counts must be a two-dimensional numeric matrix")
    if not np.isfinite(matrix).all() or np.any(matrix < 0):
        raise ProxyStop("Counts must be finite and nonnegative")
    if np.issubdtype(matrix.dtype, np.floating) and not np.equal(matrix, np.floor(matrix)).all():
        raise ProxyStop("Counts must be integers")
    return matrix


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_panel(path: Path = PANEL_PATH) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {
        "gene_symbol",
        "rho_viability_100uM",
        "p_value",
        "p_adjusted",
        "direction",
    }
    if not rows or not required.issubset(rows[0]):
        raise ProxyStop("The frozen gene panel has missing columns or no rows")
    if len(rows) != EXPECTED_PANEL_SIZE:
        raise ProxyStop(f"Expected 121 panel genes; found {len(rows)}")
    symbols = [row["gene_symbol"].strip() for row in rows]
    if any(not symbol for symbol in symbols) or len(symbols) != len(set(symbols)):
        raise ProxyStop("The frozen gene panel has empty or duplicate symbols")
    positive = negative = 0
    for row in rows:
        rho = float(row["rho_viability_100uM"])
        direction = "resistance_associated" if rho > 0 else "sensitivity_associated" if rho < 0 else "invalid"
        if direction == "invalid" or row["direction"] != direction:
            raise ProxyStop(f"Gene panel direction does not match coefficient for {row['gene_symbol']}")
        positive += rho > 0
        negative += rho < 0
    if (positive, negative) != (EXPECTED_POSITIVE, EXPECTED_NEGATIVE):
        raise ProxyStop(f"Expected 48 positive and 73 negative genes; found {positive} and {negative}")
    return rows


def verify_freeze(root: Path = PROJECT_ROOT, freeze_path: Path = FREEZE_PATH) -> dict[str, str]:
    expected_files = {
        "paper/trackA_tmz_response_121_gene_panel_2026-11-07.tsv",
        "paper/trackA_tmz_response_score_protocol_2026-11-07.md",
        "src/trackA_tmz_response_score.py",
        "tests/test_trackA_tmz_response_score.py",
        "src/189_trackA_yu2020_lodo.py",
        "paper/trackA_yu2020/PROTOCOL.md",
        "paper/trackA_yu2020/SHA256SUMS",
    }
    if not freeze_path.is_file():
        raise ProxyStop("Proxy freeze manifest is missing")
    entries: dict[str, str] = {}
    for line in freeze_path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        parts = text.split(maxsplit=1)
        if len(parts) != 2:
            raise ProxyStop("Malformed proxy freeze manifest")
        digest, relative_path = parts
        relative_path = relative_path.strip().lstrip("*")
        if relative_path in entries:
            raise ProxyStop(f"Duplicate freeze entry: {relative_path}")
        entries[relative_path] = digest.lower()
    if set(entries) != expected_files:
        missing = sorted(expected_files - set(entries))
        extra = sorted(set(entries) - expected_files)
        raise ProxyStop(f"Proxy freeze manifest mismatch; missing={missing}; extra={extra}")
    for relative_path, expected in entries.items():
        path = root / relative_path
        if not path.is_file():
            raise ProxyStop(f"Frozen file is missing: {relative_path}")
        actual = sha256_file(path)
        if actual != expected:
            raise ProxyStop(f"Frozen file hash mismatch: {relative_path}")
    return entries


def _gene_index(gene_ids: Iterable[str], panel_symbols: Iterable[str]) -> dict[str, int]:
    ids = [str(gene) for gene in gene_ids]
    if len(ids) != len(set(ids)):
        raise ProxyStop("Count matrix contains duplicate gene identifiers")
    index = {gene: i for i, gene in enumerate(ids)}
    missing = sorted(set(panel_symbols) - set(index))
    if missing:
        raise ProxyStop(f"Exact panel mapping failed; {len(missing)} source symbols are missing")
    return index


def rank_percentiles(log_expression: np.ndarray) -> np.ndarray:
    values = np.asarray(log_expression, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
        raise ProxyStop("Rank input must be a finite one-dimensional expression vector")
    return (rankdata(values, method="average") - 1.0) / (values.size - 1.0)


def score_rank_percentiles(
    rank_percentile: np.ndarray,
    gene_ids: Iterable[str],
    positive_genes: Iterable[str],
    negative_genes: Iterable[str],
) -> float:
    ranks = np.asarray(rank_percentile, dtype=np.float64)
    ids = [str(gene) for gene in gene_ids]
    if ranks.ndim != 1 or len(ids) != ranks.size or not np.isfinite(ranks).all():
        raise ProxyStop("Gene IDs and finite rank percentiles must have matching lengths")
    index = _gene_index(ids, (*positive_genes, *negative_genes))
    positive = [index[gene] for gene in positive_genes]
    negative = [index[gene] for gene in negative_genes]
    if not positive or not negative or set(positive) & set(negative):
        raise ProxyStop("Positive and negative gene groups must be nonempty and disjoint")
    return float(np.mean(ranks[positive]) - np.mean(ranks[negative]))


def exact_two_sided_sign_flip(differences: Iterable[float]) -> dict[str, Any]:
    values = np.asarray(tuple(differences), dtype=np.float64)
    if values.ndim != 1 or values.size != 6 or not np.isfinite(values).all():
        raise ProxyStop("Sign-flip test requires six finite donor-level differences")
    observed = float(values.mean())
    signs = np.asarray(
        [[1.0 if (mask >> bit) & 1 else -1.0 for bit in range(values.size)] for mask in range(1 << values.size)],
        dtype=np.float64,
    )
    null_statistics = (signs * values[None, :]).mean(axis=1)
    exceedances = int(np.count_nonzero(np.abs(null_statistics) >= abs(observed) - 1e-15))
    return {
        "method": "exact_two_sided_sign_flip",
        "analysis_unit": "donor",
        "n_donors": int(values.size),
        "sign_assignments": int(1 << values.size),
        "mean_difference": observed,
        "p_value": float(exceedances / (1 << values.size)),
    }


def _per_cell_expression_features(counts: np.ndarray, chunk_size: int = 512) -> tuple[np.ndarray, np.ndarray]:
    matrix = _validated_count_matrix(counts)
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ProxyStop("Matched-control features need a nonempty cell-by-gene count matrix")
    totals = matrix.sum(axis=1, dtype=np.float64)
    if np.any(totals <= 0):
        raise ProxyStop("A cell has zero total counts")
    mean_log_cpm = np.empty(matrix.shape[1], dtype=np.float64)
    detection = np.empty(matrix.shape[1], dtype=np.float64)
    for start in range(0, matrix.shape[1], chunk_size):
        stop = min(start + chunk_size, matrix.shape[1])
        block = np.asarray(matrix[:, start:stop], dtype=np.float64)
        log_cpm = np.log2(block / totals[:, None] * 1_000_000.0 + 1.0)
        mean_log_cpm[start:stop] = log_cpm.mean(axis=0)
        detection[start:stop] = (block > 0).mean(axis=0)
    return mean_log_cpm, detection


def generate_matched_control_panels(
    counts: np.ndarray,
    gene_ids: Iterable[str],
    panel: list[dict[str, str]],
    *,
    n_panels: int = N_CONTROL_PANELS,
    seed: int = CONTROL_SEED,
    nearest_pool: int = 20,
) -> list[tuple[tuple[str, ...], tuple[str, ...]]]:
    ids = [str(gene) for gene in gene_ids]
    if len(ids) != len(set(ids)):
        raise ProxyStop("Count matrix contains duplicate gene identifiers")
    symbol_index = _gene_index(ids, (row["gene_symbol"] for row in panel))
    source_ids = [row["gene_symbol"] for row in panel]
    source_indices = [symbol_index[gene] for gene in source_ids]
    source_set = set(source_indices)
    background = [i for i in range(len(ids)) if i not in source_set]
    if len(background) < len(panel):
        raise ProxyStop("Too few background genes for unique matched control panels")
    mean_log_cpm, detection = _per_cell_expression_features(counts)
    features = np.column_stack((mean_log_cpm, detection))
    center = features.mean(axis=0)
    scale = features.std(axis=0)
    scale[scale == 0] = 1.0
    standardized = (features - center) / scale
    background_features = standardized[background]
    nearest_by_source: list[list[int]] = []
    for source_index in source_indices:
        distance = np.sqrt(np.square(background_features - standardized[source_index]).sum(axis=1))
        order = sorted(range(len(background)), key=lambda k: (float(distance[k]), ids[background[k]]))
        nearest_by_source.append([background[k] for k in order])

    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    out: list[tuple[tuple[str, ...], tuple[str, ...]]] = []
    positive_positions = [i for i, row in enumerate(panel) if float(row["rho_viability_100uM"]) > 0]
    negative_positions = [i for i, row in enumerate(panel) if float(row["rho_viability_100uM"]) < 0]
    for _ in range(n_panels):
        used: set[int] = set()
        matched: list[str | None] = [None] * len(panel)
        for source_position in rng.permutation(len(panel)):
            available: list[int] = []
            for candidate in nearest_by_source[int(source_position)]:
                if candidate not in used:
                    available.append(candidate)
                    if len(available) == nearest_pool:
                        break
            if not available:
                raise ProxyStop("Could not form a unique expression-matched control panel")
            chosen = available[int(rng.integers(len(available)))]
            used.add(chosen)
            matched[int(source_position)] = ids[chosen]
        if any(gene is None for gene in matched) or len(set(matched)) != len(panel):
            raise ProxyStop("A matched control panel is incomplete or contains duplicates")
        control_genes = [str(gene) for gene in matched]
        out.append(
            (
                tuple(control_genes[i] for i in positive_positions),
                tuple(control_genes[i] for i in negative_positions),
            )
        )
    return out


def compute_proxy_scores(
    counts: np.ndarray,
    donors: Iterable[str],
    labels: Iterable[int],
    gene_ids: Iterable[str],
    panel: list[dict[str, str]],
    *,
    control_counts: np.ndarray | None = None,
    n_control_panels: int = N_CONTROL_PANELS,
    control_seed: int = CONTROL_SEED,
) -> dict[str, Any]:
    matrix = _validated_count_matrix(counts)
    donor_values = np.asarray(tuple(donors), dtype=str)
    label_values = np.asarray(tuple(labels), dtype=np.int64)
    ids = [str(gene) for gene in gene_ids]
    if matrix.ndim != 2 or matrix.shape != (len(donor_values), len(ids)) or len(label_values) != len(donor_values):
        raise ProxyStop("Counts, donors, labels, and genes have inconsistent shapes")
    if len(ids) != len(set(ids)):
        raise ProxyStop("Count matrix contains duplicate gene identifiers")
    if not set(np.unique(label_values)).issubset({0, 1}):
        raise ProxyStop("Labels must be 0=Periphery or 1=Core")
    positive_genes = tuple(row["gene_symbol"] for row in panel if float(row["rho_viability_100uM"]) > 0)
    negative_genes = tuple(row["gene_symbol"] for row in panel if float(row["rho_viability_100uM"]) < 0)
    _gene_index(ids, (*positive_genes, *negative_genes))
    observed_donors = set(donor_values.tolist())
    if observed_donors != set(PRIMARY_DONORS):
        raise ProxyStop(f"Expected A5 primary donor set {PRIMARY_DONORS}; found {sorted(observed_donors)}")
    for donor in PRIMARY_DONORS:
        for label in (1, 0):
            n = int(np.count_nonzero((donor_values == donor) & (label_values == label)))
            if n < MINIMUM_CELLS_PER_REGION:
                raise ProxyStop(f"Donor {donor} has only {n} QC-passed cells for label {label}")

    sample_log_cpm: dict[tuple[str, int], np.ndarray] = {}
    sample_cell_counts: dict[str, dict[str, int]] = {}
    for donor in PRIMARY_DONORS:
        sample_cell_counts[donor] = {}
        for label, region in ((1, "core"), (0, "periphery")):
            keep = (donor_values == donor) & (label_values == label)
            sample_cell_counts[donor][region] = int(keep.sum())
            pooled = matrix[keep].sum(axis=0, dtype=np.float64)
            total = float(pooled.sum())
            if total <= 0:
                raise ProxyStop(f"Donor {donor} {region} pseudobulk has no counts")
            sample_log_cpm[(donor, label)] = np.log2(pooled / total * 1_000_000.0 + 1.0)

    main_scores: dict[str, dict[str, float]] = {}
    rank_vectors: dict[tuple[str, int], np.ndarray] = {}
    for donor in PRIMARY_DONORS:
        for label in (1, 0):
            rank_vectors[(donor, label)] = rank_percentiles(sample_log_cpm[(donor, label)])
        core = score_rank_percentiles(rank_vectors[(donor, 1)], ids, positive_genes, negative_genes)
        periphery = score_rank_percentiles(rank_vectors[(donor, 0)], ids, positive_genes, negative_genes)
        main_scores[donor] = {
            "core_score": core,
            "periphery_score": periphery,
            "core_minus_periphery": core - periphery,
        }
    deltas = [main_scores[donor]["core_minus_periphery"] for donor in PRIMARY_DONORS]

    if control_counts is None:
        control_counts = matrix
    controls = generate_matched_control_panels(
        control_counts,
        ids,
        panel,
        n_panels=n_control_panels,
        seed=control_seed,
    )
    index = {gene: i for i, gene in enumerate(ids)}
    control_mean_differences = np.empty(len(controls), dtype=np.float64)
    for panel_index, (control_positive, control_negative) in enumerate(controls):
        donor_differences = []
        for donor in PRIMARY_DONORS:
            core_ranks = rank_vectors[(donor, 1)]
            peri_ranks = rank_vectors[(donor, 0)]
            control_pos = [index[gene] for gene in control_positive]
            control_neg = [index[gene] for gene in control_negative]
            core_score = float(np.mean(core_ranks[control_pos]) - np.mean(core_ranks[control_neg]))
            periphery_score = float(np.mean(peri_ranks[control_pos]) - np.mean(peri_ranks[control_neg]))
            donor_differences.append(core_score - periphery_score)
        control_mean_differences[panel_index] = float(np.mean(donor_differences))

    observed_mean = float(np.mean(deltas))
    controls_summary = {
        "n_panels": int(len(controls)),
        "seed": int(control_seed),
        "selection": "one-to-one matching on mean per-cell log2(CPM+1) and detection fraction; randomly selected among 20 nearest available genes",
        "observed_mean_delta_percentile": float(np.mean(control_mean_differences <= observed_mean)),
        "control_mean_delta_median": float(np.median(control_mean_differences)),
        "control_mean_delta_q025": float(np.quantile(control_mean_differences, 0.025)),
        "control_mean_delta_q975": float(np.quantile(control_mean_differences, 0.975)),
        "use": "descriptive matched reference; not a significance test",
    }
    signs = {"positive": int(np.sum(np.asarray(deltas) > 0)), "negative": int(np.sum(np.asarray(deltas) < 0)), "zero": int(np.sum(np.asarray(deltas) == 0))}
    return {
        "record_type": "exploratory_trackA_tmz_response_direction_score",
        "status": "complete_exploratory",
        "cohort": {"donors": list(PRIMARY_DONORS), "n_donors": len(PRIMARY_DONORS), "cell_counts_by_donor_region": sample_cell_counts},
        "panel": {"n_genes": len(panel), "positive_resistance_associated": len(positive_genes), "negative_sensitivity_associated": len(negative_genes), "exact_symbol_mapping": "passed"},
        "score": {"definition": "mean transcriptome-wide expression percentile of positive-coefficient genes minus mean percentile of negative-coefficient genes", "range": [-1.0, 1.0], "direction": "higher means expression ranks align more with higher in-vitro viability at 100 uM TMZ; not measured drug response", "paired_difference": "Core minus Periphery", "donor_results": main_scores, "mean": observed_mean, "median": float(np.median(deltas)), "minimum": float(np.min(deltas)), "maximum": float(np.max(deltas)), "direction_counts": signs},
        "exact_sign_flip": exact_two_sided_sign_flip(deltas),
        "matched_controls": controls_summary,
        "interpretation": "mixed-cell regional expression proxy only; no donor is labeled sensitive or resistant; not a validated 121-gene predictor; does not calibrate Track C or link Track A to Duke",
    }


def _load_frozen_a5_dataset() -> tuple[Any, tuple[str, ...], dict[str, Any], dict[str, str]]:
    spec = importlib.util.spec_from_file_location("trackA_a5_frozen_loader", A5_RUNNER_PATH)
    if spec is None or spec.loader is None:
        raise ProxyStop("Could not load the frozen A5 data loader")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    dataset, strict_donors, _sensitivity_donors, flow, input_hashes = module.load_dataset()
    if tuple(strict_donors) != PRIMARY_DONORS:
        raise ProxyStop(f"Frozen A5 primary donors changed: {strict_donors}")
    return dataset, tuple(strict_donors), flow, input_hashes


def run_proxy_analysis(output_path: Path = OUTPUT_PATH) -> dict[str, Any]:
    if output_path.exists():
        raise ProxyStop(f"Refusing to overwrite existing result: {output_path}")
    frozen_hashes = verify_freeze()
    panel = read_panel()
    dataset, donors, flow, input_hashes = _load_frozen_a5_dataset()
    result = compute_proxy_scores(dataset.counts, dataset.donors, dataset.labels, dataset.gene_ids, panel)
    package_versions: dict[str, str] = {}
    for package in ("numpy", "scipy", "pandas", "scikit-learn", "torch"):
        try:
            package_versions[package] = version(package)
        except PackageNotFoundError:
            package_versions[package] = "not-installed"
    result["provenance"] = {
        "git_commit": _git_commit(),
        "freeze_manifest_sha256": sha256_file(FREEZE_PATH),
        "frozen_file_sha256": frozen_hashes,
        "panel_sha256": sha256_file(PANEL_PATH),
        "protocol_sha256": sha256_file(PROTOCOL_PATH),
        "a5_input_hashes": input_hashes,
        "a5_cell_flow": {key: value for key, value in flow.items() if key != "input_hashes"},
        "source_supplement_pmcid": "PMC10575865",
        "source_supplement_sha256": SOURCE_SUPPLEMENT_SHA256,
        "source_supplement_md5": SOURCE_SUPPLEMENT_MD5,
        "python": sys.version,
        "platform": platform.platform(),
        "package_versions": package_versions,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if temporary_path.exists():
        raise ProxyStop(f"Temporary result path already exists: {temporary_path}")
    temporary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary_path.replace(output_path)
    return result


def _git_commit() -> str | None:
    import subprocess

    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=PROJECT_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Exploratory Track A 121-gene TMZ response-direction score")
    parser.add_argument("--run", action="store_true", required=True, help="Run once under the frozen protocol")
    args = parser.parse_args()
    if not args.run:
        return 2
    try:
        result = run_proxy_analysis()
    except ProxyStop as exc:
        print(f"STOP: {exc}", file=sys.stderr, flush=True)
        return 2
    print(json.dumps({"status": result["status"], "n_donors": result["cohort"]["n_donors"], "result_path": str(OUTPUT_PATH)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
