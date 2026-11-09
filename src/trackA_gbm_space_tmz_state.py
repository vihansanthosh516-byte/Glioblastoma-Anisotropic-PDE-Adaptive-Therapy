from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any, Iterable, Sequence

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata, t


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = Path.home() / "Downloads" / "GBM_space_snRNA.h5ad"
DEFAULT_OUTPUT = PROJECT_ROOT / "output" / "trackA_gbm_space_tmz_state" / "result.json"
PROTOCOL_PATH = PROJECT_ROOT / "paper" / "trackA_gbm_space_tmz_state_protocol_2026-11-08.md"
PANEL_PATH = PROJECT_ROOT / "paper" / "trackA_tmz_response_121_gene_panel_2026-11-07.tsv"
SEAL_PATH = PROJECT_ROOT / "paper" / "trackA_gbm_space_tmz_state_SHA256SUMS"
INPUT_FINGERPRINT_PATH = PROJECT_ROOT / "paper" / "trackA_gbm_space_tmz_state_input_fingerprint_2026-11-08.json"

EXPECTED_SHAPE = (1_025_329, 36_601)
EXPECTED_DONORS = (
    "AT3", "AT4", "AT5", "AT6", "AT7", "AT9", "AT10", "AT11", "AT12", "AT13", "AT14", "AT15"
)
DEVELOPMENTAL_STATES = frozenset(
    {"AC-progenitor-like", "OPC-like", "NPC-neuronal-like", "OPC-NPC-like", "OPC-neuronal-like"}
)
INJURY_HYPOXIA_STATES = frozenset({"AC-gliosis-like", "Gliosis-like", "Hypoxic"})
EXPECTED_A5_OMITTED = frozenset({"AF131215.2", "PPIAP22"})
EXPECTED_GBM_SPACE_OMITTED = frozenset({"CROCCP3", "HERC2P2", "HKR1", "LINCR-0001"})
EXPECTED_PANEL_SIZE = 121
EXPECTED_SHARED_SIZE = 115
EXPECTED_POSITIVE = 47
EXPECTED_NEGATIVE = 68
MIN_CELLS_PER_DONOR_GROUP = 200
MATCH_SAMPLE_SIZE = 25_000
CONTROL_PANEL_COUNT = 1_000
CONTROL_SEED = 20261108
CONTROL_NEAREST_POOL = 20
BLOCK_ROWS = 512


class AnalysisStop(RuntimeError):
    """Raised when a frozen input, cohort, or scoring gate fails."""


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def verify_input_fingerprint(
    input_path: Path,
    fingerprint_path: Path = INPUT_FINGERPRINT_PATH,
) -> str:
    input_path = Path(input_path)
    fingerprint = json.loads(Path(fingerprint_path).read_text(encoding="utf-8"))
    if input_path.name != fingerprint.get("filename"):
        raise AnalysisStop("Input H5AD filename differs from the frozen input fingerprint")
    if input_path.stat().st_size != fingerprint.get("bytes"):
        raise AnalysisStop("Input H5AD byte count differs from the frozen input fingerprint")
    digest = sha256_file(input_path)
    if digest != fingerprint.get("sha256"):
        raise AnalysisStop("Input H5AD SHA-256 differs from the frozen pre-score fingerprint")
    return digest


def git_commit(root: Path = PROJECT_ROOT) -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True, capture_output=True, text=True
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return result.stdout.strip()


def verify_seal(root: Path = PROJECT_ROOT, seal_path: Path = SEAL_PATH) -> dict[str, str]:
    root = Path(root)
    seal_path = Path(seal_path)
    expected_files = {
        "paper/trackA_gbm_space_tmz_state_protocol_2026-11-08.md",
        "paper/trackA_gbm_space_tmz_state_input_fingerprint_2026-11-08.json",
        "paper/trackA_tmz_response_121_gene_panel_2026-11-07.tsv",
        "paper/trackA_tmz_response_score_protocol_2026-11-07.md",
        "paper/trackA_tmz_response_119of121_amendment_2026-11-08.md",
        "src/trackA_tmz_response_score.py",
        "src/trackA_gbm_space_tmz_state.py",
        "tests/test_trackA_gbm_space_tmz_state.py",
    }
    if not seal_path.is_file():
        raise AnalysisStop("Pre-score SHA256 seal is missing")
    entries: dict[str, str] = {}
    for line in seal_path.read_text(encoding="utf-8").splitlines():
        text = line.strip()
        if not text or text.startswith("#"):
            continue
        parts = text.split(maxsplit=1)
        if len(parts) != 2:
            raise AnalysisStop("Malformed pre-score SHA256 seal")
        digest, relative_path = parts
        relative_path = relative_path.strip().lstrip("*")
        if relative_path in entries:
            raise AnalysisStop(f"Duplicate file in pre-score seal: {relative_path}")
        entries[relative_path] = digest.lower()
    if set(entries) != expected_files:
        raise AnalysisStop(
            "Pre-score seal file list mismatch; "
            f"missing={sorted(expected_files - set(entries))}; extra={sorted(set(entries) - expected_files)}"
        )
    for relative_path, expected in entries.items():
        path = root / relative_path
        if not path.is_file() or sha256_file(path) != expected:
            raise AnalysisStop(f"Pre-score seal mismatch: {relative_path}")
    return entries


def read_source_panel(path: Path = PANEL_PATH) -> list[dict[str, Any]]:
    frame = pd.read_csv(path, sep="\t")
    required = {"gene_symbol", "rho_viability_100uM"}
    if not required.issubset(frame.columns):
        raise AnalysisStop("Frozen source panel has missing required columns")
    rows = frame.to_dict(orient="records")
    symbols = [str(row["gene_symbol"]).strip() for row in rows]
    if len(rows) != EXPECTED_PANEL_SIZE or len(symbols) != len(set(symbols)):
        raise AnalysisStop("Frozen source panel must contain 121 unique symbols")
    for row, symbol in zip(rows, symbols, strict=True):
        row["gene_symbol"] = symbol
        rho = float(row["rho_viability_100uM"])
        if not np.isfinite(rho) or rho == 0:
            raise AnalysisStop("Source panel directions must be finite and nonzero")
    return rows


def map_shared_panel(
    panel: Sequence[dict[str, Any]],
    gene_symbols: Iterable[Any],
    gene_ids: Iterable[Any],
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    source_symbols = [str(row["gene_symbol"]) for row in panel]
    if len(source_symbols) != EXPECTED_PANEL_SIZE or len(set(source_symbols)) != EXPECTED_PANEL_SIZE:
        raise AnalysisStop("The unchanged 121-gene source panel is required")
    ids = [str(value) for value in gene_ids]
    if len(ids) != len(set(ids)):
        raise AnalysisStop("GBM-Space raw.var gene identifiers are not unique")
    symbol_values = ["" if pd.isna(value) else str(value).strip() for value in gene_symbols]
    if len(symbol_values) != len(ids):
        raise AnalysisStop("GBM-Space raw.var symbols and gene IDs have different lengths")
    symbol_positions: dict[str, list[int]] = {}
    for index, symbol in enumerate(symbol_values):
        if symbol:
            symbol_positions.setdefault(symbol, []).append(index)
    a5_panel = [row for row in panel if row["gene_symbol"] not in EXPECTED_A5_OMITTED]
    if len(a5_panel) != 119:
        raise AnalysisStop("A5 119-of-121 panel definition changed")
    missing_from_a5 = set(source_symbols) - {row["gene_symbol"] for row in a5_panel}
    if missing_from_a5 != EXPECTED_A5_OMITTED:
        raise AnalysisStop("The A5 source-panel omissions differ from the frozen amendment")
    shared_missing = {row["gene_symbol"] for row in a5_panel if row["gene_symbol"] not in symbol_positions}
    if shared_missing != EXPECTED_GBM_SPACE_OMITTED:
        raise AnalysisStop("GBM-Space panel omissions differ from the metadata-only mapping")
    mapped = [row for row in a5_panel if row["gene_symbol"] in symbol_positions]
    if len(mapped) != EXPECTED_SHARED_SIZE:
        raise AnalysisStop("The cross-dataset exact panel must contain 115 genes")
    symbol_to_index: dict[str, int] = {}
    for row in mapped:
        symbol = row["gene_symbol"]
        positions = symbol_positions[symbol]
        if len(positions) != 1:
            raise AnalysisStop(f"GBM-Space target symbol is not one-to-one: {symbol}")
        symbol_to_index[symbol] = positions[0]
    positive = sum(float(row["rho_viability_100uM"]) > 0 for row in mapped)
    negative = sum(float(row["rho_viability_100uM"]) < 0 for row in mapped)
    if (positive, negative) != (EXPECTED_POSITIVE, EXPECTED_NEGATIVE):
        raise AnalysisStop("Shared panel direction counts differ from the frozen 47/68 specification")
    return mapped, symbol_to_index


def classify_metadata(
    obs: pd.DataFrame,
    *,
    expected_donors: Sequence[str] = EXPECTED_DONORS,
    minimum_cells: int = MIN_CELLS_PER_DONOR_GROUP,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, dict[str, int]]]:
    required = {"donor_id", "cell_status", "annotation_coarse", "total_counts"}
    if not required.issubset(obs.columns):
        raise AnalysisStop(f"Required AnnData obs columns are missing: {sorted(required - set(obs.columns))}")
    donors = obs["donor_id"].astype("string").fillna("").to_numpy(dtype=str)
    statuses = obs["cell_status"].astype("string").fillna("").to_numpy(dtype=str)
    states = obs["annotation_coarse"].astype("string").fillna("").to_numpy(dtype=str)
    if len(donors) != len(obs) or len(statuses) != len(obs) or len(states) != len(obs):
        raise AnalysisStop("Cell metadata arrays have inconsistent lengths")
    malignant = statuses == "Malignant"
    if np.any(malignant & (donors == "")):
        raise AnalysisStop("A malignant nucleus has no donor ID")
    found_donors = tuple(sorted(set(donors[malignant])))
    if found_donors != tuple(sorted(expected_donors)):
        raise AnalysisStop(f"GBM-Space donor set changed: {found_donors}")
    donor_index = {donor: i for i, donor in enumerate(expected_donors)}
    group_codes = np.full(len(obs), -1, dtype=np.int16)
    counts: dict[str, dict[str, int]] = {donor: {"developmental_like": 0, "gliosis_hypoxia": 0} for donor in expected_donors}
    for donor in expected_donors:
        dmask = malignant & (donors == donor)
        dev_mask = dmask & np.isin(states, tuple(DEVELOPMENTAL_STATES))
        injury_mask = dmask & np.isin(states, tuple(INJURY_HYPOXIA_STATES))
        counts[donor]["developmental_like"] = int(dev_mask.sum())
        counts[donor]["gliosis_hypoxia"] = int(injury_mask.sum())
        group_codes[dev_mask] = donor_index[donor] * 2
        group_codes[injury_mask] = donor_index[donor] * 2 + 1
        if counts[donor]["developmental_like"] < minimum_cells:
            raise AnalysisStop(f"{donor} has fewer than {minimum_cells} developmental-like malignant nuclei")
        if counts[donor]["gliosis_hypoxia"] < minimum_cells:
            raise AnalysisStop(f"{donor} has fewer than {minimum_cells} gliosis/hypoxia malignant nuclei")
    total_counts = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(dtype=np.float64)
    if not np.isfinite(total_counts).all() or np.any(total_counts < 0):
        raise AnalysisStop("obs.total_counts contains invalid values")
    return group_codes, malignant, total_counts, counts


def _validate_sparse_count_block(block: Any, expected_totals: np.ndarray) -> sparse.csr_matrix:
    stored_values = np.asarray(
        block.data if sparse.issparse(block) else np.asarray(block).reshape(-1),
        dtype=np.float64,
    )
    if not np.isfinite(stored_values).all() or np.any(stored_values < 0):
        raise AnalysisStop("raw.X contains non-finite or negative stored values")
    if np.any(np.abs(stored_values - np.rint(stored_values)) > 1e-5):
        raise AnalysisStop("raw.X contains non-integer stored values; it is not a raw count matrix")
    matrix = sparse.csr_matrix(block)
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    row_sums = np.asarray(matrix.sum(axis=1), dtype=np.float64).reshape(-1)
    if row_sums.shape != expected_totals.shape or not np.allclose(
        row_sums, expected_totals, rtol=1e-6, atol=1e-3
    ):
        raise AnalysisStop("raw.X row sums do not match obs.total_counts")
    return matrix


def _accumulate_cell_features(matrix: sparse.csr_matrix) -> tuple[np.ndarray, np.ndarray]:
    if matrix.shape[0] == 0:
        return np.zeros(matrix.shape[1], dtype=np.float64), np.zeros(matrix.shape[1], dtype=np.float64)
    row_totals = np.asarray(matrix.sum(axis=1), dtype=np.float64).reshape(-1)
    if np.any(row_totals <= 0):
        raise AnalysisStop("A selected malignant nucleus has zero raw counts")
    row_index = np.repeat(np.arange(matrix.shape[0]), np.diff(matrix.indptr))
    log_cpm_nonzero = np.log2(matrix.data.astype(np.float64) / row_totals[row_index] * 1_000_000.0 + 1.0)
    mean_log_sum = np.bincount(matrix.indices, weights=log_cpm_nonzero, minlength=matrix.shape[1])
    detected = np.bincount(matrix.indices, minlength=matrix.shape[1]).astype(np.float64)
    return mean_log_sum, detected


def aggregate_backed_counts(
    raw_x: Any,
    group_codes: np.ndarray,
    malignant: np.ndarray,
    expected_totals: np.ndarray,
    match_sample_indices: Iterable[int],
    *,
    block_rows: int = BLOCK_ROWS,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    n_obs, n_genes = raw_x.shape
    if not (len(group_codes) == len(malignant) == len(expected_totals) == n_obs):
        raise AnalysisStop("raw.X and observation metadata dimensions do not match")
    if block_rows < 1:
        raise AnalysisStop("block_rows must be positive")
    sample_indices = np.asarray(tuple(match_sample_indices), dtype=np.int64)
    if sample_indices.ndim != 1 or len(sample_indices) == 0 or np.any(sample_indices < 0) or np.any(sample_indices >= n_obs):
        raise AnalysisStop("Expression-match sample indices are invalid")
    if len(np.unique(sample_indices)) != len(sample_indices):
        raise AnalysisStop("Expression-match sample indices must be unique")
    sample_lookup = np.zeros(n_obs, dtype=bool)
    sample_lookup[sample_indices] = True
    active_codes = np.unique(group_codes[group_codes >= 0])
    if active_codes.size == 0 or not np.array_equal(active_codes, np.arange(active_codes[-1] + 1)) or active_codes.size % 2:
        raise AnalysisStop("Complete paired donor-state group codes were not supplied")
    n_donors = int(active_codes.size // 2)
    pseudobulk = np.zeros((n_donors * 2, n_genes), dtype=np.int64)
    feature_log_sum = np.zeros(n_genes, dtype=np.float64)
    feature_detected = np.zeros(n_genes, dtype=np.float64)
    n_match = 0
    for start in range(0, n_obs, block_rows):
        stop = min(start + block_rows, n_obs)
        matrix = _validate_sparse_count_block(raw_x[start:stop, :], expected_totals[start:stop])
        int_matrix = matrix.copy()
        int_matrix.data = np.rint(int_matrix.data).astype(np.int64)
        local_codes = group_codes[start:stop]
        for code in np.unique(local_codes[local_codes >= 0]):
            mask = local_codes == code
            pseudobulk[int(code)] += np.asarray(int_matrix[mask].sum(axis=0), dtype=np.int64).reshape(-1)
        local_match = sample_lookup[start:stop] & malignant[start:stop]
        if np.any(local_match):
            selected = int_matrix[local_match].tocsr()
            log_sum, detected = _accumulate_cell_features(selected)
            feature_log_sum += log_sum
            feature_detected += detected
            n_match += int(selected.shape[0])
    if n_match != len(sample_indices):
        raise AnalysisStop("Expression-match sample did not contain the required malignant nuclei")
    if np.any(pseudobulk.sum(axis=1) <= 0):
        raise AnalysisStop("A donor-state pseudobulk is empty")
    return pseudobulk, feature_log_sum / n_match, feature_detected / n_match, n_match


def rank_percentiles(log_expression: np.ndarray) -> np.ndarray:
    values = np.asarray(log_expression, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
        raise AnalysisStop("Rank input must be a finite one-dimensional vector")
    return (rankdata(values, method="average") - 1.0) / (values.size - 1.0)


def score_counts(counts: np.ndarray, positive_indices: Sequence[int], negative_indices: Sequence[int]) -> float:
    counts = np.asarray(counts, dtype=np.float64)
    if counts.ndim != 1 or not np.isfinite(counts).all() or np.any(counts < 0) or counts.sum() <= 0:
        raise AnalysisStop("Pseudobulk counts must be a finite nonnegative vector with positive total")
    if not positive_indices or not negative_indices or set(positive_indices) & set(negative_indices):
        raise AnalysisStop("The two signed gene groups must be nonempty and disjoint")
    log_cpm = np.log2(counts / counts.sum() * 1_000_000.0 + 1.0)
    ranks = rank_percentiles(log_cpm)
    return float(ranks[list(positive_indices)].mean() - ranks[list(negative_indices)].mean())


def exact_sign_flip(differences: Sequence[float]) -> dict[str, Any]:
    values = np.asarray(differences, dtype=np.float64)
    if values.ndim != 1 or values.size != len(EXPECTED_DONORS) or not np.isfinite(values).all():
        raise AnalysisStop("The primary sign-flip calculation requires 12 finite donor differences")
    observed = float(values.mean())
    exceed = 0
    total = 1 << values.size
    for mask in range(total):
        statistic = float(np.mean(values * np.asarray([1.0 if (mask >> i) & 1 else -1.0 for i in range(values.size)])))
        exceed += int(abs(statistic) >= abs(observed) - 1e-15)
    return {"method": "exact_two_sided_sign_flip", "n_donors": int(values.size), "assignments": total, "mean": observed, "p_value": exceed / total}


def generate_matched_controls(
    mean_log_cpm: np.ndarray,
    detection: np.ndarray,
    gene_ids: Sequence[str],
    panel: Sequence[dict[str, Any]],
    symbol_to_index: dict[str, int],
    *,
    excluded_gene_indices: Iterable[int] = (),
    n_panels: int = CONTROL_PANEL_COUNT,
    seed: int = CONTROL_SEED,
    nearest_pool: int = CONTROL_NEAREST_POOL,
) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    means = np.asarray(mean_log_cpm, dtype=np.float64)
    detects = np.asarray(detection, dtype=np.float64)
    ids = np.asarray([str(value) for value in gene_ids], dtype=str)
    if means.ndim != 1 or detects.shape != means.shape or len(ids) != len(means):
        raise AnalysisStop("Expression-matching features and gene IDs must align")
    if not np.isfinite(means).all() or not np.isfinite(detects).all():
        raise AnalysisStop("Expression-matching features are not finite")
    target_indices = [symbol_to_index[str(row["gene_symbol"])] for row in panel]
    if len(set(target_indices)) != len(panel):
        raise AnalysisStop("Mapped panel indices are not unique")
    excluded_indices = {int(index) for index in excluded_gene_indices}
    if any(index < 0 or index >= len(ids) for index in excluded_indices):
        raise AnalysisStop("An excluded source-panel gene index is out of range")
    excluded_indices.update(target_indices)
    features = np.column_stack((means, detects))
    scale = features.std(axis=0)
    scale[scale == 0] = 1.0
    standardized = (features - features.mean(axis=0)) / scale
    background = np.asarray([i for i in range(len(ids)) if i not in excluded_indices], dtype=np.int64)
    if len(background) < len(panel):
        raise AnalysisStop("Too few background genes to build unique control panels")
    nearest_by_source: list[np.ndarray] = []
    for source_index in target_indices:
        delta = standardized[background] - standardized[source_index]
        distance = np.sqrt(np.square(delta).sum(axis=1))
        order = np.lexsort((ids[background], distance))
        nearest_by_source.append(background[order])
    rng = np.random.Generator(np.random.PCG64DXSM(seed))
    positive_positions = [i for i, row in enumerate(panel) if float(row["rho_viability_100uM"]) > 0]
    negative_positions = [i for i, row in enumerate(panel) if float(row["rho_viability_100uM"]) < 0]
    controls: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for _ in range(n_panels):
        used: set[int] = set()
        assigned = [-1] * len(panel)
        for source_position in rng.permutation(len(panel)):
            available: list[int] = []
            for candidate in nearest_by_source[int(source_position)]:
                candidate_int = int(candidate)
                if candidate_int not in used:
                    available.append(candidate_int)
                    if len(available) == nearest_pool:
                        break
            if not available:
                raise AnalysisStop("Could not form a unique expression-matched control panel")
            chosen = int(rng.choice(available))
            used.add(chosen)
            assigned[int(source_position)] = chosen
        if len(set(assigned)) != len(panel) or min(assigned) < 0:
            raise AnalysisStop("An expression-matched control panel is incomplete")
        controls.append((tuple(assigned[i] for i in positive_positions), tuple(assigned[i] for i in negative_positions)))
    return controls


def _package_versions() -> dict[str, str]:
    output: dict[str, str] = {}
    for package in ("anndata", "h5py", "numpy", "pandas", "scipy"):
        try:
            output[package] = version(package)
        except PackageNotFoundError:
            output[package] = "not-installed"
    return output


def run_analysis(input_path: Path = DEFAULT_INPUT, output_path: Path = DEFAULT_OUTPUT) -> dict[str, Any]:
    input_path = Path(input_path)
    output_path = Path(output_path)
    if output_path.exists():
        raise AnalysisStop(f"Refusing to overwrite existing result: {output_path}")
    if not input_path.is_file():
        raise AnalysisStop(f"Input H5AD does not exist: {input_path}")
    start_time = time.perf_counter()
    frozen_hashes = verify_seal()
    input_hash = verify_input_fingerprint(input_path)
    adata = ad.read_h5ad(input_path, backed="r")
    try:
        if tuple(adata.shape) != EXPECTED_SHAPE:
            raise AnalysisStop(f"H5AD shape changed: {adata.shape}")
        if adata.raw is None or tuple(adata.raw.shape) != EXPECTED_SHAPE:
            raise AnalysisStop("The frozen raw.X count matrix is missing or has the wrong shape")
        obs = adata.obs
        group_codes, malignant, expected_totals, metadata_cell_counts = classify_metadata(obs)
        symbols = adata.raw.var["SYMBOL"].astype("string").to_numpy()
        gene_ids = [str(value) for value in adata.raw.var_names]
        panel = read_source_panel()
        mapped_panel, symbol_to_index = map_shared_panel(panel, symbols, gene_ids)
        all_panel_symbols = {str(row["gene_symbol"]).strip() for row in panel}
        excluded_control_indices = [
            i for i, symbol in enumerate(symbols)
            if not pd.isna(symbol) and str(symbol).strip() in all_panel_symbols
        ]
        positive_indices = [symbol_to_index[row["gene_symbol"]] for row in mapped_panel if float(row["rho_viability_100uM"]) > 0]
        negative_indices = [symbol_to_index[row["gene_symbol"]] for row in mapped_panel if float(row["rho_viability_100uM"]) < 0]
        malignant_rows = np.flatnonzero(malignant)
        if len(malignant_rows) < MATCH_SAMPLE_SIZE:
            raise AnalysisStop("Too few malignant nuclei for the frozen expression-matching sample")
        rng = np.random.Generator(np.random.PCG64DXSM(CONTROL_SEED))
        match_sample_indices = np.sort(rng.choice(malignant_rows, size=MATCH_SAMPLE_SIZE, replace=False))
        pseudobulk, feature_means, feature_detection, match_n = aggregate_backed_counts(
            adata.raw.X,
            group_codes,
            malignant,
            expected_totals,
            match_sample_indices,
        )
        scores: dict[str, dict[str, float]] = {}
        differences: list[float] = []
        donor_rows: list[dict[str, Any]] = []
        control_panels = generate_matched_controls(
            feature_means,
            feature_detection,
            gene_ids,
            mapped_panel,
            symbol_to_index,
            excluded_gene_indices=excluded_control_indices,
        )
        control_differences = np.zeros(CONTROL_PANEL_COUNT, dtype=np.float64)
        n_donors = len(EXPECTED_DONORS)
        for donor_index, donor in enumerate(EXPECTED_DONORS):
            dev_counts = pseudobulk[donor_index * 2]
            injury_counts = pseudobulk[donor_index * 2 + 1]
            dev_score = score_counts(dev_counts, positive_indices, negative_indices)
            injury_score = score_counts(injury_counts, positive_indices, negative_indices)
            difference = injury_score - dev_score
            scores[donor] = {"developmental_like": dev_score, "gliosis_hypoxia": injury_score, "difference": difference}
            differences.append(difference)
            donor_rows.append({
                "donor_id": donor,
                "n_developmental_like": metadata_cell_counts[donor]["developmental_like"],
                "n_gliosis_hypoxia": metadata_cell_counts[donor]["gliosis_hypoxia"],
                "score_developmental_like": dev_score,
                "score_gliosis_hypoxia": injury_score,
                "difference_gliosis_hypoxia_minus_developmental": difference,
            })
            dev_log = np.log2(dev_counts / dev_counts.sum() * 1_000_000.0 + 1.0)
            injury_log = np.log2(injury_counts / injury_counts.sum() * 1_000_000.0 + 1.0)
            dev_ranks = rank_percentiles(dev_log)
            injury_ranks = rank_percentiles(injury_log)
            for panel_index, (positive_control, negative_control) in enumerate(control_panels):
                control_differences[panel_index] += (
                    float(injury_ranks[list(positive_control)].mean() - injury_ranks[list(negative_control)].mean())
                    - float(dev_ranks[list(positive_control)].mean() - dev_ranks[list(negative_control)].mean())
                ) / n_donors
        diff_array = np.asarray(differences, dtype=np.float64)
        mean_difference = float(diff_array.mean())
        standard_error = float(diff_array.std(ddof=1) / np.sqrt(n_donors))
        critical_t = float(t.ppf(0.975, df=n_donors - 1))
        ci = [mean_difference - critical_t * standard_error, mean_difference + critical_t * standard_error]
        result: dict[str, Any] = {
            "status": "exploratory_complete",
            "record_type": "post-analysis aggregate result",
            "question": "Within tumor, does the fixed 115-gene in-vitro TMZ viability-associated score differ between author-defined malignant developmental-like and gliosis/hypoxia groups?",
            "interpretation": "This is an expression association with an in-vitro viability-derived score, not a measured or predicted patient drug response or resistance label.",
            "cohort": {"n_donors": n_donors, "donor_ids": list(EXPECTED_DONORS), "n_malignant_nuclei": int(malignant.sum()), "state_group_cell_counts": metadata_cell_counts},
            "gene_panel": {"source_size": EXPECTED_PANEL_SIZE, "a5_mapped_size": 119, "shared_exact_size": len(mapped_panel), "positive_viability_associated": EXPECTED_POSITIVE, "negative_viability_associated": EXPECTED_NEGATIVE, "omitted_symbols": sorted(EXPECTED_A5_OMITTED | EXPECTED_GBM_SPACE_OMITTED), "substitutions": []},
            "primary": {"estimand": "equal-donor mean of gliosis/hypoxia minus developmental-like score", "mean_difference": mean_difference, "median_difference": float(np.median(diff_array)), "range": [float(diff_array.min()), float(diff_array.max())], "positive_donors": int((diff_array > 0).sum()), "negative_donors": int((diff_array < 0).sum()), "zero_donors": int((diff_array == 0).sum()), "95_percent_t_interval": ci, "exact_sign_flip": exact_sign_flip(diff_array), "donor_results": donor_rows},
            "matched_controls": {"n_panels": CONTROL_PANEL_COUNT, "n_sampled_malignant_nuclei_for_matching": match_n, "seed": CONTROL_SEED, "nearest_pool": CONTROL_NEAREST_POOL, "median_control_mean_difference": float(np.median(control_differences)), "control_2_5_97_5_percentiles": [float(v) for v in np.quantile(control_differences, [0.025, 0.975])], "observed_percentile": float(np.mean(control_differences <= mean_difference))},
            "provenance": {"input_path": str(input_path), "input_bytes": input_path.stat().st_size, "input_sha256": input_hash, "protocol_sha256": sha256_file(PROTOCOL_PATH), "panel_sha256": sha256_file(PANEL_PATH), "seal_sha256": sha256_file(SEAL_PATH), "sealed_file_hashes": frozen_hashes, "git_commit": git_commit(), "python": sys.version, "platform": platform.platform(), "package_versions": _package_versions(), "elapsed_seconds": time.perf_counter() - start_time},
        }
    finally:
        adata.file.close()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    if temporary_path.exists():
        raise AnalysisStop(f"Temporary output already exists: {temporary_path}")
    temporary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary_path.replace(output_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description="Exploratory GBM-Space response-associated state score")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--run", action="store_true", required=True, help="Run once under the frozen protocol")
    args = parser.parse_args()
    if not args.run:
        return 2
    try:
        result = run_analysis(args.input, args.output)
    except (AnalysisStop, OSError, ValueError) as exc:
        print(f"STOP: {exc}", file=sys.stderr, flush=True)
        return 2
    print(json.dumps({"status": result["status"], "n_donors": result["cohort"]["n_donors"], "result_path": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
