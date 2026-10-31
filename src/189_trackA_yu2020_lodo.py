from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import platform
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("OMP_NUM_THREADS", "8")
os.environ.setdefault("MKL_NUM_THREADS", "8")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "8")

import numpy as np
import pandas as pd
import sklearn
import torch
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.preprocessing import OneHotEncoder
from sklearn.svm import LinearSVC
from torch import nn
from torch.nn import functional as F


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROTOCOL_DIR = PROJECT_ROOT / "paper" / "trackA_yu2020"
HASH_FILE = PROTOCOL_DIR / "SHA256SUMS"
OUTPUT_DIR = PROJECT_ROOT / "output" / "trackA_yu2020"
MATRIX_PATH = PROJECT_ROOT / "data" / "external" / "gse117891" / "GSE117891_all_6148.umi.count.matrix.tsv.gz"
CELL_TABLE_PATH = PROJECT_ROOT / "data" / "external" / "trackA_audit_meta" / "GSE117891_Sample_barcode_cell_information.xlsx"
TABLE_S1_PATH = PROJECT_ROOT / "data" / "external" / "yu2020_supp" / "Table S1. ClinicalInformation.xlsx"
GBMAP_PATH = PROJECT_ROOT / "data" / "external" / "gbmap_obs" / "gbmap_core_obs_Yu2020.tsv"
SEED = 20261031
MIN_CELLS_PER_CLASS = 20
MIN_FEATURES = 2_000
MAX_EPOCHS = 500
PATIENCE = 25
INNER_IMPROVEMENT_TOL = 1e-12
LAMBDA_1_GRID = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
LAMBDA_2_GRID = (1e-4, 1e-3, 1e-2, 1e-1, 1.0)
LAMBDA_W_GRID = (1e-3, 1e-2, 1e-1)
BLEND_T_GRID = (0.1, 0.25, 0.5, 1.0)
BLEND_GAMMA_GRID = (0.0, 0.25, 0.5, 0.75, 1.0)
LOGISTIC_C_GRID = (0.001, 0.01, 0.1, 1.0, 10.0)
SVM_C_GRID = LOGISTIC_C_GRID
GBT_RATE_GRID = (0.05, 0.1, 0.2)
K_NEIGHBORS = 15
DEVICE = torch.device("cpu")


class StopRun(RuntimeError):
    """Raised when a frozen analysis gate fails."""


@dataclass
class Dataset:
    counts: np.ndarray
    labels: np.ndarray
    donors: np.ndarray
    cells: np.ndarray
    points: np.ndarray
    cell_types: np.ndarray
    point_sizes: np.ndarray
    gene_ids: np.ndarray
    flow: dict[str, Any]

    def subset(self, included_donors: tuple[str, ...]) -> "Dataset":
        keep = np.isin(self.donors, included_donors)
        return Dataset(
            counts=self.counts[keep],
            labels=self.labels[keep],
            donors=self.donors[keep],
            cells=self.cells[keep],
            points=self.points[keep],
            cell_types=self.cell_types[keep],
            point_sizes=self.point_sizes[keep],
            gene_ids=self.gene_ids,
            flow={**self.flow, "cohort_donors": list(included_donors)},
        )


class ExpressionPreprocessor:
    def __init__(self, selected: np.ndarray, means: np.ndarray, scales: np.ndarray, gene_ids: np.ndarray):
        self.selected = selected
        self.means = means.astype(np.float32, copy=False)
        self.scales = scales.astype(np.float32, copy=False)
        self.gene_ids = gene_ids

    @classmethod
    def fit(cls, train_counts: np.ndarray, gene_ids: np.ndarray) -> tuple["ExpressionPreprocessor", np.ndarray]:
        if train_counts.ndim != 2:
            raise StopRun("Training count matrix is not two-dimensional")
        detected = np.count_nonzero(train_counts, axis=0)
        min_detected = int(math.ceil(0.05 * train_counts.shape[0]))
        eligible = np.flatnonzero(detected >= min_detected)
        if eligible.size < MIN_FEATURES:
            raise StopRun(f"Only {eligible.size} genes pass the training-only 5% filter; need {MIN_FEATURES}")
        totals = train_counts.sum(axis=1, dtype=np.float64)
        if np.any(totals <= 0):
            raise StopRun("A QC-passed training cell has zero library size")
        log_x = np.log1p(train_counts[:, eligible].astype(np.float32) * (10_000.0 / totals[:, None]).astype(np.float32))
        variances = log_x.var(axis=0, dtype=np.float64)
        eligible_genes = gene_ids[eligible].astype(str)
        order = np.lexsort((eligible_genes, -variances))[:MIN_FEATURES]
        selected = eligible[order]
        means = log_x[:, order].mean(axis=0, dtype=np.float64).astype(np.float32)
        scales = log_x[:, order].std(axis=0, dtype=np.float64).astype(np.float32)
        scales[scales == 0] = 1.0
        prep = cls(selected=selected, means=means, scales=scales, gene_ids=gene_ids)
        return prep, prep.transform(train_counts)

    def transform(self, counts: np.ndarray) -> np.ndarray:
        totals = counts.sum(axis=1, dtype=np.float64)
        if np.any(totals <= 0):
            raise StopRun("A transformed cell has zero library size")
        log_x = np.log1p(counts[:, self.selected].astype(np.float32) * (10_000.0 / totals[:, None]).astype(np.float32))
        return np.ascontiguousarray((log_x - self.means) / self.scales, dtype=np.float32)


class GatedClassifier(nn.Module):
    def __init__(self, n_features: int):
        super().__init__()
        self.gate_logits = nn.Parameter(torch.zeros(n_features, dtype=torch.float32))
        self.weights = nn.Parameter(torch.empty((2, n_features), dtype=torch.float32))
        self.bias = nn.Parameter(torch.zeros(2, dtype=torch.float32))
        nn.init.xavier_uniform_(self.weights)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        gates = torch.sigmoid(self.gate_logits)
        logits = F.linear(x * gates, self.weights, self.bias)
        return logits, gates


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_hashes() -> dict[str, str]:
    if not HASH_FILE.is_file():
        raise StopRun(f"Missing checksum manifest: {HASH_FILE}")
    verified: dict[str, str] = {}
    for line_number, raw_line in enumerate(HASH_FILE.read_text(encoding="utf-8").splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split(maxsplit=1)
        if len(fields) != 2 or len(fields[0]) != 64:
            raise StopRun(f"Malformed checksum manifest row {line_number}")
        expected, relative = fields
        relative = relative.lstrip("*").strip()
        path = (HASH_FILE.parent / relative).resolve()
        if not path.is_file():
            raise StopRun(f"Missing sealed file: {path}")
        actual = sha256_file(path)
        if actual.lower() != expected.lower():
            raise StopRun(f"SHA-256 mismatch: {path}")
        verified[str(path.relative_to(PROJECT_ROOT))] = actual
    required = {
        str(p.relative_to(PROJECT_ROOT))
        for p in (PROTOCOL_DIR / "PROTOCOL.md", PROTOCOL_DIR / "donors.tsv", PROJECT_ROOT / "paper" / "research_plan_complete_project_appendix.md", Path(__file__).resolve(), MATRIX_PATH, CELL_TABLE_PATH, TABLE_S1_PATH, GBMAP_PATH)
    }
    absent = sorted(required - set(verified))
    if absent:
        raise StopRun(f"Checksum manifest does not seal required files: {absent}")
    return verified


def current_commit() -> str | None:
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, check=True, capture_output=True, text=True).stdout.strip()
    except Exception:
        return None


def _clean_text(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def load_donor_table() -> tuple[pd.DataFrame, tuple[str, ...], tuple[str, ...], dict[str, str]]:
    path = PROTOCOL_DIR / "donors.tsv"
    donors = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    required = {"patient", "geo_patient", "tier", "tumoral_points", "peritumoral_points"}
    if not required.issubset(donors.columns):
        raise StopRun(f"donors.tsv missing columns: {sorted(required - set(donors.columns))}")
    if donors["patient"].duplicated().any() or donors["geo_patient"].duplicated().any():
        raise StopRun("Donor table has duplicate donor identifiers")
    strict = tuple(donors.loc[donors["tier"].eq("strict"), "patient"].tolist())
    lenient = tuple(donors.loc[donors["tier"].eq("lenient"), "patient"].tolist())
    if len(strict) != 6 or len(lenient) != 2:
        raise StopRun(f"Unexpected donor tiers: {len(strict)} strict, {len(lenient)} lenient")
    geo_to_patient = dict(zip(donors["geo_patient"], donors["patient"], strict=True))
    return donors, strict, strict + lenient, geo_to_patient


def load_metadata(geo_to_patient: dict[str, str], donor_table: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    cell_table = pd.read_excel(CELL_TABLE_PATH, sheet_name="Cell_information", dtype=str)
    needed_cell_cols = {"cell", "patient", "point"}
    if not needed_cell_cols.issubset(cell_table.columns):
        raise StopRun(f"Cell table missing columns: {sorted(needed_cell_cols - set(cell_table.columns))}")
    cell_table = cell_table[["cell", "patient", "point"]].copy()
    for column in cell_table.columns:
        cell_table[column] = cell_table[column].map(_clean_text)
    if cell_table["cell"].eq("").any() or cell_table["cell"].duplicated().any():
        raise StopRun("Cell table has empty or duplicate cell identifiers")

    s1 = pd.read_excel(TABLE_S1_PATH, sheet_name="Point and Cell information", dtype=str)
    needed_s1_cols = {"Sample point", "Sample type"}
    if not needed_s1_cols.issubset(s1.columns):
        raise StopRun(f"Table S1 missing columns: {sorted(needed_s1_cols - set(s1.columns))}")
    point_labels: dict[str, str] = {}
    for _, row in s1.iterrows():
        point = _clean_text(row["Sample point"])
        sample_type = _clean_text(row["Sample type"])
        if not point:
            continue
        if sample_type not in {"Tumoral", "Peritumoral", "Normal"}:
            raise StopRun(f"Unexpected Table S1 sample type at {point}: {sample_type!r}")
        if point in point_labels and point_labels[point] != sample_type:
            raise StopRun(f"Conflicting Table S1 labels for {point}")
        point_labels[point] = sample_type

    donor_point_labels: dict[str, str] = {}
    for _, donor in donor_table.iterrows():
        patient = donor["patient"]
        for listed, expected in ((donor["tumoral_points"], "Tumoral"), (donor["peritumoral_points"], "Peritumoral")):
            for suffix in filter(None, map(str.strip, listed.split(","))):
                key = patient + suffix
                if point_labels.get(key) != expected:
                    raise StopRun(f"Table S1 and donors.tsv disagree for {key}: expected {expected}")
                if key in donor_point_labels and donor_point_labels[key] != expected:
                    raise StopRun(f"donors.tsv assigns conflicting labels for {key}")
                donor_point_labels[key] = expected

    mapped_rows: list[dict[str, Any]] = []
    reasons: dict[str, int] = {"other_patient": 0, "unmapped_point": 0, "normal_point": 0, "unlisted_target_point": 0}
    for row in cell_table.itertuples(index=False):
        cell, geo_patient, geo_point = row.cell, row.patient, row.point
        if geo_patient not in geo_to_patient:
            reasons["other_patient"] += 1
            continue
        if not geo_point.startswith(geo_patient):
            raise StopRun(f"Cell table patient/point conflict for {cell}")
        patient = geo_to_patient[geo_patient]
        point_key = patient + geo_point[len(geo_patient):]
        sample_type = point_labels.get(point_key)
        if sample_type is None:
            reasons["unmapped_point"] += 1
            continue
        if sample_type == "Normal":
            reasons["normal_point"] += 1
            continue
        if sample_type not in {"Tumoral", "Peritumoral"}:
            raise StopRun(f"Unexpected mapped label {sample_type!r} for {cell}")
        if donor_point_labels.get(point_key) != sample_type:
            reasons["unlisted_target_point"] += 1
            continue
        mapped_rows.append({"cell": cell, "donor": patient, "point": point_key, "label": int(sample_type == "Tumoral")})
    mapped = pd.DataFrame(mapped_rows)
    if mapped.empty or mapped["cell"].duplicated().any():
        raise StopRun("No uniquely mapped target cells or duplicate mapped cells")

    gbmap = pd.read_csv(GBMAP_PATH, sep="\t", dtype=str, keep_default_na=False)
    needed_gbmap_cols = {"obs_index", "annotation_level_2"}
    if not needed_gbmap_cols.issubset(gbmap.columns):
        raise StopRun(f"GBmap snapshot missing columns: {sorted(needed_gbmap_cols - set(gbmap.columns))}")
    if gbmap["obs_index"].duplicated().any():
        raise StopRun("GBmap snapshot contains duplicate exact cell identifiers")
    type_map = dict(zip(gbmap["obs_index"], gbmap["annotation_level_2"], strict=True))
    mapped["cell_type"] = mapped["cell"].map(type_map).replace("", np.nan)
    metadata_stats = {
        "cell_table_rows": int(len(cell_table)),
        "s1_labeled_points": int(len(point_labels)),
        "mapped_tumoral_or_peritumoral_cells": int(len(mapped)),
        "mapping_exclusions": reasons,
        "gbmap_snapshot_rows": int(len(gbmap)),
        "gbmap_unmapped_cells": int(mapped["cell_type"].isna().sum()),
        "gbmap_version": "CELLxGENE dataset c888b684-6c51-431f-972a-6c963044cef0; version 861acfd8-25f0-418b-a445-aa96da232827",
    }
    return mapped, metadata_stats


def read_count_matrix(mapped_cells: set[str]) -> tuple[np.ndarray, np.ndarray, list[str], dict[str, int]]:
    with gzip.open(MATRIX_PATH, "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header or len(header) < 3:
        raise StopRun("Count matrix has no usable tab-separated header")
    if len(header[1:]) != len(set(header[1:])):
        raise StopRun("Count matrix has duplicate column identifiers")
    frame = pd.read_csv(MATRIX_PATH, sep="\t", compression="gzip", index_col=0, low_memory=False)
    if frame.index.astype(str).duplicated().any():
        raise StopRun("Count matrix has duplicate row identifiers")
    if frame.columns.astype(str).duplicated().any():
        raise StopRun("Count matrix has duplicate cell or gene column identifiers")
    row_ids = set(frame.index.astype(str))
    col_ids = set(frame.columns.astype(str))
    row_overlap = len(row_ids & mapped_cells)
    col_overlap = len(col_ids & mapped_cells)
    if col_overlap >= 20 and row_overlap < 20:
        orientation = "genes_by_cells"
        selected_cells = [cell for cell in sorted(mapped_cells) if cell in col_ids]
        selected = frame.loc[:, selected_cells]
        gene_ids = frame.index.astype(str).to_numpy()
        counts = selected.to_numpy(copy=True).T
    elif row_overlap >= 20 and col_overlap < 20:
        orientation = "cells_by_genes"
        selected_cells = [cell for cell in sorted(mapped_cells) if cell in row_ids]
        selected = frame.loc[selected_cells, :]
        gene_ids = frame.columns.astype(str).to_numpy()
        counts = selected.to_numpy(copy=True)
    else:
        raise StopRun(f"Count matrix orientation is ambiguous or unmatched (row overlap={row_overlap}, column overlap={col_overlap})")
    try:
        numeric = np.asarray(counts, dtype=np.float64)
    except (TypeError, ValueError) as exc:
        raise StopRun("Count matrix contains nonnumeric entries") from exc
    if not np.isfinite(numeric).all() or np.any(numeric < 0) or not np.equal(numeric, np.floor(numeric)).all():
        raise StopRun("Count matrix must contain finite, nonnegative integer UMI counts")
    if len(gene_ids) != len(set(gene_ids.tolist())):
        raise StopRun("Count matrix has duplicate gene identifiers")
    counts = np.ascontiguousarray(numeric, dtype=np.int32)
    matrix_stats = {
        "orientation": orientation,
        "genes": int(counts.shape[1]),
        "matched_matrix_cells": int(counts.shape[0]),
        "mapped_cells_missing_from_matrix": int(len(mapped_cells - set(selected_cells))),
        "matrix_cell_ids_without_target_metadata": int(len((row_ids if orientation == "cells_by_genes" else col_ids) - mapped_cells)),
    }
    return counts, gene_ids.astype(str), selected_cells, matrix_stats


def load_dataset() -> tuple[Dataset, tuple[str, ...], tuple[str, ...], dict[str, Any], dict[str, str]]:
    hashes = verify_hashes()
    donor_table, strict_donors, sensitivity_donors, geo_to_patient = load_donor_table()
    metadata, metadata_stats = load_metadata(geo_to_patient, donor_table)
    mapped_ids = set(metadata["cell"].astype(str))
    raw_counts, gene_ids, matrix_cells, matrix_stats = read_count_matrix(mapped_ids)
    index_by_cell = {cell: index for index, cell in enumerate(matrix_cells)}
    rows = metadata[metadata["cell"].isin(index_by_cell)].copy()
    row_positions = np.array([index_by_cell[cell] for cell in rows["cell"]], dtype=np.int64)
    counts = raw_counts[row_positions]
    umis = counts.sum(axis=1, dtype=np.int64)
    detected = np.count_nonzero(counts, axis=1)
    qc = (umis >= 1_000) & (detected >= 500)
    rows = rows.loc[qc].reset_index(drop=True)
    counts = counts[qc]
    umis = umis[qc]
    detected = detected[qc]
    point_counts = rows.groupby("point")["cell"].transform("size").to_numpy(dtype=np.int32)
    cell_types = rows["cell_type"].astype(object).where(rows["cell_type"].notna(), None).to_numpy()
    dataset = Dataset(
        counts=counts,
        labels=rows["label"].to_numpy(dtype=np.int64),
        donors=rows["donor"].astype(str).to_numpy(),
        cells=rows["cell"].astype(str).to_numpy(),
        points=rows["point"].astype(str).to_numpy(),
        cell_types=cell_types,
        point_sizes=np.log1p(point_counts.astype(np.float32)),
        gene_ids=gene_ids,
        flow={
            **metadata_stats,
            **{key: value for key, value in matrix_stats.items() if not key.startswith("_")},
            "qc_pass_cells": int(qc.sum()),
            "qc_fail_cells": int((~qc).sum()),
            "qc_rules": {"minimum_umis": 1_000, "minimum_detected_genes": 500},
            "qc_by_donor_class": {},
            "input_hashes": hashes,
        },
    )
    for donor in sensitivity_donors:
        for label, label_name in ((1, "core"), (0, "periphery")):
            count = int(np.sum((dataset.donors == donor) & (dataset.labels == label)))
            dataset.flow["qc_by_donor_class"].setdefault(donor, {})[label_name] = count
    for donor in strict_donors:
        counts_for_donor = dataset.flow["qc_by_donor_class"].get(donor, {})
        if counts_for_donor.get("core", 0) < MIN_CELLS_PER_CLASS or counts_for_donor.get("periphery", 0) < MIN_CELLS_PER_CLASS:
            raise StopRun(f"Primary donor {donor} has fewer than {MIN_CELLS_PER_CLASS} QC-passed cells in one class: {counts_for_donor}")
    return dataset, strict_donors, sensitivity_donors, dataset.flow, hashes


def tie_ranks(cell_ids: np.ndarray) -> np.ndarray:
    keyed = sorted(
        range(len(cell_ids)),
        key=lambda index: hashlib.sha256(f"{SEED}\0{str(cell_ids[index])}".encode("utf-8")).digest(),
    )
    ranks = np.empty(len(cell_ids), dtype=np.int64)
    ranks[np.asarray(keyed, dtype=np.int64)] = np.arange(len(cell_ids), dtype=np.int64)
    return ranks


def normalize_rows(x: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(x, axis=1, keepdims=True)
    return np.divide(x, norms, out=np.zeros_like(x), where=norms > 0)


def graph_edges(x: np.ndarray, tie: np.ndarray, k: int = K_NEIGHBORS) -> np.ndarray:
    n = len(x)
    if n < 2:
        return np.empty((0, 2), dtype=np.int64)
    unit = normalize_rows(x)
    similarity = unit @ unit.T
    np.fill_diagonal(similarity, -np.inf)
    edges: set[tuple[int, int]] = set()
    neighbors = min(k, n - 1)
    for i in range(n):
        order = np.lexsort((tie, -similarity[i]))
        chosen = [j for j in order if j != i][:neighbors]
        for j in chosen:
            edges.add((min(i, int(j)), max(i, int(j))))
    if not edges:
        return np.empty((0, 2), dtype=np.int64)
    return np.asarray(sorted(edges), dtype=np.int64)


def query_neighbors(x_query: np.ndarray, x_train: np.ndarray, train_tie: np.ndarray, k: int = K_NEIGHBORS) -> tuple[np.ndarray, np.ndarray]:
    if x_train.shape[0] == 0:
        raise StopRun("Neighbor query has no training cells")
    sims = normalize_rows(x_query) @ normalize_rows(x_train).T
    neighbors = min(k, x_train.shape[0])
    indices = np.empty((len(x_query), neighbors), dtype=np.int64)
    values = np.empty((len(x_query), neighbors), dtype=np.float32)
    for i in range(len(x_query)):
        order = np.lexsort((train_tie, -sims[i]))[:neighbors]
        indices[i] = order
        values[i] = sims[i, order]
    return indices, values


def neighbor_probabilities(x_query: np.ndarray, x_train: np.ndarray, train_tie: np.ndarray, p_train: np.ndarray, temperatures: tuple[float, ...]) -> dict[float, np.ndarray]:
    indices, similarities = query_neighbors(x_query, x_train, train_tie)
    output: dict[float, np.ndarray] = {}
    gathered = p_train[indices]
    for temperature in temperatures:
        logits = similarities / temperature
        logits -= logits.max(axis=1, keepdims=True)
        weights = np.exp(logits)
        weights /= weights.sum(axis=1, keepdims=True)
        output[temperature] = np.sum(gathered * weights[:, :, None], axis=1)
    return output


def class_balanced_aurc(y: np.ndarray, probabilities: np.ndarray, tie: np.ndarray) -> float:
    y = np.asarray(y, dtype=np.int64)
    probabilities = checked_probabilities(probabilities, len(y))
    p_core = probabilities[:, 1]
    confidence = np.maximum(p_core, 1.0 - p_core)
    order = np.lexsort((tie, -confidence))
    accepted_y = y[order]
    predicted = (p_core[order] >= 0.5).astype(np.int64)
    errors = (predicted != accepted_y).astype(np.int64)
    n = len(y)
    if n < 2:
        raise StopRun("AURC requires at least two scored cells")
    k_min = int(math.ceil(n / 2.0))
    k_values = np.arange(k_min, n + 1, dtype=np.int64)
    risks = np.empty(len(k_values), dtype=np.float64)
    c0 = np.cumsum(accepted_y == 0)
    c1 = np.cumsum(accepted_y == 1)
    e0 = np.cumsum(errors * (accepted_y == 0))
    e1 = np.cumsum(errors * (accepted_y == 1))
    for i, k in enumerate(k_values):
        risk0 = e0[k - 1] / c0[k - 1] if c0[k - 1] else 1.0
        risk1 = e1[k - 1] / c1[k - 1] if c1[k - 1] else 1.0
        risks[i] = 0.5 * (risk0 + risk1)
    coverages = k_values.astype(np.float64) / n
    return float(np.trapezoid(risks, coverages))


def metric_record(y: np.ndarray, p: np.ndarray, tie: np.ndarray) -> dict[str, float]:
    p = checked_probabilities(p, len(y))
    p_core = p[:, 1]
    prediction = (p_core >= 0.5).astype(np.int64)
    return {
        "class_balanced_aurc": class_balanced_aurc(y, p, tie),
        "balanced_accuracy": float(balanced_accuracy_score(y, prediction)),
        "brier": float(np.mean((p_core - y) ** 2)),
    }


def subgroup_metric_record(y: np.ndarray, p: np.ndarray, tie: np.ndarray) -> dict[str, Any]:
    y = np.asarray(y, dtype=np.int64)
    p = checked_probabilities(p, len(y))
    if len(y) == 0:
        return {"n": 0, "n_core": 0, "n_periphery": 0, "class_balanced_aurc": None, "balanced_accuracy": None, "brier": None}
    p_core = p[:, 1]
    return {
        "n": int(len(y)),
        "n_core": int(np.sum(y == 1)),
        "n_periphery": int(np.sum(y == 0)),
        "class_balanced_aurc": class_balanced_aurc(y, p, tie) if len(y) >= 2 else None,
        "balanced_accuracy": float(balanced_accuracy_score(y, (p_core >= 0.5).astype(np.int64))),
        "brier": float(np.mean((p_core - y) ** 2)),
    }


def checked_probabilities(p: np.ndarray, n_rows: int) -> np.ndarray:
    p = np.asarray(p, dtype=np.float64)
    if p.shape != (n_rows, 2):
        raise StopRun(f"Probability array has shape {p.shape}; expected {(n_rows, 2)}")
    if not np.isfinite(p).all() or np.any(p < 0.0) or np.any(p > 1.0):
        raise StopRun("Model produced invalid probabilities")
    if len(p) and not np.allclose(p.sum(axis=1), 1.0, rtol=0.0, atol=1e-5):
        raise StopRun("Model probability rows do not sum to one")
    return p


def first_setting_within_tolerance(scores: np.ndarray, tolerance: float = INNER_IMPROVEMENT_TOL) -> tuple[int, ...]:
    scores = np.asarray(scores, dtype=np.float64)
    if scores.size == 0 or not np.isfinite(scores).all():
        raise StopRun("Hyperparameter search produced empty or non-finite AURC scores")
    minimum = float(np.min(scores))
    first = int(np.flatnonzero(scores.ravel() <= minimum + tolerance)[0])
    return tuple(int(index) for index in np.unravel_index(first, scores.shape))


def _to_tensor(x: np.ndarray) -> torch.Tensor:
    return torch.as_tensor(np.asarray(x, dtype=np.float32), device=DEVICE)


def train_hybrid(
    x_train: np.ndarray,
    y_train: np.ndarray,
    donor_train: np.ndarray,
    edges: np.ndarray,
    lambda_1: float,
    lambda_2: float,
    lambda_w: float,
    loss_kind: str,
    seed: int,
    x_validation: np.ndarray | None = None,
    y_validation: np.ndarray | None = None,
    fixed_epochs: int | None = None,
) -> tuple[GatedClassifier, int]:
    torch.manual_seed(seed)
    x_t = _to_tensor(x_train)
    y_t = torch.as_tensor(y_train, dtype=torch.long, device=DEVICE)
    donor_names = tuple(dict.fromkeys(map(str, donor_train)))
    donor_tensors = [torch.as_tensor(np.flatnonzero(donor_train == donor), dtype=torch.long, device=DEVICE) for donor in donor_names]
    if not donor_tensors or any(len(idx) == 0 for idx in donor_tensors):
        raise StopRun("Candidate fit has an empty donor group")
    if edges.size:
        edge_t = torch.as_tensor(edges, dtype=torch.long, device=DEVICE)
    else:
        edge_t = torch.empty((0, 2), dtype=torch.long, device=DEVICE)
    model = GatedClassifier(x_train.shape[1]).to(DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    epochs = fixed_epochs if fixed_epochs is not None else MAX_EPOCHS
    best_epoch = epochs
    best_validation = float("inf")
    best_state = None
    patience = 0
    for epoch in range(epochs):
        model.train()
        optimizer.zero_grad(set_to_none=True)
        logits, gates = model(x_t)
        donor_losses = torch.stack([F.cross_entropy(logits.index_select(0, idx), y_t.index_select(0, idx)) for idx in donor_tensors])
        if loss_kind == "smooth_worst_donor":
            data_loss = 0.1 * (torch.logsumexp(donor_losses / 0.1, dim=0) - math.log(len(donor_tensors)))
        elif loss_kind == "mean_donor":
            data_loss = donor_losses.mean()
        else:
            raise StopRun(f"Unknown candidate loss: {loss_kind}")
        ridge = 0.5 * lambda_w * torch.sum(model.weights * model.weights)
        gate_penalty = lambda_1 * torch.sum(gates)
        if lambda_2 and edge_t.shape[0]:
            probs = torch.softmax(logits, dim=1)
            edge_delta = probs[edge_t[:, 0]] - probs[edge_t[:, 1]]
            graph_penalty = lambda_2 * torch.sum(edge_delta * edge_delta) / edge_t.shape[0]
        else:
            graph_penalty = torch.zeros((), dtype=torch.float32, device=DEVICE)
        objective = data_loss + ridge + gate_penalty + graph_penalty
        if not torch.isfinite(objective):
            raise StopRun("Candidate training produced a non-finite objective")
        objective.backward()
        optimizer.step()

        if fixed_epochs is not None:
            continue
        if x_validation is None or y_validation is None:
            raise StopRun("Early-stopped candidate fit lacks an inner-validation donor")
        model.eval()
        with torch.no_grad():
            val_logits, _ = model(_to_tensor(x_validation))
            val_loss = float(F.cross_entropy(val_logits, torch.as_tensor(y_validation, dtype=torch.long, device=DEVICE)).item())
        if not math.isfinite(val_loss):
            raise StopRun("Candidate training produced a non-finite validation loss")
        if val_loss < best_validation - INNER_IMPROVEMENT_TOL:
            best_validation = val_loss
            best_epoch = epoch + 1
            best_state = {key: value.detach().cpu().clone() for key, value in model.state_dict().items()}
            patience = 0
        else:
            patience += 1
            if patience >= PATIENCE:
                break
    if fixed_epochs is None and best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    return model, int(best_epoch)


def hybrid_probabilities(model: GatedClassifier, x: np.ndarray) -> np.ndarray:
    model.eval()
    with torch.no_grad():
        probs = torch.softmax(model(_to_tensor(x))[0], dim=1).cpu().numpy()
    return probs.astype(np.float32, copy=False)


def make_blend_settings(ablation: str) -> list[tuple[float, float]]:
    if ablation == "ablation_graph_and_blend":
        return [(0.1, 0.0)]
    return [(temperature, gamma) for temperature in BLEND_T_GRID for gamma in BLEND_GAMMA_GRID]


def candidate_configs(ablation: str) -> list[tuple[float, float, float, str]]:
    if ablation == "ablation_no_gate":
        return [(0.0, l2, lw, "smooth_worst_donor") for l2 in LAMBDA_2_GRID for lw in LAMBDA_W_GRID]
    if ablation == "ablation_graph_and_blend":
        return [(l1, 0.0, lw, "smooth_worst_donor") for l1 in LAMBDA_1_GRID for lw in LAMBDA_W_GRID]
    loss_kind = "mean_donor" if ablation == "ablation_mean_donor" else "smooth_worst_donor"
    return [(l1, l2, lw, loss_kind) for l1 in LAMBDA_1_GRID for l2 in LAMBDA_2_GRID for lw in LAMBDA_W_GRID]


def tune_candidate(
    data: Dataset,
    outer_train_idx: np.ndarray,
    outer_train_donors: tuple[str, ...],
    global_tie: np.ndarray,
    ablation: str,
) -> tuple[dict[str, Any], dict[str, Any]]:
    configs = candidate_configs(ablation)
    blends = make_blend_settings(ablation)
    scores = np.zeros((len(configs), len(blends)), dtype=np.float64)
    epochs = np.zeros((len(configs), len(outer_train_donors)), dtype=np.int32)
    split_timings: list[float] = []
    for fold_index, inner_donor in enumerate(outer_train_donors):
        started = time.perf_counter()
        inner_val = outer_train_idx[data.donors[outer_train_idx] == inner_donor]
        inner_train = outer_train_idx[data.donors[outer_train_idx] != inner_donor]
        inner_train_donors = tuple(d for d in outer_train_donors if d != inner_donor)
        if set(data.labels[inner_train].tolist()) != {0, 1} or set(data.labels[inner_val].tolist()) != {0, 1}:
            raise StopRun(f"Inner split lacks both classes for {inner_donor}")
        prep, x_train = ExpressionPreprocessor.fit(data.counts[inner_train], data.gene_ids)
        x_val = prep.transform(data.counts[inner_val])
        y_train, y_val = data.labels[inner_train], data.labels[inner_val]
        d_train = data.donors[inner_train]
        tie_train, tie_val = global_tie[inner_train], global_tie[inner_val]
        edges = graph_edges(x_train, tie_train) if ablation != "ablation_graph_and_blend" else np.empty((0, 2), dtype=np.int64)
        for config_index, (lambda_1, lambda_2, lambda_w, loss_kind) in enumerate(configs):
            model, best_epoch = train_hybrid(
                x_train, y_train, d_train, edges, lambda_1, lambda_2, lambda_w, loss_kind,
                seed=SEED, x_validation=x_val, y_validation=y_val,
            )
            epochs[config_index, fold_index] = best_epoch
            p_train = hybrid_probabilities(model, x_train)
            p_val = hybrid_probabilities(model, x_val)
            if len(blends) == 1 and blends[0][1] == 0.0:
                scores[config_index, 0] += class_balanced_aurc(y_val, p_val, tie_val)
                continue
            neighbor_by_t = neighbor_probabilities(x_val, x_train, tie_train, p_train, BLEND_T_GRID)
            for blend_index, (temperature, gamma) in enumerate(blends):
                p = (1.0 - gamma) * p_val + gamma * neighbor_by_t[temperature]
                scores[config_index, blend_index] += class_balanced_aurc(y_val, p, tie_val)
        split_timings.append(time.perf_counter() - started)
    scores /= len(outer_train_donors)
    best_config, best_blend = first_setting_within_tolerance(scores)
    lambda_1, lambda_2, lambda_w, loss_kind = configs[int(best_config)]
    temperature, gamma = blends[int(best_blend)]
    selected = {
        "lambda_1": lambda_1,
        "lambda_2": lambda_2,
        "lambda_W": lambda_w,
        "loss": loss_kind,
        "temperature": temperature,
        "gamma": gamma,
        "inner_mean_aurc": float(scores[best_config, best_blend]),
        "median_best_epoch": int(np.median(epochs[best_config])),
        "best_epochs_by_inner_donor": epochs[best_config].astype(int).tolist(),
        "inner_donor_order": list(outer_train_donors),
    }
    timing = {"inner_fold_seconds": split_timings, "total_seconds": float(sum(split_timings)), "fit_count": int(len(configs) * len(outer_train_donors))}
    return selected, timing


def fit_outer_candidate(data: Dataset, train_idx: np.ndarray, test_idx: np.ndarray, tie: np.ndarray, selected: dict[str, Any], ablation: str) -> tuple[np.ndarray, float]:
    prep, x_train = ExpressionPreprocessor.fit(data.counts[train_idx], data.gene_ids)
    x_test = prep.transform(data.counts[test_idx])
    y_train = data.labels[train_idx]
    edges = graph_edges(x_train, tie[train_idx]) if ablation != "ablation_graph_and_blend" else np.empty((0, 2), dtype=np.int64)
    started = time.perf_counter()
    model, _ = train_hybrid(
        x_train, y_train, data.donors[train_idx], edges,
        float(selected["lambda_1"]), float(selected["lambda_2"]), float(selected["lambda_W"]), str(selected["loss"]),
        seed=SEED, fixed_epochs=int(selected["median_best_epoch"]),
    )
    p_train = hybrid_probabilities(model, x_train)
    p_test = hybrid_probabilities(model, x_test)
    gamma = float(selected["gamma"])
    if gamma:
        neighbor = neighbor_probabilities(x_test, x_train, tie[train_idx], p_train, (float(selected["temperature"]),))[float(selected["temperature"])]
        p_test = (1.0 - gamma) * p_test + gamma * neighbor
    return p_test, time.perf_counter() - started


def _fit_sklearn(model: Any, x: Any, y: np.ndarray, sample_weight: np.ndarray | None = None) -> Any:
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        if sample_weight is None:
            model.fit(x, y)
        else:
            model.fit(x, y, sample_weight=sample_weight)
    if any(issubclass(item.category, ConvergenceWarning) for item in caught):
        raise StopRun(f"Required estimator did not converge: {model.__class__.__name__}")
    return model


def fit_predict_logistic(x_train: np.ndarray, y_train: np.ndarray, x_eval: np.ndarray, c_value: float, balanced: bool = False) -> np.ndarray:
    model = LogisticRegression(C=c_value, penalty="l2", solver="lbfgs", class_weight="balanced" if balanced else None, max_iter=1_000, random_state=SEED)
    _fit_sklearn(model, x_train, y_train)
    return model.predict_proba(x_eval).astype(np.float32)


def _fit_svm(x: np.ndarray, y: np.ndarray, c_value: float) -> LinearSVC:
    model = LinearSVC(C=c_value, loss="squared_hinge", penalty="l2", dual=True, class_weight=None, max_iter=10_000, random_state=SEED)
    return _fit_sklearn(model, x, y)


def fit_predict_svm_calibrated(x_train: np.ndarray, y_train: np.ndarray, donor_train: np.ndarray, x_eval: np.ndarray, c_value: float) -> np.ndarray:
    donors = tuple(dict.fromkeys(map(str, donor_train)))
    if len(donors) < 3:
        raise StopRun("SVM sigmoid calibration needs at least three training donors")
    oof_margin = np.full(len(y_train), np.nan, dtype=np.float64)
    for held_donor in donors:
        fit = donor_train != held_donor
        hold = donor_train == held_donor
        if set(y_train[fit].tolist()) != {0, 1}:
            raise StopRun(f"SVM calibration fold without both classes: {held_donor}")
        fold_model = _fit_svm(x_train[fit], y_train[fit], c_value)
        oof_margin[hold] = fold_model.decision_function(x_train[hold])
    if not np.isfinite(oof_margin).all():
        raise StopRun("SVM cross-fitted margins are incomplete")
    donor_weights = np.zeros(len(y_train), dtype=np.float64)
    for donor in donors:
        mask = donor_train == donor
        donor_weights[mask] = 1.0 / (len(donors) * int(mask.sum()))
    calibrator = LogisticRegression(C=1.0, solver="lbfgs", max_iter=1_000, class_weight=None, random_state=SEED)
    _fit_sklearn(calibrator, oof_margin.reshape(-1, 1), y_train, sample_weight=donor_weights)
    final_model = _fit_svm(x_train, y_train, c_value)
    eval_margin = final_model.decision_function(x_eval).reshape(-1, 1)
    return calibrator.predict_proba(eval_margin).astype(np.float32)


def fit_predict_gbt(x_train: np.ndarray, y_train: np.ndarray, x_eval: np.ndarray, learning_rate: float) -> np.ndarray:
    model = GradientBoostingClassifier(
        loss="log_loss", learning_rate=learning_rate, n_estimators=100, max_depth=2,
        min_samples_leaf=10, subsample=1.0, random_state=SEED,
    )
    _fit_sklearn(model, x_train, y_train)
    return model.predict_proba(x_eval).astype(np.float32)


def fit_predict_knn(x_train: np.ndarray, y_train: np.ndarray, tie_train: np.ndarray, x_eval: np.ndarray) -> np.ndarray:
    indices, _ = query_neighbors(x_eval, x_train, tie_train)
    p_core = y_train[indices].mean(axis=1)
    return np.column_stack((1.0 - p_core, p_core)).astype(np.float32)


def tune_logistic(data: Dataset, train_idx: np.ndarray, train_donors: tuple[str, ...], global_tie: np.ndarray) -> tuple[float, float]:
    scores = np.zeros(len(LOGISTIC_C_GRID), dtype=np.float64)
    started = time.perf_counter()
    for inner_donor in train_donors:
        val_idx = train_idx[data.donors[train_idx] == inner_donor]
        fit_idx = train_idx[data.donors[train_idx] != inner_donor]
        prep, x_fit = ExpressionPreprocessor.fit(data.counts[fit_idx], data.gene_ids)
        x_val = prep.transform(data.counts[val_idx])
        for i, c_value in enumerate(LOGISTIC_C_GRID):
            p = fit_predict_logistic(x_fit, data.labels[fit_idx], x_val, c_value)
            scores[i] += class_balanced_aurc(data.labels[val_idx], p, global_tie[val_idx])
    scores /= len(train_donors)
    best = first_setting_within_tolerance(scores)[0]
    return LOGISTIC_C_GRID[best], time.perf_counter() - started


def tune_svm(data: Dataset, train_idx: np.ndarray, train_donors: tuple[str, ...], global_tie: np.ndarray) -> tuple[float, float]:
    scores = np.zeros(len(SVM_C_GRID), dtype=np.float64)
    started = time.perf_counter()
    for inner_donor in train_donors:
        val_idx = train_idx[data.donors[train_idx] == inner_donor]
        fit_idx = train_idx[data.donors[train_idx] != inner_donor]
        prep, x_fit = ExpressionPreprocessor.fit(data.counts[fit_idx], data.gene_ids)
        x_val = prep.transform(data.counts[val_idx])
        for i, c_value in enumerate(SVM_C_GRID):
            p = fit_predict_svm_calibrated(x_fit, data.labels[fit_idx], data.donors[fit_idx], x_val, c_value)
            scores[i] += class_balanced_aurc(data.labels[val_idx], p, global_tie[val_idx])
    scores /= len(train_donors)
    best = first_setting_within_tolerance(scores)[0]
    return SVM_C_GRID[best], time.perf_counter() - started


def tune_gbt(data: Dataset, train_idx: np.ndarray, train_donors: tuple[str, ...], global_tie: np.ndarray) -> tuple[float, float]:
    scores = np.zeros(len(GBT_RATE_GRID), dtype=np.float64)
    started = time.perf_counter()
    for inner_donor in train_donors:
        val_idx = train_idx[data.donors[train_idx] == inner_donor]
        fit_idx = train_idx[data.donors[train_idx] != inner_donor]
        prep, x_fit = ExpressionPreprocessor.fit(data.counts[fit_idx], data.gene_ids)
        x_val = prep.transform(data.counts[val_idx])
        for i, rate in enumerate(GBT_RATE_GRID):
            p = fit_predict_gbt(x_fit, data.labels[fit_idx], x_val, rate)
            scores[i] += class_balanced_aurc(data.labels[val_idx], p, global_tie[val_idx])
    scores /= len(train_donors)
    best = first_setting_within_tolerance(scores)[0]
    return GBT_RATE_GRID[best], time.perf_counter() - started


def tune_baselines(data: Dataset, train_idx: np.ndarray, train_donors: tuple[str, ...], global_tie: np.ndarray) -> tuple[dict[str, Any], dict[str, float]]:
    logistic_c, logistic_sec = tune_logistic(data, train_idx, train_donors, global_tie)
    svm_c, svm_sec = tune_svm(data, train_idx, train_donors, global_tie)
    gbt_rate, gbt_sec = tune_gbt(data, train_idx, train_donors, global_tie)
    return {"logistic_C": logistic_c, "svm_C": svm_c, "gbt_learning_rate": gbt_rate}, {"logistic_tuning_seconds": logistic_sec, "svm_tuning_seconds": svm_sec, "gbt_tuning_seconds": gbt_sec}


def fit_diagnostics(data: Dataset, train_idx: np.ndarray, test_idx: np.ndarray) -> tuple[np.ndarray | None, np.ndarray, dict[str, Any]]:
    train_types = data.cell_types[train_idx]
    test_types = data.cell_types[test_idx]
    train_mask = np.array([value is not None and not pd.isna(value) for value in train_types])
    test_mask = np.array([value is not None and not pd.isna(value) for value in test_types])
    p_type = None
    type_stats: dict[str, Any] = {"training_mapped": int(train_mask.sum()), "test_mapped": int(test_mask.sum()), "test_unmapped": int((~test_mask).sum())}
    if train_mask.any() and set(data.labels[train_idx][train_mask].tolist()) == {0, 1}:
        encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=True, dtype=np.float32)
        xtr = encoder.fit_transform(np.asarray(train_types[train_mask], dtype=str).reshape(-1, 1))
        xte = encoder.transform(np.asarray(test_types[test_mask], dtype=str).reshape(-1, 1)) if test_mask.any() else None
        model = LogisticRegression(C=1.0, penalty="l2", solver="lbfgs", class_weight="balanced", max_iter=1_000, random_state=SEED)
        _fit_sklearn(model, xtr, data.labels[train_idx][train_mask])
        if xte is not None:
            p_type = model.predict_proba(xte).astype(np.float32)
    point_model = LogisticRegression(C=1.0, penalty="l2", solver="lbfgs", class_weight="balanced", max_iter=1_000, random_state=SEED)
    _fit_sklearn(point_model, data.point_sizes[train_idx].reshape(-1, 1), data.labels[train_idx])
    p_point = point_model.predict_proba(data.point_sizes[test_idx].reshape(-1, 1)).astype(np.float32)
    return p_type, p_point, type_stats


def fit_final_baselines(data: Dataset, train_idx: np.ndarray, test_idx: np.ndarray, global_tie: np.ndarray, selected: dict[str, Any]) -> dict[str, np.ndarray]:
    prep, x_train = ExpressionPreprocessor.fit(data.counts[train_idx], data.gene_ids)
    x_test = prep.transform(data.counts[test_idx])
    y_train = data.labels[train_idx]
    outputs = {
        "logistic": fit_predict_logistic(x_train, y_train, x_test, float(selected["logistic_C"])),
        "knn15": fit_predict_knn(x_train, y_train, global_tie[train_idx], x_test),
        "svm": fit_predict_svm_calibrated(x_train, y_train, data.donors[train_idx], x_test, float(selected["svm_C"])),
        "gradient_boosting": fit_predict_gbt(x_train, y_train, x_test, float(selected["gbt_learning_rate"])),
    }
    return outputs


def fit_outer_models_for_timing(data: Dataset, train_idx: np.ndarray, global_tie: np.ndarray, selected: dict[str, Any]) -> dict[str, float]:
    prep, x_train = ExpressionPreprocessor.fit(data.counts[train_idx], data.gene_ids)
    y_train = data.labels[train_idx]
    timing: dict[str, float] = {}
    for name in ("candidate", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor"):
        params = selected[name]
        edges = graph_edges(x_train, global_tie[train_idx]) if name != "ablation_graph_and_blend" else np.empty((0, 2), dtype=np.int64)
        started = time.perf_counter()
        train_hybrid(
            x_train, y_train, data.donors[train_idx], edges,
            float(params["lambda_1"]), float(params["lambda_2"]), float(params["lambda_W"]), str(params["loss"]),
            seed=SEED, fixed_epochs=int(params["median_best_epoch"]),
        )
        timing[name] = time.perf_counter() - started

    baseline_started = time.perf_counter()
    logistic = LogisticRegression(C=float(selected["baselines"]["logistic_C"]), penalty="l2", solver="lbfgs", class_weight=None, max_iter=1_000, random_state=SEED)
    _fit_sklearn(logistic, x_train, y_train)
    fit_predict_svm_calibrated(x_train, y_train, data.donors[train_idx], x_train, float(selected["baselines"]["svm_C"]))
    gbt = GradientBoostingClassifier(
        loss="log_loss", learning_rate=float(selected["baselines"]["gbt_learning_rate"]), n_estimators=100,
        max_depth=2, min_samples_leaf=10, subsample=1.0, random_state=SEED,
    )
    _fit_sklearn(gbt, x_train, y_train)
    timing["baselines_and_calibration"] = time.perf_counter() - baseline_started

    diagnostic_started = time.perf_counter()
    fit_diagnostics(data, train_idx, train_idx)
    timing["leakage_diagnostics"] = time.perf_counter() - diagnostic_started
    return timing


def summarize_donor_metrics(per_donor: list[dict[str, Any]], field: str) -> dict[str, float | None]:
    values = [float(row[field]) for row in per_donor if row.get(field) is not None and math.isfinite(float(row[field]))]
    if not values:
        return {"mean": None, "min": None, "max": None}
    return {"mean": float(np.mean(values)), "min": float(np.min(values)), "max": float(np.max(values))}


def run_outer_fold(data: Dataset, donors: tuple[str, ...], test_donor: str, global_tie: np.ndarray, score_test: bool) -> tuple[dict[str, Any], float]:
    started = time.perf_counter()
    all_idx = np.arange(len(data.labels), dtype=np.int64)
    test_idx = all_idx[data.donors == test_donor]
    train_idx = all_idx[data.donors != test_donor]
    train_donors = tuple(donor for donor in donors if donor != test_donor)
    if set(data.labels[test_idx].tolist()) != {0, 1} or set(data.labels[train_idx].tolist()) != {0, 1}:
        raise StopRun(f"Outer fold {test_donor} lacks both classes")
    selected: dict[str, Any] = {}
    timing: dict[str, Any] = {}
    for name in ("candidate", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor"):
        params, elapsed = tune_candidate(data, train_idx, train_donors, global_tie, name)
        selected[name] = params
        timing[f"{name}_tuning_seconds"] = elapsed["total_seconds"]
        timing[f"{name}_inner_fit_count"] = elapsed["fit_count"]
    baseline_selected, baseline_timing = tune_baselines(data, train_idx, train_donors, global_tie)
    selected["baselines"] = baseline_selected
    timing.update(baseline_timing)
    if not score_test:
        # Pilot fits only on outer-training donors; it never predicts on held-out cells.
        timing["outer_training_refits_seconds"] = fit_outer_models_for_timing(data, train_idx, global_tie, selected)
        return {"timing": timing, "n_train": int(len(train_idx)), "n_test": int(len(test_idx))}, time.perf_counter() - started

    probabilities: dict[str, np.ndarray] = {}
    model_fit_sec: dict[str, float] = {}
    for name in ("candidate", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor"):
        p, elapsed = fit_outer_candidate(data, train_idx, test_idx, global_tie, selected[name], name)
        probabilities[name] = p
        model_fit_sec[name] = elapsed
    t0 = time.perf_counter()
    probabilities.update(fit_final_baselines(data, train_idx, test_idx, global_tie, selected["baselines"]))
    model_fit_sec["baselines_fit_and_calibration_sec"] = time.perf_counter() - t0
    p_type, p_point, type_stats = fit_diagnostics(data, train_idx, test_idx)
    if p_type is not None:
        probabilities["cell_type_only"] = p_type
    probabilities["point_size_only"] = p_point
    timing["final_model_fit_seconds"] = model_fit_sec

    test_y = data.labels[test_idx]
    test_tie = global_tie[test_idx]
    test_types = data.cell_types[test_idx]
    mapped_mask = np.array([value is not None and not pd.isna(value) for value in test_types])
    metrics: dict[str, Any] = {}
    for name, p in probabilities.items():
        y = test_y
        tie = test_tie
        if name == "cell_type_only":
            y, tie = test_y[mapped_mask], test_tie[mapped_mask]
            if len(y) == 0 or set(y.tolist()) != {0, 1}:
                metrics[name] = {"n": int(len(y)), "n_core": int(np.sum(y == 1)), "n_periphery": int(np.sum(y == 0)), "balanced_accuracy": None}
                continue
        record = metric_record(y, p, tie)
        record.update({"n": int(len(y)), "n_core": int(np.sum(y == 1)), "n_periphery": int(np.sum(y == 0))})
        if name in {"candidate", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor", "logistic", "knn15", "svm", "gradient_boosting"}:
            if mapped_mask.any():
                if p.shape[0] == len(test_y):
                    mapped_record = subgroup_metric_record(test_y[mapped_mask], p[mapped_mask], test_tie[mapped_mask])
                    record["mapped_subset_metrics"] = mapped_record
                    has_both_mapped_classes = set(test_y[mapped_mask].tolist()) == {0, 1}
                    record["mapped_subset_balanced_accuracy"] = mapped_record["balanced_accuracy"] if has_both_mapped_classes else None
                else:
                    record["mapped_subset_metrics"] = None
                    record["mapped_subset_balanced_accuracy"] = None
            else:
                record["mapped_subset_metrics"] = None
                record["mapped_subset_balanced_accuracy"] = None
        metrics[name] = record
    expression_models = ("candidate", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor", "logistic", "knn15", "svm", "gradient_boosting")
    within_cell_type: dict[str, Any] = {}
    valid_types = sorted({str(value) for value in test_types if value is not None and not pd.isna(value)})
    for cell_type in valid_types:
        type_mask = np.asarray([value is not None and not pd.isna(value) and str(value) == cell_type for value in test_types])
        within_cell_type[cell_type] = {
            model_name: subgroup_metric_record(test_y[type_mask], probabilities[model_name][type_mask], test_tie[type_mask])
            for model_name in expression_models
        }
    metrics["within_cell_type"] = within_cell_type
    metrics["cell_type_mapping"] = type_stats
    return {"test_donor": test_donor, "n_train": int(len(train_idx)), "n_test": int(len(test_idx)), "selected": selected, "metrics": metrics, "timing": timing}, time.perf_counter() - started


def cohort_summary(folds: list[dict[str, Any]], donors: tuple[str, ...]) -> dict[str, Any]:
    model_names = ["candidate", "logistic", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor", "knn15", "svm", "gradient_boosting", "point_size_only"]
    per_donor: dict[str, dict[str, Any]] = {name: {} for name in model_names}
    cell_type_records: list[dict[str, Any]] = []
    within_cell_type_per_donor: dict[str, Any] = {}
    for fold in folds:
        donor = fold["test_donor"]
        for model_name in model_names:
            per_donor[model_name][donor] = fold["metrics"].get(model_name, {})
        cell_type_records.append({"donor": donor, **fold["metrics"].get("cell_type_only", {})})
        within_cell_type_per_donor[donor] = fold["metrics"].get("within_cell_type", {})
    summary: dict[str, Any] = {
        "donors": list(donors),
        "per_donor": per_donor,
        "cell_type_only_per_donor": cell_type_records,
        "within_cell_type_per_donor": within_cell_type_per_donor,
        "models": {},
    }
    for model_name in model_names:
        rows = [per_donor[model_name].get(donor, {}) for donor in donors]
        summary["models"][model_name] = {
            field: summarize_donor_metrics(rows, field)
            for field in ("class_balanced_aurc", "balanced_accuracy", "brier", "mapped_subset_balanced_accuracy")
        }
    celltype_ba = [r.get("balanced_accuracy") for r in cell_type_records if r.get("balanced_accuracy") is not None]
    summary["cell_type_only_mean_balanced_accuracy"] = float(np.mean(celltype_ba)) if celltype_ba else None
    expression_names = [name for name in model_names if name != "point_size_only"]
    celltype_by_donor = {row["donor"]: row.get("balanced_accuracy") for row in cell_type_records}
    common_donors = [
        donor for donor in donors
        if celltype_by_donor.get(donor) is not None
        and all(per_donor[name][donor].get("mapped_subset_balanced_accuracy") is not None for name in expression_names)
    ]
    expression_ba = [
        (name, float(np.mean([per_donor[name][donor]["mapped_subset_balanced_accuracy"] for donor in common_donors])))
        for name in expression_names
    ]
    if expression_ba and common_donors:
        best_name, best_ba = max(expression_ba, key=lambda item: item[1])
        mean_celltype_ba = float(np.mean([celltype_by_donor[donor] for donor in common_donors]))
        summary["composition_reading"] = {
            "best_expression_model_on_uniquely_mapped_cells": best_name,
            "best_expression_mean_balanced_accuracy": best_ba,
            "cell_type_only_mean_balanced_accuracy": mean_celltype_ba,
            "common_donors_scored": common_donors,
            "within_0_05_comparable": bool(abs(best_ba - mean_celltype_ba) <= 0.05),
        }
    else:
        summary["composition_reading"] = {"status": "not possible: no sufficient uniquely mapped two-class donor subset"}
    return summary


def run_cohort(data: Dataset, donors: tuple[str, ...], global_tie: np.ndarray, score_test: bool = True) -> tuple[dict[str, Any], float]:
    if not set(donors).issubset(set(data.donors.tolist())):
        raise StopRun(f"Cohort has missing donors: {sorted(set(donors) - set(data.donors.tolist()))}")
    for donor in donors:
        for label, label_name in ((1, "core"), (0, "periphery")):
            n = int(np.sum((data.donors == donor) & (data.labels == label)))
            if n == 0:
                raise StopRun(f"Donor {donor} has no {label_name} cells after QC")
    fold_records: list[dict[str, Any]] = []
    elapsed_total = 0.0
    for index, donor in enumerate(donors, start=1):
        print(f"[{index}/{len(donors)}] {donor}: nested donor fold starting", flush=True)
        record, elapsed = run_outer_fold(data, donors, donor, global_tie, score_test=score_test)
        elapsed_total += elapsed
        if score_test:
            fold_records.append(record)
        print(f"[{index}/{len(donors)}] {donor}: {elapsed:.1f}s", flush=True)
    if not score_test:
        return {"pilot_outer_folds": len(donors), "pilot_seconds": elapsed_total}, elapsed_total
    summary = cohort_summary(fold_records, donors)
    summary["folds"] = fold_records
    summary["runtime_seconds"] = elapsed_total
    if len(donors) == 6:
        candidate = summary["models"]["candidate"]["class_balanced_aurc"]
        logistic = summary["models"]["logistic"]["class_balanced_aurc"]
        candidate_mean = candidate["mean"]
        logistic_mean = logistic["mean"]
        per_donor_lower = sum(
            summary["per_donor"]["candidate"][donor]["class_balanced_aurc"]
            < summary["per_donor"]["logistic"][donor]["class_balanced_aurc"]
            for donor in donors
        )
        summary["frozen_descriptive_reading"] = {
            "candidate_mean_aurc_lower_than_logistic": bool(candidate_mean < logistic_mean),
            "candidate_lower_in_at_least_5_of_6_donors": bool(per_donor_lower >= 5),
            "donors_with_lower_candidate_aurc": int(per_donor_lower),
            "candidate_better_descriptive": bool(candidate_mean < logistic_mean and per_donor_lower >= 5),
            "no_inferential_superiority_claim": True,
        }
    return summary, elapsed_total


def build_celltype_strata(data: Dataset, summary: dict[str, Any], donors: tuple[str, ...]) -> dict[str, Any]:
    output: dict[str, Any] = {}
    per_donor_scores = summary.get("within_cell_type_per_donor", {})
    for cell_type in sorted({str(v) for v in data.cell_types if v is not None and not pd.isna(v)}):
        eligible_donors = []
        counts: dict[str, dict[str, int]] = {}
        for donor in donors:
            mask = (data.donors == donor) & (data.cell_types.astype(str) == cell_type)
            n_core = int(np.sum(mask & (data.labels == 1)))
            n_periphery = int(np.sum(mask & (data.labels == 0)))
            counts[donor] = {"core": n_core, "periphery": n_periphery}
            if n_core >= 10 and n_periphery >= 10:
                eligible_donors.append(donor)
        if len(eligible_donors) >= 4:
            per_donor = {
                donor: per_donor_scores.get(donor, {}).get(cell_type, {})
                for donor in eligible_donors
            }
            model_names = ("candidate", "logistic", "ablation_no_gate", "ablation_graph_and_blend", "ablation_mean_donor", "knn15", "svm", "gradient_boosting")
            donor_macro: dict[str, Any] = {}
            for model_name in model_names:
                model_rows = [per_donor[donor].get(model_name, {}) for donor in eligible_donors]
                donor_macro[model_name] = {
                    field: summarize_donor_metrics(model_rows, field)
                    for field in ("balanced_accuracy", "class_balanced_aurc", "brier")
                }
            output[cell_type] = {
                "status": "scored",
                "eligible_donors": eligible_donors,
                "counts": counts,
                "per_donor_metrics": per_donor,
                "donor_macro_metrics": donor_macro,
            }
    return output or {"status": "not possible: no cell type has >=10 cells of both classes in >=4 donors"}


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def safe_commit_hash() -> str | None:
    return current_commit()


def run_pilot(
    dataset: Dataset,
    strict: tuple[str, ...],
    sensitivity: tuple[str, ...],
    hashes: dict[str, str],
    data_load_seconds: float,
) -> int:
    if (OUTPUT_DIR / "pilot_timing.json").exists():
        raise StopRun("Pilot timing file already exists; refusing to overwrite")
    counts = dataset.flow["qc_by_donor_class"]
    sensitivity_ok = all(
        counts.get(donor, {}).get("core", 0) >= MIN_CELLS_PER_CLASS
        and counts.get(donor, {}).get("periphery", 0) >= MIN_CELLS_PER_CLASS
        for donor in sensitivity
    )
    if sensitivity_ok:
        pilot_donors = sensitivity
        cohort_name = "sensitivity_8_donor"
        outer_fold_count = 14
    else:
        pilot_donors = strict
        cohort_name = "primary_6_donor_sensitivity_omitted"
        outer_fold_count = 6
    pilot_data = dataset.subset(pilot_donors)
    pilot_tie = tie_ranks(pilot_data.cells)
    donor = pilot_donors[0]
    print(f"Timing-only pilot: {cohort_name} outer fold, held-out donor {donor}; no held-out predictions or metrics", flush=True)
    record, fold_seconds = run_outer_fold(pilot_data, pilot_donors, donor, pilot_tie, score_test=False)
    projected = fold_seconds * outer_fold_count + data_load_seconds
    payload = {
        "status": "timing_only",
        "pilot_cohort": cohort_name,
        "projected_outer_fold_count": outer_fold_count,
        "held_out_donor": donor,
        "held_out_predictions_or_metrics_saved": False,
        "inner_selection_metrics_saved": False,
        "pilot_fold_seconds": fold_seconds,
        "input_load_seconds": data_load_seconds,
        "projected_total_seconds": projected,
        "projected_hours": projected / 3600.0,
        "six_hour_limit_seconds": 21_600,
        "passes_runtime_gate": bool(projected <= 21_600),
        "fit_timing": record,
        "input_hashes": hashes,
        "script_sha256": hashes.get(str(Path(__file__).resolve().relative_to(PROJECT_ROOT))),
        "code_commit_before_run": safe_commit_hash(),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "limitations": ["Projection assumes one timing fold represents the selected cohort; cell counts and fit times vary by donor."],
    }
    atomic_json(OUTPUT_DIR / "pilot_timing.json", payload)
    print(f"Projected full run: {payload['projected_hours']:.2f} hours", flush=True)
    return 0 if payload["passes_runtime_gate"] else 3


def run_full(
    dataset: Dataset,
    strict: tuple[str, ...],
    sensitivity: tuple[str, ...],
    global_tie: np.ndarray,
    flow: dict[str, Any],
    hashes: dict[str, str],
    data_load_seconds: float,
) -> None:
    output_path = OUTPUT_DIR / "result.json"
    if output_path.exists():
        raise StopRun(f"Refusing to overwrite existing result: {output_path}")
    if not (OUTPUT_DIR / "pilot_timing.json").is_file():
        raise StopRun("Run and review timing pilot before full scoring")
    pilot = json.loads((OUTPUT_DIR / "pilot_timing.json").read_text(encoding="utf-8"))
    if not pilot.get("passes_runtime_gate"):
        raise StopRun("Timing pilot projects more than six hours; write and commit an amendment before changing the grid")
    runner_key = str(Path(__file__).resolve().relative_to(PROJECT_ROOT))
    if pilot.get("script_sha256") != hashes.get(runner_key) or pilot.get("input_hashes") != hashes:
        raise StopRun("Timing pilot was made with different sealed inputs or runner; repeat the pilot")
    if pilot.get("code_commit_before_run") != safe_commit_hash():
        raise StopRun("Timing pilot was made from a different Git commit; repeat the pilot")
    primary_data = dataset.subset(strict)
    primary_summary, primary_sec = run_cohort(primary_data, strict, global_tie[np.isin(dataset.donors, strict)], score_test=True)
    sensitivity_counts = dataset.flow["qc_by_donor_class"]
    sensitivity_ok = all(
        sensitivity_counts.get(donor, {}).get("core", 0) >= MIN_CELLS_PER_CLASS
        and sensitivity_counts.get(donor, {}).get("periphery", 0) >= MIN_CELLS_PER_CLASS
        for donor in sensitivity
    )
    if sensitivity_ok:
        sensitivity_data = dataset.subset(sensitivity)
        sensitivity_tie = tie_ranks(sensitivity_data.cells)
        sensitivity_summary, sensitivity_sec = run_cohort(sensitivity_data, sensitivity, sensitivity_tie, score_test=True)
    else:
        sensitivity_summary = {"status": "omitted: at least one lenient donor has fewer than 20 QC-passed cells in a class", "donor_class_counts": {d: sensitivity_counts.get(d, {}) for d in sensitivity}}
        sensitivity_sec = 0.0
    type_strata = build_celltype_strata(primary_data, primary_summary, strict)
    result = {
        "schema_version": 1,
        "status": "completed_descriptive_nested_donor_cv",
        "question": "Does the candidate improve donor-held-out Core-versus-Periphery prediction over ordinary L2 logistic regression?",
        "primary_cohort": primary_summary,
        "sensitivity_cohort": sensitivity_summary,
        "cell_type_stratum_eligibility": type_strata,
        "cell_flow_and_qc": flow,
        "implementation_choices": {
            "final_candidate_epoch": "median best epoch for the selected setting across inner donor folds",
            "candidate_early_stopping": {
                "monitor": "inner validation cross-entropy of unblended classifier probabilities",
                "maximum_epochs": MAX_EPOCHS,
                "patience_epochs": PATIENCE,
                "equal_loss_tie": "earliest epoch",
            },
            "hyperparameter_tie_tolerance": INNER_IMPROVEMENT_TOL,
            "hyperparameter_tie_rule": "first setting in ascending grid order within tolerance of the minimum inner mean AURC",
            "svm_probability_calibration": {
                "method": "donor-held-out margins followed by logistic sigmoid calibration",
                "calibrator": "LogisticRegression(C=1, solver='lbfgs', max_iter=1000, class_weight=None)",
                "calibration_donor_weight": "equal total weight per donor",
                "c_grid": list(SVM_C_GRID),
            },
            "gradient_boosting": {
                "learning_rate_grid": list(GBT_RATE_GRID),
                "n_estimators": 100,
                "max_depth": 2,
                "min_samples_leaf": 10,
                "subsample": 1.0,
                "loss": "log_loss",
                "internal_early_stopping": False,
                "seed": SEED,
            },
            "aurc_coverage": {
                "accepted_cell_counts": "every integer k from ceil(N/2) through N",
                "integration": "trapezoid over actual k/N coverage values",
                "missing_class_risk": 1.0,
            },
        },
        "protocol": {
            "path": "paper/trackA_yu2020/PROTOCOL.md",
            "amendment": "A2",
            "protocol_sha256": hashes[str((PROTOCOL_DIR / "PROTOCOL.md").relative_to(PROJECT_ROOT))],
            "donor_table_sha256": hashes[str((PROTOCOL_DIR / "donors.tsv").relative_to(PROJECT_ROOT))],
            "appendix_sha256": hashes[str((PROJECT_ROOT / "paper" / "research_plan_complete_project_appendix.md").relative_to(PROJECT_ROOT))],
            "seed": SEED,
            "class_mapping": {"Core_Tumoral": 1, "Periphery_Peritumoral": 0},
        },
        "source_hashes": hashes,
        "code_sha256": hashes.get(str(Path(__file__).resolve().relative_to(PROJECT_ROOT))),
        "code_commit_before_run": safe_commit_hash(),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
            "pytorch": torch.__version__,
        },
        "runtime_seconds": {
            "input_load": data_load_seconds,
            "primary": primary_sec,
            "sensitivity": sensitivity_sec,
            "total": data_load_seconds + primary_sec + sensitivity_sec,
        },
        "limits": [
            "Six primary donors support descriptive estimates only; no inferential superiority claim is made.",
            "Each surgical point is one library and region is nested within point inside donor.",
            "Core/periphery labels are author-supplied MRI-marked point labels, not histology labels for each cell.",
            "GBmap annotation_level_2 is expression-derived and is used only in the composition diagnostic.",
            "No cell-level result establishes a cell-intrinsic biological effect or spatial transfer.",
        ],
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(output_path, result)
    print(f"Wrote {output_path}", flush=True)


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description="Frozen donor-held-out Track A Yu 2020 analysis")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--pilot", action="store_true", help="Run timing-only pilot; save no held-out predictions or scores")
    mode.add_argument("--run", action="store_true", help="Run primary and eligible sensitivity analyses after pilot passes")
    args = parser.parse_args()
    try:
        if args.pilot and (OUTPUT_DIR / "pilot_timing.json").exists():
            raise StopRun("Pilot timing file already exists; refusing to overwrite")
        if args.run:
            if (OUTPUT_DIR / "result.json").exists():
                raise StopRun(f"Refusing to overwrite existing result: {OUTPUT_DIR / 'result.json'}")
            pilot_path = OUTPUT_DIR / "pilot_timing.json"
            if not pilot_path.is_file():
                raise StopRun("Run and review timing pilot before full scoring")
            pilot_preflight = json.loads(pilot_path.read_text(encoding="utf-8"))
            if not pilot_preflight.get("passes_runtime_gate"):
                raise StopRun("Timing pilot projects more than six hours; write and commit an amendment before changing the grid")
        sealed_hashes = verify_hashes()
        if args.run:
            runner_key = str(Path(__file__).resolve().relative_to(PROJECT_ROOT))
            if pilot_preflight.get("script_sha256") != sealed_hashes.get(runner_key) or pilot_preflight.get("input_hashes") != sealed_hashes:
                raise StopRun("Timing pilot was made with different sealed inputs or runner; repeat the pilot")
            if pilot_preflight.get("code_commit_before_run") != safe_commit_hash():
                raise StopRun("Timing pilot was made from a different Git commit; repeat the pilot")
        n_threads = min(8, os.cpu_count() or 1)
        torch.set_num_threads(n_threads)
        torch.set_num_interop_threads(1)
        torch.use_deterministic_algorithms(True)
        load_started = time.perf_counter()
        dataset, strict, sensitivity, flow, hashes = load_dataset()
        data_load_seconds = time.perf_counter() - load_started
        global_tie = tie_ranks(dataset.cells)
        if args.pilot:
            return run_pilot(dataset, strict, sensitivity, hashes, data_load_seconds)
        run_full(dataset, strict, sensitivity, global_tie, flow, hashes, data_load_seconds)
        return 0
    except StopRun as exc:
        print(f"STOP: {exc}", file=sys.stderr, flush=True)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
