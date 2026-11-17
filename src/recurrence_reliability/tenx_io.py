"""Reading 10x per-library files (matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz) and gene modules."""
from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.io
import scipy.sparse as sp

PROJECT_ROOT = Path(__file__).resolve().parents[2]
GSM_DIR = PROJECT_ROOT / "data" / "recurrence_reliability" / "raw" / "gsm"
NEFTEL = PROJECT_ROOT / "data" / "external" / "neftel2019" / "IDHwt.GBM.MetaModules.tsv"
VERHAAK = PROJECT_ROOT / "data" / "external" / "neftel2019" / "TCGA.Verhaak.GBM.txt"


@dataclass
class Library:
    library_id: str
    counts: sp.csc_matrix          # genes x cells, integer counts
    gene_ids: np.ndarray           # Ensembl IDs
    symbols: np.ndarray
    feature_types: np.ndarray
    barcodes: np.ndarray           # GEM suffix removed
    header: str


def _lines(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return [line.rstrip("\n") for line in f if line.strip()]


def read_header(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        return f.readline().strip()


def load_library(library_id: str, gsm_dir: Path = GSM_DIR) -> Library:
    from .metadata import norm_barcode
    feats = [l.split("\t") for l in _lines(gsm_dir / f"{library_id}_features.tsv.gz")]
    widths = {len(f) for f in feats}
    if widths == {1}:
        # GSE174554 deposits symbol-only feature files (verified 2026-11-16): no Ensembl IDs or
        # feature-type column. The symbol is the only gene identifier; row order is the matrix row order.
        feats = [[f[0], f[0], "Gene Expression (assumed; symbol-only file)"] for f in feats]
    elif widths != {3}:
        raise ValueError(f"{library_id}: features file has {widths} columns")
    bcs = np.array([norm_barcode(b) for b in _lines(gsm_dir / f"{library_id}_barcodes.tsv.gz")])
    mpath = gsm_dir / f"{library_id}_matrix.mtx.gz"
    header = read_header(mpath)
    with gzip.open(mpath, "rb") as f:
        m = scipy.io.mmread(f)
    m = sp.coo_matrix(m)
    return Library(library_id, m.tocsc(), np.array([f[0] for f in feats]), np.array([f[1] for f in feats]),
                   np.array([f[2] for f in feats]), bcs, header)


def load_modules(path: Path, columns: list[str]) -> dict[str, list[str]]:
    """Gene lists from a tab file with one module per column (blank cells ignored)."""
    rows = [l.rstrip("\n").split("\t") for l in path.read_text(encoding="utf-8").splitlines() if l.strip()]
    head = rows[0]
    out = {}
    for c in columns:
        j = head.index(c)
        out[c] = [r[j].strip() for r in rows[1:] if j < len(r) and r[j].strip()]
    return out


def neftel_mes() -> dict[str, list[str]]:
    return load_modules(NEFTEL, ["MESlike1", "MESlike2"])


def verhaak_mes() -> dict[str, list[str]]:
    return load_modules(VERHAAK, ["Mesenchymal"])
