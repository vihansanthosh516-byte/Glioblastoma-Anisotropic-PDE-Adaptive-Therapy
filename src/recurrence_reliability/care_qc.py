"""Run the locked CARE count QC. This module computes no program scores."""
import csv
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import rdata

from .qc import MAX_MITO_FRAC, MIN_GENES

ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data/recurrence_reliability/raw/care"
META = ROOT / "output/recurrence_reliability_care_metadata"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def specimen_qc(path, labels):
    with gzip.open(path, "rb") as stream:
        obj = rdata.parser.parse_data(stream.read(), extension=".rds").object
    attrs = {}
    node = obj.attributes
    while node is not None and node.info.type.name == "LIST":
        tag = node.tag
        while tag.info.type.name == "REF":
            tag = tag.referenced_object
        key = tag.value.value.decode()
        attrs[key], node = node.value
    if obj.info.type.name not in {"REAL", "INT"} or not {"dim", "dimnames"} <= attrs.keys():
        raise ValueError(f"{path.name}: unsupported count matrix representation")

    def decode(item):
        while item.info.type.name == "REF":
            item = item.referenced_object
        return item.value.decode()

    genes, cells = ([decode(item) for item in axis.value] for axis in attrs["dimnames"].value)
    counts = obj.value.reshape(tuple(attrs["dim"].value), order="F")
    if counts.shape != (len(genes), len(cells)):
        raise ValueError(f"{path.name}: matrix dimensions differ from dimnames")
    normalize = lambda barcode: re.sub(r"-1$", "", barcode)
    indexes = {normalize(cell): i for i, cell in enumerate(cells)}
    keys = {normalize(cell): label for cell, label in labels.items()}
    if len(indexes) != len(cells) or len(keys) != len(labels):
        raise ValueError(f"{path.name}: ambiguous GEM suffix normalization")
    if any(re.fullmatch(r"[^-]+-[ACGT]{16}(?:-1)?", cell) is None for cell in (*cells, *labels)):
        raise ValueError(f"{path.name}: unexpected cell ID format")
    matched = indexes.keys() & keys.keys()
    mito = [i for i, gene in enumerate(genes) if gene.upper().startswith("MT-")]
    if not mito:
        raise ValueError(f"{path.name}: no mitochondrial gene symbols")
    good = 0
    for key in matched:
        if keys[key] != "Malignant":
            continue
        col = counts[:, indexes[key]]
        if not np.all(np.isfinite(col) & (col >= 0) & (col == np.floor(col))):
            raise ValueError(f"{path.name}: invalid count values")
        total = col.sum()
        if np.count_nonzero(col) >= MIN_GENES and col[mito].sum() / total <= MAX_MITO_FRAC:
            good += 1
    return dict(labels=len(labels), matched=len(matched), matrix_only=len(cells)-len(matched),
                malignant_after_qc=good, passes=len(matched)*100 >= len(labels)*95 and good >= 50,
                sha256=digest(path))


def main():
    summary = json.loads((META / "care_label_extract_summary.json").read_text(encoding="utf-8"))
    manifest = json.loads((META / "care_eligibility_manifest.json").read_text(encoding="utf-8"))
    pairs = [p for p in summary["pairs"] if p["passes_50_before_wang_filters"]]
    ids = {p[k] for p in pairs for k in ("primary", "first_recurrence")}
    rows = {s["Sample ID"]: s for s in manifest["specimens"] if s["Sample ID"] in ids}
    if len(pairs) != 33 or len(ids) != 66 or rows.keys() != ids:
        raise ValueError("CARE manifest does not match the locked cohort size")
    labels = {sid: {} for sid in ids}
    label_path = RAW / "care_celltype_4col.tsv.gz"
    with gzip.open(label_path, "rt", encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["ID"] in labels:
                table = labels[row["ID"]]
                if row["CellID"] in table:
                    raise ValueError("Duplicate CARE author cell ID")
                table[row["CellID"]] = row["CellType"]
    results, missing = {}, []
    for sid, row in rows.items():
        path = RAW / "counts" / row["geo_file"]
        if not path.is_file() or path.stat().st_size != row["geo_file_bytes"]:
            missing.append(sid)
            continue
        if len(labels[sid]) != summary["specimens"][sid]["n_cells_author_table"]:
            raise ValueError(f"{sid}: author label count differs from blinded audit")
        results[sid] = specimen_qc(path, labels[sid])
        print(sid, results[sid]["matched"], results[sid]["malignant_after_qc"], flush=True)
    report = {"specimens": results, "missing": missing, "label_sha256": digest(label_path)}
    if not missing:
        report["eligible_pairs"] = [p["patient"] for p in pairs
                                    if all(results[p[k]]["passes"] for k in ("primary", "first_recurrence"))]
    name = "care_qc_partial.json" if missing else "care_qc_report.json"
    (META / name).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
