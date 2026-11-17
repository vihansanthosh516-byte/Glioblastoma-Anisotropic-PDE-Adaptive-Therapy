"""Phase C technical QC (paper/recurrence_reliability_qc_plan.md). Computes no program scores.

    py -3 src/recurrence_reliability/qc.py

Writes output/recurrence_reliability/qc_report.json.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    __package__ = "recurrence_reliability"

from recurrence_reliability import tenx_io as rio  # noqa: E402
from recurrence_reliability.metadata import (MIN_MALIGNANT_NUCLEI, OUT, PROJECT_ROOT,  # noqa: E402
                                             load_annotation, sha256)

MIN_GENES = 200          # Wang et al. Methods
MAX_MITO_FRAC = 0.025    # Wang et al. Methods
MIN_MODULE_MAPPED = 0.80


def check_library(lib: rio.Library) -> dict:
    c = lib.counts.tocoo()
    n_genes, n_cells = lib.counts.shape
    data = c.data
    pairs = c.row.astype(np.int64) * n_cells + c.col
    out = {
        "header": lib.header,
        "shape_genes_cells": [int(n_genes), int(n_cells)],
        "dims_match_files": bool(n_genes == len(lib.gene_ids) and n_cells == len(lib.barcodes)),
        "nnz": int(c.nnz),
        "all_values_positive_integers": bool(np.all(data > 0) and np.all(np.mod(data, 1) == 0)),
        "duplicate_entries": int(len(pairs) - len(np.unique(pairs))),
        "feature_types": sorted(set(lib.feature_types.tolist())),
        "unique_gene_ids": bool(len(set(lib.gene_ids)) == len(lib.gene_ids)),
        "duplicated_symbols": int(len(lib.symbols) - len(set(lib.symbols))),
        "unique_barcodes": bool(len(set(lib.barcodes)) == len(lib.barcodes)),
        "gene_id_hash": hashlib.sha256("\n".join(lib.gene_ids).encode()).hexdigest(),
    }
    return out


def nucleus_metrics(lib: rio.Library) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    m = lib.counts
    umi = np.asarray(m.sum(axis=0)).ravel()
    genes = np.diff(m.indptr)
    mito = np.array([s.upper().startswith("MT-") for s in lib.symbols])
    mt = np.asarray(m[mito].sum(axis=0)).ravel()
    return umi, genes, np.divide(mt, umi, out=np.zeros_like(mt, dtype=float), where=umi > 0)


def module_mapping(symbols: np.ndarray, modules: dict[str, list[str]]) -> dict:
    counts: dict[str, int] = {}
    for s in symbols:
        counts[s] = counts.get(s, 0) + 1
    res = {}
    for name, genes in modules.items():
        mapped = [g for g in genes if counts.get(g, 0) == 1]
        res[name] = {"n_genes": len(genes), "n_mapped_unique": len(mapped),
                     "fraction": round(len(mapped) / len(genes), 4),
                     "dropped": [g for g in genes if counts.get(g, 0) != 1],
                     "passes": len(mapped) / len(genes) >= MIN_MODULE_MAPPED}
    return res


def main() -> None:
    manifest = json.loads((OUT / "metadata_manifest.json").read_text(encoding="utf-8"))
    ann = load_annotation()
    eligible = {sid for p in manifest["pairs"] if p["eligible"] for sid in p["specimens"].values()}
    report = {"phase": "C technical QC (no program scores computed)", "date": "2026-11-16",
              "filters": {"MIN_GENES": MIN_GENES, "MAX_MITO_FRAC": MAX_MITO_FRAC,
                          "MIN_MALIGNANT_NUCLEI": MIN_MALIGNANT_NUCLEI},
              "libraries": {}, "specimens": {}, "file_sha256": {}}
    gene_hashes, symbols_ref = set(), None
    for s in manifest["specimens"]:
        if s["specimen"] not in eligible:
            continue
        for l in s["libraries"]:
            indep = any(r["library_id"] == l["library_id"] and r["relationship_to_primary_library"] == "independent_capture"
                        for r in s["library_relationships"])
            if not (l["is_primary_library"] or indep):
                continue
            lid = l["library_id"]
            for suf in ("matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"):
                p = rio.GSM_DIR / f"{lid}_{suf}"
                report["file_sha256"][p.name] = sha256(p)
            lib = rio.load_library(lid)
            chk = check_library(lib)
            gene_hashes.add(chk["gene_id_hash"])
            symbols_ref = lib.symbols if symbols_ref is None else symbols_ref
            umi, genes, mito = nucleus_metrics(lib)
            labels = ann.get(l["annotation_key"], {}) if l["annotation_key"] else {}
            lab = np.array([labels.get(b, "NA") for b in lib.barcodes])
            keep = (genes >= MIN_GENES) & (mito <= MAX_MITO_FRAC)
            tum = lab == "Tumor"
            chk.update({
                "role": "primary" if l["is_primary_library"] else "independent_capture_secondary",
                "annotation_key": l["annotation_key"],
                "join_rate": round(float(np.mean(lab != "NA")), 4),
                "removed_by_min_genes": int(np.sum(genes < MIN_GENES)),
                "removed_by_mito": int(np.sum(mito > MAX_MITO_FRAC)),
                "n_tumor_after_filters": int(np.sum(keep & tum)),
                "n_all_after_filters": int(np.sum(keep)),
                "median_umi_tumor": float(np.median(umi[keep & tum])) if np.any(keep & tum) else None,
                "median_genes_tumor": float(np.median(genes[keep & tum])) if np.any(keep & tum) else None,
                "median_umi_all": float(np.median(umi[keep])) if np.any(keep) else None,
            })
            report["libraries"][lid] = chk
    fails = []
    for s in manifest["specimens"]:
        if s["specimen"] not in eligible:
            continue
        p = next(l for l in s["libraries"] if l["is_primary_library"])
        q = report["libraries"][p["library_id"]]
        problems = []
        # Symbol-only feature files: duplicated symbols are allowed here and handled by the module
        # mapping rule (a gene with several matching rows is dropped); feature type cannot be checked.
        if not (q["dims_match_files"] and q["all_values_positive_integers"] and q["duplicate_entries"] == 0
                and q["unique_barcodes"] and len(q["feature_types"]) == 1):
            problems.append("matrix/feature/barcode integrity check failed")
        if q["n_tumor_after_filters"] < MIN_MALIGNANT_NUCLEI:
            problems.append(f"{q['n_tumor_after_filters']} tumor nuclei after filters < {MIN_MALIGNANT_NUCLEI}")
        report["specimens"][s["specimen"]] = {"primary_library": p["library_id"], "problems": problems}
        if problems:
            fails.append(s["specimen"])
    report["feature_lists_identical"] = len(gene_hashes) == 1
    mods = {**rio.neftel_mes(), **{"Verhaak_Mesenchymal": rio.verhaak_mes()["Mesenchymal"]}}
    report["module_mapping"] = module_mapping(symbols_ref, mods)
    report["module_sources"] = {"neftel": {"path": "data/external/neftel2019/IDHwt.GBM.MetaModules.tsv", "sha256": sha256(rio.NEFTEL)},
                                "verhaak": {"path": "data/external/neftel2019/TCGA.Verhaak.GBM.txt", "sha256": sha256(rio.VERHAAK)}}
    final = []
    for p in manifest["pairs"]:
        if not p["eligible"]:
            continue
        bad = [sid for sid in p["specimens"].values() if sid in fails]
        if not bad:
            final.append({"table1_pair": p["table1_pair"], **p["specimens"]})
    report["qc_failed_specimens"] = fails
    report["final_cohort"] = final
    report["final_cohort_n_patients"] = len(final)
    report["study_stop"] = not all(v["passes"] for k, v in report["module_mapping"].items() if k.startswith("MES"))
    (OUT / "qc_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"n_patients": len(final), "qc_failed": fails, "features_identical": report["feature_lists_identical"],
                      "module_mapping": {k: v["fraction"] for k, v in report["module_mapping"].items()},
                      "study_stop": report["study_stop"]}))


if __name__ == "__main__":
    main()
