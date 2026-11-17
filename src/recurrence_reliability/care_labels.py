"""Blinded, column-restricted extraction of the CARE author cell-type table (gates C6/C7).

Authorised by the user on 2026-11-16. Reads data/recurrence_reliability/raw/care/celltype_meta_data_2025_01_08.RDS
(GBM-CARE-WT commit tree 383b0a3; Zenodo 10.5281/zenodo.14966015, CC BY 4.0).

Prints and stores ONLY: column names and storage types (schema); the four kept columns CellID, Sample,
ID, CellType; and per-specimen cell / Malignant counts for the 84 sample-level eligible specimens.
Every other column is dropped without being printed, summarised, or written. No expression file is read.

    py -3 src/recurrence_reliability/care_labels.py
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from recurrence_reliability.rds_reader import column, column_names, read_rds  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW = PROJECT_ROOT / "data" / "recurrence_reliability" / "raw" / "care"
SRC = RAW / "celltype_meta_data_2025_01_08.RDS"
OUT = PROJECT_ROOT / "output" / "recurrence_reliability_care_metadata"
KEEP = ["CellID", "Sample", "ID", "CellType"]
FLAG_PATTERN = re.compile(r"(os|surv|status|outcome|mgmt|interval|state|score|program|mp|age|sex|gender|death)", re.I)
TYPE_NAMES = {10: "logical", 13: "integer", 14: "double", 16: "character", 19: "list"}
MIN_MALIGNANT = 50


def main() -> None:
    df = read_rds(SRC)
    names = column_names(df)
    schema = []
    for n, col in zip(names, df.value):
        kind = "factor" if "levels" in col.attrs else TYPE_NAMES.get(col.rtype, str(col.rtype))
        schema.append({"column": n, "storage": kind, "name_flagged_outcome_like": bool(FLAG_PATTERN.search(n)) and n not in KEEP})
    print("class:", list(df.attrs["class"].value), "| n columns:", len(names))
    for s in schema:
        print(f"  {s['column']:<30} {s['storage']:<10} {'FLAGGED (not read)' if s['name_flagged_outcome_like'] else ''}")
    missing = [k for k in KEEP if k not in names]
    if missing:
        raise SystemExit(f"required columns missing: {missing}; stopping without reading anything else")
    kept = {k: column(df, k) for k in KEEP}
    del df   # other columns are never decoded further, printed, or written

    n = len(kept["CellID"])
    manifest = json.loads((OUT / "care_eligibility_manifest.json").read_text(encoding="utf-8"))
    st1 = {s["Sample ID"]: s for s in manifest["specimens"]}
    eligible = {x for p in manifest["pairs"] if p["status"] == "sample_level_eligible" for x in (p["primary"], p["first_recurrence"])}

    issues = []
    ids_unique = len(set(kept["CellID"])) == n
    bad_form = 0
    sample_to_id = defaultdict(set)
    barcode_re = re.compile(r"^(?P<s>[^-]+)-(?P<bc>[ACGT]{16})(-\d+)?$")
    for cid, smp, sid in zip(kept["CellID"], kept["Sample"], kept["ID"]):
        m = barcode_re.match(cid or "")
        if not m or m.group("s") != smp:
            bad_form += 1
        sample_to_id[smp].add(sid)
    if bad_form:
        issues.append(f"{bad_form} CellIDs not of the form <Sample>-<16-base barcode>[-n] with matching Sample")
    multi = {k: sorted(v) for k, v in sample_to_id.items() if len(v) != 1}
    if multi:
        issues.append(f"Samples mapping to several IDs: {multi}")
    # `Sample` is the processing lab's internal library code (e.g. TKU3074, COLUMBIA..., SMMG001A), not the
    # Table S1 Tumor ID; the specimen join key is `ID` (P#T#), which must be unique per Sample and present in
    # Table S1 / GEO. The Sample-code prefix is recorded as provenance only.
    id_to_samples = defaultdict(set)
    for smp, idset in sample_to_id.items():
        for sid in idset:
            id_to_samples[sid].add(smp)
    not_in_st1 = sorted(i for i in id_to_samples if i not in st1)
    several = {k: sorted(v) for k, v in id_to_samples.items() if len(v) != 1}
    if not_in_st1:
        issues.append(f"IDs not in Table S1/GEO: {not_in_st1}")
    if several:
        issues.append(f"IDs with several Sample codes: {several}")
    prefix_by_lab = Counter((re.match(r"[A-Za-z]+", next(iter(id_to_samples[i]))).group(0), st1[i]["Tissue source"], st1[i]["Lab processed"])
                            for i in id_to_samples if i in st1)

    per_id = defaultdict(Counter)
    for sid, ct in zip(kept["ID"], kept["CellType"]):
        if sid in eligible:
            per_id[sid][ct] += 1
    specimens = {sid: {"n_cells_author_table": sum(per_id[sid].values()), "n_malignant_author_table": per_id[sid].get("Malignant", 0)}
                 for sid in sorted(eligible, key=lambda x: (int(x[1:].split("T")[0]), x))}
    pairs = []
    for p in manifest["pairs"]:
        if p["status"] != "sample_level_eligible":
            continue
        a, b = specimens[p["primary"]], specimens[p["first_recurrence"]]
        low = [s for s in (p["primary"], p["first_recurrence"]) if specimens[s]["n_malignant_author_table"] < MIN_MALIGNANT]
        pairs.append({"patient": p["patient"], "primary": p["primary"], "first_recurrence": p["first_recurrence"],
                      "malignant_primary": a["n_malignant_author_table"], "malignant_first_recurrence": b["n_malignant_author_table"],
                      "passes_50_before_wang_filters": not low, "below_50": low})
    summary = {
        "scope": "blinded column-restricted read of the CARE author cell-type table; schema + 4 columns only; no expression read",
        "authorised": "user, 2026-11-16",
        "date": "2026-11-16",
        "source": {"file": "data/celltype_meta_data_2025_01_08.RDS", "repo_commit_tree": "383b0a320f0993a4529167cdb591dac63e103674",
                   "url": "https://raw.githubusercontent.com/dravishays/GBM-CARE-WT/383b0a320f0993a4529167cdb591dac63e103674/data/celltype_meta_data_2025_01_08.RDS",
                   "sha256": hashlib.sha256(SRC.read_bytes()).hexdigest(), "bytes": SRC.stat().st_size,
                   "license": "Zenodo 10.5281/zenodo.14966015 CC BY 4.0; GitHub MIT"},
        "schema": schema,
        "columns_kept": KEEP,
        "n_rows": n,
        "cellid_unique": ids_unique,
        "celltype_labels": sorted({c for c in kept["CellType"] if c is not None}),
        "n_samples_in_table": len(sample_to_id),
        "sample_code_prefix_by_tissue_source_and_lab": [{"prefix": a, "tissue_source": b, "lab": c, "n_specimens": n}
                                                        for (a, b, c), n in sorted(prefix_by_lab.items())],
        "eligible_specimens_present": sum(1 for s in specimens.values() if s["n_cells_author_table"] > 0),
        "issues": issues,
        "specimens": specimens,
        "pairs": pairs,
        "pairs_passing_50_before_wang_filters": sum(p["passes_50_before_wang_filters"] for p in pairs),
        "note": "Counts are after author QC (200-10,000 genes, MT <= 3%) and confidence filtering; the Wang filters (>= 200 genes, MT <= 2.5%) need the count matrices and are applied only after the external protocol is locked.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "care_label_extract_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    with gzip.open(RAW / "care_celltype_4col.tsv.gz", "wt", encoding="utf-8") as f:   # local only (git-ignored)
        f.write("\t".join(KEEP) + "\n")
        for row in zip(*(kept[k] for k in KEEP)):
            f.write("\t".join("" if v is None else str(v) for v in row) + "\n")
    print(json.dumps({k: summary[k] for k in ("n_rows", "cellid_unique", "celltype_labels", "n_samples_in_table",
                                              "eligible_specimens_present", "pairs_passing_50_before_wang_filters", "issues")}))


if __name__ == "__main__":
    main()
