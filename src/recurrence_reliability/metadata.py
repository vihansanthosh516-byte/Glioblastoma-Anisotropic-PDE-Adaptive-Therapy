"""Phase B: candidate-patient table and eligibility manifest from metadata only.

Reads the existing metadata JSON artifacts, the authors' tumor/normal annotation file, and the
per-library barcode lists. Never reads expression matrices, scores, or survival.

    py -3 src/recurrence_reliability/metadata.py

Writes output/recurrence_reliability/metadata_manifest.json and candidate_table.tsv.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
META = PROJECT_ROOT / "output" / "paper_strengthening_metadata"
RAW = PROJECT_ROOT / "data" / "recurrence_reliability" / "raw"
OUT = PROJECT_ROOT / "output" / "recurrence_reliability"

ANNOTATION = RAW / "GSE174554_Tumor_normal_metadata.txt.gz"
ANNOTATION_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE174nnn/GSE174554/suppl/GSE174554_Tumor_normal_metadata.txt.gz"
GSM_URL = "https://ftp.ncbi.nlm.nih.gov/geo/samples/{stub}nnn/{gsm}/suppl/{name}"

# Fixed before any malignant count was computed (2026-11-16). Technical floor only: the
# per-specimen mean score and its cell-bootstrap SE need enough nuclei to be non-degenerate.
MIN_MALIGNANT_NUCLEI = 50
MIN_JOIN_RATE = 0.95          # annotated barcodes / library barcodes
SAME_CAPTURE_MIN_OVERLAP = 0.90   # of the smaller library's barcodes
INDEPENDENT_MAX_OVERLAP = 0.01


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def norm_barcode(b: str) -> str:
    """Strip the 10x GEM-well suffix ('-1' in barcodes.tsv, '.1' in the annotation file)."""
    return re.sub(r"[-.]\d+$", "", b.strip())


def load_annotation(path: Path = ANNOTATION) -> dict[str, dict[str, str]]:
    """{annotation sample key: {barcode: 'Tumor'|'Normal'}}. Identical duplicate rows collapse;
    conflicting duplicates raise."""
    ann: dict[str, dict[str, str]] = defaultdict(dict)
    with gzip.open(path, "rt") as f:
        header = next(f).split()
        if header != ["Sample#_Barcode", "Tumor_Normal_annotation"]:
            raise ValueError(f"unexpected annotation header {header}")
        for line in f:
            key, label = line.split()
            sample, bc = key.split("_", 1)
            if label not in ("Tumor", "Normal"):
                raise ValueError(f"unexpected label {label!r}")
            b = norm_barcode(bc)
            if ann[sample].get(b, label) != label:
                raise ValueError(f"conflicting labels for {key}")
            ann[sample][b] = label
    return dict(ann)


def read_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt") as f:
        return [line.strip() for line in f if line.strip()]


def match_annotation(barcodes: set[str], ann: dict[str, dict[str, str]], specimen: str) -> tuple[str | None, str]:
    """Return the annotation key whose barcode set best covers this library, among the two
    source-named candidates (specimen, specimen+'v2'). The match must be unique."""
    scored = []
    for key in (specimen, specimen + "v2"):
        if key in ann:
            a = set(ann[key])
            scored.append((len(a & barcodes), len(a - barcodes), key))
    if not scored:
        return None, "no annotation sample for specimen"
    good = [s for s in scored if s[1] == 0 and s[0] >= MIN_JOIN_RATE * len(barcodes)]
    if len(good) == 1:
        return good[0][2], "unique annotation key; annotation barcodes are a subset of library barcodes"
    if len(good) > 1:
        return None, "ambiguous: more than one annotation key covers the library"
    return None, "no annotation key covers the library at the required join rate"


def library_relationship(a: set[str], b: set[str]) -> tuple[str, float]:
    frac = len(a & b) / max(1, min(len(a), len(b)))
    if frac >= SAME_CAPTURE_MIN_OVERLAP:
        return "same_capture_barcode_subset", frac
    if frac <= INDEPENDENT_MAX_OVERLAP:
        return "independent_capture", frac
    return "ambiguous_partial_overlap", frac


def build() -> dict:
    cov_path = META / "recurrence_processed_coverage.json"
    cw_path = META / "recurrence_table1_crosswalk.json"
    cov = json.loads(cov_path.read_text(encoding="utf-8"))
    cw = json.loads(cw_path.read_text(encoding="utf-8"))
    t1 = {r["ID"]: r for r in cw["snRNA_samples"]}
    ann = load_annotation()

    specimens, pairs = [], defaultdict(list)
    for s in cov["specimens"]:
        sid, gsm = s["specimen"], s["gsm"]
        row = t1[sid]
        libs = []
        for lib in s["libraries"]:
            name = lib["files"]["barcodes.tsv.gz"]["name"]
            path = RAW / "gsm" / name
            raw = read_barcodes(path)
            bset = {norm_barcode(b) for b in raw}
            key, why = match_annotation(bset, ann, sid)
            labels = ann.get(key, {}) if key else {}
            n_t = sum(1 for b in bset if labels.get(b) == "Tumor")
            n_n = sum(1 for b in bset if labels.get(b) == "Normal")
            libs.append({
                "library_id": lib["prefix"],
                "is_primary_library": lib["prefix"] == f"{gsm}_{sid}",
                "barcodes_file": name,
                "barcodes_url": GSM_URL.format(stub=gsm[:-3], gsm=gsm, name=name),
                "barcodes_sha256": sha256(path),
                "n_barcodes": len(raw),
                "n_unique_barcodes": len(bset),
                "annotation_key": key,
                "annotation_match_note": why,
                "n_tumor": n_t, "n_normal": n_n,
                "n_unannotated": len(bset) - n_t - n_n,
                "join_rate": round((n_t + n_n) / max(1, len(bset)), 4),
                "_bset": bset,
            })
        rel = []
        prim = [l for l in libs if l["is_primary_library"]]
        if len(prim) != 1:
            raise ValueError(f"{sid}: expected exactly one primary library")
        for l in libs:
            if not l["is_primary_library"]:
                kind, frac = library_relationship(prim[0]["_bset"], l["_bset"])
                rel.append({"library_id": l["library_id"], "relationship_to_primary_library": kind,
                            "overlap_fraction_of_smaller": round(frac, 4)})
        for l in libs:
            l.pop("_bset")
        p = prim[0]
        reasons = []
        if p["annotation_key"] is None:
            reasons.append(f"primary library has no unambiguous author annotation ({p['annotation_match_note']})")
        elif p["join_rate"] < MIN_JOIN_RATE:
            reasons.append(f"annotation join rate {p['join_rate']} < {MIN_JOIN_RATE}")
        elif p["n_tumor"] < MIN_MALIGNANT_NUCLEI:
            reasons.append(f"{p['n_tumor']} author-labelled tumor nuclei < {MIN_MALIGNANT_NUCLEI}")
        spec = {
            "specimen": sid, "gsm": gsm, "pair_id_geo_coverage": s["pair_id"],
            "table1_pair": row["Pair#"], "stage": row["Stage"], "geo_stage": s["stage"],
            "diagnosis": row["Diagnosis"], "idh": row["IDH"],
            "libraries": libs, "library_relationships": rel,
            "specimen_exclusion_reasons": reasons,
        }
        specimens.append(spec)
        pairs[row["Pair#"]].append(spec)

    pair_rows = []
    for pid, members in sorted(pairs.items()):
        stages = sorted(m["stage"] for m in members)
        reasons = []
        if stages != ["Primary", "Recurrent"]:
            reasons.append(f"pair does not have one primary and one recurrent specimen: {stages}")
        for m in members:
            reasons += [f"{m['specimen']}: {r}" for r in m["specimen_exclusion_reasons"]]
        pair_rows.append({"table1_pair": pid, "specimens": {m["stage"]: m["specimen"] for m in members},
                          "eligible": not reasons, "exclusion_reasons": reasons})
    return {"specimens": specimens, "pairs": pair_rows, "cov_path": cov_path, "cw_path": cw_path}


def main() -> None:
    res = build()
    OUT.mkdir(parents=True, exist_ok=True)
    eligible = [p for p in res["pairs"] if p["eligible"]]
    manifest = {
        "phase": "B metadata/annotation/batch eligibility (no expression values read)",
        "date": "2026-11-16",
        "rules": {"MIN_MALIGNANT_NUCLEI": MIN_MALIGNANT_NUCLEI, "MIN_JOIN_RATE": MIN_JOIN_RATE,
                  "SAME_CAPTURE_MIN_OVERLAP": SAME_CAPTURE_MIN_OVERLAP,
                  "INDEPENDENT_MAX_OVERLAP": INDEPENDENT_MAX_OVERLAP,
                  "primary_library_rule": "the library whose file prefix is <GSM>_<specimen> (no batch suffix)",
                  "annotation_rule": "author annotation key (specimen or specimen+'v2') whose barcodes are a subset of the library's barcodes and cover >= MIN_JOIN_RATE of them; must be unique"},
        "sources": {
            "annotation": {"url": ANNOTATION_URL, "sha256": sha256(ANNOTATION), "bytes": ANNOTATION.stat().st_size,
                           "provenance": "Wang et al. 2022 tumor/normal labels; Methods: CONICSmat v1.0 CNV calls (LRT < 0.05, BIC difference > 50, posterior > 0.5)"},
            "coverage_manifest": {"path": str(res["cov_path"].relative_to(PROJECT_ROOT)).replace("\\", "/"), "sha256": sha256(res["cov_path"])},
            "table1_crosswalk": {"path": str(res["cw_path"].relative_to(PROJECT_ROOT)).replace("\\", "/"), "sha256": sha256(res["cw_path"]),
                                 "workbook_sha256": "58e42239d4e105cfbbda6441bd640afbd64e0c55f87e07b1e0aaab0d1eff1c42"},
            "geo_series": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE174554",
        },
        "counts": {"candidate_pairs": len(res["pairs"]), "eligible_pairs": len(eligible),
                   "eligible_specimens": 2 * len(eligible)},
        "pairs": res["pairs"],
        "specimens": res["specimens"],
    }
    (OUT / "metadata_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    cols = ["table1_pair", "stage", "specimen", "gsm", "diagnosis", "idh", "primary_library", "annotation_key",
            "n_barcodes", "n_tumor", "n_normal", "join_rate", "other_libraries", "pair_eligible", "exclusion_reasons"]
    lines = ["\t".join(cols)]
    pair_by = {p["table1_pair"]: p for p in res["pairs"]}
    for s in sorted(res["specimens"], key=lambda x: (x["table1_pair"], x["stage"])):
        p = next(l for l in s["libraries"] if l["is_primary_library"])
        others = ";".join(f"{r['library_id']}={r['relationship_to_primary_library']}({r['overlap_fraction_of_smaller']})"
                          for r in s["library_relationships"])
        pr = pair_by[s["table1_pair"]]
        lines.append("\t".join(str(x) for x in [s["table1_pair"], s["stage"], s["specimen"], s["gsm"], s["diagnosis"], s["idh"],
                     p["library_id"], p["annotation_key"], p["n_barcodes"], p["n_tumor"], p["n_normal"], p["join_rate"],
                     others or "-", pr["eligible"], " | ".join(pr["exclusion_reasons"]) or "-"]))
    (OUT / "candidate_table.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(manifest["counts"]))


if __name__ == "__main__":
    main()
