"""CARE (GSE274546) sample-level eligibility gate. Metadata only.

Reads the official GEO brieftext records, the GEO supplementary file list, and SELECTED columns of
Spitzer et al. 2025 Supplementary Table 1 (identity, stage, assay flags, diagnosis, IDH, pre-surgery
radiation and alkylating agents). Survival, outcome, surgical interval, MGMT, age, sex, location and
steroid columns are never read. No expression file, cell annotation, score or state is opened.

    py -3 src/recurrence_reliability/care_metadata.py

Writes output/recurrence_reliability_care_metadata/care_eligibility_manifest.json.
"""
from __future__ import annotations

import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RAW = PROJECT_ROOT / "data" / "recurrence_reliability" / "raw" / "care"
OUT = PROJECT_ROOT / "output" / "recurrence_reliability_care_metadata"

ST1_COLUMNS = ["Sample ID", "Tumor ID", "Patient ID", "Time point", "Primary or recurrent", "Tissue source",
               "Lab processed", "10x", "Smart-seq2", "Diagnosis", "IDH mutation status",
               "Radiation before the surgery", "Alkylate agents before the surgery"]
STAGE_OF_INDEX = {"0": "Primary", "1": "1st Recurrent", "2": "2nd Recurrent"}
SOURCES = {
    "gse_gsm.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546&targ=gsm&form=text&view=brief",
    "gse_self.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274546&targ=self&form=text&view=brief",
    "gse548_self.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274548&targ=self&form=text&view=brief",
    "filelist.txt": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE274nnn/GSE274546/suppl/filelist.txt",
    "spitzer_st1.xlsx": "https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-025-02168-4/MediaObjects/41588_2025_2168_MOESM2_ESM.xlsx",
    "tree.json": "https://api.github.com/repos/dravishays/GBM-CARE-WT/git/trees/HEAD?recursive=1",
    "zen.json": "https://zenodo.org/api/records/14966015",
}


def sha256(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def clean(v) -> str:
    return re.sub(r"\s+", " ", str(v or "")).strip()


def parse_gsm(path: Path) -> dict[str, dict]:
    recs, cur = {}, None
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("^SAMPLE"):
            cur = {"gsm": line.split("=", 1)[1].strip()}
        elif cur is not None and line.startswith("!Sample_"):
            k, _, v = line[8:].partition(" = ")
            cur.setdefault(k, []).append(v)
            if k == "data_row_count":
                m = re.match(r"Patient (\d+), Timepoint (\d+)", cur["title"][0])
                _, code, stage = [x.strip() for x in cur["source_name_ch1"][0].split(",")]
                prot = [c.split(":", 1)[1].strip() for c in cur["characteristics_ch1"] if c.startswith("protocol number")]
                sid = f"P{m.group(1)}T{m.group(2)}"
                recs[sid] = {"gsm": cur["gsm"], "geo_title": cur["title"][0], "geo_code": code, "geo_stage": stage,
                             "geo_protocol": prot[0], "geo_file": cur["supplementary_file_1"][0].rsplit("/", 1)[1],
                             "organism": cur["organism_ch1"][0]}
    return recs


def parse_filelist(path: Path) -> dict[str, int]:
    out = {}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split("\t")
        if len(f) >= 4 and f[0] == "File":
            out[f[1]] = int(f[3])
    return out


def read_st1(path: Path) -> list[dict]:
    import openpyxl
    ws = openpyxl.load_workbook(path, read_only=True)["TableS1"]
    rows = ws.iter_rows(min_row=3, values_only=True)
    hdr = [clean(h) for h in next(rows)]
    idx = {c: hdr.index(c) for c in ST1_COLUMNS}      # only these columns are ever touched
    out = []
    for r in rows:
        rec = {c: clean(r[i]) for c, i in idx.items()}
        if any(rec.values()):
            out.append(rec)
    return out


def pair_rule(stages: dict[str, dict]) -> tuple[str, list[str]]:
    """Locked Wang question: untreated primary vs first recurrence after radiotherapy + TMZ only."""
    reasons = []
    if "Primary" not in stages or "1st Recurrent" not in stages:
        return "excluded", [f"no primary + 1st-recurrence pair (stages: {sorted(stages)})"]
    p, r = stages["Primary"], stages["1st Recurrent"]
    if p["Radiation before the surgery"] != "No" or p["Alkylate agents before the surgery"] != "No":
        reasons.append("primary specimen not treatment-naive")
    if r["Radiation before the surgery"] != "Yes" or r["Alkylate agents before the surgery"] != "Yes (TMZ)":
        reasons.append(f"first recurrence not after RT + TMZ only (RT={r['Radiation before the surgery']}, "
                       f"alkylating={r['Alkylate agents before the surgery']})")
    for s in (p, r):
        if s["IDH mutation status"] != "wildtype" or s["Diagnosis"] != "Glioblastoma, IDH-wildtype":
            reasons.append(f"{s['Sample ID']} not IDH-wildtype GBM")
    return ("sample_level_eligible" if not reasons else "excluded"), reasons


def build() -> dict:
    geo = parse_gsm(RAW / "gse_gsm.txt")
    sizes = parse_filelist(RAW / "filelist.txt")
    st1 = read_st1(RAW / "spitzer_st1.xlsx")
    crosswalk_issues, specimens, by_patient = [], [], defaultdict(dict)
    for s in st1:
        sid = s["Sample ID"]
        g = geo.get(sid)
        code, idx = s["Tumor ID"].split("_")
        if g is None:
            crosswalk_issues.append(f"{sid}: no GEO record")
            continue
        if code != g["geo_code"] or s["Primary or recurrent"] != g["geo_stage"] or STAGE_OF_INDEX[idx] != g["geo_stage"] \
                or f"_{sid}_" not in g["geo_file"] or s["Patient ID"] != sid.split("T")[0]:
            crosswalk_issues.append(f"{sid}: identity fields disagree between Table S1 and GEO")
        rec = {**s, **g, "geo_file_bytes": sizes.get(g["geo_file"])}
        specimens.append(rec)
        by_patient[s["Patient ID"]][s["Primary or recurrent"]] = rec
    if set(geo) - {s["Sample ID"] for s in st1}:
        crosswalk_issues.append(f"GEO records not in Table S1: {sorted(set(geo) - {s['Sample ID'] for s in st1})}")
    pairs = []
    for pid, stages in sorted(by_patient.items(), key=lambda kv: int(kv[0][1:])):
        status, reasons = pair_rule(stages)
        pairs.append({"patient": pid, "site_code": next(iter(stages.values()))["geo_code"],
                      "stages_present": sorted(stages), "status": status, "reasons": reasons,
                      "primary": stages.get("Primary", {}).get("Sample ID"),
                      "first_recurrence": stages.get("1st Recurrent", {}).get("Sample ID"),
                      "protocols": sorted({v["geo_protocol"] for v in stages.values()}),
                      "labs": sorted({v["Lab processed"] for v in stages.values()})})
    return {"specimens": specimens, "pairs": pairs, "crosswalk_issues": crosswalk_issues}


def main() -> None:
    res = build()
    elig = [p for p in res["pairs"] if p["status"] == "sample_level_eligible"]
    pr1 = [p for p in res["pairs"] if p["primary"] and p["first_recurrence"]]
    manifest = {
        "scope": "CARE GSE274546 sample-level metadata gate; no expression, cell annotation, state, score or survival read",
        "date": "2026-11-16",
        "sources": {k: {"url": u, "sha256": sha256(RAW / k), "accessed": "2026-11-16"} for k, u in SOURCES.items()},
        "st1_columns_read": ST1_COLUMNS,
        "st1_columns_not_read": ["Gender", "Age at initial diagnosis", "MGMT protomoter methylation status", "Tumor location",
                                 "Site of recurrence", "Steroid treatment before the surgery",
                                 "Number of cycles treated by adjuvant TMZ before the surgery", "Surgical interval (month)",
                                 "Overall survival (month)", "Outcome", "unnamed trailing columns"],
        "counts": {"specimens": len(res["specimens"]), "patients": len(res["pairs"]),
                   "patients_with_primary_and_first_recurrence": len(pr1),
                   "sample_level_eligible_pairs": len(elig),
                   "cell_level_confirmed_pairs": 0},
        "crosswalk_issues": res["crosswalk_issues"],
        "pairs": res["pairs"],
        "specimens": res["specimens"],
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "care_eligibility_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(manifest["counts"]), "crosswalk issues:", len(res["crosswalk_issues"]))


if __name__ == "__main__":
    main()
