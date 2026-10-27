#!/usr/bin/env python3
"""Script 151: Track A donor table (plan appendix A1). Data preparation only: no model is fit or scored.

Reads ONLY obs metadata (h5py) and small slices of X from output/01_filtered_three_class.h5ad, resolves
sample (ID1) -> donor with explicit alias rules, and writes
  data/manifests/trackA_donors.csv
  output/trackA_donor_table/summary.json

Alias rules (resolve_donors):
  R1  A sample ID ending in a known suffix (accutase, citeseq) is a candidate alias of the ID without the suffix.
  R2  Samples whose orig.ident donor part (text after the "tissue:" prefix) is shared are candidate aliases.
  R3  A candidate group is merged only if diagnosis, age, sex and MGMT agree for every member. Otherwise it is
      NOT merged and the disagreement is recorded as a conflict.
  Metadata-only collisions (same age, sex, MGMT, diagnosis but no R1/R2 link) are NOT merged in `donor_id`; they are
  listed as candidates and merged only in the extra column `donor_id_conservative` (plan: "when unsure, merge").
"""
from __future__ import annotations

import hashlib
import json
import re
import sys
from pathlib import Path as _Path

import h5py
import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
H5AD = PROJECT_ROOT / "output" / "01_filtered_three_class.h5ad"
MANIFEST = PROJECT_ROOT / "data" / "manifests" / "trackA_donors.csv"
SUMMARY = PROJECT_ROOT / "output" / "trackA_donor_table" / "summary.json"
SCRIPT = "151_trackA_donor_table"
SUFFIXES = ("accutase", "citeseq")
META_FIELDS = ("diagnosis", "age", "sex", "mgmt")
MISSING_META_TOKENS = {"", "na", "n/a", "nan", "none", "null", "<na>", "unknown", "not known",
                       "not available", "not recorded", "not reported", "missing"}
ASSAY_BY_SUFFIX = {"citeseq": "CITE-seq (inferred from ID suffix)"}
ASSAY_DEFAULT = "not recorded in file (no assay column)"
CSV_COLUMNS = ["sample_id", "donor_id", "alias_evidence", "diagnosis", "is_gbm", "source_fresh_frozen", "assay",
               "n_core", "n_periphery", "n_healthy", "age", "sex", "mgmt", "orig_ident", "donor_id_conservative",
               "metadata_only_candidates"]


# ----------------------------------------------------------------------------------------------- alias logic
def strip_suffix(sample_id: str) -> tuple[str, str | None]:
    for s in SUFFIXES:
        m = re.fullmatch(rf"(.+?){s}", sample_id, flags=re.IGNORECASE)
        if m:
            return m.group(1), s
    return sample_id, None


def donor_part(orig_ident: str) -> str:
    return orig_ident.split(":", 1)[-1]


def _is_missing_meta(value) -> bool:
    """Return true for nulls and common string encodings of missing metadata."""
    try:
        is_null = pd.isna(value)
        if isinstance(is_null, (bool, np.bool_)) and is_null:
            return True
    except (TypeError, ValueError):
        pass
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return str(value).strip().casefold() in MISSING_META_TOKENS


def _meta_issues(rows: list[pd.Series]) -> tuple[dict[str, list[str]], list[str]]:
    """Return conflicting and missing fields; missing metadata cannot confirm identity."""
    conflicting, missing = {}, []
    for field in META_FIELDS:
        values = [row[field] for row in rows]
        known = {str(value).strip() for value in values if not _is_missing_meta(value)}
        if any(_is_missing_meta(value) for value in values):
            missing.append(field)
        if len(known) > 1:
            conflicting[field] = sorted(known)
    return conflicting, missing


def _meta_str(r: pd.Series) -> str:
    return f"age={r['age']}, sex={r['sex']}, mgmt={r['mgmt']}, diagnosis={r['diagnosis']}"


def resolve_donors(samples: pd.DataFrame):
    """samples: one row per sample_id with diagnosis, age, sex, mgmt, orig_idents (iterable of strings).

    Returns (donor_id Series, evidence Series, merges list, conflicts list, links list). The returned donor_id
    is built through src.donor_cv.donor_map (refuses blank IDs, chained aliases, one sample -> two donors).
    """
    s = samples.set_index("sample_id", drop=False)
    ids = list(s.index)
    # candidate links: (frozenset of sample ids) -> reasons
    reasons: dict[frozenset, list[str]] = {}
    base_of = {i: strip_suffix(i) for i in ids}
    by_base: dict[str, list[str]] = {}
    for i in ids:
        by_base.setdefault(base_of[i][0], []).append(i)
    for base, members in by_base.items():
        if len(members) > 1:
            sufs = sorted({base_of[m][1] for m in members if base_of[m][1]})
            reasons.setdefault(frozenset(members), []).append(
                f"ID suffix rule: {sorted(members)} share base '{base}' (suffix {sufs})")
    by_part: dict[str, list[str]] = {}
    for i in ids:
        for oi in s.loc[i, "orig_idents"]:
            by_part.setdefault(donor_part(oi), []).append(i)
    for part, members in by_part.items():
        members = sorted(set(members))
        if len(members) > 1:
            reasons.setdefault(frozenset(members), []).append(
                f"shared orig.ident donor part '{part}' in {members}")
    # union-find over candidate links (so overlapping reasons form one group)
    parent = {i: i for i in ids}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for grp in reasons:
        g = sorted(grp)
        for other in g[1:]:
            parent[find(other)] = find(g[0])
    groups: dict[str, list[str]] = {}
    for i in ids:
        groups.setdefault(find(i), []).append(i)

    alias_to_canon: dict[str, str] = {}
    evidence = {i: "none (singleton, no alias link found)" for i in ids}
    merges, conflicts, links = [], [], []
    for members in groups.values():
        if len(members) == 1:
            continue
        members = sorted(members, key=lambda m: (len(m), m))
        rows = [s.loc[m] for m in members]
        why = [r for k, rs in reasons.items() if k <= set(members) for r in rs]
        links.append({"members": members, "reasons": why})
        conf, missing = _meta_issues(rows)
        if conf or missing:
            conflicts.append({"members": members, "reasons": why, "conflicting_fields": conf,
                              "missing_fields": missing})
            issue_text = []
            if conf:
                issue_text.append(f"conflicting metadata {conf}")
            if missing:
                issue_text.append(f"missing metadata {missing}")
            for m in members:
                evidence[m] = f"NOT merged: candidate alias of {members} but " + "; ".join(issue_text)
            continue
        canon = members[0]
        ev = "; ".join(why) + f"; metadata agree ({_meta_str(rows[0])})"
        for m in members:
            evidence[m] = ev if m != canon else ev + " [canonical ID]"
            if m != canon:
                alias_to_canon[m] = canon
        merges.append({"donor_id": canon, "members": members, "evidence": ev})

    from src.donor_cv import donor_map  # reuse validated mapping (refuses chained aliases / blank IDs)
    donor = donor_map(pd.DataFrame({"sample_id": ids, "donor_id": ids}), aliases=alias_to_canon)
    return donor.reindex(ids), pd.Series(evidence).reindex(ids), merges, conflicts, links


def metadata_only_candidates(samples: pd.DataFrame, donor: pd.Series):
    """Post-merge donors sharing age, sex, MGMT and diagnosis but not linked by an alias rule.
    Missing metadata is not identity evidence. Returns (candidate groups list, conservative donor Series).
    Groups with a diagnosis conflict are listed in `diag_conflicts` and not merged."""
    d = samples.assign(donor_id=donor.reindex(samples["sample_id"]).to_numpy())
    one = d.drop_duplicates("donor_id")
    complete = one[list(META_FIELDS)].apply(lambda col: col.map(lambda value: not _is_missing_meta(value)))
    one = one.loc[complete.all(axis=1)]
    groups, diag_conflicts, alias = [], [], {}
    for key, g in one.groupby(["age", "sex", "mgmt"], observed=True):
        if g["donor_id"].nunique() < 2:
            continue
        ids = sorted(g["donor_id"])
        if g["diagnosis"].nunique() > 1:
            diag_conflicts.append({"donors": ids, "age_sex_mgmt": list(map(str, key)),
                                   "diagnoses": sorted(g["diagnosis"].astype(str).unique())})
            continue
        canon = sorted(ids, key=lambda m: (len(m), m))[0]
        groups.append({"donors": ids, "age_sex_mgmt": list(map(str, key)), "diagnosis": str(g["diagnosis"].iloc[0]),
                       "source": sorted(g["source"].astype(str).unique()), "canonical_if_merged": canon})
        for i in ids:
            if i != canon:
                alias[i] = canon
    cons = donor.map(lambda x: alias.get(x, x))
    return groups, diag_conflicts, cons


# ----------------------------------------------------------------------------------------------- h5ad reading
def read_obs_col(obs: h5py.Group, name: str) -> np.ndarray:
    node = obs[name]
    if isinstance(node, h5py.Group) and "categories" in node:
        cats = np.array([c.decode() if isinstance(c, bytes) else c for c in node["categories"][:]], dtype=object)
        codes = node["codes"][:]
        out = np.empty(len(codes), dtype=object)
        ok = codes >= 0
        out[ok] = cats[codes[ok]]
        out[~ok] = None
        return out
    vals = node["values"][:] if isinstance(node, h5py.Group) else node[:]
    return np.array([v.decode() if isinstance(v, bytes) else v for v in vals], dtype=object)


def sha256_file(path: _Path, chunk: int = 1 << 24) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def raw_count_check(f: h5py.File, ncount: np.ndarray, n_head: int = 2000) -> dict:
    """Small-slice check that X holds integer raw counts. Never loads the full matrix."""
    X = f["X"]
    res = {"X_encoding": X.attrs.get("encoding-type"), "X_shape": [int(v) for v in X.attrs["shape"]],
           "X_dtype": str(X["data"].dtype), "layers_keys": list(f["layers"].keys()), "has_raw_group": "raw" in f,
           "uns_keys": list(f["uns"].keys())}
    nnz = X["data"].shape[0]
    head = X["data"][:n_head].astype(np.float64)
    res["first_values_checked"] = int(len(head))
    res["first_values_all_integer"] = bool(np.all(head == np.round(head)))
    res["first_values_min"] = float(head.min())
    res["first_values_max"] = float(head.max())
    # extra slices (middle and end) so the check is not only the first cells
    extra = []
    for name, start in (("middle", nnz // 2), ("end", max(nnz - n_head, 0))):
        v = X["data"][start:start + n_head].astype(np.float64)
        extra.append({"slice": name, "n": int(len(v)), "all_integer": bool(np.all(v == np.round(v))),
                      "max": float(v.max())})
    res["extra_slices"] = extra
    # five seeded random 100k-value slices: minimum, quantiles, share of integer-valued entries
    rng = np.random.default_rng(0)
    rand = []
    for st in rng.integers(0, max(nnz - 100000, 1), 5):
        v = X["data"][int(st):int(st) + 100000].astype(np.float64)
        rand.append({"start": int(st), "n": int(len(v)), "min": float(v.min()), "max": float(v.max()),
                     "p1_p50_p99": [float(q) for q in np.percentile(v, [1, 50, 99])],
                     "fraction_integer_valued": float(np.mean(v == np.round(v)))})
    res["random_slices"] = rand
    # row sums of first 200 cells vs obs nCount_RNA
    ip = X["indptr"][:201]
    d = X["data"][: int(ip[-1])].astype(np.float64)
    rs = np.array([d[ip[i]:ip[i + 1]].sum() for i in range(200)])
    nc = ncount[:200].astype(np.float64)
    res["first_200_cells_rowsum_equals_nCount_RNA"] = bool(np.allclose(rs, nc, rtol=0, atol=0.5))
    res["first_200_cells_rowsum_min_max"] = [float(rs.min()), float(rs.max())]
    res["first_200_cells_nCount_RNA_min_max"] = [float(nc.min()), float(nc.max())]
    ok = (res["first_values_all_integer"] and all(e["all_integer"] for e in extra)
          and res["first_200_cells_rowsum_equals_nCount_RNA"]
          and all(r["fraction_integer_valued"] == 1.0 for r in rand) and not res["layers_keys"] and not res["has_raw_group"])
    res["verdict"] = ("X appears to be integer raw counts (sampled slices only; no layers/raw alternative present)"
                      if ok else "NOT integer raw counts: X holds non-integer, already-transformed values (plan A1.1: stop and document); "
                                                "raw counts not found in X, layers or raw")
    res["raw_counts_matrix"] = "X" if ok else None
    return res


# ----------------------------------------------------------------------------------------------- main
def main() -> int:
    sys.path.insert(0, str(PROJECT_ROOT))
    with h5py.File(H5AD, "r") as f:
        obs = f["obs"]
        c = {n: read_obs_col(obs, n) for n in ["ID1", "orig.ident", "class", "Diagnosis", "source", "Gender", "Age",
                                                "MGMT_methylation_status", "nCount_RNA"]}
        has_assay = [n for n in obs.keys() if "assay" in n.lower() or "chemistry" in n.lower() or "platform" in n.lower()]
        raw = raw_count_check(f, c["nCount_RNA"].astype(float))
    df = pd.DataFrame({"sample_id": c["ID1"].astype(str), "orig_ident": c["orig.ident"].astype(str),
                       "cls": c["class"].astype(str), "diagnosis": c["Diagnosis"].astype(str),
                       "source": c["source"].astype(str), "sex": c["Gender"].astype(str),
                       "age": c["Age"].astype(str), "mgmt": c["MGMT_methylation_status"].astype(str)})
    n_cells = len(df)

    rows = []
    for sid, g in df.groupby("sample_id", sort=True):
        row = {"sample_id": sid, "orig_idents": sorted(g["orig_ident"].unique())}
        for f_ in ("diagnosis", "source", "sex", "age", "mgmt"):
            vals = sorted(g[f_].unique())
            if len(vals) != 1:
                raise ValueError(f"{sid}: field {f_} not constant within sample: {vals}")
            row[f_] = vals[0]
        for k, name in (("n_core", "Core"), ("n_periphery", "Periphery"), ("n_healthy", "Healthy")):
            row[k] = int((g["cls"] == name).sum())
        rows.append(row)
    samples = pd.DataFrame(rows)

    donor, evidence, merges, conflicts, links = resolve_donors(samples)
    cand_groups, cand_diag_conf, cons = metadata_only_candidates(samples, donor)
    cand_for = {}
    for g in cand_groups:
        for d_ in g["donors"]:
            cand_for[d_] = "; ".join(x for x in g["donors"] if x != d_)

    out = samples.copy()
    out["donor_id"] = donor.to_numpy()
    out["alias_evidence"] = evidence.to_numpy()
    out["is_gbm"] = out["diagnosis"].str.fullmatch("Glioblastoma")
    out["source_fresh_frozen"] = out["source"]
    out["assay"] = [ASSAY_BY_SUFFIX.get((strip_suffix(s)[1] or "").lower(), ASSAY_DEFAULT) for s in out["sample_id"]]
    out["orig_ident"] = out["orig_idents"].map(";".join)
    out["donor_id_conservative"] = cons.to_numpy()
    out["metadata_only_candidates"] = out["donor_id"].map(lambda d_: cand_for.get(d_, ""))
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    out[CSV_COLUMNS].to_csv(MANIFEST, index=False)

    def donor_tables(col: str):
        agg = out.groupby(col).agg(n_core=("n_core", "sum"), n_periphery=("n_periphery", "sum"),
                                   n_healthy=("n_healthy", "sum"), is_gbm=("is_gbm", "all"),
                                   members=("sample_id", lambda s: sorted(s)))
        gbm = agg[agg["is_gbm"]]
        two = gbm[(gbm["n_core"] > 0) & (gbm["n_periphery"] > 0)]
        core_only = gbm[(gbm["n_core"] > 0) & (gbm["n_periphery"] == 0)]
        rec = lambda t: {k: {"n_core": int(r.n_core), "n_periphery": int(r.n_periphery), "samples": r.members}
                         for k, r in t.iterrows()}
        return {"n_donors": int(len(agg)), "gbm_donors": sorted(gbm.index), "n_gbm_donors": int(len(gbm)),
                "two_class_gbm_donors": rec(two), "n_two_class_gbm_donors": int(len(two)),
                "core_only_gbm_donors": rec(core_only), "n_core_only_gbm_donors": int(len(core_only)),
                "non_gbm_donors": {k: {"diagnosis": out.loc[out[col] == k, "diagnosis"].iloc[0],
                                       "n_healthy": int(r.n_healthy)} for k, r in agg[~agg["is_gbm"]].iterrows()}}

    primary, conservative = donor_tables("donor_id"), donor_tables("donor_id_conservative")
    summary = {
        "script": SCRIPT, "input_file": str(H5AD.relative_to(PROJECT_ROOT)).replace("\\", "/"),
        "input_sha256": sha256_file(H5AD), "input_size_bytes": H5AD.stat().st_size,
        "n_cells": int(n_cells), "n_samples": int(len(samples)),
        "class_cell_counts": {k: int(v) for k, v in df["cls"].value_counts().items()},
        "n_donors_after_merge": primary["n_donors"],
        "n_donors_after_merge_conservative": conservative["n_donors"],
        "primary_rule_donors": primary,
        "conservative_rule_donors": conservative,
        "merges_applied": merges, "alias_links_found": links, "conflicts": conflicts,
        "metadata_only_candidates_not_merged_in_donor_id": cand_groups,
        "metadata_only_diagnosis_conflicts": cand_diag_conf,
        "assay_column_present": has_assay,
        "assay_note": "No assay column in obs. 'CITE-seq' is inferred only from the 'citeseq' ID suffix; other samples are unrecorded.",
        "raw_count_check": raw,
        "caveats": ["Raw-count check reads slices only (first/middle/end values and first 200 cells), not the full matrix.",
                    "Donor identity rests on ID suffixes, shared orig.ident and age/sex/MGMT metadata; no genotype check was done.",
                    "No classifier was fit or scored."],
    }
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
    print(f"wrote {MANIFEST} ({len(out)} rows) and {SUMMARY}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
