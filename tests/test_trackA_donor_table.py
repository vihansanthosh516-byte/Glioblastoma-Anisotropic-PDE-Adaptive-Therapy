"""Alias merging in src/151_trackA_donor_table.py on synthetic sample metadata (no real data)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

_spec = importlib.util.spec_from_file_location("trackA_donor_table", REPO / "src" / "151_trackA_donor_table.py")
dt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(dt)


def S(sid, oi=None, dx="Glioblastoma", age="64", sex="Male", mgmt="methylated", source="Frozen"):
    return {"sample_id": sid, "orig_idents": oi or [f"core:{sid}"], "diagnosis": dx, "age": age, "sex": sex,
            "mgmt": mgmt, "source": source}


def run(rows):
    return dt.resolve_donors(pd.DataFrame(rows))


def test_suffix_merge_when_metadata_agree():
    donor, ev, merges, conflicts, _ = run([S("A1"), S("A1citeseq"), S("B2", age="50")])
    assert donor["A1citeseq"] == "A1" and donor["A1"] == "A1" and donor["B2"] == "B2"
    assert len(merges) == 1 and conflicts == []
    assert "suffix" in ev["A1citeseq"] and "metadata agree" in ev["A1citeseq"]


def test_shared_orig_ident_merge():
    donor, _, merges, conflicts, _ = run([S("X7", oi=["core:X7"]), S("X7accutase", oi=["core:X7"])])
    assert donor["X7accutase"] == "X7" and len(merges) == 1 and conflicts == []


def test_refuse_merge_on_conflict():
    donor, ev, merges, conflicts, _ = run([S("A1"), S("A1citeseq", age="65")])
    assert donor["A1"] == "A1" and donor["A1citeseq"] == "A1citeseq"
    assert merges == [] and len(conflicts) == 1
    assert "age" in conflicts[0]["conflicting_fields"]
    assert ev["A1citeseq"].startswith("NOT merged")


@pytest.mark.parametrize("field,val", [("sex", "Female"), ("mgmt", "unmethylated"), ("dx", "Meningioma")])
def test_refuse_on_each_field(field, val):
    kw = {"dx": val} if field == "dx" else {field: val}
    donor, _, merges, conflicts, _ = run([S("A1"), S("A1citeseq", **kw)])
    assert donor["A1"] != donor["A1citeseq"] and merges == [] and len(conflicts) == 1


def test_lone_suffix_sample_not_renamed():
    donor, ev, merges, _, _ = run([S("C3citeseq")])
    assert donor["C3citeseq"] == "C3citeseq" and merges == [] and ev["C3citeseq"].startswith("none")


def test_metadata_only_collision_not_merged_but_flagged():
    samples = pd.DataFrame([S("A1"), S("B2")])  # same age/sex/MGMT/dx, no alias link
    donor, *_ = dt.resolve_donors(samples)
    assert donor["A1"] != donor["B2"]
    groups, diag_conf, cons = dt.metadata_only_candidates(samples, donor)
    assert len(groups) == 1 and groups[0]["donors"] == ["A1", "B2"] and diag_conf == []
    assert cons["A1"] == cons["B2"]


def test_metadata_only_diagnosis_conflict_not_merged():
    samples = pd.DataFrame([S("A1"), S("B2", dx="Meningioma")])
    donor, *_ = dt.resolve_donors(samples)
    groups, diag_conf, cons = dt.metadata_only_candidates(samples, donor)
    assert groups == [] and len(diag_conf) == 1 and cons["A1"] != cons["B2"]
