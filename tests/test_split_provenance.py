"""Split integrity and forecast-cache provenance on synthetic inputs (no patient data read or written)."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import leakage_checks as lc  # noqa: E402


def _load(name, file):
    spec = importlib.util.spec_from_file_location(name, REPO / "src" / file)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ms = _load("manifest_94_t", "94_patient_manifest.py")
pm = _load("pde_100_t", "100_pde_manifest.py")

PIDS = [f"SYN-{i:03d}" for i in range(40)]


def rule_split():
    return pd.DataFrame({"patient_id": PIDS, "fold": [ms.fold_of(p) for p in PIDS]})


def test_rule_split_is_valid():
    s = rule_split()
    lc.assert_one_fold_per_patient(s, ms.N_FOLDS)
    lc.assert_split_matches_rule(s, ms.fold_of)


@pytest.mark.parametrize("fold", [float("nan"), 1.5, -1, 5, "x"])
def test_bad_fold_values_rejected(fold):
    s = rule_split()
    s["fold"] = s["fold"].astype(object)
    s.loc[0, "fold"] = fold
    with pytest.raises(AssertionError):
        lc.assert_one_fold_per_patient(s, 5)


def test_blank_patient_and_empty_fold_rejected():
    s = rule_split()
    s.loc[0, "patient_id"] = None
    with pytest.raises(AssertionError):
        lc.assert_one_fold_per_patient(s)
    s = rule_split()
    s = s[s["fold"] != 2]
    with pytest.raises(AssertionError, match="empty folds"):
        lc.assert_one_fold_per_patient(s, 5)


def test_edited_split_rejected():
    s = rule_split()
    s.loc[3, "fold"] = (s.loc[3, "fold"] + 1) % 5
    lc.assert_one_fold_per_patient(s, 5)          # still well formed ...
    with pytest.raises(AssertionError, match="differs from its rule"):
        lc.assert_split_matches_rule(s, ms.fold_of)   # ... but not the frozen rule


def test_pair_fold_column_must_match_split():
    s = rule_split()
    pairs = pd.DataFrame({"patient_id": PIDS[:3], "fold": s["fold"].head(3).to_numpy()})
    lc.assert_pairs_follow_split(pairs, s)
    pairs.loc[1, "fold"] = (pairs.loc[1, "fold"] + 1) % 5
    with pytest.raises(AssertionError, match="pair fold differs"):
        lc.assert_pairs_follow_split(pairs, s)


def test_unmapped_patient_rejected():
    folds = lc.fold_map(rule_split(), 5)
    lc.assert_all_have_fold(PIDS[:5], folds)
    with pytest.raises(AssertionError, match="no fold"):
        lc.assert_all_have_fold(PIDS[:5] + ["NEW-1"], folds)


def test_load_split_checks_rule(tmp_path):
    f = tmp_path / "split.csv"
    rule_split().to_csv(f, index=False)
    assert pm.load_split(f).to_dict() == dict(zip(PIDS, [ms.fold_of(p) for p in PIDS]))
    s = rule_split()
    s.loc[0, "fold"] = (s.loc[0, "fold"] + 1) % 5
    s.to_csv(f, index=False)
    with pytest.raises(AssertionError):
        pm.load_split(f)


# ---- forecast cache signature ----

GRID = {"arms": [("aniso", 1.0), ("iso_homog", 1.0)], "ds": [0.01, 0.03], "rhos": [0.1, 0.2], "selected": None}
PAIR = {"patient_id": "SYN-001", "tp_in": 1, "tp_out": 2, "day_in": 30.0, "dt_days": 60.0, "fold": 3,
        "treatment": {"radiation_end_day": 40}}


def sig(pair=PAIR, h=2.0, grid=GRID, solver="TensorFK", atlas="a1"):
    return pm.cache_signature(pair, h, grid, solver, atlas)


def test_signature_stable_and_sensitive():
    assert sig() == sig(dict(PAIR))
    changed = [sig(h=1.0), sig(solver="TensorFKMonotone"), sig(atlas="a2"),
               sig(grid={**GRID, "rhos": [0.1, 0.3]}), sig(grid={**GRID, "arms": [("aniso", 10.0)]}),
               sig(pair={**PAIR, "dt_days": 61.0}), sig(pair={**PAIR, "day_in": 31.0}),
               sig(pair={**PAIR, "fold": 2}), sig(pair={**PAIR, "treatment": None})]
    assert len({sig(), *changed}) == len(changed) + 1


def test_signature_selected_mode_tracks_this_fold_cell():
    sel = {"3": {"aniso|r=1": "aniso|r=1|d=0.01|rho=0.1"}, "2": {"aniso|r=1": "aniso|r=1|d=0.03|rho=0.2"}}
    g = {**GRID, "selected": sel}
    sel2 = {"3": {"aniso|r=1": "aniso|r=1|d=0.03|rho=0.1"}, "2": sel["2"]}
    sel3 = {"3": sel["3"], "2": {"aniso|r=1": "aniso|r=1|d=0.3|rho=0.2"}}
    assert sig(grid=g) != sig(grid={**GRID, "selected": sel2})   # this fold's cell changed -> recompute
    assert sig(grid=g) == sig(grid={**GRID, "selected": sel3})   # another fold's cell changed -> reuse


def test_worker_reuses_only_matching_cache(tmp_path, monkeypatch):
    calls = []

    def fake_forecast(pair, atlas, h, grid):
        calls.append(pair["dt_days"])
        return {"patient_id": pair["patient_id"], "tp_in": pair["tp_in"], "tp_out": pair["tp_out"],
                "dt_days": pair["dt_days"], "dice": {"no_change": 0.5, "grid": {}}}

    monkeypatch.setattr(pm, "forecast_pair", fake_forecast)
    monkeypatch.setattr(pm.ria, "_ATLAS", None, raising=False)
    monkeypatch.setattr(pm.ria, "_ATLAS_SHA", "a1", raising=False)
    monkeypatch.setattr(pm, "_CACHE", tmp_path)
    monkeypatch.setattr(pm, "_H", 2.0)
    monkeypatch.setattr(pm, "_GRID", GRID)
    r1 = pm.worker(dict(PAIR))
    assert r1["cache_sig"] == sig(solver=pm.ria.TensorFK.__name__) and r1["fold"] == 3
    pm.worker(dict(PAIR))
    assert len(calls) == 1                                   # identical inputs: cache hit
    monkeypatch.setattr(pm, "_GRID", {**GRID, "rhos": [0.1, 0.3]})
    pm.worker(dict(PAIR))
    assert len(calls) == 2                                   # grid changed: recomputed
    monkeypatch.setattr(pm.ria, "_ATLAS_SHA", "a2", raising=False)
    pm.worker(dict(PAIR))
    assert len(calls) == 3                                   # atlas changed: recomputed
    f = tmp_path / "SYN-001_1_2.json"
    rec = json.loads(f.read_text())
    del rec["cache_sig"]                                     # legacy record without a signature
    f.write_text(json.dumps(rec))
    pm.worker(dict(PAIR))
    assert len(calls) == 4


# ---- fold-wise cell selection uses training patients only ----

def _write_cache(out, recs):
    c = out / "cache_h2"
    c.mkdir(parents=True)
    for i, r in enumerate(recs):
        (c / f"{r['patient_id']}_{i}.json").write_text(json.dumps(r))


def test_select_cells_uses_training_folds_only(tmp_path, monkeypatch):
    split = rule_split()
    man = tmp_path / "man"
    man.mkdir()
    split.to_csv(man / "split_mu.csv", index=False)
    monkeypatch.setattr(pm, "MAN", man)
    monkeypatch.setattr(pm, "OUT", tmp_path)
    a, b = "aniso|r=1|d=0.01|rho=0.1", "aniso|r=1|d=0.03|rho=0.1"
    # fold-0 patients strongly prefer b; all others mildly prefer a. If fold 0 leaked into its own selection,
    # the pooled mean would flip fold 0 to b. Training-only selection must pick a for fold 0.
    recs = [{"patient_id": p, "dice": {"no_change": 0.5, "grid": ({a: 0.0, b: 1.0} if f == 0 else {a: 0.6, b: 0.5})}}
            for p, f in zip(split["patient_id"], split["fold"])]
    _write_cache(tmp_path, recs)
    sel = pm.select_cells()
    assert sel["0"]["aniso|r=1"] == a
    assert set(sel) == {str(k) for k in range(pm.N_FOLDS)}


def test_select_cells_refuses_unmapped_or_mixed_grid(tmp_path, monkeypatch):
    man = tmp_path / "man"
    man.mkdir()
    rule_split().to_csv(man / "split_mu.csv", index=False)
    monkeypatch.setattr(pm, "MAN", man)
    monkeypatch.setattr(pm, "OUT", tmp_path / "o1")
    a, b = "aniso|r=1|d=0.01|rho=0.1", "aniso|r=1|d=0.03|rho=0.1"
    _write_cache(tmp_path / "o1", [{"patient_id": "UNSPLIT-1", "dice": {"no_change": 0.5, "grid": {a: 0.1, b: 0.2}}}])
    with pytest.raises(AssertionError, match="no fold"):
        pm.select_cells()
    monkeypatch.setattr(pm, "OUT", tmp_path / "o2")
    _write_cache(tmp_path / "o2", [{"patient_id": PIDS[0], "dice": {"no_change": 0.5, "grid": {a: 0.1, b: 0.2}}},
                                   {"patient_id": PIDS[1], "dice": {"no_change": 0.5, "grid": {a: 0.1}}}])
    with pytest.raises(AssertionError, match="different grid"):
        pm.select_cells()
