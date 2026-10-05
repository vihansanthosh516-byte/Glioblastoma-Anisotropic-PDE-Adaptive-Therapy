"""Leakage checks (analysis_plan_v1.md A1.1, masterplan §103).

Importable helpers. Every experiment script calls these before fitting; tests call them on the
committed manifests. Each helper raises AssertionError with the offending IDs.
"""
from __future__ import annotations

import pandas as pd


def assert_disjoint(*named_sets: tuple[str, set]) -> None:
    """No patient ID may appear in two of the named sets (train / validation / test)."""
    names = [n for n, _ in named_sets]
    for i in range(len(named_sets)):
        for j in range(i + 1, len(named_sets)):
            both = named_sets[i][1] & named_sets[j][1]
            assert not both, f"patients in both {names[i]} and {names[j]}: {sorted(both)[:5]}"


def assert_no_duplicate_scans(scans: pd.DataFrame) -> None:
    d = scans.duplicated(["patient_id", "dataset", "scan_number"])
    assert not d.any(), f"duplicate scan rows: {scans.loc[d, ['patient_id', 'scan_number']].head().values.tolist()}"


def assert_unique_patients(patients: pd.DataFrame) -> None:
    d = patients.duplicated(["dataset", "patient_id"])
    assert not d.any(), f"duplicate patient rows: {patients.loc[d, 'patient_id'].head().tolist()}"


def assert_one_fold_per_patient(split: pd.DataFrame) -> None:
    assert split["patient_id"].is_unique, "a patient has more than one fold"
    assert split["fold"].between(0, split["fold"].max()).all() and (split["fold"] >= 0).all(), "bad fold id"


def assert_pairs_follow_split(pairs: pd.DataFrame, split: pd.DataFrame) -> None:
    """Every forecast pair belongs to a patient with exactly one fold (so pairs cannot straddle folds)."""
    missing = set(pairs["patient_id"]) - set(split["patient_id"])
    assert not missing, f"pairs for patients with no fold: {sorted(missing)[:5]}"


def assert_causal_pairs(pairs: pd.DataFrame) -> None:
    """Forecast input is strictly earlier than its target; later scans never feed earlier forecasts."""
    bad = pairs[pairs["day_in"] >= pairs["day_out"]]
    assert bad.empty, f"non-causal pairs (input not before target): {bad[['patient_id', 'tp_in', 'tp_out']].head().values.tolist()}"
    assert (pairs["dt_days"] == pairs["day_out"] - pairs["day_in"]).all(), "dt_days inconsistent with day_out - day_in"
    for pid, g in pairs.groupby("patient_id"):
        g = g.sort_values("day_in")
        assert (g["day_out"].values[:-1] <= g["day_out"].values[1:]).all(), f"{pid}: pair order not temporal"


def assert_train_before_test_fit(fit_ids: set, eval_ids: set, what: str = "fit") -> None:
    """A fitted object (scaler, atlas, graph, threshold) must be built from patients disjoint from its evaluation set."""
    both = fit_ids & eval_ids
    assert not both, f"{what} built using evaluation patients: {sorted(both)[:5]}"


def fold_sets(split: pd.DataFrame, k: int) -> tuple[set, set]:
    """(train patients, test patients) for fold k of a patient-level split."""
    test = set(split.loc[split["fold"] == k, "patient_id"])
    train = set(split.loc[split["fold"] != k, "patient_id"])
    return train, test
