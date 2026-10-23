"""Leakage checks (analysis_plan_v1.md A1.1, masterplan §103).

Importable helpers; callers must invoke them explicitly before fitting. Their
presence does not prove every experiment uses them or avoids test-label access.
Tests exercise committed manifests and crafted invalid inputs. Each helper raises
AssertionError with the offending IDs.
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


def assert_one_fold_per_patient(split: pd.DataFrame, n_folds: int | None = None) -> None:
    """One integer fold per patient, no blank IDs or folds; with n_folds, every fold 0..n_folds-1 is used."""
    assert split["patient_id"].notna().all(), "blank patient ID in split"
    assert split["patient_id"].is_unique, "a patient has more than one fold"
    fold = pd.to_numeric(split["fold"], errors="coerce")
    bad = split.loc[fold.isna() | (fold != fold.round()) | (fold < 0), "patient_id"]
    assert bad.empty, f"bad fold id for patients: {bad.head().tolist()}"
    if n_folds is not None:
        out = split.loc[fold >= n_folds, "patient_id"]
        assert out.empty, f"fold >= {n_folds} for patients: {out.head().tolist()}"
        empty = sorted(set(range(n_folds)) - set(fold.astype(int)))
        assert not empty, f"empty folds: {empty}"


def assert_split_matches_rule(split: pd.DataFrame, fold_fn) -> None:
    """Every stored fold equals fold_fn(patient_id); catches a hand-edited or stale split file."""
    bad = [(p, f) for p, f in zip(split["patient_id"], split["fold"]) if fold_fn(p) != f]
    assert not bad, f"split differs from its rule (patient, stored fold): {bad[:5]}"


def fold_map(split: pd.DataFrame, n_folds: int | None = None) -> pd.Series:
    """Validated patient -> int fold Series, for callers that map folds onto pairs or records."""
    assert_one_fold_per_patient(split, n_folds)
    return split.set_index("patient_id")["fold"].astype(int)


def assert_all_have_fold(patient_ids, folds: pd.Series, what: str = "records") -> None:
    """Every patient in a pair table or result cache has a fold; an unmapped patient would count as training data in
    every fold."""
    missing = sorted(set(patient_ids) - set(folds.index))
    assert not missing, f"{what} for patients with no fold: {missing[:5]}"


def assert_pairs_follow_split(pairs: pd.DataFrame, split: pd.DataFrame) -> None:
    """Every forecast pair belongs to a patient with exactly one fold (so pairs cannot straddle folds)."""
    missing = set(pairs["patient_id"]) - set(split["patient_id"])
    assert not missing, f"pairs for patients with no fold: {sorted(missing)[:5]}"
    if "fold" in pairs.columns:
        want = pairs["patient_id"].map(split.set_index("patient_id")["fold"])
        bad = pairs.loc[want != pairs["fold"], "patient_id"]
        assert bad.empty, f"pair fold differs from split: {sorted(set(bad))[:5]}"


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
