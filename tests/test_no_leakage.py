"""Leakage and split tests on the committed manifests (analysis_plan_v1.md A1.1-A1.3)."""
from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
MAN = REPO / "data" / "manifests"

from src import leakage_checks as lc  # noqa: E402

spec = importlib.util.spec_from_file_location("manifest_script", REPO / "src" / "94_patient_manifest.py")
ms = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ms)


@pytest.fixture(scope="module")
def patients():
    return pd.read_csv(MAN / "patient_manifest.csv", low_memory=False)


@pytest.fixture(scope="module")
def scans():
    return pd.read_csv(MAN / "scan_manifest.csv", low_memory=False)


@pytest.fixture(scope="module")
def pairs():
    return pd.read_csv(MAN / "forecast_pairs_mu.csv")


@pytest.fixture(scope="module")
def split():
    return pd.read_csv(MAN / "split_mu.csv")


def test_unique_patients(patients):
    lc.assert_unique_patients(patients)


def test_no_duplicate_scans(scans):
    lc.assert_no_duplicate_scans(scans)


def test_one_fold_per_patient(split):
    lc.assert_one_fold_per_patient(split)


def test_pairs_follow_split(pairs, split):
    lc.assert_pairs_follow_split(pairs, split)


def test_every_fold_disjoint_train_test(split):
    for k in range(ms.N_FOLDS):
        train, test = lc.fold_sets(split, k)
        lc.assert_disjoint(("train", train), ("test", test))
        assert train and test


def test_causal_pairs(pairs):
    lc.assert_causal_pairs(pairs)


def test_split_is_deterministic_and_order_free(split):
    """Fold depends only on (seed, patient_id), not on row order or cohort size."""
    for pid, fold in zip(split["patient_id"].head(50), split["fold"].head(50)):
        h = int(hashlib.sha256(f"{ms.SPLIT_SEED}|{pid}".encode()).hexdigest(), 16) % ms.N_FOLDS
        assert h == fold
        assert ms.fold_of(pid) == fold


def test_folds_are_balanced(split):
    counts = split["fold"].value_counts()
    assert counts.min() >= 0.6 * counts.max()


def test_eligible_patients_have_primary_pairs(patients, pairs):
    mu = patients[patients["dataset"] == "MU"]
    prim = pairs[pairs["in_primary"]].groupby("patient_id").size()
    for pid in mu.loc[mu["eligible"] == True, "patient_id"]:  # noqa: E712
        assert prim.get(pid, 0) > 0, f"{pid} eligible with no primary pair"


def test_primary_pairs_respect_window(pairs):
    p = pairs[pairs["in_primary"]]
    assert p["dt_days"].between(ms.DT_MIN, ms.DT_MAX).all()
    assert p["input_core_nonempty"].all()


def test_exclusions_all_have_a_reason(patients):
    mu = patients[patients["dataset"] == "MU"]
    bad = mu[(mu["eligible"] == False) & (mu["reason_excluded"].isna() | (mu["reason_excluded"] == ""))]  # noqa: E712
    assert bad.empty


def test_helpers_catch_leaks():
    with pytest.raises(AssertionError):
        lc.assert_disjoint(("train", {"a", "b"}), ("test", {"b", "c"}))
    with pytest.raises(AssertionError):
        lc.assert_train_before_test_fit({"a", "b"}, {"b"}, "scaler")
    bad = pd.DataFrame({"patient_id": ["p"], "tp_in": [1], "tp_out": [2], "day_in": [10], "day_out": [5], "dt_days": [-5]})
    with pytest.raises(AssertionError):
        lc.assert_causal_pairs(bad)
    dup = pd.DataFrame({"patient_id": ["p", "p"], "dataset": ["MU", "MU"], "scan_number": [1, 1]})
    with pytest.raises(AssertionError):
        lc.assert_no_duplicate_scans(dup)
