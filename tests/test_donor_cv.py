"""Donor identity and nested donor-grouped CV (src/donor_cv.py) on synthetic data."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from src import donor_cv as dc  # noqa: E402


def acc(a, b):
    return float(np.mean(np.asarray(a) == np.asarray(b)))


def synth(n_donors=12, cells=20, seed=0):
    rng = np.random.default_rng(seed)
    groups = np.repeat([f"D{i:02d}" for i in range(n_donors)], cells)
    y = np.tile(np.arange(cells) % 2, n_donors)
    rows = np.arange(len(y), dtype=float)
    X = np.column_stack([rows, y + rng.normal(0, 0.3, len(y))])   # column 0 = row id, column 1 = signal
    return X, y, groups


def test_donor_map_resolves_aliases_and_refuses_conflicts():
    s = pd.DataFrame({"sample_id": ["s1", "s2", "s3"], "donor_id": ["P1", "P1-old", "P2"]})
    m = dc.donor_map(s, aliases={"P1-old": "P1"})
    assert m.to_dict() == {"s1": "P1", "s2": "P1", "s3": "P2"}
    bad = pd.DataFrame({"sample_id": ["s1", "s1"], "donor_id": ["P1", "P2"]})
    with pytest.raises(ValueError, match="more than one donor"):
        dc.donor_map(bad)
    with pytest.raises(ValueError, match="blank"):
        dc.donor_map(pd.DataFrame({"sample_id": ["s1"], "donor_id": [" "]}))
    with pytest.raises(ValueError, match="themselves aliases"):
        dc.donor_map(s, aliases={"a": "b", "b": "c"})
    # same sample listed twice under two aliases of one donor is fine
    dup = pd.DataFrame({"sample_id": ["s1", "s1"], "donor_id": ["P1", "P1-old"]})
    assert dc.donor_map(dup, aliases={"P1-old": "P1"}).to_dict() == {"s1": "P1"}


def test_grouped_folds_partition_donors():
    _, y, g = synth()
    folds = dc.grouped_folds(y, g, 4, seed=1)
    test_donors = [d for _, te in folds for d in np.unique(g[te])]
    assert sorted(test_donors) == sorted(np.unique(g))
    with pytest.raises(ValueError, match="cannot fill"):
        dc.grouped_folds(y, g, 13, seed=1)
    with pytest.raises(ValueError, match="missing donor"):
        dc.grouped_folds([0, 1], np.array(["a", None], dtype=object), 2, seed=1)


def test_grouped_folds_supports_leave_one_donor_out_and_checks_lengths():
    _, y, groups = synth(n_donors=6, cells=10)
    folds = dc.grouped_folds(y, groups, "lodo", seed=1)
    assert len(folds) == 6
    assert all(len(np.unique(groups[test])) == 1 for _, test in folds)
    with pytest.raises(ValueError, match="equal length"):
        dc.grouped_folds([0, 1], ["a"], 2, seed=1)


def test_nested_cv_supports_leave_one_donor_out_at_both_levels():
    X, y, groups = synth(n_donors=6, cells=20)

    def fit_predict(Xa, ya, Xb, hp, seed):
        return (Xb[:, 1] > 0.5).astype(int)

    result = dc.nested_grouped_cv(X, y, groups, [None], fit_predict, acc,
                                  outer="lodo", inner="lodo", seed=4)
    assert len(result) == 6
    assert all(len(row["test_donors"]) == 1 for row in result)


def test_nested_cv_never_fits_on_outer_test_rows():
    X, y, g = synth()
    calls = []

    def fit_predict(Xa, ya, Xb, hp, seed):
        calls.append(set(Xa[:, 0].astype(int)))
        thr = np.median(Xa[:, 1]) + hp
        return (Xb[:, 1] > thr).astype(int)

    res = dc.nested_grouped_cv(X, y, g, [0.0, 0.4], fit_predict, acc, outer=4, inner=3, seed=5)
    assert len(res) == 4
    pos = 0
    for r in res:
        test_rows = set(r["test_index"].tolist())
        n_calls = 2 * 3 + 1                   # grid x inner folds, then one refit
        for fit_rows in calls[pos:pos + n_calls]:
            assert not fit_rows & test_rows
        pos += n_calls
        assert not set(r["test_donors"]) & set(np.unique(g[np.setdiff1d(np.arange(len(y)), r["test_index"])]))
    assert pos == len(calls)


def test_nested_cv_selects_best_inner_hp():
    X, y, g = synth()

    def fit_predict(Xa, ya, Xb, hp, seed):   # hp 'flip' is always wrong
        p = (Xb[:, 1] > 0.5).astype(int)
        return 1 - p if hp == "flip" else p

    res = dc.nested_grouped_cv(X, y, g, ["flip", "ok"], fit_predict, acc, outer=3, inner=3, seed=0)
    assert all(r["selected_hp"] == "ok" for r in res)
    assert all(r["outer_score"] > 0.8 for r in res)


def test_nested_cv_can_minimize_lower_is_better_score():
    X, y, g = synth(n_donors=8, cells=20)

    def fit_predict(Xa, ya, Xb, hp, seed):
        pred = (Xb[:, 1] > 0.5).astype(int)
        return 1 - pred if hp == "bad" else pred

    def error(y_true, y_pred):
        return float(np.mean(np.asarray(y_true) != np.asarray(y_pred)))

    res = dc.nested_grouped_cv(X, y, g, ["bad", "good"], fit_predict, error,
                               outer=4, inner=3, seed=2, score_direction="minimize")
    assert all(row["selected_hp"] == "good" for row in res)
    assert all(row["score_direction"] == "minimize" for row in res)


def test_donor_macro_score_weights_donors_equally():
    y = np.array([1, 1, 0, 0, 0, 0])
    pred = np.array([0, 0, 0, 0, 0, 0])
    groups = np.array(["small", "small", "large", "large", "large", "large"])
    macro, by_donor = dc._score_fold(y, pred, groups, acc, donor_macro=True)
    assert macro == 0.5
    assert by_donor == {"large": 1.0, "small": 0.0}


def test_nested_cv_enforces_minimum_two_class_training_donors():
    donors = np.repeat([f"D{i}" for i in range(6)], 8)
    y = np.repeat(np.arange(6) % 2, 8)  # every donor has only one class
    X = np.column_stack([np.arange(len(y)), y]).astype(float)

    def fit_predict(Xa, ya, Xb, hp, seed):
        return np.zeros(len(Xb), dtype=int)

    with pytest.raises(ValueError, match="two-class donors"):
        dc.nested_grouped_cv(X, y, donors, [None], fit_predict, acc, outer=3, inner=2)


def test_donor_lookup_does_not_transfer_across_grouped_folds():
    """A model that memorises donor -> label gets near-perfect accuracy if donors leak, and chance under
    donor-grouped CV. Here each donor has one label, so donor identity fully predicts the label."""
    rng = np.random.default_rng(3)
    donors = np.repeat(np.arange(20), 10)
    y = donors % 2
    X = np.column_stack([np.arange(len(y)), donors, rng.normal(size=len(y))]).astype(float)

    def lookup(Xa, ya, Xb, hp, seed):
        table = dict(zip(Xa[:, 1], ya))
        return np.array([table.get(d, 0) for d in Xb[:, 1]])

    res = dc.nested_grouped_cv(X, y, donors.astype(str), [None], lookup, acc, outer=5, inner=2, seed=0,
                               min_train_two_class_donors=0)
    assert np.mean([r["outer_score"] for r in res]) <= 0.6


def test_per_donor_scores():
    s = dc.per_donor_scores([1, 1, 0, 0], [1, 0, 0, 0], ["a", "a", "b", "b"], acc)
    assert s.to_dict() == {"a": 0.5, "b": 1.0}
