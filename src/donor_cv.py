"""Donor identity and donor-grouped nested cross-validation for Track A (no data I/O).

Cells or spots from one donor must never sit on both sides of a fit/evaluate boundary. This module gives one
place to (1) resolve sample -> donor identity, including aliases, and refuse conflicts, and (2) run nested
donor-grouped CV: hyperparameters are chosen on inner folds built only from the outer-training donors, then
refit on all outer-training donors and scored once on the outer-test donors.

fit_predict(X_train, y_train, X_eval, hp, seed) must do all preprocessing (scaling, PCA, gene selection, graph
building) inside the call, from X_train only; this module cannot see inside it.
"""
from __future__ import annotations

from typing import Callable, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd
from sklearn.model_selection import LeaveOneGroupOut, StratifiedGroupKFold

from src import leakage_checks as lc


def donor_map(samples: pd.DataFrame, sample_col: str = "sample_id", donor_col: str = "donor_id",
              aliases: Mapping[str, str] | None = None) -> pd.Series:
    """Validated sample -> canonical donor Series.

    aliases maps alternative donor labels to one canonical label (e.g. two IDs for the same patient in two
    releases). Refuses blank IDs, a sample assigned to two donors, an alias whose target is itself an alias,
    and an alias that would merge a donor with an unrelated canonical ID listed elsewhere in aliases.
    """
    aliases = dict(aliases or {})
    chained = sorted(a for a, b in aliases.items() if b in aliases)
    if chained:
        raise ValueError(f"alias targets that are themselves aliases: {chained[:5]}")
    df = samples[[sample_col, donor_col]].copy()
    if df.isna().any().any() or (df.astype(str).apply(lambda c: c.str.strip()) == "").any().any():
        raise ValueError("blank sample or donor ID")
    df[donor_col] = df[donor_col].astype(str).map(lambda d: aliases.get(d, d))
    n = df.drop_duplicates().groupby(sample_col)[donor_col].nunique()
    bad = n[n > 1].index.tolist()
    if bad:
        raise ValueError(f"samples assigned to more than one donor: {bad[:5]}")
    return df.drop_duplicates().set_index(sample_col)[donor_col]


def grouped_folds(y: Sequence, groups: Sequence, n_splits: int | str, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """Build checked grouped folds. Use ``n_splits='lodo'`` for leave-one-donor-out,
    or an integer for shuffled StratifiedGroupKFold. Every group appears in exactly
    one test fold; train and test groups are disjoint."""
    y, groups = np.asarray(y), np.asarray(groups)
    if y.ndim != 1 or groups.ndim != 1:
        raise ValueError("labels and donor IDs must be one-dimensional")
    if len(y) != len(groups):
        raise ValueError("labels and donor IDs must have equal length")
    if len(y) == 0:
        raise ValueError("cannot split empty data")
    if pd.isna(groups).any():
        raise ValueError("missing donor ID")
    unique_groups = pd.unique(groups)
    n_groups = len(unique_groups)
    if n_splits == "lodo":
        folds = list(LeaveOneGroupOut().split(np.zeros(len(y)), y, groups))
    elif isinstance(n_splits, (int, np.integer)) and not isinstance(n_splits, (bool, np.bool_)):
        if n_splits < 2:
            raise ValueError("n_splits must be at least 2, or 'lodo'")
        if n_groups < n_splits:
            raise ValueError(f"{n_groups} donors cannot fill {n_splits} grouped folds")
        folds = list(StratifiedGroupKFold(n_splits=int(n_splits), shuffle=True, random_state=seed)
                     .split(np.zeros(len(y)), y, groups))
    else:
        raise ValueError("n_splits must be an integer or 'lodo'")
    seen: list = []
    for tr, te in folds:
        if len(tr) == 0 or len(te) == 0:
            raise ValueError("empty grouped fold")
        lc.assert_disjoint(("train", set(groups[tr])), ("test", set(groups[te])))
        seen += list(pd.unique(groups[te]))
    if len(seen) != n_groups or set(seen) != set(unique_groups):
        raise ValueError("a donor is in no test fold or in more than one")
    return folds


def _score_fold(y_true, y_pred, groups, score: Callable, donor_macro: bool) -> tuple[float, dict | None]:
    if not donor_macro:
        value = float(score(y_true, y_pred))
        if not np.isfinite(value):
            raise ValueError("scorer returned a non-finite value")
        return value, None
    by_donor = per_donor_scores(y_true, y_pred, groups, score)
    if not np.isfinite(by_donor.to_numpy(dtype=float)).all():
        raise ValueError("donor-macro scoring requires a finite score for every scored donor")
    values = {donor: float(value) for donor, value in by_donor.items()}
    return float(np.mean(list(values.values()))), values


def _two_class_group_count(y: np.ndarray, groups: np.ndarray) -> int:
    return sum(np.unique(y[groups == group]).size >= 2 for group in pd.unique(groups))


def nested_grouped_cv(X, y, groups, grid: Iterable, fit_predict: Callable, score: Callable, *,
                      outer: int | str = 5, inner: int | str = 3, seed: int = 0,
                      score_direction: str = "maximize", donor_macro: bool = False,
                      min_train_two_class_donors: int = 4) -> list[dict]:
    """Nested donor-grouped selection with explicit metric direction and weighting.

    ``outer`` or ``inner`` may be ``'lodo'``. With ``donor_macro=True``, donors receive
    equal weight across all inner validation folds and in each outer score. Otherwise,
    scores are averaged by fold as before. Metric ties resolve to the first grid item.
    """
    y, groups, X = np.asarray(y), np.asarray(groups), np.asarray(X)
    if X.ndim < 2 or len(X) != len(y) or len(y) != len(groups):
        raise ValueError("X, labels, and donor IDs must have matching row counts")
    if score_direction not in {"maximize", "minimize"}:
        raise ValueError("score_direction must be 'maximize' or 'minimize'")
    if not isinstance(min_train_two_class_donors, (int, np.integer)) or min_train_two_class_donors < 0:
        raise ValueError("min_train_two_class_donors must be a non-negative integer")
    grid = list(grid)
    if not grid:
        raise ValueError("hyperparameter grid must not be empty")
    out = []
    for k, (tr, te) in enumerate(grouped_folds(y, groups, outer, seed)):
        outer_train_groups = _two_class_group_count(y[tr], groups[tr])
        if outer_train_groups < min_train_two_class_donors:
            raise ValueError(f"outer training fold {k} has {outer_train_groups} two-class donors; "
                             f"requires {min_train_two_class_donors}")
        inner_folds = grouped_folds(y[tr], groups[tr], inner, seed + 1 + k)
        table = []
        for hp in grid:
            fold_scores, donor_scores = [], {}
            for a, b in inner_folds:
                inner_train_groups = _two_class_group_count(y[tr][a], groups[tr][a])
                if inner_train_groups < min_train_two_class_donors:
                    raise ValueError(f"inner training fold has {inner_train_groups} two-class donors; "
                                     f"requires {min_train_two_class_donors}")
                pred = fit_predict(X[tr][a], y[tr][a], X[tr][b], hp, seed)
                value, per_donor = _score_fold(y[tr][b], pred, groups[tr][b], score, donor_macro)
                fold_scores.append(value)
                if per_donor is not None:
                    if set(donor_scores) & set(per_donor):
                        raise ValueError("a donor appears in multiple inner validation folds")
                    donor_scores.update(per_donor)
            values = list(donor_scores.values()) if donor_macro else fold_scores
            table.append({"hp": hp, "inner_mean": float(np.mean(values)),
                          "inner_scores": [float(v) for v in fold_scores],
                          "inner_donor_scores": donor_scores if donor_macro else None})
        direction = 1 if score_direction == "maximize" else -1
        best = max(range(len(grid)), key=lambda i: (direction * table[i]["inner_mean"], -i))
        lc.assert_train_before_test_fit(set(groups[tr]), set(groups[te]), "outer refit")
        pred = fit_predict(X[tr], y[tr], X[te], grid[best], seed)
        outer_score, outer_donor_scores = _score_fold(y[te], pred, groups[te], score, donor_macro)
        out.append({"fold": k, "selected_hp": grid[best], "inner": table,
                    "outer_score": outer_score, "score_direction": score_direction,
                    "donor_macro": bool(donor_macro), "outer_donor_scores": outer_donor_scores,
                    "test_donors": sorted(map(str, pd.unique(groups[te]))), "n_test_rows": int(len(te)),
                    "pred": pred, "test_index": te})
    return out


def per_donor_scores(y_true, y_pred, groups, score: Callable) -> pd.Series:
    """Score each donor separately, so evaluation is reported per donor rather than pooled over cells
    (a donor with many cells cannot dominate)."""
    y_true, y_pred, groups = map(np.asarray, (y_true, y_pred, groups))
    if y_true.ndim != 1 or y_pred.ndim != 1 or groups.ndim != 1:
        raise ValueError("labels, predictions, and donor IDs must be one-dimensional")
    if not (len(y_true) == len(y_pred) == len(groups)):
        raise ValueError("labels, predictions, and donor IDs must have equal length")
    return pd.Series({g: score(y_true[groups == g], y_pred[groups == g]) for g in np.unique(groups)})
