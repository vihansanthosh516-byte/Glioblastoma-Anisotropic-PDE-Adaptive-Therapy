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
from sklearn.model_selection import StratifiedGroupKFold

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


def grouped_folds(y: Sequence, groups: Sequence, n_splits: int, seed: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """StratifiedGroupKFold folds with checks: each group is in exactly one test fold, train/test donors are
    disjoint, and no fold is empty."""
    y, groups = np.asarray(y), np.asarray(groups)
    if pd.isna(groups).any():
        raise ValueError("missing donor ID")
    n_groups = len(np.unique(groups))
    if n_groups < n_splits:
        raise ValueError(f"{n_groups} donors cannot fill {n_splits} grouped folds")
    folds = list(StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
                 .split(np.zeros(len(y)), y, groups))
    seen: list = []
    for tr, te in folds:
        if len(tr) == 0 or len(te) == 0:
            raise ValueError("empty grouped fold")
        lc.assert_disjoint(("train", set(groups[tr])), ("test", set(groups[te])))
        seen += list(np.unique(groups[te]))
    if sorted(seen) != sorted(np.unique(groups)):
        raise ValueError("a donor is in no test fold or in more than one")
    return folds


def nested_grouped_cv(X, y, groups, grid: Iterable, fit_predict: Callable, score: Callable, *,
                      outer: int = 5, inner: int = 3, seed: int = 0) -> list[dict]:
    """Per outer fold: choose hp by mean inner-fold score on outer-training donors only (ties: first in grid),
    refit on all outer-training rows, predict outer-test rows once. Returns one record per outer fold with
    the selected hp, the inner score table, the outer score and the outer-test donors."""
    y, groups = np.asarray(y), np.asarray(groups)
    grid = list(grid)
    out = []
    for k, (tr, te) in enumerate(grouped_folds(y, groups, outer, seed)):
        inner_folds = grouped_folds(y[tr], groups[tr], inner, seed + 1 + k)
        table = []
        for hp in grid:
            s = [score(y[tr][b], fit_predict(X[tr][a], y[tr][a], X[tr][b], hp, seed)) for a, b in inner_folds]
            table.append({"hp": hp, "inner_mean": float(np.mean(s)), "inner_scores": [float(v) for v in s]})
        best = max(range(len(grid)), key=lambda i: (table[i]["inner_mean"], -i))
        lc.assert_train_before_test_fit(set(groups[tr]), set(groups[te]), "outer refit")
        pred = fit_predict(X[tr], y[tr], X[te], grid[best], seed)
        out.append({"fold": k, "selected_hp": grid[best], "inner": table, "outer_score": float(score(y[te], pred)),
                    "test_donors": sorted(map(str, np.unique(groups[te]))), "n_test_rows": int(len(te)),
                    "pred": pred, "test_index": te})
    return out


def per_donor_scores(y_true, y_pred, groups, score: Callable) -> pd.Series:
    """Score each donor separately, so evaluation is reported per donor rather than pooled over cells
    (a donor with many cells cannot dominate)."""
    y_true, y_pred, groups = map(np.asarray, (y_true, y_pred, groups))
    return pd.Series({g: score(y_true[groups == g], y_pred[groups == g]) for g in np.unique(groups)})
