"""Fit and score the Track A state-consistent gene-pair method."""

from dataclasses import dataclass

import numpy as np


def _expression_array(expression, ndim):
    values = np.asarray(expression, dtype=float)
    if values.ndim != ndim or not np.isfinite(values).all():
        raise ValueError(f"Expected a finite {ndim}-dimensional expression array")
    return values


def _votes(expression, first, second):
    return np.sign(expression[..., first] - expression[..., second])


@dataclass(frozen=True)
class MarginPairModel:
    """Frozen indices and orientation in the training state and gene order."""

    first: np.ndarray
    second: np.ndarray
    directions: np.ndarray
    shape: tuple[int, int]

    def score(self, expression):
        """Score specimens shaped (specimen, state, gene)."""
        values = _expression_array(expression, 3)
        if values.shape[1:] != self.shape:
            raise ValueError("Specimen state and gene counts must match training")
        if self.directions.size == 0:
            return np.zeros(len(values))
        votes = _votes(values, self.first, self.second)
        return (votes * self.directions).mean(axis=(1, 2))


def fit_margin_pairs(core, margin, genes, count, *, weakest_state=True):
    """Fit on paired arrays shaped (donor, state, gene).

    Callers must harmonize the gene and state axes before fitting and scoring.
    The mean-change ablation uses the same eligibility and scoring rules.
    """
    core = _expression_array(core, 3)
    margin = _expression_array(margin, 3)
    genes = tuple(genes)
    if core.shape != margin.shape or 0 in core.shape:
        raise ValueError("Training arrays must match on nonempty donor, state and gene axes")
    if len(genes) != core.shape[2] or len(set(genes)) != len(genes):
        raise ValueError("Provide unique gene symbols for the gene axis")
    if not isinstance(count, int) or count < 1:
        raise ValueError("Pair count must be a positive integer")
    order = np.argsort(genes)
    left, right = np.triu_indices(len(genes), k=1)
    first, second = order[left], order[right]
    strengths = np.empty(len(first))
    orientations = np.empty(len(first))
    # Bound temporary vote arrays when the prefilter retains many genes.
    for start in range(0, len(first), 4096):
        block = slice(start, start + 4096)
        change = _votes(margin, first[block], second[block])
        change -= _votes(core, first[block], second[block])
        state_change = change.mean(axis=0)
        agreement = (state_change > 0).all(axis=0) | (state_change < 0).all(axis=0)
        absolute = np.abs(state_change)
        strength = absolute.min(axis=0) if weakest_state else absolute.mean(axis=0)
        strengths[block] = np.where(agreement, strength, 0)
        orientations[block] = np.sign(state_change[0])
    ranking = np.argsort(-strengths, kind="stable")
    selected = ranking[strengths[ranking] > 0][:count]
    return MarginPairModel(first[selected], second[selected], orientations[selected], core.shape[1:])
