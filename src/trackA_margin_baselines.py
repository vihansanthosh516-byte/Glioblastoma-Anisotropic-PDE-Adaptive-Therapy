"""Portable fixed-signature and paired-centroid baselines for Track A."""

from collections.abc import Mapping

import numpy as np


def signature_weights(genes, weights: Mapping[str, float], *, minimum_coverage: float):
    """Match frozen signed weights against an already harmonized assay universe."""
    genes = tuple(genes)
    if len(genes) != len(set(genes)):
        raise ValueError("Duplicate assay gene symbols")
    if not weights or not 0 < minimum_coverage <= 1:
        raise ValueError("Require a nonempty signature and a positive coverage threshold")
    if any(not np.isfinite(w) or w == 0 for w in weights.values()):
        raise ValueError("Signature weights must be finite and nonzero")
    vector = np.array([weights.get(gene, 0.0) for gene in genes], dtype=float)
    coverage = np.count_nonzero(vector) / len(weights)
    if coverage < minimum_coverage:
        raise ValueError("Frozen signature gene coverage was missed")
    return vector / np.abs(vector).sum(), coverage


def score_linear(expression, weights):
    """Apply fixed weights to already transformed specimen-by-feature rows."""
    expression = np.asarray(expression, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if expression.ndim != 2 or weights.shape != (expression.shape[1],):
        raise ValueError("Specimen rows and weight vector must share the feature axis")
    if not np.isfinite(expression).all() or not np.isfinite(weights).all():
        raise ValueError("Expression and weights must be finite")
    return expression @ weights


def fit_paired_centroid(core, margin):
    """Fit a linear margin score from equally weighted development-donor contrasts."""
    core, margin = np.asarray(core, dtype=float), np.asarray(margin, dtype=float)
    if core.ndim != 2 or core.shape != margin.shape or 0 in core.shape:
        raise ValueError("Require paired nonempty donor-by-feature arrays")
    if not np.isfinite(core).all() or not np.isfinite(margin).all():
        raise ValueError("Development features must be finite")
    return (margin - core).mean(axis=0)
