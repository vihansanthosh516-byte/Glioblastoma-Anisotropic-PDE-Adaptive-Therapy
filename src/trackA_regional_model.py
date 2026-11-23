"""Development model for paired regional prediction with a spatial fallback."""

from dataclasses import dataclass
from math import ceil

import numpy as np
from sklearn.linear_model import LogisticRegression


NONSPATIAL = ("composition", "density_depth", "malignant_expression")
SPATIAL = ("nominal_spatial", "alternative_spatial", "molecule_spatial")
ARMS = {
    "composition": ("composition",),
    "density_depth": ("density_depth",),
    "malignant_expression": ("malignant_expression",),
    "expression_composition": ("malignant_expression", "composition"),
    "nonspatial": NONSPATIAL,
    **{name: (*NONSPATIAL, name) for name in SPATIAL},
}


def _blocks(features):
    if set(features) != set(NONSPATIAL + SPATIAL):
        raise ValueError("Provide all six named feature blocks")
    arrays = {name: np.asarray(value, dtype=float) for name, value in features.items()}
    if any(a.ndim != 2 or 0 in a.shape for a in arrays.values()):
        raise ValueError("Each block must be a nonempty pair-by-feature matrix")
    if len({len(a) for a in arrays.values()}) != 1:
        raise ValueError("All blocks must use the same pairs in the same order")
    if any(not np.isfinite(arrays[name]).all() for name in NONSPATIAL):
        raise ValueError("Nonspatial features are required for every pair")
    if arrays["nominal_spatial"].shape != arrays["alternative_spatial"].shape:
        raise ValueError("Reconstructions must share the same feature definitions")
    return arrays


def _join(blocks, names):
    return np.column_stack([blocks[name] for name in names])


def _fit_margin(features, orientation, strength):
    # Every canonical pair is mirrored. It contributes equal weight regardless
    # of which of its specimen identifiers sorts first.
    oriented = orientation[:, None] * features
    scale = np.sqrt(np.mean(oriented**2, axis=0))
    scale[scale == 0] = 1.0
    scaled = oriented / scale
    model = LogisticRegression(C=strength, fit_intercept=False, max_iter=2000)
    model.fit(np.concatenate((scaled, -scaled)), np.repeat((1, 0), len(scaled)))
    return model.coef_[0] / scale


@dataclass(frozen=True)
class RegionalModel:
    coefficients: dict[str, np.ndarray]
    block_widths: dict[str, int]

    def predict(self, features, *, bound, audit_qualified):
        """Return signed margins and eligibility, before revealing region labels.

        A common bound is in units of the fitted shared margin and covers both
        reconstructions jointly. NaN or infinity means unavailable. Qualify its
        application separately for each pair using a prediction-time audit.
        """
        blocks = _blocks(features)
        if {name: value.shape[1] for name, value in blocks.items()} != self.block_widths:
            raise ValueError("Prediction feature blocks must match the trained widths")
        scores, available = {}, {}
        for arm, names in ARMS.items():
            margin = _join(blocks, names) @ self.coefficients[arm]
            available[arm] = np.isfinite(margin)
            scores[arm] = (
                np.where(available[arm], margin, scores["nonspatial"])
                if arm in SPATIAL else margin
            )

        reconstructed = self.reconstruction_margins(blocks)
        decisions = spatial_fallback(scores["nonspatial"], reconstructed,
                                     bound=bound, audit_qualified=audit_qualified)
        for i, name in enumerate(("shared_nominal", "shared_alternative")):
            available[name] = np.isfinite(reconstructed[:, i])
            scores[name] = np.where(available[name], reconstructed[:, i], scores["nonspatial"])
        scores.update(decisions["scores"])
        available.update(decisions["used_spatial"])
        return scores, available

    def reconstruction_margins(self, features):
        """Apply the same coefficients to each of the two reconstructions."""
        blocks = _blocks(features)
        if {name: value.shape[1] for name, value in blocks.items()} != self.block_widths:
            raise ValueError("Reconstruction feature blocks must match the trained widths")
        return np.column_stack([
            _join(blocks, (*NONSPATIAL, name)) @ self.coefficients["shared"]
            for name in SPATIAL[:2]
        ])


def fit_regional(features, orientation, *, strength):
    """Fit only on development pairs with all arms available."""
    blocks = _blocks(features)
    orientation = np.asarray(orientation, dtype=float)
    n = len(blocks["composition"])
    if orientation.shape != (n,) or not np.isin(orientation, (-1, 1)).all():
        raise ValueError("Each development pair needs a signed infiltrating/core order")
    if not np.isfinite(strength) or strength <= 0:
        raise ValueError("Require a positive finite regularization setting")
    if any(not np.isfinite(a).all() for a in blocks.values()):
        raise ValueError("All training arms require the same complete development pairs")

    weights = {
        arm: _fit_margin(_join(blocks, names), orientation, strength)
        for arm, names in ARMS.items()
    }
    both = np.concatenate([_join(blocks, (*NONSPATIAL, name)) for name in SPATIAL[:2]])
    weights["shared"] = _fit_margin(both, np.tile(orientation, 2), strength / 2)
    return RegionalModel(weights, {name: block.shape[1] for name, block in blocks.items()})


def calibrate_margin_bound(reconstructed, reference, *, alpha):
    """Calibrate from independent reference margins on disjoint development pairs.

    The reference uses the same score scale. One maximum absolute reconstruction
    error counts per donor. Split-conformal coverage requires exchangeability
    with the target condition. An infinite bound records inadequate calibration.
    """
    margins = np.asarray(reconstructed, dtype=float)
    reference = np.asarray(reference, dtype=float)
    if margins.ndim != 2 or margins.shape[1] != 2 or reference.shape != (len(margins),):
        raise ValueError("Provide one reference and two reconstructions per development pair")
    if not 0 < alpha < 1 or not np.isfinite(margins).all() or not np.isfinite(reference).all():
        raise ValueError("Require finite calibration margins and 0 < alpha < 1")
    rank = ceil((len(reference) + 1) * (1 - alpha))
    if rank > len(reference):
        return float("inf")
    errors = np.max(np.abs(margins - reference[:, None]), axis=1)
    return float(np.partition(errors, rank - 1)[rank - 1])


def spatial_fallback(nonspatial, reconstructed, *, bound, audit_qualified):
    """Use a spatial margin when both uncertainty intervals certify its sign."""
    fallback = np.asarray(nonspatial, dtype=float)
    margins = np.asarray(reconstructed, dtype=float)
    audit = np.asarray(audit_qualified)
    radius = np.asarray(bound, dtype=float)
    if fallback.ndim != 1 or margins.shape != (len(fallback), 2):
        raise ValueError("Require one nonspatial and two reconstructed margins per pair")
    if radius.ndim == 0:
        radius = np.full(len(fallback), radius)
    if (radius.shape != fallback.shape or audit.shape != fallback.shape
            or audit.dtype != np.dtype(bool) or np.any(radius < 0)
            or not np.isfinite(fallback).all()):
        raise ValueError("Require a nonnegative bound, strict audit flag and finite fallback per pair")
    valid = np.isfinite(margins).all(axis=1)
    agreement = valid & ((margins > 0).all(axis=1) | (margins < 0).all(axis=1))
    finite_bound = np.isfinite(radius)
    certified = agreement & audit & finite_bound & (np.abs(margins) > radius[:, None]).all(axis=1)
    # Both reconstructions agree where the nominal value is used, so choosing
    # the nominal margin also keeps the compare-to-nominal ablation direct.
    return {
        "scores": {
            "agreement_only": np.where(agreement, margins[:, 0], fallback),
            "bounded_spatial": np.where(certified, margins[:, 0], fallback),
        },
        "used_spatial": {"agreement_only": agreement, "bounded_spatial": certified},
    }
