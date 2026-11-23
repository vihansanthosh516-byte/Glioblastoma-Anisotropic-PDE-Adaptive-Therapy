"""Origin-level development benchmark for measured TMZ viability."""

from dataclasses import dataclass

import numpy as np
from sklearn.linear_model import Ridge
from sklearn.preprocessing import StandardScaler

from .trackA_margin_baselines import score_linear


def _matrix(values):
    result = np.asarray(values, dtype=float)
    if result.ndim != 2 or 0 in result.shape or not np.isfinite(result).all():
        raise ValueError("Require a finite, nonempty culture-by-feature matrix")
    return result


def _ids(values, n, label):
    result = np.asarray(values)
    if (result.shape != (n,) or result.dtype.kind != "U"
            or np.any(result == "") or np.any(result != np.char.strip(result))):
        raise ValueError(f"Require one nonempty canonical {label} ID per culture")
    return result


def _groups(origins):
    names, inverse, counts = np.unique(origins, return_inverse=True, return_counts=True)
    return names, inverse, counts


def _means(values, inverse, counts):
    sums = np.zeros((len(counts), *values.shape[1:]))
    np.add.at(sums, inverse, values)
    return sums / counts.reshape((-1,) + (1,) * (values.ndim - 1))


@dataclass(frozen=True)
class Cultures:
    """One row per independently derived culture after assay adjudication."""

    expression: np.ndarray
    viability: np.ndarray
    origins: np.ndarray
    culture_ids: np.ndarray

    def __post_init__(self):
        x = _matrix(self.expression)
        y = np.asarray(self.viability, dtype=float)
        origins = _ids(self.origins, len(x), "origin")
        cultures = _ids(self.culture_ids, len(x), "culture")
        if y.shape != (len(x),) or not np.isfinite(y).all():
            raise ValueError("Require matched finite viability for each culture")
        if len(np.unique(cultures)) != len(cultures):
            raise ValueError("Technical wells or duplicate cultures must be resolved before fitting")
        for name, value in zip(("expression", "viability", "origins", "culture_ids"),
                               (x, y, origins, cultures), strict=True):
            object.__setattr__(self, name, value)


def _fit_shared(x, y, true_origins, contrast_origins, penalty, within_weight):
    _, groups, counts = _groups(true_origins)
    mean_x, mean_y = _means(x, groups, counts), _means(y, groups, counts)
    center_x, center_y = mean_x.mean(axis=0), mean_y.mean()
    design, target = (mean_x - center_x) / np.sqrt(len(counts)), (mean_y - center_y) / np.sqrt(len(counts))
    if within_weight:
        _, contrast_groups, sizes = _groups(contrast_origins)
        repeated = sizes[contrast_groups] > 1
        repeated_count = np.count_nonzero(sizes > 1)
        if repeated_count == 0:
            raise ValueError("The within-origin term requires independently derived sister cultures")
        centered_x = x - _means(x, contrast_groups, sizes)[contrast_groups]
        centered_y = y - _means(y, contrast_groups, sizes)[contrast_groups]
        weight = np.sqrt(within_weight / (repeated_count * sizes[contrast_groups[repeated]]))
        design = np.vstack((design, centered_x[repeated] * weight[:, None]))
        target = np.concatenate((target, centered_y[repeated] * weight))
    fit = Ridge(alpha=penalty, fit_intercept=False).fit(design, target)
    return float(center_y - center_x @ fit.coef_), fit.coef_


def _shuffle_links(origins, strata, seed):
    # Only the contrast term uses these reassigned links. Assay and feature
    # rows stay together, and true origins still define the mean term.
    strata = _ids(strata, len(origins), "assay stratum")
    for origin in np.unique(origins):
        if np.unique(strata[origins == origin]).size != 1:
            raise ValueError("A training origin cannot span assay strata for this control")
    rng = np.random.default_rng(seed)
    shuffled = origins.copy()
    for stratum in np.unique(strata):
        mask = strata == stratum
        shuffled[mask] = rng.permutation(origins[mask])
    # Renaming whole groups does not disrupt the biological links.
    if np.unique(np.column_stack((origins, shuffled)), axis=0).shape[0] == np.unique(origins).size:
        raise ValueError("This shuffle did not break any sister links")
    return shuffled


def matched_signature_weights(fixed, bins, *, count, seed):
    """Sample disjoint genes within development-defined bins, keeping fixed weights."""
    fixed = np.asarray(fixed, dtype=float)
    if fixed.ndim != 1 or not np.isfinite(fixed).all() or not np.any(fixed):
        raise ValueError("Provide nonzero finite source-signature weights")
    bins = _ids(bins, len(fixed), "matching bin")
    if type(count) is not int or count < 1:
        raise ValueError("Require a positive number of matched controls")
    source = fixed != 0
    rng = np.random.default_rng(seed)
    result = np.zeros((count, len(fixed)))
    for bin_name in np.unique(bins[source]):
        weights = fixed[source & (bins == bin_name)]
        pool = np.flatnonzero(~source & (bins == bin_name))
        if len(pool) < len(weights):
            raise ValueError(f"Insufficient nonsignature genes in matching bin {bin_name}")
        for row in result:
            row[rng.choice(pool, size=len(weights), replace=False)] = weights
    return result


def score_origins(cultures, predictions):
    """Report primary origin-mean error and secondary within-origin contrast error."""
    names, inverse, counts = _groups(cultures.origins)
    observed = _means(cultures.viability, inverse, counts)
    scores = {}
    for arm, prediction in predictions.items():
        prediction = np.asarray(prediction, dtype=float)
        if prediction.shape != cultures.viability.shape or not np.isfinite(prediction).all():
            raise ValueError("Every arm must predict every scored culture")
        predicted = _means(prediction, inverse, counts)
        absolute = np.abs(predicted - observed)
        centered_error = prediction - cultures.viability - (predicted - observed)[inverse]
        within_rmse = np.sqrt(_means(centered_error**2, inverse, counts))
        scores[arm] = {
            "origin_mae": float(absolute.mean()),
            "mean_within_origin_rmse": float(within_rmse[counts > 1].mean()) if np.any(counts > 1) else None,
            "per_origin": {
                str(name): {
                    "cultures": int(n), "observed": float(truth), "prediction": float(estimate),
                    "absolute_error": float(error), "within_origin_rmse": float(within) if n > 1 else None,
                }
                for name, n, truth, estimate, error, within in zip(
                    names, counts, observed, predicted, absolute, within_rmse, strict=True)
            },
        }
    return scores


def predict_tmz_arms(training, expression, origins, culture_ids, *,
                     penalty, within_weight, strata, shuffle_seed, mgmt=None, signatures=None):
    """Fit development arms and predict on disjoint origins without evaluation outcomes.

    Optional MGMT values are supplied as one train/evaluation vector pair.
    Signature weights are externally fixed in the common, harmonized feature
    order. An origin-balanced linear calibration is fitted on development only.
    """
    evaluation = _matrix(expression)
    evaluation_origins = _ids(origins, len(evaluation), "origin")
    evaluation_cultures = _ids(culture_ids, len(evaluation), "culture")
    if np.unique(evaluation_cultures).size != len(evaluation):
        raise ValueError("Evaluation culture IDs must be unique")
    if (np.intersect1d(training.origins, evaluation_origins).size
            or np.intersect1d(training.culture_ids, evaluation_cultures).size):
        raise ValueError("Training and evaluation must be disjoint by culture and origin")
    if training.expression.shape[1] != evaluation.shape[1]:
        raise ValueError("Training and evaluation must share the frozen feature axis")
    if not np.isfinite([penalty, within_weight]).all() or min(penalty, within_weight) <= 0:
        raise ValueError("Require positive finite penalty and within-origin weight")

    _, inverse, counts = _groups(training.origins)
    n = len(training.origins)
    row_weight = 1.0 / (len(counts) * counts[inverse])
    scale = StandardScaler().fit(training.expression, sample_weight=row_weight)
    x, eval_x = scale.transform(training.expression), scale.transform(evaluation)
    y, origins = training.viability, training.origins
    shuffled = _shuffle_links(origins, strata, shuffle_seed)
    models = {
        "within_origin": (origins, within_weight),
        "origin_mean_only": (origins, 0),
        "shuffled_sister_links": (shuffled, within_weight),
    }
    predictions = {}
    for arm, (links, weight) in models.items():
        intercept, slope = _fit_shared(x, y, origins, links, penalty, weight)
        predictions[arm] = intercept + eval_x @ slope
    predictions["regularized_expression"] = Ridge(alpha=penalty).fit(
        x, y, sample_weight=np.full(n, 1 / n),
    ).predict(eval_x)
    predictions["origin_balanced_expression"] = Ridge(alpha=penalty).fit(
        x, y, sample_weight=row_weight,
    ).predict(eval_x)

    auxiliary = {}
    if mgmt is not None:
        train_mgmt, eval_mgmt = (np.asarray(values, dtype=float) for values in mgmt)
        if (train_mgmt.shape != y.shape or eval_mgmt.shape != (len(evaluation),)
                or not np.isin(train_mgmt, (0, 1)).all() or not np.isin(eval_mgmt, (0, 1)).all()):
            raise ValueError("MGMT comparator needs a known binary status on the common cohort")
        auxiliary["mgmt"] = train_mgmt, eval_mgmt
    for arm, weights in (signatures or {}).items():
        if arm in predictions or arm in auxiliary:
            raise ValueError(f"Signature name collides with model arm {arm}")
        auxiliary[arm] = score_linear(training.expression, weights), score_linear(evaluation, weights)
    mean_y = _means(y, inverse, counts)
    for arm, (train_score, eval_score) in auxiliary.items():
        predictions[arm] = Ridge(alpha=penalty).fit(
            _means(train_score, inverse, counts)[:, None], mean_y,
            sample_weight=np.full(len(counts), 1 / len(counts)),
        ).predict(eval_score[:, None])
    return predictions


def benchmark_tmz(training, evaluation, **options):
    """Compose the prediction and scoring stages for development runs."""
    predictions = predict_tmz_arms(
        training, evaluation.expression, evaluation.origins, evaluation.culture_ids, **options,
    )
    return score_origins(evaluation, predictions)
