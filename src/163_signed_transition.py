#!/usr/bin/env python3
"""Signed, count-aware forecasts of observed tumor-mask transitions.

This module provides the mathematical core for a Track B extension. It models
newly annotated label-1/3 voxels and disappearing label-1/3 voxels as two
distinct outcomes. These are transitions in supplied segmentations, not
verified tumor-cell birth, death, progression, or treatment response.

The legacy count component is a shared-Gamma-frailty Poisson process over at-risk
voxel-days. This is a low-hazard approximation to bounded voxel transitions;
callers must report count saturation when a forecast approaches the number
of voxels at risk. A shared frailty couples the addition and removal counts.
Spatial scores are supplied by a caller and calibrated separately to the
forecast expected counts, so no future mask is used for calibration.
The bounded model below uses ``SharedFrailtyBinomialPrior`` instead. Its
rank-constrained forecast maps count-distribution tails to voxel probabilities;
paired mask draws preserve the shared-frailty addition/removal dependence.
Optional correlated score fields model spatial location uncertainty, while
boundary envelopes show sensitivity to stated contour perturbations. Neither
component identifies biological change separately from annotation error.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterator, Sequence

import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, gaussian_filter
from scipy.optimize import minimize
from scipy.special import expit, gammaln, logsumexp, roots_genlaguerre
from scipy.stats import binom


@dataclass(frozen=True)
class SignedTransitionTargets:
    """Observed additions/removals plus their forecast-origin risk sets."""

    addition_risk: np.ndarray
    removal_risk: np.ndarray
    additions: np.ndarray
    removals: np.ndarray

    def validate(self) -> None:
        arrays = tuple(
            np.asarray(value, dtype=bool)
            for value in (self.addition_risk, self.removal_risk, self.additions, self.removals)
        )
        if not arrays or any(array.shape != arrays[0].shape for array in arrays[1:]):
            raise ValueError("signed transition arrays must share one grid")
        add_risk, remove_risk, additions, removals = arrays
        if np.any(additions & ~add_risk) or np.any(removals & ~remove_risk):
            raise ValueError("transition events must lie inside their forecast-origin risk sets")
        if np.any(additions & removals) or np.any(add_risk & remove_risk):
            raise ValueError("addition and removal risk sets must be disjoint")


def _tumor_mask(labels: np.ndarray, name: str) -> np.ndarray:
    values = np.asarray(labels)
    if values.ndim != 3:
        raise ValueError(f"{name} labels must be three-dimensional")
    if not np.issubdtype(values.dtype, np.integer):
        rounded = np.rint(values)
        if not np.all(np.isfinite(values)) or not np.allclose(values, rounded):
            raise ValueError(f"{name} labels must contain integer codes")
        values = rounded.astype(np.int16)
    unknown = np.setdiff1d(np.unique(values), np.array([0, 1, 2, 3, 4]))
    if unknown.size:
        raise ValueError(f"{name} labels contain undocumented codes: {unknown.tolist()}")
    return (values == 1) | (values == 3)


def make_signed_transition_targets(
    current_labels: np.ndarray,
    future_labels: np.ndarray,
    addition_domain: np.ndarray,
    removal_domain: np.ndarray | None = None,
) -> SignedTransitionTargets:
    """Create future-blind risk sets and signed annotation-transition targets.

    ``addition_domain`` and ``removal_domain`` must be defined using only the
    current scan and earlier information. Additions are evaluated outside the
    current label-1/3 mask. Removals are evaluated only inside that current
    mask. Label-4 cavity voxels at the current visit are excluded from the
    addition risk set. A caller may pass the full current tumor mask as the
    removal domain; omitting it uses that mask directly.
    """
    current = np.asarray(current_labels)
    future = np.asarray(future_labels)
    add_domain = np.asarray(addition_domain, dtype=bool)
    if current.ndim != 3 or future.shape != current.shape or add_domain.shape != current.shape:
        raise ValueError("current labels, future labels, and addition domain must share a 3-D grid")
    if removal_domain is None:
        remove_domain = np.ones(current.shape, dtype=bool)
    else:
        remove_domain = np.asarray(removal_domain, dtype=bool)
        if remove_domain.shape != current.shape:
            raise ValueError("removal domain must share the label grid")

    current_tumor = _tumor_mask(current, "current")
    future_tumor = _tumor_mask(future, "future")
    addition_risk = add_domain & ~current_tumor & (current != 4)
    removal_risk = remove_domain & current_tumor
    additions = future_tumor & ~current_tumor & addition_risk
    removals = current_tumor & ~future_tumor & removal_risk
    result = SignedTransitionTargets(
        addition_risk=addition_risk,
        removal_risk=removal_risk,
        additions=additions,
        removals=removals,
    )
    result.validate()
    return result


@dataclass(frozen=True)
class SignedTransitionTargetBounds:
    """Observed transition masks and boundary-perturbation sensitivity bounds."""

    nominal: SignedTransitionTargets
    addition_certain: np.ndarray
    addition_possible: np.ndarray
    removal_certain: np.ndarray
    removal_possible: np.ndarray
    perturbation_radius_mm: float

    def validate(self) -> None:
        self.nominal.validate()
        arrays = tuple(
            np.asarray(value, dtype=bool)
            for value in (
                self.addition_certain,
                self.addition_possible,
                self.removal_certain,
                self.removal_possible,
            )
        )
        if any(array.shape != self.nominal.additions.shape for array in arrays):
            raise ValueError("transition-bound arrays must share the nominal target grid")
        add_low, add_high, remove_low, remove_high = arrays
        if np.any(add_low & ~self.nominal.additions) or np.any(self.nominal.additions & ~add_high):
            raise ValueError("addition masks must satisfy certain <= nominal <= possible")
        if np.any(remove_low & ~self.nominal.removals) or np.any(self.nominal.removals & ~remove_high):
            raise ValueError("removal masks must satisfy certain <= nominal <= possible")
        if not np.isfinite(self.perturbation_radius_mm) or self.perturbation_radius_mm <= 0:
            raise ValueError("perturbation radius must be finite and positive")


def make_boundary_sensitivity_bounds(
    current_labels: np.ndarray,
    future_labels: np.ndarray,
    addition_domain: np.ndarray,
    voxel_spacing_mm: Sequence[float],
    *,
    perturbation_radius_mm: float,
    removal_domain: np.ndarray | None = None,
) -> SignedTransitionTargetBounds:
    """Bound signed labels under a stated physical boundary perturbation.

    Erosion and dilation define conservative lower and upper transition masks.
    The radius is a sensitivity setting, not an estimate of annotation or
    registration error. This cannot identify biological change from label
    error; it reports how much the endpoint depends on supplied contours.
    Possible masks may include uncertain current-boundary voxels. Addition and
    removal upper masks are componentwise bounds, not one joint scenario.
    """
    spacing = np.asarray(voxel_spacing_mm, dtype=np.float64)
    if spacing.shape != (3,) or np.any(spacing <= 0) or not np.all(np.isfinite(spacing)):
        raise ValueError("voxel_spacing_mm must contain three finite positive values")
    if not np.isfinite(perturbation_radius_mm) or perturbation_radius_mm <= 0:
        raise ValueError("perturbation_radius_mm must be finite and positive")
    current = np.asarray(current_labels)
    future = np.asarray(future_labels)
    if current.ndim != 3 or future.shape != current.shape:
        raise ValueError("current and future labels must share one 3-D grid")
    nominal = make_signed_transition_targets(
        current_labels,
        future_labels,
        addition_domain,
        removal_domain=removal_domain,
    )
    reach = np.ceil(perturbation_radius_mm / spacing).astype(int)
    axes = [np.arange(-n, n + 1) for n in reach]
    offsets = np.stack(np.meshgrid(*axes, indexing="ij"), axis=-1)
    distances_squared = np.sum((offsets * spacing) ** 2, axis=-1)
    structure = distances_squared <= perturbation_radius_mm**2 + 1e-12
    structure[tuple(reach)] = True

    current_tumor = _tumor_mask(current, "current")
    future_tumor = _tumor_mask(future, "future")
    current_inner = binary_erosion(current_tumor, structure=structure, border_value=0)
    current_outer = binary_dilation(current_tumor, structure=structure, border_value=0)
    future_inner = binary_erosion(future_tumor, structure=structure, border_value=0)
    future_outer = binary_dilation(future_tumor, structure=structure, border_value=0)
    add_scope = (np.asarray(addition_domain, dtype=bool) | current_tumor) & (current != 4)
    if removal_domain is None:
        remove_scope = np.ones(current.shape, dtype=bool)
    else:
        remove_scope = np.asarray(removal_domain, dtype=bool) | current_outer
    remove_scope &= current != 4
    bounds = SignedTransitionTargetBounds(
        nominal=nominal,
        addition_certain=future_inner & ~current_outer & nominal.addition_risk,
        addition_possible=future_outer & ~current_inner & add_scope,
        removal_certain=current_inner & ~future_outer & nominal.removal_risk,
        removal_possible=current_outer & ~future_inner & remove_scope,
        perturbation_radius_mm=float(perturbation_radius_mm),
    )
    bounds.validate()
    return bounds


@dataclass(frozen=True)
class LegacySignedCountObservation:
    """One scan interval's event counts and at-risk voxel-day exposures."""

    patient_id: str
    addition_count: int
    removal_count: int
    addition_exposure_voxel_days: float
    removal_exposure_voxel_days: float

    def validate(self) -> None:
        if not self.patient_id:
            raise ValueError("patient_id must be nonempty")
        if self.addition_count < 0 or self.removal_count < 0:
            raise ValueError("transition counts must be nonnegative")
        exposures = (self.addition_exposure_voxel_days, self.removal_exposure_voxel_days)
        if any(not np.isfinite(value) or value < 0 for value in exposures):
            raise ValueError("transition exposures must be finite and nonnegative")
        if self.addition_count and self.addition_exposure_voxel_days <= 0:
            raise ValueError("positive additions require positive addition exposure")
        if self.removal_count and self.removal_exposure_voxel_days <= 0:
            raise ValueError("positive removals require positive removal exposure")


@dataclass(frozen=True)
class SignedCountForecast:
    """Posterior predictive first and second moments for both transition counts."""

    expected_additions: float
    expected_removals: float
    variance_additions: float
    variance_removals: float
    covariance: float
    posterior_frailty_mean: float
    posterior_frailty_variance: float


@dataclass(frozen=True)
class LegacySharedGammaPoissonPrior:
    """Shared patient frailty with separate addition and removal intensities.

    Conditional on patient frailty ``Z_i``, counts are independent Poisson
    variables with means ``Z_i * rate * voxel_day_exposure``. The same
    ``Z_i`` induces positive dependence between additions and removals.
    Intensities and frailty shape are fitted on training patients only.
    """

    addition_rate_per_voxel_day: float
    removal_rate_per_voxel_day: float
    frailty_shape: float
    n_patients: int
    log_likelihood: float

    @classmethod
    def fit(cls, observations: Sequence[LegacySignedCountObservation]) -> "LegacySharedGammaPoissonPrior":
        """Fit the shared-frailty count law by patient-level marginal likelihood."""
        if not observations:
            raise ValueError("at least one training observation is required")
        grouped: dict[str, list[LegacySignedCountObservation]] = {}
        for observation in observations:
            observation.validate()
            grouped.setdefault(observation.patient_id, []).append(observation)
        if len(grouped) < 2:
            raise ValueError("at least two training patients are required to fit frailty")

        add_count = np.asarray(
            [sum(row.addition_count for row in rows) for rows in grouped.values()], dtype=np.float64
        )
        remove_count = np.asarray(
            [sum(row.removal_count for row in rows) for rows in grouped.values()], dtype=np.float64
        )
        add_exposure = np.asarray(
            [sum(row.addition_exposure_voxel_days for row in rows) for rows in grouped.values()],
            dtype=np.float64,
        )
        remove_exposure = np.asarray(
            [sum(row.removal_exposure_voxel_days for row in rows) for rows in grouped.values()],
            dtype=np.float64,
        )
        if add_count.sum() <= 0 or remove_count.sum() <= 0:
            raise ValueError("both transition types need observed events to estimate the joint count law")
        if add_exposure.sum() <= 0 or remove_exposure.sum() <= 0:
            raise ValueError("both transition types need positive training exposure")

        add_rate0 = float(add_count.sum() / add_exposure.sum())
        remove_rate0 = float(remove_count.sum() / remove_exposure.sum())

        def negative_log_likelihood(log_parameters: np.ndarray) -> float:
            log_add_rate, log_remove_rate, log_shape = np.asarray(log_parameters, dtype=np.float64)
            add_rate, remove_rate, shape = np.exp((log_add_rate, log_remove_rate, log_shape))
            add_mean = add_rate * add_exposure
            remove_mean = remove_rate * remove_exposure
            total_count = add_count + remove_count
            total_mean = add_mean + remove_mean
            if np.any((add_count > 0) & (add_mean <= 0)) or np.any((remove_count > 0) & (remove_mean <= 0)):
                return float("inf")

            log_event_terms = np.zeros_like(total_count)
            positive_add = add_count > 0
            positive_remove = remove_count > 0
            log_event_terms[positive_add] += (
                add_count[positive_add] * np.log(add_mean[positive_add])
                - gammaln(add_count[positive_add] + 1.0)
            )
            log_event_terms[positive_remove] += (
                remove_count[positive_remove] * np.log(remove_mean[positive_remove])
                - gammaln(remove_count[positive_remove] + 1.0)
            )
            log_likelihood = (
                gammaln(shape + total_count)
                - gammaln(shape)
                + shape * np.log(shape)
                - (shape + total_count) * np.log(shape + total_mean)
                + log_event_terms
            )
            if not np.all(np.isfinite(log_likelihood)):
                return float("inf")
            return -float(np.sum(log_likelihood))

        initial = np.log([add_rate0, remove_rate0, 1.0])
        rate_bounds = []
        for rate0 in (add_rate0, remove_rate0):
            center = np.log(max(rate0, 1e-12))
            rate_bounds.append((max(-60.0, center - 20.0), min(40.0, center + 20.0)))
        bounds = [*rate_bounds, (-9.0, 14.0)]
        fit = minimize(
            negative_log_likelihood,
            initial,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1e-12},
        )
        if not fit.success or not np.all(np.isfinite(fit.x)):
            raise RuntimeError(f"joint count-prior fit failed: {fit.message}")
        add_rate, remove_rate, shape = np.exp(fit.x)
        return cls(
            addition_rate_per_voxel_day=float(add_rate),
            removal_rate_per_voxel_day=float(remove_rate),
            frailty_shape=float(shape),
            n_patients=len(grouped),
            log_likelihood=-float(fit.fun),
        )

    def forecast(
        self,
        *,
        prior_additions: int,
        prior_removals: int,
        prior_addition_exposure_voxel_days: float,
        prior_removal_exposure_voxel_days: float,
        future_addition_exposure_voxel_days: float,
        future_removal_exposure_voxel_days: float,
    ) -> SignedCountForecast:
        """Update patient frailty from past transitions and forecast both counts."""
        counts = (prior_additions, prior_removals)
        exposures = (
            prior_addition_exposure_voxel_days,
            prior_removal_exposure_voxel_days,
            future_addition_exposure_voxel_days,
            future_removal_exposure_voxel_days,
        )
        if any(value < 0 for value in counts):
            raise ValueError("prior transition counts must be nonnegative")
        if any(not np.isfinite(value) or value < 0 for value in exposures):
            raise ValueError("transition exposures must be finite and nonnegative")
        if (prior_additions and prior_addition_exposure_voxel_days <= 0) or (
            prior_removals and prior_removal_exposure_voxel_days <= 0
        ):
            raise ValueError("positive prior counts require positive prior exposure")

        prior_events = float(prior_additions + prior_removals)
        prior_mean = (
            self.addition_rate_per_voxel_day * prior_addition_exposure_voxel_days
            + self.removal_rate_per_voxel_day * prior_removal_exposure_voxel_days
        )
        posterior_shape = self.frailty_shape + prior_events
        posterior_rate = self.frailty_shape + prior_mean
        frailty_mean = posterior_shape / posterior_rate
        frailty_variance = posterior_shape / (posterior_rate * posterior_rate)

        add_future_mean = self.addition_rate_per_voxel_day * future_addition_exposure_voxel_days
        remove_future_mean = self.removal_rate_per_voxel_day * future_removal_exposure_voxel_days
        expected_additions = add_future_mean * frailty_mean
        expected_removals = remove_future_mean * frailty_mean
        variance_additions = expected_additions + add_future_mean * add_future_mean * frailty_variance
        variance_removals = expected_removals + remove_future_mean * remove_future_mean * frailty_variance
        covariance = add_future_mean * remove_future_mean * frailty_variance
        return SignedCountForecast(
            expected_additions=float(expected_additions),
            expected_removals=float(expected_removals),
            variance_additions=float(variance_additions),
            variance_removals=float(variance_removals),
            covariance=float(covariance),
            posterior_frailty_mean=float(frailty_mean),
            posterior_frailty_variance=float(frailty_variance),
        )


@dataclass(frozen=True)
class BoundedSignedCountObservation:
    """One interval's observed events and finite forecast-origin risk sets."""

    patient_id: str
    addition_count: int
    removal_count: int
    addition_risk_voxels: int
    removal_risk_voxels: int
    interval_days: float

    def validate(self) -> None:
        if not self.patient_id:
            raise ValueError("patient_id must be nonempty")
        counts = (self.addition_count, self.removal_count)
        risks = (self.addition_risk_voxels, self.removal_risk_voxels)
        if any(isinstance(value, (bool, np.bool_)) or int(value) != value for value in counts + risks):
            raise ValueError("counts and risk-set sizes must be integers")
        if any(value < 0 for value in counts + risks):
            raise ValueError("counts and risk-set sizes must be nonnegative")
        if self.addition_count > self.addition_risk_voxels:
            raise ValueError("addition count exceeds its forecast-origin risk set")
        if self.removal_count > self.removal_risk_voxels:
            raise ValueError("removal count exceeds its forecast-origin risk set")
        if not np.isfinite(self.interval_days) or self.interval_days <= 0:
            raise ValueError("interval_days must be finite and positive")


def bounded_observation_from_targets(
    patient_id: str, targets: SignedTransitionTargets, interval_days: float
) -> BoundedSignedCountObservation:
    """Summarize a signed target map without losing its finite denominators."""
    targets.validate()
    observation = BoundedSignedCountObservation(
        patient_id=patient_id,
        addition_count=int(np.count_nonzero(targets.additions)),
        removal_count=int(np.count_nonzero(targets.removals)),
        addition_risk_voxels=int(np.count_nonzero(targets.addition_risk)),
        removal_risk_voxels=int(np.count_nonzero(targets.removal_risk)),
        interval_days=float(interval_days),
    )
    observation.validate()
    return observation


def _frailty_quadrature(shape: float, order: int) -> tuple[np.ndarray, np.ndarray]:
    """Return nodes and normalized weights for Z ~ Gamma(shape, rate=shape)."""
    if not np.isfinite(shape) or shape <= 0 or order < 8:
        raise ValueError("invalid Gamma-frailty quadrature settings")
    nodes, weights = roots_genlaguerre(order, shape - 1.0)
    if np.any(nodes <= 0) or not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(weights)):
        raise FloatingPointError("Gamma-frailty quadrature returned invalid nodes")
    with np.errstate(divide="ignore"):
        log_weights = np.log(weights) - gammaln(shape)
    log_weights -= logsumexp(log_weights)
    return nodes / shape, log_weights


def _conditional_binomial_log_likelihood(
    frailty: np.ndarray,
    observations: Sequence[BoundedSignedCountObservation],
    addition_rate: float,
    removal_rate: float,
) -> np.ndarray:
    """Patient log likelihood conditional on each quadrature frailty node."""
    result = np.zeros(frailty.shape, dtype=np.float64)
    for row in observations:
        for count, risk, rate in (
            (row.addition_count, row.addition_risk_voxels, addition_rate),
            (row.removal_count, row.removal_risk_voxels, removal_rate),
        ):
            if risk == 0:
                continue
            hazard = np.minimum(frailty * rate * row.interval_days, 1e300)
            probability = -np.expm1(-hazard)
            log_probability = np.log(np.maximum(probability, np.finfo(np.float64).tiny))
            log_survival = -hazard
            result += gammaln(risk + 1.0) - gammaln(count + 1.0) - gammaln(risk - count + 1.0)
            if count:
                result += count * log_probability
            if risk > count:
                result += (risk - count) * log_survival
    return result


@dataclass(frozen=True)
class BoundedSignedCountForecast:
    """Joint predictive count law and moments for one future interval."""

    expected_additions: float
    expected_removals: float
    variance_additions: float
    variance_removals: float
    covariance: float
    posterior_frailty_mean: float
    posterior_frailty_variance: float
    frailty_nodes: tuple[float, ...]
    posterior_weights: tuple[float, ...]
    addition_rate_per_voxel_day: float
    removal_rate_per_voxel_day: float
    future_addition_risk_voxels: int
    future_removal_risk_voxels: int
    future_interval_days: float

    def sample_counts(self, *, n_draws: int, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
        """Draw paired counts while preserving shared-frailty dependence."""
        if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 1:
            raise ValueError("n_draws must be positive")
        rng = np.random.default_rng(seed)
        nodes = np.asarray(self.frailty_nodes, dtype=np.float64)
        weights = np.asarray(self.posterior_weights, dtype=np.float64)
        if nodes.ndim != 1 or weights.shape != nodes.shape or nodes.size == 0:
            raise ValueError("forecast lacks a valid frailty quadrature distribution")
        if (
            np.any(nodes <= 0)
            or np.any(weights < 0)
            or not np.all(np.isfinite(nodes))
            or not np.all(np.isfinite(weights))
            or float(weights.sum()) <= 0
        ):
            raise ValueError("forecast frailty nodes or weights are invalid")
        weights = weights / weights.sum()
        selected = rng.choice(nodes.size, size=n_draws, p=weights)
        frailty = nodes[selected]
        add_probability = -np.expm1(
            -frailty * self.addition_rate_per_voxel_day * self.future_interval_days
        )
        remove_probability = -np.expm1(
            -frailty * self.removal_rate_per_voxel_day * self.future_interval_days
        )
        return (
            rng.binomial(self.future_addition_risk_voxels, add_probability),
            rng.binomial(self.future_removal_risk_voxels, remove_probability),
        )

    def log_probability(self, additions: int, removals: int) -> float:
        """Return joint predictive log probability for one observed count pair."""
        observation = BoundedSignedCountObservation(
            patient_id="predictive-score",
            addition_count=additions,
            removal_count=removals,
            addition_risk_voxels=self.future_addition_risk_voxels,
            removal_risk_voxels=self.future_removal_risk_voxels,
            interval_days=self.future_interval_days,
        )
        observation.validate()
        nodes = np.asarray(self.frailty_nodes, dtype=np.float64)
        weights = np.asarray(self.posterior_weights, dtype=np.float64)
        if nodes.ndim != 1 or weights.shape != nodes.shape or nodes.size == 0:
            raise ValueError("forecast lacks a valid frailty quadrature distribution")
        if np.any(nodes <= 0) or np.any(weights < 0) or float(weights.sum()) <= 0:
            raise ValueError("forecast frailty nodes or weights are invalid")
        if not np.all(np.isfinite(nodes)) or not np.all(np.isfinite(weights)):
            raise ValueError("forecast frailty nodes and weights must be finite")
        normalized_weights = weights / weights.sum()
        positive_weights = normalized_weights > 0
        log_weights = np.full(normalized_weights.shape, -np.inf, dtype=np.float64)
        log_weights[positive_weights] = np.log(normalized_weights[positive_weights])
        conditional = _conditional_binomial_log_likelihood(
            nodes,
            [observation],
            self.addition_rate_per_voxel_day,
            self.removal_rate_per_voxel_day,
        )
        return float(logsumexp(log_weights + conditional))


@dataclass(frozen=True)
class SharedFrailtyBinomialPrior:
    """Bounded addition/removal law with shared patient-level Gamma frailty.

    Rates are per voxel-day. Given patient frailty Z and interval length dt,
    each at-risk voxel transitions with probability 1 - exp(-Z * rate * dt).
    Additions and removals are conditionally independent given Z, but share Z.
    Fit integrates Z out by Gauss-Laguerre quadrature. This is an annotation-
    transition model, not a biological cell-birth or cell-death model.
    """

    addition_rate_per_voxel_day: float
    removal_rate_per_voxel_day: float
    frailty_shape: float
    n_patients: int
    log_likelihood: float
    quadrature_order: int = 64
    shape_at_fit_boundary: bool = False

    @classmethod
    def fit(
        cls,
        observations: Sequence[BoundedSignedCountObservation],
        *,
        quadrature_order: int = 64,
    ) -> "SharedFrailtyBinomialPrior":
        """Fit rates and frailty shape using training patients only."""
        if quadrature_order < 8:
            raise ValueError("quadrature_order must be at least 8")
        if not observations:
            raise ValueError("at least one training observation is required")
        grouped: dict[str, list[BoundedSignedCountObservation]] = {}
        for observation in observations:
            observation.validate()
            grouped.setdefault(observation.patient_id, []).append(observation)
        if len(grouped) < 2:
            raise ValueError("at least two training patients are required")

        add_events = sum(row.addition_count for row in observations)
        remove_events = sum(row.removal_count for row in observations)
        add_exposure = sum(row.addition_risk_voxels * row.interval_days for row in observations)
        remove_exposure = sum(row.removal_risk_voxels * row.interval_days for row in observations)
        if add_events <= 0 or remove_events <= 0:
            raise ValueError("both addition and removal events are required to fit the signed count law")
        if (add_events and add_exposure <= 0) or (remove_events and remove_exposure <= 0):
            raise ValueError("positive transition counts require positive risk exposure")

        active_addition = True
        active_removal = True
        initial_rates = (
            max(add_events / add_exposure, 1e-12) if active_addition else 0.0,
            max(remove_events / remove_exposure, 1e-12) if active_removal else 0.0,
        )
        active_rates = [index for index, active in enumerate((active_addition, active_removal)) if active]
        initial = [np.log(initial_rates[index]) for index in active_rates] + [0.0]
        bounds: list[tuple[float, float]] = []
        for index in active_rates:
            center = float(np.log(initial_rates[index]))
            bounds.append((max(-60.0, center - 20.0), min(40.0, center + 20.0)))
        bounds.append((np.log(0.1), np.log(148.0)))

        def objective(parameters: np.ndarray) -> float:
            rates = [0.0, 0.0]
            for index, log_rate in zip(active_rates, parameters[:-1]):
                rates[index] = float(np.exp(log_rate))
            shape = float(np.exp(parameters[-1]))
            try:
                frailty, log_weights = _frailty_quadrature(shape, quadrature_order)
                total_log_likelihood = 0.0
                for rows in grouped.values():
                    conditional = _conditional_binomial_log_likelihood(
                        frailty, rows, rates[0], rates[1]
                    )
                    total_log_likelihood += float(logsumexp(log_weights + conditional))
            except (FloatingPointError, ValueError, OverflowError):
                return float("inf")
            return -total_log_likelihood if np.isfinite(total_log_likelihood) else float("inf")

        fit = minimize(
            objective,
            np.asarray(initial, dtype=np.float64),
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 2000, "ftol": 1e-10},
        )
        if not fit.success or not np.all(np.isfinite(fit.x)):
            raise RuntimeError(f"bounded joint count-prior fit failed: {fit.message}")
        rates = [0.0, 0.0]
        for index, log_rate in zip(active_rates, fit.x[:-1]):
            rates[index] = float(np.exp(log_rate))
        shape = float(np.exp(fit.x[-1]))
        on_boundary = bool(
            np.isclose(fit.x[-1], bounds[-1][0], atol=1e-4)
            or np.isclose(fit.x[-1], bounds[-1][1], atol=1e-4)
        )
        return cls(
            addition_rate_per_voxel_day=rates[0],
            removal_rate_per_voxel_day=rates[1],
            frailty_shape=shape,
            n_patients=len(grouped),
            log_likelihood=-float(fit.fun),
            quadrature_order=quadrature_order,
            shape_at_fit_boundary=on_boundary,
        )

    def forecast(
        self,
        *,
        prior_observations: Sequence[BoundedSignedCountObservation],
        future_addition_risk_voxels: int,
        future_removal_risk_voxels: int,
        future_interval_days: float,
        future_addition_rate_multiplier: float = 1.0,
        future_removal_rate_multiplier: float = 1.0,
    ) -> BoundedSignedCountForecast:
        """Return bounded predictive moments, updating frailty from past data.

        Optional direction-specific multipliers recalibrate only the future
        hazards. The patient frailty posterior is still updated with the rates
        estimated by the training fold, so the multiplier is a held-out
        predictive calibration layer rather than a change to the history model.
        """
        for row in prior_observations:
            row.validate()
        patient_ids = {row.patient_id for row in prior_observations}
        if len(patient_ids) > 1:
            raise ValueError("frailty update must use one patient's prior intervals")
        future_risks = (future_addition_risk_voxels, future_removal_risk_voxels)
        if any(isinstance(value, (bool, np.bool_)) or int(value) != value or value < 0 for value in future_risks):
            raise ValueError("future risk-set sizes must be nonnegative integers")
        if not np.isfinite(future_interval_days) or future_interval_days <= 0:
            raise ValueError("future_interval_days must be finite and positive")
        rate_multipliers = (
            future_addition_rate_multiplier,
            future_removal_rate_multiplier,
        )
        if any(not np.isfinite(value) or value < 0 for value in rate_multipliers):
            raise ValueError("future rate multipliers must be finite and nonnegative")

        frailty, log_weights = _frailty_quadrature(self.frailty_shape, self.quadrature_order)
        if prior_observations:
            log_weights = log_weights + _conditional_binomial_log_likelihood(
                frailty,
                prior_observations,
                self.addition_rate_per_voxel_day,
                self.removal_rate_per_voxel_day,
            )
            log_weights -= logsumexp(log_weights)
        weights = np.exp(log_weights)

        future_addition_rate = self.addition_rate_per_voxel_day * future_addition_rate_multiplier
        future_removal_rate = self.removal_rate_per_voxel_day * future_removal_rate_multiplier
        add_hazard = frailty * future_addition_rate * future_interval_days
        remove_hazard = frailty * future_removal_rate * future_interval_days
        add_probability = -np.expm1(-add_hazard)
        remove_probability = -np.expm1(-remove_hazard)
        add_conditional_mean = future_addition_risk_voxels * add_probability
        remove_conditional_mean = future_removal_risk_voxels * remove_probability
        mean_add = float(np.dot(weights, add_conditional_mean))
        mean_remove = float(np.dot(weights, remove_conditional_mean))
        variance_add = float(
            np.dot(weights, future_addition_risk_voxels * add_probability * (1 - add_probability)
                   + (add_conditional_mean - mean_add) ** 2)
        )
        variance_remove = float(
            np.dot(weights, future_removal_risk_voxels * remove_probability * (1 - remove_probability)
                   + (remove_conditional_mean - mean_remove) ** 2)
        )
        covariance = float(
            np.dot(weights, (add_conditional_mean - mean_add) * (remove_conditional_mean - mean_remove))
        )
        frailty_mean = float(np.dot(weights, frailty))
        frailty_variance = float(np.dot(weights, (frailty - frailty_mean) ** 2))
        return BoundedSignedCountForecast(
            expected_additions=mean_add,
            expected_removals=mean_remove,
            variance_additions=variance_add,
            variance_removals=variance_remove,
            covariance=covariance,
            posterior_frailty_mean=frailty_mean,
            posterior_frailty_variance=frailty_variance,
            frailty_nodes=tuple(float(value) for value in frailty),
            posterior_weights=tuple(float(value) for value in weights),
            addition_rate_per_voxel_day=future_addition_rate,
            removal_rate_per_voxel_day=future_removal_rate,
            future_addition_risk_voxels=int(future_addition_risk_voxels),
            future_removal_risk_voxels=int(future_removal_risk_voxels),
            future_interval_days=float(future_interval_days),
        )


@dataclass(frozen=True)
class SignedSpatialForecast:
    """Calibrated marginal probabilities for the two disjoint transitions."""

    addition_probability: np.ndarray
    removal_probability: np.ndarray
    requested_addition_count: float
    requested_removal_count: float
    used_addition_count: float
    used_removal_count: float
    addition_count_clipped: bool
    removal_count_clipped: bool
    addition_count_variance: float = 0.0
    removal_count_variance: float = 0.0
    addition_removal_count_covariance: float = 0.0


def _calibrate_mass(scores: np.ndarray, expected_count: float) -> np.ndarray:
    values = np.asarray(scores, dtype=np.float64).reshape(-1)
    if not np.all(np.isfinite(values)):
        raise ValueError("spatial scores must be finite")
    if not np.isfinite(expected_count) or expected_count < 0:
        raise ValueError("expected count must be finite and nonnegative")
    if values.size == 0:
        if expected_count > 0:
            raise ValueError("positive expected count cannot fit an empty risk set")
        return np.empty(0, dtype=np.float64)
    if expected_count > values.size + 1e-8:
        raise ValueError("expected count exceeds finite risk-set size")
    target = min(float(expected_count), float(values.size))
    if target <= 1e-12:
        return np.zeros(values.size, dtype=np.float64)
    if target >= values.size - 1e-12:
        return np.ones(values.size, dtype=np.float64)
    centered = values - float(np.median(values))
    bound = max(40.0, float(np.max(np.abs(centered))) + 40.0)
    low, high = -bound, bound
    for _ in range(100):
        intercept = 0.5 * (low + high)
        if float(np.sum(expit(intercept + centered))) < target:
            low = intercept
        else:
            high = intercept
    return expit(0.5 * (low + high) + centered)


def calibrate_signed_spatial_forecast(
    addition_scores: np.ndarray,
    removal_scores: np.ndarray,
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    expected_additions: float,
    expected_removals: float,
) -> SignedSpatialForecast:
    """Expected-count-only calibration baseline for two future-blind score maps.

    ``addition_scores`` and ``removal_scores`` are full-grid score maps. The
    addition map is calibrated only outside the current tumor mask; the
    removal map is calibrated only inside it. Returned maps are zero outside
    their respective forecast-origin risk sets. This path does not preserve
    count variances or addition/removal covariance. Use
    ``rank_constrained_signed_forecast`` for the full bounded count law.
    """
    add_score = np.asarray(addition_scores, dtype=np.float64)
    remove_score = np.asarray(removal_scores, dtype=np.float64)
    add_risk = np.asarray(addition_risk, dtype=bool)
    remove_risk = np.asarray(removal_risk, dtype=bool)
    if add_score.ndim != 3 or any(
        array.shape != add_score.shape for array in (remove_score, add_risk, remove_risk)
    ):
        raise ValueError("signed scores and risk masks must share one 3-D grid")
    if np.any(add_risk & remove_risk):
        raise ValueError("addition and removal risk sets must be disjoint")

    n_add = int(np.count_nonzero(add_risk))
    n_remove = int(np.count_nonzero(remove_risk))
    if any(not np.isfinite(value) or value < 0 for value in (expected_additions, expected_removals)):
        raise ValueError("expected transition counts must be finite and nonnegative")
    if expected_additions > n_add + 1e-8 or expected_removals > n_remove + 1e-8:
        raise ValueError("expected transition count exceeds its finite risk set")
    add_used = min(float(expected_additions), float(n_add))
    remove_used = min(float(expected_removals), float(n_remove))
    add_probability = np.zeros(add_score.shape, dtype=np.float64)
    remove_probability = np.zeros(remove_score.shape, dtype=np.float64)
    add_probability[add_risk] = _calibrate_mass(add_score[add_risk], add_used)
    remove_probability[remove_risk] = _calibrate_mass(remove_score[remove_risk], remove_used)
    return SignedSpatialForecast(
        addition_probability=add_probability,
        removal_probability=remove_probability,
        requested_addition_count=float(expected_additions),
        requested_removal_count=float(expected_removals),
        used_addition_count=add_used,
        used_removal_count=remove_used,
        addition_count_clipped=not np.isclose(add_used, expected_additions),
        removal_count_clipped=not np.isclose(remove_used, expected_removals),
    )


def _count_tail_probabilities(
    risk_size: int,
    rate: float,
    interval_days: float,
    frailty_nodes: Sequence[float],
    posterior_weights: Sequence[float],
) -> np.ndarray:
    """Compute P(K >= rank) for a frailty-mixture binomial count law."""
    if risk_size < 0:
        raise ValueError("risk_size must be nonnegative")
    if risk_size == 0:
        return np.empty(0, dtype=np.float64)
    nodes = np.asarray(frailty_nodes, dtype=np.float64)
    weights = np.asarray(posterior_weights, dtype=np.float64)
    if nodes.ndim != 1 or weights.shape != nodes.shape or nodes.size == 0:
        raise ValueError("count forecast lacks frailty quadrature nodes and weights")
    if (
        np.any(nodes <= 0)
        or np.any(weights < 0)
        or not np.all(np.isfinite(nodes))
        or not np.all(np.isfinite(weights))
        or float(weights.sum()) <= 0
        or not np.isfinite(rate)
        or rate < 0
        or not np.isfinite(interval_days)
        or interval_days <= 0
    ):
        raise ValueError("invalid frailty quadrature distribution")
    weights = weights / weights.sum()
    ranks = np.arange(1, risk_size + 1, dtype=np.int64)
    tails = np.zeros(risk_size, dtype=np.float64)
    for node, weight in zip(nodes, weights):
        probability = float(-np.expm1(-node * rate * interval_days))
        tails += float(weight) * binom.sf(ranks - 1, risk_size, probability)
    tails = np.minimum.accumulate(np.clip(tails, 0.0, 1.0))
    return tails


def rank_constrained_signed_forecast(
    addition_scores: np.ndarray,
    removal_scores: np.ndarray,
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    count_forecast: BoundedSignedCountForecast,
) -> SignedSpatialForecast:
    """Map joint bounded count forecasts to spatial marginal probabilities.

    Each score map ranks voxels inside its own forecast-origin risk set. A
    voxel at rank r receives P(K >= r), where K follows the corresponding
    marginal of the shared-frailty binomial forecast. Exact score ties receive
    the mean probability across their tied ranks. Thus map probability mass
    equals the forecast expected count, while ``count_forecast`` retains the
    addition/removal count covariance and supports paired count sampling.
    """
    add_score = np.asarray(addition_scores, dtype=np.float64)
    remove_score = np.asarray(removal_scores, dtype=np.float64)
    add_risk = np.asarray(addition_risk, dtype=bool)
    remove_risk = np.asarray(removal_risk, dtype=bool)
    if add_score.ndim != 3 or any(
        array.shape != add_score.shape for array in (remove_score, add_risk, remove_risk)
    ):
        raise ValueError("signed scores and risk masks must share one 3-D grid")
    if not np.all(np.isfinite(add_score)) or not np.all(np.isfinite(remove_score)):
        raise ValueError("signed spatial scores must be finite")
    if np.any(add_risk & remove_risk):
        raise ValueError("addition and removal risk sets must be disjoint")
    add_indices = np.flatnonzero(add_risk)
    remove_indices = np.flatnonzero(remove_risk)
    if add_indices.size != count_forecast.future_addition_risk_voxels:
        raise ValueError("addition risk set differs from count-forecast denominator")
    if remove_indices.size != count_forecast.future_removal_risk_voxels:
        raise ValueError("removal risk set differs from count-forecast denominator")

    add_tail = _count_tail_probabilities(
        int(add_indices.size),
        count_forecast.addition_rate_per_voxel_day,
        count_forecast.future_interval_days,
        count_forecast.frailty_nodes,
        count_forecast.posterior_weights,
    )
    remove_tail = _count_tail_probabilities(
        int(remove_indices.size),
        count_forecast.removal_rate_per_voxel_day,
        count_forecast.future_interval_days,
        count_forecast.frailty_nodes,
        count_forecast.posterior_weights,
    )
    if not np.isclose(np.sum(add_tail), count_forecast.expected_additions, rtol=1e-6, atol=1e-6):
        raise ValueError("addition count moments disagree with count-distribution tails")
    if not np.isclose(np.sum(remove_tail), count_forecast.expected_removals, rtol=1e-6, atol=1e-6):
        raise ValueError("removal count moments disagree with count-distribution tails")

    def rank_map(indices: np.ndarray, score_map: np.ndarray, tails: np.ndarray) -> np.ndarray:
        probabilities = np.zeros(score_map.size, dtype=np.float64)
        if indices.size:
            order = np.argsort(-score_map.reshape(-1)[indices], kind="mergesort")
            sorted_scores = score_map.reshape(-1)[indices[order]]
            sorted_probabilities = tails.copy()
            starts = np.r_[0, np.flatnonzero(sorted_scores[1:] != sorted_scores[:-1]) + 1]
            ends = np.r_[starts[1:], sorted_scores.size]
            for start, end in zip(starts, ends):
                sorted_probabilities[start:end] = float(np.mean(tails[start:end]))
            probabilities[indices[order]] = sorted_probabilities
        return probabilities.reshape(score_map.shape)

    add_probability = rank_map(add_indices, add_score, add_tail)
    remove_probability = rank_map(remove_indices, remove_score, remove_tail)
    return SignedSpatialForecast(
        addition_probability=add_probability,
        removal_probability=remove_probability,
        requested_addition_count=count_forecast.expected_additions,
        requested_removal_count=count_forecast.expected_removals,
        used_addition_count=count_forecast.expected_additions,
        used_removal_count=count_forecast.expected_removals,
        addition_count_clipped=False,
        removal_count_clipped=False,
        addition_count_variance=count_forecast.variance_additions,
        removal_count_variance=count_forecast.variance_removals,
        addition_removal_count_covariance=count_forecast.covariance,
    )


@dataclass(frozen=True)
class CoupledSignedForecast:
    """Spatial marginal maps paired with their shared-frailty count law."""

    spatial: SignedSpatialForecast
    counts: BoundedSignedCountForecast


def couple_signed_spatial_and_count_forecast(
    addition_scores: np.ndarray,
    removal_scores: np.ndarray,
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    count_model: SharedFrailtyBinomialPrior,
    *,
    prior_observations: Sequence[BoundedSignedCountObservation],
    future_interval_days: float,
) -> CoupledSignedForecast:
    """Apply one count law to both spatial heads and retain count dependence."""
    add_risk = np.asarray(addition_risk, dtype=bool)
    remove_risk = np.asarray(removal_risk, dtype=bool)
    if add_risk.ndim != 3 or remove_risk.shape != add_risk.shape:
        raise ValueError("addition and removal risks must share one 3-D grid")
    if np.any(add_risk & remove_risk):
        raise ValueError("addition and removal risk sets must be disjoint")
    count_forecast = count_model.forecast(
        prior_observations=prior_observations,
        future_addition_risk_voxels=int(np.count_nonzero(add_risk)),
        future_removal_risk_voxels=int(np.count_nonzero(remove_risk)),
        future_interval_days=future_interval_days,
    )
    spatial_forecast = rank_constrained_signed_forecast(
        addition_scores,
        removal_scores,
        add_risk,
        remove_risk,
        count_forecast,
    )
    return CoupledSignedForecast(spatial=spatial_forecast, counts=count_forecast)


def iter_ranked_signed_mask_draws(
    addition_scores: np.ndarray,
    removal_scores: np.ndarray,
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    count_forecast: BoundedSignedCountForecast,
    *,
    n_draws: int,
    seed: int = 0,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield paired binary masks with sampled counts and spatial ranks.

    The count law supplies the number of additions and removals. Scores choose
    locations by rank. Exact score ties are broken randomly per draw. The
    iterator stores only two full-grid masks at a time.
    """
    add_score = np.asarray(addition_scores, dtype=np.float64)
    remove_score = np.asarray(removal_scores, dtype=np.float64)
    add_risk = np.asarray(addition_risk, dtype=bool)
    remove_risk = np.asarray(removal_risk, dtype=bool)
    if add_score.ndim != 3 or any(
        array.shape != add_score.shape for array in (remove_score, add_risk, remove_risk)
    ):
        raise ValueError("signed scores and risk masks must share one 3-D grid")
    if not np.all(np.isfinite(add_score)) or not np.all(np.isfinite(remove_score)):
        raise ValueError("signed spatial scores must be finite")
    if np.any(add_risk & remove_risk):
        raise ValueError("addition and removal risk sets must be disjoint")
    add_indices = np.flatnonzero(add_risk)
    remove_indices = np.flatnonzero(remove_risk)
    if add_indices.size != count_forecast.future_addition_risk_voxels:
        raise ValueError("addition risk set differs from count-forecast denominator")
    if remove_indices.size != count_forecast.future_removal_risk_voxels:
        raise ValueError("removal risk set differs from count-forecast denominator")
    add_counts, remove_counts = count_forecast.sample_counts(n_draws=n_draws, seed=seed)
    rng = np.random.default_rng(seed + 1)
    flat_add_score = add_score.reshape(-1)
    flat_remove_score = remove_score.reshape(-1)
    for add_count, remove_count in zip(add_counts, remove_counts):
        add_order = np.lexsort((rng.random(add_indices.size), -flat_add_score[add_indices]))
        remove_order = np.lexsort((rng.random(remove_indices.size), -flat_remove_score[remove_indices]))
        add_mask = np.zeros(add_score.size, dtype=bool)
        remove_mask = np.zeros(remove_score.size, dtype=bool)
        add_mask[add_indices[add_order[: int(add_count)]]] = True
        remove_mask[remove_indices[remove_order[: int(remove_count)]]] = True
        yield add_mask.reshape(add_score.shape), remove_mask.reshape(remove_score.shape)


@dataclass(frozen=True)
class CorrelatedSignedSpatialForecast:
    """Monte Carlo marginal maps from spatially correlated score fields."""

    addition_probability: np.ndarray
    removal_probability: np.ndarray
    addition_monte_carlo_se: np.ndarray
    removal_monte_carlo_se: np.ndarray
    n_draws: int
    score_noise_sd: float
    correlation_length_mm: float


def iter_correlated_signed_mask_draws(
    addition_scores: np.ndarray,
    removal_scores: np.ndarray,
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    count_forecast: BoundedSignedCountForecast,
    voxel_spacing_mm: Sequence[float],
    *,
    score_noise_sd: float,
    correlation_length_mm: float,
    n_draws: int,
    seed: int = 0,
) -> Iterator[tuple[np.ndarray, np.ndarray]]:
    """Yield exact-count masks with spatially correlated location variation.

    A Gaussian-smoothed score perturbation models residual spatial dependence.
    Its scale and correlation length are sensitivity parameters, not measured
    segmentation uncertainty. Select them inside training folds only.
    """
    add_score = np.asarray(addition_scores, dtype=np.float64)
    remove_score = np.asarray(removal_scores, dtype=np.float64)
    add_risk = np.asarray(addition_risk, dtype=bool)
    remove_risk = np.asarray(removal_risk, dtype=bool)
    spacing = np.asarray(voxel_spacing_mm, dtype=np.float64)
    if add_score.ndim != 3 or any(
        array.shape != add_score.shape for array in (remove_score, add_risk, remove_risk)
    ):
        raise ValueError("signed scores and risk masks must share one 3-D grid")
    if spacing.shape != (3,) or np.any(spacing <= 0) or not np.all(np.isfinite(spacing)):
        raise ValueError("voxel_spacing_mm must contain three finite positive values")
    if not np.isfinite(score_noise_sd) or score_noise_sd < 0:
        raise ValueError("score_noise_sd must be finite and nonnegative")
    if not np.isfinite(correlation_length_mm) or correlation_length_mm <= 0:
        raise ValueError("correlation_length_mm must be finite and positive")
    if isinstance(n_draws, (bool, np.bool_)) or int(n_draws) != n_draws or n_draws < 1:
        raise ValueError("n_draws must be a positive integer")
    if not np.all(np.isfinite(add_score)) or not np.all(np.isfinite(remove_score)):
        raise ValueError("signed spatial scores must be finite")
    if np.any(add_risk & remove_risk):
        raise ValueError("addition and removal risk sets must be disjoint")
    add_indices = np.flatnonzero(add_risk)
    remove_indices = np.flatnonzero(remove_risk)
    if add_indices.size != count_forecast.future_addition_risk_voxels:
        raise ValueError("addition risk set differs from count-forecast denominator")
    if remove_indices.size != count_forecast.future_removal_risk_voxels:
        raise ValueError("removal risk set differs from count-forecast denominator")

    add_counts, remove_counts = count_forecast.sample_counts(n_draws=n_draws, seed=seed)
    rng = np.random.default_rng(seed + 1)
    sigma = tuple(float(correlation_length_mm / value) for value in spacing)
    def perturbed_order(
        score_map: np.ndarray, risk: np.ndarray, indices: np.ndarray
    ) -> np.ndarray:
        if indices.size == 0 or score_noise_sd == 0:
            noise = np.zeros(score_map.shape, dtype=np.float32)
        else:
            noise = gaussian_filter(
                rng.standard_normal(score_map.shape, dtype=np.float32), sigma=sigma, mode="reflect"
            )
            local = noise[risk].astype(np.float64)
            local_sd = float(np.std(local))
            if local_sd > 1e-12:
                noise = ((noise - float(np.mean(local))) / local_sd).astype(np.float32)
            else:
                noise.fill(0.0)
        perturbed = score_map.reshape(-1)[indices] + score_noise_sd * noise.reshape(-1)[indices]
        return indices[np.lexsort((rng.random(indices.size), -perturbed))]

    for add_count, remove_count in zip(add_counts, remove_counts):
        add_order = perturbed_order(add_score, add_risk, add_indices)
        remove_order = perturbed_order(remove_score, remove_risk, remove_indices)
        add_mask = np.zeros(add_score.size, dtype=bool)
        remove_mask = np.zeros(remove_score.size, dtype=bool)
        add_mask[add_order[: int(add_count)]] = True
        remove_mask[remove_order[: int(remove_count)]] = True
        yield add_mask.reshape(add_score.shape), remove_mask.reshape(remove_score.shape)


def monte_carlo_correlated_signed_forecast(
    addition_scores: np.ndarray,
    removal_scores: np.ndarray,
    addition_risk: np.ndarray,
    removal_risk: np.ndarray,
    count_forecast: BoundedSignedCountForecast,
    voxel_spacing_mm: Sequence[float],
    *,
    score_noise_sd: float,
    correlation_length_mm: float,
    n_draws: int,
    seed: int = 0,
) -> CorrelatedSignedSpatialForecast:
    """Estimate voxel marginals from paired correlated, exact-count mask draws."""
    add_sum = np.zeros(np.asarray(addition_risk).shape, dtype=np.float64)
    remove_sum = np.zeros(np.asarray(removal_risk).shape, dtype=np.float64)
    for add_mask, remove_mask in iter_correlated_signed_mask_draws(
        addition_scores,
        removal_scores,
        addition_risk,
        removal_risk,
        count_forecast,
        voxel_spacing_mm,
        score_noise_sd=score_noise_sd,
        correlation_length_mm=correlation_length_mm,
        n_draws=n_draws,
        seed=seed,
    ):
        add_sum += add_mask
        remove_sum += remove_mask
    add_probability = add_sum / n_draws
    remove_probability = remove_sum / n_draws
    add_se = np.sqrt(add_probability * (1.0 - add_probability) / n_draws)
    remove_se = np.sqrt(remove_probability * (1.0 - remove_probability) / n_draws)
    return CorrelatedSignedSpatialForecast(
        addition_probability=add_probability.astype(np.float32),
        removal_probability=remove_probability.astype(np.float32),
        addition_monte_carlo_se=add_se.astype(np.float32),
        removal_monte_carlo_se=remove_se.astype(np.float32),
        n_draws=int(n_draws),
        score_noise_sd=float(score_noise_sd),
        correlation_length_mm=float(correlation_length_mm),
    )
