"""Generate assumption-only molecular and MRI signals for virtual-patient studies.

This module does not read patient data or estimate a molecular-to-resistance link.
The A5-derived ``q`` value is a scenario setting for a hypothetical report of a
Core-like/Periphery-like state. One call represents one virtual world, so ``q``
is held fixed across its ``n`` patients. It is not patient-level assay validation.
"""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral

import numpy as np
from scipy.special import ndtr, ndtri


# Envelope of the assumed f0 ranges in the calibration families.
F0_REFERENCE_BOUNDS = (0.001 * float(np.exp(-0.7)), 0.9)


@dataclass(frozen=True)
class SignalObservations:
    """Signals a policy may receive. This object contains no latent truth."""

    molecular_report: np.ndarray
    mri_score: np.ndarray


@dataclass(frozen=True)
class SignalTruth:
    """Latent variables for evaluation and generator tests only."""

    resistance_z: np.ndarray
    molecular_state_z: np.ndarray
    molecular_state: np.ndarray
    mri_latent_z: np.ndarray


@dataclass(frozen=True)
class SignalModelParameters:
    """Assumed factor coefficients known to a policy, never fitted on test outcomes.

    ``a`` and ``b`` are loadings on the resistance score, and ``d`` is the
    residual-factor coefficient. They are not generally the empirical
    correlations in a virtual cohort, especially when its resistance scores
    come from a non-normal distribution of ``f0`` values.
    """

    q: float
    a: float
    b: float
    d: float

    def __post_init__(self):
        object.__setattr__(self, "q", _bounded_probability(self.q, "q"))
        object.__setattr__(self, "a", _bounded_coefficient(self.a, "a"))
        object.__setattr__(self, "b", _bounded_coefficient(self.b, "b"))
        object.__setattr__(self, "d", _bounded_coefficient(self.d, "d"))


def standard_normal_marginal_signature(parameters: SignalModelParameters) -> tuple[float, float]:
    """Return (q, rho), which determines paired-signal marginals under normal R.

    This identity assumes a standard-normal latent resistance score and independent
    standard-normal factor errors. After resistance is marginalized, the two latent
    Gaussian scores are bivariate normal with correlation

        rho = a*b + sqrt(1-a^2)*sqrt(1-b^2)*d.

    The observed binary molecular report and continuous MRI score then depend on
    ``q`` and ``rho``, not on ``a``, ``b``, and ``d`` separately. This is a structural
    identifiability diagnostic, not an estimator or an empirical correlation.
    """
    if not all(hasattr(parameters, name) for name in ("q", "a", "b", "d")):
        raise TypeError("parameters must provide q, a, b, and d")
    parameters = SignalModelParameters(
        q=parameters.q, a=parameters.a, b=parameters.b, d=parameters.d
    )
    q, a, b, d = parameters.q, parameters.a, parameters.b, parameters.d
    rho = a * b + np.sqrt(max(0.0, 1.0 - a * a)) * np.sqrt(
        max(0.0, 1.0 - b * b)
    ) * d
    return float(q), float(np.clip(rho, -1.0, 1.0))


def _bounded_probability(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite number in [0, 1]")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite number in [0, 1]") from exc
    if not np.isfinite(result) or not 0.0 <= result <= 1.0:
        raise ValueError(f"{name} must be a finite number in [0, 1]")
    return result


def _bounded_coefficient(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite coefficient in [-1, 1]")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite coefficient in [-1, 1]") from exc
    if not np.isfinite(result) or not -1.0 <= result <= 1.0:
        raise ValueError(f"{name} must be a finite coefficient in [-1, 1]")
    return result


def _template_location_scale(template: np.ndarray) -> tuple[float, float]:
    """Return the training-only MRI template's mean and sample SD."""
    marginal = np.asarray(template, dtype=float).reshape(-1)
    if marginal.size < 2 or not np.isfinite(marginal).all():
        raise ValueError("mri_template must contain at least two finite values")
    location = float(np.mean(marginal))
    scale = float(np.std(marginal, ddof=1))
    if not np.isfinite(location) or not np.isfinite(scale) or scale <= 0.0:
        raise ValueError("mri_template must have a finite, positive sample SD")
    return location, scale


def _calibrate_mri_score(values: np.ndarray, template: np.ndarray | None) -> np.ndarray:
    """Calibrate latent MRI scores to a history-only template's mean and sample SD."""
    if template is None:
        return values.copy()
    location, scale = _template_location_scale(template)
    return location + scale * np.asarray(values, dtype=float)


def f0_to_resistance_score(f0, *, reference_bounds=F0_REFERENCE_BOUNDS) -> np.ndarray:
    """Map simulator f0 through one fixed log-uniform reference CDF to a Gaussian score."""
    values = np.asarray(f0, dtype=float)
    try:
        low, high = (float(reference_bounds[0]), float(reference_bounds[1]))
    except (TypeError, ValueError, IndexError, OverflowError) as exc:
        raise ValueError("reference_bounds must contain two finite values") from exc
    if (
        not np.isfinite(values).all()
        or np.any(values <= 0.0)
        or np.any(values >= 1.0)
        or not np.isfinite([low, high]).all()
        or low <= 0.0
        or high <= low
    ):
        raise ValueError("f0 and reference bounds must be finite, with 0 < f0 < 1 and 0 < low < high")
    probability = (np.log(values) - np.log(low)) / (np.log(high) - np.log(low))
    return ndtri(np.clip(probability, 1e-8, 1.0 - 1e-8))


def _uncalibrate_mri_score(values, template) -> np.ndarray:
    """Invert the history-only affine location/scale calibration exactly."""
    scores = np.asarray(values, dtype=float)
    if not np.isfinite(scores).all():
        raise ValueError("MRI scores must be finite")
    location, scale = _template_location_scale(template)
    return (scores - location) / scale


def generate_integrated_signals(
    n: int,
    *,
    q: float,
    a: float,
    b: float,
    d: float,
    seed: int,
    f0: np.ndarray | None = None,
    mri_template: np.ndarray | None = None,
) -> tuple[SignalObservations, SignalTruth]:
    """Draw paired synthetic signals and separate policy inputs from latent truth.

    If ``f0`` is supplied, one fixed log-uniform reference CDF maps simulator
    resistance fraction to ``Z_resistance``. If omitted, the older independent-normal
    generator remains available for isolated component tests.

    Independent standard-normal variables define:

        Z_bio = a * Z_resistance + sqrt(1-a^2) * E_bio
        Z_mri = b * Z_resistance + sqrt(1-b^2) * (d * E_bio + sqrt(1-d^2) * E_mri)

    ``a`` sets the latent Gaussian association that induces the link between the
    Core-like/Periphery-like molecular state and resistance after thresholding.
    ``b`` links the latent MRI score to resistance. ``d`` sets residual dependence.
    A binary molecular report gives the correct sign of ``Z_bio`` with probability
    ``q``. The report and MRI score are the only fields returned in observations.

    The parameters are scenario assumptions. ``mri_template`` must be built from
    training-only or otherwise permitted history data. Its mean and sample standard
    deviation calibrate the synthetic MRI-score units; the template's full empirical
    distribution is not used.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError("seed must be a non-negative integer")

    q = _bounded_probability(q, "q")
    a = _bounded_coefficient(a, "a")
    b = _bounded_coefficient(b, "b")
    d = _bounded_coefficient(d, "d")

    rng = np.random.default_rng(int(seed))
    if f0 is None:
        resistance_z = rng.standard_normal(int(n))
    else:
        f0_values = np.asarray(f0, dtype=float)
        if f0_values.shape != (int(n),):
            raise ValueError(f"f0 must have shape ({int(n)},)")
        resistance_z = f0_to_resistance_score(f0_values)
    e_bio = rng.standard_normal(int(n))
    e_mri = rng.standard_normal(int(n))

    molecular_state_z = a * resistance_z + np.sqrt(max(0.0, 1.0 - a * a)) * e_bio
    mri_latent_z = b * resistance_z + np.sqrt(max(0.0, 1.0 - b * b)) * (
        d * e_bio + np.sqrt(max(0.0, 1.0 - d * d)) * e_mri
    )
    molecular_state = molecular_state_z >= 0.0
    molecular_report = molecular_state ^ (rng.random(int(n)) >= q)
    mri_score = _calibrate_mri_score(mri_latent_z, mri_template)

    observations = SignalObservations(
        molecular_report=np.asarray(molecular_report, dtype=bool),
        mri_score=np.asarray(mri_score, dtype=float),
    )
    truth = SignalTruth(
        resistance_z=np.asarray(resistance_z, dtype=float),
        molecular_state_z=np.asarray(molecular_state_z, dtype=float),
        molecular_state=np.asarray(molecular_state, dtype=bool),
        mri_latent_z=np.asarray(mri_latent_z, dtype=float),
    )
    return observations, truth


def integrated_signal_log_likelihood(
    f0_particles,
    observations: SignalObservations,
    parameters: SignalModelParameters,
    *,
    channels=("molecular", "mri"),
    mri_template: np.ndarray | None = None,
) -> np.ndarray:
    """Return log p(observed signals | f0 particles) for particle-filter weighting.

    The joint likelihood follows the Gaussian dependence model used by the generator.
    The optional history-only MRI transform is affine. Its Jacobian is constant across
    f0 particles, so the transformed Gaussian likelihood is sufficient for comparing
    particles.
    """
    if not all(hasattr(observations, name) for name in ("molecular_report", "mri_score")):
        raise TypeError("observations must provide molecular_report and mri_score")
    if not all(hasattr(parameters, name) for name in ("q", "a", "b", "d")):
        raise TypeError("parameters must provide q, a, b, and d")
    parameters = SignalModelParameters(
        q=parameters.q, a=parameters.a, b=parameters.b, d=parameters.d
    )
    channels = frozenset(channels)
    unknown = channels - {"molecular", "mri"}
    if unknown:
        raise ValueError(f"unknown signal channels: {sorted(unknown)}")

    particles = np.asarray(f0_particles, dtype=float)
    if particles.ndim < 1 or particles.shape[0] < 1:
        raise ValueError("f0_particles must have patients on axis 0")
    if not np.isfinite(particles).all() or np.any(particles <= 0.0) or np.any(particles >= 1.0):
        raise ValueError("f0_particles must be finite and in (0, 1)")
    n_patients = particles.shape[0]
    reports = np.asarray(observations.molecular_report)
    mri_scores = np.asarray(observations.mri_score, dtype=float)
    if reports.shape != (n_patients,) or mri_scores.shape != (n_patients,):
        raise ValueError("signal arrays must have one value per patient")
    if not np.isin(reports, [False, True, 0, 1]).all():
        raise ValueError("molecular_report must contain only boolean values")
    reports = reports.astype(bool, copy=False)
    if not np.isfinite(mri_scores).all():
        raise ValueError("mri_score must be finite")

    resistance_z = f0_to_resistance_score(particles)
    result = np.zeros_like(resistance_z, dtype=float)
    q, a, b, d = parameters.q, parameters.a, parameters.b, parameters.d
    var_bio = max(0.0, 1.0 - a * a)
    var_mri = max(0.0, 1.0 - b * b)

    mri_z = None
    if "mri" in channels:
        mri_z = (
            _uncalibrate_mri_score(mri_scores, mri_template)
            if mri_template is not None
            else mri_scores
        )[:, None]
        safe_var_mri = max(var_mri, 1e-12)
        result += -0.5 * (
            np.log(2.0 * np.pi * safe_var_mri)
            + ((mri_z - b * resistance_z) ** 2) / safe_var_mri
        )

    if "molecular" in channels:
        mean_bio = a * resistance_z
        if mri_z is not None and var_mri > 1e-12:
            mean_bio = mean_bio + (
                np.sqrt(var_bio) * d / np.sqrt(var_mri)
            ) * (mri_z - b * resistance_z)
            var_conditional = var_bio * max(0.0, 1.0 - d * d)
        else:
            var_conditional = var_bio
        p_state_positive = _normal_positive_probability(mean_bio, var_conditional)
        p_report_positive = q * p_state_positive + (1.0 - q) * (1.0 - p_state_positive)
        p_report = np.where(reports[:, None], p_report_positive, 1.0 - p_report_positive)
        result += np.log(np.clip(p_report, 1e-300, 1.0))
    return result


def _normal_positive_probability(mean, variance):
    """Return P(N(mean, variance) >= 0), including deterministic limits."""
    mean = np.asarray(mean, dtype=float)
    variance = float(variance)
    if variance <= 0.0:
        return np.where(mean >= 0.0, 1.0, 0.0)
    return ndtr(mean / np.sqrt(variance))


__all__ = [
    "F0_REFERENCE_BOUNDS",
    "SignalModelParameters",
    "SignalObservations",
    "SignalTruth",
    "f0_to_resistance_score",
    "generate_integrated_signals",
    "integrated_signal_log_likelihood",
    "standard_normal_marginal_signature",
]
