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
from scipy.special import ndtr


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


def _bounded_correlation(value: float, name: str) -> float:
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be a finite number in [-1, 1]")
    try:
        result = float(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a finite number in [-1, 1]") from exc
    if not np.isfinite(result) or not -1.0 <= result <= 1.0:
        raise ValueError(f"{name} must be a finite number in [-1, 1]")
    return result


def _to_empirical_marginal(values: np.ndarray, template: np.ndarray | None) -> np.ndarray:
    """Map latent Gaussian ranks to a history-only empirical signal distribution."""
    if template is None:
        return values.copy()
    marginal = np.asarray(template, dtype=float).reshape(-1)
    if marginal.size == 0 or not np.isfinite(marginal).all():
        raise ValueError("mri_template must be non-empty and finite")
    probabilities = ndtr(values)
    return np.quantile(marginal, probabilities, method="linear")


def generate_integrated_signals(
    n: int,
    *,
    q: float,
    a: float,
    b: float,
    d: float,
    seed: int,
    mri_template: np.ndarray | None = None,
) -> tuple[SignalObservations, SignalTruth]:
    """Draw paired synthetic signals and separate policy inputs from latent truth.

    Independent standard-normal variables define:

        Z_bio = a * Z_resistance + sqrt(1-a^2) * E_bio
        Z_mri = b * Z_resistance + sqrt(1-b^2) * (d * E_bio + sqrt(1-d^2) * E_mri)

    ``a`` sets the latent Gaussian association that induces the link between the
    Core-like/Periphery-like molecular state and resistance after thresholding.
    ``b`` links the latent MRI score to resistance. ``d`` sets residual dependence.
    A binary molecular report gives the correct sign of ``Z_bio`` with probability
    ``q``. The report and MRI score are the only fields returned in observations.

    The parameters are scenario assumptions. ``mri_template`` must be built from
    training-only or otherwise permitted history data; its empirical marginal is
    applied by a Gaussian-rank transform.
    """
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n < 1:
        raise ValueError("n must be a positive integer")
    if isinstance(seed, (bool, np.bool_)) or not isinstance(seed, Integral) or seed < 0:
        raise ValueError("seed must be a non-negative integer")

    q = _bounded_probability(q, "q")
    a = _bounded_correlation(a, "a")
    b = _bounded_correlation(b, "b")
    d = _bounded_correlation(d, "d")

    rng = np.random.default_rng(int(seed))
    resistance_z = rng.standard_normal(int(n))
    e_bio = rng.standard_normal(int(n))
    e_mri = rng.standard_normal(int(n))

    molecular_state_z = a * resistance_z + np.sqrt(max(0.0, 1.0 - a * a)) * e_bio
    mri_latent_z = b * resistance_z + np.sqrt(max(0.0, 1.0 - b * b)) * (
        d * e_bio + np.sqrt(max(0.0, 1.0 - d * d)) * e_mri
    )
    molecular_state = molecular_state_z >= 0.0
    molecular_report = molecular_state ^ (rng.random(int(n)) >= q)
    mri_score = _to_empirical_marginal(mri_latent_z, mri_template)

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


__all__ = ["SignalObservations", "SignalTruth", "generate_integrated_signals"]
