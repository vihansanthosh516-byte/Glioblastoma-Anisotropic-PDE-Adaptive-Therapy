"""Leave-one-out model averaging for the Track C paired-signal belief update.

The shared signal parameters are unknown. For each virtual patient, this module
weights signal models using only the other patients' paired baseline signals, then
averages that patient's conditional signal likelihood over those weights. It does
not use latent f0 values, treatment outcomes, or patient data.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
from scipy.special import logsumexp

from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
    integrated_signal_log_likelihood,
)


def _load_sibling(name: str, filename: str):
    path = Path(__file__).with_name(filename)
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load {filename}: {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_GRID = _load_sibling("trackc_signal_fisher_grid_for_loo", "212_trackc_signal_fisher_grid.py")
_DEPENDENCE = _load_sibling("trackc_dependence_for_loo", "210_trackc_dependence_identifiability.py")
_MPC = _load_sibling("trackc_robust_belief_mpc_for_loo", "208_trackc_robust_belief_mpc.py")


DEFAULT_MODEL_GRID = tuple(
    SignalModelParameters(*point) for point in _GRID.parameter_grid()
)


def _validate_observations(observations: SignalObservations) -> SignalObservations:
    if not all(hasattr(observations, name) for name in ("molecular_report", "mri_score")):
        raise TypeError("observations must provide molecular_report and mri_score")
    reports = np.asarray(observations.molecular_report)
    scores = np.asarray(observations.mri_score, dtype=float)
    if reports.ndim != 1 or scores.shape != reports.shape or reports.size < 1:
        raise ValueError("observations must have equal, nonempty one-dimensional arrays")
    if not np.isin(reports, (False, True, 0, 1)).all():
        raise ValueError("molecular reports must be binary")
    if not np.isfinite(scores).all():
        raise ValueError("MRI scores must be finite")
    return SignalObservations(reports.astype(bool, copy=True), scores.copy())


def leave_one_out_log_model_weights(
    observations: SignalObservations,
    *,
    model_parameters=DEFAULT_MODEL_GRID,
    f0_bounds=F0_REFERENCE_BOUNDS,
    prior_weights=None,
    n_f0_nodes=256,
) -> np.ndarray:
    """Return normalized ``(patient, model)`` log posterior weights.

    For patient i, the weight of each shared signal model uses the marginal
    likelihood from every patient except i. The latent f0 values are integrated
    out under the fixed declared log-uniform population prior.
    """
    observations = _validate_observations(observations)
    models = tuple(model_parameters)
    if not models or any(not isinstance(row, SignalModelParameters) for row in models):
        raise ValueError("model_parameters must be a nonempty sequence of SignalModelParameters")
    if len(set(models)) != len(models):
        raise ValueError("model_parameters must not contain duplicates")
    if isinstance(n_f0_nodes, bool) or int(n_f0_nodes) != n_f0_nodes or n_f0_nodes < 2:
        raise ValueError("n_f0_nodes must be an integer of at least 2")

    if prior_weights is None:
        prior = np.full(len(models), 1.0 / len(models), dtype=float)
    else:
        prior = np.asarray(prior_weights, dtype=float)
        if prior.shape != (len(models),) or not np.isfinite(prior).all() or (prior < 0).any():
            raise ValueError("prior_weights must be a finite nonnegative value per model")
        if prior.sum() <= 0.0:
            raise ValueError("prior_weights must have positive total mass")
        prior = prior / prior.sum()
    log_prior = np.full(prior.shape, -np.inf, dtype=float)
    positive_prior = prior > 0.0
    log_prior[positive_prior] = np.log(prior[positive_prior])

    # Rows are candidate models; columns are patients. Marginalization integrates
    # out f0, so no hidden simulator values enter peer-based model weighting.
    log_likelihood = np.stack(
        [
            _DEPENDENCE.marginal_signal_log_density(
                observations,
                parameters,
                f0_bounds=f0_bounds,
                n_nodes=int(n_f0_nodes),
            )
            for parameters in models
        ],
        axis=0,
    )
    # Sum explicit peer columns. Subtracting each patient's likelihood from a
    # cohort total is algebraically equivalent, but can let that patient's own
    # value affect its weights through floating-point cancellation.
    peer_log_likelihood = np.stack(
        [
            np.sum(np.delete(log_likelihood, patient, axis=1), axis=1)
            for patient in range(observations.molecular_report.size)
        ],
        axis=0,
    )
    unnormalized = peer_log_likelihood + log_prior[None, :]
    return unnormalized - logsumexp(unnormalized, axis=1, keepdims=True)


def loo_model_averaged_signal_log_likelihood(
    f0_particles,
    observations: SignalObservations,
    model_parameters,
    log_model_weights,
) -> np.ndarray:
    """Mix patient-specific conditional signal likelihoods using peer-only weights."""
    observations = _validate_observations(observations)
    models = tuple(model_parameters)
    weights = np.asarray(log_model_weights, dtype=float)
    particles = np.asarray(f0_particles, dtype=float)
    n_patients = observations.molecular_report.size
    if particles.ndim != 2 or particles.shape[0] != n_patients:
        raise ValueError("f0_particles must have shape (patients, particles)")
    if not np.isfinite(particles).all() or (particles <= 0.0).any() or (particles >= 1.0).any():
        raise ValueError("f0_particles must be finite and in (0, 1)")
    if weights.shape != (n_patients, len(models)) or not np.isfinite(weights).all():
        raise ValueError("log_model_weights must have shape (patients, models) and be finite")
    if not np.allclose(logsumexp(weights, axis=1), 0.0, atol=1e-10, rtol=0.0):
        raise ValueError("each patient's log model weights must be normalized")

    conditional = np.stack(
        [
            integrated_signal_log_likelihood(
                particles,
                observations,
                parameters,
                channels=("molecular", "mri"),
            )
            for parameters in models
        ],
        axis=1,
    )
    return logsumexp(weights[:, :, None] + conditional, axis=1)


class LOOModelAveragedBeliefMPC(_MPC.RobustBeliefMPC):
    """Scan-driven belief MPC with a leave-one-out model-averaged signal update."""

    def __init__(
        self,
        observations: SignalObservations,
        *,
        model_parameters=DEFAULT_MODEL_GRID,
        f0_bounds=F0_REFERENCE_BOUNDS,
        prior_weights=None,
        n_f0_nodes=256,
        **kwargs,
    ):
        self.observations = _validate_observations(observations)
        self.model_parameters = tuple(model_parameters)
        if not self.model_parameters or any(
            not isinstance(row, SignalModelParameters) for row in self.model_parameters
        ):
            raise ValueError("model_parameters must be a nonempty sequence of SignalModelParameters")
        self.f0_bounds = tuple(map(float, f0_bounds))
        self.prior_weights = None if prior_weights is None else np.asarray(prior_weights, dtype=float).copy()
        self.n_f0_nodes = int(n_f0_nodes)
        if "ranges" in kwargs:
            self.ranges = kwargs.pop("ranges")
        else:
            self.ranges = _MPC.INTEGRATED_BELIEF_RANGES
        super().__init__(ranges=self.ranges, **kwargs)

    def reset(self, scen, rng):
        """Fit peer-only model weights, initialize particles, then mix signal likelihoods."""
        if "sigma" not in scen:
            raise ValueError("scen must provide sigma for scan likelihood")
        n_patients = len(np.asarray(scen["sigma"]))
        if self.observations.molecular_report.shape != (n_patients,):
            raise ValueError("baseline observations must match the simulated cohort size")
        self.log_model_weights = leave_one_out_log_model_weights(
            self.observations,
            model_parameters=self.model_parameters,
            f0_bounds=self.f0_bounds,
            prior_weights=self.prior_weights,
            n_f0_nodes=self.n_f0_nodes,
        )
        super().reset(scen, rng)
        mixed_log_likelihood = loo_model_averaged_signal_log_likelihood(
            self.pf.theta["f0"],
            self.observations,
            self.model_parameters,
            self.log_model_weights,
        )
        if mixed_log_likelihood.shape != self.pf.logw.shape or not np.isfinite(mixed_log_likelihood).all():
            raise ArithmeticError("model-averaged signal update is invalid")
        self.pf.logw += mixed_log_likelihood


__all__ = [
    "DEFAULT_MODEL_GRID",
    "LOOModelAveragedBeliefMPC",
    "leave_one_out_log_model_weights",
    "loo_model_averaged_signal_log_likelihood",
]
