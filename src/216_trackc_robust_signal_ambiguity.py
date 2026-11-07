"""Finite-set minimax signal MPC prototype for Track C synthetic worlds.

The controller keeps one particle filter per predeclared signal model. It scores
each short treatment pattern under every model and chooses the pattern with the
best worst-case posterior-predictive utility. The supplied model set is fixed;
this module does not estimate biology or validate the signal-to-resistance link.
"""
from __future__ import annotations

import copy
import importlib.util
import sys
from pathlib import Path

import numpy as np

import trackc_sim as tc
from trackc_integrated_signal_model import (
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


_MPC = _load_sibling("trackc_robust_signal_base_mpc", "208_trackc_robust_belief_mpc.py")


def maximin_action_index(model_action_values) -> int:
    """Choose action with greatest worst-case value; ties go to lowest index."""
    values = np.asarray(model_action_values, dtype=float)
    if values.ndim != 2 or not values.shape[0] or not values.shape[1]:
        raise ValueError("model_action_values must be a nonempty 2D array")
    if not np.isfinite(values).all():
        raise ValueError("model_action_values must contain only finite values")
    return int(np.argmax(np.min(values, axis=0)))


def keyed_stratified_uniforms(seed, patient_indices, *, cycle: int, n_draws: int) -> np.ndarray:
    """Make stable per-patient stratified uniforms, independent of model-set size/order."""
    for value, name in ((seed, "seed"), (cycle, "cycle"), (n_draws, "n_draws")):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{name} must be an integer")
    if seed < 0 or cycle < 0 or n_draws < 1:
        raise ValueError("seed and cycle must be non-negative and n_draws must be positive")
    rows = np.asarray(patient_indices)
    if rows.ndim != 1 or rows.dtype.kind not in "iu" or rows.dtype.kind == "b" or (rows < 0).any():
        raise ValueError("patient_indices must be a one-dimensional array of non-negative integers")
    uniforms = np.empty((rows.size, int(n_draws)), dtype=float)
    strata = np.arange(int(n_draws), dtype=float)
    for index, patient in enumerate(rows):
        rng = np.random.default_rng(
            np.random.SeedSequence([int(seed), int(cycle), int(patient)])
        )
        uniforms[index] = (strata + rng.random(int(n_draws))) / int(n_draws)
    return uniforms


def stratified_particle_indices(weights, uniforms) -> np.ndarray:
    """Select weighted particles with supplied common stratified uniforms."""
    probabilities = np.asarray(weights, dtype=float)
    quantiles = np.asarray(uniforms, dtype=float)
    if (
        probabilities.ndim != 2
        or not probabilities.shape[0]
        or not probabilities.shape[1]
        or not np.isfinite(probabilities).all()
        or (probabilities < 0.0).any()
        or (probabilities.sum(axis=1) <= 0.0).any()
    ):
        raise ValueError("weights must be a nonempty finite non-negative 2D array with positive row sums")
    if (
        quantiles.ndim != 2
        or quantiles.shape[0] != probabilities.shape[0]
        or not quantiles.shape[1]
        or not np.isfinite(quantiles).all()
        or (quantiles < 0.0).any()
        or (quantiles >= 1.0).any()
    ):
        raise ValueError("uniforms must be a finite 2D array with matching rows and values in [0, 1)")
    normalized = probabilities / probabilities.sum(axis=1, keepdims=True)
    cumulative = np.cumsum(normalized, axis=1)
    cumulative[:, -1] = 1.0
    indices = (quantiles[:, :, None] >= cumulative[:, None, :]).sum(axis=2)
    return np.minimum(indices, probabilities.shape[1] - 1).astype(np.intp, copy=False)


def _copy_observations(observations):
    if not all(hasattr(observations, name) for name in ("molecular_report", "mri_score")):
        raise TypeError("observations must provide molecular_report and mri_score")
    reports = np.asarray(observations.molecular_report)
    scores = np.asarray(observations.mri_score, dtype=float)
    if reports.ndim != 1 or scores.ndim != 1 or reports.shape != scores.shape or not reports.size:
        raise ValueError("signal observations must be nonempty one-dimensional arrays of equal length")
    if not np.isin(reports, (False, True, 0, 1)).all():
        raise ValueError("molecular reports must be binary")
    if not np.isfinite(scores).all():
        raise ValueError("MRI signal scores must be finite")
    return SignalObservations(reports.astype(bool, copy=True), scores.copy())


class SignalAmbiguityBeliefMPC(_MPC.RobustBeliefMPC):
    """Receding-horizon MPC robust to a finite, caller-frozen signal-model set.

    Each model has its own particle filter, conditioned on the same patient's
    paired baseline signals and serial burden scans. At each visible scan, each
    treatment pattern receives a posterior-predictive utility under every model.
    The policy applies the pattern that maximizes the minimum model-specific
    utility. No model is selected after seeing evaluation outcomes.

    This is a finite-set minimax comparator, not a learned or clinically
    calibrated model. Results depend on the supplied set and simulator assumptions.
    """

    def __init__(
        self,
        observations: SignalObservations,
        *,
        model_parameters,
        mri_template=None,
        **kwargs,
    ):
        self.observations = _copy_observations(observations)
        models = tuple(model_parameters)
        if not models or any(not isinstance(model, SignalModelParameters) for model in models):
            raise ValueError("model_parameters must be a nonempty sequence of SignalModelParameters")
        if len(set(models)) != len(models):
            raise ValueError("model_parameters must not contain duplicates")
        self.model_parameters = tuple(sorted(models, key=lambda x: (x.q, x.a, x.b, x.d)))
        self.mri_template = None if mri_template is None else np.asarray(mri_template, dtype=float).copy()
        if self.mri_template is not None and not np.isfinite(self.mri_template).all():
            raise ValueError("mri_template must contain only finite values")
        ranges = kwargs.pop("ranges", None)
        if ranges is None:
            # Match simulator f0 support by default. Other world families must pass
            # their frozen support explicitly, as the shared policy benchmark does.
            ranges = dict(_MPC.BELIEF_RANGES)
            ranges["f0"] = tuple(tc.RANGES["f0"])
        super().__init__(ranges=ranges, **kwargs)
        self.model_filters = []
        self.model_nadir_est = None

    def reset(self, scen, rng):
        """Clone one common prior particle cloud, then condition each copy on one model."""
        if "sigma" not in scen:
            raise ValueError("scen must provide sigma for scan likelihood")
        n_patients = len(np.asarray(scen["sigma"]))
        if self.observations.molecular_report.shape != (n_patients,):
            raise ValueError("signal arrays must have one value per simulated patient")
        self._prediction_seed = (
            int(rng.integers(0, 2**32, dtype=np.uint64))
            if len(self.model_parameters) > 1 else None
        )
        super().reset(scen, rng)
        common_prior = copy.deepcopy(self.pf)
        filters = []
        for model_index, parameters in enumerate(self.model_parameters):
            # Keep first filter attached to the policy RNG. This makes the one-model
            # case exactly match SignalConditionedBeliefMPC's random-stream contract.
            particle_filter = self.pf if model_index == 0 else copy.deepcopy(common_prior)
            signal_log_likelihood = integrated_signal_log_likelihood(
                particle_filter.theta["f0"],
                self.observations,
                parameters,
                channels=("molecular", "mri"),
                mri_template=self.mri_template,
            )
            if signal_log_likelihood.shape != particle_filter.logw.shape:
                raise ArithmeticError("signal likelihood has the wrong shape")
            if not np.isfinite(signal_log_likelihood).all():
                raise ArithmeticError("signal likelihood contains non-finite values")
            particle_filter.logw += signal_log_likelihood
            filters.append(particle_filter)
        self.model_filters = filters
        self.pf = self.model_filters[0]
        self.model_nadir_est = np.full((len(self.model_filters), n_patients), np.nan, dtype=float)
        self.nadir_est = self.model_nadir_est[0]

    def _update_visible(self, j, y):
        """Apply same observed burden likelihood to every model-specific filter."""
        y = np.asarray(y, dtype=float)
        if y.shape != (self.n,):
            raise ValueError(f"expected {self.n} readings, got shape {y.shape}")
        seen = np.flatnonzero(~np.isnan(y))
        if seen.size == 0:
            return seen
        if not np.isfinite(y[seen]).all():
            raise ValueError("visible readings must be finite or NaN")
        for particle_filter in self.model_filters:
            burden = np.maximum(particle_filter.S + particle_filter.R, 1e-300)
            particle_filter.logw[seen] += -0.5 * (
                (y[seen, None] - np.log(burden[seen])) / particle_filter.sigma[seen]
            ) ** 2
            weights = particle_filter.w
            effective_sample_size = 1.0 / (weights ** 2).sum(1)
            resample = seen[effective_sample_size[seen] < particle_filter.P / 2]
            if resample.size:
                particle_filter._resample_move(resample, j, weights[resample])
        self.pf = self.model_filters[0]
        return seen

    def _sample_from_filter(self, particle_filter, rows, *, uniforms=None):
        if uniforms is not None:
            weights = particle_filter.w[rows]
            pick = stratified_particle_indices(weights, uniforms)
            row_grid = rows[:, None]
            return {
                "S": particle_filter.S[row_grid, pick].copy(),
                "R": particle_filter.R[row_grid, pick].copy(),
                "theta": {
                    name: value[row_grid, pick].copy()
                    for name, value in particle_filter.theta.items()
                },
            }
        previous = self.pf
        self.pf = particle_filter
        try:
            return super()._sample_posterior(rows)
        finally:
            self.pf = previous

    def _replan(self, rows, j, nadir_prev):
        """Choose one shared plan using worst-case utility across signal models."""
        if rows.size == 0:
            return
        previous_nadir = np.asarray(nadir_prev, dtype=float)
        if previous_nadir.shape == (rows.size,):
            previous_nadir = previous_nadir[None, :]
        if previous_nadir.shape != (len(self.model_filters), rows.size):
            raise ValueError("nadir_prev must have one row per signal model and patient")

        # One-model case must reproduce the established fixed-likelihood controller.
        if len(self.model_filters) == 1:
            self.pf = self.model_filters[0]
            self.nadir_est = self.model_nadir_est[0]
            return super()._replan(rows, j, previous_nadir[0])

        scores = np.empty((len(self.model_filters), self.patterns.shape[0], rows.size), dtype=float)
        self._planning_observation_noise = None
        self._planning_scan_available = None
        if self.stop_endpoint == "observed_two_scan":
            noise = self.rng.standard_normal(
                (rows.size, self.prediction_particles, self.forecast_cycles + 1)
            )
            noise *= self.model_filters[0].sigma[rows, 0][:, None, None]
            self._planning_observation_noise = np.zeros(
                (self.n, self.prediction_particles, self.forecast_cycles + 1), dtype=float
            )
            self._planning_observation_noise[rows] = noise
            self._planning_scan_available = np.ones(
                (self.n, self.prediction_particles, self.forecast_cycles + 1), dtype=bool
            )
            probability = self.assumed_future_scan_miss_probability
            if probability >= 1.0:
                self._planning_scan_available[rows] = False
            elif probability > 0.0:
                self._planning_scan_available[rows] = self.rng.random(
                    (rows.size, self.prediction_particles, self.forecast_cycles + 1)
                ) >= probability
        uniforms = keyed_stratified_uniforms(
            self._prediction_seed,
            rows,
            cycle=j,
            n_draws=self.prediction_particles,
        )
        for model_index, particle_filter in enumerate(self.model_filters):
            self.pf = particle_filter
            state = self._sample_from_filter(particle_filter, rows, uniforms=uniforms)
            carrying_capacity = state["theta"]["B0"]
            burden = state["S"] + state["R"]
            t_base_now = np.where(burden >= tc.PROG_FACTOR * carrying_capacity, 0.0, np.nan)
            if j > 0 and j % self.scan_every == 0:
                threshold = tc.RANO_FACTOR * np.maximum(
                    previous_nadir[model_index, :, None], tc.RANO_FLOOR * carrying_capacity
                )
                has_prior_nadir = np.isfinite(previous_nadir[model_index, :, None])
                t_rano_now = np.where(has_prior_nadir & (burden >= threshold), 0.0, np.nan)
            else:
                t_rano_now = np.full(burden.shape, np.nan)
            nadir_start = self.model_nadir_est[model_index, rows].copy()
            for pattern_index, pattern in enumerate(self.patterns):
                scores[model_index, pattern_index] = self._forecast_pattern(
                    state,
                    rows,
                    j,
                    pattern,
                    nadir_start,
                    t_base_now,
                    t_rano_now,
                )

        self.pf = self.model_filters[0]
        best = np.asarray(
            [maximin_action_index(scores[:, :, patient_index]) for patient_index in range(rows.size)],
            dtype=int,
        )
        for offset in range(self.scan_every):
            cycle = j + offset
            if cycle >= tc.N_CYCLES:
                break
            self.plan[rows, cycle] = self.patterns[best, offset] & (self.used[rows] < self.budget)

    def act(self, j, y):
        """Replan only at visible scans; missing scans cause no update or treatment change."""
        y = np.asarray(y, dtype=float)
        if y.shape != (self.n,):
            raise ValueError(f"expected {self.n} readings, got shape {y.shape}")
        if not 0 <= j < tc.N_CYCLES:
            raise ValueError(f"cycle index must be in [0, {tc.N_CYCLES - 1}]")
        if j % self.scan_every == 0:
            self._update_observed_stop(y)
            visible = self._update_visible(j, y)
            if visible.size:
                nadir_prev = self.model_nadir_est[:, visible].copy()
                for model_index, particle_filter in enumerate(self.model_filters):
                    burden_mean = particle_filter.burden_mean()
                    has_prior_nadir = np.isfinite(nadir_prev[model_index])
                    self.model_nadir_est[model_index, visible] = np.where(
                        has_prior_nadir,
                        np.minimum(nadir_prev[model_index], burden_mean[visible]),
                        burden_mean[visible],
                    )
                active_mask = ~self.observed_stopped[visible]
                active = visible[active_mask]
                if active.size:
                    self._replan(active, j, nadir_prev[:, active_mask])
        action = self.plan[:, j].copy() & (self.used < self.budget)
        if self.stop_endpoint == "observed_two_scan":
            action &= ~self.observed_stopped
        return action

    def advance(self, j, treat):
        """Apply same realized action and propagate every model-specific filter."""
        treat = np.asarray(treat, dtype=bool)
        if treat.shape != (self.n,):
            raise ValueError(f"expected {self.n} treatment actions, got shape {treat.shape}")
        applied = treat & (self.used < self.budget)
        self.used += applied
        for particle_filter in self.model_filters:
            particle_filter.propagate(j, applied)
        self.pf = self.model_filters[0]


__all__ = [
    "SignalAmbiguityBeliefMPC",
    "keyed_stratified_uniforms",
    "maximin_action_index",
    "stratified_particle_indices",
]
