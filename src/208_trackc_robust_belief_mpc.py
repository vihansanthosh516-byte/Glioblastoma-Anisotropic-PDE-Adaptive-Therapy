"""Finite-horizon belief MPC for Track C (simulation only).

This policy uses noisy log-burden readings and its own particle filter. It never reads
the simulated patient's true parameter arrays. At each observed scan it compares every
binary treatment pattern over the next scan interval. It repeats each pattern over a
two-interval prediction window, then applies only the first interval and replans at the
next scan.

The planning utility is the lower of the baseline-referenced and RANO-like restricted
times, minus a fixed course penalty. The forecast is deliberately local: it has no
terminal value after 168 days, uses a broad independent prior for m, omits switching and
reversion, and estimates past RANO nadir from the filter's burden mean. These are
approximations for a bounded comparator, not a clinical controller or a novelty claim.
The outer evaluation must still score both endpoints at 365 and 730 days, actual courses,
and paired lower-tail outcomes.
"""
from __future__ import annotations

from itertools import product

import numpy as np

import trackc_sim as tc
from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
    integrated_signal_log_likelihood,
)


# The anchored development family contains both weak and strong kill values. The default
# ParticleFilter prior (m=6..20) excludes most weak-kill patients, so this policy uses a
# broad support interval. It does not encode the calibrated world's mixture structure.
BELIEF_RANGES = dict(tc.RANGES, m=(0.1, 20.0, "log"))
INTEGRATED_BELIEF_RANGES = dict(
    BELIEF_RANGES, f0=(F0_REFERENCE_BOUNDS[0], F0_REFERENCE_BOUNDS[1], "log")
)


def robust_endpoint_utility(
    ttp_baseline,
    ttp_rano,
    courses_baseline,
    courses_rano=None,
    *,
    tau=730.0,
    course_penalty=5.0,
):
    """Per-draw utility under the worse of two endpoint definitions.

    Each endpoint gets its own course count because progression stops treatment in the
    simulator. The utility is min(RMST_baseline - cost, RMST_RANO - cost). Inputs may be
    scalars or broadcastable arrays. ``tau`` is the restricted-time horizon in days;
    ``course_penalty`` is an assumed number of restricted days per course.
    """
    tau = float(tau)
    course_penalty = float(course_penalty)
    if not np.isfinite(tau) or tau <= 0:
        raise ValueError("tau must be finite and positive")
    if not np.isfinite(course_penalty) or course_penalty < 0:
        raise ValueError("course_penalty must be finite and non-negative")
    if courses_rano is None:
        courses_rano = courses_baseline
    tb, tr, cb, cr = np.broadcast_arrays(
        np.asarray(ttp_baseline, dtype=float),
        np.asarray(ttp_rano, dtype=float),
        np.asarray(courses_baseline, dtype=float),
        np.asarray(courses_rano, dtype=float),
    )
    if not all(np.isfinite(a).all() for a in (tb, tr, cb, cr)):
        raise ValueError("times and course counts must be finite")
    if np.any(tb < 0) or np.any(tr < 0) or np.any(cb < 0) or np.any(cr < 0):
        raise ValueError("times and course counts must be non-negative")
    u_baseline = np.minimum(tb, tau) - course_penalty * cb
    u_rano = np.minimum(tr, tau) - course_penalty * cr
    return np.minimum(u_baseline, u_rano)


class RobustBeliefMPC:
    """Bounded scan-driven MPC using the existing Track C policy interface.

    Parameters:
        scan_every: visible scan interval in 28-day cycles (3 gives 84-day scans).
        budget: maximum courses. Must match ``simulate(..., budget=...)`` when overridden.
        course_penalty: assumed restricted days assigned to one course.
        n_particles: filter particles per patient.
        prediction_particles: posterior draws used by the planner per patient.
        forecast_cycles: local forecast length. Must be a multiple of ``scan_every``;
            the selected scan-interval pattern repeats throughout this window.
        ranges: optional parameter prior passed to ``ParticleFilter``.

    Benchmarking: ``evaluate_robust_belief_mpc`` applies the policy's scan interval and
    course budget to ``simulate``. A harness can replay the same pre-generated ``y``
    scan stream to this policy and every comparator; the policy reads only ``sigma``
    from ``scen`` and the supplied readings. When a scheduled scan has no usable
    reading, this policy holds treatment for that patient until its next usable scan.
    """

    uses_filter = True

    def __init__(
        self,
        scan_every=3,
        budget=tc.BUDGET,
        course_penalty=5.0,
        n_particles=500,
        prediction_particles=16,
        forecast_cycles=6,
        ranges=None,
    ):
        self.scan_every = _as_int(scan_every, "scan_every")
        self.budget = _as_int(budget, "budget")
        self.course_penalty = float(course_penalty)
        self.n_particles = _as_int(n_particles, "n_particles")
        self.prediction_particles = _as_int(prediction_particles, "prediction_particles")
        self.forecast_cycles = _as_int(forecast_cycles, "forecast_cycles")
        if not 1 <= self.scan_every <= 4:
            raise ValueError("scan_every must be between 1 and 4 cycles")
        if not 1 <= self.budget <= tc.BUDGET:
            raise ValueError(f"budget must be between 1 and {tc.BUDGET}")
        if self.course_penalty < 0 or not np.isfinite(self.course_penalty):
            raise ValueError("course_penalty must be finite and non-negative")
        if self.n_particles < 2:
            raise ValueError("n_particles must be at least 2")
        if not 1 <= self.prediction_particles <= self.n_particles:
            raise ValueError("prediction_particles must be in [1, n_particles]")
        if not self.scan_every <= self.forecast_cycles <= 12:
            raise ValueError("forecast_cycles must be between scan_every and 12")
        if self.forecast_cycles % self.scan_every:
            raise ValueError("forecast_cycles must be a multiple of scan_every")
        self.ranges = BELIEF_RANGES if ranges is None else ranges
        self.patterns = np.asarray(
            sorted(product((False, True), repeat=self.scan_every), key=lambda p: (sum(p), p)),
            dtype=bool,
        )

    def reset(self, scen, rng):
        """Start a filter from scan precision only; ignore true parameter arrays in scen."""
        if "sigma" not in scen:
            raise ValueError("scen must provide sigma for the scan likelihood")
        self.rng = rng
        self.pf = tc.ParticleFilter(
            scen["sigma"], self.n_particles, rng, ranges=self.ranges
        )
        self.n = self.pf.N
        self.used = np.zeros(self.n, dtype=int)
        self.plan = np.zeros((self.n, tc.N_CYCLES), dtype=bool)
        # The RANO nadir starts at the first observed burden, not the unconditioned prior.
        # It remains NaN until a patient has a usable observation.
        self.nadir_est = np.full(self.n, np.nan, dtype=float)

    def _update_visible(self, j, y):
        """Update only rows with an actual reading; missing rows get no likelihood or resample."""
        y = np.asarray(y, dtype=float)
        seen = np.flatnonzero(~np.isnan(y))
        if seen.size == 0:
            return seen
        B = np.maximum(self.pf.S + self.pf.R, 1e-300)
        self.pf.logw[seen] += -0.5 * (
            (y[seen, None] - np.log(B[seen])) / self.pf.sigma[seen]
        ) ** 2
        w = self.pf.w
        ess = 1.0 / (w ** 2).sum(1)
        resample = seen[ess[seen] < self.pf.P / 2]
        if resample.size:
            self.pf._resample_move(resample, j, w[resample])
        return seen

    def _sample_posterior(self, rows):
        """Draw a small weighted posterior sample for bounded planning."""
        weights = self.pf.w[rows]
        pick = np.empty((rows.size, self.prediction_particles), dtype=int)
        for k in range(rows.size):
            pick[k] = self.rng.choice(
                self.pf.P,
                size=self.prediction_particles,
                replace=True,
                p=weights[k],
            )
        row_grid = rows[:, None]
        state = {
            "S": self.pf.S[row_grid, pick].copy(),
            "R": self.pf.R[row_grid, pick].copy(),
            "theta": {
                name: value[row_grid, pick].copy()
                for name, value in self.pf.theta.items()
            },
        }
        return state

    def _forecast_pattern(self, state, rows, j, pattern, nadir_start, t_base_now, t_rano_now):
        """Predict both endpoints under one repeated scan-interval pattern."""
        theta = state["theta"]
        S_base, R_base = state["S"].copy(), state["R"].copy()
        S_rano, R_rano = state["S"].copy(), state["R"].copy()
        B0 = theta["B0"]
        ttp_base = t_base_now.copy()
        ttp_rano = t_rano_now.copy()
        courses_base = np.zeros_like(ttp_base, dtype=float)
        courses_rano = np.zeros_like(ttp_rano, dtype=float)
        nadir = np.minimum(nadir_start[:, None], B0)

        for h in range(self.forecast_cycles):
            cycle = j + h
            on = bool(pattern[h % self.scan_every])
            can_dose_base = (self.used[rows, None] + courses_base < self.budget) & np.isnan(ttp_base)
            can_dose_rano = (self.used[rows, None] + courses_rano < self.budget) & np.isnan(ttp_rano)
            treat_base = can_dose_base & on
            treat_rano = can_dose_rano & on
            courses_base += treat_base
            courses_rano += treat_rano
            S_base, R_base = tc.run_cycle(
                S_base,
                R_base,
                theta,
                treat_base,
                h * tc.CYCLE_DAYS,
                prog_level=tc.PROG_FACTOR * B0,
                ttp=ttp_base,
            )
            S_rano, R_rano = tc.run_cycle(
                S_rano,
                R_rano,
                theta,
                treat_rano,
                h * tc.CYCLE_DAYS,
            )
            if (cycle + 1) % self.scan_every == 0:
                B_rano = S_rano + R_rano
                threshold = tc.RANO_FACTOR * np.maximum(nadir, tc.RANO_FLOOR * B0)
                hit = np.isnan(ttp_rano) & (B_rano >= threshold)
                ttp_rano[hit] = (h + 1) * tc.CYCLE_DAYS
                nadir = np.minimum(nadir, B_rano)

        forecast_days = self.forecast_cycles * tc.CYCLE_DAYS
        ttp_base = np.where(np.isnan(ttp_base), forecast_days, ttp_base)
        ttp_rano = np.where(np.isnan(ttp_rano), forecast_days, ttp_rano)
        utility = robust_endpoint_utility(
            ttp_base,
            ttp_rano,
            courses_base,
            courses_rano,
            tau=forecast_days,
            course_penalty=self.course_penalty,
        )
        return utility.mean(axis=1)

    def _replan(self, rows, j, nadir_prev):
        if rows.size == 0:
            return
        state = self._sample_posterior(rows)
        B0 = state["theta"]["B0"]
        burden = state["S"] + state["R"]
        t_base_now = np.where(burden >= tc.PROG_FACTOR * B0, 0.0, np.nan)
        if j > 0 and j % self.scan_every == 0:
            threshold = tc.RANO_FACTOR * np.maximum(nadir_prev[:, None], tc.RANO_FLOOR * B0)
            has_prior_nadir = np.isfinite(nadir_prev[:, None])
            t_rano_now = np.where(has_prior_nadir & (burden >= threshold), 0.0, np.nan)
        else:
            t_rano_now = np.full(burden.shape, np.nan)

        # Current scan updates the observable nadir estimate. Use that updated value for
        # future predictions; nadir_prev is only for the current scan's RANO call.
        nadir_start = self.nadir_est[rows].copy()
        scores = np.empty((self.patterns.shape[0], rows.size), dtype=float)
        for k, pattern in enumerate(self.patterns):
            scores[k] = self._forecast_pattern(
                state, rows, j, pattern, nadir_start, t_base_now, t_rano_now
            )
        best = np.argmax(scores, axis=0)
        for offset in range(self.scan_every):
            cycle = j + offset
            if cycle >= tc.N_CYCLES:
                break
            self.plan[rows, cycle] = self.patterns[best, offset] & (self.used[rows] < self.budget)

    def act(self, j, y):
        """Replan only at observed scans; repeat the stored block on unobserved cycles."""
        y = np.asarray(y, dtype=float)
        if y.shape != (self.n,):
            raise ValueError(f"expected {self.n} readings, got shape {y.shape}")
        if not 0 <= j < tc.N_CYCLES:
            raise ValueError(f"cycle index must be in [0, {tc.N_CYCLES - 1}]")
        if j % self.scan_every == 0:
            visible = self._update_visible(j, y)
            if visible.size:
                nadir_prev = self.nadir_est[visible].copy()
                burden_mean = self.pf.burden_mean()
                has_prior_nadir = np.isfinite(nadir_prev)
                self.nadir_est[visible] = np.where(
                    has_prior_nadir,
                    np.minimum(nadir_prev, burden_mean[visible]),
                    burden_mean[visible],
                )
                self._replan(visible, j, nadir_prev)
        return self.plan[:, j].copy() & (self.used < self.budget)

    def advance(self, j, treat):
        """Record applied courses and propagate the filter under known past actions."""
        treat = np.asarray(treat, dtype=bool)
        if treat.shape != (self.n,):
            raise ValueError(f"expected {self.n} treatment actions, got shape {treat.shape}")
        applied = treat & (self.used < self.budget)
        self.used += applied
        self.pf.propagate(j, applied)


class SignalConditionedBeliefMPC(RobustBeliefMPC):
    """Robust-belief MPC initialized with explicit molecular and/or MRI signal evidence.

    The policy receives only observed signals and a fixed assumed likelihood. It never
    receives the simulated patient's true f0. Use the same broad integrated prior for
    signal and no-signal comparisons.
    """

    def __init__(
        self,
        observations: SignalObservations | None,
        parameters: SignalModelParameters,
        *,
        channels=("molecular", "mri"),
        mri_template=None,
        **kwargs,
    ):
        self.signal_channels = frozenset(channels)
        if self.signal_channels - {"molecular", "mri"}:
            raise ValueError("channels may contain only 'molecular' and 'mri'")
        if self.signal_channels and not all(
            hasattr(observations, name) for name in ("molecular_report", "mri_score")
        ):
            raise TypeError("observations must provide molecular_report and mri_score")
        if not all(hasattr(parameters, name) for name in ("q", "a", "b", "d")):
            raise TypeError("parameters must provide q, a, b, and d")
        self.signal_parameters = SignalModelParameters(
            q=parameters.q, a=parameters.a, b=parameters.b, d=parameters.d
        )
        if observations is None:
            self.signal_observations = None
        else:
            self.signal_observations = SignalObservations(
                molecular_report=np.asarray(observations.molecular_report).copy(),
                mri_score=np.asarray(observations.mri_score, dtype=float).copy(),
            )
        self.mri_template = None if mri_template is None else np.asarray(mri_template, dtype=float).copy()
        ranges = kwargs.pop("ranges", None)
        super().__init__(
            ranges=INTEGRATED_BELIEF_RANGES if ranges is None else ranges,
            **kwargs,
        )

    def reset(self, scen, rng):
        """Initialize prior particles, then condition weights on supplied baseline signals."""
        n_patients = len(scen["sigma"])
        if self.signal_channels:
            if (
                self.signal_observations.molecular_report.shape != (n_patients,)
                or self.signal_observations.mri_score.shape != (n_patients,)
            ):
                raise ValueError("signal arrays must have one value per simulated patient")
        super().reset(scen, rng)
        if self.signal_channels:
            self.pf.logw += integrated_signal_log_likelihood(
                self.pf.theta["f0"],
                self.signal_observations,
                self.signal_parameters,
                channels=self.signal_channels,
                mri_template=self.mri_template,
            )


def _as_int(value, name):
    """Accept integral numeric settings without silently truncating fractions or booleans."""
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{name} must be an integer")
    try:
        numeric = float(value)
        result = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not np.isfinite(numeric) or numeric != result:
        raise ValueError(f"{name} must be an integer")
    return result


def evaluate_robust_belief_mpc(evaluation_scenarios, policy, seed, *, endpoint="rano"):
    """Evaluate a frozen policy with its configured scan interval and course budget."""
    if not isinstance(policy, RobustBeliefMPC):
        raise TypeError("policy must be a RobustBeliefMPC")
    eval_scenarios = dict(evaluation_scenarios)
    eval_scenarios["scan_every"] = policy.scan_every
    return tc.simulate(
        eval_scenarios, policy, seed, endpoint=endpoint, budget=policy.budget
    )


__all__ = [
    "BELIEF_RANGES",
    "INTEGRATED_BELIEF_RANGES",
    "RobustBeliefMPC",
    "SignalConditionedBeliefMPC",
    "evaluate_robust_belief_mpc",
    "robust_endpoint_utility",
]
