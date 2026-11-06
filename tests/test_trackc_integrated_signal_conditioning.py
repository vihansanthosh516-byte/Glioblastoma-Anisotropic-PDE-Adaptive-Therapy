from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "src"))

import trackc_sim as tc
from trackc_integrated_signal_model import (
    F0_REFERENCE_BOUNDS,
    SignalModelParameters,
    SignalObservations,
    f0_to_resistance_score,
    generate_integrated_signals,
    integrated_signal_log_likelihood,
)
from src.trackc_integrated_signal_model import (
    SignalModelParameters as NamespacedSignalParameters,
    SignalObservations as NamespacedSignalObservations,
)

_SPEC = importlib.util.spec_from_file_location(
    "integrated_signal_mpc", REPO / "src" / "208_trackc_robust_belief_mpc.py"
)
_MPC = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MPC
_SPEC.loader.exec_module(_MPC)


def test_f0_reference_score_is_monotone_and_covers_calibration_envelope():
    low, high = F0_REFERENCE_BOUNDS
    values = np.array([low, 0.001, 0.01, 0.1, 0.9])
    scores = f0_to_resistance_score(values)
    assert np.isfinite(scores).all()
    assert np.all(np.diff(scores) > 0)
    assert scores[0] < -5.0
    assert scores[-1] > 5.0
    assert low == pytest.approx(0.001 * np.exp(-0.7))
    assert high == pytest.approx(0.9)


def test_generator_uses_supplied_f0_for_both_signal_channels():
    low, high = F0_REFERENCE_BOUNDS
    f0 = np.geomspace(low, high, 1_000)
    observed, truth = generate_integrated_signals(
        f0.size, q=1.0, a=1.0, b=1.0, d=0.0, seed=20261105, f0=f0
    )
    np.testing.assert_allclose(truth.resistance_z, f0_to_resistance_score(f0))
    np.testing.assert_array_equal(observed.molecular_report, truth.molecular_state)
    np.testing.assert_allclose(observed.mri_score, truth.resistance_z)


def test_zero_signal_association_gives_flat_f0_likelihood_even_with_residual_dependence():
    f0_particles = np.array([[0.001, 0.004, 0.02, 0.2, 0.8]])
    observations = SignalObservations(
        molecular_report=np.array([True]), mri_score=np.array([0.7])
    )
    log_likelihood = integrated_signal_log_likelihood(
        f0_particles,
        observations,
        SignalModelParameters(q=0.9, a=0.0, b=0.0, d=0.7),
    )
    assert np.ptp(log_likelihood[0]) < 1e-12


@pytest.mark.parametrize(("a", "ascending"), [(0.8, True), (-0.8, False)])
def test_molecular_likelihood_order_tracks_declared_link_sign(a, ascending):
    f0_particles = np.array([[0.001, 0.8]])
    observations = SignalObservations(
        molecular_report=np.array([True]), mri_score=np.array([0.0])
    )
    log_likelihood = integrated_signal_log_likelihood(
        f0_particles,
        observations,
        SignalModelParameters(q=0.9, a=a, b=0.0, d=0.0),
        channels=("molecular",),
    )
    assert bool(log_likelihood[0, 1] > log_likelihood[0, 0]) == ascending


@pytest.mark.parametrize(("b", "ascending"), [(0.8, True), (-0.8, False)])
def test_mri_likelihood_order_tracks_declared_link_sign(b, ascending):
    f0_particles = np.array([[0.001, 0.8]])
    observations = SignalObservations(
        molecular_report=np.array([False]), mri_score=np.array([1.0])
    )
    log_likelihood = integrated_signal_log_likelihood(
        f0_particles,
        observations,
        SignalModelParameters(q=0.9, a=0.0, b=b, d=0.0),
        channels=("mri",),
    )
    assert bool(log_likelihood[0, 1] > log_likelihood[0, 0]) == ascending


def test_empirical_mri_transform_produces_finite_particle_likelihood():
    observations = SignalObservations(
        molecular_report=np.array([True]), mri_score=np.array([2.5])
    )
    log_likelihood = integrated_signal_log_likelihood(
        np.array([[0.005, 0.02, 0.2]]),
        observations,
        SignalModelParameters(q=0.75, a=0.5, b=0.5, d=0.0),
        mri_template=np.array([-2.0, -0.5, 1.0, 3.0]),
    )
    assert log_likelihood.shape == (1, 3)
    assert np.isfinite(log_likelihood).all()


def test_correct_joint_likelihood_prefers_generating_f0_on_average():
    n_per_group = 2_000
    f0 = np.r_[np.full(n_per_group, 0.002), np.full(n_per_group, 0.5)]
    parameters = SignalModelParameters(q=0.85, a=0.7, b=0.6, d=0.25)
    observations, _ = generate_integrated_signals(
        f0.size,
        q=parameters.q,
        a=parameters.a,
        b=parameters.b,
        d=parameters.d,
        seed=20261112,
        f0=f0,
    )
    candidates = np.tile(np.array([0.002, 0.5]), (f0.size, 1))
    log_likelihood = integrated_signal_log_likelihood(
        candidates, observations, parameters
    )
    difference = log_likelihood[:, 1] - log_likelihood[:, 0]
    assert difference[:n_per_group].mean() < 0.0
    assert difference[n_per_group:].mean() > 0.0


def _small_policy(observations, parameters, channels):
    return _MPC.SignalConditionedBeliefMPC(
        observations,
        parameters,
        channels=channels,
        scan_every=1,
        budget=2,
        n_particles=128,
        prediction_particles=8,
        forecast_cycles=1,
    )


def test_policy_conditions_on_observations_but_not_hidden_true_f0():
    scenarios = tc.draw_scenarios(3, 20261106)
    parameters = SignalModelParameters(q=0.9, a=0.5, b=0.5, d=0.0)
    observations, _ = generate_integrated_signals(
        3,
        q=0.8,
        a=0.7,
        b=0.7,
        d=0.2,
        seed=20261107,
        f0=scenarios["f0"],
    )
    policy_view = {"sigma": scenarios["sigma"].copy()}

    first = _small_policy(observations, parameters, ("molecular", "mri"))
    second = _small_policy(observations, parameters, ("molecular", "mri"))
    first.reset(policy_view, np.random.default_rng(91))
    second.reset(policy_view, np.random.default_rng(91))

    assert np.any(first.pf.logw != 0.0)
    np.testing.assert_array_equal(first.pf.theta["f0"], second.pf.theta["f0"])
    np.testing.assert_array_equal(first.pf.logw, second.pf.logw)
    y = np.log(np.maximum(scenarios["B0"], 1e-9))
    np.testing.assert_array_equal(first.act(0, y), second.act(0, y))


def test_no_signal_policy_matches_base_mpc_with_same_integrated_prior():
    scenarios = tc.draw_scenarios(2, 20261108)
    parameters = SignalModelParameters(q=0.8, a=0.5, b=0.5, d=0.0)
    no_signal = _small_policy(None, parameters, ())
    base = _MPC.RobustBeliefMPC(
        ranges=_MPC.INTEGRATED_BELIEF_RANGES,
        scan_every=1,
        budget=2,
        n_particles=128,
        prediction_particles=8,
        forecast_cycles=1,
    )
    no_signal.reset(scenarios, np.random.default_rng(31))
    base.reset(scenarios, np.random.default_rng(31))
    np.testing.assert_array_equal(no_signal.pf.theta["f0"], base.pf.theta["f0"])
    np.testing.assert_array_equal(no_signal.pf.logw, base.pf.logw)
    readings = np.log(np.maximum(scenarios["B0"], 1e-9))
    np.testing.assert_array_equal(no_signal.act(0, readings), base.act(0, readings))


def test_policy_accepts_namespaced_signal_types_and_runs_through_simulator():
    scenarios = tc.draw_scenarios(2, 20261109)
    namespaced_observations = NamespacedSignalObservations(
        molecular_report=np.array([True, False]), mri_score=np.array([0.5, -0.5])
    )
    namespaced_parameters = NamespacedSignalParameters(q=0.8, a=0.5, b=0.5, d=0.0)
    policy = _MPC.SignalConditionedBeliefMPC(
        namespaced_observations,
        namespaced_parameters,
        scan_every=1,
        budget=2,
        n_particles=64,
        prediction_particles=4,
        forecast_cycles=1,
    )
    result = tc.simulate(scenarios, policy, seed=20261110, endpoint="baseline", budget=2)
    assert result["ttp"].shape == (2,)
    assert result["courses"].shape == (2,)
    assert np.isfinite(result["ttp"]).all()
    assert np.all(result["courses"] <= 2)


@pytest.mark.parametrize(
    ("observations", "parameters", "channels"),
    [
        (None, SignalModelParameters(0.8, 0.5, 0.5, 0.0), ("molecular",)),
        (
            SignalObservations(np.array([True]), np.array([0.0])),
            SignalModelParameters(0.8, 0.5, 0.5, 0.0),
            ("unknown",),
        ),
    ],
)
def test_policy_rejects_missing_or_unknown_signal_channels(observations, parameters, channels):
    with pytest.raises((TypeError, ValueError)):
        _MPC.SignalConditionedBeliefMPC(observations, parameters, channels=channels)
