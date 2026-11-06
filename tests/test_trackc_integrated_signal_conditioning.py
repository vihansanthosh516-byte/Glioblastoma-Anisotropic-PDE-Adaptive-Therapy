from __future__ import annotations

import importlib.util
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import quad

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
    standard_normal_marginal_signature,
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


def test_uninformative_molecular_report_has_constant_likelihood_at_q_half():
    f0_particles = np.array([[0.001, 0.004, 0.02, 0.2, 0.8]])
    parameters = SignalModelParameters(q=0.5, a=0.9, b=-0.6, d=0.75)
    for report in (False, True):
        likelihood = integrated_signal_log_likelihood(
            f0_particles,
            SignalObservations(
                molecular_report=np.array([report]), mri_score=np.array([0.0])
            ),
            parameters,
            channels=("molecular",),
        )
        np.testing.assert_allclose(likelihood, np.log(0.5), rtol=0.0, atol=1e-12)


def test_zero_residual_dependence_factorizes_joint_likelihood():
    f0_particles = np.array([[0.001, 0.004, 0.02, 0.2, 0.8]])
    observations = SignalObservations(
        molecular_report=np.array([True]), mri_score=np.array([-0.35])
    )
    parameters = SignalModelParameters(q=0.82, a=0.45, b=-0.3, d=0.0)
    joint = integrated_signal_log_likelihood(
        f0_particles, observations, parameters
    )
    molecular = integrated_signal_log_likelihood(
        f0_particles, observations, parameters, channels=("molecular",)
    )
    mri = integrated_signal_log_likelihood(
        f0_particles, observations, parameters, channels=("mri",)
    )
    np.testing.assert_allclose(joint, molecular + mri, rtol=0.0, atol=1e-12)


def test_residual_dependence_changes_joint_but_not_single_channel_likelihoods():
    f0_particles = np.array([[0.001, 0.004, 0.02, 0.2, 0.8]])
    observations = SignalObservations(
        molecular_report=np.array([True]), mri_score=np.array([0.75])
    )
    independent = SignalModelParameters(q=0.82, a=0.45, b=-0.3, d=0.0)
    dependent = SignalModelParameters(q=0.82, a=0.45, b=-0.3, d=0.7)

    joint_independent = integrated_signal_log_likelihood(
        f0_particles, observations, independent
    )
    joint_dependent = integrated_signal_log_likelihood(
        f0_particles, observations, dependent
    )
    molecular_independent = integrated_signal_log_likelihood(
        f0_particles, observations, independent, channels=("molecular",)
    )
    molecular_dependent = integrated_signal_log_likelihood(
        f0_particles, observations, dependent, channels=("molecular",)
    )
    mri_independent = integrated_signal_log_likelihood(
        f0_particles, observations, independent, channels=("mri",)
    )
    mri_dependent = integrated_signal_log_likelihood(
        f0_particles, observations, dependent, channels=("mri",)
    )

    assert not np.allclose(joint_independent, joint_dependent)
    np.testing.assert_allclose(molecular_independent, molecular_dependent)
    np.testing.assert_allclose(mri_independent, mri_dependent)


def test_standard_normal_marginal_signals_do_not_identify_a_b_and_d_separately():
    first = SignalModelParameters(q=0.8, a=0.5, b=0.4, d=0.6)
    q, rho = standard_normal_marginal_signature(first)

    a_alt, b_alt = 0.8, 0.5
    d_alt = (rho - a_alt * b_alt) / (
        math.sqrt(1.0 - a_alt * a_alt) * math.sqrt(1.0 - b_alt * b_alt)
    )
    second = SignalModelParameters(q=q, a=a_alt, b=b_alt, d=d_alt)

    assert first.d != pytest.approx(second.d)
    assert (first.a, first.b) != (second.a, second.b)
    assert standard_normal_marginal_signature(second) == pytest.approx((q, rho))


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


@pytest.mark.parametrize("report", [False, True])
def test_joint_likelihood_matches_direct_latent_integral(report):
    parameters = SignalModelParameters(q=0.83, a=0.45, b=-0.35, d=0.6)
    f0 = 0.031
    mri_score = 0.4
    observations = SignalObservations(
        molecular_report=np.array([report]), mri_score=np.array([mri_score])
    )
    got = integrated_signal_log_likelihood(
        np.array([[f0]]), observations, parameters
    )

    z_resistance = float(f0_to_resistance_score(np.array([f0]))[0])
    sd_bio = math.sqrt(1.0 - parameters.a**2)
    sd_mri_given_bio = math.sqrt(
        (1.0 - parameters.b**2) * (1.0 - parameters.d**2)
    )

    def integrand(z_bio):
        density_bio = math.exp(
            -0.5 * ((z_bio - parameters.a * z_resistance) / sd_bio) ** 2
        ) / (math.sqrt(2.0 * math.pi) * sd_bio)
        mri_mean = (
            parameters.b * z_resistance
            + math.sqrt(1.0 - parameters.b**2)
            * parameters.d
            * (z_bio - parameters.a * z_resistance)
            / sd_bio
        )
        density_mri = math.exp(
            -0.5 * ((mri_score - mri_mean) / sd_mri_given_bio) ** 2
        ) / (math.sqrt(2.0 * math.pi) * sd_mri_given_bio)
        state_positive = z_bio >= 0.0
        report_probability = (
            parameters.q
            if report == state_positive
            else 1.0 - parameters.q
        )
        return density_bio * density_mri * report_probability

    direct = quad(integrand, -np.inf, 0.0, epsabs=1e-12)[0]
    direct += quad(integrand, 0.0, np.inf, epsabs=1e-12)[0]
    assert float(np.exp(got[0, 0])) == pytest.approx(direct, rel=1e-9, abs=1e-12)


@pytest.mark.parametrize(
    "parameters",
    [
        SignalModelParameters(q=0.83, a=0.45, b=-0.35, d=0.6),
        SignalModelParameters(q=0.67, a=-0.4, b=0.55, d=-0.7),
        SignalModelParameters(q=0.91, a=0.0, b=0.0, d=0.8),
    ],
)
def test_joint_signal_likelihood_normalizes_over_reports_and_mri(parameters):
    f0 = 0.037
    total_probability = 0.0
    for report in (False, True):
        def density(mri_score):
            observations = SignalObservations(
                molecular_report=np.array([report]),
                mri_score=np.array([mri_score]),
            )
            log_likelihood = integrated_signal_log_likelihood(
                np.array([[f0]]), observations, parameters
            )
            return float(np.exp(log_likelihood[0, 0]))

        total_probability += quad(
            density, -np.inf, np.inf, epsabs=1e-10, epsrel=1e-10
        )[0]
    assert total_probability == pytest.approx(1.0, abs=1e-9)


def test_template_calibration_preserves_joint_likelihood_ratios():
    template = np.array([-3.0, -1.0, 0.5, 2.0, 4.0])
    location = float(np.mean(template))
    scale = float(np.std(template, ddof=1))
    latent_score = -0.65
    f0_particles = np.array([[0.003, 0.02, 0.15, 0.7]])
    parameters = SignalModelParameters(q=0.82, a=0.4, b=0.5, d=-0.35)
    reports = np.array([True])

    raw_likelihood = integrated_signal_log_likelihood(
        f0_particles,
        SignalObservations(reports, np.array([latent_score])),
        parameters,
    )
    calibrated_likelihood = integrated_signal_log_likelihood(
        f0_particles,
        SignalObservations(reports, np.array([location + scale * latent_score])),
        parameters,
        mri_template=template,
    )

    np.testing.assert_allclose(
        calibrated_likelihood, raw_likelihood, rtol=0.0, atol=1e-13
    )


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


def test_serial_burden_scans_add_one_likelihood_each_after_baseline_signals():
    observations = SignalObservations(np.array([True]), np.array([0.25]))
    parameters = SignalModelParameters(q=0.85, a=0.5, b=0.4, d=0.6)
    policy = _small_policy(observations, parameters, ("molecular", "mri"))
    policy.reset({"sigma": np.array([1e6])}, np.random.default_rng(22))

    baseline_logw = policy.pf.logw.copy()
    burden = np.maximum(policy.pf.S + policy.pf.R, 1e-300)
    first_reading = np.array([np.log(np.median(burden[0]))])
    first_increment = -0.5 * (
        (first_reading[:, None] - np.log(burden)) / policy.pf.sigma
    ) ** 2
    policy._update_visible(0, first_reading)
    np.testing.assert_array_equal(policy.pf.logw, baseline_logw + first_increment)

    before_second = policy.pf.logw.copy()
    second_reading = first_reading + 0.1
    second_increment = -0.5 * (
        (second_reading[:, None] - np.log(burden)) / policy.pf.sigma
    ) ** 2
    policy._update_visible(1, second_reading)
    np.testing.assert_array_equal(policy.pf.logw, before_second + second_increment)


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
