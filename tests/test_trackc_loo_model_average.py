import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.special import logsumexp
import trackc_sim as tc

from trackc_integrated_signal_model import (
    SignalModelParameters,
    SignalObservations,
    integrated_signal_log_likelihood,
)
_SOURCE = Path(__file__).parents[1] / "src" / "214_trackc_loo_model_average.py"
_SPEC = importlib.util.spec_from_file_location("trackc_loo_model_average", _SOURCE)
assert _SPEC is not None and _SPEC.loader is not None
_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = _MODULE
_SPEC.loader.exec_module(_MODULE)
DEFAULT_MODEL_GRID = _MODULE.DEFAULT_MODEL_GRID
LOOModelAveragedBeliefMPC = _MODULE.LOOModelAveragedBeliefMPC
leave_one_out_log_model_weights = _MODULE.leave_one_out_log_model_weights
loo_model_averaged_signal_log_likelihood = _MODULE.loo_model_averaged_signal_log_likelihood
loo_signal_predictive_scores = _MODULE.loo_signal_predictive_scores
marginal_signal_log_densities = _MODULE.marginal_signal_log_densities


def _observations():
    return SignalObservations(
        molecular_report=np.asarray([0, 1, 1, 0, 1], dtype=bool),
        mri_score=np.asarray([-0.7, 0.2, 1.1, -0.1, 0.5], dtype=float),
    )


def _models():
    return (
        SignalModelParameters(q=0.65, a=0.15, b=0.35, d=0.0),
        SignalModelParameters(q=0.8, a=0.45, b=0.6, d=0.35),
        SignalModelParameters(q=0.95, a=0.7, b=0.25, d=-0.25),
    )


def test_default_model_grid_is_unique_and_has_expected_size():
    assert len(DEFAULT_MODEL_GRID) == 81
    assert len(set(DEFAULT_MODEL_GRID)) == len(DEFAULT_MODEL_GRID)


def test_loo_weights_are_normalized_and_exclude_each_patient_exactly():
    observations = _observations()
    models = _models()
    weights = leave_one_out_log_model_weights(
        observations, model_parameters=models, f0_bounds=(0.03, 0.8), n_f0_nodes=24
    )
    np.testing.assert_allclose(logsumexp(weights, axis=1), 0.0, atol=1e-12)

    changed = SignalObservations(
        molecular_report=observations.molecular_report.copy(),
        mri_score=observations.mri_score.copy(),
    )
    changed.molecular_report[2] = ~changed.molecular_report[2]
    changed.mri_score[2] += 4.0
    changed_weights = leave_one_out_log_model_weights(
        changed, model_parameters=models, f0_bounds=(0.03, 0.8), n_f0_nodes=24
    )
    np.testing.assert_array_equal(changed_weights[2], weights[2])
    assert not np.allclose(changed_weights[0], weights[0])


def test_zero_prior_weight_stays_exactly_zero_after_loo_update():
    weights = leave_one_out_log_model_weights(
        _observations(),
        model_parameters=_models(),
        prior_weights=np.asarray([0.0, 0.25, 0.75]),
        f0_bounds=(0.03, 0.8),
        n_f0_nodes=24,
    )
    assert np.isneginf(weights[:, 0]).all()
    np.testing.assert_allclose(logsumexp(weights, axis=1), 0.0, atol=1e-12)


def test_loo_predictive_scores_match_manual_peer_and_uniform_mixtures():
    observations = _observations()
    models = _models()
    result = loo_signal_predictive_scores(
        observations, model_parameters=models, f0_bounds=(0.03, 0.8), n_f0_nodes=24
    )
    log_marginal = marginal_signal_log_densities(
        observations, model_parameters=models, f0_bounds=(0.03, 0.8), n_f0_nodes=24
    )
    weights = leave_one_out_log_model_weights(
        observations, model_parameters=models, f0_bounds=(0.03, 0.8), n_f0_nodes=24
    )
    np.testing.assert_allclose(
        result["loo_log_predictive_density"], logsumexp(weights + log_marginal, axis=1)
    )
    np.testing.assert_allclose(
        result["uniform_grid_log_predictive_density"],
        logsumexp(log_marginal - np.log(len(models)), axis=1),
    )
    np.testing.assert_array_equal(result["log_model_weights"], weights)
    assert np.isfinite(result["loo_log_predictive_density"]).all()


def test_single_model_mixture_matches_direct_conditional_likelihood():
    observations = _observations()
    model = _models()[1]
    particles = np.tile(np.asarray([0.05, 0.2, 0.6, 0.85]), (5, 1))
    log_weights = np.zeros((5, 1))
    mixed = loo_model_averaged_signal_log_likelihood(
        particles, observations, (model,), log_weights
    )
    direct = integrated_signal_log_likelihood(
        particles,
        observations,
        model,
        channels=("molecular", "mri"),
    )
    np.testing.assert_allclose(mixed, direct, rtol=0.0, atol=1e-12)


def test_mpc_reset_adds_finite_peer_weighted_signal_likelihood():
    observations = _observations()
    policy = LOOModelAveragedBeliefMPC(
        observations,
        model_parameters=_models(),
        f0_bounds=(0.03, 0.8),
        n_f0_nodes=24,
        n_particles=32,
        prediction_particles=2,
    )
    policy.reset({"sigma": np.full(5, 0.2)}, np.random.default_rng(17))
    assert policy.log_model_weights.shape == (5, 3)
    assert policy.pf.logw.shape == (5, 32)
    assert np.isfinite(policy.pf.logw).all()
    np.testing.assert_allclose(logsumexp(policy.log_model_weights, axis=1), 0.0, atol=1e-12)


def test_loo_signal_mpc_runs_both_endpoints_on_one_synthetic_path():
    observations = _observations()
    scenarios = tc.draw_scenarios(5, seed=20261106, sigmas=(0.05,))
    policy = LOOModelAveragedBeliefMPC(
        observations,
        model_parameters=_models(),
        f0_bounds=(0.03, 0.8),
        n_f0_nodes=24,
        n_particles=16,
        prediction_particles=2,
        scan_every=3,
        budget=2,
        forecast_cycles=3,
        stop_endpoint="either",
    )

    outcome = tc.simulate_two_endpoints(
        scenarios, policy, seed=[20261106, 3], stop_endpoint="either", budget=2
    )

    assert outcome["stop_endpoint"] == "either"
    assert outcome["ttp_baseline"].shape == (5,)
    assert outcome["ttp_rano"].shape == (5,)
    assert outcome["courses"].shape == (5,)
    assert np.isfinite(outcome["ttp_baseline"]).all()
    assert np.isfinite(outcome["ttp_rano"]).all()
    assert np.all((outcome["courses"] >= 0) & (outcome["courses"] <= 2))


@pytest.mark.parametrize(
    "observations",
    [
        SignalObservations(np.asarray([0, 2]), np.asarray([0.0, 1.0])),
        SignalObservations(np.asarray([0, 1]), np.asarray([0.0, np.nan])),
    ],
)
def test_invalid_observations_are_rejected(observations):
    with pytest.raises(ValueError):
        leave_one_out_log_model_weights(
            observations, model_parameters=_models(), n_f0_nodes=8
        )
