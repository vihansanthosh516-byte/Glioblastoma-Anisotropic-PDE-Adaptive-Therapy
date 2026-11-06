import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest
import trackc_sim as tc

from trackc_integrated_signal_model import SignalModelParameters, SignalObservations


REPO = Path(__file__).parents[1]


def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, REPO / "src" / filename)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


_MPC = _load("trackc_mpc_for_ambiguity_test", "208_trackc_robust_belief_mpc.py")
_AMBIGUITY = _load("trackc_robust_signal_ambiguity", "216_trackc_robust_signal_ambiguity.py")
SignalAmbiguityBeliefMPC = _AMBIGUITY.SignalAmbiguityBeliefMPC
maximin_action_index = _AMBIGUITY.maximin_action_index


def _models():
    return (
        SignalModelParameters(0.65, 0.15, 0.35, 0.0),
        SignalModelParameters(0.80, 0.45, 0.60, 0.35),
        SignalModelParameters(0.95, 0.70, 0.25, -0.25),
    )


def _observations(n=2):
    return SignalObservations(
        molecular_report=np.asarray([True, False][:n]),
        mri_score=np.asarray([0.4, -0.3][:n], dtype=float),
    )


def _policy(models):
    return SignalAmbiguityBeliefMPC(
        _observations(),
        model_parameters=models,
        scan_every=1,
        budget=2,
        n_particles=64,
        prediction_particles=4,
        forecast_cycles=1,
    )


def test_maximin_action_uses_worst_model_and_first_index_for_ties():
    values = np.asarray([[4.0, 7.0, 2.0], [5.0, 3.0, 2.0]])
    assert maximin_action_index(values) == 0
    tied = np.asarray([[2.0, 2.0], [1.0, 1.0]])
    assert maximin_action_index(tied) == 0


@pytest.mark.parametrize(
    "values",
    [np.empty((0, 2)), np.empty((2, 0)), np.asarray([[1.0, np.nan]]), np.asarray([[1.0, np.inf]])],
)
def test_maximin_action_rejects_empty_or_nonfinite_values(values):
    with pytest.raises(ValueError):
        maximin_action_index(values)


def test_models_are_canonicalized_and_duplicate_models_are_rejected():
    models = _models()
    first = _policy((models[2], models[0], models[1]))
    second = _policy(models)
    assert first.model_parameters == second.model_parameters
    with pytest.raises(ValueError, match="duplicates"):
        _policy((models[0], models[0]))
    with pytest.raises(ValueError, match="nonempty"):
        _policy(())


def test_default_resistance_fraction_support_matches_simulator():
    policy = _policy(_models())
    assert policy.ranges["f0"] == tc.RANGES["f0"]


def test_each_model_gets_own_signal_conditioned_filter_and_scan_updates():
    policy = _policy(_models())
    policy.reset({"sigma": np.full(2, 0.2)}, np.random.default_rng(31))
    assert len(policy.model_filters) == 3
    assert policy.pf is policy.model_filters[0]
    before = [pf.logw.copy() for pf in policy.model_filters]
    assert any(np.any(pf.logw != 0.0) for pf in policy.model_filters)

    policy._update_visible(0, np.asarray([np.nan, np.nan]))
    for pf, original in zip(policy.model_filters, before):
        np.testing.assert_array_equal(pf.logw, original)

    policy._update_visible(0, np.asarray([0.1, -0.2]))
    for pf, original in zip(policy.model_filters, before):
        assert np.any(pf.logw != original)


def test_one_model_matches_fixed_signal_mpc_on_same_random_stream():
    model = _models()[1]
    observations = _observations()
    scenarios = tc.draw_scenarios(2, seed=20261106)
    robust = SignalAmbiguityBeliefMPC(
        observations,
        model_parameters=(model,),
        scan_every=1,
        budget=2,
        n_particles=64,
        prediction_particles=4,
        forecast_cycles=1,
        ranges=_MPC.INTEGRATED_BELIEF_RANGES,
    )
    fixed = _MPC.SignalConditionedBeliefMPC(
        observations,
        model,
        scan_every=1,
        budget=2,
        n_particles=64,
        prediction_particles=4,
        forecast_cycles=1,
    )
    robust.reset(scenarios, np.random.default_rng(91))
    fixed.reset(scenarios, np.random.default_rng(91))
    np.testing.assert_array_equal(robust.pf.theta["f0"], fixed.pf.theta["f0"])
    np.testing.assert_array_equal(robust.pf.logw, fixed.pf.logw)

    for cycle in range(3):
        reading = np.log(np.maximum(scenarios["B0"], 1e-9)) + 0.03 * cycle
        robust_action = robust.act(cycle, reading)
        fixed_action = fixed.act(cycle, reading)
        np.testing.assert_array_equal(robust_action, fixed_action)
        robust.advance(cycle, robust_action)
        fixed.advance(cycle, fixed_action)


def test_uninformative_signal_model_matches_no_signal_mpc():
    model = SignalModelParameters(q=0.5, a=0.0, b=0.0, d=0.0)
    observations = SignalObservations(
        molecular_report=np.asarray([True, False]),
        mri_score=np.asarray([1.2, -0.8]),
    )
    scenarios = tc.draw_scenarios(2, seed=20261109)
    ambiguity = SignalAmbiguityBeliefMPC(
        observations,
        model_parameters=(model,),
        scan_every=1,
        budget=2,
        n_particles=64,
        prediction_particles=4,
        forecast_cycles=1,
        ranges=_MPC.INTEGRATED_BELIEF_RANGES,
    )
    no_signal = _MPC.RobustBeliefMPC(
        scan_every=1,
        budget=2,
        n_particles=64,
        prediction_particles=4,
        forecast_cycles=1,
        ranges=_MPC.INTEGRATED_BELIEF_RANGES,
    )

    ambiguity_outcome = tc.simulate_two_endpoints(
        scenarios,
        ambiguity,
        seed=[20261109, 3],
        stop_endpoint="either",
        budget=2,
    )
    no_signal_outcome = tc.simulate_two_endpoints(
        scenarios,
        no_signal,
        seed=[20261109, 3],
        stop_endpoint="either",
        budget=2,
    )

    for name in ("ttp_baseline", "ttp_rano", "courses"):
        np.testing.assert_array_equal(ambiguity_outcome[name], no_signal_outcome[name])


def test_ambiguity_policy_runs_shared_endpoint_simulation_with_finite_budget():
    scenarios = tc.draw_scenarios(2, seed=20261107)
    policy = _policy(_models())
    outcome = tc.simulate_two_endpoints(
        scenarios,
        policy,
        seed=[20261107, 2],
        stop_endpoint="either",
        budget=2,
    )
    assert outcome["ttp_baseline"].shape == (2,)
    assert outcome["ttp_rano"].shape == (2,)
    assert np.isfinite(outcome["ttp_baseline"]).all()
    assert np.isfinite(outcome["ttp_rano"]).all()
    assert np.all((outcome["courses"] >= 0) & (outcome["courses"] <= 2))


def test_model_input_order_does_not_change_policy_or_outcomes():
    models = _models()
    scenarios = tc.draw_scenarios(2, seed=20261108)
    first = _policy(models)
    second = _policy(tuple(reversed(models)))
    first_outcome = tc.simulate_two_endpoints(
        scenarios, first, seed=[20261108, 7], stop_endpoint="either", budget=2
    )
    second_outcome = tc.simulate_two_endpoints(
        scenarios, second, seed=[20261108, 7], stop_endpoint="either", budget=2
    )
    np.testing.assert_array_equal(first_outcome["ttp_baseline"], second_outcome["ttp_baseline"])
    np.testing.assert_array_equal(first_outcome["ttp_rano"], second_outcome["ttp_rano"])
    np.testing.assert_array_equal(first_outcome["courses"], second_outcome["courses"])
