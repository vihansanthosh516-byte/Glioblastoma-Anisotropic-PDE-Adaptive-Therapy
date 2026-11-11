"""Synthetic mechanics tests for the integrated driver. No policy utility is computed."""
import functools
import importlib.util
import inspect
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import trackc_integrated_driver as drv
import trackc_seed_registry as registry
from trackc_filter_priors import PF_PRIOR_INTEGRATED_V3
from trackc_integrated_signal_model import SignalModelParameters

ROOT = Path(__file__).resolve().parent.parent
SIG = SignalModelParameters(0.7, 0.2, 0.3, 0.1)
SETTINGS = SimpleNamespace(stop_endpoint="observed_two_scan")


def entries():
    return [registry._entry("reserved", "dev", "reserved", [(1, 9999)], [], purpose="stage2_dev")]


def stub_benchmark(**kwargs):
    """Records inputs only. Simulates nothing and scores nothing."""
    ev = kwargs["evaluation_scenarios"]
    return {"world_id": kwargs["world_id"], "stub": True,
            "n_eval": len(ev["sigma"]),
            "same_patients_for_all_arms": True,
            "signals_sha": hash(kwargs["evaluation_observations"].mri_score.tobytes())}


def run(worlds, **kw):
    kw.setdefault("benchmark", stub_benchmark)
    return drv.run_driver(worlds, mode=kw.pop("mode", "mechanics"), family=kw.pop("family", "single"),
                          patients_per_cohort=3, signal_parameters=SIG, model_parameters=(SIG,),
                          settings=SETTINGS, entries=entries(), **kw)


def test_protocol_loaded_from_config():
    p = drv.load_protocol()
    assert p.prior_version == "PF_PRIOR_INTEGRATED_V3"
    assert p.precision_half_width_days == 14.0
    assert p.stop_endpoint == "observed_two_scan"
    assert p.scan_interval_days == 84 and p.scan_every_cycles == 3
    assert p.open_gates, "config still lists open evaluation gates"


def test_evaluation_mode_refused_while_gates_open():
    worlds = drv.build_world_schedule(2, scenario_base=100, policy_base=20100, prefix="w")
    with pytest.raises(drv.EvaluationBlocked, match="open gates"):
        run(worlds, mode="evaluation")


def test_evaluation_refused_if_precision_target_changed():
    p = drv.load_protocol()
    altered = drv.DriverProtocol(**{**p.__dict__, "open_gates": (), "precision_half_width_days": 30.0})
    with pytest.raises(drv.EvaluationBlocked, match="14-day"):
        drv.check_evaluation_gates(altered)


def test_mechanics_refuses_real_benchmark():
    spec = importlib.util.spec_from_file_location("bench_drv_test", ROOT / "src" / "211_trackc_policy_benchmark.py")
    mod = importlib.util.module_from_spec(spec); sys.modules["bench_drv_test"] = mod; spec.loader.exec_module(mod)
    worlds = drv.build_world_schedule(1, scenario_base=100, policy_base=20100, prefix="w")
    for fn in (mod.run_shared_policy_world, functools.partial(mod.run_shared_policy_world),
               functools.wraps(mod.run_shared_policy_world)(lambda **k: None)):
        with pytest.raises(PermissionError, match="refuses the real policy benchmark"):
            run(worlds, benchmark=fn)


def test_mechanics_run_is_complete_reproducible_and_uses_v3_and_controls():
    worlds = drv.build_world_schedule(3, scenario_base=200, policy_base=20200, prefix="m")
    a, b = run(worlds), run(worlds)
    assert a["n_worlds"] == 3 and [r["world_id"] for r in a["records"]] == ["m0000", "m0001", "m0002"]
    assert a["inputs_sha256"] == b["inputs_sha256"]
    assert a["policy_outcomes_computed"] is False
    for arm in drv.PRIMARY_ARMS + drv.REQUIRED_CONTROL_ARMS:
        assert arm in a["arm_plan"]
    assert all(r["stub"] for r in a["records"])


def test_schedule_seeds_distinct_and_validated():
    worlds = drv.build_world_schedule(4, scenario_base=500, policy_base=20500, prefix="s")
    seeds = [w[k] for w in worlds for k in w if k != "world_id"]
    assert len(seeds) == len(set(seeds)) == 44
    dup = [dict(worlds[0]), dict(worlds[0], world_id="dup")]
    with pytest.raises(ValueError):
        run(dup)
    with pytest.raises(ValueError, match="overlap"):
        drv.build_world_schedule(4, scenario_base=500, policy_base=510, prefix="s")


def test_arm_plan_check_rejects_missing_controls_or_prior():
    good = {"filter_ranges": PF_PRIOR_INTEGRATED_V3, "include_ambiguity_controller": True,
            "include_single_signal_controls": True, "primary_contrast": "minimax"}
    drv.check_arm_plan(good)
    for bad in ({**good, "include_single_signal_controls": False}, {**good, "filter_ranges": None},
                {**good, "primary_contrast": "loo"}):
        with pytest.raises(ValueError):
            drv.check_arm_plan(bad)


def test_signals_are_simulated_from_declared_assumptions_only():
    src = inspect.getsource(drv)
    for forbidden in ("Duke", "GSE117891", "output/trackA", "nibabel", "read_csv"):
        assert forbidden not in src.split('"""', 2)[-1]
    worlds = drv.build_world_schedule(1, scenario_base=900, policy_base=20900, prefix="z")
    out = run(worlds)
    assert out["records"][0]["n_eval"] == 3


def test_single_signal_arm_restricts_likelihood_channels():
    spec = importlib.util.spec_from_file_location("amb_drv_test", ROOT / "src" / "216_trackc_robust_signal_ambiguity.py")
    mod = importlib.util.module_from_spec(spec); sys.modules["amb_drv_test"] = mod; spec.loader.exec_module(mod)
    from trackc_integrated_signal_model import SignalObservations
    obs = SignalObservations(np.array([True, False]), np.array([0.1, -0.2]))
    pol = mod.SignalAmbiguityBeliefMPC(obs, model_parameters=(SIG,), channels=("mri",), ranges=PF_PRIOR_INTEGRATED_V3,
                                       n_particles=8, prediction_particles=4)
    assert pol.signal_channels == ("mri",)
    with pytest.raises(ValueError):
        mod.SignalAmbiguityBeliefMPC(obs, model_parameters=(SIG,), channels=(), ranges=PF_PRIOR_INTEGRATED_V3)
