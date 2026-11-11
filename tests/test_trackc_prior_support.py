"""Prior support and propagation checks (amendment rev 3, 2026-11-10).

Outcome-blind. No policy is simulated and no utility is computed.
"""
import importlib.util
import inspect
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

import trackc_integrated_sampler as smp
import trackc_sim as tc
from trackc_filter_priors import PF_PRIOR_INTEGRATED_V3 as V3, PF_PRIOR_NARROW_V2 as NARROW

ROOT = Path(__file__).resolve().parent.parent
CONFIG = json.loads((ROOT / "configs" / "trackABC_integrated_v2.json").read_text(encoding="utf-8"))
PRIOR = CONFIG["particle_filter_prior"]


def _load(name, filename):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, ROOT / "src" / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def benchmark():
    return _load("trackc_policy_benchmark_prior_test", "211_trackc_policy_benchmark.py")


@pytest.fixture(scope="module")
def mpc():
    return _load("trackc_robust_belief_mpc_prior_test", "208_trackc_robust_belief_mpc.py")


def _axis(family, name):
    for axis, lo, hi in smp.WORLD_AXIS_BOUNDS[family]:
        if axis == name:
            return lo, hi
    raise KeyError(name)


def generator_support():
    lo_r, hi_r = _axis("single", "log_r_shift")
    lo_f, hi_f = _axis("single", "log_f0_shift")
    lo_c, hi_c = _axis("single", "B0_center")
    m_lo = tc.RANGES["m"][0] * min(math.exp(_axis("single", "log_m_scale")[0]),
                                   math.exp(_axis("mixture", "log_m_scale_non")[0]))
    return {
        "r": (tc.RANGES["r"][0] * math.exp(lo_r), tc.RANGES["r"][1] * math.exp(hi_r)),
        "m": (m_lo, tc.RANGES["m"][1]),
        "f0": (tc.RANGES["f0"][0] * math.exp(lo_f), min(tc.RANGES["f0"][1] * math.exp(hi_f), 0.9)),
        "B0": (lo_c - smp.B0_CENTER_HALF_WIDTH, hi_c + smp.B0_CENTER_HALF_WIDTH),
        "c": tc.RANGES["c"][:2],
        "eps": tc.RANGES["eps"][:2],
    }


# ---------------------------------------------------------------- support
@pytest.mark.parametrize("key", tc.PARAMS)
def test_v3_covers_generator_support(key):
    lo, hi = generator_support()[key]
    assert V3[key][0] <= lo * (1 + 1e-12) and hi <= V3[key][1] * (1 + 1e-12)


def test_narrow_differs_from_v3_only_in_r_and_b0(mpc):
    assert {k for k in tc.PARAMS if NARROW[k] != V3[k]} == {"r", "B0"}
    assert NARROW == mpc.INTEGRATED_BELIEF_RANGES
    assert NARROW != {k: tc.RANGES[k] for k in tc.PARAMS}


def test_legacy_defaults_unchanged(mpc):
    assert tc.RANGES["r"] == (0.001, 0.05, "log") and tc.RANGES["B0"] == (0.05, 0.5, "lin")
    assert tc.RANGES["m"] == (6.0, 20.0, "log") and tc.RANGES["f0"] == (0.001, 0.1, "log")
    assert mpc.INTEGRATED_BELIEF_RANGES["r"] == (0.001, 0.05, "log")


def test_config_prior_is_v3():
    assert PRIOR["version"] == "PF_PRIOR_INTEGRATED_V3"
    assert tuple(PRIOR["r_day_minus1_log_uniform"]) == V3["r"][:2]
    assert tuple(PRIOR["B0_uniform"]) == V3["B0"][:2]


def test_sampled_patients_inside_v3():
    for family in ("single", "mixture"):
        for seed in range(40):
            world = smp.draw_world_parameters(family, 9_000_000 + seed)
            cohort = smp.draw_patient_cohort(world, 12, 9_100_000 + seed)
            for key in tc.PARAMS:
                assert V3[key][0] <= cohort[key].min() and cohort[key].max() <= V3[key][1]


def test_fraction_outside_narrow_is_recorded():
    """Observed counts for one fixed synthetic batch; not population rates."""
    out, total = {"r": 0, "B0": 0}, 0
    for family in ("single", "mixture"):
        for seed in range(100):
            world = smp.draw_world_parameters(family, 9_200_000 + seed)
            cohort = smp.draw_patient_cohort(world, 12, 9_300_000 + seed)
            total += 12
            for key in out:
                v = cohort[key]
                out[key] += int(np.sum((v < NARROW[key][0]) | (v > NARROW[key][1])))
    assert (total, out["r"], out["B0"]) == (2400, 246, 324)


# ---------------------------------------------------------------- propagation
def _filter_arms(catalog):
    for name, spec in catalog.items():
        policy = spec["factory"]()
        if getattr(policy, "uses_filter", False):
            yield name, policy


def test_every_baseline_filter_arm_gets_v3(benchmark):
    arms = list(_filter_arms(benchmark.baseline_policy_catalog(benchmark.PolicyBenchmarkSettings(), V3)))
    assert len(arms) >= 4
    for name, policy in arms:
        assert policy.ranges is V3, name


def test_baseline_filter_arms_keep_historical_prior_by_default(benchmark):
    for name, policy in _filter_arms(benchmark.baseline_policy_catalog(benchmark.PolicyBenchmarkSettings())):
        assert policy.ranges in (tc.RANGES, tc.RANGES_WIDE_M), name


def test_filter_arm_reset_builds_filter_with_v3(benchmark):
    catalog = benchmark.baseline_policy_catalog(benchmark.PolicyBenchmarkSettings(mpc_particles=8, prediction_particles=8), V3)
    for name, policy in _filter_arms(catalog):
        policy.reset({"sigma": np.array([0.1, 0.2])}, np.random.default_rng(0))
        assert policy.pf.ranges is V3, name
        for i, key in enumerate(tc.PARAMS):
            values = policy.pf.theta[key]
            assert V3[key][0] <= values.min() and values.max() <= V3[key][1], (name, key)


def test_mpc_arms_take_ranges_only_from_common_mpc(benchmark):
    source = inspect.getsource(benchmark.run_shared_policy_world)
    assert source.count('"ranges":') == 1
    assert '"ranges": mpc.INTEGRATED_BELIEF_RANGES if filter_ranges is None else filter_ranges' in source
    assert "else baseline_policy_catalog(settings, filter_ranges)" in source
    block = source[source.index("def make_signal_policies"):source.index("return policies_for_signals")]
    constructors = block.count("MPC(")
    assert constructors == 6
    assert block.count("**deepcopy(common_mpc)") == constructors


def test_mpc_classes_store_supplied_ranges(mpc):
    ambiguity = _load("trackc_robust_signal_ambiguity_prior_test", "216_trackc_robust_signal_ambiguity.py")
    assert mpc.RobustBeliefMPC(ranges=V3).ranges is V3
    assert "ranges" in inspect.signature(ambiguity.SignalAmbiguityBeliefMPC.__init__).parameters \
        or "kwargs" in str(inspect.signature(ambiguity.SignalAmbiguityBeliefMPC.__init__))


def test_world_runner_passes_v3():
    from trackc_world_runner import prepare_synthetic_world
    from trackc_integrated_signal_model import SignalModelParameters
    import trackc_seed_registry as registry
    import trackc_world_schedule as schedule
    from types import SimpleNamespace
    w = {role: 1 + i for i, role in enumerate(schedule.SCENARIO_ROLES)}
    w.update({role: 1002 + i for i, role in enumerate(schedule.POLICY_ROLES)})
    w["world_id"] = "a"
    prepared = prepare_synthetic_world(
        w, family="single", patients_per_cohort=2,
        signal_parameters=SignalModelParameters(0.7, 0.2, 0.3, 0.1),
        model_parameters=(SignalModelParameters(0.7, 0.2, 0.3, 0.1),),
        settings=SimpleNamespace(stop_endpoint="observed_two_scan"),
        entries=[registry._entry("reserved", "dev", "reserved", [(1, 1000)], [], purpose="stage2_dev")],
    )
    assert prepared["benchmark_kwargs"]["filter_ranges"] is V3
