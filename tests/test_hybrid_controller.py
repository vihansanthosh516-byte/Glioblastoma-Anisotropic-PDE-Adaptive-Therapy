import numpy as np

from src.trust_signal import RobustMPCController, compute_trust_signal, benefit_distribution
from src.hybrid_controller import (
    MPCAsPolicy,
    HybridGateA,
    HybridGateB,
    TrustConditionedPolicy,
    rl_choose_action,
)


def test_benefit_distribution_shape():
    ctrl = RobustMPCController(seed=0)
    b = benefit_distribution(ctrl, 5000.0, 0.02, 0.013, 600.0, step=100, n_samples=32)
    assert b.shape == (32,)


def test_trust_bounded_0_to_1():
    ctrl = RobustMPCController(seed=0)
    sig = compute_trust_signal(ctrl, 5000.0, 0.02, 0.013, 600.0, step=100)
    assert 0.0 <= sig["trust"] <= 1.0
    assert sig["mpc_decision"] in (0.0, 1.0)


def test_high_rho_more_uncertain_lower_trust():
    ctrl = RobustMPCController(seed=0)
    lo = compute_trust_signal(ctrl, 5000.0, 0.005, 0.013, 600.0, step=100, n_samples=32)
    hi = compute_trust_signal(ctrl, 5000.0, 0.05, 0.013, 600.0, step=100, n_samples=32)
    assert lo["trust"] > hi["trust"]


def _make_solver_env():
    import importlib.util
    import sys
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    spec = importlib.util.spec_from_file_location(
        "vcs", ROOT / "src" / "64_virtual_cohort_simulation.py"
    )
    vcs = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(vcs)
    solver = vcs.FastPDESolver(grid_size=(16, 16, 16), is_training=False)
    env = vcs.GbmTherapyEnv(solver, vcs.RL_REWARD_WEIGHTS)
    obs, _ = env.reset()  # reset initializes solver.initial_volume
    return env


def test_mpc_policy_action_mapping():
    env = _make_solver_env()
    mpc = MPCAsPolicy(rho=0.02, D=0.013, target_volume_mm3=env.solver.initial_volume * 0.12)
    obs, _ = env.reset()
    action, diag = mpc.choose_action(env, obs)
    assert action in (0, 3)


def test_gate_a_always_returns_action():
    env = _make_solver_env()
    gate = HybridGateA(rho=0.02, D=0.013, target_volume_mm3=env.solver.initial_volume * 0.12, tau=0.05)
    obs, _ = env.reset()
    for _ in range(3):
        action, diag = gate.choose_action(env, obs)
        assert action in (0, 1, 2, 3)
        assert "trust" in diag
        obs, _, term, trunc, _ = env.step(action)
        if term or trunc:
            break
    assert len(gate.trust_history) >= 1
    assert all(0 <= t <= 1 for t in gate.trust_history)


def test_gate_b_obs_dim_and_action():
    env = _make_solver_env()
    gate = HybridGateB(rho=0.02, D=0.013, target_volume_mm3=env.solver.initial_volume * 0.12)
    obs, _ = env.reset()
    aug = gate._augment_obs(obs, 0.7)
    assert aug.shape == (6,)
    action, diag = gate.choose_action(env, obs)
    assert action in (0, 1, 2, 3)
    assert "trust" in diag


def test_policy_fallback_when_no_torch():
    # rl_choose_action with None policy must not crash (falls back to combo=3)
    action = rl_choose_action(None, np.zeros(6, dtype=np.float32))
    assert action == 3