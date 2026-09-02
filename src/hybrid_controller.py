#!/usr/bin/env python3
"""Trust-gated MPC-RL hybrid therapy controller.

Two architectures, both operating as policies inside the 3D anisotropic PDE
environment (FastPDESolver + GbmTherapyEnv from 64_virtual_cohort_simulation):

  Architecture A (hard gate):
      |benefit| > tau  ->  trust the robust-MPC decision
      |benefit| <= tau ->  delegate to the RL policy

  Architecture B (trust-conditioned RL):
      The RL policy's observation is augmented with the normalized trust score
      from the robust MPC. The policy learns when to emulate MPC (high trust)
      and when to deviate (low trust).

The MPC controller used for trust extraction is the corrected
RobustMPCController (see trust_signal module). The RL policy is a small MLP
(REINFORCE), compatible with the env's action space {0: Rest, 1: TMZ, 2: RT,
3: Combo}. MPC's binary decision is mapped to {Rest, Combo} for a fair match
with the Stupp active phase.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import trust_signal  # noqa: E402
from trust_signal import (  # noqa: E402
    CalibrationTrust,
    RobustMPCController,
    compute_trust_signal,
)

try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.distributions import Categorical

    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False
    torch = nn = optim = Categorical = None


# --------------------------------------------------------------------------- #
# RL policy (shared by Arch A delegate and Arch B)
# --------------------------------------------------------------------------- #
class TrustConditionedPolicy(nn.Module):
    """MLP policy. obs_dim = 5 (env) or 6 (env + trust)."""

    def __init__(self, obs_dim: int = 6, n_actions: int = 4, hidden: int = 64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(obs_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, n_actions),
        )
        self._init_action_bias()

    def _init_action_bias(self):
        with torch.no_grad():
            torch.nn.init.constant_(self.net[-1].bias, 0.0)
            self.net[-1].bias.data[2] = 0.5  # RT preference
            self.net[-1].bias.data[3] = 1.0  # Combo preference

    def forward(self, x):
        logits = self.net(x)
        return Categorical(logits=logits)


def train_reinforce(
    env,
    obs_dim: int = 6,
    episodes: int = 60,
    lr: float = 1e-2,
    gamma: float = 0.99,
    entropy_coef: float = 0.01,
    seed: Optional[int] = 0,
) -> Optional[TrustConditionedPolicy]:
    """Train the trust-conditioned policy via REINFORCE.

    Keeps the policy that achieves the best mean episode return so training
    noise doesn't leave a degenerate (always-rest) policy.
    """
    if not HAS_TORCH:
        return None
    if seed is not None:
        torch.manual_seed(seed)
        np.random.seed(seed)

    policy = TrustConditionedPolicy(obs_dim=obs_dim)
    optimizer = optim.Adam(policy.parameters(), lr=lr)
    best_policy_state = None
    best_mean_return = -float("inf")

    for ep in range(episodes):
        obs, _ = env.reset()
        log_probs, rewards, entropies = [], [], []
        for _ in range(env.max_steps):
            obs_t = torch.FloatTensor(np.asarray(obs, dtype=np.float32)).unsqueeze(0)
            dist = policy(obs_t)
            action = dist.sample()
            log_probs.append(dist.log_prob(action))
            entropies.append(dist.entropy())
            obs, reward, terminated, truncated, _ = env.step(action.item())
            rewards.append(reward)
            if terminated or truncated:
                break
        returns = []
        R = 0.0
        for r in reversed(rewards):
            R = r + gamma * R
            returns.insert(0, R)
        returns = torch.tensor(returns, dtype=torch.float32)
        returns = (returns - returns.mean()) / (returns.std() + 1e-8)
        loss = -torch.stack(log_probs).sum() * returns.sum() - entropy_coef * torch.stack(entropies).sum()
        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(policy.parameters(), max_norm=1.0)
        optimizer.step()

        ep_return = sum(rewards)
        if ep_return > best_mean_return:
            best_mean_return = ep_return
            best_policy_state = {k: v.clone() for k, v in policy.state_dict().items()}

    if best_policy_state is not None:
        policy.load_state_dict(best_policy_state)
    return policy


def rl_choose_action(policy, obs: np.ndarray) -> int:
    """Deterministic argmax action from a trained policy."""
    if policy is None:
        return 3  # fallback: combo
    obs_t = torch.FloatTensor(np.asarray(obs, dtype=np.float32)).unsqueeze(0)
    with torch.no_grad():
        return int(policy(obs_t).probs.argmax().item())


# --------------------------------------------------------------------------- #
# MPC-as-policy wrapper inside the 3D env
# --------------------------------------------------------------------------- #
class MPCAsPolicy:
    """Wraps RobustMPCController to act in the 3D PDE env.

    MPC reads the current scalar volume from the env, computes its decision,
    and maps binary dose -> action (Combo when dosing, Rest when holding).
    """

    ACTION_DOSE = 3  # Combo (TMZ + RT), matches Stupp active phase
    ACTION_HOLD = 0  # Rest

    def __init__(
        self,
        rho: float,
        D: float,
        target_volume_mm3: float,
        controller: Optional[RobustMPCController] = None,
    ):
        self.rho = rho
        self.D = D
        self.target_volume = target_volume_mm3
        self.controller = controller or RobustMPCController(seed=0)
        self.step = 0
        self.last_predicted_volume: Optional[float] = None

    def reset(self):
        self.step = 0
        self.last_predicted_volume = None

    def choose_action(self, env, obs: np.ndarray) -> Tuple[int, Dict[str, Any]]:
        volume = float(env.solver.u.sum() * env.solver.dx**3)
        dose, _ = self.controller.optimize_control(
            current_volume_mm3=volume,
            rho_nominal=self.rho,
            D_nominal=self.D,
            target_volume=self.target_volume,
            step=self.step,
        )
        # Predicted volume at the end of the current horizon under the chosen
        # decision (used for calibration trust).
        self.last_predicted_volume = trust_signal._predict_volume_corrected(
            volume, self.rho, self.D, self.controller.horizon, self.step,
            drug_on=dose >= 0.5, target_volume=self.target_volume,
            w_tumor=self.controller.w_tumor, w_drug=self.controller.w_drug,
        )
        self.step += 1
        action = self.ACTION_DOSE if dose >= 0.5 else self.ACTION_HOLD
        return action, {"mpc_dose": dose, "predicted_volume": self.last_predicted_volume}


# --------------------------------------------------------------------------- #
# Architecture A: hard-gate hybrid
# --------------------------------------------------------------------------- #
class HybridGateA:
    """Hard-switch between MPC (confident) and RL (uncertain region).

    trust_mode:
      - "margin": |benefit_of_dosing| > tau -> MPC, else RL (parameter-ensemble)
      - "calibration": MPC calibration error < threshold -> MPC, else RL
    """

    def __init__(
        self,
        rho: float,
        D: float,
        target_volume_mm3: float,
        policy=None,
        tau: float = 0.05,
        n_samples: int = trust_signal.DEFAULT_N_SAMPLES,
        seed: Optional[int] = None,
        trust_mode: str = "margin",
        cal_threshold: float = 0.5,
    ):
        self.rho = rho
        self.D = D
        self.target_volume = target_volume_mm3
        self.policy = policy
        self.tau = tau
        self.n_samples = n_samples
        self.trust_mode = trust_mode
        self.cal_threshold = cal_threshold
        self.controller = RobustMPCController(seed=seed or 0)
        self.mpc = MPCAsPolicy(rho=rho, D=D, target_volume_mm3=target_volume_mm3, controller=self.controller)
        self.calibration = CalibrationTrust()
        self.step = 0
        self.trust_history: list[float] = []

    def reset(self):
        self.step = 0
        self.trust_history = []
        self.mpc.reset()
        self.calibration = CalibrationTrust()

    def choose_action(self, env, obs: np.ndarray) -> Tuple[int, Dict[str, Any]]:
        volume = float(env.solver.u.sum() * env.solver.dx**3)
        # Always compute MPC decision + prediction so calibration can be tracked.
        mpc_action, mpc_diag = self.mpc.choose_action(env, obs)
        pred = mpc_diag.get("predicted_volume")
        if pred is not None and self.calibration.n_observations > 0:
            # Record previous prediction against current observed volume.
            pass

        if self.trust_mode == "calibration":
            # Record prediction from the PREVIOUS step against the CURRENT volume.
            # mpc.last_predicted_volume is from the step just completed.
            if self.mpc.last_predicted_volume is not None:
                self.calibration.record(self.mpc.last_predicted_volume, volume)
            trust = self.calibration.trust()
            confident = trust >= self.cal_threshold
        else:
            sig = compute_trust_signal(
                self.controller,
                current_volume_mm3=volume,
                rho_nominal=self.rho,
                D_nominal=self.D,
                target_volume=self.target_volume,
                step=self.step,
                n_samples=self.n_samples,
                seed=self.step % (2**31 - 1),
            )
            trust = sig["trust"]
            confident = abs(sig["benefit_mean"]) > self.tau

        self.step += 1
        self.trust_history.append(trust)

        if confident:
            action = mpc_action
            delegate = False
        else:
            action = rl_choose_action(self.policy, np.asarray(obs, dtype=np.float32))
            delegate = True
        return action, {"trust": trust, "delegate": delegate, "mpc_action": mpc_action}


# --------------------------------------------------------------------------- #
# Architecture B: trust-conditioned RL
# --------------------------------------------------------------------------- #
class HybridGateB:
    """RL policy conditioned on the MPC trust score as an extra input feature."""

    OBS_DIM = 6  # 5 env obs + 1 trust

    def __init__(
        self,
        rho: float,
        D: float,
        target_volume_mm3: float,
        policy: Optional[TrustConditionedPolicy] = None,
        n_samples: int = trust_signal.DEFAULT_N_SAMPLES,
        seed: Optional[int] = None,
        trust_mode: str = "margin",
    ):
        self.rho = rho
        self.D = D
        self.target_volume = target_volume_mm3
        self.policy = policy
        self.n_samples = n_samples
        self.trust_mode = trust_mode
        self.controller = RobustMPCController(seed=seed or 0)
        self.mpc = MPCAsPolicy(rho=rho, D=D, target_volume_mm3=target_volume_mm3, controller=self.controller)
        self.calibration = CalibrationTrust()
        self.step = 0
        self.trust_history: list[float] = []

    def reset(self):
        self.step = 0
        self.trust_history = []
        self.mpc.reset()
        self.calibration = CalibrationTrust()

    def _augment_obs(self, obs: np.ndarray, trust: float) -> np.ndarray:
        return np.concatenate([np.asarray(obs, dtype=np.float32), [np.clip(trust, 0, 1)]]).astype(np.float32)

    def choose_action(self, env, obs: np.ndarray) -> Tuple[int, Dict[str, Any]]:
        volume = float(env.solver.u.sum() * env.solver.dx**3)
        _, mpc_diag = self.mpc.choose_action(env, obs)

        if self.trust_mode == "calibration":
            if self.mpc.last_predicted_volume is not None:
                self.calibration.record(self.mpc.last_predicted_volume, volume)
            trust = self.calibration.trust()
        else:
            sig = compute_trust_signal(
                self.controller,
                current_volume_mm3=volume,
                rho_nominal=self.rho,
                D_nominal=self.D,
                target_volume=self.target_volume,
                step=self.step,
                n_samples=self.n_samples,
                seed=self.step % (2**31 - 1),
            )
            trust = sig["trust"]

        self.step += 1
        self.trust_history.append(trust)
        aug = self._augment_obs(obs, trust)
        action = rl_choose_action(self.policy, aug)
        return action, {"trust": trust, "benefit_mean": mpc_diag.get("benefit_mean")}


# --------------------------------------------------------------------------- #
# Env observation wrapper that injects trust during RL training
# --------------------------------------------------------------------------- #
class TrustObsWrapper:
    """Wraps GbmTherapyEnv to append a trust feature to observations."""

    def __init__(self, env, trust_fn):
        self.env = env
        self.trust_fn = trust_fn
        self.max_steps = env.max_steps

    def reset(self, seed=None, options=None):
        obs, info = self.env.reset(seed=seed, options=options)
        return self._wrap(obs), info

    def step(self, action):
        obs, reward, terminated, truncated, info = self.env.step(action)
        return self._wrap(obs), reward, terminated, truncated, info

    def _wrap(self, obs):
        trust = self.trust_fn()
        return np.concatenate([np.asarray(obs, dtype=np.float32), [np.clip(trust, 0, 1)]]).astype(np.float32)


def make_trust_fn(controller, rho, D, target_volume, get_volume, n_samples=trust_signal.DEFAULT_N_SAMPLES, step=0):
    def _trust():
        nonlocal step
        sig = compute_trust_signal(
            controller, get_volume(), rho, D, target_volume, step=step,
            n_samples=n_samples, seed=step % (2**31 - 1),
        )
        step += 1
        return sig["trust"]
    return _trust


if __name__ == "__main__":
    print("hybrid_controller imported OK")