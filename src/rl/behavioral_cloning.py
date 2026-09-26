#!/usr/bin/env python3
"""
src/rl/behavioral_cloning.py
============================
Train a policy network via supervised imitation learning (behavioral cloning)
from the model-based OptimalController (Phase 2).

Literature basis:
- arXiv:2609.06667: Behavioral cloning from a Pontryagin-derived teacher achieves
  100% cure rate on adaptive tumor therapy, whereas actor-critic exploration collapses.
- Crucial finding: Use the BC policy directly; do NOT fine-tune with SAC/critic-driven RL.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
from torch.distributions import Categorical

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "src" / "rl"))

from optimal_controller import OptimalController


# --------------------------------------------------------------------------- #
# Policy Network Architectures
# --------------------------------------------------------------------------- #
class PolicyNet(nn.Module):
    """Binary action policy network for 2-state Pontryagin ODE [M_s, M_r]."""

    def __init__(self, state_dim: int = 2, hidden: int = 64):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, 1),
            nn.Sigmoid(),
        )

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        return self.net(state)


class PolicyNetwork(nn.Module):
    """
    4-action multi-modality policy network directly compatible with
    GbmTherapyEnv and script 59 sensitivity analysis.
    
    Observation (5): [norm_vol, u_max, day_frac, chemo_tox, rad_tox]
    Actions (4): 0=Holiday, 1=TMZ, 2=RT, 3=Combo
    """

    def __init__(self, obs_dim: int = 5, n_actions: int = 4, hidden: int = 64):
        super().__init__()
        self.obs_dim = obs_dim
        self.n_actions = n_actions
        self.net = nn.Sequential(
            nn.Linear(obs_dim, hidden),
            nn.ReLU(),
            nn.Linear(hidden, hidden),
            nn.ReLU(),
            nn.Linear(hidden, n_actions),
        )

    def forward(self, x: torch.Tensor) -> Categorical:
        logits = self.net(x)
        return Categorical(logits=logits)


# --------------------------------------------------------------------------- #
# Demonstration Generation
# --------------------------------------------------------------------------- #
def generate_ode_demonstrations(
    rho_s: float = 0.0195,
    n_trajectories: int = 200,
    n_steps: int = 90,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """Generate (state, action) pairs from OptimalController on ODE system."""
    rng = np.random.default_rng(seed)
    controller = OptimalController(rho_s=rho_s)
    states = []
    actions = []

    for _ in range(n_trajectories):
        M_s0 = float(rng.uniform(0.05, 0.35))
        M_r0 = float(rng.uniform(1e-5, 5e-3))
        ms_traj, mr_traj, u_traj = controller.rollout(M_s0, M_r0, n_steps=n_steps)

        for i in range(len(u_traj)):
            states.append([ms_traj[i], mr_traj[i]])
            actions.append([u_traj[i]])

    return np.array(states, dtype=np.float32), np.array(actions, dtype=np.float32)


def generate_env_demonstrations(
    rho_list: Optional[List[float]] = None,
    n_samples_per_rho: int = 1500,
    seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate diverse (obs, action) pairs covering GbmTherapyEnv observation space:
    obs = [norm_vol, u_max, day_frac, chemo_tox, rad_tox]
    action in {0, 1, 2, 3}
    
    Generates physically consistent sequential rollouts using FastPDESolver and OptimalController.
    """
    import importlib
    m = importlib.import_module("59_sensitivity_analysis")
    FastPDESolver = m.FastPDESolver
    GbmTherapyEnv = m.GbmTherapyEnv

    rng = np.random.default_rng(seed)
    if rho_list is None:
        rho_list = [0.0195]

    states = []
    actions = []

    # Number of 90-day trajectories per rho
    n_trajs_per_rho = max(5, int(np.ceil(n_samples_per_rho / 90)))

    for rho in rho_list:
        controller = OptimalController(rho_s=rho)
        for _ in range(n_trajs_per_rho):
            D_w = float(rng.uniform(0.001, 0.008))
            alpha = float(rng.uniform(0.5, 1.5))

            solver = FastPDESolver(
                grid_size=(16, 16, 16),
                dt_pde=0.5,
                rho=rho,
                D_white=D_w,
                alpha_sens=alpha,
                is_training=False,
            )
            env = GbmTherapyEnv(solver)
            obs, _ = env.reset()
            env.solver.u *= 0.1
            env.solver.initial_volume = float(env.solver.u.sum() * env.solver.dx**3)
            min_vol = env.solver.initial_volume

            for _ in range(env.max_steps):
                current_vol = float(env.solver.u.sum() * env.solver.dx**3)
                if current_vol < min_vol:
                    min_vol = current_vol
                min_norm = min_vol / max(env.solver.initial_volume, 1e-6)

                act = controller.get_env_action(obs, rho=rho, min_norm_vol=min_norm)
                states.append(obs.copy())
                actions.append(act)

                obs, _, term, _, _ = env.step(act)
                if term:
                    break

    return np.array(states, dtype=np.float32), np.array(actions, dtype=np.int64)


# --------------------------------------------------------------------------- #
# Behavioral Cloning Training Loops
# --------------------------------------------------------------------------- #
def train_bc_ode(
    rho_s: float = 0.0195,
    epochs: int = 200,
    lr: float = 1e-3,
    output_path: Optional[Path] = None,
) -> PolicyNet:
    """Train binary PolicyNet on 2-state ODE demonstrations."""
    states, actions = generate_ode_demonstrations(rho_s=rho_s, n_trajectories=200)
    print(f"[BC ODE] Generated {len(states)} demonstration pairs for rho={rho_s:.4f}")

    states_t = torch.FloatTensor(states)
    actions_t = torch.FloatTensor(actions)

    policy = PolicyNet(state_dim=2, hidden=64)
    optimizer = torch.optim.Adam(policy.parameters(), lr=lr)
    criterion = nn.BCELoss()

    for epoch in range(epochs):
        optimizer.zero_grad()
        pred = policy(states_t)
        loss = criterion(pred, actions_t)
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 50 == 0 or epoch == epochs - 1:
            acc = ((pred > 0.5).float() == actions_t).float().mean()
            print(f"  Epoch {epoch+1:3d}/{epochs}: loss={loss.item():.4f}, accuracy={acc.item()*100:.2f}%")

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(policy.state_dict(), str(output_path))
        print(f"[BC ODE] Saved model -> {output_path}")

    return policy


def train_bc_env(
    rho_list: Optional[List[float]] = None,
    epochs: int = 250,
    lr: float = 2e-3,
    output_path: Optional[Path] = None,
) -> PolicyNetwork:
    """Train 4-action PolicyNetwork on GbmTherapyEnv state demonstrations."""
    if rho_list is None:
        rho_list = [0.0195]

    states, actions = generate_env_demonstrations(rho_list=rho_list, n_samples_per_rho=2000)
    print(f"[BC Env] Generated {len(states)} state-action pairs across {len(rho_list)} rho values")

    states_t = torch.FloatTensor(states)
    actions_t = torch.LongTensor(actions)

    policy = PolicyNetwork(obs_dim=5, n_actions=4, hidden=64)
    optimizer = torch.optim.Adam(policy.parameters(), lr=lr, weight_decay=1e-5)
    criterion = nn.CrossEntropyLoss()

    batch_size = 256
    n_batches = int(np.ceil(len(states) / batch_size))

    for epoch in range(epochs):
        perm = torch.randperm(len(states_t))
        epoch_loss = 0.0

        for b in range(n_batches):
            idx = perm[b * batch_size : (b + 1) * batch_size]
            b_states = states_t[idx]
            b_actions = actions_t[idx]

            optimizer.zero_grad()
            logits = policy.net(b_states)
            loss = criterion(logits, b_actions)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        if (epoch + 1) % 50 == 0 or epoch == epochs - 1:
            with torch.no_grad():
                logits = policy.net(states_t)
                preds = logits.argmax(dim=-1)
                acc = (preds == actions_t).float().mean()
            print(f"  Epoch {epoch+1:3d}/{epochs}: loss={epoch_loss/n_batches:.4f}, acc={acc.item()*100:.2f}%")

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        torch.save(policy.state_dict(), str(output_path))
        print(f"[BC Env] Saved model -> {output_path}")

    return policy


def main():
    parser = argparse.ArgumentParser(description="Train Behavioral Cloning Policy from Optimal Controller")
    parser.add_argument("--mode", choices=["env", "ode"], default="env")
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--output", type=str, default="output/bc_policy.pt")
    args = parser.parse_args()

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = PROJECT_ROOT / out_path

    if args.mode == "env":
        policy = train_bc_env(epochs=args.epochs, output_path=out_path)
    else:
        policy = train_bc_ode(epochs=args.epochs, output_path=out_path)


if __name__ == "__main__":
    main()
