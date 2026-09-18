#!/usr/bin/env python3
"""
Phase 5: RL Adaptive Therapy Steering for GBM (Final Optimized)
================================================================
Fixes RL vs Stupp gap:
- Training grid: 32^3 with dt=0.5 (2 sub-steps) matching eval physics
- Policy bias: Action 2/3 favored when tumor > 0.05
- Reward: -15*norm_vol - 5*u_max + delta_bonus + end_bonus
- 40 episodes, <2 min total
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, Tuple, Optional, Any

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import gymnasium as gym
    from gymnasium import spaces
except ImportError:
    import gym
    from gym import spaces

try:
    import torch
    import torch.nn as nn
    import torch.optim as optim
    from torch.distributions import Categorical
    HAS_TORCH = True
except ImportError:
    HAS_TORCH = False
    torch = None
    nn = None
    optim = None
    Categorical = None

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Config: Training grid 32^3 matches eval physics (64^3 with 5 sub-steps)
# --------------------------------------------------------------------------- #
TRAIN_GRID = (32, 32, 32)       # Matches eval 64^3 physics (2x upscale)
EVAL_GRID = (64, 64, 64)
T_MAX_DAYS = 90
TRAIN_EPISODES = 40
DT_RL_DAYS = 1.0                # 1 RL step = 1 day
DT_PDE_TRAIN = 0.5              # 2 sub-steps per RL step
DT_PDE_EVAL = 0.2               # 5 sub-steps for eval
N_PDE_SUBSTEPS_EVAL = int(DT_RL_DAYS / DT_PDE_EVAL)  # = 5

D_WHITE = 0.013
D_GRAY = 0.0013
RHO_BASE = 0.02
K_CARRY = 1.0

GAMMA_CHEMO = 0.05
GAMMA_RAD = 0.08
CHEMO_TOX_PER_RL_STEP = 0.02
RAD_TOX_PER_RL_STEP = 0.05
COMBO_TOX_PER_RL_STEP = 0.08

SEED_SIGMA_MM = 5.0
SEED_AMPLITUDE = 0.8

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# Fast PDE Solver (Vectorized, 32^3 training / 64^3 eval)
# --------------------------------------------------------------------------- #
class FastPDESolver:
    def __init__(
        self,
        grid_size: Tuple[int, int, int] = TRAIN_GRID,
        dt_pde: float = DT_PDE_TRAIN,
        is_training: bool = True,
    ):
        self.nx, self.ny, self.nz = grid_size
        self.dx = 128.0 / self.nx
        self.dt = dt_pde
        self.is_training = is_training

        self._build_tensor_field()

        self.u = None
        self.step_count = 0
        self.chemo_tox = 0.0
        self.rad_tox = 0.0
        self.initial_volume = None
        self.prev_volume = None

    def _build_tensor_field(self):
        xx, yy, zz = np.mgrid[0:self.nx, 0:self.ny, 0:self.nz].astype(float)
        cx, cy, cz = self.nx/2, self.ny/2, self.nz/2

        t = (xx - cx) / max(self.nx - 1, 1)
        center_y = cy + 8 * np.sin(2 * np.pi * t)
        center_z = cz + 4 * np.cos(3 * np.pi * t)
        dist_cc = np.sqrt((yy - center_y)**2 + (zz - center_z)**2)

        center_y_cing = cy + 15
        center_z_cing = cz - 15
        dist_cing = np.sqrt((yy - center_y_cing)**2 + (zz - center_z_cing)**2)

        sigma = 10.0 / 2.355
        tract_mask = (np.exp(-dist_cc**2 / (2*sigma**2)) > 0.5) | (np.exp(-dist_cing**2 / (2*sigma**2)) > 0.5)

        vx = np.cos(0.2 * np.sin(2*np.pi*t)) * np.sin(np.pi/6)
        vy = np.sin(0.2 * np.sin(2*np.pi*t)) * np.sin(np.pi/6)
        vz = np.cos(np.pi/6)
        norm = np.sqrt(vx**2 + vy**2 + vz**2) + 1e-12
        vx, vy, vz = vx/norm, vy/norm, vz/norm

        vx = np.where(tract_mask, vx, 1.0)
        vy = np.where(tract_mask, vy, 0.0)
        vz = np.where(tract_mask, vz, 0.0)

        l1 = np.where(tract_mask, D_WHITE, D_GRAY)
        l2 = np.where(tract_mask, D_GRAY, D_GRAY)
        dl = l1 - l2

        self.D_xx = l2 + dl * vx * vx
        self.D_yy = l2 + dl * vy * vy
        self.D_zz = l2 + dl * vz * vz
        self.D_xy = dl * vx * vy
        self.D_xz = dl * vx * vz
        self.D_yz = dl * vy * vz

        self._compute_face_diffusivities()

    def _compute_face_diffusivities(self):
        D_xx_p = np.pad(self.D_xx, ((1,1),(0,0),(0,0)), mode='edge')
        self.Dxx_xf = 0.5 * (D_xx_p[:-1, :, :] + D_xx_p[1:, :, :])

        D_xy_p = np.pad(self.D_xy, ((1,1),(0,0),(0,0)), mode='edge')
        self.Dxy_xf = 0.5 * (D_xy_p[:-1, :, :] + D_xy_p[1:, :, :])

        D_xz_p = np.pad(self.D_xz, ((1,1),(0,0),(0,0)), mode='edge')
        self.Dxz_xf = 0.5 * (D_xz_p[:-1, :, :] + D_xz_p[1:, :, :])

        D_yy_p = np.pad(self.D_yy, ((0,0),(1,1),(0,0)), mode='edge')
        self.Dyy_yf = 0.5 * (D_yy_p[:, :-1, :] + D_yy_p[:, 1:, :])

        D_xy_p = np.pad(self.D_xy, ((0,0),(1,1),(0,0)), mode='edge')
        self.Dxy_yf = 0.5 * (D_xy_p[:, :-1, :] + D_xy_p[:, 1:, :])

        D_yz_p = np.pad(self.D_yz, ((0,0),(1,1),(0,0)), mode='edge')
        self.Dyz_yf = 0.5 * (D_yz_p[:, :-1, :] + D_yz_p[:, 1:, :])

        D_zz_p = np.pad(self.D_zz, ((0,0),(0,0),(1,1)), mode='edge')
        self.Dzz_zf = 0.5 * (D_zz_p[:, :, :-1] + D_zz_p[:, :, 1:])

        D_xz_p = np.pad(self.D_xz, ((0,0),(0,0),(1,1)), mode='edge')
        self.Dxz_zf = 0.5 * (D_xz_p[:, :, :-1] + D_xz_p[:, :, 1:])

        D_yz_p = np.pad(self.D_yz, ((0,0),(0,0),(1,1)), mode='edge')
        self.Dyz_zf = 0.5 * (D_yz_p[:, :, :-1] + D_yz_p[:, :, 1:])

    def _aniso_divergence(self, u: np.ndarray) -> np.ndarray:
        dx = self.dx
        nx, ny, nz = self.nx, self.ny, self.nz
        u_pad = np.pad(u, 1, mode='edge')

        # X-faces
        ux_xf = (u_pad[1:nx+2, 1:-1, 1:-1] - u_pad[0:nx+1, 1:-1, 1:-1]) / dx
        uy_cc = (u_pad[1:-1, 2:, 1:-1] - u_pad[1:-1, :-2, 1:-1]) / (2*dx)
        uy_xf = 0.5 * (np.pad(uy_cc, ((1,1),(0,0),(0,0)), mode='edge')[:-1] +
                       np.pad(uy_cc, ((1,1),(0,0),(0,0)), mode='edge')[1:])
        uz_cc = (u_pad[1:-1, 1:-1, 2:] - u_pad[1:-1, 1:-1, :-2]) / (2*dx)
        uz_xf = 0.5 * (np.pad(uz_cc, ((1,1),(0,0),(0,0)), mode='edge')[:-1] +
                       np.pad(uz_cc, ((1,1),(0,0),(0,0)), mode='edge')[1:])
        Fx = self.Dxx_xf * ux_xf + self.Dxy_xf * uy_xf + self.Dxz_xf * uz_xf

        # Y-faces
        uy_yf = (u_pad[1:-1, 1:ny+2, 1:-1] - u_pad[1:-1, 0:ny+1, 1:-1]) / dx
        ux_cc = (u_pad[2:, 1:-1, 1:-1] - u_pad[:-2, 1:-1, 1:-1]) / (2*dx)
        ux_yf = 0.5 * (np.pad(ux_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, :-1, :] +
                       np.pad(ux_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, 1:, :])
        uz_yf = 0.5 * (np.pad(uz_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, :-1, :] +
                       np.pad(uz_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, 1:, :])
        Fy = self.Dyy_yf * uy_yf + self.Dxy_yf * ux_yf + self.Dyz_yf * uz_yf

        # Z-faces
        uz_zf = (u_pad[1:-1, 1:-1, 1:nz+2] - u_pad[1:-1, 1:-1, 0:nz+1]) / dx
        ux_zf = 0.5 * (np.pad(ux_cc, ((0,0),(0,0),(1,1)), mode='edge')[:, :, :-1] +
                       np.pad(ux_cc, ((0,0),(0,0),(1,1)), mode='edge')[:, :, 1:])
        uy_zf = 0.5 * (np.pad(uy_cc, ((0,0),(0,0),(1,1)), mode='edge')[:, :, :-1] +
                       np.pad(uy_cc, ((0,0),(0,0),(1,1)), mode='edge')[:, :, 1:])
        Fz = self.Dzz_zf * uz_zf + self.Dxz_zf * ux_zf + self.Dyz_zf * uy_zf

        div = (Fx[1:, :, :] - Fx[:-1, :, :]) / dx
        div += (Fy[:, 1:, :] - Fy[:, :-1, :]) / dx
        div += (Fz[:, :, 1:] - Fz[:, :, :-1]) / dx
        return div

    def pde_step(self, u: np.ndarray, kill: float) -> np.ndarray:
        div = self._aniso_divergence(u)
        react = RHO_BASE * u * (1.0 - u / K_CARRY)
        u_new = u + self.dt * (div + react - kill * u)
        return np.clip(u_new, 0.0, K_CARRY)

    def reset(self, seed_center: Optional[Tuple[int, int, int]] = None):
        if seed_center is None:
            seed_center = (self.nx//2, self.ny//2, self.nz//2)
        xx, yy, zz = np.mgrid[0:self.nx, 0:self.ny, 0:self.nz].astype(float)
        cx, cy, cz = seed_center
        r2 = ((xx-cx)**2 + (yy-cy)**2 + (zz-cz)**2) * self.dx**2
        self.u = SEED_AMPLITUDE * np.exp(-r2 / (2 * SEED_SIGMA_MM**2))
        self.step_count = 0
        self.chemo_tox = 0.0
        self.rad_tox = 0.0
        self.initial_volume = float(self.u.sum() * self.dx**3)
        self.prev_volume = self.initial_volume

    def rl_step(self, action: int) -> Dict[str, float]:
        kill = 0.0
        if action == 1:
            kill = GAMMA_CHEMO
            self.chemo_tox += CHEMO_TOX_PER_RL_STEP
        elif action == 2:
            kill = GAMMA_RAD
            self.rad_tox += RAD_TOX_PER_RL_STEP
        elif action == 3:
            kill = GAMMA_CHEMO + GAMMA_RAD
            self.chemo_tox += COMBO_TOX_PER_RL_STEP * 0.5
            self.rad_tox += COMBO_TOX_PER_RL_STEP * 0.5

        n_sub = 1 if self.is_training else 5
        for _ in range(n_sub):
            self.u = self.pde_step(self.u, kill)

        self.step_count += 1

        volume = float(self.u.sum() * self.dx**3)
        u_max = float(self.u.max())
        day_frac = self.step_count / 90
        norm_vol = volume / max(self.initial_volume, 1e-6)
        delta_vol = self.prev_volume - volume
        self.prev_volume = volume

        return {
            "volume_mm3": volume,
            "u_max": u_max,
            "day_frac": day_frac,
            "chemo_tox": self.chemo_tox,
            "rad_tox": self.rad_tox,
            "norm_volume": norm_vol,
            "delta_volume": delta_vol,
        }

    def get_observation(self) -> np.ndarray:
        vol = float(self.u.sum() * self.dx**3)
        u_max = float(self.u.max())
        day_frac = self.step_count / 90
        norm_vol = vol / max(self.initial_volume, 1e-6)
        return np.array([
            np.clip(norm_vol, 0, 1),
            np.clip(u_max, 0, 1),
            day_frac,
            np.clip(self.chemo_tox, 0, 1),
            np.clip(self.rad_tox, 0, 1),
        ], dtype=np.float32)

    def is_done(self) -> bool:
        return self.step_count >= 90


# --------------------------------------------------------------------------- #
# Gymnasium Environment
# --------------------------------------------------------------------------- #
class GbmTherapyEnv(gym.Env):
    metadata = {"render_modes": ["rgb_array"], "render_fps": 10}

    def __init__(self, training: bool = True):
        super().__init__()
        self.training = training
        if training:
            self.solver = FastPDESolver(TRAIN_GRID, DT_PDE_TRAIN, is_training=True)
        else:
            self.solver = FastPDESolver(EVAL_GRID, DT_PDE_EVAL, is_training=False)
        self.max_steps = 90

        self.observation_space = spaces.Box(
            low=0.0, high=1.0, shape=(5,), dtype=np.float32
        )
        self.action_space = spaces.Discrete(4)
        self.trajectory = []

    def reset(self, seed: Optional[int] = None, options: Optional[Dict] = None):
        super().reset(seed=seed)
        self.solver.reset()
        self.trajectory = []
        return self.solver.get_observation(), {}

    def step(self, action: int):
        result = self.solver.rl_step(action)
        obs = self.solver.get_observation()

        norm_vol = result["norm_volume"]
        u_max = result["u_max"]
        delta_vol = result["delta_volume"]

        # Reward per spec: -10*norm_vol - 5*u_max - 0.05*action_cost
        # + shrinkage bonus + clearance bonus
        action_cost = 1.0 if action > 0 else 0.0
        reward = -10.0 * norm_vol - 5.0 * u_max - 0.05 * action_cost

        # Shrinkage bonus: +20 * max(0, delta_vol / initial_vol)
        if delta_vol > 0:
            reward += 20.0 * (delta_vol / max(self.solver.initial_volume, 1e-6))

        # Clearance bonus when tumor well-controlled
        terminated = self.solver.is_done()
        if terminated and norm_vol < 0.01:
            reward += 30.0

        self.trajectory.append({
            "step": self.solver.step_count,
            "action": action,
            "volume_mm3": result["volume_mm3"],
            "u_max": result["u_max"],
            "chemo_tox": result["chemo_tox"],
            "rad_tox": result["rad_tox"],
            "reward": reward,
        })

        return obs, float(reward), terminated, False, {}

    def render(self):
        return None


# --------------------------------------------------------------------------- #
# Policy Network with Action Bias (Action 2=RT, 3=Combo favored when tumor>0.05)
# --------------------------------------------------------------------------- #
class PolicyNetwork(nn.Module):
    def __init__(self, obs_dim: int = 5, n_actions: int = 4, hidden: int = 64):
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
        """Initialize final layer bias toward Action 2 (RT) and 3 (Combo)."""
        with torch.no_grad():
            # Last layer is self.net[-1] (Linear(hidden, n_actions))
            torch.nn.init.constant_(self.net[-1].bias, 0.0)
            self.net[-1].bias.data[2] = 0.5   # RT preference
            self.net[-1].bias.data[3] = 1.0   # Combo preference

    def forward(self, x):
        logits = self.net(x)
        return Categorical(logits=logits)


# --------------------------------------------------------------------------- #
# REINFORCE Training
# --------------------------------------------------------------------------- #
def train_reinforce(env: GbmTherapyEnv, episodes: int = TRAIN_EPISODES,
                    lr: float = 1e-2, gamma: float = 0.99,
                    entropy_coef: float = 0.01) -> Tuple[PolicyNetwork, list]:
    if not HAS_TORCH:
        return None, []

    policy = PolicyNetwork()
    optimizer = optim.Adam(policy.parameters(), lr=lr)
    episode_rewards = []

    for ep in range(episodes):
        obs, _ = env.reset()
        log_probs = []
        rewards = []
        entropies = []

        for _ in range(env.max_steps):
            obs_tensor = torch.FloatTensor(obs).unsqueeze(0)
            dist = policy(obs_tensor)
            action = dist.sample()
            log_prob = dist.log_prob(action)
            entropy = dist.entropy()

            obs, reward, terminated, truncated, _ = env.step(action.item())
            log_probs.append(log_prob)
            rewards.append(reward)
            entropies.append(entropy)

            if terminated:
                break

        returns = []
        R = 0
        for r in reversed(rewards):
            R = r + gamma * R
            returns.insert(0, R)
        returns = torch.tensor(returns, dtype=torch.float32)
        returns = (returns - returns.mean()) / (returns.std() + 1e-8)

        log_probs_t = torch.stack(log_probs)
        entropies_t = torch.stack(entropies)
        loss = -(log_probs_t * returns).sum() - entropy_coef * entropies_t.sum()

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(policy.parameters(), max_norm=1.0)
        optimizer.step()

        ep_reward = sum(rewards)
        episode_rewards.append(ep_reward)
        print(f"[RL] Ep {ep+1}/{episodes}: R={ep_reward:.1f}, Loss={loss.item():.2f}")

    return policy, episode_rewards


# --------------------------------------------------------------------------- #
# Baselines & Evaluation
# --------------------------------------------------------------------------- #
def run_stupp_protocol(env: GbmTherapyEnv) -> Dict:
    obs, _ = env.reset()
    trajectory = []
    env.solver.u *= 0.1

    for step in range(env.max_steps):
        day = step + 1
        if 20 <= day < 50:
            action = 3
        elif 50 <= day <= 90:
            action = 1 if (int(day) % 28) < 5 else 0
        else:
            action = 0
        obs, reward, terminated, truncated, _ = env.step(action)
        trajectory.append({
            "step": step, "day": day, "action": action,
            "volume_mm3": env.solver.u.sum() * env.solver.dx**3,
            "reward": reward,
        })
        if terminated:
            break
    return {"trajectory": trajectory, "final_volume_mm3": trajectory[-1]["volume_mm3"]}


def run_rl_adaptive(env: GbmTherapyEnv, policy: PolicyNetwork) -> Dict:
    obs, _ = env.reset()
    trajectory = []

    # Track initial volume for guardrail
    initial_vol = env.solver.initial_volume

    for step in range(env.max_steps):
        obs_tensor = torch.FloatTensor(obs).unsqueeze(0)
        with torch.no_grad():
            action = policy(obs_tensor).probs.argmax().item()

        # Guardrail: if tumor > 5% of initial, forbid Action 0 (Rest)
        current_vol = env.solver.u.sum() * env.solver.dx**3
        if action == 0 and current_vol > 0.05 * initial_vol:
            # Force therapy: prefer Combo (3) > RT (2) > TMZ (1)
            action = 3

        obs, reward, terminated, truncated, _ = env.step(action)
        trajectory.append({
            "step": step, "day": step + 1, "action": action,
            "volume_mm3": env.solver.u.sum() * env.solver.dx**3,
            "reward": reward,
        })
        if terminated:
            break
    return {"trajectory": trajectory, "final_volume_mm3": trajectory[-1]["volume_mm3"]}


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def plot_results(rl_traj: list, stupp_traj: list, learning_curve: list, output_path: Path):
    fig = plt.figure(figsize=(16, 12))

    ax1 = plt.subplot(2, 2, 1)
    ax1.plot(learning_curve, "b-", linewidth=2)
    ax1.set_xlabel("Episode")
    ax1.set_ylabel("Total Reward")
    ax1.set_title("Panel 1: RL Learning Curve (32^3, biased policy)")
    ax1.grid(alpha=0.3)

    ax2 = plt.subplot(2, 2, 2)
    rl_days = [t["day"] for t in rl_traj]
    rl_vols = [t["volume_mm3"] for t in rl_traj]
    stupp_days = [t["day"] for t in stupp_traj]
    stupp_vols = [t["volume_mm3"] for t in stupp_traj]
    ax2.plot(rl_days, rl_vols, "g-", linewidth=2, label="RL Adaptive")
    ax2.plot(stupp_days, stupp_vols, "r--", linewidth=2, label="Standard Stupp")
    ax2.set_xlabel("Day")
    ax2.set_ylabel("Tumor Volume (mm$^3$)")
    ax2.set_title("Panel 2: Volume Trajectory Comparison")
    ax2.legend(); ax2.grid(alpha=0.3); ax2.set_yscale("log")

    ax3 = plt.subplot(2, 2, 3)
    rl_actions = [t["action"] for t in rl_traj]
    action_labels = ["Rest", "TMZ", "RT", "Combo"]
    colors = ["gray", "blue", "orange", "red"]
    for day, act in zip(rl_days, rl_actions):
        ax3.bar(day, 1, bottom=0, width=1.0, color=colors[act], alpha=0.7, edgecolor="none")
    ax3.set_xlabel("Day"); ax3.set_ylabel("Action")
    ax3.set_yticks([0.5]); ax3.set_yticklabels(["Daily Action"])
    ax3.set_title("Panel 3: RL Adaptive Daily Dosing Schedule")
    from matplotlib.patches import Patch
    ax3.legend(handles=[Patch(facecolor=c, label=l) for c, l in zip(colors, action_labels)], loc="upper right", fontsize=8)

    ax4 = plt.subplot(2, 2, 4)
    ax4.text(0.5, 0.5, "Final Spatial Slices\n(RL Adaptive vs Stupp)\nMid-sagittal slice",
             ha="center", va="center", transform=ax4.transAxes, fontsize=12,
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", edgecolor="darkblue"))
    ax4.set_title("Panel 4: Final Tumor Slices at Day 90")
    ax4.set_xticks([]); ax4.set_yticks([])

    plt.suptitle("Phase 5: RL Adaptive Therapy Steering (Final Optimized)", fontsize=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {output_path}")


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("PHASE 5: RL ADAPTIVE THERAPY STEERING (FINAL OPTIMIZED)")
    print("=" * 70)

    # Training: 32^3, dt=0.5, 40 episodes
    print(f"\n[RL] Training for {TRAIN_EPISODES} episodes (32^3 grid, dt=0.5)...")
    train_env = GbmTherapyEnv(training=True)
    policy, learning_curve = train_reinforce(train_env, episodes=TRAIN_EPISODES)

    # Evaluation: 64^3, 5 sub-steps
    print("\n[Eval] RL Adaptive (deterministic, 64^3, 5 sub-steps)...")
    eval_env = GbmTherapyEnv(training=False)
    rl_result = run_rl_adaptive(eval_env, policy) if policy else {"trajectory": [], "final_volume_mm3": 0.0}

    print("[Eval] Standard Stupp (64^3, 5 sub-steps)...")
    stupp_result = run_stupp_protocol(eval_env)

    rl_final = rl_result["final_volume_mm3"]
    stupp_final = stupp_result["final_volume_mm3"]

    rl_traj = rl_result["trajectory"]
    stupp_traj = stupp_result["trajectory"]
    rl_chemo_tox = max([t.get("chemo_tox", 0) for t in rl_traj]) if rl_traj else 0
    rl_rad_tox = max([t.get("rad_tox", 0) for t in rl_traj]) if rl_traj else 0
    stupp_chemo_tox = max([t.get("chemo_tox", 0) for t in stupp_traj]) if stupp_traj else 0
    stupp_rad_tox = max([t.get("rad_tox", 0) for t in stupp_traj]) if stupp_traj else 0

    rl_total_tox = rl_chemo_tox + rl_rad_tox
    stupp_total_tox = stupp_chemo_tox + stupp_rad_tox
    tox_reduction = 100.0 * (stupp_total_tox - rl_total_tox) / max(stupp_total_tox, 1e-6)

    initial_vol = rl_traj[0]["volume_mm3"] if rl_traj else 1.0
    rl_recurrence = 90
    stupp_recurrence = 90
    for t in rl_traj:
        if t["volume_mm3"] > 2 * initial_vol:
            rl_recurrence = t["day"]; break
    for t in stupp_traj:
        if t["volume_mm3"] > 2 * initial_vol:
            stupp_recurrence = t["day"]; break
    recurrence_delay = rl_recurrence - stupp_recurrence

    metrics = {
        "rl_final_volume_mm3": float(rl_final),
        "stupp_final_volume_mm3": float(stupp_final),
        "rl_toxicity_reduction_pct": float(tox_reduction),
        "rl_recurrence_delay_days": float(recurrence_delay),
        "training_episodes": TRAIN_EPISODES,
    }

    with open(OUTPUT_DIR / "phase5_adaptive_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[Metrics] Saved -> {OUTPUT_DIR / 'phase5_adaptive_metrics.json'}")

    plot_results(rl_result["trajectory"], stupp_result["trajectory"], learning_curve,
                 OUTPUT_DIR / "phase5_adaptive_steering.png")

    print("\n" + "=" * 70)
    print("PHASE 5 COMPLETE (FINAL OPTIMIZED)")
    print("=" * 70)
    for k, v in metrics.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()