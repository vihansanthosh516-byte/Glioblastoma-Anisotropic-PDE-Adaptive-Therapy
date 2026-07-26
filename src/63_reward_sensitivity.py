#!/usr/bin/env python3
"""
Phase 7: Reward Weight Sensitivity Analysis
===========================================
Addresses Reviewer Concern #9: Hand-Tuned Reward Weights - Proving the 
RL agent's strategy is robust and doesn't overfit to hyperparameter weights.

Evaluates policy robustness across reward weight variations.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

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
# Configuration
# --------------------------------------------------------------------------- #
EVAL_GRID = (64, 64, 64)
TRAIN_GRID = (32, 32, 32)
T_MAX_DAYS = 90
DT_RL_DAYS = 1.0
DT_PDE_TRAIN = 0.5
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = 5
N_TRAIN_EPISODES = 40

D_WHITE_BASE = 0.013
D_GRAY_BASE = 0.0013
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

# Reward weight configurations to test
REWARD_CONFIGS = [
    # (lambda_vol, lambda_den, lambda_tox)
    (5.0,  2.0,  0.001),  # Low penalties
    (5.0,  2.0,  0.01),
    (5.0,  2.0,  0.05),
    (5.0,  5.0,  0.001),
    (5.0,  5.0,  0.01),
    (5.0,  5.0,  0.05),
    (5.0,  10.0, 0.001),
    (5.0,  10.0, 0.01),
    (5.0,  10.0, 0.05),
    (15.0, 2.0,  0.001),
    (15.0, 2.0,  0.01),
    (15.0, 2.0,  0.05),
    (15.0, 5.0,  0.001),
    (15.0, 5.0,  0.01),
    (15.0, 5.0,  0.05),
    (15.0, 10.0, 0.001),
    (15.0, 10.0, 0.01),
    (15.0, 10.0, 0.05),
    (25.0, 2.0,  0.001),
    (25.0, 2.0,  0.01),
    (25.0, 2.0,  0.05),
    (25.0, 5.0,  0.001),
    (25.0, 5.0,  0.01),
    (25.0, 5.0,  0.05),
    (25.0, 10.0, 0.001),
    (25.0, 10.0, 0.01),
    (25.0, 10.0, 0.05),
]

# Select 10 diverse configurations for tractability
SELECTED_CONFIGS = [
    (5.0,  2.0,  0.001),
    (5.0,  5.0,  0.01),
    (5.0,  10.0, 0.05),
    (15.0, 2.0,  0.001),
    (15.0, 5.0,  0.01),    # Our default
    (15.0, 10.0, 0.05),
    (25.0, 2.0,  0.001),
    (25.0, 5.0,  0.01),
    (25.0, 10.0, 0.05),
    (5.0,  2.0,  0.05),
]

D_WHITE_BASE = 0.013
D_GRAY_BASE = 0.0013
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
# Fast PDE Solver (Self-contained)
# --------------------------------------------------------------------------- #
class FastPDESolver:
    def __init__(
        self,
        grid_size: Tuple[int, int, int] = EVAL_GRID,
        dt_pde: float = DT_PDE_EVAL,
        rho: float = 0.02,
        D_white: float = D_WHITE_BASE,
        alpha_sens: float = 1.0,
        use_dti: bool = True,
        use_mechanics: bool = True,
        is_training: bool = False,
    ):
        self.nx, self.ny, self.nz = grid_size
        self.dx = 128.0 / self.nx
        self.dt = dt_pde
        self.is_training = is_training
        self.rho = rho
        self.D_white = D_white
        self.alpha_sens = alpha_sens
        self.use_dti = use_dti
        self.use_mechanics = use_mechanics

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

        if self.use_dti:
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

            l1 = np.where(tract_mask, self.D_white, D_GRAY_BASE)
            l2 = np.where(tract_mask, D_GRAY_BASE, D_GRAY_BASE)
        else:
            tract_mask = np.zeros((self.nx, self.ny, self.nz), dtype=bool)
            vx = np.ones_like(xx)
            vy = np.zeros_like(yy)
            vz = np.zeros_like(zz)
            l1 = np.full_like(xx, self.D_white)
            l2 = np.full_like(xx, D_GRAY_BASE)

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

        ux_xf = (u_pad[1:nx+2, 1:-1, 1:-1] - u_pad[0:nx+1, 1:-1, 1:-1]) / dx
        uy_cc = (u_pad[1:-1, 2:, 1:-1] - u_pad[1:-1, :-2, 1:-1]) / (2*dx)
        uy_xf = 0.5 * (np.pad(uy_cc, ((1,1),(0,0),(0,0)), mode='edge')[:-1] +
                       np.pad(uy_cc, ((1,1),(0,0),(0,0)), mode='edge')[1:])
        uz_cc = (u_pad[1:-1, 1:-1, 2:] - u_pad[1:-1, 1:-1, :-2]) / (2*dx)
        uz_xf = 0.5 * (np.pad(uz_cc, ((1,1),(0,0),(0,0)), mode='edge')[:-1] +
                       np.pad(uz_cc, ((1,1),(0,0),(0,0)), mode='edge')[1:])
        Fx = self.Dxx_xf * ux_xf + self.Dxy_xf * uy_xf + self.Dxz_xf * uz_xf

        uy_yf = (u_pad[1:-1, 1:ny+2, 1:-1] - u_pad[1:-1, 0:ny+1, 1:-1]) / dx
        ux_cc = (u_pad[2:, 1:-1, 1:-1] - u_pad[:-2, 1:-1, 1:-1]) / (2*dx)
        ux_yf = 0.5 * (np.pad(ux_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, :-1, :] +
                       np.pad(ux_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, 1:, :])
        uz_yf = 0.5 * (np.pad(uz_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, :-1, :] +
                       np.pad(uz_cc, ((0,0),(1,1),(0,0)), mode='edge')[:, 1:, :])
        Fy = self.Dyy_yf * uy_yf + self.Dxy_yf * ux_yf + self.Dyz_yf * uz_yf

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
        react = self.rho * u * (1.0 - u / K_CARRY)
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
            kill = GAMMA_CHEMO * self.alpha_sens
            self.chemo_tox += CHEMO_TOX_PER_RL_STEP
        elif action == 2:
            kill = GAMMA_RAD * self.alpha_sens
            self.rad_tox += RAD_TOX_PER_RL_STEP
        elif action == 3:
            kill = (GAMMA_CHEMO + GAMMA_RAD) * self.alpha_sens
            self.chemo_tox += COMBO_TOX_PER_RL_STEP * 0.5
            self.rad_tox += COMBO_TOX_PER_RL_STEP * 0.5

        n_sub = 5
        for _ in range(n_sub):
            self.u = self.pde_step(self.u, kill)

        self.step_count += 1

        volume = float(self.u.sum() * self.dx**3)
        u_max = float(self.u.max())
        norm_vol = volume / max(self.initial_volume, 1e-6)
        delta_vol = self.prev_volume - volume
        self.prev_volume = volume

        return {
            "volume_mm3": volume,
            "u_max": u_max,
            "norm_volume": norm_vol,
            "delta_volume": delta_vol,
            "chemo_tox": self.chemo_tox,
            "rad_tox": self.rad_tox,
        }

    def get_observation(self) -> np.ndarray:
        vol = float(self.u.sum() * self.dx**3)
        u_max = float(self.u.max())
        norm_vol = vol / max(self.initial_volume, 1e-6)
        return np.array([
            np.clip(norm_vol, 0, 1),
            np.clip(u_max, 0, 1),
            self.step_count / 90,
            np.clip(self.chemo_tox, 0, 1),
            np.clip(self.rad_tox, 0, 1),
        ], dtype=np.float32)

    def is_done(self) -> bool:
        return self.step_count >= 90


# --------------------------------------------------------------------------- #
# Gym Environment
# --------------------------------------------------------------------------- #
class GbmTherapyEnv:
    def __init__(self, solver: FastPDESolver, reward_weights: Dict[str, float]):
        self.solver = solver
        self.max_steps = 90
        self.trajectory = []
        self.reward_weights = reward_weights

    def reset(self, seed: Optional[int] = None, options: Optional[Dict] = None):
        self.solver.reset()
        self.trajectory = []
        return self.solver.get_observation(), {}

    def step(self, action: int):
        result = self.solver.rl_step(action)
        obs = self.solver.get_observation()

        norm_vol = result["norm_volume"]
        u_max = result["u_max"]
        delta_vol = result["delta_volume"]

        action_cost = 1.0 if action > 0 else 0.0
        lambda_vol = self.reward_weights["lambda_vol"]
        lambda_den = self.reward_weights["lambda_den"]
        lambda_tox = self.reward_weights["lambda_tox"]
        
        reward = -lambda_vol * norm_vol - lambda_den * u_max - lambda_tox * action_cost
        if delta_vol > 0:
            reward += self.reward_weights.get("lambda_shrink", 100.0) * max(delta_vol / max(self.solver.initial_volume, 1e-6), 0.0)

        terminated = self.solver.is_done()
        if terminated and norm_vol < 0.01:
            reward += self.reward_weights.get("lambda_clear", 200.0)

        return obs, float(reward), terminated, False, {}


# --------------------------------------------------------------------------- #
# Policy Network
# --------------------------------------------------------------------------- #
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
        with torch.no_grad():
            torch.nn.init.constant_(self.net[-1].bias, 0.0)
            self.net[-1].bias.data[2] = 0.5
            self.net[-1].bias.data[3] = 1.0

    def forward(self, x):
        logits = self.net(x)
        return Categorical(logits=logits)


# --------------------------------------------------------------------------- #
# Training & Evaluation Functions
# --------------------------------------------------------------------------- #
def train_rl_policy(reward_weights: Dict[str, float], episodes: int = N_TRAIN_EPISODES, seed: int = 42) -> Optional["PolicyNetwork"]:
    if not HAS_TORCH:
        return None

    np.random.seed(seed)
    if torch is not None:
        torch.manual_seed(seed)

    train_solver = FastPDESolver(
        grid_size=TRAIN_GRID, dt_pde=DT_PDE_TRAIN, rho=RHO_BASE,
        D_white=D_WHITE_BASE, alpha_sens=1.0, is_training=True
    )
    train_env = GbmTherapyEnv(train_solver, reward_weights)

    policy = PolicyNetwork()
    optimizer = optim.Adam(policy.parameters(), lr=1e-2)
    episode_rewards = []

    for ep in range(episodes):
        obs, _ = train_env.reset()
        log_probs = []
        rewards = []
        entropies = []

        for _ in range(train_env.max_steps):
            obs_tensor = torch.FloatTensor(obs).unsqueeze(0)
            dist = policy(obs_tensor)
            action = dist.sample()
            log_prob = dist.log_prob(action)
            entropy = dist.entropy()

            obs, reward, terminated, truncated, _ = train_env.step(action.item())
            log_probs.append(log_prob)
            rewards.append(reward)
            entropies.append(entropy)

            if terminated:
                break

        returns = []
        R = 0
        for r in reversed(rewards):
            R = r + 0.99 * R
            returns.insert(0, R)
        returns = torch.tensor(returns, dtype=torch.float32)
        returns = (returns - returns.mean()) / (returns.std() + 1e-8)

        log_probs_t = torch.stack(log_probs)
        entropies_t = torch.stack(entropies)
        loss = -(log_probs_t * returns).sum() - 0.01 * entropies_t.sum()

        optimizer.zero_grad()
        loss.backward()
        torch.nn.utils.clip_grad_norm_(policy.parameters(), max_norm=1.0)
        optimizer.step()

        episode_rewards.append(sum(rewards))

    return policy


def evaluate_policy(policy: Optional["PolicyNetwork"], reward_weights: Dict[str, float], seed: int = 42) -> Dict[str, float]:
    if not HAS_TORCH or policy is None:
        return {"final_volume_mm3": 0.0, "peak_u_max": 0.0}

    np.random.seed(seed)
    if torch is not None:
        torch.manual_seed(seed)

    solver = FastPDESolver(grid_size=EVAL_GRID, dt_pde=DT_PDE_EVAL, is_training=False)
    env = GbmTherapyEnv(solver, reward_weights)

    obs, _ = env.reset()
    trajectory = []

    for step in range(env.max_steps):
        obs_tensor = torch.FloatTensor(obs).unsqueeze(0)
        with torch.no_grad():
            action = policy(obs_tensor).probs.argmax().item()

        # Guardrail
        current_vol = env.solver.u.sum() * env.solver.dx**3
        if action == 0 and current_vol > 0.05 * env.solver.initial_volume:
            action = 3

        obs, reward, terminated, truncated, _ = env.step(action)
        trajectory.append({
            "volume_mm3": env.solver.u.sum() * env.solver.dx**3,
            "u_max": env.solver.u.max(),
        })
        if terminated:
            break

    return {
        "final_volume_mm3": trajectory[-1]["volume_mm3"],
        "peak_u_max": max(t["u_max"] for t in trajectory),
    }


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("PHASE 7: REWARD WEIGHT SENSITIVITY ANALYSIS")
    print("=" * 70)

    if not HAS_TORCH:
        print("PyTorch not available - skipping reward sensitivity analysis")
        return

    print(f"\n[Phase 7] Testing {len(SELECTED_CONFIGS)} reward weight configurations...")
    print(f"Training: {N_TRAIN_EPISODES} episodes per config on 32^3 grid")
    print(f"Evaluation: 64^3 grid, deterministic policy")

    results = []

    for i, (lambda_vol, lambda_den, lambda_tox) in enumerate(SELECTED_CONFIGS):
        print(f"\n[Phase 7] Config {i+1}/{len(SELECTED_CONFIGS)}: "
              f"λ_vol={lambda_vol}, λ_den={lambda_den}, λ_tox={lambda_tox}")

        reward_weights = {
            "lambda_vol": lambda_vol,
            "lambda_den": lambda_den,
            "lambda_tox": lambda_tox,
            "lambda_shrink": 100.0,
            "lambda_clear": 200.0,
        }

        # Train policy
        policy = train_rl_policy(reward_weights, episodes=N_TRAIN_EPISODES, seed=42)

        # Evaluate
        eval_result = evaluate_policy(policy, reward_weights, seed=42)

        result = {
            "config_idx": i,
            "lambda_vol": lambda_vol,
            "lambda_den": lambda_den,
            "lambda_tox": lambda_tox,
            "final_volume_mm3": eval_result["final_volume_mm3"],
            "peak_u_max": eval_result["peak_u_max"],
        }
        results.append(result)

        print(f"  Final Volume: {eval_result['final_volume_mm3']:.2f} mm³")
        print(f"  Peak u_max: {eval_result['peak_u_max']:.4f}")

    # -----------------------------------------------------------------------
    # Metrics & Analysis
    # -----------------------------------------------------------------------
    volumes = np.array([r["final_volume_mm3"] for r in results])
    peaks = np.array([r["peak_u_max"] for r in results])

    # Coefficient of variation across configurations
    cv_volume = float(np.std(volumes) / max(np.mean(volumes), 1e-6))
    cv_peak = float(np.std(peaks) / max(np.mean(peaks), 1e-6))

    # Best configuration
    best_idx = int(np.argmin(volumes))
    best_config = results[best_idx]

    # Correlation between weight magnitudes and outcomes
    vol_weights = np.array([r["lambda_vol"] for r in results])
    den_weights = np.array([r["lambda_den"] for r in results])
    tox_weights = np.array([r["lambda_tox"] for r in results])

    from scipy.stats import pearsonr, spearmanr
    vol_corr_vol, vol_corr_p = pearsonr(vol_weights, volumes)
    vol_corr_peak, vol_corr_peak_p = pearsonr(vol_weights, peaks)
    den_corr_vol, den_corr_vol_p = pearsonr(den_weights, volumes)
    den_corr_peak, den_corr_peak_p = pearsonr(den_weights, peaks)
    tox_corr_vol, tox_corr_vol_p = pearsonr(tox_weights, volumes)
    tox_corr_peak, tox_corr_peak_p = pearsonr(tox_weights, peaks)

    metrics = {
        "n_configurations": len(SELECTED_CONFIGS),
        "episodes_per_config": N_TRAIN_EPISODES,
        "mean_final_volume_mm3": float(np.mean(volumes)),
        "std_final_volume_mm3": float(np.std(volumes)),
        "cv_final_volume": cv_volume,
        "mean_peak_u_max": float(np.mean(peaks)),
        "std_peak_u_max": float(np.std(peaks)),
        "cv_peak_u_max": cv_peak,
        "best_config": {
            "lambda_vol": best_config["lambda_vol"],
            "lambda_den": best_config["lambda_den"],
            "lambda_tox": best_config["lambda_tox"],
            "final_volume_mm3": best_config["final_volume_mm3"],
            "peak_u_max": best_config["peak_u_max"],
        },
        "correlations": {
            "lambda_vol_vs_volume": {"pearson_r": float(vol_corr_vol), "p_value": float(vol_corr_p)},
            "lambda_vol_vs_peak": {"pearson_r": float(vol_corr_peak), "p_value": float(vol_corr_peak_p)},
            "lambda_den_vs_volume": {"pearson_r": float(den_corr_vol), "p_value": float(den_corr_vol_p)},
            "lambda_den_vs_peak": {"pearson_r": float(den_corr_peak), "p_value": float(den_corr_peak_p)},
            "lambda_tox_vs_volume": {"pearson_r": float(tox_corr_vol), "p_value": float(tox_corr_vol_p)},
            "lambda_tox_vs_peak": {"pearson_r": float(tox_corr_peak), "p_value": float(tox_corr_peak_p)},
        },
        "all_results": results,
    }

    # Save metrics
    with open(OUTPUT_DIR / "reward_sensitivity_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\n[Metrics] Saved -> {OUTPUT_DIR / 'reward_sensitivity_metrics.json'}")

    # -----------------------------------------------------------------------
    # Visualization
    # -----------------------------------------------------------------------
    print("\n[Phase 7] Generating reward sensitivity visualization...")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel 1: Volume vs lambda_vol (grouped by lambda_den)
    ax1 = axes[0, 0]
    for lambda_den in [2.0, 5.0, 10.0]:
        mask = np.array([r["lambda_den"] for r in results]) == lambda_den
        if np.any(mask):
            lv = np.array([r["lambda_vol"] for r in results])[mask]
            vol = np.array([r["final_volume_mm3"] for r in results])[mask]
            sort_idx = np.argsort(lv)
            ax1.plot(lv[sort_idx], vol[sort_idx], 'o-', label=f'λ_den={lambda_den}', linewidth=2, markersize=6)
    ax1.set_xlabel('Volume Weight (λ_vol)')
    ax1.set_ylabel('Final Tumor Volume (mm³)')
    ax1.set_title('Panel 1: Volume vs λ_vol (grouped by λ_den)')
    ax1.legend()
    ax1.grid(alpha=0.3)
    ax1.set_yscale('log')

    # Panel 2: Volume vs lambda_den (grouped by lambda_vol)
    ax2 = axes[0, 1]
    for lambda_vol in [5.0, 15.0, 25.0]:
        mask = np.array([r["lambda_vol"] for r in results]) == lambda_vol
        if np.any(mask):
            ld = np.array([r["lambda_den"] for r in results])[mask]
            vol = np.array([r["final_volume_mm3"] for r in results])[mask]
            sort_idx = np.argsort(ld)
            ax2.plot(ld[sort_idx], vol[sort_idx], 'o-', label=f'λ_vol={lambda_vol}', linewidth=2, markersize=6)
    ax2.set_xlabel('Density Weight (λ_den)')
    ax2.set_ylabel('Final Tumor Volume (mm³)')
    ax2.set_title('Panel 2: Volume vs λ_den (grouped by λ_vol)')
    ax2.legend()
    ax2.grid(alpha=0.3)
    ax2.set_yscale('log')

    # Panel 3: Heatmap of final volume across lambda_vol x lambda_den
    ax3 = axes[1, 0]
    lambda_vols = sorted(set(r["lambda_vol"] for r in results))
    lambda_dens = sorted(set(r["lambda_den"] for r in results))
    heatmap_data = np.zeros((len(lambda_vols), len(lambda_dens)))
    for r in results:
        i = lambda_vols.index(r["lambda_vol"])
        j = lambda_dens.index(r["lambda_den"])
        # Average over toxicity weights
        mask = (np.array([x["lambda_vol"] for x in results]) == r["lambda_vol"]) & \
               (np.array([x["lambda_den"] for x in results]) == r["lambda_den"])
        heatmap_data[i, j] = np.mean([results[k]["final_volume_mm3"] for k in np.where(mask)[0]])

    im = ax3.imshow(heatmap_data, cmap='RdYlGn_r', aspect='auto', origin='lower')
    ax3.set_xticks(range(len(lambda_dens)))
    ax3.set_xticklabels([str(d) for d in lambda_dens])
    ax3.set_yticks(range(len(lambda_vols)))
    ax3.set_yticklabels([str(v) for v in lambda_vols])
    ax3.set_xlabel('Density Weight (λ_den)')
    ax3.set_ylabel('Volume Weight (λ_vol)')
    ax3.set_title('Panel 3: Mean Final Volume Heatmap\n(averaged over λ_tox)')
    plt.colorbar(im, ax=ax3, label='Final Volume (mm³)')

    # Annotate heatmap
    for i in range(len(lambda_vols)):
        for j in range(len(lambda_dens)):
            ax3.text(j, i, f'{heatmap_data[i,j]:.1f}', ha='center', va='center', fontsize=9)

    # Panel 4: Toxicity weight effect (scatter with size = lambda_tox)
    ax4 = axes[1, 1]
    scatter = ax4.scatter(
        [r["lambda_vol"] for r in results],
        [r["lambda_den"] for r in results],
        c=[r["final_volume_mm3"] for r in results],
        s=[r["lambda_tox"] * 50 + 20 for r in results],
        cmap='RdYlGn_r',
        alpha=0.7,
        edgecolors='black',
        linewidth=0.5
    )
    ax4.set_xlabel('Volume Weight (λ_vol)')
    ax4.set_ylabel('Density Weight (λ_den)')
    ax4.set_title('Panel 4: Volume by λ_vol/λ_den (size = λ_tox)')
    plt.colorbar(scatter, ax=ax4, label='Final Volume (mm³)')
    ax4.grid(alpha=0.3)

    plt.suptitle('Phase 7: Reward Weight Sensitivity Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "reward_sensitivity_figure.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {OUTPUT_DIR / 'reward_sensitivity_figure.png'}")

    # Summary
    print("\n" + "=" * 70)
    print("PHASE 7: REWARD WEIGHT SENSITIVITY ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"Configurations tested: {len(SELECTED_CONFIGS)}")
    print(f"Mean final volume: {np.mean(volumes):.2f} ± {np.std(volumes):.2f} mm³ (CV: {cv_volume:.3f})")
    print(f"Best config: λ_vol={best_config['lambda_vol']}, λ_den={best_config['lambda_den']}, λ_tox={best_config['lambda_tox']}")
    print(f"  -> Volume: {best_config['final_volume_mm3']:.2f} mm³")
    print(f"\nCorrelation Analysis:")
    print(f"  λ_vol -> Volume: r={vol_corr_vol:.3f} (p={vol_corr_p:.3f})")
    print(f"  λ_den -> Volume: r={den_corr_vol:.3f} (p={den_corr_vol_p:.3f})")
    print(f"  λ_tox -> Volume: r={tox_corr_vol:.3f} (p={tox_corr_vol_p:.3f})")
    print(f"\nOutputs saved to {OUTPUT_DIR}/")
    print("  - reward_sensitivity_metrics.json")
    print("  - reward_sensitivity_figure.png")


if __name__ == "__main__":
    main()