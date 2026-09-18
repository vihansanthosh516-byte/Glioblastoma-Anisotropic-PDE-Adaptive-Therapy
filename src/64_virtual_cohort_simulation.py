#!/usr/bin/env python3
"""
Phase 8: Multi-Patient Virtual Cohort Validation
==================================================
Validates RL Adaptive Therapy against Standard Stupp Protocol across
a virtual cohort of 20 patients with clinically sampled parameters.

Outputs:
  - output/phase8_cohort_analysis.png (4-panel figure)
  - output/phase8_cohort_metrics.json (statistical metrics)
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
from scipy import stats
from scipy.stats import norm

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
N_PATIENTS = 20
EVAL_GRID = (64, 64, 64)
T_MAX_DAYS = 90
DT_RL_DAYS = 1.0
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = 5

# Clinical parameter distributions
RHO_MEAN, RHO_STD = 0.025, 0.005        # 1/day
D_WHITE_MEAN, D_WHITE_STD = 0.0012, 0.0003  # cm^2/day
GAMMA_CHEMO_MIN, GAMMA_CHEMO_MAX = 0.02, 0.08
ALPHA_RT_MIN, ALPHA_RT_MAX = 0.015, 0.045

# Fixed parameters
D_GRAY_BASE = 0.0013
RHO_BASE = 0.02
K_CARRY = 1.0
GAMMA_RAD_BASE = 0.08
CHEMO_TOX_PER_RL_STEP = 0.02
RAD_TOX_PER_RL_STEP = 0.05
COMBO_TOX_PER_RL_STEP = 0.08
SEED_SIGMA_MM = 5.0
SEED_AMPLITUDE = 0.8
DT_RL_DAYS = 1.0
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = 5

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# RL Reward Weights (tuned from Phase 7)
RL_REWARD_WEIGHTS = {
    "lambda_vol": 15.0,
    "lambda_den": 5.0,
    "lambda_tox": 0.01,
    "lambda_shrink": 100.0,
    "lambda_clear": 200.0,
}

# --------------------------------------------------------------------------- #
# Fast PDE Solver (Self-contained)
# --------------------------------------------------------------------------- #
class FastPDESolver:
    def __init__(
        self,
        grid_size: Tuple[int, int, int] = EVAL_GRID,
        dt_pde: float = DT_PDE_EVAL,
        rho: float = RHO_BASE,
        D_white: float = D_WHITE_MEAN,
        alpha_sens: float = 1.0,
        gamma_chemo: float = GAMMA_CHEMO_MAX,
        alpha_rt: float = ALPHA_RT_MAX,
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
        self.gamma_chemo = gamma_chemo
        self.alpha_rt = alpha_rt
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
            kill = self.gamma_chemo * self.alpha_sens
            self.chemo_tox += CHEMO_TOX_PER_RL_STEP
        elif action == 2:
            kill = self.alpha_rt * self.alpha_sens
            self.rad_tox += RAD_TOX_PER_RL_STEP
        elif action == 3:
            kill = (self.gamma_chemo + self.alpha_rt) * self.alpha_sens
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
            "u_max": float(self.u.max()),
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
# Policy Network (Heuristic fallback - same as Phase 5/6/7)
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
# Baseline Protocols
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
            "u_max": float(env.solver.u.max()),
            "reward": reward,
        })
        if terminated:
            break
    return {
        "trajectory": trajectory,
        "final_volume_mm3": trajectory[-1]["volume_mm3"],
        "peak_u_max": max(t.get("u_max", 0) for t in trajectory) if trajectory else 0,
        "progressed": trajectory[-1]["volume_mm3"] > 500,
    }


def run_rl_adaptive(env: GbmTherapyEnv, policy=None) -> Dict:
    obs, _ = env.reset()
    trajectory = []

    initial_vol = env.solver.initial_volume

    for step in range(env.max_steps):
        if policy is not None and HAS_TORCH:
            obs_tensor = torch.FloatTensor(obs).unsqueeze(0)
            with torch.no_grad():
                action = policy(obs_tensor).probs.argmax().item()
        else:
            # Heuristic RL policy (same as Phase 5/6/7)
            current_vol = env.solver.u.sum() * env.solver.dx**3
            norm_vol = current_vol / max(env.solver.initial_volume, 1e-6)

            if norm_vol > 0.05:
                action = 3
            elif norm_vol > 0.01:
                action = 2
            else:
                action = 0

        # Guardrail
        current_vol = env.solver.u.sum() * env.solver.dx**3
        if action == 0 and current_vol > 0.05 * env.solver.initial_volume:
            action = 3

        obs, reward, terminated, truncated, _ = env.step(action)
        u_max = float(env.solver.u.max())
        trajectory.append({
            "step": step, "day": step + 1, "action": action,
            "volume_mm3": env.solver.u.sum() * env.solver.dx**3,
            "u_max": u_max,
            "reward": reward,
        })
        if terminated:
            break
    return {
        "trajectory": trajectory,
        "final_volume_mm3": trajectory[-1]["volume_mm3"],
        "peak_u_max": max(t.get("u_max", 0) for t in trajectory) if trajectory else 0,
        "progressed": trajectory[-1]["volume_mm3"] > 500,
    }


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def generate_virtual_cohort(n_patients: int = N_PATIENTS, seed: int = 12345) -> List[Dict[str, float]]:
    """Generate virtual patient cohort with clinically sampled parameters."""
    np.random.seed(seed)
    
    cohort = []
    for i in range(n_patients):
        # Sample from clinical distributions
        rho = max(0.005, min(0.05, np.random.normal(RHO_MEAN, RHO_STD)))
        D_white = max(0.0005, min(0.003, np.random.normal(D_WHITE_MEAN, D_WHITE_STD)))
        gamma_chemo = np.random.uniform(GAMMA_CHEMO_MIN, GAMMA_CHEMO_MAX)
        alpha_rt = np.random.uniform(ALPHA_RT_MIN, ALPHA_RT_MAX)
        
        cohort.append({
            "patient_id": i,
            "rho": float(rho),
            "D_white": float(D_white),
            "gamma_chemo": float(gamma_chemo),
            "alpha_rt": float(alpha_rt),
        })
    return cohort


def evaluate_patient(patient_params: Dict[str, float], rl_policy=None, patient_id: int = 0) -> Dict[str, Any]:
    """Evaluate a single patient with both Stupp and RL Adaptive protocols."""
    solver = FastPDESolver(
        grid_size=EVAL_GRID,
        dt_pde=DT_PDE_EVAL,
        rho=patient_params["rho"],
        D_white=patient_params["D_white"],
        alpha_sens=1.0,
        gamma_chemo=patient_params["gamma_chemo"],
        alpha_rt=patient_params["alpha_rt"],
        is_training=False,
    )
    env = GbmTherapyEnv(solver, RL_REWARD_WEIGHTS)

    # Evaluate Stupp
    stupp_result = run_stupp_protocol(env)
    
    # Recreate environment for RL evaluation (fresh solver state)
    solver_rl = FastPDESolver(
        grid_size=EVAL_GRID,
        dt_pde=DT_PDE_EVAL,
        rho=patient_params["rho"],
        D_white=patient_params["D_white"],
        alpha_sens=1.0,
        gamma_chemo=patient_params["gamma_chemo"],
        alpha_rt=patient_params["alpha_rt"],
        is_training=False,
    )
    env_rl = GbmTherapyEnv(solver_rl, RL_REWARD_WEIGHTS)
    rl_result = run_rl_adaptive(env_rl, rl_policy)

    return {
        "patient_id": patient_id,
        "params": {
            "rho": patient_params["rho"],
            "D_white": patient_params["D_white"],
            "gamma_chemo": patient_params["gamma_chemo"],
            "alpha_rt": patient_params["alpha_rt"],
        },
        "stupp": {
            "final_volume_mm3": stupp_result["final_volume_mm3"],
            "peak_u_max": stupp_result["peak_u_max"],
            "progressed": stupp_result["progressed"],
            "trajectory": stupp_result["trajectory"],
        },
        "rl_adaptive": {
            "final_volume_mm3": rl_result["final_volume_mm3"],
            "peak_u_max": rl_result["peak_u_max"],
            "progressed": rl_result["progressed"],
            "trajectory": rl_result["trajectory"],
        },
    }


# --------------------------------------------------------------------------- #
# Statistical Analysis
# --------------------------------------------------------------------------- #
def compute_statistics(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Compute statistical metrics from cohort results."""
    stupp_volumes = np.array([r["stupp"]["final_volume_mm3"] for r in results])
    rl_volumes = np.array([r["rl_adaptive"]["final_volume_mm3"] for r in results])
    stupp_progressed = np.array([r["stupp"]["progressed"] for r in results])
    rl_progressed = np.array([r["rl_adaptive"]["progressed"] for r in results])
    
    # Paired t-test
    t_stat, p_value = stats.ttest_rel(rl_volumes, stupp_volumes)
    
    # Wilcoxon signed-rank (non-parametric)
    wilcoxon_stat, wilcoxon_p = stats.wilcoxon(rl_volumes, stupp_volumes)
    
    # Progression-free rates
    stupp_pf_rate = float(np.mean(~stupp_progressed))
    rl_pf_rate = float(np.mean(~rl_progressed))
    
    # McNemar's test for progression rates
    both_progressed = np.sum(stupp_progressed & rl_progressed)
    stupp_only = np.sum(stupp_progressed & ~rl_progressed)
    rl_only = np.sum(rl_progressed & ~stupp_progressed)
    neither = np.sum(~stupp_progressed & ~rl_progressed)
    
    # Effect size (Cohen's d)
    diff = rl_volumes - stupp_volumes
    cohens_d = float(np.mean(diff) / np.std(diff)) if np.std(diff) > 0 else 0.0
    
    # Cohen's kappa for progression agreement
    # Agreement matrix:
    # both_progressed, stupp_only, rl_only, neither
    both_prog = both_progressed
    stupp_only_prog = stupp_only
    rl_only_prog = rl_only
    neither_prog = neither
    total = len(results)
    
    po = (both_prog + neither_prog) / total  # observed agreement
    pe = ((both_prog + stupp_only_prog) * (both_prog + rl_only_prog) + 
          (stupp_only_prog + neither_prog) * (rl_only_prog + neither_prog)) / (total * total)  # expected agreement
    cohens_kappa = float((po - pe) / max(1 - pe, 1e-10)) if pe < 1 else 1.0

    return {
        "cohort_size": len(results),
        "rl_mean_final_volume_mm3": float(np.mean(rl_volumes)),
        "rl_std_final_volume_mm3": float(np.std(rl_volumes)),
        "stupp_mean_final_volume_mm3": float(np.mean(stupp_volumes)),
        "stupp_std_final_volume_mm3": float(np.std(stupp_volumes)),
        "paired_t_test_p_value": float(p_value),
        "wilcoxon_p_value": float(wilcoxon_p),
        "cohens_d": cohens_d,
        "rl_progression_free_rate": rl_pf_rate,
        "stupp_progression_free_rate": stupp_pf_rate,
        "progression_difference": float(rl_pf_rate - stupp_pf_rate),
        "mcnemar_chi2": float((abs(stupp_only - rl_only) - 1)**2 / max(stupp_only + rl_only, 1)) if (stupp_only + rl_only) > 0 else 0.0,
        "mcnemar_p_value": float(stats.chi2.sf((abs(stupp_only - rl_only) - 1)**2 / max(stupp_only + rl_only, 1), 1)) if (stupp_only + rl_only) > 0 else 1.0,
        "cohens_kappa": cohens_kappa,
    }


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def create_visualization(results: List[Dict[str, Any]], stats: Dict[str, Any], output_path: Path):
    """Create 4-panel cohort analysis figure."""
    fig = plt.figure(figsize=(16, 12))

    # Extract trajectory data
    days = np.arange(1, 91)
    
    # Collect trajectories
    stupp_trajs = []
    rl_trajs = []
    for r in results:
        stupp_traj = [t["volume_mm3"] for t in r["stupp"]["trajectory"]]
        rl_traj = [t["volume_mm3"] for t in r["rl_adaptive"]["trajectory"]]
        if len(stupp_traj) == 90:
            stupp_trajs.append(stupp_traj)
        if len(rl_traj) == 90:
            rl_trajs.append(rl_traj)
    
    stupp_trajs = np.array(stupp_trajs)  # (N, 90)
    rl_trajs = np.array(rl_trajs)
    
    # Panel 1: Population Volume Trajectories (Mean ± 95% CI)
    ax1 = plt.subplot(2, 2, 1)
    if len(stupp_trajs) > 0 and len(rl_trajs) > 0:
        days = np.arange(1, 91)
        
        stupp_mean = np.mean(stupp_trajs, axis=0)
        stupp_std = np.std(stupp_trajs, axis=0)
        stupp_ci_lower = stupp_mean - 1.96 * stupp_std / np.sqrt(len(stupp_trajs))
        stupp_ci_upper = stupp_mean + 1.96 * stupp_std / np.sqrt(len(stupp_trajs))
        
        rl_mean = np.mean(rl_trajs, axis=0)
        rl_std = np.std(rl_trajs, axis=0)
        rl_ci_lower = rl_mean - 1.96 * rl_std / np.sqrt(len(rl_trajs))
        rl_ci_upper = rl_mean + 1.96 * rl_std / np.sqrt(len(rl_trajs))
        
        ax1.plot(days, stupp_mean, 'r-', linewidth=2, label='Stupp (Mean)')
        ax1.fill_between(days, stupp_ci_lower, stupp_ci_upper, color='red', alpha=0.2, label='Stupp 95% CI')
        
        ax1.plot(days, rl_mean, 'b-', linewidth=2, label='RL Adaptive (Mean)')
        ax1.fill_between(days, rl_ci_lower, rl_ci_upper, color='blue', alpha=0.2, label='RL 95% CI')
        
        ax1.set_xlabel('Day')
        ax1.set_ylabel('Tumor Volume (mm³)')
        ax1.set_title('Panel 1: Population Volume Trajectories\n(Mean ± 95% CI)')
        ax1.set_yscale('log')
        ax1.set_ylim(0.5, 5000)
        ax1.legend(loc='upper right')
        ax1.grid(alpha=0.3)
    
    # Panel 2: Patient-Level Final Volume Comparison (Paired Scatter/Box Plot)
    ax2 = plt.subplot(2, 2, 2)
    stupp_finals = np.array([r["stupp"]["final_volume_mm3"] for r in results])
    rl_finals = np.array([r["rl_adaptive"]["final_volume_mm3"] for r in results])
    
    # Scatter with identity line
    max_vol = max(np.max(stupp_finals), np.max(rl_finals)) * 1.1
    ax2.plot([0, max_vol], [0, max_vol], 'k--', alpha=0.5, label='Identity')
    ax2.scatter(stupp_finals, rl_finals, c='purple', s=80, alpha=0.7, edgecolors='black', linewidth=0.5)
    
    # Highlight patients who progressed on Stupp but not RL
    for i, r in enumerate(results):
        if r["stupp"]["progressed"] and not r["rl_adaptive"]["progressed"]:
            ax2.scatter(r["stupp"]["final_volume_mm3"], r["rl_adaptive"]["final_volume_mm3"], 
                       c='green', s=120, marker='*', edgecolors='black', linewidth=1.5, zorder=10)
        elif r["rl_adaptive"]["progressed"] and not r["stupp"]["progressed"]:
            ax2.scatter(r["stupp"]["final_volume_mm3"], r["rl_adaptive"]["final_volume_mm3"], 
                       c='orange', s=120, marker='*', edgecolors='black', linewidth=1.5, zorder=10)
    
    ax2.set_xlabel('Stupp Final Volume (mm³)')
    ax2.set_ylabel('RL Adaptive Final Volume (mm³)')
    ax2.set_title('Panel 2: Patient-Level Final Volume\n(RL vs Stupp)')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend(loc='upper left')
    ax2.grid(alpha=0.3)
    
    # Add box plot inset
    ax2_inset = ax2.inset_axes([0.6, 0.1, 0.35, 0.35])
    box_data = [stupp_finals, rl_finals]
    bp = ax2_inset.boxplot(box_data, labels=['Stupp', 'RL'], patch_artist=True,
                           showfliers=True, widths=0.5)
    colors = ['red', 'blue']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
    ax2_inset.set_yscale('log')
    ax2_inset.set_ylabel('Volume (mm³)')
    ax2_inset.set_title('Distribution')
    
    # Panel 3: Toxicity vs Efficacy Trade-off
    ax3 = plt.subplot(2, 2, 3)
    for r in results:
        stupp_result = r["stupp"]
        rl_result = r["rl_adaptive"]
        # Compute cumulative toxicity
        stupp_traj = stupp_result["trajectory"]
        rl_traj = rl_result["trajectory"]
        
        stupp_tox = sum(1 for t in stupp_traj if t["action"] > 0) / 90.0
        rl_tox = sum(1 for t in rl_traj if t["action"] > 0) / 90.0
        
        stupp_vol = stupp_result["final_volume_mm3"]
        rl_vol = rl_result["final_volume_mm3"]
        
        ax3.scatter(stupp_tox, stupp_vol, c='red', s=60, alpha=0.6, edgecolors='black', marker='o')
        ax3.scatter(rl_tox, rl_vol, c='blue', s=60, alpha=0.6, edgecolors='black', marker='s')
        
        # Connect paired points
        ax3.plot([stupp_tox, rl_tox], [stupp_vol, rl_vol], 'k-', alpha=0.3, linewidth=0.5)
    
    ax3.set_xlabel('Cumulative Toxicity Exposure (fraction of days with therapy)')
    ax3.set_ylabel('Final Tumor Volume (mm³)')
    ax3.set_title('Panel 3: Toxicity vs Efficacy Trade-off\n(Red=Stupp, Blue=RL Adaptive)')
    ax3.set_yscale('log')
    ax3.legend(['Stupp', 'RL Adaptive'], loc='upper right')
    ax3.grid(alpha=0.3)
    
    # Panel 4: Kaplan-Meier Progression-Free Survival
    ax4 = plt.subplot(2, 2, 4)
    
    # Progression = volume > 500 mm³
    # Compute time-to-progression for each patient
    def get_ttp(trajectory):
        for t in trajectory:
            if t["volume_mm3"] > 500:
                return t["day"]
        return 90  # censored at 90 days
    
    stupp_ttp = np.array([get_ttp(r["stupp"]["trajectory"]) for r in results])
    rl_ttp = np.array([get_ttp(r["rl_adaptive"]["trajectory"]) for r in results])
    
    # Kaplan-Meier curves
    for label, ttp_data, color in [("Stupp", stupp_ttp, 'red'), ("RL Adaptive", rl_ttp, 'blue')]:
        # Sort by time
        sorted_idx = np.argsort(ttp_data)
        sorted_ttp = ttp_data[sorted_idx]
        
        # Kaplan-Meier estimator
        n_at_risk = len(ttp_data)
        survival = np.ones(n_at_risk)
        for i, t in enumerate(sorted_ttp):
            if t < 90:  # progression event
                survival[i:] = survival[i:] * (n_at_risk - i - 1) / max(n_at_risk - i, 1)
        
        ax4.step(sorted_ttp, survival, where='post', color=color, linewidth=2, label=label)
    
    ax4.set_xlabel('Day')
    ax4.set_ylabel('Progression-Free Probability')
    ax4.set_title('Panel 4: Kaplan-Meier Progression-Free Survival\n(Progression = Volume > 500 mm³)')
    ax4.set_xlim(0, 90)
    ax4.set_ylim(0, 1.05)
    ax4.legend(loc='lower left')
    ax4.grid(alpha=0.3)
    
    plt.suptitle('Phase 8: Multi-Patient Virtual Cohort Validation\n(20 Patients, 90-Day Horizon)', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "phase8_cohort_analysis.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {output_path}")


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("PHASE 8: MULTI-PATIENT VIRTUAL COHORT VALIDATION")
    print("=" * 70)

    # 1. Generate virtual cohort
    print(f"\n[Phase 8] Generating virtual cohort of {N_PATIENTS} patients...")
    cohort = generate_virtual_cohort(N_PATIENTS, seed=42)
    for p in cohort:
        print(f"  Patient {p['patient_id']}: rho={p['rho']:.4f}, D_white={p['D_white']:.4f}, "
              f"gamma_chemo={p['gamma_chemo']:.3f}, alpha_rt={p['alpha_rt']:.3f}")

    # 2. Evaluate each patient
    print(f"\n[Phase 8] Evaluating cohort (Stupp vs RL Adaptive)...")
    results = []
    for patient in cohort:
        print(f"  Patient {patient['patient_id']:2d}: rho={patient['rho']:.4f}, "
              f"D_white={patient['D_white']:.4f}...", end=" ", flush=True)
        result = evaluate_patient(patient, rl_policy=None, patient_id=patient["patient_id"])
        results.append(result)
        winner = "RL" if result["rl_adaptive"]["final_volume_mm3"] < result["stupp"]["final_volume_mm3"] else "Stupp"
        print(f"RL: {result['rl_adaptive']['final_volume_mm3']:.2f} mm³, "
              f"Stupp: {result['stupp']['final_volume_mm3']:.2f} mm³ -> {winner} wins")

    # 3. Statistics
    print("\n[Phase 8] Computing statistics...")
    stats = compute_statistics(results)

    # 4. Save metrics
    metrics = {
        "cohort_size": N_PATIENTS,
        "rl_mean_final_volume_mm3": float(stats["rl_mean_final_volume_mm3"]),
        "stupp_mean_final_volume_mm3": float(stats["stupp_mean_final_volume_mm3"]),
        "rl_std_final_volume_mm3": float(stats["rl_std_final_volume_mm3"]),
        "stupp_std_final_volume_mm3": float(stats["stupp_std_final_volume_mm3"]),
        "paired_t_test_p_value": float(stats["paired_t_test_p_value"]),
        "wilcoxon_p_value": float(stats["wilcoxon_p_value"]),
        "cohens_d": float(stats["cohens_d"]),
        "rl_progression_free_rate": float(stats["rl_progression_free_rate"]),
        "stupp_progression_free_rate": float(stats["stupp_progression_free_rate"]),
        "progression_difference": float(stats["progression_difference"]),
        "mcnemar_p_value": float(stats["mcnemar_p_value"]),
        "cohens_kappa": float(stats["cohens_kappa"]),
        "patient_details": [
            {
                "patient_id": int(r["patient_id"]),
                "params": {k: float(v) for k, v in r["params"].items()},
                "stupp_final_volume_mm3": float(r["stupp"]["final_volume_mm3"]),
                "rl_final_volume_mm3": float(r["rl_adaptive"]["final_volume_mm3"]),
                "stupp_progressed": bool(r["stupp"]["progressed"]),
                "rl_progressed": bool(r["rl_adaptive"]["progressed"]),
            }
            for r in results
        ],
    }

    with open(OUTPUT_DIR / "phase8_cohort_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"[Metrics] Saved -> {OUTPUT_DIR / 'phase8_cohort_metrics.json'}")

    # 5. Visualization
    print("\n[Phase 8] Generating cohort analysis visualization...")
    create_visualization(results, stats, OUTPUT_DIR / "phase8_cohort_analysis.png")

    # Summary
    print("\n" + "=" * 70)
    print("PHASE 8 COMPLETE")
    print("=" * 70)
    print(f"Cohort Size: {N_PATIENTS}")
    print(f"RL Mean Volume: {stats['rl_mean_final_volume_mm3']:.2f} ± {stats['rl_std_final_volume_mm3']:.2f} mm³")
    print(f"Stupp Mean Volume: {stats['stupp_mean_final_volume_mm3']:.2f} ± {stats['stupp_std_final_volume_mm3']:.2f} mm³")
    print(f"Paired t-test p-value: {stats['paired_t_test_p_value']:.4f}")
    print(f"Wilcoxon p-value: {stats['wilcoxon_p_value']:.4f}")
    print(f"Cohen's d: {stats['cohens_d']:.3f}")
    print(f"RL Progression-Free Rate: {stats['rl_progression_free_rate']:.1%}")
    print(f"Stupp Progression-Free Rate: {stats['stupp_progression_free_rate']:.1%}")
    print(f"McNemar's p-value: {stats['mcnemar_p_value']:.4f}")
    print(f"\nOutputs saved to {OUTPUT_DIR}/")
    print("  - phase8_cohort_metrics.json")
    print("  - phase8_cohort_analysis.png")


if __name__ == "__main__":
    main()