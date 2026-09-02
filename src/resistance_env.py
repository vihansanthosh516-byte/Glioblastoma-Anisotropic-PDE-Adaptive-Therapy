#!/usr/bin/env python3
"""Resistance-aware 3D anisotropic PDE environment for adaptive therapy.

Extends the FastPDESolver from 64_virtual_cohort_simulation with a resistant
subclone (Lotka-Volterra competition + mutation), following the validated
dual-clone model in Track B (44_adaptive_therapy.py):

    du_s/dt = div(D grad u_s) + rho_s u_s (1 - u_s - u_r) - kill(C) u_s - mutation
    du_r/dt = div(D grad u_r) + rho_r u_r (1 - u_s - u_r) + mutation

Resistant cells are NOT killed by chemotherapy (partial) and carry a fitness
cost (rho_r < rho_s). Over-dosing selects for resistance -> late escape. This
makes dose-sparing (adaptive therapy) valuable, which is where the trust-gated
hybrid is expected to win.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))

_spec = importlib.util.spec_from_file_location("vcs", ROOT / "src" / "64_virtual_cohort_simulation.py")
_vcs = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(_vcs)

FastPDESolver = _vcs.FastPDESolver
GbmTherapyEnv = _vcs.GbmTherapyEnv

# Resistance parameters (from Track B 44_adaptive_therapy.py)
MUTATION_RATE = 1e-4        # per division (raised so resistance manifests in-trial)
RHO_RESISTANT = 0.020        # fitness cost of resistance (slightly below rho)
RESISTANCE_KILL_FRACTION = 0.0  # resistant cells unaffected by chemo/RT


class ResistancePDESolver(FastPDESolver):
    """Dual-clone PDE solver: sensitive (u_s) + resistant (u_r)."""

    def __init__(
        self,
        grid_size=(16, 16, 16),
        dt_pde: float = 0.2,
        rho: float = 0.025,
        D_white: float = 0.0012,
        alpha_sens: float = 1.0,
        gamma_chemo: float = 0.08,
        alpha_rt: float = 0.045,
        mutation_rate: float = MUTATION_RATE,
        rho_resistant: float = RHO_RESISTANT,
        resistant_kill_fraction: float = RESISTANCE_KILL_FRACTION,
        initial_resistant_fraction: float = 1e-4,
        is_training: bool = False,
    ):
        super().__init__(
            grid_size=grid_size, dt_pde=dt_pde, rho=rho, D_white=D_white,
            alpha_sens=alpha_sens, gamma_chemo=gamma_chemo, alpha_rt=alpha_rt,
            is_training=is_training,
        )
        self.mutation_rate = mutation_rate
        self.rho_resistant = rho_resistant
        self.resistant_kill_fraction = resistant_kill_fraction
        self.initial_resistant_fraction = initial_resistant_fraction
        self.u_s = None
        self.u_r = None

    def reset(self, seed_center: Optional[Tuple[int, int, int]] = None):
        if seed_center is None:
            seed_center = (self.nx // 2, self.ny // 2, self.nz // 2)
        xx, yy, zz = np.mgrid[0:self.nx, 0:self.ny, 0:self.nz].astype(float)
        cx, cy, cz = seed_center
        r2 = ((xx - cx) ** 2 + (yy - cy) ** 2 + (zz - cz) ** 2) * self.dx ** 2
        base = _vcs.SEED_AMPLITUDE * np.exp(-r2 / (2 * _vcs.SEED_SIGMA_MM ** 2))
        self.u_s = base * (1.0 - self.initial_resistant_fraction)
        self.u_r = base * self.initial_resistant_fraction
        self.u = self.u_s + self.u_r
        self.step_count = 0
        self.chemo_tox = 0.0
        self.rad_tox = 0.0
        self.initial_volume = float(self.u.sum() * self.dx ** 3)
        self.prev_volume = self.initial_volume

    def pde_step_dual(self, kill_s: float) -> None:
        """Advance both clones one PDE step (shared anisotropic diffusion)."""
        div_s = self._aniso_divergence(self.u_s)
        div_r = self._aniso_divergence(self.u_r)
        total = self.u_s + self.u_r

        growth_s = self.rho * self.u_s * (1.0 - total)
        growth_r = self.rho_resistant * self.u_r * (1.0 - total)
        mutation = self.mutation_rate * np.maximum(growth_s, 0.0)

        kill_r = kill_s * self.resistant_kill_fraction

        u_s_new = self.u_s + self.dt * (div_s + growth_s - kill_s * self.u_s - mutation)
        u_r_new = self.u_r + self.dt * (div_r + growth_r + mutation - kill_r * self.u_r)

        self.u_s = np.clip(u_s_new, 0.0, _vcs.K_CARRY)
        self.u_r = np.clip(u_r_new, 0.0, _vcs.K_CARRY)
        self.u = self.u_s + self.u_r

    def rl_step(self, action: int) -> Dict[str, float]:
        kill = 0.0
        if action == 1:
            kill = self.gamma_chemo * self.alpha_sens
            self.chemo_tox += _vcs.CHEMO_TOX_PER_RL_STEP
        elif action == 2:
            kill = self.alpha_rt * self.alpha_sens
            self.rad_tox += _vcs.RAD_TOX_PER_RL_STEP
        elif action == 3:
            kill = (self.gamma_chemo + self.alpha_rt) * self.alpha_sens
            self.chemo_tox += _vcs.COMBO_TOX_PER_RL_STEP * 0.5
            self.rad_tox += _vcs.COMBO_TOX_PER_RL_STEP * 0.5

        n_sub = 1 if self.is_training else 5
        for _ in range(n_sub):
            self.pde_step_dual(kill)

        self.step_count += 1
        volume = float(self.u.sum() * self.dx ** 3)
        self.prev_volume = volume
        return {
            "volume_mm3": volume,
            "u_max": float(self.u.max()),
            "norm_volume": volume / max(self.initial_volume, 1e-6),
            "delta_volume": volume - self.prev_volume,
            "chemo_tox": self.chemo_tox,
            "rad_tox": self.rad_tox,
            "resistant_fraction": float(self.u_r.sum() / max(self.u.sum(), 1e-12)),
        }


class ResistanceGbmEnv:
    """Gymnasium-like env for the dual-clone solver (same interface as GbmTherapyEnv)."""

    def __init__(self, solver: ResistancePDESolver, reward_weights: Dict[str, float], max_steps: int = 90):
        self.solver = solver
        self.max_steps = max_steps
        self.trajectory = []
        self.reward_weights = reward_weights

    def reset(self, seed=None, options=None):
        self.solver.reset()
        self.trajectory = []
        return self.solver.get_observation(), {}

    def step(self, action: int):
        result = self.solver.rl_step(action)
        obs = self.solver.get_observation()
        norm_vol = result["norm_volume"]
        u_max = result["u_max"]

        # Per-step reward: volume + density control, with a light resistance term.
        lambda_vol = self.reward_weights["lambda_vol"]
        lambda_den = self.reward_weights["lambda_den"]
        lambda_tox = self.reward_weights["lambda_tox"]
        lambda_res = self.reward_weights.get("lambda_res", 8.0)

        reward = (
            -lambda_vol * norm_vol
            - lambda_den * u_max
            - lambda_tox * (1.0 if action > 0 else 0.0)
            - lambda_res * result["resistant_fraction"]
        )
        if result["delta_volume"] > 0:
            reward += self.reward_weights.get("lambda_shrink", 100.0) * max(
                result["delta_volume"] / max(self.solver.initial_volume, 1e-6), 0.0
            )
        terminated = self.solver.step_count >= self.max_steps
        if terminated:
            # Endpoint reward: penalize residual volume AND resistance burden.
            # This is the biologically honest trial endpoint: a tumor that is
            # 99% resistant is a recurrence waiting to happen, even if small.
            res_frac = result["resistant_fraction"]
            reward -= self.reward_weights.get("lambda_res_end", 30.0) * res_frac
            if norm_vol < 0.01:
                reward += self.reward_weights.get("lambda_clear", 200.0)

        self.trajectory.append({
            "step": self.solver.step_count,
            "action": int(action),
            "volume_mm3": result["volume_mm3"],
            "u_max": u_max,
            "reward": reward,
            "resistant_fraction": result["resistant_fraction"],
        })
        return obs, float(reward), terminated, False, {}


def make_resistance_env(rho: float, D_white: float, grid_size=(16, 16, 16)) -> ResistanceGbmEnv:
    solver = ResistancePDESolver(grid_size=grid_size, dt_pde=0.2, rho=rho, D_white=D_white, is_training=False)
    env = ResistanceGbmEnv(solver, _vcs.RL_REWARD_WEIGHTS)
    env.reset()
    return env


if __name__ == "__main__":
    env = make_resistance_env(0.025, 0.0012, (16, 16, 16))
    # Run Stupp-like protocol to verify resistance develops under continuous dosing
    for step in range(90):
        day = step + 1
        action = 3 if 20 <= day < 50 else (1 if (day % 28) < 5 and day >= 50 else 0)
        obs, reward, term, trunc, _ = env.step(action)
        if term or trunc:
            break
    print("Final volume:", env.solver.u.sum() * env.solver.dx ** 3)
    print("Resistant fraction:", env.solver.u_r.sum() / max(env.solver.u.sum(), 1e-12))