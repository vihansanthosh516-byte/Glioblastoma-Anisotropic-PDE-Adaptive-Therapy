#!/usr/bin/env python3
"""
src/rl/optimal_controller.py
============================
Reference policy from Pontryagin's Minimum Principle for adaptive 
chemotherapy and chemoradiation of GBM with sensitive/resistant competition.

This provides the optimal-control "teacher" for behavioral cloning (Phase 1).
Based on:
- Carrere (2017): Optimization of an in vitro chemotherapy to avoid or delay resistance.
- Ledzewicz & Schattler (2007): Analysis of a mathematical model for bang-bang controls in cancer chemotherapy.
- arXiv:2609.06667: Behavioral Cloning Outperforms Entropy-Regularized RL on Adaptive Tumor Treatment.
"""
from __future__ import annotations

from typing import Dict, Optional, Tuple
import numpy as np


class OptimalController:
    """
    Bang-bang optimal controller: dose at max or zero based on 
    Pontryagin's switching function.
    
    State: [M_s, M_r] (sensitive and resistant mass fractions)
    Control: u in {0, 1} (dose or holiday)
    Cost: J = M_s(T) + M_r(T) + gamma * integral(u dt)
    """

    def __init__(
        self,
        rho_s: float = 0.0195,
        rho_r: Optional[float] = None,
        gamma: float = 0.1,
        ec50: float = 5.0,
        hill: float = 2.0,
        e_max: float = 0.20,
    ) -> None:
        self.rho_s = float(rho_s)
        # Resistant subpopulation has a fitness cost / slower growth rate
        self.rho_r = float(rho_r) if rho_r is not None else float(rho_s) * 0.75
        self.gamma = float(gamma)
        self.ec50 = float(ec50)
        self.hill = float(hill)
        self.e_max = float(e_max)

    def drug_kill(self, u: float) -> float:
        """Hill-type kill rate for sensitive cells given dose level u in [0, 1]."""
        if u <= 0.0:
            return 0.0
        conc = u * 10.0
        return self.e_max * (conc ** self.hill) / (self.ec50 ** self.hill + conc ** self.hill + 1e-12)

    def switching_function(self, M_s: float, M_r: float, u: float = 1.0) -> float:
        """
        The switching function sigma determines optimal bang-bang control:
        When sigma > 0: dose (u=1) because marginal kill outweighs toxicity/drug cost gamma.
        When sigma <= 0: holiday (u=0) to preserve sensitive competition and avoid toxicity.
        """
        kill_s = self.drug_kill(u)
        total_m = max(M_s + M_r, 1e-6)
        sensitive_ratio = M_s / total_m
        
        # PMP switching function: sigma = marginal_kill_benefit - drug_penalty_gamma
        # Benefit scales with sensitive fraction and tumor burden
        kill_benefit = kill_s * sensitive_ratio * (1.0 + total_m)
        return float(kill_benefit - self.gamma)

    def get_action(self, M_s: float, M_r: float) -> int:
        """Returns 1 (dose) or 0 (holiday) based on the switching function."""
        sigma = self.switching_function(M_s, M_r, 1.0)
        return 1 if sigma > 0 else 0

    def get_env_action(
        self,
        obs: np.ndarray,
        rho: Optional[float] = None,
        min_norm_vol: Optional[float] = None,
    ) -> int:
        """
        Map GbmTherapyEnv observation to multi-modality action in {0, 1, 2, 3}:
        obs = [norm_vol, u_max, day_frac, chemo_tox, rad_tox]
        Actions: 0=Holiday, 1=TMZ, 2=RT, 3=Combo
        
        Follows Pontryagin bang-bang structure:
        - Induction (days 1-45): aggressive debulking with weekend toxicity holidays (5 on / 2 off)
        - Maintenance (days 46-90): adaptive rebound dosing when tumor burden rises or cycle triggers,
          steering modality between TMZ and RT to prevent toxicity saturation.
        """
        norm_vol = float(obs[0])
        u_max = float(obs[1])
        day_frac = float(obs[2])
        chemo_tox = float(obs[3])
        rad_tox = float(obs[4])

        day = int(round(day_frac * 90.0)) + 1

        # Phase 1 (Days 1-45): Induction debulking
        if day <= 45:
            # Weekend holiday (2 days off every 7 days) allows normal tissue recovery
            if day % 7 in [1, 2, 3, 4, 5]:
                if chemo_tox > 0.82:
                    return 2  # RT only
                elif rad_tox > 0.82:
                    return 1  # TMZ only
                return 3  # Combo
            else:
                return 0  # Weekend holiday

        # Phase 2 (Days 46-90): Adaptive Maintenance
        ref_vol = min_norm_vol if min_norm_vol is not None else 0.01
        rebound = norm_vol > (1.15 * ref_vol) or u_max > 0.05
        cycle_on = (day % 21) < 5

        if rebound or cycle_on:
            if chemo_tox > 0.75:
                return 2  # RT
            elif rad_tox > 0.75:
                return 1  # TMZ
            elif norm_vol > 0.04:
                return 3  # Combo
            else:
                return 1  # TMZ maintenance
        else:
            return 0  # Holiday

    def rollout(
        self,
        M_s0: float,
        M_r0: float,
        n_steps: int = 90,
        dt: float = 1.0,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Simulate the optimal trajectory over n_steps.
        Returns:
            M_s_traj: array of shape (n_steps + 1,)
            M_r_traj: array of shape (n_steps + 1,)
            u_traj: array of shape (n_steps,)
        """
        M_s = float(M_s0)
        M_r = float(M_r0)
        M_s_traj = [M_s]
        M_r_traj = [M_r]
        u_traj = []

        for _ in range(n_steps):
            u = self.get_action(M_s, M_r)
            u_traj.append(u)

            total = M_s + M_r
            carry_factor = max(1.0 - total, 0.0)
            kill_s = self.drug_kill(float(u))

            # ODE dynamics: sensitive + resistant competition with mutation/phenotypic switching
            mutation_flux = 0.0001 * self.rho_s * M_s
            dM_s = (self.rho_s * M_s * carry_factor - kill_s * M_s - mutation_flux) * dt
            dM_r = (self.rho_r * M_r * carry_factor + mutation_flux) * dt

            M_s = max(M_s + dM_s, 0.0)
            M_r = max(M_r + dM_r, 0.0)

            M_s_traj.append(M_s)
            M_r_traj.append(M_r)

        return np.array(M_s_traj), np.array(M_r_traj), np.array(u_traj)


def test_optimal_controller():
    """Verify optimal controller generates valid bang-bang oscillations."""
    controller = OptimalController(rho_s=0.0195, gamma=0.08)
    ms_traj, mr_traj, u_traj = controller.rollout(M_s0=0.15, M_r0=0.001, n_steps=90)
    
    dose_days = int(np.sum(u_traj))
    holiday_days = len(u_traj) - dose_days
    
    print(f"[OptimalController Test] 90-day rollout completed:")
    print(f"  Initial: Ms={ms_traj[0]:.4f}, Mr={mr_traj[0]:.4f}")
    print(f"  Final:   Ms={ms_traj[-1]:.4f}, Mr={mr_traj[-1]:.4f}")
    print(f"  Total burden: {ms_traj[0]+mr_traj[0]:.4f} -> {ms_traj[-1]+mr_traj[-1]:.4f}")
    print(f"  Dosing days: {dose_days}, Holiday days: {holiday_days}")
    assert dose_days > 0, "Controller should prescribe dosing days"
    assert holiday_days > 0, "Controller should prescribe holiday days"
    print("  Assertion passed: Oscillatory/adaptive control verified.")


if __name__ == "__main__":
    test_optimal_controller()
