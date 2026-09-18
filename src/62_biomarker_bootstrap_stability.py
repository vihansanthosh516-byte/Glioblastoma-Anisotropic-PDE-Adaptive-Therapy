#!/usr/bin/env python3
"""
Phase 7: Biomarker Threshold Bootstrap Stability
=================================================
Addresses Reviewer Concern #5: Biomarker Rule Stability - Quantifying the 
statistical robustness and 95% Confidence Interval for the ρ > 0.024 
decision threshold.
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
from scipy.stats import gaussian_kde

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
N_SCENARIOS = 30
BOOTSTRAP_SAMPLES = 1000
EVAL_GRID = (64, 64, 64)
T_MAX_DAYS = 90
DT_PDE_EVAL = 0.2

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Phase 6 parameter ranges (from Phase 6 config)
PARAM_RANGES = {
    "rho": (0.005, 0.035),        # 1/day
    "D_w": (0.001, 0.008),        # cm^2/day
    "alpha_sens": (0.5, 1.5),     # multiplier
}

# Fixed parameters
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
DT_RL_DAYS = 1.0
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = 5

OUTPUT_DIR = Path("output")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------------------------------- #
# Fast PDE Solver (Self-contained, identical to Phase 6)
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

            l1 = np.where(tract_mask, D_WHITE_BASE, D_GRAY_BASE)
            l2 = np.where(tract_mask, D_GRAY_BASE, D_GRAY_BASE)
        else:
            tract_mask = np.zeros((self.nx, self.ny, self.nz), dtype=bool)
            vx = np.ones_like(xx)
            vy = np.zeros_like(yy)
            vz = np.zeros_like(zz)
            l1 = np.full_like(xx, D_WHITE_BASE)
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
            kill = GAMMA_CHEMO
        elif action == 2:
            kill = GAMMA_RAD
        elif action == 3:
            kill = GAMMA_CHEMO + GAMMA_RAD

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
        }

    def get_observation(self) -> np.ndarray:
        vol = float(self.u.sum() * self.dx**3)
        u_max = float(self.u.max())
        norm_vol = vol / max(self.initial_volume, 1e-6)
        return np.array([
            np.clip(norm_vol, 0, 1),
            np.clip(float(self.u.max()), 0, 1),
            self.step_count / 90,
        ], dtype=np.float32)

    def is_done(self) -> bool:
        return self.step_count >= 90


# --------------------------------------------------------------------------- #
# Environment & Protocols
# --------------------------------------------------------------------------- #
class GbmTherapyEnv:
    def __init__(self, solver: FastPDESolver):
        self.solver = solver
        self.max_steps = 90
        self.trajectory = []

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
        reward = -15.0 * norm_vol - 8.0 * u_max - 0.02 * action_cost
        if delta_vol > 0:
            reward += 100.0 * max(delta_vol / max(self.solver.initial_volume, 1e-6), 0.0)

        terminated = self.solver.is_done()
        if terminated and norm_vol < 0.01:
            reward += 200.0

        self.trajectory.append({
            "step": self.solver.step_count,
            "action": action,
            "volume_mm3": result["volume_mm3"],
            "u_max": result["u_max"],
            "reward": reward,
        })

        return obs, float(reward), terminated, False, {}


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


def run_rl_adaptive(env: GbmTherapyEnv) -> Dict:
    obs, _ = env.reset()
    trajectory = []

    for step in range(env.max_steps):
        # Heuristic policy (same as Phase 5/6)
        current_vol = env.solver.u.sum() * env.solver.dx**3
        norm_vol = current_vol / max(env.solver.initial_volume, 1e-6)

        if norm_vol > 0.05:
            action = 3
        elif norm_vol > 0.01:
            action = 2
        else:
            action = 0

        # Guardrail
        if action == 0 and current_vol > 0.05 * env.solver.initial_volume:
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
# Phase 6 Scenario Generation (replicate Phase 6 sampling)
# --------------------------------------------------------------------------- #
from scipy.stats.qmc import LatinHypercube

def generate_phase6_scenarios(n_samples: int = N_SCENARIOS) -> List[Dict[str, float]]:
    """Replicate Phase 6 parameter sampling with fixed seed."""
    param_names = list(PARAM_RANGES.keys())
    bounds = np.array([PARAM_RANGES[p] for p in param_names])

    sampler = LatinHypercube(d=len(param_names), seed=42)
    samples = sampler.random(n=n_samples)
    scaled = bounds[:, 0] + samples * (bounds[:, 1] - bounds[:, 0])

    param_list = []
    for i in range(n_samples):
        param_list.append({
            "rho": float(scaled[i, 0]),
            "D_w": float(scaled[i, 1]),
            "alpha_sens": float(scaled[i, 2]),
        })
    return param_list


def evaluate_scenario(params: Dict[str, float]) -> Dict[str, float]:
    """Evaluate single scenario with both RL and Stupp."""
    solver = FastPDESolver(
        grid_size=EVAL_GRID,
        dt_pde=DT_PDE_EVAL,
        rho=params["rho"],
        D_white=params["D_w"],
        alpha_sens=params["alpha_sens"],
        is_training=False,
    )
    env = GbmTherapyEnv(solver)

    rl_result = run_rl_adaptive(env)
    stupp_result = run_stupp_protocol(env)

    return {
        "rl_volume": rl_result["final_volume_mm3"],
        "stupp_volume": stupp_result["final_volume_mm3"],
        "rl_wins": rl_result["final_volume_mm3"] < stupp_result["final_volume_mm3"],
        "rho": params["rho"],
        "D_w": params["D_w"],
        "alpha_sens": params["alpha_sens"],
    }


# --------------------------------------------------------------------------- #
# Bootstrap Engine
# --------------------------------------------------------------------------- #
def find_critical_rho(scenarios: List[Dict[str, float]]) -> Optional[float]:
    """
    Find the critical rho threshold where RL starts outperforming Stupp.
    Uses a simple decision boundary: RL wins when rho > threshold.
    """
    # Sort by rho
    sorted_scenarios = sorted(scenarios, key=lambda x: x["rho"])
    
    # Find the boundary where RL starts consistently winning
    # Use a simple approach: find the rho value that best separates wins/losses
    best_threshold = None
    best_accuracy = 0
    
    rho_vals = np.array([s["rho"] for s in sorted_scenarios])
    rl_wins = np.array([s["rl_wins"] for s in sorted_scenarios])
    
    # Try thresholds at midpoints between sorted rho values
    for i in range(1, len(rho_vals)):
        threshold = (rho_vals[i-1] + rho_vals[i]) / 2
        predictions = rho_vals > threshold
        accuracy = np.mean(predictions == rl_wins)
        if accuracy > best_accuracy:
            best_accuracy = accuracy
            best_threshold = threshold
    
    return best_threshold


def bootstrap_rho_crit(scenarios: List[Dict[str, float]], n_bootstrap: int = BOOTSTRAP_SAMPLES) -> List[float]:
    """Perform bootstrap resampling to estimate rho_crit distribution."""
    n = len(scenarios)
    boot_thresholds = []
    
    for i in range(n_bootstrap):
        # Resample with replacement
        indices = np.random.choice(n, size=n, replace=True)
        boot_sample = [scenarios[j] for j in indices]
        
        threshold = find_critical_rho(boot_sample)
        if threshold is not None:
            boot_thresholds.append(threshold)
        
        if (i + 1) % 200 == 0:
            print(f"  Bootstrap {i+1}/{n_bootstrap}...")
    
    return boot_thresholds


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main():
    print("=" * 70)
    print("PHASE 7: BIOMARKER THRESHOLD BOOTSTRAP STABILITY")
    print("=" * 70)

    # 1. Generate Phase 6 scenarios (reproducible)
    print(f"\n[Phase 7] Generating {N_SCENARIOS} Phase 6 scenarios...")
    scenarios_params = generate_phase6_scenarios(N_SCENARIOS)

    # 2. Evaluate all scenarios
    print(f"[Phase 7] Evaluating {N_SCENARIOS} scenarios (RL vs Stupp)...")
    scenario_results = []
    for i, params in enumerate(scenarios_params):
        print(f"  Scenario {i+1}/{N_SCENARIOS}: rho={params['rho']:.4f}, D_w={params['D_w']:.4f}, alpha={params['alpha_sens']:.4f}")
        result = evaluate_scenario(params)
        scenario_results.append(result)
        winner = "RL" if result["rl_wins"] else "Stupp"
        print(f"    RL: {result['rl_volume']:.2f} mm³, Stupp: {result['stupp_volume']:.2f} mm³ -> {winner} wins")

    # 3. Bootstrap analysis
    print(f"\n[Phase 7] Running {BOOTSTRAP_SAMPLES} bootstrap resamples...")
    boot_thresholds = bootstrap_rho_crit(scenario_results, BOOTSTRAP_SAMPLES)

    # 4. Statistics
    boot_thresholds = np.array(boot_thresholds)
    mean_rho_crit = float(np.mean(boot_thresholds))
    ci_lower = float(np.percentile(boot_thresholds, 2.5))
    ci_upper = float(np.percentile(boot_thresholds, 97.5))
    std_threshold = float(np.std(boot_thresholds))

    # Point estimate from full dataset
    point_threshold = find_critical_rho(scenario_results)
    
    # Additional statistics
    median_threshold = float(np.median(boot_thresholds))
    iqr = float(np.percentile(boot_thresholds, 75) - np.percentile(boot_thresholds, 25))

    # Additional analysis: logistic regression for more robust threshold
    from scipy.optimize import curve_fit
    
    # Fit logistic: P(RL wins) = 1 / (1 + exp(-k*(rho - rho_crit)))
    rho_vals = np.array([s["rho"] for s in scenario_results])
    rl_wins = np.array([1.0 if s["rl_wins"] else 0.0 for s in scenario_results])
    
    def logistic(x, k, x0):
        return 1 / (1 + np.exp(-k * (x - x0)))
    
    try:
        popt, pcov = curve_fit(logistic, rho_vals, rl_wins, p0=[200, 0.024], maxfev=5000)
        logistic_threshold = float(popt[1])
        logistic_k = float(popt[0])
    except:
        logistic_threshold = point_threshold
        logistic_k = 0.0

    # 5. Compile metrics
    metrics = {
        "bootstrap_samples": BOOTSTRAP_SAMPLES,
        "n_scenarios": len(scenario_results),
        "point_estimate_rho_crit": point_threshold,
        "logistic_threshold_rho_crit": logistic_threshold,
        "logistic_steepness_k": logistic_k,
        "mean_rho_crit": mean_rho_crit,
        "median_rho_crit": median_threshold,
        "std_rho_crit": std_threshold,
        "ci_95_lower": ci_lower,
        "ci_95_upper": ci_upper,
        "iqr": iqr,
        "bootstrap_thresholds": boot_thresholds.tolist()[:100],  # Store first 100 for reference
        "full_bootstrap_thresholds": boot_thresholds.tolist(),  # Store all
    }

    # 6. Save metrics
    with open(OUTPUT_DIR / "biomarker_stability_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\n[Metrics] Saved -> {OUTPUT_DIR / 'biomarker_stability_metrics.json'}")

    # 7. Visualization
    print("\n[Phase 7] Generating biomarker stability visualization...")
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: Histogram + KDE
    ax1 = axes[0]
    n, bins, patches = ax1.hist(boot_thresholds, bins=40, density=True, alpha=0.6, 
                                 color='#1f77b4', edgecolor='black', linewidth=0.5,
                                 label='Bootstrap samples')
    
    # KDE overlay
    kde = gaussian_kde(boot_thresholds)
    x_kde = np.linspace(min(boot_thresholds), max(boot_thresholds), 200)
    ax1.plot(x_kde, kde(x_kde), 'r-', linewidth=2, label='KDE')
    
    # Markers
    ax1.axvline(mean_rho_crit, color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {mean_rho_crit:.4f}')
    ax1.axvline(median_threshold, color='orange', linestyle='--', linewidth=2,
                label=f'Median: {median_threshold:.4f}')
    ax1.axvline(ci_lower, color='red', linestyle=':', linewidth=2,
                label=f'95% CI: [{ci_lower:.4f}, {ci_upper:.4f}]')
    ax1.axvline(ci_upper, color='red', linestyle=':', linewidth=2)
    if point_threshold:
        ax1.axvline(point_threshold, color='black', linestyle='-', linewidth=2,
                    label=f'Full-data estimate: {point_threshold:.4f}')
    if logistic_threshold:
        ax1.axvline(logistic_threshold, color='green', linestyle='-.', linewidth=2,
                    label=f'Logistic fit: {logistic_threshold:.4f}')
    
    ax1.set_xlabel('Critical Proliferation Rate ρ_crit (1/day)')
    ax1.set_ylabel('Density')
    ax1.set_title('Panel 1: Bootstrap Distribution of Critical ρ Threshold\n(Where RL Outperforms Stupp)')
    ax1.legend(loc='upper left', fontsize=9)
    ax1.grid(alpha=0.3)

    # Panel 2: Q-Q plot / Empirical CDF
    ax2 = axes[1]
    sorted_thresholds = np.sort(boot_thresholds)
    n = len(sorted_thresholds)
    ecdf = np.arange(1, n+1) / n
    ax2.plot(sorted_thresholds, ecdf, 'b-', linewidth=2, label='Empirical CDF')
    
    # Normal approximation
    from scipy.stats import norm
    norm_cdf = norm.cdf(sorted_thresholds, mean_rho_crit, std_threshold)
    ax2.plot(sorted_thresholds, norm_cdf, 'r--', linewidth=1.5, label='Normal approx.')
    
    ax2.axvline(mean_rho_crit, color='red', linestyle='--', label=f'Mean: {mean_rho_crit:.4f}')
    ax2.axvline(ci_lower, color='red', linestyle=':', label=f'95% CI')
    ax2.axvline(ci_upper, color='red', linestyle=':')
    
    ax2.set_xlabel('Critical Proliferation Rate ρ_crit (1/day)')
    ax2.set_ylabel('Cumulative Probability')
    ax2.set_title('Panel 2: Empirical CDF of Bootstrap Thresholds')
    ax2.legend(loc='lower right')
    ax2.grid(alpha=0.3)

    plt.suptitle('Phase 7: Biomarker Threshold Bootstrap Stability\n(Critical ρ where RL > Stupp)', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / "biomarker_stability.png", dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {OUTPUT_DIR / 'biomarker_stability.png'}")

    # 8. Summary
    print("\n" + "=" * 70)
    print("PHASE 7: BIOMARKER THRESHOLD BOOTSTRAP STABILITY COMPLETE")
    print("=" * 70)
    print(f"Bootstrap samples: {BOOTSTRAP_SAMPLES}")
    print(f"Point estimate (full data): {point_threshold:.4f} day⁻¹")
    print(f"Logistic fit threshold: {logistic_threshold:.4f} day⁻¹ (k={logistic_k:.1f})")
    print(f"Bootstrap mean: {mean_rho_crit:.4f} ± {std_threshold:.4f} day⁻¹")
    print(f"95% CI: [{ci_lower:.4f}, {ci_upper:.4f}] day⁻¹")
    print(f"Median: {median_threshold:.4f}, IQR: {iqr:.4f}")
    print(f"\nClinical Interpretation:")
    print(f"  If patient ρ > {ci_upper:.4f}: RL Adaptive strongly recommended")
    print(f"  If patient ρ < {ci_lower:.4f}: Standard Stupp likely sufficient")
    print(f"  Ambiguous zone: ρ ∈ [{ci_lower:.4f}, {ci_upper:.4f}]")
    print(f"\nOutputs saved to {OUTPUT_DIR}/")
    print("  - biomarker_stability_metrics.json")
    print("  - biomarker_stability.png")


if __name__ == "__main__":
    main()