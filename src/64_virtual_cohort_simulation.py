#!/usr/bin/env python3
"""
Script 64: 20-Patient Virtual Cohort at an EQUAL DRUG BUDGET
============================================================
Every arm gets the same total drug AUC (Stupp's) and the same seed tumour.
Arms come from src/rl/equal_budget_arms.py (shared with 60, 62, 65):
  Stupp | budgeted 59 heuristic | RL = PPO policy (script 66) | DAgger oracle (script 66)

Physics is this script's own FastPDESolver. Virtual patients differ in rho,
D_white and in their own kill rates (gamma_chemo, alpha_rt), so the most
drug-efficient action differs between patients. Neither script-66 policy saw
patient-specific kill rates in training.

Outputs:
  - output/phase8_cohort_metrics.json
  - output/phase8_cohort_analysis.png
"""
from __future__ import annotations

import json
import sys
import time
import warnings
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from rl.equal_budget_arms import ARMS, STUPP_BUDGET, evaluate_arms, s66, validate_against_numpy  # noqa: E402

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
N_PATIENTS = 20
COHORT_SEED = 42
EVAL_GRID = (64, 64, 64)
T_MAX_DAYS = 90
DT_RL_DAYS = 1.0
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = 5
PROGRESSION_THRESHOLD_MM3 = 500.0
BOOTSTRAP_SAMPLES = 2000

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

OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

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


# --------------------------------------------------------------------------- #
# Equal-budget cohort evaluation
# --------------------------------------------------------------------------- #
ARM_LABELS = {"stupp": "Stupp", "heuristic_budgeted": "Budgeted heuristic (59)",
              "ppo": "RL: PPO (66)", "dagger_oracle": "DAgger oracle (66)"}
ACTION_NAMES = ["off", "chemo", "rad", "combo"]


def kill_row(p: Dict[str, float]) -> List[float]:
    """Kill rate per action exactly as FastPDESolver.rl_step computes it (alpha_sens = 1)."""
    return [0.0, p["gamma_chemo"], p["alpha_rt"], p["gamma_chemo"] + p["alpha_rt"]]


def make_solver(p: Dict[str, float]) -> "FastPDESolver":
    return FastPDESolver(grid_size=EVAL_GRID, dt_pde=DT_PDE_EVAL, rho=p["rho"], D_white=p["D_white"],
                         alpha_sens=1.0, gamma_chemo=p["gamma_chemo"], alpha_rt=p["alpha_rt"])


def arm_stats(arm: Dict[str, list], ref: Dict[str, list], rng: np.random.Generator) -> Dict[str, Any]:
    f, fs = np.array(arm["final_volume_mm3"]), np.array(ref["final_volume_mm3"])
    b, bs = np.array(arm["mean_burden_mm3"]), np.array(ref["mean_burden_mm3"])
    v0 = np.array(arm["initial_volume_mm3"])
    lr = np.log(f / fs)
    boot = np.array([lr[rng.integers(0, len(lr), len(lr))].mean() for _ in range(BOOTSTRAP_SAMPLES)])
    out = {
        "mean_final_volume_mm3": float(f.mean()), "std_final_volume_mm3": float(f.std(ddof=1)),
        "median_final_volume_mm3": float(np.median(f)),
        "mean_drug_auc": float(np.mean(arm["drug_auc"])),
        "drug_auc_ratio_vs_stupp": float(np.mean(arm["drug_auc"]) / np.mean(ref["drug_auc"])),
        "win_rate_vs_stupp_pct": float((f < fs).mean() * 100),
        "mean_log_ratio_final_vs_stupp": float(lr.mean()),
        "mean_log_ratio_final_vs_stupp_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
        "burden_win_rate_vs_stupp_pct": float((b < bs).mean() * 100),
        "mean_log_ratio_burden_vs_stupp": float(np.log(b / bs).mean()),
        "n_progressed_v90_gt_500mm3": int((f > PROGRESSION_THRESHOLD_MM3).sum()),
        "n_grew_over_90_days": int((f > v0).sum()),
        "wilcoxon_p_value": None, "paired_t_p_value_log": None, "cohens_d_log": None,
    }
    if np.any(lr != 0):
        out["wilcoxon_p_value"] = float(stats.wilcoxon(np.log(f), np.log(fs)).pvalue)
        out["paired_t_p_value_log"] = float(stats.ttest_rel(np.log(f), np.log(fs)).pvalue)
        if lr.std(ddof=1) > 0:
            out["cohens_d_log"] = float(lr.mean() / lr.std(ddof=1))
    return out


def create_visualization(cohort, res, arm_summary, output_path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    stupp_f = np.array(res["stupp"]["final_volume_mm3"])
    ids = np.arange(len(cohort))

    ax = axes[0, 0]
    for i, arm in enumerate(res):
        ax.bar(ids + (i - 1.5) * 0.2, res[arm]["final_volume_mm3"], 0.2, label=ARM_LABELS[arm])
    ax.set_xlabel("patient")
    ax.set_ylabel("final volume (mm^3)")
    ax.set_title("Panel 1: Per-patient day-90 volume (equal drug budget)")
    ax.set_xticks(ids)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[0, 1]
    arms = [a for a in res if a != "stupp"]
    for i, arm in enumerate(arms):
        lr = np.log(np.array(res[arm]["final_volume_mm3"]) / stupp_f)
        ax.scatter(np.full(len(lr), i) + np.random.default_rng(i).uniform(-0.12, 0.12, len(lr)), lr, s=30,
                   edgecolors="black", alpha=0.8)
        lo, hi = arm_summary[arm]["mean_log_ratio_final_vs_stupp_ci95"]
        ax.errorbar(i + 0.3, lr.mean(), yerr=[[lr.mean() - lo], [hi - lr.mean()]], fmt="D", color="black", capsize=5)
        ax.text(i, max(lr) + 0.02, f"win {arm_summary[arm]['win_rate_vs_stupp_pct']:.0f}%",
                ha="center", fontsize=9)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels([ARM_LABELS[a] for a in arms])
    ax.set_ylabel("log(V_arm / V_Stupp) at day 90   (< 0 beats Stupp)")
    ax.set_title("Panel 2: Paired effect vs Stupp (mean + bootstrap 95% CI)")
    ax.grid(alpha=0.3, axis="y")

    ax = axes[1, 0]
    eff = np.array([[k / i for k, i in zip(kill_row(p)[1:], [0.5, 0.75, 1.0])] for p in cohort])
    best = eff.argmax(1) + 1
    lr_ppo = np.log(np.array(res["ppo"]["final_volume_mm3"]) / stupp_f)
    for a, col in ((1, "#fdae61"), (2, "#abd9e9"), (3, "#2c7bb6")):
        m = best == a
        if m.any():
            ax.scatter(np.array([p["rho"] for p in cohort])[m], lr_ppo[m], color=col, s=60, edgecolors="black",
                       label=f"most kill per unit drug: {ACTION_NAMES[a]} (n={m.sum()})")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("PPO: log(V_PPO / V_Stupp)")
    ax.set_title("Panel 3: PPO effect by patient drug-efficiency profile")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    days = np.arange(1, T_MAX_DAYS + 1)
    for arm in res:
        V = np.array(res[arm]["volume_history_mm3"])
        ax.plot(days, np.exp(np.log(V).mean(0)), lw=2, label=ARM_LABELS[arm])
    ax.set_yscale("log")
    ax.set_xlabel("day")
    ax.set_ylabel("geometric-mean volume (mm^3)")
    ax.set_title("Panel 4: Cohort trajectories")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, which="both")

    plt.suptitle(f"Script 64: {len(cohort)}-patient virtual cohort at equal drug budget (AUC = {STUPP_BUDGET:g})",
                 fontsize=15, fontweight="bold")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {output_path}")


def main():
    print("=" * 70)
    print("SCRIPT 64: VIRTUAL COHORT AT EQUAL DRUG BUDGET")
    print("=" * 70)
    cohort = generate_virtual_cohort(N_PATIENTS, seed=COHORT_SEED)
    kills = [kill_row(p) for p in cohort]
    val_err = validate_against_numpy(make_solver(cohort[0]), kills[0])
    assert val_err < 1e-9, f"batched solver disagrees with numpy solver: {val_err}"

    t0 = time.time()
    res = evaluate_arms([lambda p=p: make_solver(p) for p in cohort], kills, chunk=N_PATIENTS)
    print(f"  evaluated {len(ARMS)} arms x {len(cohort)} patients in {time.time() - t0:.0f}s (validation err {val_err:.1e})")

    rng = np.random.default_rng(64)
    arm_summary = {arm: arm_stats(r, res["stupp"], rng) for arm, r in res.items()}
    patients = []
    for i, p in enumerate(cohort):
        eff = {ACTION_NAMES[a]: kills[i][a] / s66.s59.ACTION_INTENSITY[a] for a in (1, 2, 3)}
        row = {"patient_id": int(p["patient_id"]), "params": {k: float(v) for k, v in p.items() if k != "patient_id"},
               "kill_per_unit_drug": eff, "most_efficient_action": max(eff, key=eff.get)}
        for arm, r in res.items():
            acts = r["action_history"][i]
            row[arm] = {"final_volume_mm3": r["final_volume_mm3"][i], "drug_auc": r["drug_auc"][i],
                        "mean_burden_mm3": r["mean_burden_mm3"][i],
                        "action_days": {ACTION_NAMES[a]: int(acts.count(a)) for a in range(4)},
                        "beats_stupp": bool(r["final_volume_mm3"][i] < res["stupp"]["final_volume_mm3"][i])}
        patients.append(row)

    ppo = arm_summary["ppo"]
    metrics = {
        "framework": "equal drug budget; arms from src/rl/equal_budget_arms.py",
        "budget_drug_auc": STUPP_BUDGET,
        "cohort_size": len(cohort),
        "cohort_seed": COHORT_SEED,
        "rl_column": "ppo",
        "validation_max_rel_err_vs_numpy_solver": val_err,
        "arms": arm_summary,
        # flat keys kept for script 65
        "rl_mean_final_volume_mm3": ppo["mean_final_volume_mm3"],
        "rl_std_final_volume_mm3": ppo["std_final_volume_mm3"],
        "stupp_mean_final_volume_mm3": arm_summary["stupp"]["mean_final_volume_mm3"],
        "stupp_std_final_volume_mm3": arm_summary["stupp"]["std_final_volume_mm3"],
        "rl_win_rate_pct": ppo["win_rate_vs_stupp_pct"],
        "rl_drug_auc_ratio_vs_stupp": ppo["drug_auc_ratio_vs_stupp"],
        "wilcoxon_p_value": ppo["wilcoxon_p_value"],
        "paired_t_test_p_value": ppo["paired_t_p_value_log"],
        "cohens_d": ppo["cohens_d_log"],
        "patient_details": patients,
        "notes": [
            f"Progression (V90 > {PROGRESSION_THRESHOLD_MM3:g} mm^3, the original definition) is reported as counts;"
            " McNemar/kappa on progression were dropped.",
            "Statistics compare log volumes (paired); win = lower day-90 volume than Stupp for the same patient.",
            "Kill rates are patient-specific here; script-66 policies were trained with fixed kill rates.",
        ],
    }
    path = OUTPUT_DIR / "phase8_cohort_metrics.json"
    path.write_text(json.dumps(metrics, indent=2))
    print(f"[Metrics] Saved -> {path}")
    create_visualization(cohort, res, arm_summary, OUTPUT_DIR / "phase8_cohort_analysis.png")

    print("\n" + "=" * 70)
    for arm, s in arm_summary.items():
        print(f"  {ARM_LABELS[arm]:26s} mean {s['mean_final_volume_mm3']:7.2f} +/- {s['std_final_volume_mm3']:6.2f}  "
              f"win {s['win_rate_vs_stupp_pct']:5.1f}%  log-ratio {s['mean_log_ratio_final_vs_stupp']:+.3f} "
              f"CI {['%+.3f' % x for x in s['mean_log_ratio_final_vs_stupp_ci95']]}  "
              f"AUC x{s['drug_auc_ratio_vs_stupp']:.3f}  p_wilcoxon {s['wilcoxon_p_value']}")


if __name__ == "__main__":
    main()
