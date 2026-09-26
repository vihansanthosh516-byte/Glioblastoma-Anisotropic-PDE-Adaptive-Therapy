#!/usr/bin/env python3
"""
Script 60: Baselines & Ablation under an EQUAL DRUG BUDGET
==========================================================
Every arm gets the same total drug AUC (Stupp's) and the same seed tumour.
Arms come from src/rl/equal_budget_arms.py (shared with 62, 64, 65):
  Stupp | budgeted 59 heuristic | RL = PPO policy (script 66) | DAgger oracle (script 66)

Physics comes from this script's own FastPDESolver, which carries the ablation
flags. Scenarios are the 30 LHS scenarios of script 59; alpha_sens scales every
kill rate in this solver (unlike script 59, where it is unused).

Ablations: Full model, No DTI, No Mechanics, Pure Reaction-Diffusion.
Two properties of this solver are reported, not changed:
  - use_mechanics is stored but never enters pde_step, so "No Mechanics" is
    identical to the full model by construction.
  - use_dti=False is not isotropic: D_xx = D_white, D_yy = D_zz = D_gray
    everywhere (a uniform x-aligned tensor).
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

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from rl.equal_budget_arms import (STUPP_BUDGET, evaluate_arms, s66,  # noqa: E402
                                  validate_against_numpy)

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
EVAL_GRID = (64, 64, 64)
T_MAX_DAYS = 90
DT_RL_DAYS = 1.0
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = int(DT_RL_DAYS / DT_PDE_EVAL)

D_WHITE_BASE = 0.013
D_GRAY_BASE = 0.0013
K_CARRY = 1.0
GAMMA_CHEMO = 0.05
GAMMA_RAD = 0.08
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
        grid_size: Tuple[int, int, int] = (64, 64, 64),
        dt_pde: float = DT_PDE_EVAL,
        rho: float = 0.02,
        D_white: float = D_WHITE_BASE,
        alpha_sens: float = 1.0,
        use_dti: bool = True,
        use_mechanics: bool = True,
        mechanics_eta: float = 0.1,
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
        self.mechanics_eta = mechanics_eta if use_mechanics else 0.0

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
            # Isotropic scalar diffusion
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
# Equal-budget evaluation
# --------------------------------------------------------------------------- #
FULL = "Full Model (DTI + Mechanics)"
ABLATIONS = {
    FULL: {"use_dti": True, "use_mechanics": True},
    "No DTI (x-aligned tensor)": {"use_dti": False, "use_mechanics": True},
    "No Mechanics": {"use_dti": True, "use_mechanics": False},
    "Pure Reaction-Diffusion": {"use_dti": False, "use_mechanics": False},
}
ARM_LABELS = {"stupp": "Stupp", "heuristic_budgeted": "Budgeted heuristic (59)",
              "ppo": "RL: PPO (66)", "dagger_oracle": "DAgger oracle (66)"}


def kill_row(alpha_sens: float) -> List[float]:
    """Kill rate per action exactly as FastPDESolver.rl_step computes it."""
    return [0.0, GAMMA_CHEMO * alpha_sens, GAMMA_RAD * alpha_sens, (GAMMA_CHEMO + GAMMA_RAD) * alpha_sens]


def arm_stats(res: Dict[str, list], ref: Dict[str, list]) -> Dict[str, Any]:
    f, fs = np.array(res["final_volume_mm3"]), np.array(ref["final_volume_mm3"])
    b, bs = np.array(res["mean_burden_mm3"]), np.array(ref["mean_burden_mm3"])
    auc, aucs = np.array(res["drug_auc"]), np.array(ref["drug_auc"])
    return {
        "mean_final_volume_mm3": float(f.mean()),
        "median_final_volume_mm3": float(np.median(f)),
        "mean_drug_auc": float(auc.mean()),
        "drug_auc_ratio_vs_stupp": float(auc.mean() / aucs.mean()),
        "win_rate_vs_stupp_pct": float((f < fs).mean() * 100),
        "mean_log_ratio_final_vs_stupp": float(np.log(f / fs).mean()),
        "mean_burden_mm3": float(b.mean()),
        "burden_win_rate_vs_stupp_pct": float((b < bs).mean() * 100),
        "mean_log_ratio_burden_vs_stupp": float(np.log(b / bs).mean()),
        "final_volume_mm3": f.tolist(),
        "drug_auc": auc.tolist(),
    }


def main():
    print("=" * 70)
    print("SCRIPT 60: BASELINES & ABLATION AT EQUAL DRUG BUDGET")
    print("=" * 70)
    scenarios = s66.lhs_scenarios()
    kills = [kill_row(p["alpha_sens"]) for p in scenarios]
    print(f"  {len(scenarios)} LHS scenarios, budget = Stupp drug AUC = {STUPP_BUDGET:g}")

    validation, ablations, actions0 = {}, {}, {}
    for name, flags in ABLATIONS.items():
        makers = [lambda p=p, flags=flags: FastPDESolver(grid_size=EVAL_GRID, dt_pde=DT_PDE_EVAL, rho=p["rho"],
                                                          D_white=p["D_w"], alpha_sens=p["alpha_sens"], **flags)
                  for p in scenarios]
        validation[name] = validate_against_numpy(makers[0](), kills[0])
        assert validation[name] < 1e-9, f"batched solver disagrees with numpy solver ({name}): {validation[name]}"
        t0 = time.time()
        res = evaluate_arms(makers, kills)
        ablations[name] = {arm: arm_stats(r, res["stupp"]) for arm, r in res.items()}
        actions0[name] = {arm: r["action_history"][0] for arm, r in res.items()}
        print(f"\n  [{name}]  ({time.time() - t0:.0f}s, validation err {validation[name]:.1e})")
        for arm, s in ablations[name].items():
            print(f"    {ARM_LABELS[arm]:26s} mean final {s['mean_final_volume_mm3']:8.3f}  "
                  f"win {s['win_rate_vs_stupp_pct']:5.1f}%  log-ratio {s['mean_log_ratio_final_vs_stupp']:+.3f}  "
                  f"AUC {s['mean_drug_auc']:.2f} (x{s['drug_auc_ratio_vs_stupp']:.3f})")

    full = ablations[FULL]

    def impact(name):
        return (ablations[name]["ppo"]["mean_final_volume_mm3"] - full["ppo"]["mean_final_volume_mm3"]) \
            / full["ppo"]["mean_final_volume_mm3"] * 100

    metrics = {
        "framework": "equal drug budget; arms from src/rl/equal_budget_arms.py",
        "budget_drug_auc": STUPP_BUDGET,
        "n_scenarios": len(scenarios),
        "scenario_source": "script 59 LHS (rho, D_w, alpha_sens); alpha_sens scales all kill rates here",
        "rl_column": "ppo",
        "validation_max_rel_err_vs_numpy_solver": validation,
        "baselines": {"model": FULL, "arms": full},
        "ablations": ablations,
        "actions_scenario0": actions0,
        "summary": {
            "best_arm_full_model": min(full, key=lambda a: full[a]["mean_log_ratio_final_vs_stupp"]),
            "rl_win_rate_pct": full["ppo"]["win_rate_vs_stupp_pct"],
            "rl_vs_stupp_improvement_pct": (full["stupp"]["mean_final_volume_mm3"] - full["ppo"]["mean_final_volume_mm3"])
            / full["stupp"]["mean_final_volume_mm3"] * 100,
            "rl_drug_auc_ratio_vs_stupp": full["ppo"]["drug_auc_ratio_vs_stupp"],
            "rl_win_rate_by_ablation_pct": {n: a["ppo"]["win_rate_vs_stupp_pct"] for n, a in ablations.items()},
            "ablation_impact_no_dti_pct": impact("No DTI (x-aligned tensor)"),
            "ablation_impact_no_mechanics_pct": impact("No Mechanics"),
            "ablation_impact_pure_rd_pct": impact("Pure Reaction-Diffusion"),
        },
        "notes": [
            "use_mechanics never enters pde_step: 'No Mechanics' equals the full model by construction.",
            "use_dti=False gives D_xx=D_white, D_yy=D_zz=D_gray everywhere (uniform x-aligned), not isotropic.",
            "The old threshold-adaptive baseline and 40-episode REINFORCE policy were removed; they were unbudgeted.",
        ],
    }
    path = OUTPUT_DIR / "ablation_and_baselines_metrics.json"
    path.write_text(json.dumps(metrics, indent=2))
    print(f"\n[Metrics] Saved -> {path}")

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    ax = axes[0]
    arms = [a for a in full if a != "stupp"]
    vals = [full[a]["mean_log_ratio_final_vs_stupp"] for a in arms]
    bars = ax.barh([ARM_LABELS[a] for a in arms], vals, color=["green" if v < 0 else "firebrick" for v in vals])
    for bar, a in zip(bars, arms):
        ax.text(0.005, bar.get_y() + bar.get_height() / 2,
                f"win {full[a]['win_rate_vs_stupp_pct']:.0f}%   AUC x{full[a]['drug_auc_ratio_vs_stupp']:.2f}",
                va="center", fontsize=9)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("mean log(final volume / Stupp final volume)   (< 0 beats Stupp)")
    ax.set_title(f"Panel 1: Equal-budget arms vs Stupp\n({FULL}, {len(scenarios)} LHS scenarios, AUC = {STUPP_BUDGET:g})")
    ax.grid(alpha=0.3, axis="x")

    ax = axes[1]
    names = list(ABLATIONS)
    x = np.arange(len(names))
    w = 0.2
    for i, a in enumerate(full):
        ax.bar(x + (i - 1.5) * w, [ablations[n][a]["mean_final_volume_mm3"] for n in names], w, label=ARM_LABELS[a])
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=12, ha="right")
    ax.set_ylabel("mean final tumour volume (mm^3)")
    ax.set_title("Panel 2: Ablation (every arm at equal drug budget)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")
    plt.suptitle("Script 60: Baselines & Ablation at Equal Drug Budget", fontsize=15, fontweight="bold")
    plt.tight_layout()
    fig_path = OUTPUT_DIR / "ablation_study_figure.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {fig_path}")

    s = metrics["summary"]
    print("\n" + "=" * 70)
    print(f"RL (PPO) win rate vs Stupp: {s['rl_win_rate_pct']:.1f}%   mean-volume reduction "
          f"{s['rl_vs_stupp_improvement_pct']:.1f}%   AUC ratio {s['rl_drug_auc_ratio_vs_stupp']:.3f}")
    print(f"Ablation impact on PPO mean final volume: no DTI {s['ablation_impact_no_dti_pct']:+.2f}%  "
          f"no mechanics {s['ablation_impact_no_mechanics_pct']:+.2f}%  pure RD {s['ablation_impact_pure_rd_pct']:+.2f}%")


if __name__ == "__main__":
    main()
