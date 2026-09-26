#!/usr/bin/env python3
"""
Script 62: Who benefits from starting treatment early at equal drug?
====================================================================
Reframed after script 66. At an equal drug budget (Stupp's AUC), compare
  early start  = PPO arm from script 66 (combo from day 1 until the budget is spent)
  delayed start = Stupp (nothing until day 20)
Arms come from src/rl/equal_budget_arms.py; physics is this script's FastPDESolver.

benefit_i = log(V_Stupp,i / V_early,i) at day 90  (> 0: early start wins)

Cohort: MU-Glioma growing tumours with fit R^2 >= 0.5
(load_real_params_for_track_bc(min_r2=0.5, require_growing=True)): 64 patients.

Threshold rho*: benefit crosses zero. Reported three ways:
  - model rho*: dense rho grid on the same PDE (the ground truth of this model)
  - cohort rho*: estimated from the 64 patients alone
  - bootstrap 95% CI of the cohort estimate (1000 resamples of patients)
In this model a patient's outcome is a deterministic function of rho (diffusion
is negligible, same seed tumour), so the bootstrap measures how well the
cohort's rho values pin the threshold down, not biological noise.
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
from rl.equal_budget_arms import STUPP_BUDGET, ppo, stupp, evaluate_arms, validate_against_numpy  # noqa: E402
from mu_glioma_loader import load_real_params_for_track_bc  # noqa: E402

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
BOOTSTRAP_SAMPLES = 1000
BOOTSTRAP_SEED = 62
EVAL_GRID = (64, 64, 64)
T_MAX_DAYS = 90
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
DT_RL_DAYS = 1.0
DT_PDE_EVAL = 0.2
N_PDE_SUBSTEPS_EVAL = 5
KILL_ROW = [0.0, GAMMA_CHEMO, GAMMA_RAD, GAMMA_CHEMO + GAMMA_RAD]
EARLY_VS_DELAYED = {"stupp": stupp, "ppo": ppo}

OUTPUT_DIR = PROJECT_ROOT / "output"
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
# Early vs delayed start at equal drug
# --------------------------------------------------------------------------- #
def run_patients(params: List[Dict[str, float]]) -> Dict[str, np.ndarray]:
    makers = [lambda p=p: FastPDESolver(grid_size=EVAL_GRID, dt_pde=DT_PDE_EVAL, rho=p["rho"], D_white=p["D_w"])
              for p in params]
    res = evaluate_arms(makers, [KILL_ROW] * len(params), arms=EARLY_VS_DELAYED)
    s, e = res["stupp"], res["ppo"]
    arr = lambda d, k: np.array(d[k])
    return {"benefit_final": np.log(arr(s, "final_volume_mm3") / arr(e, "final_volume_mm3")),
            "benefit_burden": np.log(arr(s, "mean_burden_mm3") / arr(e, "mean_burden_mm3")),
            "stupp_final": arr(s, "final_volume_mm3"), "early_final": arr(e, "final_volume_mm3"),
            "stupp_auc": arr(s, "drug_auc"), "early_auc": arr(e, "drug_auc")}


def interp_crossing(rho: np.ndarray, benefit: np.ndarray) -> float:
    """Linear interpolation between the last sorted rho with benefit > 0 and the
    first with benefit <= 0. NaN when the sample has no patient on one side."""
    order = np.argsort(rho, kind="stable")
    r, b = rho[order], benefit[order]
    neg = np.flatnonzero(b <= 0)
    if len(neg) == 0 or neg[0] == 0:
        return float("nan")
    j = neg[0]
    i = j - 1
    return float(r[i] + (r[j] - r[i]) * b[i] / (b[i] - b[j]))


def quadratic_crossing(rho: np.ndarray, benefit: np.ndarray, lo: float = 0.0, hi: float = 0.3) -> float:
    """Root of a least-squares quadratic benefit(rho) where the fit crosses downward.
    Can extrapolate past the largest sampled rho."""
    c = np.polyfit(rho, benefit, 2)
    roots = [r.real for r in np.roots(c) if abs(r.imag) < 1e-12 and lo < r.real < hi
             and np.polyval(np.polyder(c), r.real) < 0]
    return float(min(roots)) if roots else float("nan")


def ci(samples: np.ndarray) -> Dict[str, Any]:
    ok = samples[np.isfinite(samples)]
    if len(ok) == 0:
        return {"identifiable_fraction": 0.0, "ci95": None, "median": None}
    return {"identifiable_fraction": float(len(ok) / len(samples)),
            "ci95": [float(np.percentile(ok, 2.5)), float(np.percentile(ok, 97.5))],
            "median": float(np.median(ok))}


def main():
    print("=" * 70)
    print("SCRIPT 62: EARLY vs DELAYED START AT EQUAL DRUG - rho THRESHOLD")
    print("=" * 70)
    rows = load_real_params_for_track_bc(min_r2=0.5, require_growing=True)
    rows = sorted(rows, key=lambda r: r["patient_id"])
    cohort = [{"id": r["patient_id"], "rho": float(r["rho_per_day"]), "D_w": float(r["D_mm2_per_day"])} for r in rows]
    rho = np.array([p["rho"] for p in cohort])
    print(f"  cohort: {len(cohort)} growing MU-Glioma patients (R^2 >= 0.5), rho {rho.min():.4f}-{rho.max():.4f}/day")

    val_err = validate_against_numpy(FastPDESolver(rho=cohort[0]["rho"], D_white=cohort[0]["D_w"]), KILL_ROW)
    assert val_err < 1e-9, f"batched solver disagrees with numpy solver: {val_err}"

    # Model ground truth: dense rho grid (D_w = rho/10, the cohort's convention), then refine
    t0 = time.time()
    grid = np.linspace(0.002, 0.12, 30)
    g = run_patients([{"rho": r, "D_w": r / 10} for r in grid])
    coarse = interp_crossing(grid, g["benefit_final"])
    fine_grid = np.linspace(coarse - 0.004, coarse + 0.004, 17)
    gf = run_patients([{"rho": r, "D_w": r / 10} for r in fine_grid])
    model_rho_star = interp_crossing(fine_grid, gf["benefit_final"])
    curve_rho = np.concatenate([grid, fine_grid])
    order = np.argsort(curve_rho)
    curve_b = np.concatenate([g["benefit_final"], gf["benefit_final"]])[order]
    curve_bb = np.concatenate([g["benefit_burden"], gf["benefit_burden"]])[order]
    curve_rho = curve_rho[order]
    print(f"  model rho* = {model_rho_star:.5f}/day  (grid runs {time.time() - t0:.0f}s)")
    print(f"  burden benefit over grid: min {curve_bb.min():+.3f}  max {curve_bb.max():+.3f}")

    t0 = time.time()
    c = run_patients(cohort)
    b = c["benefit_final"]
    print(f"  cohort runs {time.time() - t0:.0f}s; early start better (final) for {(b > 0).sum()}/{len(b)}")

    cohort_interp = interp_crossing(rho, b)
    cohort_quad = quadratic_crossing(rho, b)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    boot_i, boot_q = np.empty(BOOTSTRAP_SAMPLES), np.empty(BOOTSTRAP_SAMPLES)
    for k in range(BOOTSTRAP_SAMPLES):
        idx = rng.integers(0, len(rho), len(rho))
        boot_i[k] = interp_crossing(rho[idx], b[idx])
        boot_q[k] = quadratic_crossing(rho[idx], b[idx])
    bi, bq = ci(boot_i), ci(boot_q)
    above = [p for p, bf in zip(cohort, b) if p["rho"] > model_rho_star]

    metrics = {
        "question": "At equal drug (Stupp AUC), which patients benefit from starting treatment early?",
        "framework": "equal drug budget; arms from src/rl/equal_budget_arms.py",
        "budget_drug_auc": STUPP_BUDGET,
        "early_start_arm": "ppo (script 66: combo from day 1 until budget spent)",
        "delayed_start_arm": "stupp (first treatment day 20)",
        "benefit_definition": "log(V_stupp / V_early) at day 90; > 0 means early start wins",
        "cohort": {"source": "output/mu_glioma_params_real.csv", "filter": "trajectory == growing and r_squared >= 0.5",
                   "n_patients": len(cohort), "rho_range": [float(rho.min()), float(rho.max())],
                   "note": "No filter of the real cohort yields 61 patients; this is the loader's default quality filter."},
        "validation_max_rel_err_vs_numpy_solver": val_err,
        "model_threshold": {"rho_star_per_day": model_rho_star, "resolution_per_day": float(fine_grid[1] - fine_grid[0]),
                            "benefit_curve": {"rho": curve_rho.tolist(), "benefit_final": curve_b.tolist(),
                                              "benefit_burden": curve_bb.tolist()}},
        "cohort_threshold": {"interpolated_rho_star": cohort_interp, "quadratic_fit_rho_star": cohort_quad},
        "bootstrap": {"n_resamples": BOOTSTRAP_SAMPLES, "seed": BOOTSTRAP_SEED,
                      "interpolated": bi, "quadratic_fit": bq},
        "patients": [{"id": p["id"], "rho": p["rho"], "benefit_final": float(bf), "benefit_burden": float(bb),
                      "early_start_better_final": bool(bf > 0), "stupp_final_mm3": float(sf), "early_final_mm3": float(ef),
                      "stupp_auc": float(sa), "early_auc": float(ea)}
                     for p, bf, bb, sf, ef, sa, ea in zip(cohort, b, c["benefit_burden"], c["stupp_final"],
                                                          c["early_final"], c["stupp_auc"], c["early_auc"])],
        "summary": {
            "rho_star_model_per_day": model_rho_star,
            "rho_star_cohort_per_day": cohort_interp,
            "rho_star_ci95_per_day": bi["ci95"],
            "ci_identifiable_fraction": bi["identifiable_fraction"],
            "n_patients_above_rho_star": len(above),
            "patients_above_rho_star": [p["id"] for p in above],
            "n_early_start_better_final": int((b > 0).sum()),
            "pct_early_start_better_final": float((b > 0).mean() * 100),
            "pct_early_start_better_burden": float((c["benefit_burden"] > 0).mean() * 100),
            "mean_benefit_final": float(b.mean()),
            "mean_benefit_burden": float(c["benefit_burden"].mean()),
            "drug_auc_early_mean": float(c["early_auc"].mean()), "drug_auc_stupp_mean": float(c["stupp_auc"].mean()),
        },
        "notes": [
            "Outcome is a deterministic function of rho in this model (diffusion negligible, identical seed tumour);"
            " the bootstrap CI reflects how densely the cohort samples rho near the threshold.",
            "Above rho*, delaying treatment lowers day-90 volume because logistic saturation slows untreated growth;"
            " early start still lowers mean 90-day burden for every patient.",
            "Real rho values are net growth rates fitted to clinically treated tumours, used here as intrinsic rho.",
        ],
    }
    path = OUTPUT_DIR / "biomarker_stability_metrics.json"
    path.write_text(json.dumps(metrics, indent=2))
    print(f"[Metrics] Saved -> {path}")

    fig, axes = plt.subplots(1, 3, figsize=(20, 6))
    ax = axes[0]
    ax.plot(curve_rho, curve_b, "k-", lw=1.5, label="model (dense rho grid)")
    ax.scatter(rho, b, c=np.where(b > 0, "green", "firebrick"), s=35, edgecolors="black", zorder=3,
               label=f"patients (n={len(rho)})")
    ax.axhline(0, color="gray", lw=0.8)
    ax.axvline(model_rho_star, color="purple", ls="--", label=f"model rho* = {model_rho_star:.4f}")
    if bi["ci95"]:
        ax.axvspan(*bi["ci95"], color="purple", alpha=0.15, label="bootstrap 95% CI")
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("log(V_Stupp / V_early) at day 90")
    ax.set_title("Panel 1: Day-90 benefit of starting early (equal drug)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.plot(curve_rho, curve_bb, "k-", lw=1.5, label="model")
    ax.scatter(rho, c["benefit_burden"], color="green", s=35, edgecolors="black", zorder=3, label="patients")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("log(burden_Stupp / burden_early)")
    ax.set_title("Panel 2: Mean 90-day burden benefit of starting early")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[2]
    for s, lab, col in ((boot_i, "interpolated", "#1f77b4"), (boot_q, "quadratic fit", "#ff7f0e")):
        ok = s[np.isfinite(s)]
        if len(ok):
            ax.hist(ok, bins=40, alpha=0.6, color=col, label=f"{lab} ({len(ok)}/{len(s)} identifiable)")
    ax.axvline(model_rho_star, color="purple", ls="--", label="model rho*")
    ax.set_xlabel("rho* estimate (1/day)")
    ax.set_ylabel("bootstrap count")
    ax.set_title(f"Panel 3: Bootstrap of rho* ({BOOTSTRAP_SAMPLES} resamples)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)
    plt.suptitle("Script 62: Who benefits from early treatment at equal drug?", fontsize=15, fontweight="bold")
    plt.tight_layout()
    fig_path = OUTPUT_DIR / "biomarker_stability.png"
    plt.savefig(fig_path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {fig_path}")

    s = metrics["summary"]
    print("\n" + "=" * 70)
    print(f"model rho* {s['rho_star_model_per_day']:.5f}  cohort rho* {s['rho_star_cohort_per_day']:.5f}  "
          f"95% CI {s['rho_star_ci95_per_day']}  identifiable {s['ci_identifiable_fraction']:.2f}")
    print(f"early start better (final) {s['n_early_start_better_final']}/{len(cohort)}; "
          f"better on burden {s['pct_early_start_better_burden']:.0f}%; above rho*: {s['patients_above_rho_star']}")


if __name__ == "__main__":
    main()
