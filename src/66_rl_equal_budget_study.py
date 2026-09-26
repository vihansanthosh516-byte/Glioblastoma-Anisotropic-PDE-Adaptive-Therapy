#!/usr/bin/env python3
"""
Script 66: Can a properly trained RL agent beat Stupp at EQUAL drug budget?
===========================================================================
Follow-up to script 59. At equal drug AUC the 59 heuristic wins 10% vs Stupp.
This study asks whether that is a policy failure or a limit of the model.

Stages (run all with no flags, or one with --stage):
  validate  Batched torch port of the 59 PDE reproduces the numpy solver.
  theory    Closed-form per-voxel logistic/log-kill prediction vs PDE;
            timing-vs-composition decomposition of the Stupp gap.
  oracle    Black-box schedule search (cross-entropy method) for the best
            equal-budget open-loop schedule -- a ceiling for any policy.
  train     PPO (budget-aware obs, action masking) on the 91-patient
            MU-Glioma growing cohort; reward variants x seeds; BC from oracle.
  evaluate  All policies on the 30 LHS scenarios of script 59, 64^3 grid,
            same seed tumour, same budget (= Stupp AUC) for every arm.

Physics, initial condition, action kill rates, intensity scale and budget
are imported from script 59 unchanged.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

import numpy as np
import torch
import torch.nn.functional as F

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUTPUT_DIR = PROJECT_ROOT / "output" / "rl_equal_budget"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

_spec = spec_from_file_location("s59", PROJECT_ROOT / "src" / "59_sensitivity_analysis.py")
s59 = module_from_spec(_spec)
_spec.loader.exec_module(s59)

torch.set_default_dtype(torch.float64)
N_DAYS = 90
N_SUB = 5
DT_SUB = s59.DT_PDE_EVAL
KILL = torch.tensor([0.0, s59.GAMMA_CHEMO, s59.GAMMA_RAD, s59.GAMMA_CHEMO + s59.GAMMA_RAD])
INTENSITY = torch.tensor([s59.ACTION_INTENSITY[a] for a in range(4)])
N_ACTIONS = 4


def stupp_action(day: int) -> int:
    if 20 <= day < 50:
        return 3
    if 50 <= day <= 90:
        return 1 if (day % 28) < 5 else 0
    return 0


STUPP_SCHEDULE = [stupp_action(d) for d in range(1, N_DAYS + 1)]
BUDGET = float(sum(s59.ACTION_INTENSITY[a] for a in STUPP_SCHEDULE) * s59.DT_RL_DAYS)


# --------------------------------------------------------------------------- #
# Batched port of s59.FastPDESolver (same tensor field, fluxes, seed, Euler)
# --------------------------------------------------------------------------- #
class BatchedPDE:
    diffuse = True

    def __init__(self, grid: int, rho: Sequence[float], D_w: Sequence[float]):
        ref = s59.FastPDESolver(grid_size=(grid,) * 3, dt_pde=DT_SUB, rho=0.01, D_white=1.0)
        self.n = grid
        self.dx = ref.dx
        self.B = len(rho)
        self.rho = torch.tensor(rho).view(-1, 1, 1, 1)
        # D is affine in D_white: D = D_gray*I + (D_w - D_gray)*M, M from the unit-D_w field
        ref_g = s59.FastPDESolver(grid_size=(grid,) * 3, dt_pde=DT_SUB, rho=0.01, D_white=s59.D_GRAY_BASE)
        scale = (torch.tensor(D_w) - s59.D_GRAY_BASE).view(-1, 1, 1, 1) / (1.0 - s59.D_GRAY_BASE)
        self.face = {}
        for name in ["Dxx_xf", "Dxy_xf", "Dxz_xf", "Dyy_yf", "Dxy_yf", "Dyz_yf", "Dzz_zf", "Dxz_zf", "Dyz_zf"]:
            base = torch.from_numpy(getattr(ref_g, name))
            unit = torch.from_numpy(getattr(ref, name))
            self.face[name] = base.unsqueeze(0) + scale * (unit - base).unsqueeze(0)
        xx, yy, zz = np.mgrid[0:grid, 0:grid, 0:grid].astype(float)
        c = grid // 2
        r2 = ((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2) * self.dx ** 2
        seed = s59.SEED_AMPLITUDE * np.exp(-r2 / (2 * s59.SEED_SIGMA_MM ** 2))
        # both protocols in s59 scale the seed by 0.1 after reset
        self.u0 = torch.from_numpy(seed * 0.1).unsqueeze(0).repeat(self.B, 1, 1, 1)
        self.reset()

    def reset(self):
        self.u = self.u0.clone()
        self.v0 = self.volume()
        return self.v0

    def volume(self) -> torch.Tensor:
        return self.u.sum(dim=(1, 2, 3)) * self.dx ** 3

    @staticmethod
    def _pad(t, axis):
        pad = [0, 0, 0, 0, 0, 0]
        pad[2 * (2 - axis)] = pad[2 * (2 - axis) + 1] = 1
        return F.pad(t.unsqueeze(1), pad, mode="replicate").squeeze(1)

    def _face_avg(self, t, axis):
        p = self._pad(t, axis)
        n = p.shape[axis + 1]
        return 0.5 * (p.narrow(axis + 1, 0, n - 1) + p.narrow(axis + 1, 1, n - 1))

    def _div(self, u):
        dx, f = self.dx, self.face
        up = F.pad(u.unsqueeze(1), [1] * 6, mode="replicate").squeeze(1)
        ux_cc = (up[:, 2:, 1:-1, 1:-1] - up[:, :-2, 1:-1, 1:-1]) / (2 * dx)
        uy_cc = (up[:, 1:-1, 2:, 1:-1] - up[:, 1:-1, :-2, 1:-1]) / (2 * dx)
        uz_cc = (up[:, 1:-1, 1:-1, 2:] - up[:, 1:-1, 1:-1, :-2]) / (2 * dx)
        ux_xf = (up[:, 1:, 1:-1, 1:-1] - up[:, :-1, 1:-1, 1:-1]) / dx
        uy_yf = (up[:, 1:-1, 1:, 1:-1] - up[:, 1:-1, :-1, 1:-1]) / dx
        uz_zf = (up[:, 1:-1, 1:-1, 1:] - up[:, 1:-1, 1:-1, :-1]) / dx
        Fx = f["Dxx_xf"] * ux_xf + f["Dxy_xf"] * self._face_avg(uy_cc, 0) + f["Dxz_xf"] * self._face_avg(uz_cc, 0)
        Fy = f["Dyy_yf"] * uy_yf + f["Dxy_yf"] * self._face_avg(ux_cc, 1) + f["Dyz_yf"] * self._face_avg(uz_cc, 1)
        Fz = f["Dzz_zf"] * uz_zf + f["Dxz_zf"] * self._face_avg(ux_cc, 2) + f["Dyz_zf"] * self._face_avg(uy_cc, 2)
        return ((Fx[:, 1:] - Fx[:, :-1]) + (Fy[:, :, 1:] - Fy[:, :, :-1]) + (Fz[:, :, :, 1:] - Fz[:, :, :, :-1])) / dx

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        kill = KILL[actions].view(-1, 1, 1, 1)
        for _ in range(N_SUB):
            u = self.u
            div = self._div(u) if self.diffuse else 0.0
            self.u = torch.clamp(u + DT_SUB * (div + self.rho * u * (1.0 - u / s59.K_CARRY) - kill * u),
                                 0.0, s59.K_CARRY)
        return self.volume()

    def u_max(self) -> torch.Tensor:
        return self.u.amax(dim=(1, 2, 3))


class ReactionSim:
    """The 64^3 model with the diffusion term dropped. Voxels then evolve
    independently, and voxels sharing a seed value evolve identically, so the
    grid collapses exactly onto its unique seed values weighted by counts."""

    def __init__(self, rho: Sequence[float], grid: int = 64):
        dx = 128.0 / grid
        xx, yy, zz = np.mgrid[0:grid, 0:grid, 0:grid].astype(float)
        c = grid // 2
        r2 = ((xx - c) ** 2 + (yy - c) ** 2 + (zz - c) ** 2) * dx ** 2
        seed = (s59.SEED_AMPLITUDE * np.exp(-r2 / (2 * s59.SEED_SIGMA_MM ** 2)) * 0.1).ravel()
        vals, counts = np.unique(seed, return_counts=True)
        self.dx = dx
        self.B = len(rho)
        self.rho = torch.tensor(rho).view(-1, 1)
        self.w = torch.from_numpy(counts.astype(float)) * dx ** 3
        self.u0 = torch.from_numpy(vals).unsqueeze(0).repeat(self.B, 1)
        self.reset()

    def reset(self):
        self.u = self.u0.clone()
        self.v0 = self.volume()
        return self.v0

    def volume(self) -> torch.Tensor:
        return self.u @ self.w

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        kill = KILL[actions].view(-1, 1)
        for _ in range(N_SUB):
            u = self.u
            self.u = torch.clamp(u + DT_SUB * (self.rho * u * (1.0 - u / s59.K_CARRY) - kill * u), 0.0, s59.K_CARRY)
        return self.volume()

    def u_max(self) -> torch.Tensor:
        return self.u.amax(dim=1)


def run_schedules(pde: BatchedPDE, schedules: torch.Tensor) -> torch.Tensor:
    """schedules: (B, 90) int actions. Returns (B, 90) volume after each day."""
    pde.reset()
    vols = []
    for d in range(N_DAYS):
        vols.append(pde.step(schedules[:, d]))
    return torch.stack(vols, dim=1)


def rollout(pde, policy: Callable[[Dict], torch.Tensor], budget: float = BUDGET,
            intensity: torch.Tensor = INTENSITY) -> Dict:
    """Closed-loop rollout with the s59 budget rule: an action that would push
    cumulative AUC past the budget is replaced by a drug holiday."""
    v0 = pde.reset()
    B = pde.B
    auc = torch.zeros(B)
    vol, prev = v0.clone(), v0.clone()
    vols, acts, umax = [], [], []
    for d in range(N_DAYS):
        obs = {"day": d, "vol": vol, "prev": prev, "v0": v0, "auc": auc, "budget": budget,
               "u_max": pde.u_max(), "intensity": intensity}
        a = policy(obs).long()
        a = torch.where(auc + intensity[a] * s59.DT_RL_DAYS > budget + 1e-9, torch.zeros_like(a), a)
        auc = auc + intensity[a] * s59.DT_RL_DAYS
        prev = vol
        vol = pde.step(a)
        vols.append(vol)
        acts.append(a)
        umax.append(pde.u_max())
    V = torch.stack(vols, 1)
    return {"vols": V, "actions": torch.stack(acts, 1), "auc": auc, "final": V[:, -1],
            "mean_burden": V.mean(1), "v0": v0, "u_max": torch.stack(umax, 1)}


def schedule_policy(sched: Sequence[int]) -> Callable:
    s = torch.tensor(list(sched))
    return lambda obs: s[obs["day"]].expand(obs["vol"].shape[0])


def heuristic59_policy(obs: Dict) -> torch.Tensor:
    """Exact batched copy of s59.run_rl_adaptive's rule (policy=None branch + guardrail)."""
    nv = obs["vol"] / obs["v0"]
    a = torch.where(nv > 0.05, 3, torch.where(nv > 0.01, 2, 0))
    return torch.where((a == 0) & (obs["vol"] > 0.05 * obs["v0"]), 3, a)


def block_schedule(action: int, start_day: int, n: int) -> List[int]:
    return [action if start_day <= d < start_day + n else 0 for d in range(1, N_DAYS + 1)]


def lhs_scenarios() -> List[Dict[str, float]]:
    return s59.generate_parameter_samples(s59.N_SCENARIOS, method="lhs")


def save_json(obj, name: str):
    path = OUTPUT_DIR / name
    path.write_text(json.dumps(obj, indent=2))
    print(f"[saved] {path}")


# --------------------------------------------------------------------------- #
# Stage: validate
# --------------------------------------------------------------------------- #
def stage_validate():
    sc = lhs_scenarios()
    idx = [0, 3, 18]  # low / high / highest-rho scenarios
    ref = []
    t0 = time.time()
    for i in idx:
        p = sc[i]
        solver = s59.FastPDESolver(grid_size=s59.EVAL_GRID, dt_pde=s59.DT_PDE_EVAL, rho=p["rho"],
                                   D_white=p["D_w"], alpha_sens=p["alpha_sens"])
        env = s59.GbmTherapyEnv(solver)
        r = s59.run_stupp_protocol(env)
        ref.append([t["volume_mm3"] for t in r["trajectory"]])
    t_np = time.time() - t0
    t0 = time.time()
    pde = BatchedPDE(64, [sc[i]["rho"] for i in idx], [sc[i]["D_w"] for i in idx])
    sched = torch.tensor([STUPP_SCHEDULE] * len(idx))
    vols = run_schedules(pde, sched).numpy()
    t_torch = time.time() - t0
    ref = np.array(ref)
    rel = float(np.max(np.abs(vols - ref) / ref))
    out = {"scenarios": idx, "max_rel_err_trajectory": rel, "numpy_final": ref[:, -1].tolist(),
           "torch_final": vols[:, -1].tolist(), "numpy_seconds": t_np, "torch_seconds_batched": t_torch,
           "budget_stupp_auc": BUDGET}
    print(json.dumps(out, indent=2))
    assert rel < 1e-9, f"torch port disagrees with numpy solver: {rel}"

    # ReactionSim must equal the 64^3 PDE with the diffusion term switched off
    pde.diffuse = False
    nd = run_schedules(pde, sched)
    rs = ReactionSim([sc[i]["rho"] for i in idx])
    rv = run_schedules(rs, sched)
    out["reaction_sim_vs_nodiffusion_pde_max_rel_err"] = float(((rv - nd).abs() / nd).max())
    out["reaction_sim_n_voxel_classes"] = int(rs.u0.shape[1])
    full = torch.from_numpy(vols[:, -1])
    out["diffusion_effect_stupp_final_max_rel"] = float(((nd[:, -1] - full).abs() / full).max())
    print(json.dumps({k: out[k] for k in list(out)[-3:]}, indent=2))
    assert out["reaction_sim_vs_nodiffusion_pde_max_rel_err"] < 1e-9
    save_json(out, "validate.json")


# --------------------------------------------------------------------------- #
# Stage: theory
# --------------------------------------------------------------------------- #
def closed_form_volumes(pde: BatchedPDE, sched: Sequence[int]) -> torch.Tensor:
    """Per-voxel logistic + log-kill (Bernoulli ODE), no diffusion, exact in time:
    u(t) = u0 e^{G(t)} / (1 + u0 (rho/K) J(t)),  G = int (rho - k),  J = int e^{G}.
    Returns (B, 90) daily volumes."""
    k = KILL[torch.tensor(list(sched))]
    rho = pde.rho.view(-1, 1)
    g = rho - k.view(1, -1)                                  # (B, 90)
    G = torch.cumsum(g, 1)
    G_prev = G - g
    J = torch.cumsum(torch.exp(G_prev) * torch.expm1(g) / g, 1)
    u0 = pde.u0.flatten(1)                                   # (B, V)
    vols = [(u0 * torch.exp(G[:, d:d + 1]) / (1 + u0 * rho / s59.K_CARRY * J[:, d:d + 1])).sum(1)
            for d in range(N_DAYS)]
    return torch.stack(vols, 1) * pde.dx ** 3


THEORY_SCHEDULES = {
    "stupp": STUPP_SCHEDULE,
    "stupp_mix_early": [3] * 30 + [1] * 10 + [0] * 50,
    "stupp_mix_late": [0] * 50 + [1] * 10 + [3] * 30,
    "combo35_early": block_schedule(3, 1, 35),
    "combo35_stupp_aligned": block_schedule(3, 20, 35),
    "combo35_late": block_schedule(3, 56, 35),
}


def summarize(res: Dict, ref: Optional[Dict] = None) -> Dict:
    out = {"mean_final_mm3": float(res["final"].mean()), "mean_burden_mm3": float(res["mean_burden"].mean()),
           "mean_auc": float(res["auc"].mean()), "min_auc": float(res["auc"].min()),
           "max_auc": float(res["auc"].max()), "final_mm3": res["final"].tolist()}
    if ref is not None:
        out["win_rate_final_pct"] = float((res["final"] < ref["final"]).double().mean() * 100)
        out["win_rate_burden_pct"] = float((res["mean_burden"] < ref["mean_burden"]).double().mean() * 100)
        out["mean_log_ratio_final_vs_stupp"] = float(torch.log(res["final"] / ref["final"]).mean())
        out["mean_log_ratio_burden_vs_stupp"] = float(torch.log(res["mean_burden"] / ref["mean_burden"]).mean())
    return out


def stage_theory():
    sc = lhs_scenarios()
    pde = BatchedPDE(64, [p["rho"] for p in sc], [p["D_w"] for p in sc])
    results, summary = {}, {}
    for name, sched in THEORY_SCHEDULES.items():
        t0 = time.time()
        results[name] = rollout(pde, schedule_policy(sched))
        print(f"  {name:24s} mean final {results[name]['final'].mean():7.3f}  ({time.time()-t0:.0f}s)")
    t0 = time.time()
    results["heuristic59"] = rollout(pde, heuristic59_policy)
    print(f"  heuristic59              mean final {results['heuristic59']['final'].mean():7.3f}  ({time.time()-t0:.0f}s)")

    ref = results["stupp"]
    for name, res in results.items():
        summary[name] = summarize(res, ref)

    # How well does the diffusion-free closed form explain the PDE?
    fit = {}
    for name, sched in THEORY_SCHEDULES.items():
        V = closed_form_volumes(pde, sched)
        rel = (V - results[name]["vols"]).abs() / results[name]["vols"]
        lin = results[name]["v0"] * torch.exp(pde.rho.view(-1) * N_DAYS - KILL[torch.tensor(sched)].sum())
        fit[name] = {"closed_form_max_rel_err_any_day": float(rel.max()),
                     "closed_form_mean_rel_err_final": float(rel[:, -1].mean()),
                     "pure_exponential_mean_rel_err_final": float(((lin - results[name]["final"]).abs()
                                                                   / results[name]["final"]).mean())}
    pde.diffuse = False
    nd = rollout(pde, schedule_policy(STUPP_SCHEDULE))
    pde.diffuse = True
    fit["stupp_pde_without_diffusion_mean_rel_change_final"] = float(
        ((nd["final"] - ref["final"]).abs() / ref["final"]).mean())

    total_kill = {n: float(KILL[torch.tensor(s)].sum()) for n, s in THEORY_SCHEDULES.items()}
    out = {"budget": BUDGET, "grid": 64, "n_scenarios": len(sc),
           "kill_per_unit_auc": {a: float(KILL[a] / INTENSITY[a]) for a in (1, 2, 3)},
           "total_log_kill": total_kill,
           "composition_ceiling_log_gain_vs_stupp": total_kill["combo35_late"] - total_kill["stupp"],
           "closed_form_fit": fit, "schedules": summary,
           "heuristic59_actions_scenario0": results["heuristic59"]["actions"][0].tolist(),
           "heuristic59_auc_per_scenario": results["heuristic59"]["auc"].tolist()}
    save_json(out, "theory.json")
    torch.save({k: {kk: vv for kk, vv in v.items()} for k, v in results.items()}, OUTPUT_DIR / "theory_rollouts.pt")


# --------------------------------------------------------------------------- #
# Cohorts
# --------------------------------------------------------------------------- #
def real_cohort(n_test: int = 21, seed: int = 0) -> Dict[str, List]:
    """91 growing MU-Glioma patients, fixed train/test split."""
    from mu_glioma_loader import load_real_params_for_track_bc
    rows = load_real_params_for_track_bc(min_r2=-np.inf, require_growing=True)
    rows = sorted(rows, key=lambda r: r["patient_id"])
    perm = np.random.default_rng(seed).permutation(len(rows))
    test = set(perm[:n_test].tolist())

    def d_w(r):  # D is not used by the reaction-only training sim; full-PDE eval needs a value
        d = r.get("D_mm2_per_day")
        return float(d) if d is not None and np.isfinite(d) else s59.D_WHITE_BASE

    split = {"train": [], "test": []}
    for i, r in enumerate(rows):
        split["test" if i in test else "train"].append(
            {"id": r["patient_id"], "rho": float(r["rho_per_day"]), "D_w": d_w(r)})
    return split


def sched_tensor_policy(S: torch.Tensor) -> Callable:
    return lambda obs: S[:, obs["day"]]


# --------------------------------------------------------------------------- #
# Stage: oracle  (privileged: knows rho and the model; no feedback needed)
# --------------------------------------------------------------------------- #
def candidate_schedules(intensity: torch.Tensor, budget: float) -> List[Dict]:
    """Single-drug blocks that spend the budget, at every start day, plus a
    one-day fill with the largest action that fits the leftover budget."""
    cands = []
    for a in (1, 2, 3):
        n = int(math.floor(budget / float(intensity[a]) + 1e-9))
        n = min(n, N_DAYS)
        left = budget - n * float(intensity[a])
        fill = max([b for b in (1, 2, 3) if float(intensity[b]) <= left + 1e-9], default=0,
                   key=lambda b: float(intensity[b]))
        for s in range(1, N_DAYS - n + 2):
            sched = block_schedule(a, s, n)
            if fill and n < N_DAYS:
                j = s - 2 if s >= 2 else s - 1 + n  # 0-based day just before (or after) the block
                sched[j] = fill
            cands.append({"action": a, "start": s, "n": n, "fill": fill, "sched": sched})
    return cands


def structured_oracle(rhos: Sequence[float], intensity=INTENSITY, budget=BUDGET, chunk=3000) -> Dict:
    cands = candidate_schedules(intensity, budget)
    K, P = len(cands), len(rhos)
    S_all = torch.tensor([c["sched"] for c in cands])
    final = torch.empty(P, K)
    burden = torch.empty(P, K)
    pairs = [(p, k) for p in range(P) for k in range(K)]
    for i in range(0, len(pairs), chunk):
        pk = pairs[i:i + chunk]
        sim = ReactionSim([rhos[p] for p, _ in pk])
        res = rollout(sim, sched_tensor_policy(S_all[[k for _, k in pk]]), budget, intensity)
        for j, (p, k) in enumerate(pk):
            final[p, k], burden[p, k] = res["final"][j], res["mean_burden"][j]
    out = {"candidates": [{k: v for k, v in c.items() if k != "sched"} for c in cands]}
    for obj, M in (("final", final), ("burden", burden)):
        best = M.argmin(1)
        out[obj] = {"best_idx": best.tolist(), "best_value": M.min(1).values.tolist(),
                    "best_sched": [cands[b]["sched"] for b in best.tolist()],
                    "best_desc": [{k: v for k, v in cands[b].items() if k != "sched"} for b in best.tolist()]}
    out["_final"], out["_burden"] = final, burden
    return out


def cem_search(rhos: Sequence[float], objective: str, intensity=INTENSITY, budget=BUDGET,
               n=256, iters=40, elite=0.1, alpha=0.7, seed=0) -> Dict:
    """Unstructured cross-entropy search over all 4^90 daily schedules."""
    g = torch.Generator().manual_seed(seed)
    P = len(rhos)
    probs = torch.full((P, N_DAYS, N_ACTIONS), 1.0 / N_ACTIONS)
    best_val = torch.full((P,), float("inf"))
    best_sched = torch.zeros(P, N_DAYS, dtype=torch.long)
    n_el = max(1, int(n * elite))
    history = []
    sim = ReactionSim([r for r in rhos for _ in range(n)])
    for it in range(iters):
        S = torch.multinomial(probs.repeat_interleave(n, 0).view(-1, N_ACTIONS), 1, generator=g).view(P * n, N_DAYS)
        res = rollout(sim, sched_tensor_policy(S), budget, intensity)
        val = (res["final"] if objective == "final" else res["mean_burden"]).view(P, n)
        exe = res["actions"].view(P, n, N_DAYS)
        order = val.argsort(1)[:, :n_el]
        el = torch.gather(exe, 1, order.unsqueeze(2).expand(-1, -1, N_DAYS))
        freq = F.one_hot(el, N_ACTIONS).double().mean(1)
        probs = alpha * freq + (1 - alpha) * probs
        cur = val.min(1)
        better = cur.values < best_val
        best_val = torch.where(better, cur.values, best_val)
        best_sched[better] = exe[better, cur.indices[better]]
        history.append(best_val.tolist())
    return {"best_value": best_val.tolist(), "best_sched": best_sched.tolist(), "history": history}


def stage_oracle():
    sc = lhs_scenarios()
    rhos = [p["rho"] for p in sc]
    t0 = time.time()
    so = structured_oracle(rhos)
    print(f"  structured oracle on 30 LHS: {time.time()-t0:.0f}s")
    stupp = rollout(ReactionSim(rhos), schedule_policy(STUPP_SCHEDULE))
    combo_idx = [i for i, c in enumerate(so["candidates"]) if c["action"] == 3]
    out = {"budget": BUDGET, "sim": "reaction-only 64^3 (exact voxel classes)",
           "structured": {obj: {k: v for k, v in so[obj].items()} for obj in ("final", "burden")},
           "candidates": so["candidates"],
           "combo_start_curve_final_log_ratio": (torch.log(so["_final"][:, combo_idx]
                                                           / stupp["final"].view(-1, 1))).tolist(),
           "combo_start_curve_burden_log_ratio": (torch.log(so["_burden"][:, combo_idx]
                                                            / stupp["mean_burden"].view(-1, 1))).tolist(),
           "stupp_final": stupp["final"].tolist(), "stupp_burden": stupp["mean_burden"].tolist()}

    # Unstructured CEM check on six scenarios spanning rho: can it beat the structured optimum?
    pick = sorted(range(len(sc)), key=lambda i: rhos[i])[::6][:5] + [max(range(len(sc)), key=lambda i: rhos[i])]
    out["cem"] = {"scenarios": pick}
    for obj in ("final", "burden"):
        t0 = time.time()
        c = cem_search([rhos[i] for i in pick], obj)
        s_best = [so[obj]["best_value"][i] for i in pick]
        out["cem"][obj] = {"best_value": c["best_value"], "structured_best_value": s_best,
                           "cem_minus_structured_rel": [(a - b) / b for a, b in zip(c["best_value"], s_best)],
                           "best_sched": c["best_sched"], "history": c["history"]}
        print(f"  CEM {obj}: {time.time()-t0:.0f}s  rel diff vs structured "
              f"{['%+.4f' % x for x in out['cem'][obj]['cem_minus_structured_rel']]}")

    # Sensitivity of the ceiling to the intensity scale: make every action equally
    # kill-efficient per unit AUC (intensity proportional to kill). Same scale for both arms.
    prop = KILL / KILL[3]
    b_prop = float(sum(float(prop[a]) for a in STUPP_SCHEDULE))
    sp = structured_oracle(rhos, prop, b_prop)
    stupp_p = rollout(ReactionSim(rhos), schedule_policy(STUPP_SCHEDULE), b_prop, prop)
    out["kill_proportional_intensity"] = {
        "intensity": prop.tolist(), "budget": b_prop,
        "best_final_log_ratio_vs_stupp": (torch.log(torch.tensor(sp["final"]["best_value"]) / stupp_p["final"])).tolist(),
        "best_burden_log_ratio_vs_stupp": (torch.log(torch.tensor(sp["burden"]["best_value"]) / stupp_p["mean_burden"])).tolist(),
        "best_desc_final": sp["final"]["best_desc"], "best_desc_burden": sp["burden"]["best_desc"]}

    # Oracle for the real cohort (BC teacher + real-patient ceiling)
    rc = real_cohort()
    for part in ("train", "test"):
        r = structured_oracle([p["rho"] for p in rc[part]])
        out[f"real_{part}"] = {"ids": [p["id"] for p in rc[part]],
                               "final_best_sched": r["final"]["best_sched"],
                               "final_best_desc": r["final"]["best_desc"],
                               "burden_best_sched": r["burden"]["best_sched"],
                               "burden_best_desc": r["burden"]["best_desc"]}
    save_json(out, "oracle.json")


# --------------------------------------------------------------------------- #
# Stage: train  (PPO with budget-aware observations and action masking; BC)
# --------------------------------------------------------------------------- #
class ActorCritic(torch.nn.Module):
    def __init__(self, obs_dim: int = 5, hidden: int = 64):
        super().__init__()
        nn = torch.nn
        self.pi = nn.Sequential(nn.Linear(obs_dim, hidden), nn.Tanh(), nn.Linear(hidden, hidden), nn.Tanh(),
                                nn.Linear(hidden, N_ACTIONS))
        self.v = nn.Sequential(nn.Linear(obs_dim, hidden), nn.Tanh(), nn.Linear(hidden, hidden), nn.Tanh(),
                               nn.Linear(hidden, 1))


def features(obs: Dict) -> torch.Tensor:
    v = obs["vol"].clamp_min(1e-12)
    p = obs["prev"].clamp_min(1e-12)
    day = torch.full_like(v, obs["day"] / N_DAYS)
    return torch.stack([torch.log(v / obs["v0"]) / 5, (torch.log(v) - torch.log(p)) * 10, day,
                        (obs["budget"] - obs["auc"]) / obs["budget"], obs["u_max"] * 10], 1)


def action_mask(obs: Dict) -> torch.Tensor:
    return obs["auc"].unsqueeze(1) + obs["intensity"].unsqueeze(0) * s59.DT_RL_DAYS <= obs["budget"] + 1e-9


def greedy_policy(net: ActorCritic) -> Callable:
    def pol(obs):
        with torch.no_grad():
            return net.pi(features(obs)).masked_fill(~action_mask(obs), -1e9).argmax(1)
    return pol


def compute_rewards(cfg: Dict, res: Dict, ref: Dict) -> torch.Tensor:
    B = res["vols"].shape[0]
    R = torch.zeros(B, N_DAYS)
    if cfg["reward"] == "final_rel":
        R[:, -1] = -torch.log(res["final"] / ref["final"])
    elif cfg["reward"] == "burden_rel":
        R[:, -1] = -torch.log(res["mean_burden"] / ref["mean_burden"])
    elif cfg["reward"] == "legacy59":  # GbmTherapyEnv.step reward from script 59
        v0 = res["v0"].view(-1, 1)
        nv = res["vols"] / v0
        prev = torch.cat([v0, res["vols"][:, :-1]], 1)
        dv = prev - res["vols"]
        R = -10 * nv - 5 * res["u_max"] - 0.05 * (res["actions"] > 0).double() + 20 * dv.clamp_min(0) / v0
        R[:, -1] += 30.0 * (nv[:, -1] < 0.01).double()
    else:
        raise ValueError(cfg["reward"])
    return R - cfg.get("drug_penalty", 0.0) * INTENSITY[res["actions"]]


def eval_reaction(pol: Callable, rhos: Sequence[float], ref: Dict) -> Dict:
    res = rollout(ReactionSim(rhos), pol)
    lf = torch.log(res["final"] / ref["final"])
    lb = torch.log(res["mean_burden"] / ref["mean_burden"])
    return {"final_log_ratio": float(lf.mean()), "burden_log_ratio": float(lb.mean()),
            "win_final_pct": float((lf < 0).double().mean() * 100), "mean_auc": float(res["auc"].mean())}


def train_ppo(cfg: Dict, train_rho: np.ndarray, val: Dict, seed: int, iters: int, batch: int,
              eval_every: int = 10) -> (ActorCritic, Dict):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    net = ActorCritic()
    opt_pi = torch.optim.Adam(net.pi.parameters(), lr=3e-4)
    opt_v = torch.optim.Adam(net.v.parameters(), lr=1e-3)
    curve, stats = [], []
    for it in range(iters + 1):
        if it % eval_every == 0 or it == iters:
            row = {"iter": it, "episodes": it * batch}
            for name, (rhos, ref) in val.items():
                row[name] = eval_reaction(greedy_policy(net), rhos, ref)
            curve.append(row)
        if it == iters:
            break
        sim = ReactionSim(train_rho[rng.integers(len(train_rho), size=batch)].tolist())
        ref = rollout(sim, schedule_policy(STUPP_SCHEDULE))
        buf = {"f": [], "m": [], "a": [], "logp": [], "v": []}

        def pol(obs):
            f, m = features(obs), action_mask(obs)
            with torch.no_grad():
                dist = torch.distributions.Categorical(logits=net.pi(f).masked_fill(~m, -1e9))
                a = dist.sample()
                buf["f"].append(f); buf["m"].append(m); buf["a"].append(a)
                buf["logp"].append(dist.log_prob(a)); buf["v"].append(net.v(f).squeeze(1))
            return a

        res = rollout(sim, pol)
        assert torch.equal(res["actions"], torch.stack(buf["a"], 1)), "mask/budget rule mismatch"
        R = compute_rewards(cfg, res, ref)
        Vt = torch.stack(buf["v"], 1)
        adv = torch.zeros_like(R)
        last = torch.zeros(R.shape[0])
        for t in reversed(range(N_DAYS)):
            nxt = Vt[:, t + 1] if t + 1 < N_DAYS else torch.zeros_like(last)
            last = R[:, t] + nxt - Vt[:, t] + 0.95 * last
            adv[:, t] = last
        ret = (adv + Vt).flatten()
        f = torch.stack(buf["f"], 1).flatten(0, 1)
        m = torch.stack(buf["m"], 1).flatten(0, 1)
        a = torch.stack(buf["a"], 1).flatten()
        old = torch.stack(buf["logp"], 1).flatten()
        advf = adv.flatten()
        advf = (advf - advf.mean()) / (advf.std() + 1e-8)
        kls, clips, ents = [], [], []
        for _ in range(4):
            for mb in torch.randperm(len(a)).split(1024):
                dist = torch.distributions.Categorical(logits=net.pi(f[mb]).masked_fill(~m[mb], -1e9))
                lp = dist.log_prob(a[mb])
                ratio = torch.exp(lp - old[mb])
                pl = -torch.min(ratio * advf[mb], ratio.clamp(0.8, 1.2) * advf[mb]).mean()
                ent = dist.entropy().mean()
                opt_pi.zero_grad()
                (pl - 0.01 * ent).backward()
                torch.nn.utils.clip_grad_norm_(net.pi.parameters(), 0.5)
                opt_pi.step()
                vl = F.mse_loss(net.v(f[mb]).squeeze(1), ret[mb])
                opt_v.zero_grad()
                vl.backward()
                torch.nn.utils.clip_grad_norm_(net.v.parameters(), 0.5)
                opt_v.step()
                with torch.no_grad():
                    kls.append(float((old[mb] - lp).mean()))
                    clips.append(float(((ratio - 1).abs() > 0.2).double().mean()))
                    ents.append(float(ent))
        stats.append({"iter": it, "train_final_log_ratio": float(torch.log(res["final"] / ref["final"]).mean()),
                      "train_return": float(R.sum(1).mean()), "entropy": float(np.mean(ents)),
                      "approx_kl": float(np.mean(kls)), "clipfrac": float(np.mean(clips)),
                      "mean_auc": float(res["auc"].mean())})
    return net, {"curve": curve, "stats": stats}


def train_bc(teacher_sched: List[List[int]], train_rho: List[float], seed: int, epochs: int = 400) -> ActorCritic:
    torch.manual_seed(seed)
    net = ActorCritic()
    S = torch.tensor(teacher_sched)
    buf = {"f": [], "m": []}

    def pol(obs):
        buf["f"].append(features(obs)); buf["m"].append(action_mask(obs))
        return S[:, obs["day"]]

    res = rollout(ReactionSim(train_rho), pol)
    f = torch.stack(buf["f"], 1).flatten(0, 1)
    m = torch.stack(buf["m"], 1).flatten(0, 1)
    y = res["actions"].flatten()
    opt = torch.optim.Adam(net.pi.parameters(), lr=3e-3)
    for _ in range(epochs):
        loss = F.cross_entropy(net.pi(f).masked_fill(~m, -1e9), y)
        opt.zero_grad(); loss.backward(); opt.step()
    net.bc_train_acc = float((net.pi(f).masked_fill(~m, -1e9).argmax(1) == y).double().mean())
    return net


PPO_RUNS = [
    ("ppo_final_s0", {"reward": "final_rel"}, 0),
    ("ppo_final_s1", {"reward": "final_rel"}, 1),
    ("ppo_final_s2", {"reward": "final_rel"}, 2),
    ("ppo_burden_s0", {"reward": "burden_rel"}, 0),
    ("ppo_final_drug0.05_s0", {"reward": "final_rel", "drug_penalty": 0.05}, 0),
    ("ppo_final_drug0.20_s0", {"reward": "final_rel", "drug_penalty": 0.20}, 0),
    ("ppo_legacy59_s0", {"reward": "legacy59"}, 0),
]


def stage_train(iters: int = 150, batch: int = 128, runs: Optional[List[str]] = None):
    rc = real_cohort()
    train_rho = np.array([p["rho"] for p in rc["train"]])
    lhs_rho = [p["rho"] for p in lhs_scenarios()]
    test_rho = [p["rho"] for p in rc["test"]]
    val = {}
    for name, rhos in (("real_test", test_rho), ("lhs30", lhs_rho)):
        val[name] = (rhos, rollout(ReactionSim(rhos), schedule_policy(STUPP_SCHEDULE)))
    out_path = OUTPUT_DIR / "train.json"
    out = json.loads(out_path.read_text()) if out_path.exists() else {}
    out.update({"iters": iters, "batch": batch, "n_train_patients": len(train_rho),
                "train_rho_range": [float(train_rho.min()), float(train_rho.max())],
                "train_sim": "reaction-only 64^3 (exact voxel classes); eval uses full PDE"})
    out.setdefault("runs", {})
    for name, cfg, seed in PPO_RUNS:
        if runs and name not in runs:
            continue
        t0 = time.time()
        net, log = train_ppo(cfg, train_rho, val, seed, iters, batch)
        torch.save(net.state_dict(), OUTPUT_DIR / f"{name}.pt")
        log.update({"cfg": cfg, "seed": seed, "seconds": time.time() - t0})
        out["runs"][name] = log
        last = log["curve"][-1]
        print(f"  {name:24s} {time.time()-t0:5.0f}s  lhs30 final log-ratio {last['lhs30']['final_log_ratio']:+.4f} "
              f"win {last['lhs30']['win_final_pct']:5.1f}%  real_test {last['real_test']['final_log_ratio']:+.4f}")
        out_path.write_text(json.dumps(out, indent=2))

    if not runs or "bc_oracle" in runs:
        oracle = json.loads((OUTPUT_DIR / "oracle.json").read_text())
        net = train_bc(oracle["real_train"]["final_best_sched"], train_rho.tolist(), seed=0)
        torch.save(net.state_dict(), OUTPUT_DIR / "bc_oracle.pt")
        out["runs"]["bc_oracle"] = {"train_acc": net.bc_train_acc,
                                    "curve": [{"iter": 0, **{k: eval_reaction(greedy_policy(net), *v)
                                                             for k, v in val.items()}}]}
        print(f"  bc_oracle train acc {net.bc_train_acc:.3f}  {out['runs']['bc_oracle']['curve'][0]}")
    out_path.write_text(json.dumps(out, indent=2))
    print(f"[saved] {out_path}")


# --------------------------------------------------------------------------- #
# Stage: evaluate  (full 64^3 PDE with diffusion, equal budget for every arm)
# --------------------------------------------------------------------------- #
def load_net(name: str) -> ActorCritic:
    net = ActorCritic()
    net.load_state_dict(torch.load(OUTPUT_DIR / f"{name}.pt"))
    return net


def stage_evaluate():
    oracle = json.loads((OUTPUT_DIR / "oracle.json").read_text())
    train = json.loads((OUTPUT_DIR / "train.json").read_text())
    sc = lhs_scenarios()
    rc = real_cohort()
    sets = {"lhs30": (sc, oracle["structured"]), "real_test": (rc["test"], None)}
    out = {"budget": BUDGET, "grid": 64, "physics": "full s59 PDE incl. diffusion", "sets": {}}
    for set_name, (params, so) in sets.items():
        pde = BatchedPDE(64, [p["rho"] for p in params], [p["D_w"] for p in params])
        arms = {"stupp": schedule_policy(STUPP_SCHEDULE), "heuristic59": heuristic59_policy,
                "combo35_early": schedule_policy(block_schedule(3, 1, 35)),
                "combo35_late": schedule_policy(block_schedule(3, 56, 35))}
        if so is None:
            so = {obj: {"best_sched": oracle["real_test"][f"{obj}_best_sched"]} for obj in ("final", "burden")}
        arms["oracle_final"] = sched_tensor_policy(torch.tensor(so["final"]["best_sched"]))
        arms["oracle_burden"] = sched_tensor_policy(torch.tensor(so["burden"]["best_sched"]))
        for name in train["runs"]:
            arms[name] = greedy_policy(load_net(name))
        res = {}
        for name, pol in arms.items():
            t0 = time.time()
            res[name] = rollout(pde, pol)
            print(f"  [{set_name}] {name:24s} mean final {res[name]['final'].mean():8.3f}  "
                  f"auc {res[name]['auc'].mean():5.2f}  ({time.time()-t0:.0f}s)")
        ref = res["stupp"]
        block = {}
        for name, r in res.items():
            s = summarize(r, ref)
            s["drug_auc_ratio_vs_stupp"] = float(r["auc"].mean() / ref["auc"].mean())
            s["actions_scenario0_days1_28"] = r["actions"][0, :28].tolist()
            s["actions_scenario0_all"] = r["actions"][0].tolist()
            block[name] = s
        out["sets"][set_name] = block
        torch.save({k: {"vols": v["vols"], "actions": v["actions"]} for k, v in res.items()},
                   OUTPUT_DIR / f"eval_rollouts_{set_name}.pt")
    save_json(out, "evaluate.json")


# --------------------------------------------------------------------------- #
# Stage: report
# --------------------------------------------------------------------------- #
def stage_report():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    theory = json.loads((OUTPUT_DIR / "theory.json").read_text())
    oracle = json.loads((OUTPUT_DIR / "oracle.json").read_text())
    train = json.loads((OUTPUT_DIR / "train.json").read_text())
    ev = json.loads((OUTPUT_DIR / "evaluate.json").read_text())
    lhs = ev["sets"]["lhs30"]

    fig, ax = plt.subplots(2, 2, figsize=(16, 12))
    # A: when to spend a fixed budget
    a = ax[0, 0]
    starts = [c["start"] for c in oracle["candidates"] if c["action"] == 3]
    for key, col, lab in (("combo_start_curve_final_log_ratio", "#1f77b4", "final volume"),
                          ("combo_start_curve_burden_log_ratio", "#d62728", "mean 90-day burden")):
        M = np.array(oracle[key])
        a.plot(starts, M.mean(0), color=col, lw=2, label=f"{lab} (mean of 30)")
        a.fill_between(starts, M.min(0), M.max(0), color=col, alpha=0.2)
    a.axhline(0, color="k", lw=0.8)
    a.set_xlabel("start day of the 35-day combo block (budget = Stupp AUC)")
    a.set_ylabel("log(policy / Stupp)   (< 0 is better)")
    a.set_title("A. Same drug, different timing: the objective decides the schedule")
    a.legend()
    a.grid(alpha=0.3)

    # B: learning curves
    b = ax[0, 1]
    ceiling = lhs["oracle_final"]["mean_log_ratio_final_vs_stupp"]
    for name, run in train["runs"].items():
        if name.startswith("ppo"):
            xs = [r["episodes"] for r in run["curve"]]
            ys = [r["lhs30"]["final_log_ratio"] for r in run["curve"]]
            b.plot(xs, ys, lw=2 if "final_s" in name else 1, label=name)
    b.axhline(0, color="k", lw=0.8, label="Stupp")
    b.axhline(ceiling, color="green", ls="--", label="oracle ceiling (full PDE)")
    b.axhline(lhs["heuristic59"]["mean_log_ratio_final_vs_stupp"], color="gray", ls=":", label="59 heuristic")
    b.set_xlabel("training episodes (real-cohort patients)")
    b.set_ylabel("held-out LHS-30 mean log(final / Stupp final)")
    b.set_title("B. PPO learning curves (greedy policy, reaction-only sim)")
    b.set_ylim(-0.5, 1.0)
    b.legend(fontsize=7)
    b.grid(alpha=0.3)

    # C: equal-budget results on the full PDE
    c = ax[1, 0]
    names = [n for n in lhs if n != "stupp"]
    vals = [lhs[n]["mean_log_ratio_final_vs_stupp"] for n in names]
    bars = c.barh(names, vals, color=["green" if v < 0 else "firebrick" for v in vals])
    for bar, n in zip(bars, names):
        c.text(bar.get_width(), bar.get_y() + bar.get_height() / 2,
               f"  win {lhs[n]['win_rate_final_pct']:.0f}%  AUC x{lhs[n]['drug_auc_ratio_vs_stupp']:.2f}",
               va="center", fontsize=8)
    c.axvline(0, color="k", lw=0.8)
    c.axvline(theory["composition_ceiling_log_gain_vs_stupp"] * -1, color="purple", ls="--",
              label="composition-only bound (total log-kill)")
    c.set_xlabel("mean log(final / Stupp final), LHS-30, full 64^3 PDE")
    c.set_title("C. Equal drug budget: every arm vs Stupp")
    c.legend(fontsize=8)
    c.grid(alpha=0.3, axis="x")

    # D: action rasters, scenario 0
    d = ax[1, 1]
    show = [n for n in ("stupp", "heuristic59", "ppo_final_s0", "ppo_burden_s0", "oracle_final",
                        "oracle_burden", "bc_oracle") if n in lhs]
    R = np.array([lhs[n]["actions_scenario0_all"] for n in show])
    cmap = matplotlib.colors.ListedColormap(["#f0f0f0", "#fdae61", "#abd9e9", "#2c7bb6"])
    d.imshow(R, aspect="auto", cmap=cmap, vmin=-0.5, vmax=3.5, interpolation="nearest",
             extent=[0.5, N_DAYS + 0.5, len(show) - 0.5, -0.5])
    d.set_yticks(range(len(show)))
    d.set_yticklabels(show)
    d.set_xlabel("day")
    d.set_title("D. Scenario 0 schedules (grey off, orange chemo, light-blue RT, blue combo)")

    plt.suptitle(f"Script 66: RL vs Stupp at equal drug budget (AUC = {BUDGET:g})", fontsize=15, fontweight="bold")
    plt.tight_layout()
    path = OUTPUT_DIR / "rl_equal_budget_study.png"
    plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"[saved] {path}")


STAGES: Dict[str, Callable] = {"validate": stage_validate, "theory": stage_theory, "oracle": stage_oracle,
                               "train": stage_train, "evaluate": stage_evaluate, "report": stage_report}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=list(STAGES), default=None)
    ap.add_argument("--iters", type=int, default=150)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--runs", nargs="*", default=None, help="subset of PPO_RUNS names (and/or bc_oracle)")
    args = ap.parse_args()
    torch.set_num_threads(8)
    for name, fn in STAGES.items():
        if args.stage is None or args.stage == name:
            print(f"\n===== stage: {name} =====")
            if name == "train":
                fn(iters=args.iters, batch=args.batch, runs=args.runs)
            else:
                fn()


if __name__ == "__main__":
    main()
