#!/usr/bin/env python3
"""
Script 76: Does adaptive therapy gain anything once resistance is in the model?
===============================================================================
Track C negative 3: the script 59-68 model is one homogeneous population with
log-kill, so there is no clonal competition for adaptive timing to exploit
(Gatenby 2009; Leder 2014). Here the same model is split into drug-sensitive
and drug-resistant cells competing for the same carrying capacity
(Lotka-Volterra, as in Strobl et al. 2021 and Gallagher et al. 2024):

    du_s/dt = rho          u_s (1 - (u_s+u_r)/K) - k_a       u_s
    du_r/dt = rho (1 - c)  u_r (1 - (u_s+u_r)/K) - k_a * eps u_r

k_a is each evaluation set's own kill table (scripts 64/60/66/67), rho its own
growth rate, the seed is the script-59 seed split as (1-f_r0, f_r0), and the
drug budget is the Stupp AUC with the script-59 budget rule. With f_r0 = 0 the
sim reduces exactly to s67.KillReactionSim (asserted).

Literature parameters (scenario grid, every cell reported):
  f_r0  initial resistant fraction  0.001, 0.01, 0.1  (Strobl 2021 uses 1%)
  c     fitness cost of resistance  0, 0.25           (script 44 uses 25%)
  eps   residual drug effect on R   0, 0.2            (0 = fully resistant)
  occupancy  script-59 seed (peak u = 0.08 K) or near capacity (peak 0.8 K);
             Strobl 2021 predicts adaptive gains need a cost OR a tumour near
             its carrying capacity, so this is the axis that should matter
  window     treatment allowed on days 1-90 (script 68) or days 1-365

Arms (all capped at the Stupp budget):
  stupp        script-66 Stupp schedule (days 1-90)
  mtd_front    combo every day from day 1 until the budget is spent
  AT50 / AT80  Zhang et al. 2017 adaptive rule: combo while V >= V_ref, stop
               once V <= 0.5 V_ref (AT50) or 0.8 V_ref (AT80), restart when
               V >= V_ref again; V_ref = baseline volume
  best_fixed   per-patient best of the script-67 fixed schedules (privileged,
               TTP-optimal over 90-day schedules)
  paced        non-adaptive control for the window: the 35 combo days spread
               evenly over the treatment window (days 1-90 or 1-365)

A first run showed AT80 with a 365-day window beating Stupp by ~80 days even
with f_r0 = 0.001 and no fitness cost, i.e. the gain came from spreading the
same budget over a year, not from adaptivity. So the adaptive advantage is
measured against non-adaptive arms with the SAME window:
  primary   AT vs paced (same window), primary cell (f_r0 0.01, cost 0.25,
            eps 0, script-59 seed, 365-day window). The window was set to 365
            after the confounded first run, because adaptive therapy is a
            continuous on/off strategy (Zhang 2017); the 90-day cells are
            reported in full alongside it.
  also      AT vs best_nonadaptive = per-patient max TTP over
            {stupp, mtd_front, best_fixed, paced} (privileged, stricter)
  control   f_r0 = 0 cells (single population, cost/eps irrelevant). If the
            AT - paced difference is the same with and without resistant
            cells, the gain is not clonal competition.

Endpoint: RANO 2.0 TTP (>= 40% above nadir) to day 365, as script 68.
Question: in which cells does an adaptive rule beat the same-window
non-adaptive arms, and by more than in the f_r0 = 0 control?
Run: shards compute cells into output/resistance_adaptive/cells/ (resumable),
then --merge builds the summary:
    python src/76_resistance_adaptive_equal_budget.py --shard 0 --nshards 4   (x4)
    python src/76_resistance_adaptive_equal_budget.py --merge
Output: output/resistance_adaptive/evaluate.json
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable, Dict

import numpy as np
import torch
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUTPUT_DIR = PROJECT_ROOT / "output" / "resistance_adaptive"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

_spec = spec_from_file_location("s68", PROJECT_ROOT / "src" / "68_time_to_progression_equal_budget.py")
s68 = module_from_spec(_spec)
_spec.loader.exec_module(s68)
s67, s66 = s68.s67, s68.s66
s59 = s66.s59

BUDGET, STUPP_SCHEDULE, INTENSITY = s66.BUDGET, s66.STUPP_SCHEDULE, s66.INTENSITY
N_SUB, DT_SUB, K = s66.N_SUB, s66.DT_SUB, s59.K_CARRY
HORIZON = s68.HORIZON

GRID = {"f_r0": (0.001, 0.01, 0.1), "cost": (0.0, 0.25), "eps": (0.0, 0.2),
        "occupancy": ("seed59", "near_capacity"), "window": (90, 365)}
PRIMARY = {"f_r0": 0.01, "cost": 0.25, "eps": 0.0, "occupancy": "seed59", "window": 365}


class TwoPopSim(s66.ReactionSim):
    """s67.KillReactionSim split into sensitive / resistant voxel-class densities."""

    def __init__(self, rho, kill, f_r0: float, cost: float, eps: float, occupancy: str):
        super().__init__(list(rho))
        self.kill = torch.as_tensor(kill, dtype=torch.float64).view(self.B, 4)
        if occupancy == "near_capacity":
            self.u0 = self.u0 * 10.0          # undo the script-59 x0.1 seed scaling: peak 0.8 K
        self.f_r0, self.cost, self.eps = f_r0, cost, eps
        self.reset()

    def reset(self):
        self.us = self.u0 * (1.0 - getattr(self, "f_r0", 0.0))
        self.ur = self.u0 * getattr(self, "f_r0", 0.0)
        self.u = self.us + self.ur
        self.v0 = self.volume()
        return self.v0

    def step(self, actions: torch.Tensor) -> torch.Tensor:
        k = self.kill.gather(1, actions.view(-1, 1))
        for _ in range(N_SUB):
            tot = self.us + self.ur
            gs = self.rho * self.us * (1.0 - tot / K) - k * self.us
            gr = self.rho * (1.0 - self.cost) * self.ur * (1.0 - tot / K) - k * self.eps * self.ur
            self.us = torch.clamp(self.us + DT_SUB * gs, 0.0, K)
            self.ur = torch.clamp(self.ur + DT_SUB * gr, 0.0, K)
        self.u = self.us + self.ur
        return self.volume()

    def resistant_fraction(self) -> torch.Tensor:
        return (self.ur @ self.w) / (self.u @ self.w).clamp_min(1e-12)


def rollout(sim, policy: Callable, window: int) -> Dict[str, torch.Tensor]:
    """Daily loop to HORIZON. Policy acts on days < window with the script-59
    budget rule (an action that would exceed the budget becomes a holiday)."""
    v0 = sim.reset()
    B = sim.B
    auc = torch.zeros(B)
    vol, prev = v0.clone(), v0.clone()
    vols, acts = [v0], []
    for d in range(HORIZON):
        if d < window:
            obs = {"day": d, "vol": vol, "prev": prev, "v0": v0, "auc": auc, "budget": BUDGET,
                   "u_max": sim.u_max(), "intensity": INTENSITY}
            a = policy(obs).long()
            a = torch.where(auc + INTENSITY[a] * s59.DT_RL_DAYS > BUDGET + 1e-9, torch.zeros_like(a), a)
        else:
            a = torch.zeros(B, dtype=torch.long)
        auc = auc + INTENSITY[a] * s59.DT_RL_DAYS
        prev, vol = vol, sim.step(a)
        vols.append(vol)
        acts.append(a)
    V = torch.stack(vols, 1)
    return {**s68.ttp_endpoints(V), "auc": auc, "actions": torch.stack(acts, 1),
            "resistant_fraction_end": sim.resistant_fraction()}


def at_policy(stop_frac: float) -> Callable:
    """Zhang 2017 adaptive rule, batched and stateful."""
    state = {}

    def pol(obs):
        if obs["day"] == 0:
            state["on"] = torch.ones_like(obs["vol"], dtype=torch.bool)
        on = state["on"]
        on = torch.where(obs["vol"] <= stop_frac * obs["v0"], torch.zeros_like(on), on)
        on = torch.where(obs["vol"] >= obs["v0"], torch.ones_like(on), on)
        state["on"] = on
        return torch.where(on, 3, 0)
    return pol


def paced_policy(window: int) -> Callable:
    """35 combo days (the budget) spread evenly over the window."""
    n = int(BUDGET / float(INTENSITY[3]) + 1e-9)

    def pol(obs):
        d = obs["day"]
        on = (d + 1) * n // window > d * n // window
        return torch.full((obs["vol"].shape[0],), 3 if on else 0, dtype=torch.long)
    return pol


def schedule(S: torch.Tensor) -> Callable:
    """(B, 90) fixed schedule; off after day 90."""
    return lambda obs: S[:, obs["day"]] if obs["day"] < S.shape[1] else torch.zeros(S.shape[0], dtype=torch.long)


# Candidate subset for the grid: single-action budget-spending blocks starting every
# 5 days (plus the latest start). The full 417-candidate script-67 set is run on
# the primary cell only, and the gap between the two is reported.
SUBSET_IDX = [k for k, c in enumerate(s67.CANDS) if c["kind"] == "block"
              and ((c["start"] - 1) % 5 == 0 or c["start"] == max(d["start"] for d in s67.CANDS
                                                                   if d["kind"] == "block"
                                                                   and d["action"] == c["action"]))]


def best_fixed(rho, kill, cfg, full: bool = False) -> torch.Tensor:
    """Per-patient TTP-best of the script-67 candidate schedules (ties: smaller day-365 volume)."""
    C = s67.CAND_SCHED if full else s67.CAND_SCHED[SUBSET_IDX]
    P, Kc = len(rho), len(C)
    ip, ik = torch.arange(P).repeat_interleave(Kc), torch.arange(Kc).repeat(P)
    T, VH = torch.empty(P * Kc), torch.empty(P * Kc)
    chunk = 4000
    for i in range(0, P * Kc, chunk):
        sl = slice(i, i + chunk)
        sim = TwoPopSim(rho[ip[sl]].tolist(), kill[ip[sl]], cfg["f_r0"], cfg["cost"], cfg["eps"], cfg["occupancy"])
        r = rollout(sim, schedule(C[ik[sl]]), 90)
        T[sl], VH[sl] = r["ttp"], r["v_horizon"]
    T, VH = T.view(P, Kc), VH.view(P, Kc)
    top = T == T.max(1, keepdim=True).values
    return C[torch.where(top, VH, torch.full_like(VH, float("inf"))).argmin(1)]


def paired(t: np.ndarray, ref: np.ndarray, rng) -> Dict:
    d = t - ref
    boot = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)])
    return {"mean_diff_days": float(d.mean()),
            "ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            "n_longer_equal_shorter": [int((d > 0).sum()), int((d == 0).sum()), int((d < 0).sum())],
            "wilcoxon_p": float(stats.wilcoxon(t, ref).pvalue) if np.any(d != 0) else None}


def validate() -> float:
    st = s68.reaction_sets()["cohort64"]
    a = s68.rollout_ttp(s67.KillReactionSim(st["rho"].tolist(), st["kill"]), s66.schedule_policy(STUPP_SCHEDULE))
    b = rollout(TwoPopSim(st["rho"].tolist(), st["kill"], 0.0, 0.25, 0.0, "seed59"),
                schedule(torch.tensor([STUPP_SCHEDULE] * len(st["rho"]))), 90)
    err = float((a["ttp"] - b["ttp"]).abs().max())
    return err


def all_configs():
    keys = list(GRID)
    cfgs = [dict(zip(keys, v)) for v in itertools.product(*[GRID[k] for k in keys])]
    cfgs += [{"f_r0": 0.0, "cost": 0.0, "eps": 0.0, "occupancy": o, "window": w}
             for o in GRID["occupancy"] for w in GRID["window"]]
    return cfgs


def cell_path(cfg) -> Path:
    tag = "_".join(f"{k}{cfg[k]}" for k in GRID)
    return OUTPUT_DIR / "cells" / f"{tag}.json"


def run_cell(idx: int, cfg: Dict, sets: Dict, fixed_cache: Dict) -> Dict:
    rng = np.random.default_rng(7600 + idx)   # per-cell seed: shards reproduce a serial run
    cell = {"config": cfg, "control": cfg["f_r0"] == 0.0, "sets": {}}
    for set_name, st in sets.items():
        rho, kill = st["rho"], st["kill"]
        fkey = (set_name, cfg["f_r0"], cfg["cost"], cfg["eps"], cfg["occupancy"])
        if fkey not in fixed_cache:   # fixed schedules act on days 1-90 only: window-independent
            fixed_cache[fkey] = best_fixed(rho, kill, cfg)
        mk = lambda: TwoPopSim(rho.tolist(), kill, cfg["f_r0"], cfg["cost"], cfg["eps"], cfg["occupancy"])  # noqa
        B = len(rho)
        arms = {"stupp": (schedule(torch.tensor([STUPP_SCHEDULE] * B)), 90),
                "mtd_front": (lambda obs: torch.full((obs["vol"].shape[0],), 3, dtype=torch.long), cfg["window"]),
                "paced": (paced_policy(cfg["window"]), cfg["window"]),
                "AT50": (at_policy(0.5), cfg["window"]),
                "AT80": (at_policy(0.8), cfg["window"]),
                "best_fixed": (schedule(fixed_cache[fkey]), 90)}
        if cfg == PRIMARY:
            arms["best_fixed_full417"] = (schedule(best_fixed(rho, kill, cfg, full=True)), 90)
        res = {n: rollout(mk(), p, w) for n, (p, w) in arms.items()}
        T = {n: r["ttp"].numpy() for n, r in res.items()}
        best_na = np.max(np.stack([T[n] for n in ("stupp", "mtd_front", "best_fixed", "paced")]), 0)
        cell["sets"][set_name] = {"n": B, "best_nonadaptive_rmst_days": float(best_na.mean()), "arms": {
            n: {"rmst_days": float(r["ttp"].mean()), "censored_pct": float(r["ttp_censored"].double().mean() * 100),
                "mean_auc": float(r["auc"].mean()),
                "resistant_fraction_day365_median": float(r["resistant_fraction_end"].median()),
                "vs_stupp": paired(T[n], T["stupp"], rng) if n != "stupp" else None,
                "vs_paced": paired(T[n], T["paced"], rng) if n.startswith("AT") else None,
                "vs_best_nonadaptive": paired(T[n], best_na, rng) if n.startswith("AT") else None,
                "vs_best_fixed": paired(T[n], T["best_fixed"], rng) if n == "best_fixed_full417" else None}
            for n, r in res.items()}}
    return cell


def compute(shard: int, nshards: int, reverse: bool = False, only=None):
    err = validate()
    assert err == 0.0, f"TwoPopSim(f_r0=0) TTP differs from KillReactionSim by {err} days"
    (OUTPUT_DIR / "cells").mkdir(exist_ok=True)
    sets = s68.reaction_sets()
    fixed_cache = {}
    t_all = time.time()
    order = list(enumerate(all_configs()))
    for idx, cfg in (reversed(order) if reverse else order):   # reverse: a helper process for a slow shard
        if (only is not None and idx not in only) or (only is None and idx % nshards != shard)                 or cell_path(cfg).exists():
            continue
        cell = run_cell(idx, cfg, sets, fixed_cache)
        cell["validation_f_r0_0_max_ttp_diff_days"] = err
        cell_path(cfg).write_text(json.dumps(cell, indent=1))
        rs = cell["sets"]["real_test"]["arms"]
        print(f"[{idx}] {cfg}  real_test RMST stupp {rs['stupp']['rmst_days']:.1f} paced {rs['paced']['rmst_days']:.1f} "
              f"AT50 {rs['AT50']['rmst_days']:.1f} AT80 {rs['AT80']['rmst_days']:.1f} "
              f"best_fixed {rs['best_fixed']['rmst_days']:.1f} ({time.time() - t_all:.0f}s)")


def merge():
    cfgs = all_configs()
    missing = [c for c in cfgs if not cell_path(c).exists()]
    assert not missing, f"{len(missing)} cells not computed yet, e.g. {missing[0]}"
    cells = [json.loads(cell_path(c).read_text()) for c in cfgs]
    sets = list(cells[0]["sets"])
    out = {"budget": BUDGET, "horizon_day": HORIZON, "grid": GRID, "primary": PRIMARY,
           "validation_f_r0_0_max_ttp_diff_days": cells[0]["validation_f_r0_0_max_ttp_diff_days"],
           "best_fixed_candidates": {"grid_subset": len(SUBSET_IDX), "primary_full": len(s67.CANDS)},
           "cells": cells}
    # Summary per set: adaptive vs same-window paced, vs best non-adaptive, and the
    # same difference in the matching f_r0 = 0 control cell (difference-in-differences)
    ctrl = {(c["config"]["occupancy"], c["config"]["window"]): c for c in out["cells"] if c["control"]}
    summ = {}
    for set_name in sets:
        rows = []
        for c in out["cells"]:
            if c["control"]:
                continue
            k = (c["config"]["occupancy"], c["config"]["window"])
            for at in ("AT50", "AT80"):
                a = c["sets"][set_name]["arms"][at]
                a0 = ctrl[k]["sets"][set_name]["arms"][at]
                rows.append({**c["config"], "arm": at,
                             "vs_stupp": a["vs_stupp"]["mean_diff_days"], "vs_stupp_ci95": a["vs_stupp"]["ci95"],
                             "vs_paced": a["vs_paced"]["mean_diff_days"], "vs_paced_ci95": a["vs_paced"]["ci95"],
                             "vs_paced_in_control": a0["vs_paced"]["mean_diff_days"],
                             "resistance_attributable_days": a["vs_paced"]["mean_diff_days"]
                             - a0["vs_paced"]["mean_diff_days"],
                             "vs_best_nonadaptive": a["vs_best_nonadaptive"]["mean_diff_days"],
                             "vs_best_nonadaptive_ci95": a["vs_best_nonadaptive"]["ci95"]})
        summ[set_name] = {
            "n_cells_x_arms": len(rows),
            "adaptive_beats_stupp_ci": sum(r["vs_stupp_ci95"][0] > 0 for r in rows),
            "adaptive_beats_paced_ci": sum(r["vs_paced_ci95"][0] > 0 for r in rows),
            "adaptive_loses_to_paced_ci": sum(r["vs_paced_ci95"][1] < 0 for r in rows),
            "adaptive_beats_best_nonadaptive_ci": sum(r["vs_best_nonadaptive_ci95"][0] > 0 for r in rows),
            "max_resistance_attributable_days": max(r["resistance_attributable_days"] for r in rows),
            "rows": rows}
    out["summary"] = summ
    path = OUTPUT_DIR / "evaluate.json"
    path.write_text(json.dumps(out, indent=2))
    print(json.dumps({k: {kk: v for kk, v in s.items() if kk != "rows"} for k, s in summ.items()}, indent=2))
    print(f"[saved] {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshards", type=int, default=1)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--merge", action="store_true")
    ap.add_argument("--reverse", action="store_true")
    ap.add_argument("--only", type=int, nargs="*", help="cell indices to compute (helper processes)")
    a = ap.parse_args()
    if a.threads:
        torch.set_num_threads(a.threads)
    if a.merge:
        merge()
    else:
        compute(a.shard, a.nshards, a.reverse, a.only)
        if a.nshards == 1 and a.only is None:
            merge()


if __name__ == "__main__":
    main()
