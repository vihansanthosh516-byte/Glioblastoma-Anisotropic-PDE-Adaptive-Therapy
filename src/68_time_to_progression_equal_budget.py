#!/usr/bin/env python3
"""
Script 68: Time to progression at an equal drug budget
======================================================
Scripts 59-67 score every arm by tumour volume on day 90, the last treatment
day. That endpoint rewards spending the whole budget at the end (the day-90
oracle is always a late block), and above rho* ~0.075/day it rewards delaying
treatment while the tumour grows untreated (script 62). Neither is a clinical
endpoint.

Here the endpoint is time to progression (TTP) under RANO 2.0 volumetric
criteria: progression on the first day the volume is >= 40% above its nadir
(the smallest volume so far, baseline included; Wen et al. 2023). Treatment is
confined to days 1-90 at the same drug budget (Stupp's AUC, script-59 budget
rule); every arm is then followed untreated to day 365, where TTP is censored.
Reported as restricted mean TTP to day 365 and as paired differences vs Stupp.
Model-based TTP optimisation at equal total dose has found gains of about a
week (Chaudhuri et al. 2023), so effects of that size are what to expect.

Stages (run all with no flags, or one with --stage):
  oracle    per-patient TTP-optimal schedule over the script-67 candidates
            (single blocks at every start day, two-action packed late blocks),
            CEM check; which schedules win and why
  train     DAgger policy (kill-conditioned, script-67 features) from the TTP oracle
  evaluate  all arms on the script-67 evaluation sets (reaction-only sim, exact
            voxel classes) + full-PDE check on the script-64 cohort
  report    figure

Outputs: output/ttp_equal_budget/{oracle,train,evaluate}.json,
         output/ttp_equal_budget/ttp_equal_budget.png
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from typing import Callable, Dict, List, Optional

import numpy as np
import torch
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
OUTPUT_DIR = PROJECT_ROOT / "output" / "ttp_equal_budget"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

_spec = spec_from_file_location("s67", PROJECT_ROOT / "src" / "67_kill_conditioned_policy.py")
s67 = module_from_spec(_spec)
_spec.loader.exec_module(s67)
s66 = s67.s66

N_DAYS, BUDGET, STUPP_SCHEDULE = s66.N_DAYS, s66.BUDGET, s66.STUPP_SCHEDULE
HORIZON = 365          # follow-up end (day); TTP censored here
PROG_FACTOR = 1.4      # RANO 2.0 volumetric progression: >= 40% above nadir
N_TTP_TRAIN, N_TTP_VAL = 500, 200


# --------------------------------------------------------------------------- #
# Rollout with untreated follow-up and the TTP endpoint
# --------------------------------------------------------------------------- #
def rollout_ttp(sim, policy: Callable) -> Dict[str, torch.Tensor]:
    """s66.rollout for days 1-90 (budget rule), then untreated to HORIZON."""
    res = s66.rollout(sim, policy)
    off = torch.zeros(sim.B, dtype=torch.long)
    extra = [sim.step(off) for _ in range(HORIZON - N_DAYS)]
    V = torch.cat([res["v0"].view(-1, 1), res["vols"], torch.stack(extra, 1)], 1)  # V[:, d] = volume on day d
    return {**res, **ttp_endpoints(V), "V": V}


def first_day(mask: torch.Tensor) -> (torch.Tensor, torch.Tensor):
    hit = mask.any(1)
    day = torch.where(hit, mask.int().argmax(1), torch.full((mask.shape[0],), mask.shape[1] - 1))
    return day.double(), ~hit


def ttp_endpoints(V: torch.Tensor) -> Dict[str, torch.Tensor]:
    nadir = torch.cummin(V, 1).values
    prog = V >= PROG_FACTOR * nadir
    prog[:, 0] = False
    ttp, cens = first_day(prog)
    base = V >= PROG_FACTOR * V[:, :1]
    base[:, 0] = False
    ttb, cens_b = first_day(base)
    return {"ttp": ttp, "ttp_censored": cens, "time_to_1p4x_baseline": ttb, "ttb_censored": cens_b,
            "v_horizon": V[:, -1], "nadir": nadir[:, -1], "burden_horizon": V.mean(1)}


def ttp_oracle(rho: torch.Tensor, kill: torch.Tensor, chunk: int = 4000) -> Dict:
    """Per patient: the candidate with the longest TTP; ties (e.g. all censored) are
    broken by the smallest volume at the horizon."""
    C = s67.CAND_SCHED
    P, K = len(rho), len(C)
    ip, ik = torch.arange(P).repeat_interleave(K), torch.arange(K).repeat(P)
    T, VH = torch.empty(P * K), torch.empty(P * K)
    for i in range(0, P * K, chunk):
        sl = slice(i, i + chunk)
        r = rollout_ttp(s67.KillReactionSim(rho[ip[sl]].tolist(), kill[ip[sl]]), s66.sched_tensor_policy(C[ik[sl]]))
        T[sl], VH[sl] = r["ttp"], r["v_horizon"]
    T, VH = T.view(P, K), VH.view(P, K)
    top = T == T.max(1, keepdim=True).values
    best = torch.where(top, VH, torch.full_like(VH, float("inf"))).argmin(1)
    return {"best": best, "sched": C[best], "ttp": T, "v_horizon": VH}


def save_json(obj, name: str):
    path = OUTPUT_DIR / name
    path.write_text(json.dumps(obj, indent=2))
    print(f"[saved] {path}")


def describe(best: torch.Tensor) -> List[Dict]:
    return [{k: v for k, v in s67.CANDS[b].items() if k != "sched"} for b in best.tolist()]


# --------------------------------------------------------------------------- #
# Stage: oracle
# --------------------------------------------------------------------------- #
def reaction_sets() -> Dict[str, Dict]:
    """rho and kill tables of the script-67 evaluation sets (no solvers built)."""
    return {name: {"rho": torch.tensor(st["rho"], dtype=torch.float64), "kill": st["kill"], "ids": st["ids"],
                   "source": st["source"]}
            for name, st in s67.eval_sets().items()}


def stage_oracle():
    out = {"budget": BUDGET, "horizon_day": HORIZON, "progression_factor": PROG_FACTOR,
           "n_candidates": len(s67.CANDS), "sim": "reaction-only 64^3 (exact voxel classes)", "sets": {}}
    store = {}
    pools = torch.load(s67.OUTPUT_DIR / "oracle_pools.pt")
    sets = reaction_sets()
    sets["ttp_train"] = {"rho": pools["train"]["rho"][:N_TTP_TRAIN], "kill": pools["train"]["kill"][:N_TTP_TRAIN]}
    sets["ttp_val"] = {"rho": pools["val"]["rho"][:N_TTP_VAL], "kill": pools["val"]["kill"][:N_TTP_VAL]}
    for name, st in sets.items():
        t0 = time.time()
        o = ttp_oracle(st["rho"], st["kill"])
        sim = lambda: s67.KillReactionSim(st["rho"].tolist(), st["kill"])  # noqa: E731
        stupp = rollout_ttp(sim(), s66.schedule_policy(STUPP_SCHEDULE))
        d90 = rollout_ttp(sim(), s66.sched_tensor_policy(s67.kill_oracle(st["rho"], st["kill"])["sched"]))
        best_ttp = o["ttp"].gather(1, o["best"].view(-1, 1)).squeeze(1)
        desc = describe(o["best"])
        late = {a: max(c["start"] for c in s67.CANDS if c["kind"] == "block" and c["action"] == a) for a in (1, 2, 3)}
        blocks = [d for d in desc if d["kind"] == "block"]
        out["sets"][name] = {
            "n": len(st["rho"]), "seconds": time.time() - t0,
            "oracle_kind_counts": {k: sum(d["kind"] == k for d in desc) for k in ("block", "mixed")},
            "oracle_main_action_counts": {s67.ACTION_NAMES[a]: sum(d["action"] == a for d in desc) for a in (1, 2, 3)},
            "oracle_block_start_days": [d["start"] for d in blocks],
            "oracle_block_start_is_latest_pct": 100 * float(np.mean([d["start"] == late[d["action"]] for d in blocks]))
            if blocks else None,
            "restricted_mean_ttp": {"stupp": float(stupp["ttp"].mean()), "day90_oracle": float(d90["ttp"].mean()),
                                    "ttp_oracle": float(best_ttp.mean())},
            "censored_pct": {"stupp": float(stupp["ttp_censored"].double().mean() * 100),
                             "day90_oracle": float(d90["ttp_censored"].double().mean() * 100)},
            "ttp_oracle_minus_stupp_days_mean": float((best_ttp - stupp["ttp"]).mean()),
            "day90_oracle_minus_stupp_days_mean": float((d90["ttp"] - stupp["ttp"]).mean()),
            "rho": st["rho"].tolist(), "oracle_desc": desc, "stupp_ttp": stupp["ttp"].tolist(),
            "day90_oracle_ttp": d90["ttp"].tolist(), "ttp_oracle_ttp": best_ttp.tolist(),
        }
        store[name] = {"rho": st["rho"], "kill": st["kill"], "sched": o["sched"], "ttp": best_ttp,
                       "stupp_ttp": stupp["ttp"]}
        s = out["sets"][name]
        print(f"  [{name}] n={s['n']} ({s['seconds']:.0f}s)  RMST stupp {s['restricted_mean_ttp']['stupp']:.1f}  "
              f"day90-oracle {s['restricted_mean_ttp']['day90_oracle']:.1f}  "
              f"TTP-oracle {s['restricted_mean_ttp']['ttp_oracle']:.1f}  kinds {s['oracle_kind_counts']}  "
              f"late {s['oracle_block_start_is_latest_pct']}")
    torch.save(store, OUTPUT_DIR / "ttp_oracle_pools.pt")

    # CEM check on nine validation patients spanning rho; objective TTP, tie-break volume at horizon
    v = store["ttp_val"]
    pick = v["rho"].argsort()[torch.linspace(0, len(v["rho"]) - 1, 9).long()].tolist()
    out["cem_check"] = {"val_indices": pick, "rho": v["rho"][pick].tolist()}
    for mode, init in (("cold", None), ("warm_from_oracle", v["sched"][pick])):
        t0 = time.time()
        c = cem_ttp(v["rho"][pick].tolist(), v["kill"][pick], init=init)
        diff = (c["best_ttp"] - v["ttp"][pick]).tolist()
        out["cem_check"][mode] = {"cem_minus_oracle_ttp_days": diff, "cem_best_sched": c["best_sched"].tolist(),
                                  "seconds": time.time() - t0}
        print(f"  CEM {mode} ({time.time() - t0:.0f}s): TTP diff {diff}")
    save_json(out, "oracle.json")


def cem_ttp(rho: List[float], kill: torch.Tensor, n: int = 256, iters: int = 30, elite: float = 0.1,
            alpha: float = 0.7, seed: int = 0, init: Optional[torch.Tensor] = None) -> Dict:
    """s67.cem_kill with the TTP objective (longest TTP, then smallest horizon volume)."""
    F = torch.nn.functional
    g = torch.Generator().manual_seed(seed)
    P = len(rho)
    probs = torch.full((P, N_DAYS, 4), 0.25)
    if init is not None:
        probs = 0.5 * F.one_hot(init.long(), 4).double() + 0.5 * probs
    best_ttp = torch.full((P,), -1.0)
    best_vh = torch.full((P,), float("inf"))
    best_sched = torch.zeros(P, N_DAYS, dtype=torch.long)
    n_el = max(1, int(n * elite))
    for _ in range(iters):
        S = torch.multinomial(probs.repeat_interleave(n, 0).view(-1, 4), 1, generator=g).view(P * n, N_DAYS)
        r = rollout_ttp(s67.KillReactionSim([x for x in rho for _ in range(n)], kill.repeat_interleave(n, 0)),
                        s66.sched_tensor_policy(S))
        T, VH = r["ttp"].view(P, n), r["v_horizon"].view(P, n)
        score = T - 1e-3 * torch.log(VH / VH.min(1, keepdim=True).values)   # TTP is in whole days
        exe = r["actions"].view(P, n, N_DAYS)
        order = score.argsort(1, descending=True)[:, :n_el]
        el = torch.gather(exe, 1, order.unsqueeze(2).expand(-1, -1, N_DAYS))
        probs = alpha * F.one_hot(el, 4).double().mean(1) + (1 - alpha) * probs
        j = order[:, 0]
        t_b, v_b = T.gather(1, j.view(-1, 1)).squeeze(1), VH.gather(1, j.view(-1, 1)).squeeze(1)
        better = (t_b > best_ttp) | ((t_b == best_ttp) & (v_b < best_vh))
        best_ttp = torch.where(better, t_b, best_ttp)
        best_vh = torch.where(better, v_b, best_vh)
        best_sched[better] = exe[better, j[better]]
    return {"best_ttp": best_ttp, "best_v_horizon": best_vh, "best_sched": best_sched}


# --------------------------------------------------------------------------- #
# Stage: train  (DAgger from the TTP oracle, kill-conditioned features)
# --------------------------------------------------------------------------- #
def eval_ttp_pool(pol_factory: Callable, pool: Dict) -> Dict:
    r = rollout_ttp(s67.KillReactionSim(pool["rho"].tolist(), pool["kill"]), pol_factory(pool["kill"]))
    d = r["ttp"] - pool["stupp_ttp"]
    return {"rmst": float(r["ttp"].mean()), "ttp_minus_stupp_days": float(d.mean()),
            "ttp_minus_oracle_days": float((r["ttp"] - pool["ttp"]).mean()),
            "win_vs_stupp_pct": float((d > 0).double().mean() * 100), "mean_auc": float(r["auc"].mean())}


def stage_train():
    store = torch.load(OUTPUT_DIR / "ttp_oracle_pools.pt")
    train, val = store["ttp_train"], store["ttp_val"]
    ref = {"stupp": eval_ttp_pool(lambda k: s66.schedule_policy(STUPP_SCHEDULE), val),
           "efficiency_rule": eval_ttp_pool(lambda k: s66.sched_tensor_policy(s67.efficiency_rule_sched(k)), val),
           "dagger_cond_day90": eval_ttp_pool(lambda k: s67.greedy(s67.load_run("dagger_cond"), k, True), val)}
    print(f"  val reference: {ref}")
    t0 = time.time()
    net, log = s67.train_dagger(True, train, val,
                                val_eval=lambda n: eval_ttp_pool(lambda k: s67.greedy(n, k, True), val))
    torch.save(net.state_dict(), OUTPUT_DIR / "dagger_ttp.pt")
    save_json({"val_reference": ref, "dagger_ttp": log, "seconds": time.time() - t0,
               "n_train": len(train["rho"]), "n_val": len(val["rho"])}, "train.json")


# --------------------------------------------------------------------------- #
# Stage: evaluate
# --------------------------------------------------------------------------- #
def arms_for(kill: torch.Tensor, rho: torch.Tensor, oracle_sched: torch.Tensor) -> Dict[str, Callable]:
    ttp_net = s67.make_net(True)
    ttp_net.load_state_dict(torch.load(OUTPUT_DIR / "dagger_ttp.pt"))
    return {"stupp": s66.schedule_policy(STUPP_SCHEDULE), "heuristic59": s66.heuristic59_policy,
            "ppo66": s66.greedy_policy(s66.load_net("ppo_final_s0")),
            "dagger66": s66.greedy_policy(s66.load_net("bc_dagger_final")),
            "efficiency_rule": s66.sched_tensor_policy(s67.efficiency_rule_sched(kill)),
            "day90_oracle": s66.sched_tensor_policy(s67.kill_oracle(rho, kill)["sched"]),
            "dagger_cond_day90": s67.greedy(s67.load_run("dagger_cond"), kill, True),
            "ppo_cond_ft_day90": s67.greedy(s67.load_run("ppo_cond_ft"), kill, True),
            "ttp_oracle": s66.sched_tensor_policy(oracle_sched),
            "dagger_ttp": s67.greedy(ttp_net, kill, True)}


def arm_stats(r: Dict, ref: Dict, rng: np.random.Generator) -> Dict:
    t, ts = r["ttp"].numpy(), ref["ttp"].numpy()
    d = t - ts
    boot = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(2000)])
    out = {"rmst_days": float(t.mean()), "median_ttp_days": float(np.median(t)),
           "censored_pct": float(r["ttp_censored"].double().mean() * 100),
           "ttp_minus_stupp_days_mean": float(d.mean()),
           "ttp_minus_stupp_days_ci95": [float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
           "n_longer_ttp": int((d > 0).sum()), "n_equal_ttp": int((d == 0).sum()), "n_shorter_ttp": int((d < 0).sum()),
           "wilcoxon_p_vs_stupp": None,
           "time_to_1p4x_baseline_mean": float(r["time_to_1p4x_baseline"].mean()),
           "log_ratio_volume_horizon_vs_stupp": float(torch.log(r["v_horizon"] / ref["v_horizon"]).mean()),
           "log_ratio_volume_day90_vs_stupp": float(torch.log(r["final"] / ref["final"]).mean()),
           "mean_drug_auc": float(r["auc"].mean()), "max_drug_auc": float(r["auc"].max()),
           "ttp_days": t.tolist()}
    if np.any(d != 0):
        out["wilcoxon_p_vs_stupp"] = float(stats.wilcoxon(t, ts).pvalue)
    return out


def stage_evaluate(full_pde: bool = True):
    store = torch.load(OUTPUT_DIR / "ttp_oracle_pools.pt")
    rng = np.random.default_rng(68)
    out = {"budget": BUDGET, "horizon_day": HORIZON, "progression_factor": PROG_FACTOR,
           "sim": "reaction-only 64^3 (exact voxel classes) for all sets; full PDE check on cohort64", "sets": {}}
    for name, st in reaction_sets().items():
        t0 = time.time()
        arms = arms_for(st["kill"], st["rho"], store[name]["sched"])
        res = {a: rollout_ttp(s67.KillReactionSim(st["rho"].tolist(), st["kill"]), pol) for a, pol in arms.items()}
        blk = {a: arm_stats(r, res["stupp"], rng) for a, r in res.items()}
        out["sets"][name] = {"source": st["source"], "n": len(st["rho"]), "ids": st["ids"], "rho": st["rho"].tolist(),
                             "arms": blk, "actions": {a: r["actions"].tolist() for a, r in res.items()}}
        print(f"\n  [{name}] n={len(st['rho'])} ({time.time() - t0:.0f}s)")
        for a, s in blk.items():
            print(f"    {a:18s} RMST {s['rmst_days']:6.1f}  vs Stupp {s['ttp_minus_stupp_days_mean']:+6.1f} d "
                  f"[{s['ttp_minus_stupp_days_ci95'][0]:+.1f}, {s['ttp_minus_stupp_days_ci95'][1]:+.1f}]  "
                  f"+{s['n_longer_ttp']}/={s['n_equal_ttp']}/-{s['n_shorter_ttp']}  cens {s['censored_pct']:.0f}%  "
                  f"AUC max {s['max_drug_auc']:.2f}  p {s['wilcoxon_p_vs_stupp']}")
        save_json(out, "evaluate.json")

    if full_pde:
        # Diffusion conserves tumour mass, so the reaction-only sim should match the full
        # PDE closely; check that over 365 days on the script-64 cohort with its own solver.
        st = s67.eval_sets()["cohort64"]
        err = st["validate"]()
        assert err < 1e-9, f"cohort64 batched solver disagrees with numpy solver ({err})"
        arms = arms_for(st["kill"], torch.tensor(st["rho"]), store["cohort64"]["sched"])
        keep = ["stupp", "heuristic59", "efficiency_rule", "dagger_cond_day90", "ppo_cond_ft_day90",
                "ttp_oracle", "dagger_ttp"]
        t0 = time.time()
        solver = st["solver"]()
        full = {a: rollout_ttp(solver, arms[a]) for a in keep}
        red = out["sets"]["cohort64"]["arms"]
        out["full_pde_check_cohort64"] = {
            "validation_max_rel_err_vs_numpy_solver": err, "seconds": time.time() - t0,
            "arms": {a: {**arm_stats(full[a], full["stupp"], rng),
                         "ttp_max_abs_diff_vs_reaction_days": float(np.max(np.abs(full[a]["ttp"].numpy()
                                                                                  - np.array(red[a]["ttp_days"])))),
                         } for a in keep}}
        for a, s in out["full_pde_check_cohort64"]["arms"].items():
            print(f"  [full PDE cohort64] {a:18s} RMST {s['rmst_days']:6.1f}  vs Stupp "
                  f"{s['ttp_minus_stupp_days_mean']:+6.1f} d  max |TTP full - reaction| "
                  f"{s['ttp_max_abs_diff_vs_reaction_days']:.0f} d")
        save_json(out, "evaluate.json")


# --------------------------------------------------------------------------- #
# Stage: report
# --------------------------------------------------------------------------- #
def stage_report():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ev = json.loads((OUTPUT_DIR / "evaluate.json").read_text())
    orc = json.loads((OUTPUT_DIR / "oracle.json").read_text())
    arms = ["heuristic59", "ppo66", "dagger66", "efficiency_rule", "day90_oracle", "dagger_cond_day90",
            "ppo_cond_ft_day90", "dagger_ttp", "ttp_oracle"]
    sets = list(ev["sets"])
    fig, axes = plt.subplots(2, 2, figsize=(17, 12))

    ax = axes[0, 0]
    w = 0.8 / len(sets)
    for j, s in enumerate(sets):
        m = np.array([ev["sets"][s]["arms"][a]["ttp_minus_stupp_days_mean"] for a in arms])
        ci = np.array([ev["sets"][s]["arms"][a]["ttp_minus_stupp_days_ci95"] for a in arms])
        x = np.arange(len(arms)) + (j - (len(sets) - 1) / 2) * w
        ax.bar(x, m, w, label=f"{s} (n={ev['sets'][s]['n']})")
        ax.errorbar(x, m, yerr=[m - ci[:, 0], ci[:, 1] - m], fmt="none", ecolor="k", lw=0.8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels(arms, rotation=35, ha="right")
    ax.set_ylabel("mean TTP - TTP_Stupp (days)")
    ax.set_title(f"A. Time to progression vs Stupp, equal drug (RANO +40% over nadir, follow-up to day {HORIZON})")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3, axis="y")

    ax = axes[0, 1]
    s = orc["sets"]["ttp_val"]
    rho = np.array(s["rho"])
    starts = np.array([d["start"] for d in s["oracle_desc"]])
    cols = {1: "#fdae61", 2: "#abd9e9", 3: "#2c7bb6"}
    ax.scatter(rho, starts, c=[cols[d["action"]] for d in s["oracle_desc"]], edgecolors="k", s=30)
    ax.set_xscale("log")
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("oracle treatment start day")
    ax.set_title("B. TTP-optimal start day by growth rate (colour = main action: chemo/rad/combo)")
    ax.grid(alpha=0.3)

    ax = axes[1, 0]
    ax.scatter(rho, np.array(s["ttp_oracle_ttp"]) - np.array(s["stupp_ttp"]), s=20, label="TTP oracle - Stupp")
    ax.scatter(rho, np.array(s["day90_oracle_ttp"]) - np.array(s["stupp_ttp"]), s=20, marker="x",
               label="day-90 oracle - Stupp")
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xscale("log")
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("TTP difference (days)")
    ax.set_title("C. Per-patient TTP gain at equal drug (validation pool)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    days = np.arange(HORIZON + 1)
    for a in ("stupp", "efficiency_rule", "day90_oracle", "ttp_oracle", "dagger_ttp"):
        t = np.array(ev["sets"]["cohort64"]["arms"][a]["ttp_days"])
        ax.step(days, [(t > d).mean() for d in days], where="post", label=a)
    ax.set_xlabel("day")
    ax.set_ylabel("fraction progression-free")
    ax.set_title("D. Progression-free fraction, script-64 cohort (n=20)")
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    plt.suptitle(f"Script 68: time to progression at equal drug budget (AUC {BUDGET:g})", fontsize=14,
                 fontweight="bold")
    plt.tight_layout()
    path = OUTPUT_DIR / "ttp_equal_budget.png"
    plt.savefig(path, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"[saved] {path}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["oracle", "train", "evaluate", "report"])
    ap.add_argument("--no-full-pde", action="store_true")
    args = ap.parse_args()
    print("=" * 70)
    print("SCRIPT 68: TIME TO PROGRESSION AT EQUAL DRUG BUDGET")
    print("=" * 70)
    stages = [args.stage] if args.stage else ["oracle", "train", "evaluate", "report"]
    for st in stages:
        t0 = time.time()
        print(f"\n--- stage {st} ---")
        {"oracle": stage_oracle, "train": stage_train,
         "evaluate": lambda: stage_evaluate(not args.no_full_pde), "report": stage_report}[st]()
        print(f"--- stage {st} done in {time.time() - t0:.0f}s ---")


if __name__ == "__main__":
    main()
