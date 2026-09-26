#!/usr/bin/env python3
"""
Script 65: Track C Final Report (equal-drug-budget framework)
=============================================================
Aggregates the Track C outputs into
  output/final_executive_summary.json
  output/65_master_summary_figure.png

Inputs (a missing or unreadable input is recorded as such, never defaulted):
  51 inverse_est_metrics.json            inverse parameter estimation
  52 robust_mpc_benchmark.json           robust MPC benchmark
  59 phase6_sensitivity_metrics.json     drug-budget sensitivity sweep
  60 ablation_and_baselines_metrics.json equal-budget baselines + ablation
  62 biomarker_stability_metrics.json    early-vs-delayed start rho threshold
  64 phase8_cohort_metrics.json          equal-budget virtual cohort
  66 rl_equal_budget/evaluate.json       RL equal-budget study

Every treatment comparison must use the same drug budget as the shared arms
module (src/rl/equal_budget_arms.py); the report checks this.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
from rl.equal_budget_arms import ARMS, STUPP_BUDGET  # noqa: E402

OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INPUTS = {
    "51_inverse_estimation": "inverse_est_metrics.json",
    "52_robust_mpc": "robust_mpc_benchmark.json",
    "59_drug_budget_sensitivity": "phase6_sensitivity_metrics.json",
    "60_ablation_and_baselines": "ablation_and_baselines_metrics.json",
    "62_early_start_threshold": "biomarker_stability_metrics.json",
    "64_virtual_cohort": "phase8_cohort_metrics.json",
    "66_rl_equal_budget": "rl_equal_budget/evaluate.json",
}
ARM_LABELS = {"heuristic_budgeted": "Budgeted heuristic (59)", "ppo": "RL: PPO (66)",
              "dagger_oracle": "DAgger oracle (66)"}
S66_ARM_NAMES = {"heuristic_budgeted": "heuristic59", "ppo": "ppo_final_s0", "dagger_oracle": "bc_dagger_final"}


def load_inputs() -> (Dict[str, Any], Dict[str, str]):
    data, status = {}, {}
    for key, rel in INPUTS.items():
        path = OUTPUT_DIR / rel
        if not path.exists():
            status[key] = "missing"
        elif path.stat().st_size == 0:
            status[key] = "empty"
        else:
            try:
                data[key] = json.loads(path.read_text())
                status[key] = "ok"
            except json.JSONDecodeError as e:
                status[key] = f"unreadable: {e}"
        print(f"  {key:28s} {rel:40s} {status[key]}")
    return data, status


def arm_row(win: float, log_ratio: float, auc_ratio: float) -> Dict[str, float]:
    return {"win_rate_vs_stupp_pct": win, "mean_log_ratio_final_vs_stupp": log_ratio,
            "drug_auc_ratio_vs_stupp": auc_ratio}


def build_summary(d: Dict[str, Any], status: Dict[str, str]) -> Dict[str, Any]:
    s: Dict[str, Any] = {
        "report_framework": "equal drug budget (Stupp AUC) for every treatment arm",
        "budget_drug_auc": STUPP_BUDGET,
        "arms": list(ARMS),
        "input_status": status,
        "missing_inputs": [k for k, v in status.items() if v != "ok"],
        "pipeline_status": "COMPLETE" if all(v == "ok" for v in status.values()) else "INCOMPLETE",
    }

    budgets = {}
    if "59_drug_budget_sensitivity" in d:
        m = d["59_drug_budget_sensitivity"]
        budgets["59"] = m.get("mean_stupp_drug_auc")
        s["59_drug_budget_sensitivity"] = {k: m.get(k) for k in (
            "rl_win_rate_pct", "mean_rl_volume_mm3", "mean_stupp_volume_mm3", "mean_rl_drug_auc",
            "mean_stupp_drug_auc", "drug_auc_ratio", "top_sensitive_parameter")}
    if "60_ablation_and_baselines" in d:
        m = d["60_ablation_and_baselines"]
        budgets["60"] = m.get("budget_drug_auc")
        arms = m["baselines"]["arms"]
        s["60_ablation_and_baselines"] = {
            "n_scenarios": m["n_scenarios"],
            "arms_full_model": {a: arm_row(arms[a]["win_rate_vs_stupp_pct"], arms[a]["mean_log_ratio_final_vs_stupp"],
                                           arms[a]["drug_auc_ratio_vs_stupp"]) for a in ARM_LABELS},
            "summary": m["summary"], "notes": m.get("notes", [])}
    if "62_early_start_threshold" in d:
        m = d["62_early_start_threshold"]
        budgets["62"] = m.get("budget_drug_auc")
        s["62_early_start_threshold"] = {"cohort_n": m["cohort"]["n_patients"], **m["summary"]}
    if "64_virtual_cohort" in d:
        m = d["64_virtual_cohort"]
        budgets["64"] = m.get("budget_drug_auc")
        s["64_virtual_cohort"] = {
            "cohort_size": m["cohort_size"],
            "arms": {a: {**arm_row(m["arms"][a]["win_rate_vs_stupp_pct"], m["arms"][a]["mean_log_ratio_final_vs_stupp"],
                                   m["arms"][a]["drug_auc_ratio_vs_stupp"]),
                         "mean_log_ratio_ci95": m["arms"][a]["mean_log_ratio_final_vs_stupp_ci95"],
                         "wilcoxon_p_value": m["arms"][a]["wilcoxon_p_value"]} for a in ARM_LABELS}}
    if "66_rl_equal_budget" in d:
        m = d["66_rl_equal_budget"]
        budgets["66"] = m.get("budget")
        s["66_rl_equal_budget"] = {
            set_name: {a: arm_row(blk[n]["win_rate_final_pct"], blk[n]["mean_log_ratio_final_vs_stupp"],
                                  blk[n]["drug_auc_ratio_vs_stupp"]) for a, n in S66_ARM_NAMES.items()}
            for set_name, blk in m["sets"].items()}
    if "52_robust_mpc" in d:
        s["52_robust_mpc"] = d["52_robust_mpc"]
    if "51_inverse_estimation" in d:
        s["51_inverse_estimation"] = d["51_inverse_estimation"]

    s["budget_consistency"] = {
        "per_study_budget": budgets,
        "all_equal_to_arms_module": bool(budgets) and all(
            b is not None and abs(float(b) - STUPP_BUDGET) < 1e-9 for b in budgets.values()),
    }
    return s


def create_master_figure(d: Dict[str, Any], s: Dict[str, Any], path: Path):
    fig, axes = plt.subplots(2, 2, figsize=(18, 13))
    studies = []
    if "66_rl_equal_budget" in s:
        studies += [("66 LHS-30", s["66_rl_equal_budget"]["lhs30"]), ("66 real test", s["66_rl_equal_budget"]["real_test"])]
    if "60_ablation_and_baselines" in s:
        studies.append(("60 full model", s["60_ablation_and_baselines"]["arms_full_model"]))
    if "64_virtual_cohort" in s:
        studies.append(("64 virtual cohort", s["64_virtual_cohort"]["arms"]))

    for ax, key, ylabel, title in (
            (axes[0, 0], "win_rate_vs_stupp_pct", "win rate vs Stupp (%)", "A. Win rate vs Stupp at equal drug"),
            (axes[0, 1], "mean_log_ratio_final_vs_stupp", "mean log(V_arm / V_Stupp)  (< 0 better)",
             "B. Effect size vs Stupp at equal drug")):
        x = np.arange(len(studies))
        for i, a in enumerate(ARM_LABELS):
            ax.bar(x + (i - 1) * 0.27, [st[a][key] for _, st in studies], 0.27, label=ARM_LABELS[a])
        ax.set_xticks(x)
        ax.set_xticklabels([n for n, _ in studies])
        ax.axhline(0, color="k", lw=0.8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend(fontsize=9)
        ax.grid(alpha=0.3, axis="y")
    axes[0, 0].axhline(50, color="gray", ls=":", lw=1)

    ax = axes[1, 0]
    if "62_early_start_threshold" in d:
        m = d["62_early_start_threshold"]
        cur = m["model_threshold"]["benefit_curve"]
        ax.plot(cur["rho"], cur["benefit_final"], "k-", label="model")
        pr = [p["rho"] for p in m["patients"]]
        pb = [p["benefit_final"] for p in m["patients"]]
        ax.scatter(pr, pb, c=["green" if b > 0 else "firebrick" for b in pb], edgecolors="black", s=30, zorder=3,
                   label=f"MU-Glioma patients (n={len(pr)})")
        rs = m["summary"]["rho_star_model_per_day"]
        ax.axvline(rs, color="purple", ls="--", label=f"rho* = {rs:.4f}/day")
        if m["summary"]["rho_star_ci95_per_day"]:
            ax.axvspan(*m["summary"]["rho_star_ci95_per_day"], color="purple", alpha=0.15, label="bootstrap 95% CI")
        ax.axhline(0, color="gray", lw=0.8)
        ax.legend(fontsize=9)
    else:
        ax.text(0.5, 0.5, "script 62 output missing", ha="center", transform=ax.transAxes)
    ax.set_xlabel("rho (1/day)")
    ax.set_ylabel("log(V_Stupp / V_early) at day 90")
    ax.set_title("C. Script 62: when does starting early win at equal drug?")
    ax.grid(alpha=0.3)

    ax = axes[1, 1]
    ax.axis("off")
    lines = [f"Budget (every arm): drug AUC = {STUPP_BUDGET:g}",
             f"Budget consistent across studies: {s['budget_consistency']['all_equal_to_arms_module']}",
             f"Pipeline status: {s['pipeline_status']}"]
    if s["missing_inputs"]:
        lines.append("Missing / unreadable inputs: " + ", ".join(s["missing_inputs"]))
    if "59_drug_budget_sensitivity" in s:
        m = s["59_drug_budget_sensitivity"]
        lines.append(f"59: heuristic win {m['rl_win_rate_pct']:.0f}% at AUC ratio {m['drug_auc_ratio']:.3f}")
    if "60_ablation_and_baselines" in s:
        m = s["60_ablation_and_baselines"]["summary"]
        lines.append(f"60: PPO win {m['rl_win_rate_pct']:.0f}%, ablation impact no-DTI {m['ablation_impact_no_dti_pct']:+.2f}%,"
                     f" pure RD {m['ablation_impact_pure_rd_pct']:+.2f}%")
    if "62_early_start_threshold" in s:
        m = s["62_early_start_threshold"]
        lines.append(f"62: rho* {m['rho_star_model_per_day']:.4f}/day, early start better for "
                     f"{m['n_early_start_better_final']}/{m['cohort_n']} patients (day 90), "
                     f"{m['pct_early_start_better_burden']:.0f}% (burden)")
    if "64_virtual_cohort" in s:
        m = s["64_virtual_cohort"]["arms"]["ppo"]
        lines.append(f"64: PPO win {m['win_rate_vs_stupp_pct']:.0f}%, log-ratio {m['mean_log_ratio_final_vs_stupp']:+.3f}"
                     f" CI {['%+.3f' % v for v in m['mean_log_ratio_ci95']]}, Wilcoxon p {m['wilcoxon_p_value']:.2g}")
    if "52_robust_mpc" in s:
        m = s["52_robust_mpc"]
        lines.append(f"52: robust vs standard MPC final volume {m['robust']['final_volume_mean_std'][0]:.4f} vs "
                     f"{m['standard']['final_volume_mean_std'][0]:.4f}; variance reduction {m['variance_reduction_pct']:.1f}%")
    ax.text(0.0, 1.0, "D. Summary\n\n" + "\n\n".join(lines), va="top", fontsize=10, family="monospace", wrap=True)

    plt.suptitle("Track C Final Report: every arm at equal drug budget", fontsize=16, fontweight="bold")
    plt.tight_layout()
    plt.savefig(path, dpi=200, bbox_inches="tight")
    plt.close()
    print(f"[Plot] Saved -> {path}")


def main():
    print("=" * 70)
    print("SCRIPT 65: TRACK C FINAL REPORT (EQUAL DRUG BUDGET)")
    print("=" * 70)
    data, status = load_inputs()
    summary = build_summary(data, status)
    path = OUTPUT_DIR / "final_executive_summary.json"
    path.write_text(json.dumps(summary, indent=2))
    print(f"[Summary] Saved -> {path}")
    create_master_figure(data, summary, OUTPUT_DIR / "65_master_summary_figure.png")
    print(f"\nPipeline status: {summary['pipeline_status']}   missing: {summary['missing_inputs']}")
    print(f"Budget consistency: {summary['budget_consistency']}")


if __name__ == "__main__":
    main()
