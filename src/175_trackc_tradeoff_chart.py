#!/usr/bin/env python3
"""
Script 175: Track C mean-vs-tail trade-off chart (descriptive, post hoc).
========================================================================
SIMULATION ONLY. Reads the sealed results of scripts 171 and 173 and the 173 tuning file; runs no
simulation. For each suite and each arm it plots mean TTP (x) against CVaR_0.10 scenario mean (y).
The tuning panel shows all 10 script-173 configs on the tuning draw. Both sealed H1 tests failed;
the chart states this. No pass/fail claim is made here.

Run:  python src/175_trackc_tradeoff_chart.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
SEALED = {"171": OUT / "trackc_sealed_tail_hybrid" / "sealed_results.json",
          "173": OUT / "trackc_sealed_tail_trigger" / "sealed_results.json"}
TUNE = OUT / "trackc_sealed_tail_trigger" / "frozen_trigger.json"
OUT_DIR = OUT / "trackc_tradeoff"
ARM_NAME = {"hybrid": "front-load 4 + contain", "trigger": "front-load 4 + contain",
            "wide_containment": "wide containment", "comparator_raw_single": "comparator (raw 0.6)",
            "front_loaded": "front-loaded"}
COLOR = {"front-load 4 + contain": "#1f77b4", "wide containment": "#2ca02c",
         "comparator (raw 0.6)": "#7f7f7f", "front-loaded": "#d62728"}
MARK = {"171": "o", "173": "s"}


def main():
    res = {k: json.loads(p.read_text()) for k, p in SEALED.items()}
    suites = list(res["171"]["suites"])
    points = []
    for script, r in res.items():
        for nm in suites:
            for arm in r["arms"]:
                v = r["suites"][nm][arm]
                points.append({"script": script, "suite": nm, "arm": ARM_NAME[arm],
                               "mean_ttp": v["mean"], "cvar10": v["cvar10_scenario_mean"]})
    tune = json.loads(TUNE.read_text())
    tpts = [{"key": k, "mean_ttp": tune["tuning_mean_ttp"][k],
             "cvar10": tune["tuning_cvar10_scenario_mean"][k]} for k in tune["tuning_mean_ttp"]]

    fig, axes = plt.subplots(1, len(suites) + 1, figsize=(4.2 * (len(suites) + 1), 4.2))
    for ax, nm in zip(axes, suites):
        for p in (q for q in points if q["suite"] == nm):
            ax.scatter(p["mean_ttp"], p["cvar10"], c=COLOR[p["arm"]], marker=MARK[p["script"]], s=50,
                       edgecolors="k", linewidths=0.4)
        h1 = [res[s]["suites"][nm][[k for k in res[s]["suites"][nm] if k.startswith("H1_")][0]] for s in res]
        ax.set_title(f"{nm}\nH1 diff 171: {h1[0]['diff']:+.1f} [{h1[0]['ci95'][0]:.1f}, {h1[0]['ci95'][1]:.1f}]"
                     f"\n173: {h1[1]['diff']:+.1f} [{h1[1]['ci95'][0]:.1f}, {h1[1]['ci95'][1]:.1f}]", fontsize=8)
        ax.set_xlabel("mean TTP (days)")
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("worst-10% CVaR TTP (days)")
    ax = axes[-1]
    for t in tpts:
        hz = t["key"].split("[")[1].split(",")[2].strip()
        ax.scatter(t["mean_ttp"], t["cvar10"], marker="^" if t["key"].startswith("nominal") else "v",
                   c={"2": "#1f77b4", "3": "#ff7f0e", "4": "#9467bd"}[hz], s=40, edgecolors="k", linewidths=0.4)
    ax.set_title("173 tuning draw (not test)\nhorizon 2 blue, 3 orange, 4 purple\n^ nominal, v state-linked",
                 fontsize=8)
    ax.set_xlabel("mean TTP (days)")
    ax.grid(alpha=0.3)
    for arm, c in COLOR.items():
        axes[0].scatter([], [], c=c, label=arm, edgecolors="k", linewidths=0.4)
    for s, m in MARK.items():
        axes[0].scatter([], [], c="w", marker=m, label=f"script {s} seeds", edgecolors="k")
    axes[0].legend(fontsize=7, loc="center left")
    fig.suptitle("SIMULATION, post hoc, descriptive. Mean vs worst-10% trade-off. Sealed H1 "
                 "(tail vs front-loaded, margin -7 d) FAILED in both scripts 171 and 173.", fontsize=9)
    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_DIR / "tradeoff_mean_vs_cvar.png", dpi=150)
    (OUT_DIR / "tradeoff_points.json").write_text(json.dumps(
        {"script": "175_trackc_tradeoff_chart", "label": "SIMULATION; post hoc; descriptive; no new runs",
         "sources": {k: str(p.relative_to(PROJECT_ROOT)) for k, p in {**SEALED, "tune": TUNE}.items()},
         "sealed_points": points, "tuning_points": tpts}, indent=2))
    print("wrote", OUT_DIR)


if __name__ == "__main__":
    main()
