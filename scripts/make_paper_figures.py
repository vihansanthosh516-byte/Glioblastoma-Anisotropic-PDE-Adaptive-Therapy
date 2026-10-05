"""Paper figures, all drawn from JSON/CSV in output/. No number is typed by hand.

Writes paper/figures/fig{1,2,3}_*.png (multi-panel) and one PNG per panel.
Colour: categorical slots 1-3 of the reference palette (validated all-pairs), plus neutral grey.
Run: python scripts/make_paper_figures.py
"""
import json
from pathlib import Path as _Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "output"
FIG = PROJECT_ROOT / "paper" / "figures"
FIG.mkdir(parents=True, exist_ok=True)

BLUE, ORANGE, AQUA, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8b8a85"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e4e3df"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.titleweight": "bold",
    "axes.labelsize": 8, "axes.edgecolor": INK2, "axes.labelcolor": INK2, "xtick.color": INK2,
    "ytick.color": INK2, "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.dpi": 100, "savefig.dpi": 300, "legend.frameon": False,
})


def J(p):
    return json.loads((OUT / p).read_text())


def tag(ax, letter):
    ax.text(-0.02, 1.08, letter, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


# ---------------------------------------------------------------- Track A
def panel_a1(ax):
    gat = J("cgat/gat_metrics.json")["accuracy"]
    scvi, nmf = J("scvi_metrics.json")["accuracy"], J("nmf_metrics.json")["accuracy"]
    bench = pd.read_csv(OUT / "benchmark_comparison.tsv", sep="\t").set_index("Method")["Accuracy"]
    rows = [("C-GAT", gat), ("scVI", scvi),
            ("Random forest", bench["Random Forest (Classical)"]),
            ("Logistic regression", bench["Logistic Regression (Classical)"]), ("NMF", nmf),
            ("Hybrid transformer", bench["Hybrid (LR-prior Transformer)"]),
            ("Transformer", bench["Transformer (Deep Learning)"])]
    rows.sort(key=lambda r: r[1])
    y = np.arange(len(rows))
    ax.barh(y, [r[1] * 100 for r in rows], color=[BLUE if r[0] == "C-GAT" else GREY for r in rows], height=0.62)
    for yi, r in zip(y, rows):
        ax.text(r[1] * 100 + 0.8, yi, f"{r[1] * 100:.1f}", va="center", fontsize=7.5)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows])
    ax.set_xlim(0, 92)
    ax.set_xlabel("Test accuracy, 3 tumour zones (%)")
    ax.set_title("Zone classification")
    ax.grid(axis="y", visible=False)


def panel_a2(ax):
    m = J("spib_saddle_point_metrics.json")
    att = m["attractors"]
    names = ["Periphery", "Core", "Healthy"]
    keys = ["Periphery_Attractor", "Core_Attractor", "Healthy_Attractor"]
    e = [att[k]["energy"] for k in keys] + [m["saddle"]["energy"]]
    x = np.arange(4)
    ax.bar(x, e, color=[BLUE, BLUE, BLUE, ORANGE], width=0.6)
    for xi, v in zip(x, e):
        ax.text(xi, v + 0.04, f"{v:.3f}", ha="center", fontsize=7.5)
    s = m["saddle"]
    ax.text(2.55, 1.2, f"saddle: index-1\neigenvalues +{s['positive_eigenvalues']} / −{s['negative_eigenvalues']}",
            ha="right", va="center", fontsize=7.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(names + ["Transition\nsaddle"])
    ax.set_ylim(0, 2.35)
    ax.set_ylabel("Landscape energy (a.u.)")
    ax.set_title("SPIB energy landscape")
    ax.grid(axis="x", visible=False)


def panel_a3(ax):
    s = J("survival_stats_summary.json")["univariate"]
    km = J("tcga_km_summary.json")
    lab = {"age_at_diagnosis": "Age (per year)", "gender": "Gender (M vs F)", "molecular_subtype": "Molecular subtype"}
    y = np.arange(len(s))[::-1]
    for yi, r in zip(y, s):
        sig = r["cox_p"] < 0.05
        c = BLUE if sig else GREY
        ax.plot([r["cox_ci_lower"], r["cox_ci_upper"]], [yi, yi], color=c, lw=1.8)
        ax.plot(r["cox_hr"], yi, "o", color=c, ms=6, mfc=c, mec="white", mew=1)
        ax.text(1.75, yi, f"HR {r['cox_hr']:.2f}, p = {r['cox_p']:.2g}", va="center", fontsize=7.5)
    ax.axvline(1, color=INK2, lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([lab[r["covariate"]] for r in s])
    ax.set_xscale("log")
    ax.set_xlim(0.7, 3.6)
    ax.minorticks_off()
    ax.set_xticks([0.8, 1, 1.25, 1.5])
    ax.set_xticklabels(["0.8", "1", "1.25", "1.5"])
    ax.set_xlabel("Cox hazard ratio (95% CI)\n" + f"n = {km['n_analysed']}, {km['n_events']} events\nKM median OS "
                  f"{km['km_median_days']:.0f} d ({km['km_median_ci95_days'][0]:.0f}–{km['km_median_ci95_days'][1]:.0f})")
    ax.set_title("TCGA-GBM survival")
    ax.grid(axis="y", visible=False)


# ---------------------------------------------------------------- Track B
def _adaptive():
    a = J("adaptive_geometry_metrics.json")
    d = pd.DataFrame(a)
    d["earlier"] = d["ttp_adaptive"] < d["ttp_mtd"]
    d["mass_ratio"] = d["final_mass_adaptive"] / d["final_mass_mtd"]
    return d


def panel_b1a(ax):
    d = _adaptive()
    rng = np.random.default_rng(0)
    xj = d["resistant_fraction_mtd"] + rng.normal(0, 0.004, len(d))
    ax.plot([0, 1.02], [0, 1.02], color=INK2, lw=0.8, ls="--")
    ax.scatter(xj[~d.earlier], d.resistant_fraction_adaptive[~d.earlier], s=22, color=BLUE, edgecolor="white", lw=0.6, label="TTP not earlier")
    ax.scatter(xj[d.earlier], d.resistant_fraction_adaptive[d.earlier], s=22, color=ORANGE, edgecolor="white", lw=0.6, label="Adaptive progressed earlier")
    below = int((d.resistant_fraction_mtd > d.resistant_fraction_adaptive).sum())
    above = int((d.resistant_fraction_adaptive > d.resistant_fraction_mtd).sum())
    ax.text(0.03, 0.97, f"{below}/{len(d)} below diagonal\n{above}/{len(d)} above", transform=ax.transAxes, va="top", fontsize=7.5)
    ax.set_xlim(0, 1.03)
    ax.set_ylim(0, 1.03)
    ax.set_xlabel("Resistant fraction, MTD")
    ax.set_ylabel("Resistant fraction, adaptive")
    ax.set_title("Resistance selection (n = 61)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 0.83), fontsize=7, handletextpad=0.3)


def panel_b1b(ax):
    d = _adaptive()
    ax.axhline(1, color=INK2, lw=0.8, ls="--")
    ax.axvline(0.02, color=GREY, lw=0.8, ls=":")
    ax.scatter(d.rho_per_day[~d.earlier], d.mass_ratio[~d.earlier], s=22, color=BLUE, edgecolor="white", lw=0.6)
    ax.scatter(d.rho_per_day[d.earlier], d.mass_ratio[d.earlier], s=22, color=ORANGE, edgecolor="white", lw=0.6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_yticks([0.8, 1, 2, 4])
    ax.set_yticklabels(["0.8", "1", "2", "4"])
    ax.minorticks_off()
    ax.text(0.0215, ax.get_ylim()[0] * 1.2, "rho = 0.02 /day", fontsize=7, color=INK2, rotation=90, va="bottom")
    ax.set_xlabel("Fitted growth rate rho (per day)")
    ax.set_ylabel("Final tumour mass, adaptive / MTD")
    ax.set_title("Final tumour mass ratio")


def panel_b2(ax):
    r = J("forecast_labels_core/results.json")
    groups = [("All\n(n = %d)" % r["n_scored"], r["dice_out_of_fold_all"]),
              ("Grew\n(n = %d)" % r["n_grew"], r["descriptive_grew"]["dice"]),
              ("Shrank / same\n(n = %d)" % r["n_shrank_or_same"], r["descriptive_shrank_or_same"]["dice"])]
    arms = [("no_change", "No-change", GREY), ("iso_homog", "Isotropic PDE", AQUA), ("anisotropic", "Anisotropic PDE", BLUE)]
    w = 0.25
    x = np.arange(len(groups))
    for i, (k, lab, c) in enumerate(arms):
        m = np.array([g[1][k]["mean"] for g in groups])
        lo = np.array([g[1][k]["mean_ci95"][0] for g in groups])
        hi = np.array([g[1][k]["mean_ci95"][1] for g in groups])
        ax.bar(x + (i - 1) * w, m, w * 0.92, color=c, label=lab, yerr=[m - lo, hi - m],
               error_kw=dict(ecolor=INK2, lw=0.8, capsize=2))
    t = r["primary_tests_all"]["anisotropic_vs_no_change"]
    ax.text(0.02, 0.86, f"All: +{t['mean_diff']:.3f} Dice vs no-change, Holm p = {t['wilcoxon_p_holm']:.4f}",
            transform=ax.transAxes, va="top", fontsize=7.5)
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in groups])
    ax.set_ylim(0, 0.55)
    ax.set_ylabel("Out-of-fold Dice, tumour core")
    ax.set_title("Forecast on the cellular core")
    ax.legend(loc="upper left", fontsize=7, ncol=3, columnspacing=0.8, handlelength=1, bbox_to_anchor=(0.0, 1.0))
    ax.grid(axis="x", visible=False)


def panel_b3(ax):
    f = J("fractal_aniso_vs_iso.json")
    p = pd.DataFrame(f["patients"])
    for a, b in zip(p.h150_iso_elong, p.h150_aniso_elong):
        ax.plot([0, 1], [a, b], color=GREY, alpha=0.35, lw=0.7)
    ax.plot([0, 1], [p.h150_iso_elong.mean(), p.h150_aniso_elong.mean()], color=BLUE, lw=2.4, marker="o", ms=6, mec="white")
    e = f["secondary"]["elongation_150d"]
    ax.text(0.5, 0.97, f"{e['n_aniso_greater']}/{e['n']} patients higher, p = {e['wilcoxon_p']:.1e}",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Isotropic", "Anisotropic"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylim(0.95, 1.7)
    ax.set_ylabel("Tumour elongation at 150 d")
    ax.set_title("Simulated shape")
    ax.grid(axis="x", visible=False)


# ---------------------------------------------------------------- Track C
SETS = [("cohort64", "cohort64\n(synthetic)"), ("lhs60", "lhs60\n(synthetic)"), ("synth_test", "synth_test\n(synthetic)"),
        ("real_test", "real_test\n(real params)")]


def panel_c1(ax):
    e = J("rl_kill_conditioned/evaluate.json")["sets"]
    arms = [("ppo66", "Blind PPO", GREY), ("ppo_cond_ft", "Kill-conditioned PPO", BLUE), ("efficiency_rule", "One-line rule", AQUA)]
    w = 0.26
    x = np.arange(len(SETS))
    for i, (k, lab, c) in enumerate(arms):
        v = [e[s]["arms"][k]["win_rate_vs_stupp_pct"] for s, _ in SETS]
        ax.bar(x + (i - 1) * w, v, w * 0.92, color=c, label=lab)
        for xi, vi in zip(x + (i - 1) * w, v):
            ax.text(xi, vi + 1.5, f"{vi:.0f}", ha="center", fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels([l for _, l in SETS], fontsize=7)
    ax.set_ylim(0, 132)
    ax.set_ylabel("Win rate vs Stupp, day-90 volume (%)")
    ax.set_title("Conditioning on kill rates")
    ax.legend(loc="upper left", fontsize=7, ncol=3, columnspacing=0.8, handlelength=1, bbox_to_anchor=(0, 1.0))
    ax.grid(axis="x", visible=False)


def panel_c2(ax):
    arms = J("probe_paced_policies/evaluate.json")["sets"]["cohort64"]["arms"]
    sig = [0.0, 0.02, 0.05, 0.10]

    def series(prefix):
        m, lo, hi = [], [], []
        for s in sig:
            seeds = [0] if s == 0 else [0, 1, 2]
            v = [arms[f"{prefix}_L1_sig{s:.2f}_s{k}"]["day90_win_rate_vs_stupp_pct"] for k in seeds]
            m.append(np.mean(v)); lo.append(min(v)); hi.append(max(v))
        return np.array(m), np.array(lo), np.array(hi)
    blind = arms["ppo66_blind"]["day90_win_rate_vs_stupp_pct"]
    ax.axhline(blind, color=GREY, lw=1.2, ls="--")
    ax.text(0.101, blind + 2.5, f"blind PPO {blind:.0f}%", ha="right", fontsize=7, color=INK2)
    for prefix, lab, c in [("probe_rule", "Probe + rule", BLUE), ("probe_dagger_cond", "Probe + conditioned net", AQUA)]:
        m, lo, hi = series(prefix)
        ax.errorbar(sig, m, yerr=[m - lo, hi - m], color=c, marker="o", ms=5, lw=1.8, capsize=2.5, label=lab, mec="white")
    ax.set_xlim(-0.006, 0.106)
    ax.set_ylim(0, 108)
    ax.set_xlabel("Noise in log-volume (sigma)\n3-seed mean, bars = seed range (sigma 0: one seed)")
    ax.set_ylabel("Win rate vs Stupp, day 90 (%)")
    ax.set_title("Probe-then-commit (cohort64)")
    ax.legend(loc="upper right", fontsize=7)


def panel_c3(ax):
    sets = J("probe_paced_policies/evaluate.json")["sets"]
    y = np.arange(len(SETS))[::-1]
    for yi, (s, lab) in zip(y, SETS):
        a = sets[s]["arms"]["paced_heuristic59"]
        c = ORANGE if s == "real_test" else BLUE
        lo, hi = a["ttp_minus_stupp_days_ci95"]
        m = a["ttp_minus_stupp_days_mean"]
        ax.plot([lo, hi], [yi, yi], color=c, lw=2)
        ax.plot(m, yi, "o", color=c, ms=7, mec="white", mew=1)
        L, E, S = a["ttp_n_longer_equal_shorter"]
        ax.text(hi + 3, yi, f"+{m:.1f} d   ({L} longer / {E} equal / {S} shorter)", va="center", fontsize=7)
    ax.axvline(0, color=INK2, lw=0.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels([l for _, l in SETS], fontsize=7)
    ax.set_xlim(-5, 150)
    ax.set_xlabel("Time to progression minus Stupp (days, 95% CI)")
    ax.set_title("Pacing at equal drug budget")
    ax.grid(axis="y", visible=False)


# ---------------------------------------------------------------- assemble
def save_panels(name, fn):
    f, ax = plt.subplots(figsize=(3.6, 2.9), constrained_layout=True)
    fn(ax)
    f.savefig(FIG / f"{name}.png")
    plt.close(f)


def combo(fname, panels, shape, size):
    f, axes = plt.subplots(*shape, figsize=size, constrained_layout=True)
    for ax, (letter, fn) in zip(np.atleast_1d(axes).ravel(), panels):
        fn(ax)
        tag(ax, letter)
    f.savefig(FIG / fname)
    plt.close(f)


def main():
    p1 = [("A", panel_a1), ("B", panel_a2), ("C", panel_a3)]
    p2 = [("A", panel_b1a), ("B", panel_b1b), ("C", panel_b2), ("D", panel_b3)]
    p3 = [("A", panel_c1), ("B", panel_c2), ("C", panel_c3)]
    combo("fig1_track_a.png", p1, (1, 3), (10.5, 3.3))
    combo("fig2_track_b.png", p2, (2, 2), (7.6, 6.2))
    combo("fig3_track_c.png", p3, (1, 3), (11.5, 3.4))
    for name, fn in [("fig1a_classification", panel_a1), ("fig1b_spib_saddle", panel_a2), ("fig1c_tcga_survival", panel_a3),
                     ("fig2a_resistance", panel_b1a), ("fig2b_final_mass", panel_b1b), ("fig2c_forecast", panel_b2),
                     ("fig2d_elongation", panel_b3), ("fig3a_conditioned_policy", panel_c1),
                     ("fig3b_probe_commit", panel_c2), ("fig3c_paced_ttp", panel_c3)]:
        save_panels(name, fn)
    print("wrote", len(list(FIG.glob("*.png"))), "figures to", FIG)


if __name__ == "__main__":
    main()
