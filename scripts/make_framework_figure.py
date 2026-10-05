"""Framework schematic: three tracks, their data, methods, headline results, and limits.

Text is hand-written from vault/PAPER_FINDINGS_LEDGER.md. Run: python scripts/make_framework_figure.py
Writes paper/figures/fig0_framework.png
"""
from pathlib import Path as _Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUT = PROJECT_ROOT / "paper" / "figures" / "fig0_framework.png"

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, PALE = "#0b0b0b", "#52514e", "#f3f2ee"
plt.rcParams.update({"font.family": "DejaVu Sans", "savefig.dpi": 300})

TRACKS = [
    dict(c=BLUE, name="Track A", sub="Single cell to clinical  (01-41)",
         data="multiomic-gbm atlas\n140,355 cells\nTCGA-GBM, n = 518",
         meth="C-GAT zone classifier\nSPIB energy landscape\nTransfer-entropy network\nInvasion models, virtual\nknockout screen, Cox models",
         res="C-GAT 78.7% vs scVI 73.0%\nIndex-1 saddle (2 states)\nAge HR 1.031, p = 8.9e-16",
         lim="Drug screen is in silico only;\nno forecast of real scans"),
    dict(c=ORANGE, name="Track B", sub="PDE and adaptive therapy  (42-49, 72-81)",
         data="MU-Glioma-Post MRI\n154 fitted, 133 forecast\nUCSF-PDGM DTI, n = 62",
         meth="Anisotropic Fisher-Kolmogorov PDE\nStromal coupling\nAdaptive dosing vs MTD\nCore-target forecast\nSobol sensitivity",
         res="54/61 less resistance\nForecast +0.026 Dice vs no-change\nElongation 1.21 vs 1.03",
         lim="DTI orientation adds nothing;\nadaptive ends with more tumour"),
    dict(c=AQUA, name="Track C", sub="Digital twin and RL  (50-68, 75-83)",
         data="Fitted growth rates (real)\nAssumed kill rates\n3 synthetic sets + real_test (21)",
         meth="Inverse estimation, robust MPC\nPPO / DAgger at equal budget\nProbe-then-commit\nPaced heuristic\nResistance models",
         res="Conditioned policy 100% (day 90)\nPacing +2.1 d real, +52.5 d synthetic\nProbe-then-commit 100% to 62%",
         lim="A one-line rule matches RL;\nresistance gain fails on real data"),
]
X0, W, GAP = 0.03, 0.29, 0.045
ROWS = {"data": (0.70, 0.13), "meth": (0.46, 0.20), "res": (0.27, 0.15), "lim": (0.13, 0.11)}
ROW_LAB = {"data": "Data", "meth": "Methods", "res": "Headline results", "lim": "Limits"}


def box(ax, x, y, w, h, c, text, fc=PALE, ec=None, fs=6.8, bold=False, tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=0.012",
                                fc=fc, ec=ec or c, lw=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=tc,
            fontweight="bold" if bold else "normal", linespacing=1.35)


def arrow(ax, p, q, c=INK2, ls="-", lw=1.1, rad=0.0):
    ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, color=c, lw=lw, ls=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0))


def main():
    f, ax = plt.subplots(figsize=(8.2, 5.2))
    ax.set_xlim(-0.01, 1.02)
    ax.set_ylim(0, 1)
    ax.axis("off")
    xs = []
    for i, t in enumerate(TRACKS):
        x = X0 + i * (W + GAP)
        xs.append(x)
        ax.add_patch(FancyBboxPatch((x, 0.88), W, 0.095, boxstyle="round,pad=0.004,rounding_size=0.012", fc=t["c"], ec=t["c"]))
        ax.text(x + W / 2, 0.955, t["name"], ha="center", va="center", color="white", fontsize=9.5, fontweight="bold")
        ax.text(x + W / 2, 0.905, t["sub"], ha="center", va="center", color="white", fontsize=6.3)
        for k, key in enumerate(("data", "meth", "res", "lim")):
            y, h = ROWS[key]
            box(ax, x, y, W, h, t["c"], t[key], fs=6.7, fc="white" if key in ("res",) else PALE)
        for k in range(3):
            y0 = ROWS[("data", "meth", "res", "lim")[k]][0]
            y1 = ROWS[("data", "meth", "res", "lim")[k + 1]]
            arrow(ax, (x + W / 2, y0 - 0.002), (x + W / 2, y1[0] + y1[1] + 0.004), c=INK2, lw=0.9)
    # row labels on the far left margin
    for key, (y, h) in ROWS.items():
        ax.text(0.003, y + h / 2, ROW_LAB[key], rotation=90, ha="center", va="center", fontsize=6.3, color=INK2, fontweight="bold")
    # cross-track links
    ym = ROWS["meth"][0] + ROWS["meth"][1] / 2
    ax.text(xs[0] + W + GAP / 2, ym + 0.05, "no\nlink", ha="center", va="center", fontsize=5.4, color=INK2)
    arrow(ax, (xs[1] + W + 0.002, ym + 0.05), (xs[2] - 0.002, ym + 0.05), c=INK, ls="-", lw=1.5)
    ax.text(xs[1] + W + GAP / 2, ym + 0.075, "PDE\nsolver", ha="center", va="bottom", fontsize=5.2, color=INK)
    ax.text(0.5, 0.03, "Track A stands alone: the inflammation score in scripts 43-44 is 1.0 for every patient (no ID overlap with TCGA).\n"
            "Negative results reported in full: D_f invalid  |  MGMT does not predict time to progression\n"
            "forecast loses on shrinking tumours  |  resistance-driven adaptive gain fails on real data",
            ha="center", va="center", fontsize=6.2, color=INK2, linespacing=1.5)
    f.savefig(OUT, bbox_inches="tight", pad_inches=0.08, facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
