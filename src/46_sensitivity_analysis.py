"""Phase 2b: Global Sobol Sensitivity Analysis.

Reduced-ODE evaluator + SALib Sobol indices + publication-ready tornado plot.
N=500 base samples -> ~3500 model evaluations.
"""
import json
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from SALib.analyze import sobol as sobol_analyze
from SALib.sample import sobol as sobol_sample

warnings.filterwarnings("ignore")

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Load real MU-Glioma parameter distributions
from mu_glioma_loader import load_mu_glioma_params, real_cohort_stats

_params = load_mu_glioma_params()
_rhos = [p["rho_per_day"] for p in _params.values() if p["rho_per_day"] > 0]
_Ds = [p["D_mm2_per_day"] for p in _params.values() if p["D_mm2_per_day"]]

print(f"[Cohort] Loaded {len(_params)} real MU-Glioma patients")
print(f"[Cohort] rho_s range: [{min(_rhos):.5f}, {max(_rhos):.5f}] /day")
print(f"[Cohort] D_white range: [{min(_Ds):.7f}, {max(_Ds):.7f}] mm^2/day")

_cohort = real_cohort_stats()

# Physical constants (Phase 1, mm / days / ug/mL)
DT = 0.1
SIM_DAYS = 360
N_STEPS = int(SIM_DAYS / DT)

RHO_R_BASE = 0.015
K = 1.0
K_EL = np.log(2) / 0.075
C_PEAK = 10.0
HILL_COEFF = 2.0

# Same calibration as script 44: E_MAX_RATIO = 1000
E_MAX_RATIO = 1000.0
E_MAX = _cohort["rho_median"] * E_MAX_RATIO

print(f"[Params] E_MAX = {E_MAX:.4f} /day (from E_MAX_RATIO=1000, rho_median={_cohort['rho_median']:.6f})")

# TTP threshold: total cell mass fraction
TTP_FRACTION = 0.4

problem = {
    "num_vars": 5,
    "names": ["rho_s", "aniso_ratio", "mu_r", "EC50", "D_white"],
    "bounds": [
        [min(_rhos), max(_rhos)],       # rho_s from real MU-Glioma cohort
        [8.0, 12.0],                    # aniso_ratio: literature
        [8e-6, 1.2e-5],                 # mu_r: literature
        [4.0, 6.0],                     # EC50: literature (TMZ)
        [min(_Ds), max(_Ds)],           # D_white from real cohort
    ],
}


def reduced_ode(params: np.ndarray) -> float:
    """Reduced spatially-averaged ODE returning TTP (days)."""
    rho_s, aniso_ratio, mu_r, ec50, d_white = params
    rho_r = RHO_R_BASE
    mu = mu_r

    k_diff = 15.0
    k_aniso = 0.2
    eff_rho_s = rho_s * (1.0 + k_diff * d_white) * (1.0 + k_aniso * (aniso_ratio - 1.0))

    M_s = 0.05
    M_r = 1e-4
    C = 0.0

    days_on = 5
    cycle_days = 28

    for step in range(N_STEPS):
        t = step * DT
        day_in_cycle = int(t) % cycle_days

        if day_in_cycle < days_on:
            C = C_PEAK
        else:
            C *= np.exp(-K_EL * DT)

        kill = E_MAX * (C ** HILL_COEFF) / (ec50 ** HILL_COEFF + C ** HILL_COEFF + 1e-12)

        total = M_s + M_r
        dMs = (eff_rho_s * M_s * (1.0 - total) - kill * M_s) * DT
        dMr = (rho_r * M_r * (1.0 - total) + mu * 5e4 * eff_rho_s * M_s) * DT

        M_s = max(M_s + dMs, 0.0)
        M_r = max(M_r + dMr, 0.0)

        if M_s + M_r >= TTP_FRACTION:
            return t + DT

    return float(SIM_DAYS)


def run_simulation_and_get_ttp(params: np.ndarray) -> float:
    return reduced_ode(params)


def plot_tornado(Si: dict, path: Path) -> None:
    """Publication-ready tornado plot."""
    names = list(problem["names"])
    st_vals = np.asarray(Si["ST"])
    s1_vals = np.asarray(Si["S1"])
    st_conf = np.asarray(Si["ST_conf"])

    order = np.argsort(st_vals)[::-1]
    names_sorted = [names[i] for i in order]
    st_sorted = st_vals[order]
    s1_sorted = s1_vals[order]
    conf_sorted = st_conf[order]

    y = np.arange(len(names_sorted))
    fig, ax = plt.subplots(figsize=(8, 4.5))

    ax.barh(y, st_sorted, xerr=conf_sorted, color="steelblue", alpha=0.85,
            label="Total-effect ST", capsize=3, height=0.5)
    ax.scatter(s1_sorted, y, color="crimson", s=50, zorder=5,
               label="First-order S1")

    ax.set_yticks(y)
    ax.set_yticklabels(names_sorted, fontsize=10)
    ax.set_xlabel("Sobol index", fontsize=11)
    ax.set_title("Sobol Sensitivity: TTP variance decomposition", fontsize=12,
                 fontweight="bold")
    ax.legend(fontsize=9, loc="lower right")
    ax.set_xlim(-0.05, 1.05)
    ax.grid(alpha=0.3, axis="x")
    fig.tight_layout()
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[Tornado] Saved -> {path}")


def main() -> None:
    print("=" * 70)
    print("Phase 2b: Global Sobol Sensitivity Analysis")
    print("=" * 70)
    np.random.seed(42)

    N = 500
    n_evals = N * (2 * problem["num_vars"] + 2)
    print(f"Generating {N} Saltelli samples ({n_evals} evaluations)...")
    param_values = sobol_sample.sample(problem, N)
    print(f"  -> {len(param_values)} parameter combinations")

    print("Running reduced-ODE simulations...")
    Y = np.array([run_simulation_and_get_ttp(p) for p in param_values])
    print(f"  TTP range: [{Y.min():.1f}, {Y.max():.1f}] days, "
          f"mean={Y.mean():.1f} +/- {Y.std():.1f}")

    print("Computing Sobol indices...")
    Si = sobol_analyze.analyze(problem, Y, print_to_console=False)

    results = {
        "S1": {name: float(Si["S1"][i]) for i, name in enumerate(problem["names"])},
        "ST": {name: float(Si["ST"][i]) for i, name in enumerate(problem["names"])},
        "S1_conf": {name: float(Si["S1_conf"][i]) for i, name in enumerate(problem["names"])},
        "ST_conf": {name: float(Si["ST_conf"][i]) for i, name in enumerate(problem["names"])},
        "n_samples": N,
        "ttp_mean": float(Y.mean()),
        "ttp_std": float(Y.std()),
        "E_MAX": float(E_MAX),
        "E_MAX_RATIO": float(E_MAX_RATIO),
        "cohort": {
            "n_patients": len(_params),
            "rho_min": float(min(_rhos)),
            "rho_max": float(max(_rhos)),
            "D_min": float(min(_Ds)),
            "D_max": float(max(_Ds)),
        },
    }

    out_path = OUTPUT_DIR / "sobol_sensitivity_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"Results saved -> {out_path}")

    plot_tornado(Si, OUTPUT_DIR / "sobol_tornado_plot.png")

    print("\nPhase 2b complete.")


if __name__ == "__main__":
    main()