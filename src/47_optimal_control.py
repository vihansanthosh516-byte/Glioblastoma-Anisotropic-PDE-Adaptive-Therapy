"""Phase 3: MPC Optimal Control & Dual-Drug Protocol.
Evaluated on 8 real MU-Glioma patients selected from the 154-patient
cohort (positive rho, R² > 0.7, spread across rho range).
Per-patient rho drives the ODE; MPC optimizes dosing schedule.
"""
import json
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

warnings.filterwarnings("ignore")

from pathlib import Path as _Path
PROJECT_ROOT = _Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------------------- #
# Phase 1 physical constants + real MU-Glioma cohort
# --------------------------------------------------------------------------- #
from mu_glioma_loader import load_mu_glioma_params, real_cohort_stats

_params = load_mu_glioma_params()
_cohort = real_cohort_stats()

# Select 8 real patients: positive rho, R² > 0.7, spread across
# the rho range using linspace
_positive = {pid: p for pid, p in _params.items()
             if p["rho_per_day"] > 0 and p["r_squared"] > 0.7}
_sorted = sorted(_positive.items(), key=lambda x: x[1]["rho_per_day"])
if len(_sorted) < 8:
    raise RuntimeError(f"Only {len(_sorted)} eligible patients; need 8")
_indices = np.linspace(0, len(_sorted) - 1, 8).astype(int)
COHORT_PATIENTS = [_sorted[i][0] for i in _indices]

print(f"[Cohort] Selected {len(COHORT_PATIENTS)} real patients:")
for pid in COHORT_PATIENTS:
    p = _params[pid]
    print(f"  {pid}: rho={p['rho_per_day']:.6f}, V0={p['V0_mm3']:.0f}, R²={p['r_squared']:.3f}")

RHO_S = _cohort["rho_median"]                # real MU-Glioma median (E_MAX calibration)
RHO_R = RHO_S * 0.75                          # resistance fitness cost
K = 1.0
K_EL = np.log(2) / 0.075    # ~9.24 /day, TMZ
K_EL2 = np.log(2) / 0.15    # stromal inhibitor (~2x longer half-life)
C_PEAK = 10.0               # ug/mL TMZ
C2_PEAK = 5.0               # ug/mL stromal inhibitor
EC50 = 5.0                  # TMZ EC50 ug/mL
EC50_2 = 3.0                # stromal inhibitor EC50 ug/mL
HILL_COEFF = 2.0
E_MAX = RHO_S * 1000.0                        # calibrated to scripts 44 and 46
print(f"[Params] RHO_S = {RHO_S:.6f} /day, RHO_R = {RHO_R:.6f}, E_MAX = {E_MAX:.4f}")
print(f"[Params] cohort n={_cohort['n_total']}, rho_median={_cohort['rho_median']:.6f}")
GAMMA_R = 0.4              # secondary drug kill rate on resistant cells
MUTATION_RATE = 1e-5
INFUSION_RATE = C_PEAK

# MPC
MPC_HORIZON_DAYS = 14
DT = 1.0                     # 1-day steps for MPC (ODE integrated daily)
SIM_DAYS = 360
N_SIM_STEPS = int(SIM_DAYS / DT)
W_TUMOR = 1.0
W_DRUG_BASELINE = 0.1


# --------------------------------------------------------------------------- #
# Reduced ODE (1-day step)
# --------------------------------------------------------------------------- #
def step_reduced_ode(state: np.ndarray, a: float, a2: float,
                     dt: float = DT,
                     rho_s: float = RHO_S, rho_r: float = RHO_R) -> np.ndarray:
    """Single-day reduced ODE step.

    state = [M_s, M_r, C, C2]
    a in [0,1]: TMZ dose rate (fraction of C_PEAK bolus)
    a2 in [0,1]: stromal inhibitor dose rate

    dM_s/dt = rho_s * M_s * (1 - M_s - M_r) - kill_tmz * M_s
    dM_r/dt = rho_r * M_r * (1 - M_s - M_r) + mu * rho_s * M_s - kill_stromal * M_r
    dC/dt   = -k_el * C + a * infusion_rate
    dC2/dt  = -k_el2 * C2 + a2 * C2_peak
    """
    M_s, M_r, C, C2 = state
    total = M_s + M_r

    kill_tmz = E_MAX * (C ** HILL_COEFF) / (EC50 ** HILL_COEFF + C ** HILL_COEFF + 1e-12)
    # Secondary drug only kills resistant during TMZ holidays
    kill_stromal = 0.0 if a >= 0.01 else GAMMA_R * (C2 ** HILL_COEFF) / (
        EC50_2 ** HILL_COEFF + C2 ** HILL_COEFF + 1e-12)

    dMs = (rho_s * M_s * (1.0 - total) - kill_tmz * M_s) * dt
    dMr = (rho_r * M_r * (1.0 - total) + MUTATION_RATE * rho_s * M_s * 1e4
           - kill_stromal * M_r) * dt
    dC = (-K_EL * C + a * INFUSION_RATE) * dt
    dC2 = (-K_EL2 * C2 + a2 * C2_PEAK) * dt

    new = np.array(state) + np.array([dMs, dMr, dC, dC2])
    return np.maximum(new, 0.0)


# --------------------------------------------------------------------------- #
# MPC cost function: simulate 14-day horizon with control sequence
# --------------------------------------------------------------------------- #
def mpc_cost(control_seq: np.ndarray, state0: np.ndarray,
              w_tumor: float, w_drug: float,
              rho_s: float = RHO_S, rho_r: float = RHO_R) -> float:
    """Cost = expected sum over horizon over ensemble."""
    # +/- 15% uncertainty ensemble relative to patient rho
    ensemble = [
        (rho_s, rho_r),
        (rho_s * 1.15, rho_r * 1.15),
        (rho_s * 0.85, rho_r * 0.85)
    ]

    total_expected_cost = 0.0
    for rho_s_ens, rho_r_ens in ensemble:
        state = state0.copy()
        cost = 0.0
        for i in range(len(control_seq)):
            a = float(np.clip(control_seq[i], 0.0, 1.0))
            state = step_reduced_ode(state, a, 0.0, dt=DT, rho_s=rho_s_ens, rho_r=rho_r_ens)
            cost += w_tumor * (state[0] + state[1]) + w_drug * a * a
        total_expected_cost += cost

    return float(total_expected_cost / len(ensemble))


def solve_mpc_horizon(state0: np.ndarray, w_tumor: float, w_drug: float,
                      rho_s: float = RHO_S, rho_r: float = RHO_R
                      ) -> np.ndarray:
    """Solve 14-day MPC horizon, return optimal daily controls [0,1]^14."""
    res = minimize(
        mpc_cost,
        x0=np.full(MPC_HORIZON_DAYS, 0.5),
        args=(state0, w_tumor, w_drug, rho_s, rho_r),
        method="L-BFGS-B",
        bounds=[(0.0, 1.0)] * MPC_HORIZON_DAYS,
        options={"maxiter": 200, "ftol": 1e-6},
    )
    return np.clip(res.x, 0.0, 1.0)


# --------------------------------------------------------------------------- #
# Arm 1: MTD (continuous TMZ 5-on/23-off)
# --------------------------------------------------------------------------- #
def run_mtd(state0: np.ndarray,
            rho_s: float = RHO_S, rho_r: float = RHO_R) -> dict:
    state = state0.copy()
    n = N_SIM_STEPS
    Ms, Mr, C, C2, a_hist = (np.zeros(n) for _ in range(5))
    for step in range(n):
        t = step * DT
        day_in_cycle = int(t) % 28
        a = 1.0 if day_in_cycle < 5 else 0.0
        state = step_reduced_ode(state, a, 0.0, rho_s=rho_s, rho_r=rho_r)
        Ms[step], Mr[step], C[step], C2[step], a_hist[step] = (
            state[0], state[1], state[2], state[3], a)
    return {"M_s": Ms, "M_r": Mr, "C": C, "C2": C2, "a": a_hist, "a2": np.zeros(n)}


# --------------------------------------------------------------------------- #
# Arm 2: Single-agent MPC adaptive (TMZ only)
# --------------------------------------------------------------------------- #
def run_single_adaptive(state0: np.ndarray, w_tumor: float = W_TUMOR,
                         w_drug: float = W_DRUG_BASELINE,
                         rho_s: float = RHO_S, rho_r: float = RHO_R) -> dict:
    state = state0.copy()
    n = N_SIM_STEPS
    Ms, Mr, C, C2, a_hist = (np.zeros(n) for _ in range(5))
    daily_controls = None
    ctrl_idx = 0
    for step in range(n):
        if step % MPC_HORIZON_DAYS == 0:
            daily_controls = solve_mpc_horizon(state, w_tumor, w_drug,
                                               rho_s=rho_s, rho_r=rho_r)
            ctrl_idx = 0
        a = float(daily_controls[ctrl_idx])
        state = step_reduced_ode(state, a, 0.0, rho_s=rho_s, rho_r=rho_r)
        Ms[step], Mr[step], C[step], C2[step], a_hist[step] = (
            state[0], state[1], state[2], state[3], a)
        ctrl_idx += 1
    return {"M_s": Ms, "M_r": Mr, "C": C, "C2": C2, "a": a_hist, "a2": np.zeros(n)}


# --------------------------------------------------------------------------- #
# Arm 3: Dual-agent MPC adaptive (TMZ + stromal inhibitor on holidays)
# --------------------------------------------------------------------------- #
def run_dual_adaptive(state0: np.ndarray, w_tumor: float = W_TUMOR,
                       w_drug: float = W_DRUG_BASELINE,
                       rho_s: float = RHO_S, rho_r: float = RHO_R) -> dict:
    state = state0.copy()
    n = N_SIM_STEPS
    Ms, Mr, C, C2, a_hist, a2_hist = (np.zeros(n) for _ in range(6))
    daily_controls = None
    ctrl_idx = 0
    for step in range(n):
        if step % MPC_HORIZON_DAYS == 0:
            daily_controls = solve_mpc_horizon(state, w_tumor, w_drug,
                                               rho_s=rho_s, rho_r=rho_r)
            ctrl_idx = 0
        a = float(daily_controls[ctrl_idx])
        a2 = 1.0 if a < 0.01 else 0.0
        state = step_reduced_ode(state, a, a2, rho_s=rho_s, rho_r=rho_r)
        Ms[step], Mr[step], C[step], C2[step] = state[0], state[1], state[2], state[3]
        a_hist[step], a2_hist[step] = a, a2
        ctrl_idx += 1
    return {"M_s": Ms, "M_r": Mr, "C": C, "C2": C2, "a": a_hist, "a2": a2_hist}


# --------------------------------------------------------------------------- #
# Metrics
# --------------------------------------------------------------------------- #
def compute_ttp(total_mass_trajectory, nadir_fraction=1.2):
    """First time after nadir where mass exceeds 1.2x the minimum."""
    arr = np.asarray(total_mass_trajectory)
    nadir_idx = int(np.argmin(arr))
    threshold = arr[nadir_idx] * nadir_fraction
    after = arr[nadir_idx:]
    over = np.where(after > threshold)[0]
    if len(over) > 0:
        return float(nadir_idx + over[0])
    return float(len(arr))  # censored


def compute_auc(trajectory):
    """Total drug exposure: TMZ (C) + stromal inhibitor (C2)."""
    C = trajectory[:, 2]
    C2 = trajectory[:, 3]
    return float(np.trapezoid(C + C2, dx=DT))


def resistant_fraction(Ms: np.ndarray, Mr: np.ndarray) -> float:
    return float(Mr[-1] / (Ms[-1] + Mr[-1] + 1e-12))


def _trajectory_from_result(r: dict) -> np.ndarray:
    """Build [M_s, M_r, C, C2] trajectory for metric helpers."""
    return np.column_stack([r["M_s"], r["M_r"], r["C"], r["C2"]])


# --------------------------------------------------------------------------- #
# Patient initial conditions (real MU-Glioma V0)
# --------------------------------------------------------------------------- #
def get_patient_state(pid: str) -> np.ndarray:
    """Initial state from real MU-Glioma patient parameters.

    M_s0 scaled from V0_mm3: assumes V0_median ≈ 50000 mm³ maps to
    M_s0 = 0.10 (calibrated so mass range stays in [0.05, 0.30]).
    """
    p = _params[pid]
    M_s0 = float(np.clip(p["V0_mm3"] / 500000.0, 0.05, 0.30))
    M_r0 = 1e-4   # seed resistant clone
    return np.array([M_s0, M_r0, 0.0, 0.0])


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    print("=" * 70)
    print("Phase 3: MPC Optimal Control & Dual-Drug Protocol")
    print("=" * 70)
    print(f"MPC horizon: {MPC_HORIZON_DAYS} days | "
          f"baseline w_tumor={W_TUMOR}, w_drug={W_DRUG_BASELINE}")
    print(f"Trial horizon: {SIM_DAYS} days\n")

    all_results = []

    for pid in COHORT_PATIENTS:
        print(f"--- {pid} ---")
        state0 = get_patient_state(pid)
        rho_s_patient = _params[pid]["rho_per_day"]
        rho_r_patient = rho_s_patient * 0.75
        print(f"  rho_s={rho_s_patient:.6f}, rho_r={rho_r_patient:.6f}, "
              f"M_s0={state0[0]:.4f}")

        print("  Arm 1: MTD (continuous TMZ 5/23)...")
        r_mtd = run_mtd(state0, rho_s=rho_s_patient, rho_r=rho_r_patient)

        print("  Arm 2: Single-agent MPC adaptive...")
        r_sa = run_single_adaptive(state0, rho_s=rho_s_patient, rho_r=rho_r_patient)

        print("  Arm 3: Dual-agent MPC adaptive...")
        r_da = run_dual_adaptive(state0, rho_s=rho_s_patient, rho_r=rho_r_patient)

        ttp_m = compute_ttp(r_mtd["M_s"] + r_mtd["M_r"])
        ttp_s = compute_ttp(r_sa["M_s"] + r_sa["M_r"])
        ttp_d = compute_ttp(r_da["M_s"] + r_da["M_r"])
        auc_m = compute_auc(_trajectory_from_result(r_mtd))
        auc_s = compute_auc(_trajectory_from_result(r_sa))
        auc_d = compute_auc(_trajectory_from_result(r_da))
        rf_m = resistant_fraction(r_mtd["M_s"], r_mtd["M_r"])
        rf_s = resistant_fraction(r_sa["M_s"], r_sa["M_r"])
        rf_d = resistant_fraction(r_da["M_s"], r_da["M_r"])

        all_results.append({
            "patient_id": pid,
            "rho_per_day": rho_s_patient,
            "V0_mm3": _params[pid]["V0_mm3"],
            "r_squared": _params[pid]["r_squared"],
            "mtd": {"ttp_days": ttp_m, "auc": auc_m, "resistant_fraction": rf_m},
            "single_adaptive": {"ttp_days": ttp_s, "auc": auc_s, "resistant_fraction": rf_s},
            "dual_adaptive": {"ttp_days": ttp_d, "auc": auc_d, "resistant_fraction": rf_d},
        })
        print(f"    MTD    : TTP={ttp_m:6.1f}d  AUC={auc_m:8.1f}  Rfrac={rf_m:.3f}")
        print(f"    Single : TTP={ttp_s:6.1f}d  AUC={auc_s:8.1f}  Rfrac={rf_s:.3f}")
        print(f"    Dual   : TTP={ttp_d:6.1f}d  AUC={auc_d:8.1f}  Rfrac={rf_d:.3f}\n")

    # Save JSON
    out = {
        "parameters": {
            "MPC_HORIZON_DAYS": MPC_HORIZON_DAYS,
            "W_TUMOR": W_TUMOR,
            "W_DRUG_BASELINE": W_DRUG_BASELINE,
            "SIM_DAYS": SIM_DAYS,
            "DT": DT,
            "cohort_patients": COHORT_PATIENTS,
        },
        "patients": all_results,
    }
    json_path = OUTPUT_DIR / "dual_drug_comparison.json"
    with open(json_path, "w") as f:
        json.dump(out, f, indent=2)
    print(f"Results saved -> {json_path}")

    print("\n" + "=" * 70)
    print("3-Arm Benchmark Summary")
    print("=" * 70)
    hdr = (f"{'Patient (real)':<14} {'MTD TTP':>8} {'SA TTP':>8} {'DA TTP':>8} "
           f"{'MTD AUC':>10} {'SA AUC':>10} {'DA AUC':>10} "
           f"{'MTD Rf':>7} {'SA Rf':>7} {'DA Rf':>7}")
    print(hdr)
    print("-" * len(hdr))
    for e in all_results:
        m, s, d = e["mtd"], e["single_adaptive"], e["dual_adaptive"]
        print(f"{e['patient_id']:<14} {m['ttp_days']:>8.1f} {s['ttp_days']:>8.1f} "
              f"{d['ttp_days']:>8.1f} {m['auc']:>10.1f} {s['auc']:>10.1f} "
              f"{d['auc']:>10.1f} {m['resistant_fraction']:>7.3f} "
              f"{s['resistant_fraction']:>7.3f} {d['resistant_fraction']:>7.3f}")

    print("\nPhase 3 complete.")


if __name__ == "__main__":
    main()
