"""
Script 74 -- Anisotropic vs matched-isotropic tumor geometry (Track B negative 3).

Question: once the base=2 box-counting fix is in (script 42, commit b2a395a),
is the anisotropic fractal dimension D_f different from the isotropic one?

Why a new comparison instead of script 45's isotropic cache: script 45's
isotropic baseline uses its own solver settings (ISO_D=0.02, 500 steps of
dt=0.02, a 40x40 seed) and only the 8 legacy PAT_xxxx patients, so it differs
from script 42 in more than the tensor. Here the ONLY difference between the
two arms is the tensor orientation:

    aniso : script 42's per-patient tensor D(x)            (CohortSimulator)
    iso   : D_iso(x) = trace(D(x))/2 * I                   (same mean diffusivity)

Same patients (all 103 in spatial_recurrence_profiles.npz), same per-patient
rho field, same seed centre, same dt and step count, same estimator.

Endpoints (fixed before running):
  primary    D_f, script 42's box-counting estimator, script 42 horizon (150 d)
  secondary  elongation (sqrt of PCA eigenvalue ratio of the mask; 1 = round)
             tract alignment fraction (script 42's estimator)
             D_f at a 4x longer horizon (600 d), where the blob spans more box sizes
Test: paired Wilcoxon signed-rank (two-sided) and paired t-test, aniso - iso.

Output: output/fractal_aniso_vs_iso.json
"""
from __future__ import annotations

import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUT_JSON = OUTPUT_DIR / "fractal_aniso_vs_iso.json"

_spec = spec_from_file_location("s42", PROJECT_ROOT / "src" / "42_anisotropic_pde.py")
s42 = module_from_spec(_spec)
_spec.loader.exec_module(s42)

THRESHOLD = 0.1
LONG_FACTOR = 4


def elongation(u: np.ndarray, threshold: float = THRESHOLD) -> float:
    ys, xs = np.nonzero(u > threshold)
    if len(xs) < 3:
        return float("nan")
    ev = np.linalg.eigvalsh(np.cov(np.vstack([xs, ys]).astype(float)))
    return float(np.sqrt(ev[1] / max(ev[0], 1e-12)))


def simulate(D_xx, D_xy, D_yy, rho, center, n_steps: int, dt: float) -> np.ndarray:
    solver = s42.AnisotropicFKSolver(D_xx=D_xx, D_xy=D_xy, D_yy=D_yy, rho=rho, dt=dt)
    u = s42.AnisotropicFKSolver.initial_gaussian_seed(
        (solver.H, solver.W), center=center, sigma=3.0, amplitude=0.8)
    for _ in range(n_steps):
        u = solver.step(u)
    return u


def paired(a, b) -> dict:
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    a, b = a[ok], b[ok]
    d = a - b
    out = {"n": int(ok.sum()), "mean_aniso": float(a.mean()), "mean_iso": float(b.mean()),
           "mean_diff": float(d.mean()), "median_diff": float(np.median(d)),
           "sd_diff": float(d.std(ddof=1)), "n_aniso_greater": int((d > 0).sum()),
           "n_equal": int((d == 0).sum())}
    out["cohens_dz"] = float(d.mean() / d.std(ddof=1)) if d.std(ddof=1) > 0 else float("nan")
    out["wilcoxon_p"] = float(stats.wilcoxon(a, b).pvalue) if np.any(d != 0) else 1.0
    out["paired_t_p"] = float(stats.ttest_rel(a, b).pvalue) if d.std(ddof=1) > 0 else 1.0
    return out


def main() -> None:
    builder = s42.TensorFieldBuilder(
        grid_size=s42.GRID_SIZE, d_parallel=s42.D_PARALLEL_DEFAULT,
        d_perpendicular=s42.D_PERPENDICULAR_DEFAULT, d_base=s42.D_BASE,
        tract_angle_deg=45.0, tract_width=15,
        atlas_mode="atlas", atlas_lam1=0.013, atlas_lam2=0.0013)
    builder.build_tract_mask()
    builder.build_orientation_field()
    builder.build_tensor_field()
    mapper = s42.PatientParameterMapper(cohort_npz=OUTPUT_DIR / "spatial_recurrence_profiles.npz",
                                        zone_csv_root=OUTPUT_DIR)
    mapper.load()
    sim = s42.CohortSimulator(mapper=mapper, n_steps=s42.N_PATIENT_STEPS,
                              save_interval=s42.PATIENT_SAVE_INTERVAL, dt=s42.DT_DEFAULT)
    viz = s42.AnisotropicVisualizer(builder)
    theta = builder.theta_field

    rows = []
    pids = [str(p) for p in mapper.data["patient_ids"]]
    for i, pid in enumerate(pids):
        p = mapper.build_patient_tensor_and_rho(pid, builder)
        md = 0.5 * (p["D_xx"] + p["D_yy"])
        zero = np.zeros_like(md)
        center = sim.patient_seed_center(pid, N=md.shape[0])
        row = {"patient_id": pid}
        for horizon, n_steps in (("h150", s42.N_PATIENT_STEPS), ("h600", LONG_FACTOR * s42.N_PATIENT_STEPS)):
            ua = simulate(p["D_xx"], p["D_xy"], p["D_yy"], p["rho"], center, n_steps, s42.DT_DEFAULT)
            ui = simulate(md, zero, md, p["rho"], center, n_steps, s42.DT_DEFAULT)
            for arm, u in (("aniso", ua), ("iso", ui)):
                row[f"{horizon}_{arm}_fd"] = viz.fractal_dimension(u, THRESHOLD)
                row[f"{horizon}_{arm}_elong"] = elongation(u)
                row[f"{horizon}_{arm}_align"] = viz.orientation_alignment(u, theta)
                row[f"{horizon}_{arm}_area_px"] = int((u > THRESHOLD).sum())
        rows.append(row)
        print(f"[{i+1}/{len(pids)}] {pid}: D_f aniso {row['h150_aniso_fd']:.3f} iso {row['h150_iso_fd']:.3f} | "
              f"elong {row['h150_aniso_elong']:.2f}/{row['h150_iso_elong']:.2f} | "
              f"600d D_f {row['h600_aniso_fd']:.3f}/{row['h600_iso_fd']:.3f}")

    col = lambda k: [r[k] for r in rows]  # noqa: E731
    res = {
        "script": "74_fractal_aniso_vs_iso",
        "design": "paired, same patient/rho/seed/dt/steps; iso tensor = trace(D)/2 * I",
        "n_patients": len(rows),
        "threshold": THRESHOLD,
        "primary_fd_150d": paired(col("h150_aniso_fd"), col("h150_iso_fd")),
        "secondary": {
            "elongation_150d": paired(col("h150_aniso_elong"), col("h150_iso_elong")),
            "tract_alignment_150d": paired(col("h150_aniso_align"), col("h150_iso_align")),
            "fd_600d": paired(col("h600_aniso_fd"), col("h600_iso_fd")),
            "elongation_600d": paired(col("h600_aniso_elong"), col("h600_iso_elong")),
        },
        "fd_in_1_2_range": {
            k: int(sum(1.0 <= r[k] <= 2.0 for r in rows))
            for k in ("h150_aniso_fd", "h150_iso_fd", "h600_aniso_fd", "h600_iso_fd")},
        "median_area_px": {k: float(np.median(col(k))) for k in
                           ("h150_aniso_area_px", "h150_iso_area_px", "h600_aniso_area_px", "h600_iso_area_px")},
        "script45_iso_cache_fd_for_reference": json.loads(
            (OUTPUT_DIR / "isotropic_baseline_metrics.json").read_text())["patients"][0]["fractal_dimension"],
        "patients": rows,
    }
    OUT_JSON.write_text(json.dumps(res, indent=2))
    print(json.dumps({k: v for k, v in res.items() if k != "patients"}, indent=2))
    print(f"[saved] {OUT_JSON}")


if __name__ == "__main__":
    sys.exit(main())
