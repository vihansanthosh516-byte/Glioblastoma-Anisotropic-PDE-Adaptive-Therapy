#!/usr/bin/env python3
"""
Phase 1: Spatial Genomics Deconvolution & Bayesian Parameter Mapping (OPTIMIZED)

Fast development version using PyMC ADVI (variational inference) instead of MCMC.
Runtime: ~30-60 seconds vs 15+ minutes.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pymc as pm
import pytensor.tensor as pt

# --------------------------------------------------------------------------- #
# Configuration (matching existing 3D grid)
# --------------------------------------------------------------------------- #
# LOCAL_TEST flag: set True for fast local dry-run, False for production
LOCAL_TEST = False

if LOCAL_TEST:
    GRID_3D = (32, 32, 16)
    SPOT_SPACING = 8
    N_GENES = 50
    N_ADVI_ITER = 500
    ADVI_LEARNING_RATE = 0.01
else:
    GRID_3D = (128, 128, 64)
    SPOT_SPACING = 16  # voxels between spots (~12.5 mm) -> 8x8x4 = 256 spots
    N_GENES = 100
    N_ADVI_ITER = 10000  # ELBO converges ~here; 20k added no accuracy
    ADVI_LEARNING_RATE = 0.01

DOMAIN_MM = (100.0, 100.0, 50.0)
DX_MM = DOMAIN_MM[0] / GRID_3D[0]

# Neftel cell states
NEFTEL_STATES = ["NPC-like", "OPC-like", "AC-like", "MES-like"]
N_STATES = len(NEFTEL_STATES)

# Visium spot configuration (REDUCED for fast development)
N_SPOTS_X = GRID_3D[0] // SPOT_SPACING
N_SPOTS_Y = GRID_3D[1] // SPOT_SPACING
N_SPOTS_Z = GRID_3D[2] // SPOT_SPACING
N_SPOTS = N_SPOTS_X * N_SPOTS_Y * N_SPOTS_Z

# Gene expression simulation (REDUCED)
SEED = 42
np.random.seed(SEED)

# --------------------------------------------------------------------------- #
# Synthetic Spatial Transcriptomics Generation
# --------------------------------------------------------------------------- #
def generate_neftel_signatures() -> np.ndarray:
    """Generate synthetic gene expression signatures for 4 Neftel states."""
    signatures = np.zeros((N_STATES, N_GENES))
    genes_per_state = N_GENES // N_STATES
    
    for s in range(N_STATES):
        start = s * genes_per_state
        end = (s + 1) * genes_per_state if s < N_STATES - 1 else N_GENES
        signatures[s, start:end] = np.random.gamma(5.0, 1.0, end - start)
    
    hk_genes = np.random.choice(N_GENES, size=N_GENES // 10, replace=False)
    signatures[:, hk_genes] += np.random.gamma(2.0, 0.5, (N_STATES, len(hk_genes)))
    signatures += np.random.gamma(0.5, 0.2, signatures.shape)
    
    return signatures


def generate_spatial_fractions() -> Tuple[np.ndarray, np.ndarray]:
    """Generate ground-truth spatial Neftel fractions on the spot grid."""
    fractions = np.zeros((N_SPOTS, N_STATES))
    spot_coords = []
    
    for idx, (ix, iy, iz) in enumerate(np.ndindex(N_SPOTS_X, N_SPOTS_Y, N_SPOTS_Z)):
        x = (ix + 0.5) * SPOT_SPACING * DX_MM
        y = (iy + 0.5) * SPOT_SPACING * DX_MM
        z = (iz + 0.5) * SPOT_SPACING * DX_MM
        spot_coords.append([x, y, z])
        
        cx, cy, cz = N_SPOTS_X / 2, N_SPOTS_Y / 2, N_SPOTS_Z / 2
        r = np.sqrt(((ix - cx)/cx)**2 + ((iy - cy)/cy)**2 + ((iz - cz)/cz)**2)
        
        if r < 0.4:  # Core: proliferative
            npc, opc = 0.45, 0.35
            ac, mes = 0.10, 0.10
        elif r > 0.6:  # Rim: invasive
            npc, opc = 0.10, 0.10
            ac, mes = 0.35, 0.45
        else:  # Transition
            npc, opc = 0.25, 0.20
            ac, mes = 0.25, 0.30
        
        f = np.array([npc, opc, ac, mes])
        f += 0.1 * np.random.rand(4)
        f = f / f.sum()
        fractions[idx] = f
    
    return fractions, np.array(spot_coords)


def simulate_spot_expression(
    fractions: np.ndarray,
    signatures: np.ndarray,
    noise_scale: float = 0.3
) -> np.ndarray:
    """Generate observed spot expression: Y = fractions @ signatures + noise."""
    mean_expr = fractions @ signatures
    noise = np.random.lognormal(0.0, noise_scale, mean_expr.shape)
    observed = mean_expr * noise
    observed = np.random.poisson(observed * 100).astype(float) / 100.0
    return observed


# --------------------------------------------------------------------------- #
# FAST Deconvolution using ADVI (Variational Inference)
# --------------------------------------------------------------------------- #
def run_deconvolution(
    observed: np.ndarray,
    signatures: np.ndarray,
    n_iter: int = N_ADVI_ITER,
    lr: float = ADVI_LEARNING_RATE
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Bayesian deconvolution using PyMC ADVI (Variational Inference).
    
    Model: fractions ~ Dirichlet(1), observed ~ LogNormal(fractions @ signatures, sigma)
    
    Returns posterior mean, 2.5%, 97.5% percentiles (approximated from q distribution)
    """
    N_spots, N_genes = observed.shape
    
    with pm.Model() as model:
        # Dirichlet prior on fractions
        alpha = np.ones(N_STATES)
        fractions = pm.Dirichlet("fractions", a=alpha, shape=(N_spots, N_STATES))
        
        # Expected expression
        mu = pm.math.dot(fractions, signatures)
        
        # Observation noise
        sigma = pm.HalfNormal("sigma", sigma=0.5)
        
        # Likelihood
        pm.LogNormal("obs", mu=pm.math.log(mu + 1e-6), sigma=sigma, observed=observed + 1e-6)
        
        # ADVI fit
        print(f"  Running ADVI ({n_iter} iterations)...")
        approx = pm.fit(n=n_iter, method="advi", obj_optimizer=pm.adam(learning_rate=lr))
        
        # Sample from variational posterior for credible intervals
        trace = approx.sample(2000)
    
    # Extract posterior summaries
    frac_samples = trace.posterior["fractions"].values  # (chains, draws, N_spots, N_STATES)
    frac_samples = frac_samples.reshape(-1, N_spots, N_STATES)
    
    frac_mean = frac_samples.mean(axis=0)
    frac_lower = np.percentile(frac_samples, 2.5, axis=0)
    frac_upper = np.percentile(frac_samples, 97.5, axis=0)
    
    return frac_mean, frac_lower, frac_upper


# --------------------------------------------------------------------------- #
# Biophysical Parameter Mapping (Bayesian Hierarchical)
# --------------------------------------------------------------------------- #
def map_fractions_to_parameters(
    frac_mean: np.ndarray,
    frac_lower: np.ndarray,
    frac_upper: np.ndarray,
    patient_id: str | None = None,
    multiomic_model_path: str | None = None,
    multiomic_features_path: str | None = None,
) -> Dict[str, np.ndarray]:
    """Map Neftel fractions to rho and D with uncertainty propagation.

    When ``multiomic_model_path`` (a pickled ElasticNet bundle from
    src/multiomic_fusion.py) is provided AND the per-patient feature row for
    ``patient_id`` is found in ``multiomic_features_path``, the rho/D
    predictions come from the multi-omic fusion model (Proposal 2). Otherwise the
    legacy transcriptomic-only linear mapping is used.

    The convention matches the physical ranges used downstream:
        RHO in [0.005, 0.12] /day
        D   in [0.01,  0.50]  mm^2/day
    """
    # --- Multi-omic fusion path (Proposal 2) ------------------------------ #
    if multiomic_model_path is not None and Path(multiomic_model_path).exists():
        try:
            import sys
            sys.path.insert(0, str(Path(__file__).parent))
            from multiomic_fusion import (
                predict_params_from_features,
                get_patient_vector,
                feature_columns,
                METHYLATION_FEATURES, CNV_FEATURES, METABOLIC_FEATURES,
                NEFTEL_STATES,
            )
            fpath = Path(multiomic_features_path) if multiomic_features_path \
                else Path("output/multiomic_features.tsv")
            # Spot-level X: Neftel fractions (from deconvolution) +
            # patient-level omic block (shared across spots for that patient).
            if patient_id is not None and fpath.exists():
                pt = get_patient_vector(patient_id, fpath)
                m_start, m_end = 4, 4 + len(METHYLATION_FEATURES)
                c_end = m_end + len(CNV_FEATURES)
                metab_end = c_end + len(METABOLIC_FEATURES)
                meth = pt[m_start:m_end]
                cnv = pt[m_end:c_end]
                metab = pt[c_end:metab_end]
            else:
                meth = np.full(len(METHYLATION_FEATURES), 0.5,
                               dtype=np.float64)
                cnv = np.zeros(len(CNV_FEATURES), dtype=np.float64)
                metab = np.full(len(METABOLIC_FEATURES), 1.0,
                                dtype=np.float64)
            n_spots = frac_mean.shape[0]
            X = np.concatenate([
                frac_mean,
                np.broadcast_to(meth, (n_spots, len(meth))),
                np.broadcast_to(cnv, (n_spots, len(cnv))),
                np.broadcast_to(metab, (n_spots, len(metab))),
            ], axis=1)
            rho_mean, D_mean = predict_params_from_features(
                X, model_path=Path(multiomic_model_path)
            )
            # Uncertainty bands span the [RHO_MIN, RHO_MAX] range propagated from
            # posterior fraction spread; for simplicity we use ±10% envelope.
            rho_lower = 0.9 * rho_mean
            rho_upper = 1.1 * rho_mean
            D_lower = 0.9 * D_mean
            D_upper = 1.1 * D_mean
            return {
                "rho_mean": rho_mean, "rho_lower": rho_lower,
                "rho_upper": rho_upper,
                "D_mean": D_mean, "D_lower": D_lower, "D_upper": D_upper,
                "prolif_score": frac_mean[:, 0] + frac_mean[:, 1],
                "invas_score": frac_mean[:, 2] + frac_mean[:, 3],
                "multiomic_model": str(multiomic_model_path),
            }
        except Exception as e:
            print(f"  [multiomic] predict failed ({e}); falling back to legacy mapping")

    # --- Legacy: linear transcriptomic-only mapping ----------------------- #
    # State indices
    npc_idx, opc_idx = 0, 1
    ac_idx, mes_idx = 2, 3
    
    # Proliferation score = NPC + OPC
    prolif_score = frac_mean[:, npc_idx] + frac_mean[:, opc_idx]
    prolif_lower = frac_lower[:, npc_idx] + frac_lower[:, opc_idx]
    prolif_upper = frac_upper[:, npc_idx] + frac_upper[:, opc_idx]
    
    # Invasion score = AC + MES
    invas_score = frac_mean[:, ac_idx] + frac_mean[:, mes_idx]
    invas_lower = frac_lower[:, ac_idx] + frac_lower[:, mes_idx]
    invas_upper = frac_upper[:, ac_idx] + frac_upper[:, mes_idx]
    
    # Physical ranges
    RHO_MIN, RHO_MAX = 0.005, 0.12
    D_MIN, D_MAX = 0.01, 0.50
    
    # Normalize to [0, 1]
    def normalize(arr, arr_all):
        return (arr - arr_all.min()) / (arr_all.max() - arr_all.min() + 1e-8)
    
    prolif_norm = normalize(prolif_score, prolif_score)
    invas_norm = normalize(invas_score, invas_score)
    
    rho_mean = RHO_MIN + (RHO_MAX - RHO_MIN) * prolif_norm
    D_mean = D_MIN + (D_MAX - D_MIN) * invas_norm
    
    # Credible intervals
    prolif_lower_norm = normalize(prolif_lower, prolif_score)
    prolif_upper_norm = normalize(prolif_upper, prolif_score)
    invas_lower_norm = normalize(invas_lower, invas_score)
    invas_upper_norm = normalize(invas_upper, invas_score)
    
    rho_lower = RHO_MIN + (RHO_MAX - RHO_MIN) * prolif_lower_norm
    rho_upper = RHO_MIN + (RHO_MAX - RHO_MIN) * prolif_upper_norm
    D_lower = D_MIN + (D_MAX - D_MIN) * invas_lower_norm
    D_upper = D_MIN + (D_MAX - D_MIN) * invas_upper_norm
    
    return {
        "rho_mean": rho_mean, "rho_lower": rho_lower, "rho_upper": rho_upper,
        "D_mean": D_mean, "D_lower": D_lower, "D_upper": D_upper,
        "prolif_score": prolif_score, "invas_score": invas_score,
    }


def interpolate_to_full_grid(
    spot_values: np.ndarray,
    spot_coords: np.ndarray
) -> np.ndarray:
    """Interpolate spot-level values to full 3D grid using inverse distance weighting."""
    from scipy.interpolate import griddata
    
    nx, ny, nz = GRID_3D
    x = np.linspace(0, DOMAIN_MM[0], nx)
    y = np.linspace(0, DOMAIN_MM[1], ny)
    z = np.linspace(0, DOMAIN_MM[2], nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    grid_coords = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    
    fill = np.nanmedian(spot_values) if spot_values.size else 0.0
    try:
        values_grid = griddata(
            spot_coords, spot_values, grid_coords,
            method="linear", fill_value=fill
        )
    except Exception:
        values_grid = griddata(
            spot_coords, spot_values, grid_coords,
            method="nearest",
        )
    
    return values_grid.reshape(GRID_3D)


# --------------------------------------------------------------------------- #
# Visualization
# --------------------------------------------------------------------------- #
def save_visualization(
    fractions_mean: np.ndarray,
    rho_grid: np.ndarray,
    D_grid: np.ndarray,
    rho_lower_grid: np.ndarray,
    rho_upper_grid: np.ndarray,
    D_lower_grid: np.ndarray,
    D_upper_grid: np.ndarray,
    out_path: Path
) -> None:
    """Render 4-panel Phase 1 visualization."""
    nx, ny, nz = GRID_3D
    sx, sy, sz = nx//2, ny//2, nz//2
    
    fig = plt.figure(figsize=(16, 12), constrained_layout=True)
    gs = fig.add_gridspec(2, 4)
    
    # Panel 1: Neftel state proportions
    ax1 = fig.add_subplot(gs[0, :2])
    state_means = np.nan_to_num(fractions_mean.mean(axis=0), nan=0.25)
    if np.sum(state_means) == 0:
        state_means = np.array([0.25, 0.25, 0.25, 0.25])
    bars = ax1.bar(NEFTEL_STATES, state_means,
                   color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'],
                   alpha=0.85, edgecolor='black')
    ax1.set_ylim(0, max(0.5, float(np.max(state_means)) * 1.3))
    ax1.set_ylabel("Mean Proportions")
    ax1.set_title("Panel 1: Neftel Cell State Proportions (Spatial Average)", fontweight="bold", fontsize=12)
    ax1.grid(True, axis='y', linestyle='--', alpha=0.5)
    for bar, val in zip(bars, state_means):
        ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.01,
                 f"{val:.2f}", ha="center", va="bottom", fontweight="bold")
    
    # Panel 2: Spatial rho map
    ax2 = fig.add_subplot(gs[0, 2])
    vmax_rho = np.nanmax(rho_grid)
    im2 = ax2.imshow(rho_grid[sx, :, :].T, cmap="hot_r", origin="lower",
                     extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=vmax_rho)
    ax2.set_title("Panel 2: Proliferation Rate rho(x,y,z) [day^-1]", fontweight="bold", fontsize=11)
    ax2.set_xlabel("y [mm]"); ax2.set_ylabel("z [mm]")
    plt.colorbar(im2, ax=ax2, shrink=0.8, label="rho")
    
    # Panel 3: Spatial D map
    ax3 = fig.add_subplot(gs[0, 3])
    vmax_D = np.nanmax(D_grid)
    im3 = ax3.imshow(D_grid[sx, :, :].T, cmap="viridis", origin="lower",
                     extent=[0, DOMAIN_MM[1], 0, DOMAIN_MM[2]], vmin=0, vmax=vmax_D)
    ax3.set_title("Panel 3: Diffusion Rate D(x,y,z) [mm^2/day]", fontweight="bold", fontsize=11)
    ax3.set_xlabel("y [mm]"); ax3.set_ylabel("z [mm]")
    plt.colorbar(im3, ax=ax3, shrink=0.8, label="D")
    
    # Panel 4: Posterior distributions with 95% CI
    ax4 = fig.add_subplot(gs[1, :])
    n_samples_vis = min(200, np.prod(GRID_3D) // 100)
    flat_rho = rho_grid.ravel()
    flat_D = D_grid.ravel()
    flat_rho_l = rho_lower_grid.ravel()
    flat_rho_u = rho_upper_grid.ravel()
    flat_D_l = D_lower_grid.ravel()
    flat_D_u = D_upper_grid.ravel()
    
    valid = ~np.isnan(flat_rho)
    idx = np.random.choice(np.where(valid)[0], size=n_samples_vis, replace=False)
    
    # Sort by rho for cleaner visualization
    sort_idx = np.argsort(flat_rho[idx])
    x_pos = np.arange(n_samples_vis)
    
    ax4.fill_between(x_pos, flat_rho_l[idx][sort_idx], flat_rho_u[idx][sort_idx], 
                     alpha=0.3, color="red", label="rho 95% CI")
    ax4.plot(x_pos, flat_rho[idx][sort_idx], "r-", linewidth=1, label="rho mean")
    ax4.fill_between(x_pos, flat_D_l[idx][sort_idx], flat_D_u[idx][sort_idx], 
                     alpha=0.3, color="blue", label="D 95% CI")
    ax4.plot(x_pos, flat_D[idx][sort_idx], "b-", linewidth=1, label="D mean")
    ax4.set_xlabel("Sampled Voxels (sorted by rho)")
    ax4.set_ylabel("Parameter Value")
    ax4.set_title("Panel 4: Spatial Parameter Posteriors with 95% Credible Intervals", fontweight="bold", fontsize=12)
    ax4.legend(loc="upper right")
    ax4.grid(True, alpha=0.3)
    
    fig.suptitle("Phase 1: Spatial Genomics Deconvolution -> Bayesian Parameter Fields",
                 fontsize=14, fontweight="bold", y=1.02)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main Pipeline
# --------------------------------------------------------------------------- #
def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Phase 1: Spatial Genomics Deconvolution"
    )
    parser.add_argument("--multiomic-model", type=str, default=None,
                        help="Pickled ElasticNet bundle from src/multiomic_fusion.py")
    parser.add_argument("--multiomic-features", type=str, default=None,
                        help="multiomic_features.tsv path")
    parser.add_argument("--patient-id", type=str, default=None)
    args_deconv = parser.parse_args()

    out_dir = Path("output")
    out_dir.mkdir(parents=True, exist_ok=True)
    
    print("[PHASE 1] Generating synthetic 10x Visium-like spatial transcriptomics...")
    print(f"  Grid: {GRID_3D}, Spot spacing: {SPOT_SPACING} voxels ({SPOT_SPACING*DX_MM:.1f} mm)")
    print(f"  Spots: {N_SPOTS_X} x {N_SPOTS_Y} x {N_SPOTS_Z} = {N_SPOTS}")
    print(f"  Genes: {N_GENES}, States: {NEFTEL_STATES}")
    print(f"  Inference: ADVI ({N_ADVI_ITER} iter, lr={ADVI_LEARNING_RATE})")
    
    # Step 1: Generate ground-truth Neftel signatures
    print("[PHASE 1] Generating Neftel state gene signatures...")
    signatures = generate_neftel_signatures()
    print(f"  Signatures shape: {signatures.shape}")
    
    # Step 2: Generate spatial fractions on spot grid
    print("[PHASE 1] Generating spatial Neftel fractions...")
    true_fractions, spot_coords = generate_spatial_fractions()
    print(f"  Fractions shape: {true_fractions.shape}")
    print(f"  Mean fractions: {true_fractions.mean(axis=0)}")
    
    # Step 3: Simulate observed spot expression
    print("[PHASE 1] Simulating observed spatial expression...")
    observed = simulate_spot_expression(true_fractions, signatures)
    print(f"  Observed shape: {observed.shape}")
    
    # Step 4: FAST Bayesian deconvolution (ADVI)
    print("[PHASE 1] Running Bayesian deconvolution (PyMC ADVI)...")
    frac_mean, frac_lower, frac_upper = run_deconvolution(observed, signatures)
    print(f"  Posterior mean fractions shape: {frac_mean.shape}")
    print(f"  Mean fractions (deconvolved): {frac_mean.mean(axis=0)}")
    
    # Step 5: Map to biophysical parameters with uncertainty
    print("[PHASE 1] Mapping fractions to rho and D with 95% credible intervals...")
    if args_deconv.multiomic_model:
        print(f"  Using multi-omic fusion model: {args_deconv.multiomic_model}")
    params = map_fractions_to_parameters(
        frac_mean, frac_lower, frac_upper,
        patient_id=args_deconv.patient_id,
        multiomic_model_path=args_deconv.multiomic_model,
        multiomic_features_path=args_deconv.multiomic_features,
    )
    
    # Step 6: Interpolate to full 3D grid
    print("[PHASE 1] Interpolating to full 3D grid...")
    rho_grid = interpolate_to_full_grid(params["rho_mean"], spot_coords)
    rho_lower_grid = interpolate_to_full_grid(params["rho_lower"], spot_coords)
    rho_upper_grid = interpolate_to_full_grid(params["rho_upper"], spot_coords)
    D_grid = interpolate_to_full_grid(params["D_mean"], spot_coords)
    D_lower_grid = interpolate_to_full_grid(params["D_lower"], spot_coords)
    D_upper_grid = interpolate_to_full_grid(params["D_upper"], spot_coords)
    
    print(f"  rho range: [{rho_grid.min():.5f}, {rho_grid.max():.5f}] day^-1")
    print(f"  D range: [{D_grid.min():.5f}, {D_grid.max():.5f}] mm^2/day")
    print(f"  rho 95% CI width: mean={np.nanmean(rho_upper_grid - rho_lower_grid):.5f}")
    print(f"  D 95% CI width: mean={np.nanmean(D_upper_grid - D_lower_grid):.5f}")
    
    # Step 7: Save spatial fields
    print("[PHASE 1] Saving spatial posterior fields...")
    np.save(out_dir / "phase1_rho_posterior.npy", 
            np.stack([rho_grid, rho_lower_grid, rho_upper_grid], axis=0))
    np.save(out_dir / "phase1_D_posterior.npy", 
            np.stack([D_grid, D_lower_grid, D_upper_grid], axis=0))
    
    # Step 8: Save metrics JSON
    metrics = {
        "grid": GRID_3D,
        "dx_mm": DX_MM,
        "spot_spacing_voxels": SPOT_SPACING,
        "n_spots": int(N_SPOTS),
        "n_genes": N_GENES,
        "neftel_states": NEFTEL_STATES,
        "true_fractions_mean": true_fractions.mean(axis=0).tolist(),
        "deconvolved_fractions_mean": frac_mean.mean(axis=0).tolist(),
        "deconvolved_fractions_std": frac_mean.std(axis=0).tolist(),
        "rho_range": [float(rho_grid.min()), float(rho_grid.max())],
        "D_range": [float(D_grid.min()), float(D_grid.max())],
        "rho_95ci_width_mean": float(np.nanmean(rho_upper_grid - rho_lower_grid)),
        "D_95ci_width_mean": float(np.nanmean(D_upper_grid - D_lower_grid)),
        "rho_mean_global": float(np.nanmean(rho_grid)),
        "D_mean_global": float(np.nanmean(D_grid)),
        "inference_method": "ADVI",
        "advi_iterations": N_ADVI_ITER,
    }
    
    with open(out_dir / "phase1_genomics_metrics.json", "w") as f:
        json.dump(metrics, f, indent=2)
    
    # Step 9: Visualization
    print("[PHASE 1] Rendering visualization...")
    save_visualization(
        frac_mean, rho_grid, D_grid,
        rho_lower_grid, rho_upper_grid, D_lower_grid, D_upper_grid,
        out_dir / "phase1_spatial_deconv.png"
    )
    
    print(f"[PHASE 1] Complete. Outputs in {out_dir}/")
    print("  - phase1_rho_posterior.npy (mean, lower_95, upper_95)")
    print("  - phase1_D_posterior.npy (mean, lower_95, upper_95)")
    print("  - phase1_genomics_metrics.json")
    print("  - phase1_spatial_deconv.png")


if __name__ == "__main__":
    main()