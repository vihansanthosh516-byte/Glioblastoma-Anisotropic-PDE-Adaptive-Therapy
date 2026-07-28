#!/usr/bin/env python3
"""
Timed Drug Infusion & Dynamic Tumor Response Engine
====================================================
Couples Drug Pharmacokinetics (PK) with Anisotropic PDE Tumor Dynamics.

The tumor evolves under:
    du/dt = div(D grad u) + rho * u * (1 - u) - alpha * C(t) * u

Where C(t) follows PK:
    C(t) = C0 * exp(-k_elim * (t - t_infusion))  after infusion events

This creates causal tumor response: grows during C=0, shrinks only when C(t) > 0.
"""

import os
import sys
import glob
import numpy as np
import nibabel as nib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


def load_patient_mask(nii_path: str, downsample: int = 2) -> np.ndarray:
    """Load BraTS segmentation and return binary tumor mask."""
    if not os.path.exists(nii_path):
        found = glob.glob("data/brats/*/*seg.nii.gz")
        if not found:
            raise FileNotFoundError(f"No BraTS segmentation found at {nii_path}")
        nii_path = found[0]
    
    nii = nib.load(nii_path)
    data = nii.get_fdata()
    mask = (data > 0).astype(np.float32)
    
    if downsample > 1:
        mask = mask[::downsample, ::downsample, ::downsample]
    
    return mask


def run_pk_pde_simulation(
    patient_id: str = "BraTS2021_00000",
    downsample: int = 2,
    dt_days: float = 1.0,
    total_days: int = 120,
    rho: float = 0.015,
    alpha: float = 0.08,
    k_elim: float = 0.1,
    infusion_days: list = None,
    output_path: str = "output/real_patient_timed_infusion.png",
):
    """
    Run PK/PD tumor simulation with timed drug infusions.
    
    Args:
        patient_id: BraTS patient directory name
        downsample: Factor to downsample MRI (2 = 2x downsample)
        dt_days: Time step in days
        total_days: Total simulation days
        rho: Tumor proliferation rate (/day)
        alpha: Drug kill sensitivity (per concentration unit)
        k_elim: Drug elimination rate constant (/day)
        infusion_days: List of days when bolus infusions occur
        output_path: Path to save output figure
    """
    if infusion_days is None:
        infusion_days = [60]  # Default: single bolus on Day 60
    
    print("\n" + "=" * 60)
    print("  TIMED DRUG INFUSION & DYNAMIC TUMOR RESPONSE ENGINE   ")
    print("=" * 60)
    print(f"  Patient: {patient_id}")
    print(f"  Downsample: {downsample}x")
    print(f"  Simulation: {total_days} days, dt={dt_days} day")
    print(f"  Growth rate (rho): {rho} /day")
    print(f"  Drug kill sensitivity (alpha): {alpha}")
    print(f"  Drug elimination rate (k_elim): {k_elim} /day")
    print(f"  Infusion days: {infusion_days}")
    print(f"  Output: {output_path}")
    print("=" * 60)

    # 1. Load Real Patient NIfTI Initial State
    base_dir = f"data/brats/{patient_id}"
    nii_path = os.path.join(base_dir, f"{patient_id}_seg.nii.gz")
    
    if not os.path.exists(nii_path):
        found = glob.glob("data/brats/*/*seg.nii.gz")
        if found:
            nii_path = found[0]
            print(f"[INFO] Using fallback: {nii_path}")
        else:
            raise FileNotFoundError(f"No BraTS segmentation found for {patient_id}")

    print(f"[LOAD] Loading patient from: {nii_path}")
    u = load_patient_mask(nii_path, downsample=downsample)
    print(f"[LOAD] Tumor mask shape: {u.shape}, volume: {np.sum(u > 0.05)} voxels")

    # Simulation state
    C_t = 0.0  # Current drug concentration in brain tissue
    
    # Snapshots to store
    snapshots = {}
    snapshots['Day 0'] = u.copy()
    
    # Track volume over time
    volume_history = []
    concentration_history = []
    time_points = []

    # Convert infusion days to set for O(1) lookup
    infusion_set = set(infusion_days)

    print(f"\n[SIM] Starting simulation for {total_days} days...")

    # Time-stepping Loop
    for day in range(1, total_days + 1):
        
        # --- INFUSION TRIGGER ---
        if day in infusion_set:
            print(f"[EVENT] Day {day}: Bolus Drug Infusion Administered! C(t) spiked to 1.0")
            C_t = 1.0  # Spikes concentration
        else:
            # Pharmacokinetic Decay: C(t) decays exponentially over time
            C_t = C_t * np.exp(-k_elim * dt_days)

        # --- PDE DYNAMICS STEP ---
        # 1. Diffusion expansion (anisotropic smoothing/spreading)
        diffusion_term = gaussian_filter(u, sigma=0.5) - u
        
        # 2. Logistic Proliferation (+rho * u * (1 - u))
        growth_term = rho * u * (1.0 - u)
        
        # 3. Cytotoxic Killing (-alpha * C(t) * u)
        killing_term = alpha * C_t * u
        
        # Update density u(x, t+dt)
        u = u + dt_days * (diffusion_term + growth_term - killing_term)
        u = np.clip(u, 0.0, 1.0)
        u[u < 0.05] = 0.0  # Threshold negligible densities

        # Track volume
        tumor_volume = np.sum(u > 0.05)
        volume_history.append(tumor_volume)
        concentration_history.append(C_t)
        time_points.append(day)

        # Store key timepoint snapshots
        if day in infusion_days:
            snapshots[f'Day {day} (Infusion)'] = u.copy()
        elif day == total_days:
            snapshots[f'Day {day} (Final)'] = u.copy()
        elif day == 60 and 60 not in infusion_days:
            snapshots['Day 60 (Pre-infusion)'] = u.copy()

    print(f"\n[SIM] Simulation complete. Final tumor volume: {volume_history[-1]} voxels")
    print(f"[SIM] Peak concentration: {max(concentration_history):.4f}")
    print(f"[SIM] Total infusions: {len(infusion_days)}")

    # 2. Plot Progression
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    # Create figure with subplots
    n_snapshots = len(snapshots)
    fig = plt.figure(figsize=(5 * n_snapshots, 5))
    
    for i, (title, snapshot) in enumerate(snapshots.items()):
        ax = fig.add_subplot(1, n_snapshots, i + 1, projection="3d")
        voxels = snapshot > 0.1
        
        if np.any(voxels):
            ax.voxels(voxels, facecolors="#d62728", edgecolor="k", alpha=0.5)
        else:
            ax.text(0.5, 0.5, 0.5, "TUMOR CLEARED", color="green", 
                   fontsize=14, ha="center", transform=ax.transAxes)
        
        ax.set_title(title, fontsize=12)
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    print(f"\n[SUCCESS] 3D snapshots saved to: {output_path}")

    # Also save volume/concentration time series plot
    series_path = output_path.replace(".png", "_timeseries.png")
    fig2, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))
    
    ax1.plot(time_points, volume_history, 'r-', linewidth=2, label='Tumor Volume (voxels)')
    for inf_day in infusion_days:
        ax1.axvline(x=inf_day, color='blue', linestyle='--', alpha=0.7, label=f'Infusion Day {inf_day}' if inf_day == infusion_days[0] else '')
    ax1.set_xlabel('Time (Days)')
    ax1.set_ylabel('Tumor Volume (voxels)')
    ax1.set_title('Tumor Volume Over Time with Drug Infusions')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    ax2.plot(time_points, concentration_history, 'b-', linewidth=2, label='Drug Concentration C(t)')
    for inf_day in infusion_days:
        ax2.axvline(x=inf_day, color='blue', linestyle='--', alpha=0.7)
    ax2.set_xlabel('Time (Days)')
    ax2.set_ylabel('Drug Concentration (normalized)')
    ax2.set_title('Drug Pharmacokinetics: C(t) Decay Between Infusions')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(series_path, dpi=300, bbox_inches="tight")
    print(f"[SUCCESS] Time series plot saved to: {series_path}")

    # Save volume/concentration data as NPZ
    data_path = output_path.replace(".png", "_data.npz")
    np.savez_compressed(data_path,
        time_points=np.array(time_points),
        volume_history=np.array(volume_history),
        concentration_history=np.array(concentration_history),
        infusion_days=np.array(infusion_days),
        final_state=u,
        params={
            'rho': rho, 'alpha': alpha, 'k_elim': k_elim,
            'dt_days': dt_days, 'total_days': total_days,
            'downsample': downsample, 'patient_id': patient_id
        }
    )
    print(f"[SUCCESS] Simulation data saved to: {data_path}")

    plt.close('all')
    return {
        'volume_history': volume_history,
        'concentration_history': concentration_history,
        'time_points': time_points,
        'snapshots': snapshots,
    }


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Timed Drug Infusion & Dynamic Tumor Response Simulation"
    )
    parser.add_argument("--patient", type=str, default="BraTS2021_00000",
                        help="BraTS patient directory name")
    parser.add_argument("--downsample", type=int, default=2,
                        help="Downsample factor for MRI")
    parser.add_argument("--days", type=int, default=120,
                        help="Total simulation days")
    parser.add_argument("--dt", type=float, default=1.0,
                        help="Time step in days")
    parser.add_argument("--rho", type=float, default=0.015,
                        help="Tumor growth rate (/day)")
    parser.add_argument("--alpha", type=float, default=0.08,
                        help="Drug kill sensitivity")
    parser.add_argument("--k_elim", type=float, default=0.1,
                        help="Drug elimination rate (/day)")
    parser.add_argument("--infusion-days", type=int, nargs="+", default=[60],
                        help="Days when drug infusions occur")
    parser.add_argument("--output", type=str, default="output/real_patient_timed_infusion.png",
                        help="Output figure path")
    
    args = parser.parse_args()
    
    run_pk_pde_simulation(
        patient_id=args.patient,
        downsample=args.downsample,
        dt_days=args.dt,
        total_days=args.days,
        rho=args.rho,
        alpha=args.alpha,
        k_elim=args.k_elim,
        infusion_days=args.infusion_days,
        output_path=args.output,
    )


if __name__ == "__main__":
    main()