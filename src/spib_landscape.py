#!/usr/bin/env python3
"""
Phase 3: Construct energy landscape in SPIB 3D space.
"""

import numpy as np
import os
from scipy.stats import gaussian_kde
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def construct_landscape(rc, bandwidth=0.3):
    """
    Construct energy landscape U(x) = -log P(x) in SPIB 3D space.
    
    Args:
        rc: (N, 3) SPIB reaction coordinates
        bandwidth: KDE bandwidth
    
    Returns:
        energy: (N,) energy values
    """
    kde = gaussian_kde(rc.T, bw_method=bandwidth)
    log_prob = kde.logpdf(rc.T)
    energy = -log_prob
    
    # Normalize to [0, 10]
    energy = (energy - energy.min()) / (energy.max() - energy.min()) * 10
    
    return energy

def main():
    print("[PHASE 3] Constructing energy landscape in SPIB space...")
    
    # Load SPIB reaction coordinates
    rc = np.load("output/spib_rc_3d.npy")  # (15000, 3)
    labels = np.load("output/nn_y.npy")  # (15000,)
    
    print(f"  RC shape: {rc.shape}")
    print(f"  RC ranges: RC1=[{rc[:,0].min():.3f}, {rc[:,0].max():.3f}], "
          f"RC2=[{rc[:,1].min():.3f}, {rc[:,1].max():.3f}], "
          f"RC3=[{rc[:,2].min():.3f}, {rc[:,2].max():.3f}]")
    
    # Construct energy landscape
    bandwidth = 0.3
    energy = construct_landscape(rc, bandwidth=bandwidth)
    
    print(f"  Energy range: [{energy.min():.3f}, {energy.max():.3f}]")
    print(f"  Energy mean: {energy.mean():.3f}, std: {energy.std():.3f}")
    
    # Print zone energy stats
    zone_names = ['Core', 'Periphery', 'Healthy']
    for z, name in enumerate(zone_names):
        mask = labels == z
        if mask.any():
            zone_energy = energy[mask]
            print(f"  {name}: mean={zone_energy.mean():.3f}, "
                  f"range=[{zone_energy.min():.3f}, {zone_energy.max():.3f}]")
    
    # Save
    np.save("output/spib_landscape.npy", energy.astype(np.float32))
    print(f"  Saved: output/spib_landscape.npy")
    
    # Visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    colors = ['#DC143C', '#FF8C00', '#2E8B57']
    
    # RC1 vs RC2 colored by energy
    ax = axes[0, 0]
    sc = ax.scatter(rc[:, 0], rc[:, 1], c=energy, cmap='RdYlBu_r', s=3, alpha=0.7)
    plt.colorbar(sc, ax=ax, label='Energy')
    ax.set_xlabel('RC1')
    ax.set_ylabel('RC2')
    ax.set_title(f'SPIB RC Space - Energy Landscape (bandwidth={bandwidth})')
    
    # RC1 vs RC3 colored by energy
    ax = axes[0, 1]
    sc = ax.scatter(rc[:, 0], rc[:, 2], c=energy, cmap='RdYlBu_r', s=3, alpha=0.7)
    plt.colorbar(sc, ax=ax, label='Energy')
    ax.set_xlabel('RC1')
    ax.set_ylabel('RC3')
    ax.set_title('SPIB RC Space - Energy Landscape')
    
    # Energy distribution by zone
    ax = axes[1, 0]
    for z, (name, color) in enumerate(zip(zone_names, colors)):
        mask = labels == z
        if mask.any():
            ax.hist(energy[mask], bins=50, alpha=0.6, label=name, color=color, density=True)
    ax.set_xlabel('Energy')
    ax.set_ylabel('Density')
    ax.set_title('Energy Distribution by Zone')
    ax.legend()
    
    # Zone centroids in RC space
    ax = axes[1, 1]
    for z, (name, color) in enumerate(zip(zone_names, colors)):
        mask = labels == z
        if mask.any():
            centroid = rc[mask].mean(axis=0)
            ax.scatter(centroid[0], centroid[1], c=color, s=200, label=f'{name} centroid', edgecolors='black')
            ax.scatter(rc[mask, 0], rc[mask, 1], c=color, s=1, alpha=0.3)
    ax.set_xlabel('RC1')
    ax.set_ylabel('RC2')
    ax.set_title('Zone Centroids in SPIB Space')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("output/spib_landscape_validation.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Validation plot: output/spib_landscape_validation.png")
    print("[PHASE 3] Complete.")

if __name__ == "__main__":
    main()