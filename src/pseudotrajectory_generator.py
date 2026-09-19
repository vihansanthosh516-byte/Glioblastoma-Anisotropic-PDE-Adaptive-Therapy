#!/usr/bin/env python3
"""
Phase 1: Generate pseudo-trajectories from phenotypic velocity field.
SPIB requires time-resolved trajectory data (X_t, X_{t+dt}) pairs.
"""

import numpy as np
import os

def generate_pseudotrajectories(latent, velocity, labels, dt=0.05, n_steps=3):
    """
    Generate (past, present, future) triples from velocity field.
    
    Args:
        latent: (N, D) latent coordinates
        velocity: (N, D) phenotypic velocity vectors
        labels: (N,) zone labels
        dt: time step for pseudo-trajectory
        n_steps: number of forward steps to generate
    
    Returns:
        traj_data: (N * n_steps, D) concatenated trajectory frames
        traj_labels: (N * n_steps,) zone labels
    """
    N, D = latent.shape
    traj_frames = []
    traj_labels = []
    
    # Ensure labels are integer
    labels_int = labels.astype(np.int64)
    
    for step in range(n_steps):
        # Forward step
        latent_fwd = latent + velocity * dt * (step + 1)
        traj_frames.append(latent_fwd)
        traj_labels.append(labels_int)
    
    # Concatenate into trajectory format
    traj_data = np.vstack(traj_frames)
    traj_labels = np.concatenate(traj_labels)
    
    return traj_data, traj_labels


def main():
    print("[PHASE 1] Generating pseudo-trajectories for SPIB...")
    
    # Load data
    latent = np.load("output/scvi_latent.npy")
    velocity = np.load("output/phenotypic_velocity.npy")
    labels = np.load("output/nn_y.npy")
    
    print(f"  Latent shape: {latent.shape}")
    print(f"  Velocity shape: {velocity.shape}")
    print(f"  Labels shape: {labels.shape}")
    print(f"  Unique labels: {np.unique(labels)}")
    
    # Generate pseudo-trajectories with longer lag time
    dt = 0.2
    n_steps = 10
    traj_data, traj_labels = generate_pseudotrajectories(latent, velocity, labels, dt=dt, n_steps=n_steps)
    
    print(f"  Trajectory data shape: {traj_data.shape}")  # (45000, 32)
    print(f"  Trajectory labels shape: {traj_labels.shape}")
    
    # Save
    os.makedirs("output", exist_ok=True)
    np.save("output/spib_traj_data.npy", traj_data)
    np.save("output/spib_traj_labels.npy", traj_labels)
    
    print(f"  Saved: output/spib_traj_data.npy")
    print(f"  Saved: output/spib_traj_labels.npy")
    
    # Validation: check zone structure preservation
    from sklearn.decomposition import PCA
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    pca = PCA(n_components=2, random_state=42)
    traj_2d = pca.fit_transform(traj_data)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    zone_names = ['Core', 'Periphery', 'Healthy']
    colors = ['#DC143C', '#FF8C00', '#2E8B57']
    
    for z, (name, color) in enumerate(zip(zone_names, colors)):
        mask = traj_labels == z
        if mask.any():
            ax.scatter(traj_2d[mask, 0], traj_2d[mask, 1], 
                       c=color, s=1, alpha=0.5, label=name, rasterized=True)
    
    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title(f'Pseudo-trajectories in PCA space (dt={dt}, steps={n_steps})')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig("output/spib_traj_validation.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Validation plot: output/spib_traj_validation.png")
    print("[PHASE 1] Complete.")

if __name__ == "__main__":
    main()