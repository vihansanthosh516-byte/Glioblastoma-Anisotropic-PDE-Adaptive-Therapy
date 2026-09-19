#!/usr/bin/env python3
"""
Phase 2: Train SPIB to learn 3D reaction coordinate.
"""

import os
import numpy as np
import torch
import spib
from spib.spib import SPIB
from spib.utils import TimeLaggedDataset, DataNormalize

def main():
    print("[PHASE 2] Training SPIB...")
    
    # Load pseudo-trajectories
    traj_data = np.load("output/spib_traj_data.npy")  # (45000, 32)
    traj_labels = np.load("output/spib_traj_labels.npy")  # (45000,)
    
    print(f"  Trajectory data shape: {traj_data.shape}")
    print(f"  Trajectory labels shape: {traj_labels.shape}")
    
    # Convert to list of trajectories (each cell has n_steps frames)
    n_cells = 15000
    n_steps = 10
    
    data_list = []
    label_list = []
    
    for i in range(n_cells):
        data_list.append(traj_data[i * n_steps:(i + 1) * n_steps])
        label_list.append(traj_labels[i * n_steps:(i + 1) * n_steps])
    
    print(f"  Number of trajectories: {len(data_list)}")
    print(f"  Trajectory lengths: {[len(d) for d in data_list[:5]]}")
    
    # Create datasets
    # With dt=0.2, n_steps=10, we have 9 steps per trajectory.
    # Use lagtime=5 to look further into the future for better kinetic separation.
    lagtime = 5
    output_dim = 3  # 3 metastable states: Core, Periphery, Healthy
    
    # Normalize data
    all_data = np.vstack(data_list)
    data_mean = all_data.mean(axis=0)
    data_std = all_data.std(axis=0) + 1e-8
    
    normalizer = DataNormalize(mean=data_mean, std=data_std)
    
    # Apply normalization to each trajectory
    data_list_norm = []
    for traj in data_list:
        data_list_norm.append((traj - data_mean) / data_std)
    
    # Create train/test split (80/20)
    n_train = int(0.8 * n_cells)
    train_data = data_list_norm[:n_train]
    train_labels = label_list[:n_train]
    test_data = data_list_norm[n_train:]
    test_labels = label_list[n_train:]
    
    train_dataset = TimeLaggedDataset(
        train_data, train_labels, 
        lagtime=lagtime, 
        output_dim=output_dim,
        device=torch.device('cpu')
    )
    
    test_dataset = TimeLaggedDataset(
        test_data, test_labels,
        lagtime=lagtime,
        output_dim=output_dim,
        device=torch.device('cpu')
    )
    
    print(f"  Train dataset: {train_dataset.past_data.shape[0]} pairs")
    print(f"  Test dataset: {test_dataset.past_data.shape[0]} pairs")
    
    # Initialize SPIB
    spib_model = SPIB(
        output_dim=output_dim,
        data_shape=data_list[0][0].shape,  # (32,)
        encoder_type='Nonlinear',
        z_dim=3,  # 3D reaction coordinate
        lagtime=lagtime,
        beta=0.01,  # IB trade-off parameter
        learning_rate=0.001,
        neuron_num1=64,
        neuron_num2=64,
        UpdateLabel=True,
        device=torch.device('cpu'),
        path='output/spib_model'
    )
    
    print("  SPIB model initialized.")
    
    # Train
    print("  Starting training...")
    spib_model.fit(
        train_dataset=train_dataset,
        test_dataset=test_dataset,
        batch_size=2048,
        tolerance=0.01,
        patience=2,
        refinements=8,
        mask_threshold=0
    )
    
    print("[PHASE 2] Training complete!")
    
    # Save the normalizer parameters
    np.save("output/spib_normalizer_mean.npy", data_mean)
    np.save("output/spib_normalizer_std.npy", data_std)
    
    # Extract reaction coordinates for all cells
    spib_model.eval()
    with torch.no_grad():
        # Encode all normalized data
        all_data_norm = np.vstack(data_list_norm)
        all_data_tensor = torch.from_numpy(all_data_norm).float().to(torch.device('cpu'))
        
        z, _ = spib_model.encode(all_data_tensor)
        rc_all = z.cpu().numpy()  # (45000, 3)
        
        # First 15000 are original cells (step 0)
        rc_original = rc_all[:n_cells]
        
    np.save("output/spib_rc_3d.npy", rc_original)
    print(f"  Saved SPIB RC: output/spib_rc_3d.npy (shape: {rc_original.shape})")
    
    # Also save the full RC for analysis
    np.save("output/spib_rc_all.npy", rc_all)
    print(f"  Saved all SPIB RC: output/spib_rc_all.npy (shape: {rc_all.shape})")
    
    # Validation plot
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    
    labels = np.load("output/nn_y.npy")
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    zone_names = ['Core', 'Periphery', 'Healthy']
    colors = ['#DC143C', '#FF8C00', '#2E8B57']
    
    # Plot RC1 vs RC2
    ax = axes[0]
    for z, (name, color) in enumerate(zip(zone_names, colors)):
        mask = labels == z
        if mask.any():
            ax.scatter(rc_original[mask, 0], rc_original[mask, 1], 
                       c=color, s=3, alpha=0.6, label=name, rasterized=True)
    ax.set_xlabel('RC1')
    ax.set_ylabel('RC2')
    ax.set_title('SPIB Reaction Coordinates (RC1 vs RC2)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Plot RC1 vs RC3
    ax = axes[1]
    for z, (name, color) in enumerate(zip(zone_names, colors)):
        mask = labels == z
        if mask.any():
            ax.scatter(rc_original[mask, 0], rc_original[mask, 2], 
                       c=color, s=3, alpha=0.6, label=name, rasterized=True)
    ax.set_xlabel('RC1')
    ax.set_ylabel('RC3')
    ax.set_title('SPIB Reaction Coordinates (RC1 vs RC3)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("output/spib_rc_validation.png", dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"  Validation plot: output/spib_rc_validation.png")
    
    # Also compute state assignments
    with torch.no_grad():
        # Get decoder output (state probabilities)
        log_probs = spib_model.decode(z)
        probs = torch.exp(log_probs).cpu().numpy()
        state_assignments = probs.argmax(axis=1)  # (45000,)
        
    np.save("output/spib_state_assignments.npy", state_assignments)
    print(f"  State assignments: {np.bincount(state_assignments[:n_cells])}")

if __name__ == "__main__":
    main()