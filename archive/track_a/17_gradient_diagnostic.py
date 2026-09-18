import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import umap

def calculate_shannon_entropy(fractions):
    eps = 1e-10
    probs = np.clip(fractions, eps, 1.0)
    probs = probs / probs.sum(axis=1, keepdims=True)
    return -np.sum(probs * np.log2(probs), axis=1)

def main():
    print("[LOAD] Loading data artifacts...")
    y = np.load("output/nn_y.npy")
    scvi_latent = np.load("output/scvi_latent.npy")
    nmf_fractions = np.load("output/nmf_fractions.npy")
    
    zone_names = {0: "Core", 1: "Periphery", 2: "Healthy"}
    meta_module_names = ["AC", "MES", "NPC", "OPC"]
    
    print(f"   Labels: {y.shape}")
    print(f"   scVI Latent: {scvi_latent.shape}")
    print(f"   NMF Fractions: {nmf_fractions.shape}")
    
    print("\n[ENTROPY] Calculating Shannon Entropy per cell...")
    entropy = calculate_shannon_entropy(nmf_fractions)
    
    print("\n=== ENTROPY BY SPATIAL ZONE ===")
    for zone_id in range(3):
        mask = y == zone_id
        zone_entropy = entropy[mask]
        print(f"   {zone_names[zone_id]}: Mean={zone_entropy.mean():.4f}, Std={zone_entropy.std():.4f}, N={mask.sum()}")
    
    print("\n[VISUALIZATION] Computing UMAP on scVI latent space...")
    umap_model = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
    latent_2d = umap_model.fit_transform(scvi_latent)
    
    print("[PLOT] Generating diagnostic figure...")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    colors = ['#e74c3c', '#f39c12', '#2ecc71']
    for zone_id in range(3):
        mask = y == zone_id
        axes[0].scatter(
            latent_2d[mask, 0], latent_2d[mask, 1],
            c=colors[zone_id], label=zone_names[zone_id],
            alpha=0.6, s=8, edgecolors='none'
        )
    axes[0].set_title("scVI Latent Space (UMAP) — Spatial Zones", fontsize=14, fontweight='bold')
    axes[0].set_xlabel("UMAP 1")
    axes[0].set_ylabel("UMAP 2")
    axes[0].legend(title="Zone", frameon=True, fontsize=11)
    axes[0].grid(True, alpha=0.3)
    
    avg_proportions = np.zeros((3, 4))
    for zone_id in range(3):
        mask = y == zone_id
        avg_proportions[zone_id] = nmf_fractions[mask].mean(axis=0)
    
    bottom = np.zeros(3)
    zone_labels = [zone_names[i] for i in range(3)]
    for k in range(4):
        axes[1].bar(
            zone_labels, avg_proportions[:, k],
            bottom=bottom, label=meta_module_names[k],
            color=plt.cm.Set2(k), edgecolor='white', linewidth=0.5
        )
        bottom += avg_proportions[:, k]
    
    axes[1].set_title("Average NMF Meta-Module Proportions by Zone", fontsize=14, fontweight='bold')
    axes[1].set_ylabel("Fractional Abundance")
    axes[1].set_ylim(0, 1.05)
    axes[1].legend(title="Meta-Module", frameon=True, fontsize=11, loc='upper right')
    axes[1].grid(True, alpha=0.3, axis='y')
    
    for i, zone_id in enumerate(range(3)):
        axes[1].text(i, 1.02, f"H={entropy[y==zone_id].mean():.3f}", 
                     ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    os.makedirs("output", exist_ok=True)
    fig.savefig("output/gradient_failure_analysis.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("[OK] Diagnostic plot saved to output/gradient_failure_analysis.png")
    
    print("\n=== AVERAGE NMF PROPORTIONS BY ZONE ===")
    for zone_id in range(3):
        props = avg_proportions[zone_id]
        print(f"   {zone_names[zone_id]}: AC={props[0]:.3f}, MES={props[1]:.3f}, NPC={props[2]:.3f}, OPC={props[3]:.3f}")

if __name__ == "__main__":
    main()