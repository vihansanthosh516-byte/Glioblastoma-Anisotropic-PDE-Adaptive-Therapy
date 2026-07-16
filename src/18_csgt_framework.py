import os
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import umap
from scipy.stats import kruskal
from scipy.spatial.distance import cdist

def calculate_shannon_entropy(fractions):
    eps = 1e-10
    probs = np.clip(fractions, eps, 1.0)
    probs = probs / probs.sum(axis=1, keepdims=True)
    return -np.sum(probs * np.log2(probs), axis=1)

def main():
    print("[CSGT] Loading precomputed latent spaces...")
    scvi_latent = np.load("output/scvi_latent.npy")
    nmf_fractions = np.load("output/nmf_fractions.npy")
    y = np.load("output/nn_y.npy")
    
    zone_names = {0: "Core", 1: "Periphery", 2: "Healthy"}
    zone_order = [2, 1, 0]
    
    print(f"   scVI Latent: {scvi_latent.shape}")
    print(f"   NMF Fractions: {nmf_fractions.shape}")
    print(f"   Labels: {y.shape}")
    
    print("\n[CSGT] Computing Healthy Centroid in VAE space...")
    healthy_mask = y == 2
    z_healthy_centroid = scvi_latent[healthy_mask].mean(axis=0)
    print(f"   Centroid computed from {healthy_mask.sum()} healthy cells")
    
    print("\n[CSGT] Executing Microenvironmental Transition Score equation...")
    alpha = 0.5
    beta = 0.5
    sigma = 1.0
    
    entropy = calculate_shannon_entropy(nmf_fractions)
    distances = cdist(scvi_latent, z_healthy_centroid.reshape(1, -1), metric='euclidean').flatten()
    spatial_term = np.exp(-distances**2 / (2 * sigma**2))
    
    T_scores = alpha * entropy + beta * spatial_term
    
    print(f"   T-score range: [{T_scores.min():.4f}, {T_scores.max():.4f}]")
    print(f"   T-score mean: {T_scores.mean():.4f} ± {T_scores.std():.4f}")
    
    print("\n[CSGT] Statistical Gradient Analysis (Kruskal-Wallis)...")
    zone_data = [T_scores[y == z] for z in zone_order]
    stat, p_value = kruskal(*zone_data)
    is_significant = bool(p_value < 0.001)  # Force raw Python boolean
    
    print(f"   H-statistic: {stat:.4f}")
    print(f"   p-value: {p_value:.2e}")
    print(f"   Significant gradient: {'YES' if is_significant else 'NO'}")
    
    for z in zone_order:
        zone_T = T_scores[y == z]
        print(f"   {zone_names[z]}: Mean={zone_T.mean():.4f}, Median={np.median(zone_T):.4f}, N={len(zone_T)}")
    
    print("\n[CSGT] Computing UMAP for visualization...")
    umap_model = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.1)
    latent_2d = umap_model.fit_transform(scvi_latent)
    
    print("[CSGT] Generating publication-grade visual proof...")
    fig = plt.figure(figsize=(18, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1])
    
    ax1 = fig.add_subplot(gs[0])
    scatter = ax1.scatter(
        latent_2d[:, 0], latent_2d[:, 1],
        c=T_scores, cmap='plasma', s=6, alpha=0.7, edgecolors='none'
    )
    cbar = plt.colorbar(scatter, ax=ax1, shrink=0.8, pad=0.01)
    cbar.set_label(r'Microenvironmental Transition Score $\mathcal{T}_i$', fontsize=12, fontweight='bold')
    
    # Raw string prefix r'' prevents LaTeX escaping syntax bugs
    ax1.set_title(r"Panel A: Continuous CSGT Gradient Wave" + "\n" + r"(UMAP colored by $\mathcal{T}_i$)", fontsize=14, fontweight='bold')
    ax1.set_xlabel("UMAP 1")
    ax1.set_ylabel("UMAP 2")
    ax1.grid(True, alpha=0.2)
    
    ax2 = fig.add_subplot(gs[1])
    violin_labels = [zone_names[z] for z in zone_order]
    colors = ['#2ecc71', '#f39c12', '#e74c3c']
    
    parts = ax2.violinplot(zone_data, positions=[0, 1, 2], widths=0.7, showmeans=True, showmedians=True, showextrema=True)
    
    for i, pc in enumerate(parts['bodies']):
        pc.set_facecolor(colors[i])
        pc.set_edgecolor('black')
        pc.set_alpha(0.7)
        pc.set_linewidth(1)
        
    for part_name in ('cmeans', 'cmedians', 'cbars', 'cmins', 'cmaxes'):
        parts[part_name].set_edgecolor('black')
        parts[part_name].set_linewidth(1.5)
        
    ax2.set_xticks([0, 1, 2])
    ax2.set_xticklabels(violin_labels, fontsize=12)
    ax2.set_ylabel(r'$\mathcal{T}_i$ Score', fontsize=12, fontweight='bold')
    ax2.set_title("Panel B: Violin Distribution by Zone\n(Periphery as Intermediate State)", fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    for i, (data, color) in enumerate(zip(zone_data, colors)):
        ax2.text(i, data.max() + 0.02, f"n={len(data)}\nμ={data.mean():.3f}\nH={stat:.1f}",
                 ha='center', va='bottom', fontsize=9, fontweight='bold', color=color,
                 bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=color, alpha=0.9))
        
    fig.suptitle("Continuous State-Gradient Trajectory (CSGT) Model — Mathematical Proof",
                 fontsize=16, fontweight='bold', y=0.98)
    
    plt.tight_layout()
    os.makedirs("output", exist_ok=True)
    fig.savefig("output/csgt_mathematical_proof.png", dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("[OK] Mathematical proof saved to output/csgt_mathematical_proof.png")
    
    np.save("output/csgt_transition_scores.npy", T_scores)
    
    metrics = {
        "alpha": alpha,
        "beta": beta,
        "sigma": sigma,
        "healthy_centroid": z_healthy_centroid.tolist(),
        "t_score_stats": {
            "min": float(T_scores.min()),
            "max": float(T_scores.max()),
            "mean": float(T_scores.mean()),
            "std": float(T_scores.std())
        },
        "kruskal_wallis": {
            "h_statistic": float(stat),
            "p_value": float(p_value),
            "significant": is_significant  # Now clean JSON standard type
        },
        "zone_stats": {
            zone_names[z]: {
                "mean": float(T_scores[y == z].mean()),
                "median": float(np.median(T_scores[y == z])),
                "std": float(T_scores[y == z].std()),
                "n": int((y == z).sum())
            }
            for z in zone_order
        }
    }
    
    with open("output/csgt_metrics.json", "w") as f:
        json.dump(metrics, f, indent=4)
    print("[OK] CSGT metrics saved to output/csgt_metrics.json")
    print("[OK] Transition scores saved to output/csgt_transition_scores.npy")

if __name__ == "__main__":
    main()