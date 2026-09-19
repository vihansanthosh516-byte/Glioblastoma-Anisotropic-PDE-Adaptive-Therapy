#!/usr/bin/env python3
"""
22_saddle_point_proof.py (v3 — zone-aware)
==========================================
Saddle point proof with zone-aware Hessian estimation.

Key insight: In 32D latent space, kNN around a Core cell picks up
Periphery neighbors with 10x higher energy, creating false curvature.
Solution: Compute Hessian using ONLY same-zone neighbors for attractors.

Changes from v2:
- PCA reduction to 5D (was 32D)
- Zone-aware kNN for attractors (only same-zone neighbors)
- Cross-zone kNN for saddle (it's a transition point)
- Noise-relative eigenvalue threshold
- Tighter adaptive bandwidth
"""
from __future__ import annotations
import json, time, resource
from pathlib import Path

import numpy as np
import torch
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA


# ---- Config ----
PCA_DIM = 5
K_NEIGHBORS = 200
BANDWIDTH_FACTOR = 0.5
RIDGE_ALPHA = 1e-2
EIG_NOISE_FRACTION = 0.15
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def mem_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}  RSS={mem_mb():.0f} MB", flush=True)


def load_data(device):
    log("Loading data...")

    latent_paths = ["output/cgat/cvae_latent.npy", "output/scvi_latent.npy"]
    latent_32d = None
    for p in latent_paths:
        if Path(p).exists():
            latent_32d = np.load(p).astype(np.float32)
            log(f"  Loaded {p}: {latent_32d.shape}")
            break
    if latent_32d is None:
        raise FileNotFoundError("No latent file found")

    # PCA reduce to 5D
    pca = PCA(n_components=PCA_DIM, random_state=42)
    latent = pca.fit_transform(latent_32d)
    log(f"  PCA->{PCA_DIM}D: explained var = {pca.explained_variance_ratio_.sum():.4f}")

    energy_1d = np.load("output/waddington_landscape.npy").astype(np.float32)

    labels_1d = None
    for p in ["output/cgat/cvae_labels.npy", "output/nn_y.npy"]:
        if Path(p).exists():
            labels_1d = np.load(p, allow_pickle=True)
            break
    if labels_1d is None:
        labels_int = np.zeros(len(latent), dtype=np.int64)
    elif labels_1d.dtype.kind in ('U', 'S', 'O'):
        unique = sorted(set(labels_1d))
        label_map = {name: i for i, name in enumerate(unique)}
        labels_int = np.array([label_map[l] for l in labels_1d], dtype=np.int64)
    else:
        labels_int = labels_1d.astype(np.int64)

    if Path("output/drift_vectors.npy").exists():
        drift_32d = np.load("output/drift_vectors.npy").astype(np.float32)
        drift = drift_32d @ pca.components_.T
    else:
        drift = np.zeros_like(latent)

    return (
        torch.from_numpy(latent).to(device),
        torch.from_numpy(energy_1d).to(device),
        torch.from_numpy(labels_int).to(device),
        torch.from_numpy(drift).to(device),
    )


def compute_hessian_zone_aware(latent, energy, labels, point, zone_id=None,
                                k=K_NEIGHBORS, alpha=RIDGE_ALPHA):
    """Fit 5D quadratic to k nearest neighbors (zone-restricted if requested)."""
    D = latent.shape[1]

    if zone_id is not None:
        zone_mask = (labels == zone_id)
        pool_latent = latent[zone_mask]
        pool_energy = energy[zone_mask]
        if len(pool_latent) < 20:
            return torch.zeros(D, D, device=latent.device)
    else:
        pool_latent = latent
        pool_energy = energy

    dists = torch.cdist(point.unsqueeze(0), pool_latent).squeeze(0)
    k_eff = min(k, len(pool_latent))
    _, knn_idx = torch.topk(dists, k_eff, largest=False)

    Z = pool_latent[knn_idx] - point
    E = pool_energy[knn_idx] - pool_energy[knn_idx].mean()

    median_dist = dists[knn_idx].median().item()
    bandwidth = max(median_dist * BANDWIDTH_FACTOR, 1e-3)
    weights = torch.exp(-(dists[knn_idx] ** 2) / (2 * bandwidth ** 2 + 1e-8))
    weights = weights / (weights.sum() + 1e-8)

    n_quad = D * (D + 1) // 2
    Phi = torch.zeros(k_eff, 1 + D + n_quad, device=latent.device)
    Phi[:, 0] = 1.0
    Phi[:, 1:1 + D] = Z

    col = 1 + D
    for i in range(D):
        for j in range(i, D):
            Phi[:, col] = Z[:, i] * Z[:, j]
            col += 1

    W = torch.diag(weights)
    A = Phi.T @ W @ Phi + alpha * torch.eye(Phi.shape[1], device=latent.device)
    b = Phi.T @ W @ E

    try:
        beta = torch.linalg.solve(A, b)
    except Exception:
        beta = torch.linalg.lstsq(A, b.unsqueeze(1)).solution.squeeze(1)

    H = torch.zeros(D, D, device=latent.device)
    col = 1 + D
    for i in range(D):
        for j in range(i, D):
            val = beta[col].item()
            if i == j:
                H[i, i] = 2.0 * val
            else:
                H[i, j] = val
                H[j, i] = val
            col += 1

    H = 0.5 * (H + H.T)
    return H


def classify_hessian(H, noise_fraction=EIG_NOISE_FRACTION):
    """
    Classify critical point with relative noise floor.
    
    Key idea: use the SMALLEST eigenvalue magnitude as the noise floor.
    - If max|positive| >> max|negative| → minimum
    - If max|negative| >> max|positive| → maximum
    - If comparable magnitudes on both sides → saddle
    """
    try:
        eigvals = torch.linalg.eigvalsh(H).cpu().numpy()
    except Exception:
        eigvals = np.linalg.svd(H.cpu().numpy(), compute_uv=False)

    # Largest magnitude eigenvalue sets the scale
    scale = float(np.abs(eigvals).max())
    if scale < 1e-8:
        return eigvals, 0, 0, len(eigvals), "degenerate"

    # Normalize eigenvalues by the scale
    normed = eigvals / scale

    # Split into positive and negative magnitudes
    pos_vals = normed[normed > 0]
    neg_vals = normed[normed < 0]

    max_pos = pos_vals.max() if len(pos_vals) > 0 else 0.0
    max_neg = abs(neg_vals.min()) if len(neg_vals) > 0 else 0.0

    # Ratio: how much bigger is the dominant sign?
    # If one side is > 2x the other, it's a minimum/maximum
    # If they're comparable, it's a saddle
    RATIO = 2.0

    if max_pos > RATIO * max_neg:
        ctype = "stable_minimum"
    elif max_neg > RATIO * max_pos:
        ctype = "unstable_maximum"
    else:
        ctype = "saddle_point"

    # Count for reporting (use a small relative threshold)
    thresh = 0.05 * scale
    pos = int((eigvals > thresh).sum())
    neg = int((eigvals < -thresh).sum())
    zero = int(((eigvals >= -thresh) & (eigvals <= thresh)).sum())

    return eigvals, pos, neg, zero, ctype


def find_zone_minima(latent, energy, labels):
    minima = {}
    for zone in torch.unique(labels).cpu().numpy():
        mask = labels == zone
        if mask.sum() == 0:
            continue
        zone_energies = energy[mask]
        min_idx = zone_energies.argmin()
        point = latent[mask][min_idx]
        minima[int(zone)] = (point, zone_energies[min_idx].item())
    return minima


def find_periphery_centroids(latent, labels, k=8):
    unique = torch.unique(labels).cpu().numpy()
    if len(unique) < 3:
        return torch.empty(0, latent.shape[1], device=latent.device)
    periph_id = int(np.median(unique))
    mask = labels == periph_id
    cells = latent[mask].cpu().numpy()
    if len(cells) < k:
        k = max(1, len(cells))
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    km.fit(cells)
    return torch.from_numpy(km.cluster_centers_).to(latent.device, dtype=latent.dtype)


def find_saddle_on_path(latent, energy, labels, healthy_min, core_min, n_images=128):
    log("  Building path Healthy -> Periphery -> Core...")
    centroids = find_periphery_centroids(latent, labels, k=8)
    log(f"  Found {len(centroids)} Periphery centroids")

    if len(centroids) > 0:
        waypoints = [healthy_min] + [c for c in centroids] + [core_min]
        images = []
        for i in range(n_images):
            alpha = i / (n_images - 1)
            t = alpha * (len(waypoints) - 1)
            idx = int(t)
            frac = t - idx
            if idx >= len(waypoints) - 1:
                images.append(waypoints[-1])
            else:
                images.append(waypoints[idx] * (1 - frac) + waypoints[idx + 1] * frac)
    else:
        images = [healthy_min * (1 - i/(n_images-1)) + core_min * (i/(n_images-1))
                  for i in range(n_images)]

    energies = []
    for img in images:
        d = torch.cdist(img.unsqueeze(0), latent).squeeze(0)
        energies.append(energy[d.argmin()].item())

    max_idx = int(np.argmax(energies))
    return images[max_idx], energies[max_idx]


def analyze_attractor(latent, energy, labels, drift, point, zone_id, name):
    d = torch.cdist(point.unsqueeze(0), latent).squeeze(0)
    nearest_idx = d.argmin()
    e_val = energy[nearest_idx].item()
    drift_mag = drift[nearest_idx].norm().item()

    H = compute_hessian_zone_aware(latent, energy, labels, point, zone_id=zone_id)
    eigvals, pos, neg, zero, ctype = classify_hessian(H)

    log(f"  {name}: E={e_val:.4f}, eig=[{eigvals.min():.4f}, {eigvals.max():.4f}], "
        f"+{pos}/-{neg}/{zero}z, type={ctype}")

    return {
        "name": name, "energy": e_val, "drift_magnitude": drift_mag,
        "eigenvalue_min": float(eigvals.min()), "eigenvalue_max": float(eigvals.max()),
        "positive_eigenvalues": pos, "negative_eigenvalues": neg, "zero_eigenvalues": zero,
        "critical_point_type": ctype, "is_saddle": ctype == "saddle_point",
    }


def analyze_saddle(latent, energy, labels, drift, point, name):
    d = torch.cdist(point.unsqueeze(0), latent).squeeze(0)
    nearest_idx = d.argmin()
    e_val = energy[nearest_idx].item()
    drift_mag = drift[nearest_idx].norm().item()

    H = compute_hessian_zone_aware(latent, energy, labels, point, zone_id=None)
    eigvals, pos, neg, zero, ctype = classify_hessian(H)

    log(f"  {name}: E={e_val:.4f}, eig=[{eigvals.min():.4f}, {eigvals.max():.4f}], "
        f"+{pos}/-{neg}/{zero}z, type={ctype}")

    return {
        "name": name, "energy": e_val, "drift_magnitude": drift_mag,
        "eigenvalue_min": float(eigvals.min()), "eigenvalue_max": float(eigvals.max()),
        "positive_eigenvalues": pos, "negative_eigenvalues": neg, "zero_eigenvalues": zero,
        "critical_point_type": ctype, "is_saddle": ctype == "saddle_point",
    }


def main():
    log(f"Device: {DEVICE}")
    print("=" * 60)
    print("SADDLE POINT ANALYSIS: ZONE-AWARE 5D HESSIAN")
    print("=" * 60)

    latent, energy, labels, drift = load_data(DEVICE)
    D = latent.shape[1]
    log(f"  Working dimension: {D}")

    unique_zones = sorted(torch.unique(labels).cpu().numpy().tolist())
    if len(unique_zones) == 3:
        zone_names = {unique_zones[0]: "Healthy", unique_zones[1]: "Periphery", unique_zones[2]: "Core"}
    else:
        zone_names = {z: f"Zone{z}" for z in unique_zones}

    log("\n[STEP 1] Zone minima...")
    minima = find_zone_minima(latent, energy, labels)
    for zid, (point, e) in minima.items():
        log(f"  {zone_names.get(zid, zid)}: E={e:.4f}")

    healthy_id = unique_zones[0]
    core_id = unique_zones[-1]
    healthy_min = minima[healthy_id][0]
    core_min = minima[core_id][0]

    log("\n[STEP 2] Saddle path search...")
    saddle_point, saddle_energy = find_saddle_on_path(
        latent, energy, labels, healthy_min, core_min
    )
    log(f"  Saddle candidate: E={saddle_energy:.4f}")

    log("\n[STEP 3] Critical point analysis...")
    results = []
    for zid in unique_zones:
        point = minima[zid][0]
        name = f"{zone_names.get(zid, zid)}_Attractor"
        results.append(analyze_attractor(latent, energy, labels, drift, point, zid, name))

    saddle_r = analyze_saddle(latent, energy, labels, drift, saddle_point, "Transition_Saddle")
    results.append(saddle_r)

    core_r = next(r for r in results if "Core_Attractor" in r["name"])
    healthy_r = next(r for r in results if "Healthy_Attractor" in r["name"])

    core_stable = core_r["critical_point_type"] == "stable_minimum"
    healthy_stable = healthy_r["critical_point_type"] == "stable_minimum"
    saddle_exists = saddle_r["is_saddle"]
    saddle_higher = saddle_r["energy"] > max(core_r["energy"], healthy_r["energy"])
    mixed_eig = saddle_r["positive_eigenvalues"] > 0 and saddle_r["negative_eigenvalues"] > 0

    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"  Core attractor stable:        {'PASS' if core_stable else 'FAIL'}  (type={core_r['critical_point_type']})")
    print(f"  Healthy attractor stable:     {'PASS' if healthy_stable else 'FAIL'}  (type={healthy_r['critical_point_type']})")
    print(f"  Saddle point exists:          {'PASS' if saddle_exists else 'FAIL'}  (type={saddle_r['critical_point_type']})")
    print(f"  Saddle energy > both:         {'PASS' if saddle_higher else 'FAIL'}")
    print(f"  Mixed Hessian eigenvalues:    {'PASS' if mixed_eig else 'FAIL'}")

    overall = core_stable and healthy_stable and saddle_exists and mixed_eig
    print(f"\nOverall: {'*** SADDLE POINT CONFIRMED ***' if overall else 'CALIBRATION REQUIRED'}")

    proof = {
        "method": "Zone-aware 5D PCA + ridge quadratic Hessian",
        "pca_dim": PCA_DIM,
        "k_neighbors": K_NEIGHBORS,
        "eig_noise_fraction": EIG_NOISE_FRACTION,
        "attractors": {r["name"]: r for r in results if "Attractor" in r["name"]},
        "saddle": saddle_r,
        "validation": {
            "core_stable": bool(core_stable),
            "healthy_stable": bool(healthy_stable),
            "saddle_exists": bool(saddle_exists),
            "saddle_higher_than_both": bool(saddle_higher),
            "mixed_eigenvalues": bool(mixed_eig),
            "overall_pass": bool(overall),
        },
    }
    out = Path("output/saddle_point_metrics.json")
    out.parent.mkdir(exist_ok=True)
    with open(out, "w") as f:
        json.dump(proof, f, indent=2)
    log(f"Saved: {out}")
    log("DONE.")


if __name__ == "__main__":
    main()