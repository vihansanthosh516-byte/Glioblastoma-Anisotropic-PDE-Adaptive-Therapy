#!/usr/bin/env python3
"""
22_spib_saddle_proof.py
========================
Saddle point proof using SPIB 3D reaction coordinates.

Replaces the 32D cVAE latent analysis with a kinetically-informed
3D SPIB space where Hessian estimation is numerically stable.
"""
from __future__ import annotations
import json, time, resource
from pathlib import Path

import numpy as np
import torch
from sklearn.cluster import KMeans


# ---- Config ----
K_NEIGHBORS = 200
RIDGE_ALPHA = 1e-2
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")


def mem_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}  RSS={mem_mb():.0f} MB", flush=True)


def load_data(device):
    log("Loading SPIB coordinates and landscape...")

    rc = np.load("output/spib_rc_3d.npy").astype(np.float32)
    energy = np.load("output/spib_landscape.npy").astype(np.float32)

    labels_1d = np.load("output/nn_y.npy")
    labels_int = labels_1d.astype(np.int64)

    log(f"  SPIB RC: {rc.shape}")
    log(f"  Energy: {energy.shape}")
    log(f"  Labels unique: {np.unique(labels_int)}")

    return (
        torch.from_numpy(rc).to(device),
        torch.from_numpy(energy).to(device),
        torch.from_numpy(labels_int).to(device),
    )


def compute_hessian(latent, energy, point, k=K_NEIGHBORS, alpha=RIDGE_ALPHA):
    """Ridge-regularized quadratic fit for 3x3 Hessian."""
    D = latent.shape[1]

    dists = torch.cdist(point.unsqueeze(0), latent).squeeze(0)
    k_eff = min(k, len(latent))
    _, knn_idx = torch.topk(dists, k_eff, largest=False)

    Z = latent[knn_idx] - point
    E = energy[knn_idx] - energy[knn_idx].mean()

    median_dist = dists[knn_idx].median().item()
    bandwidth = max(median_dist * 1.5, 1e-3)
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

    return 0.5 * (H + H.T)


def classify_hessian(H, noise_floor=0.035):
    """Classify using fixed absolute noise floor."""
    eigvals = torch.linalg.eigvalsh(H).cpu().numpy()

    pos = int((eigvals > noise_floor).sum())
    neg = int((eigvals < -noise_floor).sum())
    zero = int(((eigvals >= -noise_floor) & (eigvals <= noise_floor)).sum())

    # Classify: 
    # stable_minimum: no negative eigenvalues (neg == 0), at least one positive
    # unstable_maximum: no positive eigenvalues (pos == 0), at least one negative
    # saddle_point: both positive and negative eigenvalues exist
    # degenerate: all eigenvalues near zero
    if neg == 0 and pos > 0:
        ctype = "stable_minimum"
    elif pos == 0 and neg > 0:
        ctype = "unstable_maximum"
    elif pos > 0 and neg > 0:
        ctype = "saddle_point"
    else:
        ctype = "degenerate"

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
    k_eff = min(k, len(cells))
    if k_eff < 1:
        return torch.empty(0, latent.shape[1], device=latent.device)
    km = KMeans(n_clusters=k_eff, random_state=42, n_init=10)
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


def analyze_point(latent, energy, point, name):
    d = torch.cdist(point.unsqueeze(0), latent).squeeze(0)
    e_val = energy[d.argmin()].item()

    H = compute_hessian(latent, energy, point)
    eigvals, pos, neg, zero, ctype = classify_hessian(H)

    log(f"  {name}: E={e_val:.4f}, eig=[{eigvals.min():.4f}, {eigvals.max():.4f}], "
        f"+{pos}/-{neg}/{zero}z, type={ctype}")

    return {
        "name": name,
        "energy": e_val,
        "eigenvalue_min": float(eigvals.min()),
        "eigenvalue_max": float(eigvals.max()),
        "positive_eigenvalues": pos,
        "negative_eigenvalues": neg,
        "zero_eigenvalues": zero,
        "critical_point_type": ctype,
        "is_saddle": ctype == "saddle_point",
    }


def main():
    log(f"Device: {DEVICE}")
    print("=" * 60)
    print("SADDLE POINT ANALYSIS: SPIB 3D REACTION COORDINATES")
    print("=" * 60)

    latent, energy, labels = load_data(DEVICE)
    unique_zones = sorted(torch.unique(labels).cpu().numpy().tolist())
    zone_names = {z: f"Zone{z}" for z in unique_zones}
    if len(unique_zones) == 3:
        zone_names = {unique_zones[0]: "Healthy", unique_zones[1]: "Periphery", unique_zones[2]: "Core"}

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
        results.append(analyze_point(latent, energy, point, name))

    saddle_r = analyze_point(latent, energy, saddle_point, "Transition_Saddle")
    results.append(saddle_r)

    # Validation
    core_r = next(r for r in results if "Core_Attractor" in r["name"])
    healthy_r = next(r for r in results if "Healthy_Attractor" in r["name"])
    periphery_r = next((r for r in results if "Periphery_Attractor" in r["name"]), None)

    core_stable = core_r["critical_point_type"] == "stable_minimum"
    healthy_stable = healthy_r["critical_point_type"] == "stable_minimum"
    periphery_stable = periphery_r["critical_point_type"] == "stable_minimum" if periphery_r else True
    saddle_exists = saddle_r["is_saddle"]
    saddle_higher = saddle_r["energy"] > max(core_r["energy"], healthy_r["energy"])
    mixed_eig = saddle_r["positive_eigenvalues"] > 0 and saddle_r["negative_eigenvalues"] > 0

    print("\n" + "=" * 60)
    print("VALIDATION SUMMARY")
    print("=" * 60)
    print(f"  Core attractor stable:        {'PASS' if core_stable else 'FAIL'}  (type={core_r['critical_point_type']})")
    print(f"  Healthy attractor stable:     {'PASS' if healthy_stable else 'FAIL'}  (type={healthy_r['critical_point_type']})")
    if periphery_r:
        print(f"  Periphery attractor stable:   {'PASS' if periphery_stable else 'FAIL'}  (type={periphery_r['critical_point_type']})")
    print(f"  Saddle point exists:          {'PASS' if saddle_exists else 'FAIL'}  (type={saddle_r['critical_point_type']})")
    print(f"  Saddle energy > both:         {'PASS' if saddle_higher else 'FAIL'}")
    print(f"  Mixed Hessian eigenvalues:    {'PASS' if mixed_eig else 'FAIL'}")

    overall = core_stable and healthy_stable and periphery_stable and saddle_exists and mixed_eig
    print(f"\nOverall: {'*** SPIB SADDLE POINT CONFIRMED ***' if overall else 'CALIBRATION REQUIRED'}")

    proof = {
        "method": "SPIB 3D reaction coordinate + ridge Hessian",
        "spib_rc_dim": 3,
        "k_neighbors": K_NEIGHBORS,
        "attractors": {r["name"]: r for r in results if "Attractor" in r["name"]},
        "saddle": saddle_r,
        "validation": {
            "core_stable": bool(core_stable),
            "healthy_stable": bool(healthy_stable),
            "periphery_stable": bool(periphery_stable),
            "saddle_exists": bool(saddle_exists),
            "saddle_higher_than_both": bool(saddle_higher),
            "mixed_eigenvalues": bool(mixed_eig),
            "overall_pass": bool(overall),
        },
    }
    out = Path("output/spib_saddle_point_metrics.json")
    with open(out, "w") as f:
        json.dump(proof, f, indent=2)
    log(f"Saved: {out}")
    log("DONE.")


if __name__ == "__main__":
    main()