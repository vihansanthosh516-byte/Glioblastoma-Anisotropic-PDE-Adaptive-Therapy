#!/usr/bin/env python3
"""
src/rl/train_cohort.py
======================
Train behavioral cloning policy across the real MU-Glioma cohort's
rho distribution (Phase 3).

Literature basis:
- arXiv:2609.12264: Reinforcement learning and adaptive control policies developed
  for a single nominal patient model fail under parametric heterogeneity.
  Training across cohort-level heterogeneity is required for robust generalization.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
from typing import Dict, List, Tuple

import numpy as np
import torch

PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "src" / "rl"))

from mu_glioma_loader import load_mu_glioma_params, real_cohort_stats
from behavioral_cloning import PolicyNetwork, generate_env_demonstrations


def get_cohort_rho_distribution() -> List[float]:
    """Load valid positive rho values from real MU-Glioma estimation."""
    params = load_mu_glioma_params()
    rhos = [
        p["rho_per_day"] for p in params.values()
        if p.get("rho_per_day") is not None and p["rho_per_day"] > 0
    ]
    return sorted(rhos)


def train_cohort_policy(
    epochs: int = 250,
    lr: float = 2e-3,
    output_path: Optional[Path] = None,
) -> PolicyNetwork:
    """
    Train a single unified PolicyNetwork across the real MU-Glioma cohort distribution.
    """
    rhos = get_cohort_rho_distribution()
    stats = real_cohort_stats()
    
    print("=" * 70)
    print("PHASE 3: COHORT-LEVEL BEHAVIORAL CLONING TRAINING")
    print("=" * 70)
    print(f"Cohort size: {len(rhos)} positive patients")
    print(f"Rho range:  [{min(rhos):.6f}, {max(rhos):.6f}] day^-1")
    print(f"Rho median: {np.median(rhos):.6f} day^-1 (Cohort summary: {stats['rho_median']:.6f})")

    # Sample representative percentiles across the real cohort to ensure robust coverage
    percentiles = [5, 15, 25, 35, 50, 65, 75, 85, 95]
    representative_rhos = [float(np.percentile(rhos, p)) for p in percentiles]
    print(f"Representative rho percentiles: {[round(r, 5) for r in representative_rhos]}")

    # Generate demonstrations across the cohort distribution
    states, actions = generate_env_demonstrations(
        rho_list=representative_rhos,
        n_samples_per_rho=1000,
        seed=42,
    )
    print(f"Generated {len(states)} demonstration pairs across cohort heterogeneity")

    states_t = torch.FloatTensor(states)
    actions_t = torch.LongTensor(actions)

    policy = PolicyNetwork(obs_dim=5, n_actions=4, hidden=64)
    optimizer = torch.optim.Adam(policy.parameters(), lr=lr, weight_decay=1e-5)
    criterion = torch.nn.CrossEntropyLoss()

    batch_size = 256
    n_batches = int(np.ceil(len(states) / batch_size))

    for epoch in range(epochs):
        perm = torch.randperm(len(states_t))
        epoch_loss = 0.0

        for b in range(n_batches):
            idx = perm[b * batch_size : (b + 1) * batch_size]
            b_states = states_t[idx]
            b_actions = actions_t[idx]

            optimizer.zero_grad()
            logits = policy.net(b_states)
            loss = criterion(logits, b_actions)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()

        if (epoch + 1) % 50 == 0 or epoch == epochs - 1:
            with torch.no_grad():
                logits = policy.net(states_t)
                preds = logits.argmax(dim=-1)
                acc = (preds == actions_t).float().mean()
            print(f"  Epoch {epoch+1:3d}/{epochs}: loss={epoch_loss/n_batches:.4f}, acc={acc.item()*100:.2f}%")

    if output_path is None:
        output_path = PROJECT_ROOT / "output" / "bc_policy_cohort.pt"

    output_path.parent.mkdir(parents=True, exist_ok=True)
    torch.save(policy.state_dict(), str(output_path))
    print(f"\n[Saved] Cohort policy model -> {output_path}")

    # Also save bucketed models for stratified analysis if requested
    buckets = np.percentile(rhos, [0, 20, 40, 60, 80, 100])
    for i in range(5):
        bucket_rhos = [r for r in rhos if buckets[i] <= r <= buckets[i+1]]
        if bucket_rhos:
            bucket_mid = float(np.median(bucket_rhos))
            b_states, b_actions = generate_env_demonstrations(
                rho_list=[bucket_mid],
                n_samples_per_rho=800,
                seed=42 + i,
            )
            b_states_t = torch.FloatTensor(b_states)
            b_actions_t = torch.LongTensor(b_actions)
            b_policy = PolicyNetwork(obs_dim=5, n_actions=4, hidden=64)
            b_opt = torch.optim.Adam(b_policy.parameters(), lr=lr)
            for _ in range(100):
                b_opt.zero_grad()
                b_loss = criterion(b_policy.net(b_states_t), b_actions_t)
                b_loss.backward()
                b_opt.step()
            b_path = PROJECT_ROOT / "output" / f"bc_policy_bucket_{i}.pt"
            torch.save(b_policy.state_dict(), str(b_path))

    print(f"[Complete] Saved unified cohort policy and 5 bucketed models to output/")
    return policy


def main():
    parser = argparse.ArgumentParser(description="Train Cohort-Level Behavioral Cloning Policy")
    parser.add_argument("--epochs", type=int, default=250)
    parser.add_argument("--output", type=str, default="output/bc_policy_cohort.pt")
    args = parser.parse_args()

    out_path = Path(args.output)
    if not out_path.is_absolute():
        out_path = PROJECT_ROOT / out_path

    train_cohort_policy(epochs=args.epochs, output_path=out_path)


if __name__ == "__main__":
    main()
