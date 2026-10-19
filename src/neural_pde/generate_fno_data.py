"""Seeded synthetic exponential-growth pairs; targets do not solve a diffusion PDE.

Historical saved datasets are unchanged. New generated pairs are toy training data,
not measured longitudinal tumor states or evidence of PDE surrogate accuracy.
"""
import os
import sys
from pathlib import Path
import numpy as np
import torch

# Ensure src is in import path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

PROJECT_ROOT = Path(__file__).resolve().parents[2]


def generate_dataset(num_samples=100, grid_size=32, steps=10,
                     output_path=PROJECT_ROOT / "output" / "fno_dataset.pt", seed=42):
    """
    Generates training pairs for FNO.
    Inputs: u0 (initial tumor state), rho (proliferation)
    Targets: clipped exponential growth without diffusion (synthetic toy model).
    """
    output_path = Path(output_path)
    if not output_path.is_absolute():
        output_path = PROJECT_ROOT / output_path
    if output_path.exists():
        raise FileExistsError(f"Refusing to overwrite historical dataset: {output_path}")
    print(f"[Phase 10 DataGen] Generating {num_samples} synthetic exponential-growth samples on a {grid_size}^3 grid...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    
    inputs = []
    targets = []
    
    for i in range(num_samples):
        # Sample random physical parameters
        rho = rng.uniform(0.01, 0.05)
        
        # Random initial tumor sphere seed position
        cx, cy, cz = rng.integers(8, grid_size - 8, size=3)
        grid = np.zeros((grid_size, grid_size, grid_size), dtype=np.float32)
        
        z, y, x = np.ogrid[:grid_size, :grid_size, :grid_size]
        mask = (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= rng.uniform(3, 6)**2
        grid[mask] = rng.uniform(0.5, 1.0)
        
        u0 = grid.copy()
        
        # Toy reaction-only target; no spatial diffusion or PDE solver is used.
        u_next = u0 * np.exp(rho * steps * 0.5)
        u_next = np.clip(u_next, 0.0, 1.0)
        
        # Pack input features: [u0, rho_field]
        rho_field = np.full_like(u0, rho)
        inp_tensor = np.stack([u0, rho_field], axis=0)  # Shape: (2, grid_size, grid_size, grid_size)
        
        inputs.append(inp_tensor)
        targets.append(u_next[np.newaxis, ...])  # Shape: (1, grid_size, grid_size, grid_size)
        
        if (i + 1) % 25 == 0 or (i + 1) == num_samples:
            print(f"[Phase 10 DataGen] Generated {i + 1}/{num_samples} samples")

    x_data = torch.tensor(np.array(inputs), dtype=torch.float32)
    y_data = torch.tensor(np.array(targets), dtype=torch.float32)
    
    torch.save({"x": x_data, "y": y_data, "metadata": {
        "seed": seed, "synthetic": True, "grid_size": grid_size, "steps": steps,
        "target_method": "clipped exponential growth; no diffusion PDE",
    }}, output_path)
    if output_path.stat().st_size == 0:
        raise RuntimeError(f"Empty dataset: {output_path}")
    print(f"[Phase 10 DataGen] Dataset saved to {output_path} (X: {x_data.shape}, Y: {y_data.shape})")

if __name__ == "__main__":
    generate_dataset(num_samples=100, grid_size=32)
