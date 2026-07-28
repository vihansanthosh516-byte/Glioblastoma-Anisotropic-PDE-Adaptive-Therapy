import os
import sys
import numpy as np
import torch

# Ensure src is in import path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

def generate_dataset(num_samples=100, grid_size=32, steps=10, output_path="output/fno_dataset.pt"):
    """
    Generates training pairs for FNO.
    Inputs: u0 (initial tumor state), rho (proliferation)
    Targets: u_next (tumor state at t + dt * steps)
    """
    print(f"[Phase 10 DataGen] Generating {num_samples} PDE trajectory samples on a {grid_size}^3 grid...")
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    
    inputs = []
    targets = []
    
    for i in range(num_samples):
        # Sample random physical parameters
        rho = np.random.uniform(0.01, 0.05)
        
        # Random initial tumor sphere seed position
        cx, cy, cz = np.random.randint(8, grid_size - 8, size=3)
        grid = np.zeros((grid_size, grid_size, grid_size), dtype=np.float32)
        
        z, y, x = np.ogrid[:grid_size, :grid_size, :grid_size]
        mask = (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= np.random.uniform(3, 6)**2
        grid[mask] = np.random.uniform(0.5, 1.0)
        
        u0 = grid.copy()
        
        # Run short PDE trajectory forward (simplified exponential growth for data gen)
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
    
    torch.save({"x": x_data, "y": y_data}, output_path)
    print(f"[Phase 10 DataGen] Dataset saved to {output_path} (X: {x_data.shape}, Y: {y_data.shape})")

if __name__ == "__main__":
    generate_dataset(num_samples=100, grid_size=32)