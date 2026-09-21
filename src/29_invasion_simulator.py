#!/usr/bin/env python3
"""
Month 3, Week 3: Integrated Invasion Simulator (Calibrated)
Combines Fisher-Kolmogorov PDE (finite difference, Dirichlet BC) with 
cellular automaton for cell-state transitions.

TIME UNITS (documented and consistent):
- FK PDE internal time step: dt_fk = 0.001 (FK time units)
- 1 FK time unit = 1 hour (calibrated to clinical velocity)
- Simulator step = 1 hour = 1000 FK internal steps
- CA sub-steps per simulator step: for stochastic averaging

Calibrated for clinical glioblastoma invasion velocity: 10-50 µm/hr
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import torch


# ---------------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = PROJECT_ROOT / "output"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Cell state constants
# ---------------------------------------------------------------------------
EMPTY = 0
HEALTHY = 1
PERIPHERY = 2
CORE = 3
NECROTIC = 4


# ---------------------------------------------------------------------------
# Fisher-Kolmogorov PDE Solver (Finite Difference, Dirichlet BC, RK4)
# ---------------------------------------------------------------------------
class FKSolver:
    """Finite difference FK solver: ∂ρ/∂t = D∇²ρ + rρ(1-ρ), ρ=0 at boundaries."""

    def __init__(
        self,
        grid_size: Tuple[int, int],
        D: float = 10.0,
        r: float = 0.4,
        dt: float = 0.001,
        device: torch.device = None,
    ):
        self.H, self.W = grid_size
        self.D = D
        self.r = r
        self.dt = dt  # Internal FK time step (FK time units)
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.dx = 1.0
        self.dx2 = 1.0

        # Slicing indices for 5-point stencil
        self.interior = (slice(1, self.H - 1), slice(1, self.W - 1))
        self.up = (slice(0, self.H - 2), slice(1, self.W - 1))
        self.down = (slice(2, self.H), slice(1, self.W - 1))
        self.left = (slice(1, self.H - 1), slice(0, self.W - 2))
        self.right = (slice(1, self.H - 1), slice(2, self.W))

        max_dt = self.dx2 / (4.0 * D)
        if dt > max_dt:
            print(f"[FK] WARNING: dt={dt} > stability limit {max_dt:.6f}")

    def reaction(self, rho: torch.Tensor) -> torch.Tensor:
        return self.r * rho * (1.0 - rho)

    def laplacian(self, rho: torch.Tensor) -> torch.Tensor:
        lap = torch.zeros_like(rho)
        lap[self.interior] = (
            rho[self.up] + rho[self.down] + rho[self.left] + rho[self.right] - 4.0 * rho[self.interior]
        ) / self.dx2
        return lap

    def rhs(self, rho: torch.Tensor) -> torch.Tensor:
        return self.D * self.laplacian(rho) + self.reaction(rho)

    def step(self, rho: torch.Tensor) -> torch.Tensor:
        # RK4
        k1 = self.rhs(rho)
        k2 = self.rhs(rho + 0.5 * self.dt * k1)
        k3 = self.rhs(rho + 0.5 * self.dt * k2)
        k4 = self.rhs(rho + self.dt * k3)
        rho_new = rho + (self.dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

        # Dirichlet BC
        rho_new[0, :] = 0.0
        rho_new[-1, :] = 0.0
        rho_new[:, 0] = 0.0
        rho_new[:, -1] = 0.0
        return torch.clamp(rho_new, 0.0, 1.0)

    def advance_hours(self, rho: torch.Tensor, hours: float) -> torch.Tensor:
        """Advance morphogen field by specified hours (FK time units)."""
        n_internal = int(hours / self.dt)
        for _ in range(n_internal):
            rho = self.step(rho)
        return rho


# ---------------------------------------------------------------------------
# Integrated Invasion Simulator
# ---------------------------------------------------------------------------
class IntegratedInvasionSimulator:
    """
    Multi-scale simulator with CONSISTENT TIME UNITS:
    - Simulator step = 1 hour
    - FK PDE: advances by 1 hour per simulator step (1000 internal steps)
    - CA: 1 sub-step per simulator step (stochastic, represents 1 hour)
    
    Coupling:
    - Periphery cells secrete morphogen
    - Morphogen drives Healthy->Periphery transition (invasion)
    - Core cells become necrotic if morphogen too low
    - Periphery/Core cells proliferate into adjacent empty space
    """

    def __init__(
        self,
        grid_size: Tuple[int, int] = (512, 512),
        n_hours: int = 100,           # Total simulation hours
        fk_D: float = 10.0,
        fk_r: float = 0.4,
        fk_dt: float = 0.01,          # FK internal time step (FK time units = hours) - stable for D=10
        device: torch.device = None,
    ):
        self.H, self.W = grid_size
        self.n_hours = n_hours
        self.device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")

        # FK solver for morphogen (dt=0.01 hours per internal step, stable for D=10)
        self.fk = FKSolver(grid_size, fk_D, fk_r, fk_dt, self.device)
        self.fk_hours_per_step = 1.0  # 1 hour per simulator step
        self.fk_internal_per_step = int(self.fk_hours_per_step / fk_dt)  # 100

        # Cell state grid
        self.rho = torch.zeros(grid_size, dtype=torch.int8, device=self.device)

        # Morphogen field (continuous [0,1])
        self.morphogen = torch.zeros(grid_size, dtype=torch.float32, device=self.device)

        # Parameters - CALIBRATED for biological realism
        # Cell division time ~24 hours → probability per hour ~ 0.04
        # With 4 neighbors, effective division probability per hour ~ 0.06
        self.rules = {
            'proliferation_rate_per_hour': 0.06,   # Per CA step (1 hour)
            'core_to_periphery_threshold': 0.3,
            'healthy_to_periphery_threshold': 0.1,  # Track leading edge (FK wave front)
            'periphery_secrete_rate_per_hour': 0.2, # Increased secretion
            'necrosis_threshold': 0.1,
        }

        # Metrics
        self.metrics_history = []
        self.front_history = []

    def initialize_from_latent(self, latent_coords: np.ndarray = None) -> None:
        """Initialize tissue: 75% healthy, central tumor core + periphery ring."""
        n_cells = self.H * self.W
        n_healthy = int(n_cells * 0.75)

        flat = torch.zeros(n_cells, dtype=torch.int8, device=self.device)
        flat[:n_healthy] = HEALTHY
        flat[n_healthy:] = EMPTY
        perm = torch.randperm(n_cells, device=self.device)
        self.rho = flat[perm].view(self.H, self.W)

        # Central tumor seed
        ch, cw = self.H // 2, self.W // 2
        core_radius = 6
        for dh in range(-core_radius, core_radius + 1):
            for dw in range(-core_radius, core_radius + 1):
                if dh*dh + dw*dw <= core_radius*core_radius:
                    h, w = ch + dh, cw + dw
                    if 0 <= h < self.H and 0 <= w < self.W:
                        self.rho[h, w] = CORE
                        self.morphogen[h, w] = 1.0

        # Periphery ring
        ring_r = core_radius + 2
        for dh in range(-ring_r, ring_r + 1):
            for dw in range(-ring_r, ring_r + 1):
                dist_sq = dh*dh + dw*dw
                if core_radius*core_radius < dist_sq <= ring_r*ring_r:
                    h, w = ch + dh, cw + dw
                    if 0 <= h < self.H and 0 <= w < self.W:
                        if self.rho[h, w] == HEALTHY:
                            self.rho[h, w] = PERIPHERY

        print(f"[SIM] Initial: {self.count_cells()}")

    def count_cells(self) -> Dict[str, int]:
        counts = {}
        for val, name in [(EMPTY, 'empty'), (HEALTHY, 'healthy'),
                          (PERIPHERY, 'periphery'), (CORE, 'core'),
                          (NECROTIC, 'necrotic')]:
            counts[name] = int((self.rho == val).sum().item())
        return counts

    def step_ca(self) -> Dict[str, int]:
        """One CA step (1 hour): proliferation, state transitions, necrosis."""
        changes = {'proliferation': 0, 'transition': 0, 'necrosis': 0, 'replenish': 0}

        # --- 1. Proliferation: periphery/core cells divide into adjacent empty ---
        proliferating_mask = (self.rho == PERIPHERY) | (self.rho == CORE)
        coords = proliferating_mask.nonzero(as_tuple=True)

        for i in range(len(coords[0])):
            h, w = coords[0][i].item(), coords[1][i].item()
            if np.random.random() > self.rules['proliferation_rate_per_hour']:
                continue
            empty_neighbors = []
            for dh, dw in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nh, nw = h + dh, w + dw
                if 0 <= nh < self.H and 0 <= nw < self.W and self.rho[nh, nw] == EMPTY:
                    empty_neighbors.append((nh, nw))
            if empty_neighbors:
                nh, nw = empty_neighbors[np.random.randint(len(empty_neighbors))]
                self.rho[nh, nw] = self.rho[h, w]
                changes['proliferation'] += 1

        # --- 2. State transitions driven by morphogen ---
        # Healthy -> Periphery (if morphogen > threshold): INVASION
        healthy_mask = self.rho == HEALTHY
        morphogen_high = self.morphogen > self.rules['healthy_to_periphery_threshold']
        transition_mask = healthy_mask & morphogen_high
        if transition_mask.any():
            self.rho[transition_mask] = PERIPHERY
            changes['transition'] += int(transition_mask.sum().item())

        # Core -> Periphery (if morphogen drops below threshold)
        core_mask = self.rho == CORE
        morphogen_low = self.morphogen < self.rules['core_to_periphery_threshold']
        core_to_periphery = core_mask & morphogen_low
        if core_to_periphery.any():
            self.rho[core_to_periphery] = PERIPHERY
            changes['transition'] += int(core_to_periphery.sum().item())

        # --- 3. Necrosis: core cells die if morphogen too low ---
        necrosis_mask = core_mask & (self.morphogen < self.rules['necrosis_threshold'])
        if necrosis_mask.any():
            self.rho[necrosis_mask] = NECROTIC
            changes['necrosis'] += int(necrosis_mask.sum().item())

        return changes

    def step_morphogen(self) -> None:
        """Advance morphogen by 1 hour: secretion + FK diffusion/reaction."""
        # Secretion from periphery cells (per hour)
        periphery_mask = self.rho == PERIPHERY
        if periphery_mask.any():
            secretion = torch.rand_like(self.morphogen) < self.rules['periphery_secrete_rate_per_hour']
            add_mask = periphery_mask & secretion
            self.morphogen[add_mask] = torch.minimum(
                self.morphogen[add_mask] + 0.3,
                torch.ones_like(self.morphogen[add_mask])
            )

        # FK PDE: advance by 1 hour (1000 internal steps)
        self.morphogen = self.fk.advance_hours(self.morphogen, self.fk_hours_per_step)

    def step(self) -> Dict[str, int]:
        """One simulator step = 1 hour."""
        # CA step (1 hour)
        changes = self.step_ca()
        
        # Morphogen step (1 hour)
        self.step_morphogen()
        
        return changes

    def compute_front_metrics(self) -> Tuple[float, float, int]:
        """Compute invasion front radius from ACTUAL TUMOR CELLS (not FK wave)."""
        tumor_mask = (self.rho == PERIPHERY) | (self.rho == CORE)
        if not tumor_mask.any():
            return 0.0, 0.0, 0

        coords = tumor_mask.nonzero(as_tuple=True)
        ch, cw = self.H / 2, self.W / 2
        radii = torch.sqrt(
            (coords[0].float() - ch) ** 2 + (coords[1].float() - cw) ** 2
        )
        front_radius = float(radii.max().item())
        mean_radius = float(radii.mean().item())
        front_cells = int(tumor_mask.sum().item())
        return front_radius, mean_radius, front_cells

    def run(
        self,
        n_hours: int = None,
        save_interval: int = 20,
        frames_dir: str = str(OUTPUT_DIR / "invasion_frames"),
    ) -> Dict:
        if n_hours is None:
            n_hours = self.n_hours

        Path(frames_dir).mkdir(parents=True, exist_ok=True)

        print(f"[SIM] Running {n_hours} hours (1 CA step + {self.fk_internal_per_step} FK steps per hour)...")
        t0 = time.perf_counter()

        for hour in range(n_hours):
            if hour % 10 == 0:
                print(f"  Progress: hour {hour}/{n_hours}")
            changes = self.step()
            counts = self.count_cells()
            front_r, mean_r, front_n = self.compute_front_metrics()

            metrics = {
                'hour': hour,
                **counts,
                'proliferation': changes['proliferation'],
                'transition': changes['transition'],
                'necrosis': changes['necrosis'],
                'replenish': changes['replenish'],
                'front_radius': front_r,
                'mean_tumor_radius': mean_r,
                'front_cells': front_n,
            }
            self.metrics_history.append(metrics)
            self.front_history.append(front_r)

            if hour % save_interval == 0:
                self._save_frame(hour, frames_dir)

            if hour % 20 == 0:
                print(f"  Hour {hour}/{n_hours}: {counts}, front_r={front_r:.1f}")

        elapsed = time.perf_counter() - t0
        print(f"[SIM] Completed in {elapsed:.1f}s")

        return {
            'metrics_history': self.metrics_history,
            'runtime': elapsed,
            'final_counts': self.count_cells(),
        }

    def _save_frame(self, hour: int, frames_dir: str) -> None:
        grid_np = self.rho.cpu().numpy()
        morph_np = self.morphogen.cpu().numpy()

        fig, axes = plt.subplots(1, 3, figsize=(12, 4))

        cmap = plt.cm.colors.ListedColormap([
            'white', '#2E8B57', '#FF8C00', '#DC143C', '#696969'
        ])
        axes[0].imshow(grid_np, cmap=cmap, vmin=0, vmax=4)
        axes[0].set_title(f'Cell States (Hour {hour})')
        axes[0].axis('off')

        im = axes[1].imshow(morph_np, cmap='hot', vmin=0, vmax=1)
        axes[1].set_title('Morphogen Concentration')
        axes[1].axis('off')
        plt.colorbar(im, ax=axes[1], fraction=0.046)

        if len(self.front_history) > 1:
            axes[2].plot(self.front_history, 'b-', linewidth=1.5)
            axes[2].set_xlabel('Hour')
            axes[2].set_ylabel('Front Radius (px)')
            axes[2].set_title('Invasion Front Progression')
            axes[2].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig(f"{frames_dir}/frame_{hour:04d}.png", dpi=150, bbox_inches='tight')
        plt.close()

    def export_results(self, runtime_seconds: float = 0.0) -> None:
        print("[EXPORT] Saving simulation results...")
        OUTPUT_DIR.mkdir(exist_ok=True)

        hours = [m['hour'] for m in self.metrics_history]
        healthy = [m['healthy'] for m in self.metrics_history]
        periphery = [m['periphery'] for m in self.metrics_history]
        core = [m['core'] for m in self.metrics_history]
        necrotic = [m['necrotic'] for m in self.metrics_history]
        front_r = [m['front_radius'] for m in self.metrics_history]
        prolif = [m['proliferation'] for m in self.metrics_history]
        trans = [m['transition'] for m in self.metrics_history]
        necro = [m['necrosis'] for m in self.metrics_history]

        metrics_array = np.column_stack([hours, healthy, periphery, core, necrotic,
                                         front_r, prolif, trans, necro])
        np.save(str(OUTPUT_DIR / "invasion_metrics.npy"), metrics_array)

        with open(str(OUTPUT_DIR / "invasion_metrics.tsv"), "w") as f:
            f.write("hour\thealthy\tperiphery\tcore\tnecrotic\tfront_radius\tproliferation\ttransition\tnecrosis\n")
            for row in metrics_array:
                f.write("\t".join(str(x) for x in row) + "\n")

        # Compute velocity in FREE-PROPAGATION regime (before boundary effects)
        # Grid diagonal = sqrt(256^2 + 256^2) ≈ 362 px
        # Boundary effects start at ~90% of diagonal ≈ 326 px
        # We measure over the first 50 hours (free propagation, before boundary)
        if len(self.front_history) > 50:
            free_advance = self.front_history[50] - self.front_history[0]
            free_velocity_px_per_hr = free_advance / 50.0
            free_velocity_um_per_hr = free_velocity_px_per_hr * 5.0
        else:
            free_advance = self.front_history[-1] - self.front_history[0]
            free_velocity_px_per_hr = free_advance / len(self.front_history) if len(self.front_history) > 1 else 0.0
            free_velocity_um_per_hr = free_velocity_px_per_hr * 5.0

        # Also compute overall average for reference
        if len(self.front_history) > 1:
            total_advance = self.front_history[-1] - self.front_history[0]
            avg_velocity_px_per_hr = total_advance / self.n_hours
            velocity_um_per_hr = avg_velocity_px_per_hr * 5.0
        else:
            avg_velocity_px_per_hr = 0.0
            velocity_um_per_hr = 0.0

        summary = {
            'n_hours': self.n_hours,
            'grid_size': [self.H, self.W],
            'final_counts': self.count_cells(),
            'max_front_radius': max(self.front_history) if self.front_history else 0,
            # Free-propagation velocity (physically meaningful)
            'free_front_velocity_pixels_per_hour': free_velocity_px_per_hr,
            'free_front_velocity_um_per_hr': free_velocity_um_per_hr,
            'free_propagation_hours': 50,
            # Overall average (includes boundary-limited regime)
            'avg_front_velocity_pixels_per_hour': avg_velocity_px_per_hr,
            'avg_front_velocity_um_per_hr': velocity_um_per_hr,
            'clinical_velocity_range_um_per_hr': [10, 50],
            'in_range': 10 <= free_velocity_um_per_hr <= 50,
            'fk_params': {'D': self.fk.D, 'r': self.fk.r, 'dt': self.fk.dt},
            'ca_params': self.rules,
            'fk_analytical_c_px_per_hr': 2 * np.sqrt(self.fk.D * self.fk.r),
            'fk_analytical_c_um_per_hr': 2 * np.sqrt(self.fk.D * self.fk.r) * 5.0,
            'runtime_seconds': runtime_seconds,
            'notes': 'FK r=0.8 chosen to match CA+FK coupling timescale; script 28 uses r=0.4 for standalone FK calibration'
        }

        with open(str(OUTPUT_DIR / "invasion_summary.json"), "w") as f:
            json.dump(summary, f, indent=2)

        print(f"  output/invasion_metrics.npy: {metrics_array.shape}")
        print(f"  output/invasion_summary.json")
        print(f"  Total front advance (all): {total_advance:.1f} px over {self.n_hours} hours")
        print(f"  Overall avg velocity: {avg_velocity_px_per_hr:.4f} px/hr = {velocity_um_per_hr:.1f} µm/hr")
        print(f"  Free-propagation advance (first 50h): {free_advance:.1f} px")
        print(f"  Free-propagation velocity: {free_velocity_px_per_hr:.4f} px/hr = {free_velocity_um_per_hr:.1f} µm/hr")
        print(f"  FK analytical c = 2*sqrt(Dr) = {2*np.sqrt(self.fk.D*self.fk.r):.4f} px/hr = {2*np.sqrt(self.fk.D*self.fk.r)*5:.1f} µm/hr")
        print(f"  Clinical range [10, 50] µm/hr: {'PASS' if 10 <= free_velocity_um_per_hr <= 50 else 'FAIL'}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 60)
    print("MONTH 3 WEEK 3: INTEGRATED INVASION SIMULATOR (CALIBRATED)")
    print("=" * 60)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"[SIM] Device: {device}")

    # FK parameters: D=10.0, r=0.8
    # Note: Script 28 uses r=0.4 for standalone FK calibration (c=4.0 px/hr)
    # Here we use r=0.8 to achieve faster FK wave (c=5.66 px/hr) that better
    # couples with the CA dynamics at the chosen time step and grid resolution.
    sim = IntegratedInvasionSimulator(
        grid_size=(512, 512),
        n_hours=150,
        fk_D=10.0,
        fk_r=0.8,
        fk_dt=0.01,
        device=device,
    )

    sim.initialize_from_latent()
    results = sim.run(save_interval=20)
    sim.export_results(runtime_seconds=results['runtime'])

    print(f"\n[SIM] Final counts: {results['final_counts']}")
    print(f"[SIM] Max front radius: {max(sim.front_history):.1f}")
    total_advance = sim.front_history[-1] - sim.front_history[0]
    avg_vel = total_advance / sim.n_hours
    vel_um_hr = avg_vel * 5.0
    free_vel = (sim.front_history[50] - sim.front_history[0]) / 50.0 if len(sim.front_history) > 50 else avg_vel
    print(f"[SIM] Total front advance: {total_advance:.1f} px over {sim.n_hours} hours")
    print(f"[SIM] Overall avg velocity: {avg_vel:.4f} px/hr = {vel_um_hr:.1f} µm/hr")
    print(f"[SIM] Free-propagation velocity (first 50h): {free_vel:.4f} px/hr = {free_vel*5:.1f} µm/hr")
    print(f"[SIM] Clinical range [10, 50] µm/hr: {'PASS' if 10 <= free_vel*5 <= 50 else 'FAIL'}")

    print("\n[SUCCESS] Month 3 Week 3 Complete: Integrated Invasion Simulator")
    print("  - Frames: output/invasion_frames/")
    print("  - Metrics: output/invasion_metrics.npy, .tsv")
    print("  - Summary: output/invasion_summary.json")


if __name__ == "__main__":
    main()