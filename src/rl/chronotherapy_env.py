#!/usr/bin/env python3
"""
Phase 12: Closed-Loop Chronotherapy RL Environment
===================================================
Gymnasium environment for real-time adaptive therapy with:
- FNO-accelerated 3D tumor growth PDE
- Virtual biosensor suite (ctDNA, IFP, BBB, pH, O2)
- Circadian-aware micro-dosing control
- Safety constraints and toxicity tracking
"""
import gymnasium as gym
from gymnasium import spaces
import numpy as np
import torch
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent))

from neural_pde.fno_solver import FNO3d
from sensing.virtual_sensor import VirtualBiosensorSuite, SensorConfig, CircadianBiosensorSuite


class ChronotherapyEnv(gym.Env):
    """
    Closed-loop chronotherapy environment for real-time adaptive therapy.
    
    State: [tumor_volume, ctdna, ifp, bbb, ph, pO2, cum_toxicity, 
            circadian_phase, drug_concentrations, time_since_last_dose]
    Action: [dose_TMZ, dose_inhibitor, dose_radiation, timing_offset]
    
    Reward: -tumor_volume - toxicity_penalty - resistance_penalty + shrinkage_bonus
    """
    
    metadata = {"render_modes": ["human", "rgb_array"], "render_fps": 4}
    
    def __init__(
        self,
        grid_size: int = 64,
        dt_hours: float = 1.0,
        max_episode_hours: int = 240,  # 10 days
        fno_model_path: str = "output/fno_model.pth",
        patient_id: str = "BraTS2021_00000",
        circadian: bool = True,
        seed: int = None,
    ):
        super().__init__()
        
        self.grid_size = grid_size
        self.dt = dt_hours  # hours per RL step
        self.max_steps = int(max_episode_hours / dt_hours)
        self.circadian = circadian
        self.patient_id = patient_id
        
        # Set random seed
        if seed is not None:
            np.random.seed(seed)
            torch.manual_seed(seed)
        
        # Load FNO model
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.fno = FNO3d(in_channels=2, out_channels=1, modes=8, width=20).to(self.device)
        
        try:
            checkpoint = torch.load(fno_model_path, map_location=self.device)
            if "model_state" in checkpoint:
                state_dict = checkpoint["model_state"]
            else:
                state_dict = checkpoint
            self.fno.load_state_dict(state_dict)
            self.fno.eval()
            print(f"[ChronotherapyEnv] FNO loaded from {fno_model_path}")
        except Exception as e:
            print(f"[ChronotherapyEnv] Warning: FNO load failed: {e}. Using untrained FNO.")
        
        # Initialize tumor state
        self._init_tumor_state()
        
        # Initialize biosensors
        if circadian:
            self.sensors = CircadianBiosensorSuite()
        else:
            self.sensors = VirtualBiosensorSuite()
        
        # Drug PK/PD parameters
        self.pk_params = {
            "TMZ": {"half_life": 1.8, "kill_rate": 0.08, "toxicity_per_dose": 0.02},
            "Inhibitor": {"half_life": 12.0, "kill_rate": 0.05, "toxicity_per_dose": 0.01},
            "Radiation": {"half_life": 0.1, "kill_rate": 0.15, "toxicity_per_dose": 0.05},
        }
        
        # Drug concentrations (normalized 0-1)
        self.drug_concs = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        self.last_dose_time = {"TMZ": -100, "Inhibitor": -100, "Radiation": -100}
        
        # Cumulative toxicity
        self.cum_toxicity = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        self.max_toxicity = {"TMZ": 1.0, "Inhibitor": 0.8, "Radiation": 0.6}
        
        # Resistance tracking
        self.resistant_fraction = 0.01  # Initial resistant fraction
        self.resistance_growth_rate = 0.001
        
        # Define action space: [TMZ_dose, Inhibitor_dose, Radiation_dose, timing_offset_hours]
        # All in [0, 1], normalized
        self.action_space = spaces.Box(
            low=np.array([0.0, 0.0, 0.0, -12.0]),
            high=np.array([1.0, 1.0, 1.0, 12.0]),
            dtype=np.float32,
        )
        
        # Define observation space
        # [tumor_vol_norm, ctdna_norm, ifp_norm, bbb_norm, ph_norm, o2_norm,
        #  cum_tox_TMZ, cum_tox_inh, cum_tox_rad,
        #  circ_phase_sin, circ_phase_cos,
        #  conc_TMZ, conc_inh, conc_rad,
        #  time_since_TMZ, time_since_inh, time_since_rad,
        #  resistant_frac]
        self.observation_space = spaces.Box(
            low=-np.inf, high=np.inf, shape=(18,), dtype=np.float32
        )
        
        # Episode tracking
        self.step_count = 0
        self.episode_history = []
        
        # Render
        self.render_mode = None
    
    def _init_tumor_state(self):
        """Initialize tumor with small spherical seed."""
        gs = self.grid_size
        self.u = np.zeros((gs, gs, gs), dtype=np.float32)
        self.rho = np.full((gs, gs, gs), 0.02, dtype=np.float32)
        
        # Spherical seed at center with slight noise
        cx, cy, cz = self.grid_size // 2, self.grid_size // 2, self.grid_size // 2
        z, y, x = np.ogrid[:gs, :gs, :gs]
        mask = (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= 5**2
        self.u[mask] = np.random.uniform(0.5, 1.0, size=mask.sum())
        self.rho[mask] = np.random.uniform(0.015, 0.03, size=mask.sum())
        
        # Diffusion tensor (synthetic DTI-like)
        self.D_tensor = self._create_synthetic_tensor()
    
    def _create_synthetic_tensor(self):
        """Create synthetic DTI tensor field with corpus callosum tract."""
        gs = self.grid_size
        D = np.zeros((3, 3, gs, gs, gs), dtype=np.float32)
        
        # Corpus callosum: horizontal tract through center
        # Use mgrid for full 3D coordinate grids
        z, y, x = np.mgrid[:gs, :gs, :gs]
        cy = gs // 2
        dist_to_tract = np.sqrt((y - cy)**2 + (z - gs//2)**2)
        in_tract = dist_to_tract < gs // 6
        
        assert in_tract.shape == (gs, gs, gs), f"in_tract shape {in_tract.shape} != ({gs}, {gs}, {gs})"
        
        # Base isotropic diffusivity (gray matter)
        D_iso = 0.0013
        D_parallel = 0.013
        D_perp = 0.0013
        
        for i in range(3):
            for j in range(3):
                D[i, j] = D_iso
        
        # Tract direction: x-axis (left-right)
        n = np.array([1.0, 0.0, 0.0])
        for i in range(3):
            for j in range(3):
                D[i, j] = np.where(in_tract, 
                                   D_perp + (D_parallel - D_perp) * n[i] * n[j],
                                   D[i, j])
        
        return D
    
    def _fno_step(self, u0: np.ndarray, rho: float, dt: float) -> np.ndarray:
        """Step tumor state forward using FNO."""
        gs = self.grid_size
        
        # Prepare input: [u0, rho_field]
        rho_field = np.full_like(self.u, rho)
        inp = np.stack([self.u, rho_field], axis=0).astype(np.float32)  # (2, gs, gs, gs)
        inp_tensor = torch.tensor(inp, dtype=torch.float32).unsqueeze(0).to(self.device)
        
        with torch.no_grad():
            out = self.fno(inp_tensor)  # (1, 1, gs, gs, gs)
        
        u_next = out.squeeze().cpu().numpy()
        u_next = np.clip(u_next, 0.0, 1.0)
        
        return u_next
    
    def _pk_step(self, drug: str, dose: float, dt: float):
        """Update drug concentration using PK model."""
        params = self.pk_params[drug]
        k_elim = np.log(2) / params["half_life"]
        
        # Elimination
        self.drug_concs[drug] *= np.exp(-k_elim * dt)
        
        # Add dose if given
        if dose > 0:
            self.drug_concs[drug] = min(1.0, self.drug_concs[drug] + dose)
            self.last_dose_time[drug] = self.step_count
    
    def _compute_reward(self, prev_vol: float, curr_vol: float) -> float:
        """Compute reward: tumor shrinkage - toxicity - resistance."""
        # Volume change
        vol_change = (prev_vol - curr_vol) / max(prev_vol, 1e-6)
        
        # Toxicity penalty
        tox_penalty = sum(self.cum_toxicity.values()) / sum(self.max_toxicity.values())
        
        # Resistance penalty
        resistance_penalty = self.resistant_fraction * 10
        
        # Shrinkage bonus
        shrinkage_bonus = max(0, vol_change) * 10
        
        # Volume penalty
        vol_penalty = curr_vol / 100000  # normalize
        
        reward = shrinkage_bonus - vol_penalty - 5 * tox_penalty - resistance_penalty
        
        return float(reward)
    
    def _get_obs(self) -> np.ndarray:
        """Construct observation vector."""
        # Normalize volume
        vol = np.sum(self.u > 0.1) * 1.0  # mm^3
        vol_norm = vol / 100000.0
        
        # Sensor readings (already normalized in 0-1 range)
        sensor_readings = self.sensors.read_all(self.u, self.rho, t=self.step_count * self.dt)
        
        ctdna_norm = sensor_readings['ctdna_fragments_per_ml'] / 1000.0
        ifp_norm = sensor_readings['ifp_mmHg'] / 50.0
        bbb_norm = sensor_readings['bbb_permeability']
        ph_norm = (sensor_readings['pH'] - 7.0) / 1.0  # centered
        o2_norm = sensor_readings['pO2_mmHg'] / 100.0
        
        # Toxicity
        tox_norm = np.array([
            self.cum_toxicity["TMZ"] / self.max_toxicity["TMZ"],
            self.cum_toxicity["Inhibitor"] / self.max_toxicity["Inhibitor"],
            self.cum_toxicity["Radiation"] / self.max_toxicity["Radiation"],
        ])
        
        # Circadian phase
        t_hours = self.step_count * self.dt
        circ_sin = np.sin(2 * np.pi * t_hours / 24)
        circ_cos = np.cos(2 * np.pi * t_hours / 24)
        
        # Drug concentrations
        concs = np.array([
            self.drug_concs["TMZ"],
            self.drug_concs["Inhibitor"],
            self.drug_concs["Radiation"],
        ])
        
        # Time since last dose
        time_since = np.array([
            (self.step_count - self.last_dose_time["TMZ"]) * self.dt,
            (self.step_count - self.last_dose_time["Inhibitor"]) * self.dt,
            (self.step_count - self.last_dose_time["Radiation"]) * self.dt,
        ])
        time_since_norm = np.clip(time_since / 48.0, 0, 1)
        
        # Resistant fraction
        res_norm = self.resistant_fraction
        
        obs = np.concatenate([
            [vol_norm], [ctdna_norm], [ifp_norm], [bbb_norm], [ph_norm], [o2_norm],
            tox_norm,
            [circ_sin], [circ_cos],
            concs,
            time_since_norm,
            [res_norm],
        ]).astype(np.float32)
        
        return obs
    
    def reset(self, seed=None, options=None):
        super().reset(seed=seed)
        
        self._init_tumor_state()
        self.drug_concs = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        self.last_dose_time = {"TMZ": -100, "Inhibitor": -100, "Radiation": -100}
        self.cum_toxicity = {"TMZ": 0.0, "Inhibitor": 0.0, "Radiation": 0.0}
        self.resistant_fraction = 0.01
        self.step_count = 0
        self.episode_history = []
        
        return self._get_obs(), {}
    
    def step(self, action):
        """Execute one RL step."""
        # Parse action
        dose_tmz = float(action[0])
        dose_inh = float(action[1])
        dose_rad = float(action[2])
        timing_offset = float(action[3])  # hours offset from current time
        
        # Apply timing offset (if within bounds)
        if abs(timing_offset) < 0.5:
            timing_offset = 0  # execute now
        
        prev_vol = np.sum(self.u > 0.1)
        
        # Update drug PK
        dt = self.dt
        for drug, dose in [("TMZ", dose_tmz), ("Inhibitor", dose_inh), ("Radiation", dose_rad)]:
            self._pk_step(drug, dose, dt)
            if dose > 0:
                self.cum_toxicity[drug] += self.pk_params[drug]["toxicity_per_dose"] * dose
        
        # Update drug concentrations array
        self.drug_concs = {
            "TMZ": self.drug_concs["TMZ"],
            "Inhibitor": self.drug_concs["Inhibitor"],
            "Radiation": self.drug_concs["Radiation"],
        }
        
        # FNO PDE step
        rho_eff = self.rho.mean() * (1 - self.resistant_fraction)
        self.u = self._fno_step(self.u, rho_eff, dt)
        
        # Apply drug killing
        kill_tmz = self.pk_params["TMZ"]["kill_rate"] * self.drug_concs["TMZ"]
        kill_inh = self.pk_params["Inhibitor"]["kill_rate"] * self.drug_concs["Inhibitor"]
        kill_rad = self.pk_params["Radiation"]["kill_rate"] * self.drug_concs["Radiation"]
        
        total_kill = kill_tmz + kill_inh + kill_rad
        if total_kill > 0:
            self.u = self.u * (1 - total_kill)
            self.u = np.clip(self.u, 0.0, 1.0)
            
            # Resistance evolution
            self.resistant_fraction += self.resistance_growth_rate * (1 - self.resistant_fraction)
            self.resistant_fraction = min(self.resistant_fraction, 0.5)
        
        # Circadian phase update
        t_hours = self.step_count * self.dt
        
        # Compute reward
        curr_vol = np.sum(self.u > 0.1)
        reward = self._compute_reward(prev_vol, curr_vol)
        
        # Update step count
        self.step_count += 1
        
        # Check termination
        terminated = False
        truncated = False
        
        # Tumor cleared
        if curr_vol < 10:
            terminated = True
            reward += 100  # Clearance bonus
        
        # Toxicity limit exceeded
        if any(self.cum_toxicity[d] >= self.max_toxicity[d] for d in self.cum_toxicity):
            terminated = True
            reward -= 50
        
        # Max steps
        if self.step_count >= self.max_steps:
            truncated = True
        
        # Tumor exploded
        if curr_vol > 200000:
            terminated = True
            reward -= 100
        
        # Record history
        self.episode_history.append({
            "step": self.step_count,
            "volume": curr_vol,
            "reward": reward,
            "action": action,
            "toxicity": dict(self.cum_toxicity),
        })
        
        return self._get_obs(), reward, terminated, truncated, {}
    
    def render(self):
        pass
    
    def close(self):
        pass


# Quick validation
if __name__ == "__main__":
    print("=" * 60)
    print("Phase 12: Closed-Loop Chronotherapy RL Environment Test")
    print("=" * 60)
    
    env = ChronotherapyEnv(grid_size=32, max_episode_hours=48, dt_hours=1.0)
    
    print(f"\n[Test] Observation space: {env.observation_space.shape}")
    print(f"[Test] Action space: {env.action_space.shape}")
    
    # Test reset
    obs, info = env.reset()
    print(f"[Test] Reset obs shape: {obs.shape}, values: {obs[:6]}")
    
    # Test a few steps
    print("\n[Test] Running 10 steps...")
    total_reward = 0
    for i in range(10):
        action = np.array([0.1, 0.0, 0.0, 0.0], dtype=np.float32)  # Small TMZ dose
        obs, reward, terminated, truncated, info = env.step(action)
        total_reward += reward
        print(f"  Step {i+1}: vol={np.sum(env.u > 0.1):.0f}, reward={reward:.3f}, tox={sum(env.cum_toxicity.values()):.3f}")
        
        if terminated or truncated:
            break
    
    print(f"\n[Test] Total reward: {total_reward:.3f}")
    print(f"[Test] Final volume: {np.sum(env.u > 0.1):.0f}")
    print("\n[SUCCESS] Phase 12 Closed-Loop RL Environment operational!")