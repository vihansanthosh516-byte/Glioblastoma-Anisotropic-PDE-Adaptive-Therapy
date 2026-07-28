#!/usr/bin/env python3
"""
Phase 11: Virtual Biosensor Suite
==================================
Simulates continuous noisy measurements from implantable nanopore/microfluidic sensors
for the RL agent's observation space.

Sensors simulated:
- ctDNA Shedding Rate: ∝ tumor surface area × local proliferation (ρ)
- Interstitial Fluid Pressure (IFP): From solid tumor stress + fluid accumulation
- BBB Permeability: Scaled by local tumor density u(x)
- pH & Oxygenation (O₂): Micro-environmental hypoxia/acidosis indicators
"""
import numpy as np
from scipy.ndimage import gaussian_filter, distance_transform_edt
from scipy.signal import savgol_filter
from typing import Dict, Tuple, Optional, List
from dataclasses import dataclass, field


@dataclass
class SensorConfig:
    """Configuration for virtual biosensor noise and sampling characteristics."""
    # ctDNA sensor (nanopore sequencing)
    ctdna_base_rate: float = 0.01          # baseline shedding rate (fragments/day/mm²)
    ctdna_noise_cv: float = 0.15           # coefficient of variation (15% CV)
    ctdna_sampling_hz: float = 1.0/24      # once per day (clinical draw)
    
    # IFP sensor (microneedle pressure transducer)
    ifp_base_pressure: float = 5.0         # mmHg baseline
    ifp_noise_mmHg: float = 1.0            # mmHg std dev
    ifp_sampling_hz: float = 1.0/6         # every 6 hours
    
    # BBB permeability (contrast-enhanced MRI derived)
    bbb_base_perm: float = 0.15            # baseline permeability
    bbb_noise_cv: float = 0.1              # 10% CV
    bbb_sampling_hz: float = 1.0/72        # weekly MRI
    
    # pH sensor (optical fiber)
    ph_base: float = 7.4                   # physiological pH
    ph_noise: float = 0.02                 # pH units std dev
    ph_sampling_hz: float = 1.0/12         # twice daily
    
    # O2 sensor (optical quenching)
    o2_base: float = 40.0                  # mmHg baseline
    o2_noise_cv: float = 0.1               # 10% CV
    o2_sampling_hz: float = 1.0/6          # every 6 hours
    
    # Temporal smoothing (moving average window in samples)
    smoothing_window: int = 3


class VirtualBiosensorSuite:
    """
    Virtual biosensor suite simulating continuous in vivo measurements.
    
    Given a 3D tumor state u(x,y,z) and biophysical parameters,
    generates noisy continuous sensor readings with realistic noise models.
    """
    
    def __init__(self, config: SensorConfig = None, dx: float = 1.0):
        self.config = config or SensorConfig()
        self.dx = dx  # mm per voxel
        
        # Sensor state for temporal correlations
        self._ctdna_state = None
        self._ifp_state = None
        self._bbb_state = None
        self._ph_state = None
        self._o2_state = None
        
        # History buffers for temporal smoothing
        self._history = {
            'ctdna': [],
            'ifp': [],
            'bbb': [],
            'ph': [],
            'o2': []
        }
    
    def compute_tumor_surface_area(self, u: np.ndarray, threshold: float = 0.1) -> float:
        """
        Compute tumor surface area using marching cubes approximation.
        Surface area ≈ boundary voxels * dx²
        """
        mask = (u > threshold).astype(np.uint8)
        # Find boundary voxels (have at least one neighbor outside)
        from scipy.ndimage import binary_erosion
        eroded = binary_erosion(mask)
        boundary = mask - eroded
        surface_area = np.sum(boundary) * (self.dx ** 2)
        return surface_area
    
    def compute_tumor_volume(self, u: np.ndarray, threshold: float = 0.1) -> float:
        """Compute tumor volume in mm³."""
        return np.sum(u > threshold) * (self.dx ** 3)
    
    def compute_ctdna_shedding(self, u: np.ndarray, rho: np.ndarray, 
                                threshold: float = 0.1) -> float:
        """
        ctDNA shedding rate ∝ tumor surface area × local proliferation.
        
        ctDNA fragments shed from tumor cells at boundary into bloodstream.
        """
        surface_area = self.compute_tumor_surface_area(u, threshold)
        
        # Weight by proliferation at boundary
        mask = (u > threshold).astype(np.float32)
        from scipy.ndimage import binary_erosion
        boundary_mask = mask - binary_erosion(mask.astype(bool)).astype(np.float32)
        
        if np.sum(boundary_mask) > 0:
            avg_rho_boundary = np.sum(rho * boundary_mask) / np.sum(boundary_mask)
        else:
            avg_rho_boundary = np.mean(rho)
        
        # Shedding rate: fragments per day
        shedding_rate = self.config.ctdna_base_rate * surface_area * avg_rho_boundary * 1000  # scale
        
        return max(shedding_rate, 0.0)
    
    def compute_ifp(self, u: np.ndarray, young_modulus: float = 1000.0, 
                    poisson_ratio: float = 0.45) -> float:
        """
        Interstitial Fluid Pressure (IFP) from solid tumor stress.
        
        IFP ≈ solid stress = E * ε / (1 - ν) where ε is volumetric strain.
        Simplified: IFP ∝ tumor volume fraction * stiffness
        """
        volume_fraction = np.mean(u)
        
        # Solid stress approximation
        solid_stress = young_modulus * volume_fraction / (1 - poisson_ratio)
        
        # Convert to mmHg (1 Pa = 0.0075 mmHg)
        ifp_mmHg = self.config.ifp_base_pressure + solid_stress * 0.0075
        
        return ifp_mmHg
    
    def compute_bbb_permeability(self, u: np.ndarray, 
                                  threshold: float = 0.3) -> np.ndarray:
        """
        Blood-Brain Barrier (BBB) permeability scaled by local tumor density.
        
        BBB breakdown occurs where tumor density > threshold.
        Permeability = base + max_breakdown * sigmoid(u - threshold)
        """
        # Sigmoid response centered at threshold
        steepness = 10.0
        permeability = self.config.bbb_base_perm + \
            0.85 * (1.0 / (1.0 + np.exp(-steepness * (u - threshold))))
        
        return np.clip(permeability, 0.0, 1.0)
    
    def compute_ph(self, u: np.ndarray, rho: np.ndarray) -> float:
        """
        Micro-environmental pH from Warburg effect (glycolysis → lactate → acidosis).
        
        pH = 7.4 - k * metabolic_rate * density
        """
        metabolic_rate = np.mean(rho * u)  # proliferation-weighted density
        ph_drop = 0.3 * metabolic_rate * 10  # scale
        ph = self.config.ph_base - ph_drop
        
        return np.clip(ph, 6.5, 7.4)
    
    def compute_o2(self, u: np.ndarray, consumption_rate: float = 0.5) -> float:
        """
        Oxygen tension (pO₂) from diffusion-limited consumption.
        
        pO₂ = pO₂_boundary - consumption * tumor_radius / D
        """
        volume = self.compute_tumor_volume(u)
        radius = (3 * volume / (4 * np.pi)) ** (1/3) if volume > 0 else 0
        
        # Oxygen diffusion from boundary
        pO2_boundary = self.config.o2_base
        pO2_drop = consumption_rate * radius * 2  # mmHg
        
        pO2 = pO2_boundary - pO2_drop
        
        return max(pO2, 1.0)
    
    def _add_noise(self, value: float, cv: float = 0.0, 
                   additive_std: float = 0.0) -> float:
        """Add multiplicative and additive noise."""
        if cv > 0:
            value *= np.random.lognormal(mean=0, sigma=cv)
        if additive_std > 0:
            value += np.random.normal(0, additive_std)
        return value
    
    def _smooth(self, history: List[float], new_value: float, 
                window: int) -> float:
        """Temporal smoothing via moving average."""
        history.append(new_value)
        if len(history) > window:
            history.pop(0)
        return np.mean(history)
    
    def read_all(self, u: np.ndarray, rho: np.ndarray, 
                 t: float = 0.0) -> Dict[str, float]:
        """
        Read all virtual sensors at current time step.
        
        Returns dict with noisy sensor readings.
        """
        # Compute ground truth values
        ctdna_true = self.compute_ctdna_shedding(u, rho)
        ifp_true = self.compute_ifp(u)
        bbb_true = self.compute_bbb_permeability(u)
        ph_true = self.compute_ph(u, rho)
        o2_true = self.compute_o2(u)
        
        # Add noise and temporal correlation
        ctdna = self._add_noise(ctdna_true, cv=self.config.ctdna_noise_cv)
        ctdna = self._smooth(self._history['ctdna'], ctdna, self.config.smoothing_window)
        
        ifp = self._add_noise(ifp_true, additive_std=self.config.ifp_noise_mmHg)
        ifp = self._smooth(self._history['ifp'], ifp, self.config.smoothing_window)
        
        bbb_mean = np.mean(bbb_true)
        bbb = self._add_noise(bbb_mean, cv=self.config.bbb_noise_cv)
        bbb = self._smooth(self._history['bbb'], bbb, self.config.smoothing_window)
        
        ph = self._add_noise(ph_true, additive_std=self.config.ph_noise)
        ph = self._smooth(self._history['ph'], ph, self.config.smoothing_window)
        
        o2 = self._add_noise(o2_true, cv=self.config.o2_noise_cv)
        o2 = self._smooth(self._history['o2'], o2, self.config.smoothing_window)
        
        return {
            'ctdna_fragments_per_ml': ctdna,
            'ifp_mmHg': ifp,
            'bbb_permeability': bbb,
            'pH': ph,
            'pO2_mmHg': o2,
            # Ground truth for validation
            '_ctdna_true': ctdna_true,
            '_ifp_true': ifp_true,
            '_bbb_map': bbb_true,
            '_ph_true': ph_true,
            '_pO2_true': o2_true,
        }
    
    def get_sensor_names(self) -> List[str]:
        """Return list of sensor output keys (excluding ground truth)."""
        return ['ctdna_fragments_per_ml', 'ifp_mmHg', 'bbb_permeability', 
                'pH', 'pO2_mmHg']


class CircadianBiosensorSuite(VirtualBiosensorSuite):
    """
    Extended biosensor suite with circadian rhythm modulation.
    
    Circadian rhythms modulate:
    - ctDNA shedding (higher at night)
    - BBB permeability (lower at night)
    - Metabolic rate (pH, O₂)
    """
    
    def __init__(self, config: SensorConfig = None, dx: float = 1.0,
                 circadian_period: float = 24.0):
        super().__init__(config, dx)
        self.circadian_period = circadian_period  # hours
    
    def _circadian_modulation(self, t: float, amplitude: float, 
                               phase: float = 0.0, baseline: float = 1.0) -> float:
        """Circadian modulation factor [0, 1]."""
        return baseline + amplitude * np.sin(2 * np.pi * t / self.circadian_period + phase)
    
    def read_all(self, u: np.ndarray, rho: np.ndarray, t: float = 0.0) -> Dict[str, float]:
        """Read sensors with circadian modulation."""
        # Get base readings
        readings = super().read_all(u, rho, t)
        
        # Apply circadian modulation
        circ_ctdna = self._circadian_modulation(t, 0.3, phase=-np.pi/2)  # peak at night
        circ_bbb = self._circadian_modulation(t, 0.2, phase=np.pi/2)     # peak during day
        circ_metab = self._circadian_modulation(t, 0.15, phase=0)        # metabolic rhythm
        
        readings['ctdna_fragments_per_ml'] *= circ_ctdna
        readings['bbb_permeability'] *= circ_bbb
        readings['pH'] = readings['pH'] * circ_metab
        readings['pO2_mmHg'] = readings['pO2_mmHg'] * circ_metab
        
        return readings


# Factory function for easy instantiation
def create_biosensor_suite(config: SensorConfig = None, 
                           circadian: bool = False,
                           dx: float = 1.0) -> VirtualBiosensorSuite:
    """Factory for creating biosensor suite instances."""
    if circadian:
        return CircadianBiosensorSuite(config, dx)
    return VirtualBiosensorSuite(config, dx)


# Quick validation test
if __name__ == "__main__":
    print("=" * 60)
    print("Phase 11: Virtual Biosensor Suite - Validation Test")
    print("=" * 60)
    
    # Create synthetic tumor
    grid_size = 64
    u = np.zeros((grid_size, grid_size, grid_size))
    rho = np.full((grid_size, grid_size, grid_size), 0.02)
    
    # Spherical tumor
    cx, cy, cz = grid_size // 2, grid_size // 2, grid_size // 2
    z, y, x = np.ogrid[:grid_size, :grid_size, :grid_size]
    mask = (x - cx)**2 + (y - cy)**2 + (z - cz)**2 <= 10**2
    u[mask] = np.random.uniform(0.5, 1.0, size=mask.sum())
    rho[mask] = np.random.uniform(0.015, 0.03, size=mask.sum())
    
    # Create sensor suite
    config = SensorConfig()
    sensors = VirtualBiosensorSuite(config)
    
    # Read sensors
    print("\n[Test] Reading virtual biosensors...")
    readings = sensors.read_all(u, rho, t=0.0)
    
    for key, val in readings.items():
        if not key.startswith('_'):
            print(f"  {key}: {val:.4f}")
    
    # Test circadian version
    print("\n[Test] Circadian modulation at t=0, 6, 12, 18h...")
    circ_sensors = CircadianBiosensorSuite(config)
    for t in [0, 6, 12, 18]:
        r = circ_sensors.read_all(u, rho, t=t)
        print(f"  t={t}h: ctDNA={r['ctdna_fragments_per_ml']:.1f}, "
              f"BBB={r['bbb_permeability']:.3f}, "
              f"pH={r['pH']:.3f}, pO2={r['pO2_mmHg']:.1f}")
    
    print("\n[SUCCESS] Virtual Biosensor Suite operational!")