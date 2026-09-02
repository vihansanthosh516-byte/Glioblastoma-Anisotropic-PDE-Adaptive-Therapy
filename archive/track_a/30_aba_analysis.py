#!/usr/bin/env python3
"""
Month 3, Week 4: ABA Analysis & Clinical Correlation
Analyzes invasion simulation metrics and correlates with empirical data.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def load_simulation_data() -> Tuple[np.ndarray, dict]:
    """Load invasion simulation metrics and summary."""
    metrics = np.load("output/invasion_metrics.npy", allow_pickle=True)
    if metrics.dtype == object:
        # Convert structured data
        steps = metrics[:, 0].astype(int)
        healthy = metrics[:, 1].astype(int)
        periphery = metrics[:, 2].astype(int)
        core = metrics[:, 3].astype(int)
        necrotic = metrics[:, 4].astype(int)
        front_r = metrics[:, 5].astype(float)
        prolif = metrics[:, 6].astype(int)
        trans = metrics[:, 7].astype(int)
        necro = metrics[:, 8].astype(int)
    else:
        steps = metrics[:, 0]
        healthy = metrics[:, 1]
        periphery = metrics[:, 2]
        core = metrics[:, 3]
        necrotic = metrics[:, 4]
        front_r = metrics[:, 5]
        prolif = metrics[:, 6]
        trans = metrics[:, 7]
        necro = metrics[:, 8]
    
    with open("output/invasion_summary.json") as f:
        summary = json.load(f)
    
    return {
        'steps': steps,
        'healthy': healthy,
        'periphery': periphery,
        'core': core,
        'necrotic': necrotic,
        'front_radius': front_r,
        'proliferation': prolif,
        'transition': trans,
        'necrosis': necro,
    }, summary


def compute_invasion_kinetics(data: dict) -> Dict:
    """Compute invasion kinetic parameters from simulation."""
    steps = data['steps']
    front_r = data['front_radius']
    core = data['core']
    periphery = data['periphery']
    necrotic = data['necrotic']
    
    # Wave speed: linear fit to front radius
    if len(front_r) > 10:
        coeffs = np.polyfit(steps, front_r, 1)
        wave_speed = coeffs[0]  # pixels/step
    else:
        wave_speed = 0.0
    
    # Core growth rate
    if len(core) > 10:
        core_coeffs = np.polyfit(steps, core, 1)
        core_growth = core_coeffs[0]
    else:
        core_growth = 0.0
    
    # Periphery dynamics
    max_periphery = max(periphery)
    final_periphery = periphery[-1]
    
    # Necrosis accumulation rate
    if len(necrotic) > 10:
        nec_coeffs = np.polyfit(steps, necrotic, 1)
        necrosis_rate = nec_coeffs[0]
    else:
        necrosis_rate = 0.0
    
    # Transition zone width (periphery + core interface)
    # Approximate as where both periphery and core are present
    transition_metric = np.array(periphery) * np.array(core)
    
    kinetics = {
        'wave_speed_pixels_per_step': float(wave_speed),
        'wave_speed_um_per_hour': float(wave_speed * 10),  # Assuming 10um/pixel, 1 step = 1 hour
        'core_growth_rate': float(core_growth),
        'max_periphery_cells': int(max_periphery),
        'final_periphery_cells': int(final_periphery),
        'necrosis_accumulation_rate': float(necrosis_rate),
        'final_necrotic_fraction': float(necrotic[-1] / sum([data[k][-1] for k in ['healthy', 'periphery', 'core', 'necrotic']])),
        'total_steps': int(steps[-1]),
    }
    
    return kinetics


def clinical_correlation_analysis(data: dict, kinetics: dict) -> Dict:
    """
    Correlate simulation metrics with clinical expectations.
    Based on glioblastoma invasion literature.
    """
    print("[CLINICAL] Performing clinical correlation analysis...")
    
    # Literature values for GBM
    # Wave speed: ~0.5-2 mm/day = ~5-20 um/hour (our scale)
    # Our wave speed in um/hour
    wave_speed_um_hr = kinetics['wave_speed_um_per_hour']
    
    # Core growth doubling time
    core_final = data['core'][-1]
    core_initial = data['core'][0]
    if core_initial > 0 and len(data['core']) > 0:
        total_time = data['steps'][-1]  # hours
        if total_time > 0:
            doubling_time_hr = total_time * np.log(2) / np.log(core_final / core_initial) if core_final > core_initial else float('inf')
        else:
            doubling_time_hr = float('inf')
    else:
        doubling_time_hr = float('inf')
    
    # Periphery as invasive front marker
    periphery_fraction = data['periphery'][-1] / (data['periphery'][-1] + data['core'][-1])
    
    # Necrotic fraction (GBM hallmark)
    total_tumor = data['core'][-1] + data['periphery'][-1] + data['necrotic'][-1]
    necrotic_fraction = data['necrotic'][-1] / total_tumor if total_tumor > 0 else 0
    
    correlation = {
        'wave_speed_um_per_hour': wave_speed_um_hr,
        'literature_range_um_per_hour': [5, 20],
        'wave_speed_in_range': 5 <= wave_speed_um_hr <= 20,
        'core_doubling_time_hours': float(doubling_time_hr),
        'literature_doubling_time_days': [7, 30],  # GBM doubling time
        'doubling_time_in_range': 7*24 <= doubling_time_hr <= 30*24 if doubling_time_hr != float('inf') else False,
        'invasive_front_periphery_fraction': float(periphery_fraction),
        'necrotic_fraction': float(necrotic_fraction),
        'literature_necrotic_fraction_range': [0.1, 0.4],
        'necrotic_fraction_in_range': 0.1 <= necrotic_fraction <= 0.4,
        'histological_pattern': 'infiltrative' if wave_speed_um_hr > 8 else 'circumscribed',
    }
    
    return correlation


def plot_invasion_dynamics(data: dict, kinetics: dict, correlation: dict) -> None:
    """Generate publication figures for invasion dynamics."""
    steps = data['steps']
    
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # 1. Cell population dynamics
    ax1 = fig.add_subplot(gs[0, :2])
    ax1.plot(steps, data['healthy'], label='Healthy', color='#2E8B57', linewidth=2)
    ax1.plot(steps, data['periphery'], label='Periphery', color='#FF8C00', linewidth=2)
    ax1.plot(steps, data['core'], label='Core', color='#DC143C', linewidth=2)
    ax1.plot(steps, data['necrotic'], label='Necrotic', color='#696969', linewidth=2)
    ax1.set_xlabel('Simulation Step (hours)')
    ax1.set_ylabel('Cell Count')
    ax1.set_title('Tumor Cell Population Dynamics', fontweight='bold', fontsize=12)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # 2. Invasion front radius
    ax2 = fig.add_subplot(gs[0, 2])
    ax2.plot(steps, data['front_radius'], color='#4B0082', linewidth=2)
    # Linear fit
    coeffs = np.polyfit(steps, data['front_radius'], 1)
    fit = np.polyval(coeffs, steps)
    ax2.plot(steps, fit, '--', color='red', alpha=0.7, label=f'Fit: {coeffs[0]:.3f} px/step')
    ax2.set_xlabel('Step')
    ax2.set_ylabel('Front Radius (pixels)')
    ax2.set_title('Invasion Front Propagation', fontweight='bold', fontsize=12)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    # 3. Periphery/Core ratio (invasive potential)
    ax3 = fig.add_subplot(gs[1, 0])
    ratio = np.array(data['periphery']) / (np.array(data['periphery']) + np.array(data['core']) + 1)
    ax3.plot(steps, ratio, color='#FF8C00', linewidth=2)
    ax3.set_xlabel('Step')
    ax3.set_ylabel('Periphery / (Periphery + Core)')
    ax3.set_title('Invasive Front Fraction', fontweight='bold', fontsize=12)
    ax3.grid(True, alpha=0.3)
    
    # 4. Event rates
    ax4 = fig.add_subplot(gs[1, 1])
    ax4.plot(steps, data['proliferation'], label='Proliferation', color='#2E8B57', alpha=0.7)
    ax4.plot(steps, data['transition'], label='Transition (P→C)', color='#FF8C00', alpha=0.7)
    ax4.plot(steps, data['necrosis'], label='Necrosis', color='#696969', alpha=0.7)
    ax4.set_xlabel('Step')
    ax4.set_ylabel('Event Count / Step')
    ax4.set_title('Cellular Event Rates', fontweight='bold', fontsize=12)
    ax4.legend(fontsize=9)
    ax4.grid(True, alpha=0.3)
    
    # 5. Phase portrait: Core vs Periphery
    ax5 = fig.add_subplot(gs[1, 2])
    ax5.plot(data['core'], data['periphery'], color='#8B0000', linewidth=1, alpha=0.7)
    ax5.scatter(data['core'][0], data['periphery'][0], color='green', s=100, label='Start', zorder=5)
    ax5.scatter(data['core'][-1], data['periphery'][-1], color='red', s=100, label='End', zorder=5)
    ax5.set_xlabel('Core Cells')
    ax5.set_ylabel('Periphery Cells')
    ax5.set_title('Phase Portrait: Core vs Periphery', fontweight='bold', fontsize=12)
    ax5.legend(fontsize=9)
    ax5.grid(True, alpha=0.3)
    
    # 6. Necrotic accumulation
    ax6 = fig.add_subplot(gs[2, 0])
    ax6.plot(steps, data['necrotic'], color='#696969', linewidth=2)
    ax6.fill_between(steps, 0, data['necrotic'], alpha=0.3, color='gray')
    ax6.set_xlabel('Step')
    ax6.set_ylabel('Necrotic Cells')
    ax6.set_title('Necrosis Accumulation', fontweight='bold', fontsize=12)
    ax6.grid(True, alpha=0.3)
    
    # 7. Velocity profile (derivative of front)
    ax7 = fig.add_subplot(gs[2, 1])
    if len(steps) > 1:
        velocity = np.diff(data['front_radius']) / np.diff(steps)
        v_steps = steps[:-1] + 0.5
        ax7.plot(v_steps, velocity, color='#4B0082', alpha=0.7, linewidth=1)
        ax7.axhline(np.mean(velocity), color='red', linestyle='--', label=f'Mean: {np.mean(velocity):.3f}')
        ax7.set_xlabel('Step')
        ax7.set_ylabel('Front Velocity (px/step)')
        ax7.set_title('Instantaneous Invasion Velocity', fontweight='bold', fontsize=12)
        ax7.legend(fontsize=9)
        ax7.grid(True, alpha=0.3)
    
    # 8. Clinical correlation summary
    ax8 = fig.add_subplot(gs[2, 2])
    ax8.axis('off')
    
    summary_text = "CLINICAL CORRELATION SUMMARY\n" + "="*30 + "\n\n"
    summary_text += f"Wave Speed: {correlation['wave_speed_um_per_hour']:.1f} µm/hr\n"
    summary_text += f"Literature: {correlation['literature_range_um_per_hour'][0]}-{correlation['literature_range_um_per_hour'][1]} µm/hr\n"
    summary_text += f"✓ In Range: {correlation['wave_speed_in_range']}\n\n"
    
    if correlation['core_doubling_time_hours'] != float('inf'):
        summary_text += f"Core Doubling: {correlation['core_doubling_time_hours']/24:.1f} days\n"
        summary_text += f"Literature: {correlation['literature_doubling_time_days'][0]}-{correlation['literature_doubling_time_days'][1]} days\n"
        summary_text += f"✓ In Range: {correlation['doubling_time_in_range']}\n\n"
    else:
        summary_text += "Core Doubling: N/A\n\n"
    
    summary_text += f"Necrotic Fraction: {correlation['necrotic_fraction']:.2f}\n"
    summary_text += f"Literature: {correlation['literature_necrotic_fraction_range'][0]}-{correlation['literature_necrotic_fraction_range'][1]}\n"
    summary_text += f"✓ In Range: {correlation['necrotic_fraction_in_range']}\n\n"
    
    summary_text += f"Pattern: {correlation['histological_pattern']}\n"
    summary_text += f"Periphery Fraction: {correlation['invasive_front_periphery_fraction']:.2f}"
    
    ax8.text(0.05, 0.95, summary_text, transform=ax8.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.suptitle('Month 3: Stochastic Agent-Based Tumor Invasion Analysis', 
                 fontsize=16, fontweight='bold', y=0.98)
    plt.savefig('output/invasion_dynamics_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("[PLOT] Saved invasion_dynamics_analysis.png")


def export_analysis_results(kinetics: dict, correlation: dict) -> None:
    """Export all analysis results."""
    # Convert numpy types to Python types
    def convert(obj):
        if isinstance(obj, (np.integer, np.floating)):
            return obj.item()
        elif isinstance(obj, np.bool_):
            return bool(obj)
        elif isinstance(obj, dict):
            return {k: convert(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert(v) for v in obj]
        return obj
    
    results = {
        'invasion_kinetics': convert(kinetics),
        'clinical_correlation': convert(correlation),
        'timestamp': time.strftime('%Y-%m-%d %H:%M:%S'),
    }
    
    with open("output/aba_analysis_results.json", "w") as f:
        json.dump(results, f, indent=2)
    
    # TSV summary
    with open("output/aba_kinetics_summary.tsv", "w") as f:
        f.write("parameter\tvalue\tunits\n")
        for k, v in kinetics.items():
            f.write(f"{k}\t{v}\t\n")
        f.write("\n# Clinical Correlation\n")
        for k, v in correlation.items():
            f.write(f"{k}\t{v}\t\n")
    
    print("[EXPORT] Saved aba_analysis_results.json and aba_kinetics_summary.tsv")


def main():
    print("=" * 60)
    print("MONTH 3 WEEK 4: ABA ANALYSIS & CLINICAL CORRELATION")
    print("=" * 60)
    
    # Load data
    data, summary = load_simulation_data()
    print(f"[LOAD] Loaded {len(data['steps'])} time points")
    
    # Compute kinetics
    kinetics = compute_invasion_kinetics(data)
    print(f"[KINETICS] Wave speed: {kinetics['wave_speed_um_per_hour']:.2f} µm/hr")
    print(f"[KINETICS] Core growth rate: {kinetics['core_growth_rate']:.2f} cells/step")
    print(f"[KINETICS] Necrosis rate: {kinetics['necrosis_accumulation_rate']:.2f} cells/step")
    
    # Clinical correlation
    correlation = clinical_correlation_analysis(data, kinetics)
    print(f"[CLINICAL] Wave speed in range: {correlation['wave_speed_in_range']}")
    print(f"[CLINICAL] Necrotic fraction in range: {correlation['necrotic_fraction_in_range']}")
    print(f"[CLINICAL] Pattern: {correlation['histological_pattern']}")
    
    # Visualization
    plot_invasion_dynamics(data, kinetics, correlation)
    
    # Export
    export_analysis_results(kinetics, correlation)
    
    print("\n[SUCCESS] Month 3 Week 4 Complete: ABA Analysis & Clinical Correlation")
    print("  - output/aba_analysis_results.json")
    print("  - output/aba_kinetics_summary.tsv")
    print("  - output/invasion_dynamics_analysis.png")


if __name__ == "__main__":
    main()