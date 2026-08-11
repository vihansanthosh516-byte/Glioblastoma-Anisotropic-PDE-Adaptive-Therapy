#!/usr/bin/env python3
"""
Master Cohort Dashboard - Aggregates All 183 Patient Digital Twins
==================================================================
Creates a single interactive Plotly HTML dashboard showing:
- All patients' tumor volume trajectories
- Personalized ρ/D parameter distributions
- Treatment outcomes by clinical event type
- Survival/progression correlations
- Interactive filtering by patient subgroups
"""

import os
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import glob


BASE_OUTPUT = Path("output/batch_digital_twins")


def load_all_patient_data():
    """Load data from all successful patient runs."""
    patients = []
    
    for patient_dir in BASE_OUTPUT.glob("PatientID_*"):
        if not patient_dir.is_dir():
            continue
        
        # Load inverse estimation
        inv_file = patient_dir / "inverse_estimation.json"
        rho, D, rho_ci, D_ci = None, None, None, None
        if inv_file.exists():
            try:
                with open(inv_file) as f:
                    data = json.load(f)
                rho = data.get('rho')
                D = data.get('D')
                rho_ci = data.get('rho_ci')
                D_ci = data.get('D_ci')
            except:
                pass
        
        # Load treatment schedule
        sched_file = patient_dir / "treatment_schedule.json"
        events = []
        if sched_file.exists():
            try:
                with open(sched_file) as f:
                    sched = json.load(f)
                events = sched.get('events', [])
            except:
                pass
        
        # Load simulation data
        sim_file = patient_dir / "clinical_simulation_data.npz"
        volume_history = None
        final_volume = None
        if sim_file.exists():
            try:
                data = np.load(sim_file, allow_pickle=True)
                volume_history = data.get('volume_history')
                if volume_history is not None and len(volume_history) > 0:
                    final_volume = float(volume_history[-1])
            except:
                pass
        
        # Determine treatment types
        treatment_types = set()
        has_surgery = False
        has_radiation = False
        has_chemo = False
        has_immuno = False
        for ev in events:
            t = ev.get('type', '')
            if t == 'surgery':
                has_surgery = True
            elif t == 'radiation':
                has_radiation = True
            elif t == 'chemotherapy':
                has_chemo = True
            elif t == 'immunotherapy':
                has_immuno = True
            treatment_types.add(t)
        
        # Load clinical timing
        clinical_file = Path(f"data/tcia/MU-Glioma-Post/{patient_dir.name}")
        clinical_data = {}
        if clinical_file.exists():
            try:
                import pandas as pd
                df = pd.read_excel("data/tcia/MU-Glioma-Post_ClinicalData-July2025.xlsx", sheet_name="MU Glioma Post")
                row = df[df['Patient_ID'] == patient_dir.name]
                if len(row) > 0:
                    r = row.iloc[0]
                    clinical_data = {
                        'time_to_progression': r.get('Time to First Progression (Days)'),
                        'overall_survival': r.get('Number of days from Diagnosis to death (Days)'),
                        'progression': r.get('Progression'),
                        'grade': r.get('Grade of Primary Brain Tumor'),
                        'diagnosis': r.get('Primary Diagnosis'),
                    }
            except:
                pass
        
        patients.append({
            'patient_id': patient_dir.name,
            'rho': rho,
            'D': D,
            'rho_ci': rho_ci,
            'D_ci': D_ci,
            'final_volume': final_volume,
            'volume_history': volume_history,
            'events': events,
            'has_surgery': has_surgery,
            'has_radiation': has_radiation,
            'has_chemo': has_chemo,
            'has_immuno': has_immuno,
            'treatment_types': list(treatment_types),
            'clinical': clinical_data,
        })
    
    return patients


def create_master_dashboard(patients):
    """Create comprehensive master cohort dashboard."""
    
    # Convert to DataFrame for easier plotting
    df = pd.DataFrame(patients)
    success_df = df[df['rho'].notna()].copy()
    
    # Create subplot figure with multiple panels
    fig = make_subplots(
        rows=3, cols=3,
        subplot_titles=[
            'Tumor Volume Trajectories (All Patients)',
            'Personalized Parameters: ρ vs D',
            'Final Tumor Volume Distribution',
            'Treatment Events Timeline',
            'Parameter Distributions by Treatment',
            'Clinical Outcomes vs Parameters',
            'Volume vs ρ (colored by D)',
            'Survival/Progression Analysis',
            'Cohort Summary Statistics'
        ],
        specs=[
            [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}],
            [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}],
            [{"secondary_y": False}, {"secondary_y": False}, {"secondary_y": False}],
        ],
        vertical_spacing=0.08,
        horizontal_spacing=0.08
    )
    
    # Color scheme
    colors = {
        'surgery_only': '#1f77b4',
        'surgery+rt': '#ff7f0e',
        'surgery+chemo': '#2ca02c',
        'surgery+rt+chemo': '#d62728',
        'with_immuno': '#9467bd',
        'other': '#8c564b'
    }
    
    # === PANEL 1: Volume Trajectories ===
    for _, patient in df.iterrows():
        if patient['volume_history'] is not None and len(patient['volume_history']) > 0:
            vol = patient['volume_history']
            days = np.arange(1, len(vol) + 1)
            
            # Determine color by treatment
            color = colors['other']
            if patient['has_immuno']:
                color = colors['with_immuno']
            elif patient['has_chemo'] and patient['has_radiation']:
                color = colors['surgery+rt+chemo']
            elif patient['has_radiation']:
                color = colors['surgery+rt']
            elif patient['has_chemo']:
                color = colors['surgery+chemo']
            else:
                color = colors['surgery_only']
            
            fig.add_trace(
                go.Scatter(
                    x=days, y=vol,
                    mode='lines',
                    line=dict(color=color, width=1),
                    opacity=0.3,
                    name=patient['patient_id'],
                    showlegend=False,
                    hovertext=patient['patient_id'],
                    hoverinfo='text+x+y'
                ),
                row=1, col=1
            )
    
    # Add median trajectory
    all_volumes = [p['volume_history'] for p in patients if p['volume_history'] is not None]
    if all_volumes:
        min_len = min(len(v) for v in all_volumes)
        median_vol = np.median([v[:min_len] for v in all_volumes], axis=0)
        fig.add_trace(
            go.Scatter(
                x=np.arange(1, min_len + 1),
                y=median_vol,
                mode='lines',
                line=dict(color='black', width=3),
                name='Median',
                showlegend=True,
                legendgroup='median'
            ),
            row=1, col=1
        )
    
    # === PANEL 2: ρ vs D Scatter ===
    personalized = success_df[success_df['rho'] != 0.005] if 'success_df' in locals() else df[df['rho'].notna() & (df['rho'] != 0.005)]
    default = df[(df['rho'] == 0.005) | df['rho'].isna()]
    
    if len(personalized) > 0:
        fig.add_trace(
            go.Scatter(
                x=personalized['rho'],
                y=personalized['D'],
                mode='markers',
                marker=dict(
                    size=10,
                    color='red',
                    symbol='circle',
                    line=dict(width=1, color='black')
                ),
                name='Personalized (n=' + str(len(personalized)) + ')',
                text=personalized['patient_id'],
                hovertemplate='%{text}<br>ρ=%{x:.6f}<br>D=%{y:.6f}<extra></extra>'
            ),
            row=1, col=2
        )
    
    if len(default) > 0:
        fig.add_trace(
            go.Scatter(
                x=default['rho'],
                y=default['D'],
                mode='markers',
                marker=dict(
                    size=8,
                    color='blue',
                    symbol='square',
                    opacity=0.5
                ),
                name='Population Default (n=' + str(len(default)) + ')',
                text=default['patient_id'],
                hovertemplate='%{text}<br>ρ=%{x:.6f}<br>D=%{y:.6f}<extra></extra>'
            ),
            row=1, col=2
        )
    
    # === PANEL 3: Final Volume Distribution ===
    volumes = df[df['final_volume'].notna()]['final_volume']
    if len(volumes) > 0:
        fig.add_trace(
            go.Histogram(
                x=volumes,
                nbinsx=30,
                marker_color='lightblue',
                marker_line=dict(color='black', width=1),
                name='Final Volume',
                showlegend=False
            ),
            row=1, col=3
        )
    
    # === PANEL 4: Treatment Timeline ===
    for i, (_, patient) in enumerate(df.iterrows()):
        y_pos = i
        for ev in patient['events']:
            start = ev.get('start_day', 0)
            end = ev.get('end_day', start)
            ev_type = ev.get('type', '')
            
            color_map = {
                'surgery': 'black',
                'radiation': 'orange',
                'chemotherapy': 'green',
                'immunotherapy': 'purple',
                'brachytherapy': 'brown',
                'other': 'gray'
            }
            
            fig.add_trace(
                go.Scatter(
                    x=[start, end],
                    y=[y_pos, y_pos],
                    mode='lines+markers',
                    line=dict(color=color_map.get(ev_type, 'gray'), width=8),
                    marker=dict(size=6),
                    name=ev_type,
                    showlegend=(i == 0),
                    legendgroup=ev_type,
                    hovertext=f"{patient['patient_id']}: {ev.get('description', ev_type)} ({start}-{end})",
                    hoverinfo='text'
                ),
                row=2, col=1
            )
    
    # === PANEL 5: Parameter Distributions by Treatment ===
    treatment_groups = []
    if len(df[df['has_radiation']]) > 0:
        treatment_groups.append(('Radiation', df[df['has_radiation']]['rho']))
    if len(df[df['has_chemo']]) > 0:
        treatment_groups.append(('Chemo', df[df['has_chemo']]['rho']))
    if len(df[df['has_immuno']]) > 0:
        treatment_groups.append(('Immuno', df[df['has_immuno']]['rho']))
    if len(df[~df['has_radiation'] & ~df['has_chemo'] & ~df['has_immuno']]) > 0:
        treatment_groups.append(('Surgery Only', df[~df['has_radiation'] & ~df['has_chemo'] & ~df['has_immuno']]['rho']))
    
    for name, vals in treatment_groups:
        vals_clean = vals.dropna()
        if len(vals_clean) > 0:
            fig.add_trace(
                go.Box(
                    y=vals_clean,
                    name=name,
                    boxpoints='outliers',
                    marker_color=colors.get(name.lower().replace(' ', '').replace('+', ''), 'gray'),
                    showlegend=False
                ),
                row=2, col=2
            )
    
    # === PANEL 6: Clinical Outcomes vs Parameters ===
    if 'clinical' in df.columns:
        clinical_df = df[df['clinical'].apply(lambda x: isinstance(x, dict) and x.get('time_to_progression') is not None)].copy()
        if len(clinical_df) > 0:
            clinical_df['ttp'] = clinical_df['clinical'].apply(lambda x: x.get('time_to_progression'))
            clinical_df = clinical_df[clinical_df['ttp'].notna()]
            
            if len(clinical_df) > 0:
                fig.add_trace(
                    go.Scatter(
                        x=clinical_df['rho'],
                        y=clinical_df['ttp'],
                        mode='markers',
                        marker=dict(
                            size=10,
                            color=clinical_df['D'],
                            colorscale='Viridis',
                            showscale=True,
                            colorbar=dict(title="D", thickness=10, len=0.3, x=0.33, y=0.35)
                        ),
                        text=clinical_df['patient_id'],
                        hovertemplate='%{text}<br>ρ=%{x:.6f}<br>TTP=%{y:.0f} days<extra></extra>',
                        name='TTP vs ρ',
                        showlegend=False
                    ),
                    row=2, col=3
                )
    
    # === PANEL 7: Volume vs ρ colored by D ===
    vol_rho_df = df[df['final_volume'].notna() & df['rho'].notna() & df['D'].notna()]
    if len(vol_rho_df) > 0:
        fig.add_trace(
            go.Scatter(
                x=vol_rho_df['rho'],
                y=vol_rho_df['final_volume'],
                mode='markers',
                marker=dict(
                    size=10,
                    color=vol_rho_df['D'],
                    colorscale='Plasma',
                    showscale=True,
                    colorbar=dict(title="D (mm²/day)", thickness=10, len=0.3, x=0.66, y=0.35)
                ),
                text=vol_rho_df['patient_id'],
                hovertemplate='%{text}<br>ρ=%{x:.6f}<br>Vol=%{y:.0f}<br>D=%{marker.color:.6f}<extra></extra>',
                name='Vol vs ρ',
                showlegend=False
            ),
            row=3, col=1
        )
    
    # === PANEL 8: Survival Analysis ===
    if 'clinical' in df.columns:
        surv_df = df[df['clinical'].apply(lambda x: isinstance(x, dict) and x.get('overall_survival') is not None)].copy()
        if len(surv_df) > 0:
            surv_df['os'] = surv_df['clinical'].apply(lambda x: x.get('overall_survival'))
            surv_df = surv_df[surv_df['os'].notna()]
            
            if len(surv_df) > 0:
                # KM-style: sort by survival time
                surv_df = surv_df.sort_values('os')
                surv_df['survival_prob'] = 1 - np.arange(len(surv_df)) / len(surv_df)
                
                fig.add_trace(
                    go.Scatter(
                        x=surv_df['os'],
                        y=surv_df['survival_prob'] * 100,
                        mode='lines+markers',
                        line=dict(color='red', width=2),
                        name='Overall Survival',
                        text=surv_df['patient_id'],
                        hovertemplate='%{text}<br>OS=%{x:.0f} days<br>S=%{y:.1f}%<extra></extra>'
                    ),
                    row=3, col=2
                )
    
    # === PANEL 9: Summary Statistics Table ===
    stats_text = f"""
    <b>COHORT SUMMARY (N={len(df)} patients)</b><br><br>
    <b>Parameter Estimation:</b><br>
    • Personalized: {len(personalized)} ({len(personalized)/len(df)*100:.1f}%)<br>
    • Population default: {len(default)} ({len(default)/len(df)*100:.1f}%)<br><br>
    <b>ρ Range:</b> {df['rho'].min():.6f} – {df['rho'].max():.6f} /day<br>
    <b>D Range:</b> {df['D'].min():.6f} – {df['D'].max():.6f} mm²/day<br><br>
    <b>Treatments:</b><br>
    • Surgery: {df['has_surgery'].sum()}<br>
    • Radiation: {df['has_radiation'].sum()}<br>
    • Chemotherapy: {df['has_chemo'].sum()}<br>
    • Immunotherapy: {df['has_immuno'].sum()}<br><br>
    <b>Final Volume:</b><br>
    • Mean: {df['final_volume'].mean():.1f} voxels<br>
    • Median: {df['final_volume'].median():.1f} voxels<br>
    • Min: {df['final_volume'].min():.1f}<br>
    • Max: {df['final_volume'].max():.1f}
    """
    
    fig.add_annotation(
        text=stats_text,
        xref="x9", yref="y9",
        x=0.5, y=0.5,
        showarrow=False,
        align="left",
        font=dict(size=11, family="monospace"),
        bgcolor="rgba(240,240,240,0.9)",
        bordercolor="black",
        borderwidth=1,
        row=3, col=3
    )
    
    # Update layout
    fig.update_layout(
        title={
            'text': 'GBM Digital Twin Cohort Master Dashboard (183 Patients)',
            'x': 0.5,
            'font': {'size': 20}
        },
        height=1600,
        width=1800,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.01,
            font=dict(size=9)
        ),
        template="plotly_white"
    )
    
    # Update axis labels
    fig.update_xaxes(title_text="Day", row=1, col=1)
    fig.update_yaxes(title_text="Tumor Volume (voxels)", row=1, col=1)
    
    fig.update_xaxes(title_text="ρ (day⁻¹)", row=1, col=2)
    fig.update_yaxes(title_text="D (mm²/day)", row=1, col=2)
    
    fig.update_xaxes(title_text="Final Volume (voxels)", row=1, col=3)
    fig.update_yaxes(title_text="Count", row=1, col=3)
    
    fig.update_xaxes(title_text="Day from Diagnosis", row=2, col=1)
    fig.update_yaxes(title_text="Patient Index", row=2, col=1)
    
    fig.update_xaxes(title_text="Treatment Group", row=2, col=2)
    fig.update_yaxes(title_text="ρ (day⁻¹)", row=2, col=2)
    
    fig.update_xaxes(title_text="ρ (day⁻¹)", row=2, col=3)
    fig.update_yaxes(title_text="Time to Progression (days)", row=2, col=3)
    
    fig.update_xaxes(title_text="ρ (day⁻¹)", row=3, col=1)
    fig.update_yaxes(title_text="Final Volume (voxels)", row=3, col=1)
    
    fig.update_xaxes(title_text="Overall Survival (days)", row=3, col=2)
    fig.update_yaxes(title_text="Survival Probability (%)", row=3, col=2)
    
    # Hide axes for stats panel
    fig.update_xaxes(visible=False, row=3, col=3)
    fig.update_yaxes(visible=False, row=3, col=3)
    
    return fig


def main():
    print("=" * 60)
    print("  MASTER COHORT DASHBOARD GENERATOR")
    print("=" * 60)
    
    print("\n[LOAD] Loading all patient data...")
    patients = load_all_patient_data()
    print(f"[LOAD] Loaded {len(patients)} patients")
    
    print("\n[BUILD] Creating master dashboard...")
    fig = create_master_dashboard(patients)
    
    output_path = BASE_OUTPUT / "master_cohort_dashboard.html"
    fig.write_html(str(output_path), include_plotlyjs='cdn')
    print(f"[SUCCESS] Master dashboard saved to: {output_path}")
    
    # Also save as static image
    try:
        img_path = BASE_OUTPUT / "master_cohort_dashboard.png"
        fig.write_image(str(img_path), width=1800, height=1600, scale=2)
        print(f"[SUCCESS] Static image saved to: {img_path}")
    except Exception as e:
        print(f"[WARNING] Could not save static image: {e}")
    
    # Print summary
    print("\n" + "=" * 60)
    print("COHORT SUMMARY")
    print("=" * 60)
    df = pd.DataFrame(patients)
    personalized = df[df['rho'].notna() & (df['rho'] != 0.005)]
    default = df[(df['rho'] == 0.005) | df['rho'].isna()]
    
    print(f"Total patients: {len(df)}")
    print(f"Personalized: {len(personalized)} ({len(personalized)/len(df)*100:.1f}%)")
    print(f"Default params: {len(default)} ({len(default)/len(df)*100:.1f}%)")
    print(f"ρ range: {df['rho'].min():.6f} – {df['rho'].max():.6f} /day")
    print(f"D range: {df['D'].min():.6f} – {df['D'].max():.6f} mm²/day")
    print(f"Surgery: {df['has_surgery'].sum()}")
    print(f"Radiation: {df['has_radiation'].sum()}")
    print(f"Chemo: {df['has_chemo'].sum()}")
    print(f"Immuno: {df['has_immuno'].sum()}")
    print(f"Final volume - Mean: {df['final_volume'].mean():.1f}, Median: {df['final_volume'].median():.1f}")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    import pandas as pd
    import sys
    sys.exit(main())