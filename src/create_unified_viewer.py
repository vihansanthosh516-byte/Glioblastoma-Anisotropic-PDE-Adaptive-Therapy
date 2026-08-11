#!/usr/bin/env python3
"""
Unified Patient 4D Viewer - Switch Between All 183 Patient 4D Simulations
=========================================================================
Creates a SINGLE HTML file with:
- Dropdown/search to select any patient
- Each patient's 4D tumor evolution (time slider + play/pause)
- Their personalized treatment schedule overlay
- Their clinical parameters (rho, D, treatment events)
"""

import os
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from pathlib import Path
import glob


BASE_OUTPUT = Path("output/batch_digital_twins")


def load_all_patient_4d_data():
    """Load 4D simulation data for all patients."""
    patients = {}
    
    for patient_dir in BASE_OUTPUT.glob("PatientID_*"):
        if not patient_dir.is_dir():
            continue
        
        pid = patient_dir.name
        
        sim_file = patient_dir / "clinical_simulation_data.npz"
        schedule_file = patient_dir / "treatment_schedule.json"
        inverse_file = patient_dir / "inverse_estimation.json"
        
        if not sim_file.exists() or not schedule_file.exists():
            continue
        
        try:
            # Load simulation data
            sim_data = np.load(sim_file, allow_pickle=True)
            volume_history = sim_data.get('volume_history')
            drug_conc = sim_data.get('drug_concentration')
            rt_dose = sim_data.get('radiation_dose')
            final_state = sim_data.get('final_state')
            
            # Load treatment schedule
            with open(schedule_file) as f:
                schedule = json.load(f)
            events = schedule.get('events', [])
            
            # Load inverse estimation
            rho, D = 0.015, 0.01
            if inverse_file.exists():
                with open(inverse_file) as f:
                    inv = json.load(f)
                rho = float(inv.get('rho', 0.015))
                D = float(inv.get('D', 0.01))
            
            # Convert numpy arrays to Python lists with native types
            vol_hist = volume_history.tolist() if volume_history is not None else []
            drug_list = drug_conc.tolist() if drug_conc is not None else []
            rt_list = rt_dose.tolist() if rt_dose is not None else []
            
            # Ensure native Python types
            vol_hist = [float(v) for v in vol_hist]
            drug_list = [float(v) for v in drug_list]
            rt_list = [float(v) for v in rt_list]
            
            init_vol = float(vol_hist[0]) if vol_hist else 0.0
            final_vol = float(vol_hist[-1]) if vol_hist else 0.0
            
            patients[pid] = {
                'patient_id': pid,
                'volume_history': vol_hist,
                'drug_concentration': drug_list,
                'radiation_dose': rt_list,
                'events': events,
                'rho': float(rho),
                'D': float(D),
                'initial_volume': init_vol,
                'final_volume': final_vol,
            }
            
        except Exception as e:
            print(f"Error loading {pid}: {e}")
            continue
    
    return patients


def create_unified_viewer(patients):
    """Create single HTML with dropdown to switch between patient 4D views."""
    
    patient_ids = sorted(patients.keys())
    n_patients = len(patient_ids)
    
    # Get max simulation days
    max_days = max(len(p['volume_history']) for p in patients.values() if p['volume_history'])
    
    # Prepare patient data for JavaScript
    patient_data = {}
    for pid in patient_ids:
        p = patients[pid]
        patient_data[pid] = {
            'patient_id': pid,
            'volume_history': p['volume_history'],
            'drug_concentration': p['drug_concentration'],
            'radiation_dose': p['radiation_dose'],
            'events': p['events'],
            'rho': p['rho'],
            'D': p['D'],
            'initial_volume': p['initial_volume'],
            'final_volume': p['final_volume'],
        }
    
def create_unified_viewer(patients):
    """Create single HTML with dropdown to switch between patient 4D views."""
    
    patient_ids = sorted(patients.keys())
    n_patients = len(patient_ids)
    
    # Prepare patient data for JavaScript
    patient_data = {}
    for pid in patient_ids:
        p = patients[pid]
        patient_data[pid] = {
            'patient_id': pid,
            'volume_history': p['volume_history'],
            'drug_concentration': p['drug_concentration'],
            'radiation_dose': p['radiation_dose'],
            'events': p['events'],
            'rho': p['rho'],
            'D': p['D'],
            'initial_volume': p['initial_volume'],
            'final_volume': p['final_volume'],
        }
    
    # Use double curly braces for JavaScript template literals in f-string
    js_template = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>GBM Digital Twin - Unified 4D Patient Viewer</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background: #f5f5f5;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            padding: 20px;
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            margin-bottom: 10px;
        }}
        .subtitle {{
            text-align: center;
            color: #7f8c8d;
            margin-bottom: 30px;
        }}
        .patient-selector {{
            background: #f8f9fa;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
            border: 1px solid #e9ecef;
        }}
        .selector-row {{
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            align-items: center;
        }}
        .search-box {{
            flex: 1;
            min-width: 250px;
        }}
        .search-box input {{
            width: 100%;
            padding: 10px 15px;
            border: 2px solid #ddd;
            border-radius: 6px;
            font-size: 16px;
        }}
        .search-box input:focus {{
            border-color: #3498db;
            outline: none;
        }}
        .dropdown-wrapper {{
            min-width: 200px;
        }}
        .dropdown-wrapper select {{
            width: 100%;
            padding: 10px 15px;
            border: 2px solid #ddd;
            border-radius: 6px;
            font-size: 16px;
            background: white;
        }}
        .patient-info {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 15px;
            margin-top: 15px;
            padding: 15px;
            background: #e8f4fd;
            border-radius: 6px;
        }}
        .info-box {{
            text-align: center;
        }}
        .info-label {{
            font-size: 12px;
            color: #7f8c8d;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .info-value {{
            font-size: 18px;
            font-weight: 600;
            color: #2c3e50;
        }}
        .charts-row {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-top: 20px;
        }}
        .chart-container {{
            background: #fafafa;
            border-radius: 8px;
            padding: 15px;
            border: 1px solid #eee;
        }}
        .chart-title {{
            font-size: 14px;
            font-weight: 600;
            color: #2c3e50;
            margin-bottom: 10px;
            padding-bottom: 8px;
            border-bottom: 2px solid #3498db;
        }}
        .timeline-container {{
            grid-column: 1 / -1;
            margin-top: 20px;
        }}
        .treatment-timeline {{
            height: 200px;
        }}
        @media (max-width: 1000px) {{
            .charts-row {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>GBM Digital Twin - Unified 4D Patient Viewer</h1>
        <p class="subtitle">183 Patients | Interactive 4D Tumor Evolution | Clinical Treatment Timelines</p>
        
        <div class="patient-selector">
            <div class="selector-row">
                <div class="dropdown-wrapper">
                    <select id="patientDropdown">
                        <option value="">-- Select Patient --</option>
                        {patient_options}
                    </select>
                </div>
            </div>
            
            <div id="patientInfo" class="patient-info" style="display:none;">
                <div class="info-box">
                    <div class="info-label">Patient ID</div>
                    <div class="info-value" id="infoPatientId">-</div>
                </div>
                <div class="info-box">
                    <div class="info-label">rho (Proliferation)</div>
                    <div class="info-value" id="infoRho">-</div>
                </div>
                <div class="info-box">
                    <div class="info-label">D (Diffusivity)</div>
                    <div class="info-value" id="infoD">-</div>
                </div>
                <div class="info-box">
                    <div class="info-label">Initial Volume</div>
                    <div class="info-value" id="infoVol0">-</div>
                </div>
                <div class="info-box">
                    <div class="info-label">Final Volume</div>
                    <div class="info-value" id="infoVolF">-</div>
                </div>
                <div class="info-box">
                    <div class="info-label">Volume Change</div>
                    <div class="info-value" id="infoVolChange">-</div>
                </div>
            </div>
        </div>
        
        <div class="charts-row">
            <div class="chart-container">
                <div class="chart-title">[Chart] Tumor Volume Evolution</div>
                <div id="volumePlot" style="height: 400px;"></div>
            </div>
            <div class="chart-container">
                <div class="chart-title">[Schedule] Treatment Schedule</div>
                <div id="treatmentPlot" style="height: 400px;"></div>
            </div>
        </div>
        
        <div class="timeline-container chart-container">
            <div class="chart-title">[Timeline] Clinical Treatment Timeline</div>
            <div id="timelinePlot" class="treatment-timeline"></div>
        </div>
    </div>

    <script>
        // Embedded patient data
        const patientData = {patient_data_json};
        const patientIds = {patient_ids_json};
        
        let currentPatient = null;
        
        // DOM elements
        const dropdown = document.getElementById('patientDropdown');
        const infoPanel = document.getElementById('patientInfo');
        
        // Dropdown change
        dropdown.addEventListener('change', function() {{
            if (this.value) {{
                loadPatient(this.value);
            }} else {{
                hideInfo();
            }}
        }});
        
        // Load patient data and render plots
        function loadPatient(pid) {{
            console.log('Loading patient:', pid);
            currentPatient = pid;
            const data = patientData[pid];
            if (!data) {
                console.error('No data for patient:', pid);
                return;
            }
            
            // Update info panel
            document.getElementById('infoPatientId').textContent = pid;
            document.getElementById('infoRho').textContent = data.rho.toFixed(6) + ' /day';
            document.getElementById('infoD').textContent = data.D.toFixed(6) + ' mm2/day';
            document.getElementById('infoVol0').textContent = data.initial_volume.toLocaleString() + ' voxels';
            document.getElementById('infoVolF').textContent = data.final_volume.toLocaleString() + ' voxels';
            
            const change = data.final_volume - data.initial_volume;
            const pct = data.initial_volume > 0 ? ((change / data.initial_volume) * 100).toFixed(1) : 0;
            const changeEl = document.getElementById('infoVolChange');
            changeEl.textContent = (change >= 0 ? '+' : '') + change.toLocaleString() + ' (' + pct + '%)';
            changeEl.style.color = change <= 0 ? '#27ae60' : '#e74c3c';
            
            infoPanel.style.display = 'grid';
            
            // Render plots
            console.log('Rendering plots for:', pid);
            renderVolumePlot(data);
            renderTreatmentPlot(data);
            renderTimelinePlot(data);
        }}
        
        function hideInfo() {{
            infoPanel.style.display = 'none';
            // Clear plots
            Plotly.purge('volumePlot');
            Plotly.purge('treatmentPlot');
            Plotly.purge('timelinePlot');
        }}
        
        function renderVolumePlot(data) {{
            const days = Array.from({{length: data.volume_history.length}}, (_, i) => i + 1);
            const volumes = data.volume_history;
            
            // Find treatment event days for annotations
            const eventMarkers = {{x: [], y: [], text: [], type: []}};
            
            data.events.forEach(ev => {{
                const day = ev.start_day;
                if (day < days.length) {{
                    eventMarkers.x.push(day);
                    eventMarkers.y.push(volumes[day] || 0);
                    eventMarkers.text.push(ev.description + ' (Day ' + day + ')');
                    eventMarkers.type.push(ev.type);
                }}
            }});
            
            const traces = [
                {{
                    x: days,
                    y: volumes,
                    mode: 'lines',
                    name: 'Tumor Volume',
                    line: {{color: '#2c3e50', width: 2}},
                    hovertemplate: 'Day %{{x}}<br>Volume: %{{y:,}} voxels<extra></extra>'
                }}
            ];
            
            // Add event markers
            if (eventMarkers.x.length > 0) {{
                const colors = {{
                    'surgery': '#000000',
                    'radiation': '#e67e22',
                    'chemotherapy': '#27ae60',
                    'immunotherapy': '#8e44ad',
                    'brachytherapy': '#8b4513',
                    'other': '#95a5a6'
                }};
                
                // Group by type for legend
                const types = [...new Set(eventMarkers.type)];
                types.forEach(type => {{
                    const indices = eventMarkers.type.map((t, i) => t === type ? i : -1).filter(i => i !== -1);
                    if (indices.length > 0) {{
                        traces.push({{
                            x: indices.map(i => eventMarkers.x[i]),
                            y: indices.map(i => eventMarkers.y[i]),
                            mode: 'markers',
                            name: type.charAt(0).toUpperCase() + type.slice(1),
                            marker: {{
                                size: 10,
                                color: colors[type] || '#3498db',
                                symbol: 'diamond',
                                line: {{width: 1, color: 'white'}}
                            }},
                            text: indices.map(i => eventMarkers.text[i]),
                            hovertemplate: '%{{text}}<br>Day %{{x}}<br>Volume: %{{y:,}}<extra></extra>',
                            showlegend: true
                        }});
                    }}
                }});
            }}
            
            const layout = {{
                title: {{
                    text: 'Tumor Volume Evolution (voxels)',
                    font: {{size: 16}}
                }},
                xaxis: {{title: 'Day', gridcolor: '#eee'}},
                yaxis: {{title: 'Volume (voxels)', gridcolor: '#eee', rangemode: 'tozero'}},
                hovermode: 'x unified',
                plot_bgcolor: '#fafafa',
                paper_bgcolor: 'white',
                legend: {{orientation: 'h', y: -0.2, x: 0.5, xanchor: 'center'}},
                margin: {{l: 60, r: 20, t: 40, b: 60}},
                height: 380
            }};
            
            try {
                Plotly.newPlot('volumePlot', traces, layout, {{responsive: true, displayModeBar: true}});
                console.log('Volume plot rendered successfully');
            } catch (err) {
                console.error('Error rendering volume plot:', err);
            }
        }}
        
        function renderTreatmentPlot(data) {{
            const days = data.drug_concentration.length;
            const x = Array.from({{length: days}}, (_, i) => i);
            
            const traces = [
                {{
                    x: x,
                    y: data.drug_concentration,
                    mode: 'lines',
                    name: 'Drug Concentration',
                    line: {{color: '#27ae60', width: 2}},
                    fill: 'tozeroy',
                    fillcolor: 'rgba(39, 174, 96, 0.15)',
                    hovertemplate: 'Day %{{x}}<br>Drug Conc.: %{{y:.2f}}<extra></extra>',
                    yaxis: 'y'
                }},
                {{
                    x: x,
                    y: data.radiation_dose.map(d => d * 5), // Scale for visibility
                    mode: 'lines',
                    name: 'Radiation Dose (×5)',
                    line: {{color: '#e67e22', width: 2, dash: 'dot'}},
                    hovertemplate: 'Day %{{x}}<br>RT Dose: %{{y:.1f}} Gy<extra></extra>',
                    yaxis: 'y2'
                }}
            ];
            
            const layout = {{
                title: {{
                    text: 'Treatment Intensity Over Time',
                    font: {{size: 16}}
                }},
                xaxis: {{title: 'Day', gridcolor: '#eee'}},
                yaxis: {{
                    title: 'Drug Concentration',
                    side: 'left',
                    gridcolor: '#eee'
                }},
                yaxis2: {{
                    title: 'Radiation Dose (Gy ×5)',
                    side: 'right',
                    overlaying: 'y',
                    showgrid: false
                }},
                hovermode: 'x unified',
                plot_bgcolor: '#fafafa',
                paper_bgcolor: 'white',
                legend: {{orientation: 'h', y: -0.2, x: 0.5, xanchor: 'center'}},
                margin: {{l: 60, r: 60, t: 40, b: 60}},
                height: 380
            };
            
            Plotly.newPlot('treatmentPlot', traces, layout, {{responsive: true, displayModeBar: true}});
            console.log('Treatment plot rendered successfully');
        } catch (err) {
            console.error('Error rendering treatment plot:', err);
        }}
        
        function renderTimelinePlot(data) {{
            const traces = [];
            const yPositions = {{}};
            let yIdx = 0;
            
            // Define event types and colors
            const colors = {{
                'surgery': '#000000',
                'radiation': '#e67e22',
                'chemotherapy': '#27ae60',
                'immunotherapy': '#8e44ad',
                'brachytherapy': '#8b4513',
                'other': '#95a5a6'
            }};
            const labels = {{
                'surgery': '[Surgery] Surgery',
                'radiation': '[Radiation] Radiation',
                'chemotherapy': '[Chemo] Chemotherapy',
                'immunotherapy': '[Immuno] Immunotherapy',
                'brachytherapy': '[Brachy] Brachytherapy',
                'other': '[Other] Other'
            }};
            
            data.events.forEach((ev, idx) => {{
                const type = ev.type;
                const start = ev.start_day;
                const end = ev.end_day || start;
                const desc = ev.description || type;
                
                if (!yPositions[type]) {{
                    yPositions[type] = yIdx++;
                }}
                
                traces.push({{
                    x: [start, end],
                    y: [yPositions[type], yPositions[type]],
                    mode: 'lines+markers',
                    name: labels[type],
                    legendgroup: type,
                    showlegend: idx === data.events.findIndex(e => e.type === type),
                    line: {{
                        color: colors[type] || '#3498db',
                        width: 12
                    }},
                    marker: {{
                        size: 8,
                        color: colors[type] || '#3498db',
                        symbol: 'square'
                    }},
                    hovertemplate: 
                        '<b>' + labels[type] + '</b><br>' +
                        'Description: ' + desc + '<br>' +
                        'Start: Day ' + start + '<br>' +
                        'End: Day ' + end + '<br>' +
                        'Duration: ' + (end - start) + ' days<extra></extra>',
                    hoverlabel: {{bgcolor: colors[type] || '#3498db'}}
                }});
            }});
            
            const layout = {{
                title: {{
                    text: 'Clinical Treatment Timeline',
                    font: {{size: 16}}
                }},
                xaxis: {{
                    title: 'Days from Diagnosis',
                    gridcolor: '#eee',
                    range: [0, 180]
                }},
                yaxis: {{
                    title: '',
                    showticklabels: false,
                    showgrid: false,
                    range: [-0.5, Object.keys(yPositions).length - 0.5]
                }},
                hovermode: 'closest',
                plot_bgcolor: '#fafafa',
                paper_bgcolor: 'white',
                legend: {{orientation: 'h', y: -0.15, x: 0.5, xanchor: 'center'}},
                margin: {{l: 40, r: 20, t: 40, b: 60}},
                height: 180
            }};
            
            Plotly.newPlot('timelinePlot', traces, layout, {{responsive: true, displayModeBar: true}});
            console.log('Timeline plot rendered successfully');
        } catch (err) {
            console.error('Error rendering timeline plot:', err);
        }}
        
        // Initialize: load first patient when DOM is ready
        function init() {
            console.log('Initializing viewer...');
            console.log('Available patients:', patientIds.length);
            console.log('First patient:', patientIds[0]);
            
            if (patientIds.length > 0) {
                dropdown.value = patientIds[0];
                loadPatient(patientIds[0]);
            } else {
                console.error('No patients available!');
            }
        }
        
        // Run init when DOM is ready
        if (document.readyState === 'loading') {
            document.addEventListener('DOMContentLoaded', init);
        } else {
            init();
        }
    </script>
</body>
</html>
"""
    
    # Build patient options for dropdown
    patient_options = ''.join([f'<option value="{pid}">{pid}</option>' for pid in patient_ids])
    
    # Serialize data for JavaScript
    patient_data_json = json.dumps(patient_data)
    patient_ids_json = json.dumps(patient_ids)
    
    # Replace placeholders in template
    html_content = js_template.replace('{patient_options}', patient_options)
    html_content = html_content.replace('{patient_data_json}', patient_data_json)
    html_content = html_content.replace('{patient_ids_json}', patient_ids_json)
    
    return html_content


def main():
    print("=" * 60)
    print("  UNIFIED 4D PATIENT VIEWER GENERATOR")
    print("=" * 60)
    
    print("\n[LOAD] Loading all patient 4D data...")
    patients = load_all_patient_4d_data()
    print(f"[LOAD] Loaded {len(patients)} patients with 4D simulation data")
    
    if not patients:
        print("[ERROR] No patient data found!")
        return 1
    
    print("\n[BUILD] Creating unified viewer HTML...")
    html_content = create_unified_viewer(patients)
    
    output_path = BASE_OUTPUT / "unified_4d_patient_viewer.html"
    with open(output_path, 'w') as f:
        f.write(html_content)
    
    print(f"[SUCCESS] Unified viewer saved to: {output_path}")
    print(f"[SUCCESS] {len(patients)} patients available in dropdown")
    
    # Open in browser
    import webbrowser
    webbrowser.open(f"file://{output_path.absolute()}")
    
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())