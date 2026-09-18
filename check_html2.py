with open('output/batch_digital_twins/unified_4d_patient_viewer.html', 'r') as f:
    c = f.read()

# Check dropdown options
import re
dropdown_section = re.search(r'<select id="patientDropdown">(.*?)</select>', c, re.DOTALL)
if dropdown_section:
    options = re.findall(r'<option value="([^"]*)">([^<]*)</option>', dropdown_section.group(1))
    print(f'Number of dropdown options: {len(options)}')
    if options:
        print(f'First 5 options: {options[:5]}')

# Check patientData structure
patient_data_match = re.search(r'const patientData = ({.*?});', c, re.DOTALL)
if patient_data_match:
    pd_str = patient_data_match.group(1)[:500]
    print(f'patientData starts with: {pd_str}')

# Check patient_ids
patient_ids_match = re.search(r'const patientIds = (\[.*?\]);', c, re.DOTALL)
if patient_ids_match:
    pi_str = patient_ids_match.group(1)[:200]
    print(f'patientIds starts with: {pi_str}')

# Check loadPatient function for issues
load_match = re.search(r'function loadPatient\(pid\) \{(.*?)^\s*\}', c, re.DOTALL | re.MULTILINE)
if load_match:
    func_body = load_match.group(1)[:500]
    print(f'loadPatient body starts with: {func_body}')

# Check if infoPanel is being set to grid
info_panel = re.search(r'infoPanel\.style\.display\s*=\s*[\'"]grid[\'"]', c)
print(f'infoPanel.style.display = grid: {bool(info_panel)}')

# Check renderVolumePlot call
render_vol = re.search(r'renderVolumePlot\(data\)', c)
print(f'renderVolumePlot called: {bool(render_vol)}')

# Check Plotly.newPlot for volumePlot
vol_newplot = re.search(r"Plotly\.newPlot\('volumePlot'", c)
print(f"Plotly.newPlot('volumePlot'): {bool(vol_newplot)}")