with open('output/batch_digital_twins/4d_3d_viewer_PatientID_0041.html', 'r') as f:
    c = f.read()
print('File size:', len(c))
print('Has plotly:', 'plotly-2.27.0.min.js' in c)
print('Has plot div:', 'id="plot"' in c)
print('Has snapshot files:', 'snapshotFiles' in c)
print('Has dayLabels:', 'dayLabels' in c)
print('Has frames:', 'framesLoaded' in c)