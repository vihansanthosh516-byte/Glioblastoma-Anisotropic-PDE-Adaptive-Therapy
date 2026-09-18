import json
with open('output/batch_digital_twins/PatientID_0041/treatment_schedule.json') as f:
    s = json.load(f)
print('Events:', len(s.get('events', [])))
for e in s.get('events', []):
    print(f'  {e["type"]}: {e["start_day"]}-{e["end_day"]} {e.get("description","")}')