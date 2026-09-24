import json

with open('output/mu_glioma_cohort.json') as f:
    cohort = json.load(f)

print(f'Total patients: {len(cohort)}')
print(f'Longitudinal (>=2 timepoints): {sum(1 for p in cohort if p["has_longitudinal_pair"])}')
print(f'With treatment timing: {sum(1 for p in cohort if p["treatment_schedule"].get("surgery_day") is not None)}')
print()

print('Sample patient 1:')
p = cohort[0]
print(f'  {p["patient_id"]}')
for tp in p['timepoints']:
    print(f'    TP{tp["number"]}: day={tp["day_from_diagnosis"]}, vol={tp["volume_mm3"]}')
print(f'  Treatment: {p["treatment_schedule"]}')
print()

print('Sample patient 2:')
p = cohort[1]
print(f'  {p["patient_id"]}')
for tp in p['timepoints']:
    print(f'    TP{tp["number"]}: day={tp["day_from_diagnosis"]}, vol={tp["volume_mm3"]}')
print(f'  Treatment: {p["treatment_schedule"]}')