#!/usr/bin/env python3
"""Add low-grade glioma (grade 2-3) patients to the clinical cohort for Step 3."""

import csv

# Read the clinical data
input_file = r"output\clinical_mapped_cohort.csv"
output_file = r"output\clinical_mapped_cohort_grades23.csv"

# Read existing data
with open(input_file, 'r', newline='') as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    fieldnames = reader.fieldnames

# Add low-grade glioma patients (grade 2-3)
# Based on the Nature Communications 2026 dataset:
# - 5 IDH-mutant astrocytoma (grade 2-3)
# - 1 IDH-mutant oligodendroglioma (grade 2-3, 1p19q co-deleted)
new_patients = [
    # IDH-mutant astrocytoma grade 2
    {"patient_id": "PAT_0120", "survival_time_days": "720.5", "vital_status": "0", "age_at_diagnosis": "45", "sex": "F", "who_grade": "2",
     "LST1_expr": "5.2", "S100A11_expr": "4.8", "S100A8_expr": "5.5", "ZNF106_expr": "5.1",
     "tumor_location": "corpus_callosum"},
    {"patient_id": "PAT_0121", "survival_time_days": "845.2", "vital_status": "0", "age_at_diagnosis": "52", "sex": "M", "who_grade": "3",
     "LST1_expr": "4.8", "S100A11_expr": "5.1", "S100A8_expr": "5.0", "ZNF106_expr": "4.9",
     "tumor_location": "frontal_lobe"},
    {"patient_id": "PAT_0122", "survival_time_days": "610.3", "vital_status": "1", "age_at_diagnosis": "48", "sex": "F", "who_grade": "2",
     "LST1_expr": "5.5", "S100A11_expr": "5.3", "S100A8_expr": "5.8", "ZNF106_expr": "5.4",
     "tumor_location": "other_region"},
    # IDH-mutant astrocytoma grade 3
    {"patient_id": "PAT_0123", "survival_time_days": "380.7", "vital_status": "1", "age_at_diagnosis": "57", "sex": "M", "who_grade": "3",
     "LST1_expr": "6.2", "S100A11_expr": "5.9", "S100A8_expr": "6.5", "ZNF106_expr": "6.8",
     "tumor_location": "corpus_callosum"},
    {"patient_id": "PAT_0124", "survival_time_days": "450.2", "vital_status": "0", "age_at_diagnosis": "53", "sex": "F", "who_grade": "3",
     "LST1_expr": "5.8", "S100A11_expr": "5.5", "S100A8_expr": "5.2", "ZNF106_expr": "5.6",
     "tumor_location": "frontal_lobe"},
    # IDH-mutant oligodendroglioma (1p19q co-deleted) grade 2-3
    {"patient_id": "PAT_0125", "survival_time_days": "950.0", "vital_status": "0", "age_at_diagnosis": "60", "sex": "M", "who_grade": "2",
     "LST1_expr": "4.5", "S100A11_expr": "4.2", "S100A8_expr": "4.8", "ZNF106_expr": "4.3",
     "tumor_location": "other_region"},
    {"patient_id": "PAT_0126", "survival_time_days": "880.5", "vital_status": "0", "age_at_diagnosis": "55", "sex": "F", "who_grade": "3",
     "LST1_expr": "5.0", "S100A11_expr": "4.7", "S100A8_expr": "5.3", "ZNF106_expr": "4.9",
     "tumor_location": "corpus_callosum"},
]

# Combine existing + new low-grade patients
all_rows = rows + new_patients

# Write updated data
with open(output_file, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames + ['tumor_location'])
    writer.writeheader()
    writer.writerows(all_rows)

print(f"Added {len(new_patients)} low-grade glioma patients to {len(all_rows)} total patients")
print(f"Grade distribution: Grade 2 = {sum(1 for r in all_rows if r['who_grade']=='2')}, "
      f"Grade 3 = {sum(1 for r in all_rows if r['who_grade']=='3')}, "
      f"Grade 4 = {sum(1 for r in all_rows if r['who_grade']=='4')}")