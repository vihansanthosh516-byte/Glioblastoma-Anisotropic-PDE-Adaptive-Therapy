import pandas as pd
import numpy as np
from scipy.stats import wilcoxon
import json

output_dir = "C:/Users/vihan/Downloads/output"

# Load all results
standard = pd.read_csv(f"{output_dir}/results_standard.csv")
atlas = pd.read_csv(f"{output_dir}/results_atlas.csv")
dki = pd.read_csv(f"{output_dir}/results_dki.csv")

standard['mode'] = 'Standard DTI'
atlas['mode'] = 'Atlas DTI'
dki['mode'] = 'DKI'

df = pd.concat([standard, atlas, dki], ignore_index=True)
df = df[df['status'] == 'success']

# Set random seed for reproducibility
np.random.seed(42)

# Get unique patients
patients = df['patient'].unique()
print(f"Total patients: {len(patients)}")

# 42/20 split
n_train = 42
n_test = 20

# Random shuffle
shuffled = np.random.permutation(patients)
train_patients = shuffled[:n_train]
test_patients = shuffled[n_train:n_train+n_test]

print(f"Train patients ({len(train_patients)}): {sorted(train_patients)}")
print(f"Test patients ({len(test_patients)}): {sorted(test_patients)}")

# Split data
train_df = df[df['patient'].isin(train_patients)]
test_df = df[df['patient'].isin(test_patients)]

# Analyze train set
print("\n" + "="*60)
print("TRAIN SET (n=42)")
print("="*60)
for mode in ['Standard DTI', 'Atlas DTI', 'DKI']:
    subset = train_df[train_df['mode'] == mode]
    aniso_mean = subset['dsc_aniso'].mean()
    iso_mean = subset['dsc_iso'].mean()
    delta_mean = subset['delta'].mean()
    win_pct = (subset['delta'] > 0).sum() / len(subset) * 100
    try:
        w, p = wilcoxon(subset['dsc_aniso'], subset['dsc_iso'])
    except:
        w, p = None, None
    print(f"  {mode}: aniso={aniso_mean:.4f}, iso={iso_mean:.4f}, delta={delta_mean:.4f}, wins={win_pct:.0f}%, p={p:.6f}")

# Analyze test set
print("\n" + "="*60)
print("TEST SET (n=20)")
print("="*60)
for mode in ['Standard DTI', 'Atlas DTI', 'DKI']:
    subset = test_df[test_df['mode'] == mode]
    aniso_mean = subset['dsc_aniso'].mean()
    iso_mean = subset['dsc_iso'].mean()
    delta_mean = subset['delta'].mean()
    win_pct = (subset['delta'] > 0).sum() / len(subset) * 100
    try:
        w, p = wilcoxon(subset['dsc_aniso'], subset['dsc_iso'])
    except:
        w, p = None, None
    print(f"  {mode}: aniso={aniso_mean:.4f}, iso={iso_mean:.4f}, delta={delta_mean:.4f}, wins={win_pct:.0f}%, p={p:.6f}")

# Check consistency (do train/test agree on which mode wins?)
print("\n" + "="*60)
print("CONSISTENCY CHECK")
print("="*60)
for mode in ['Standard DTI', 'Atlas DTI', 'DKI']:
    train_subset = train_df[train_df['mode'] == mode]
    test_subset = test_df[test_df['mode'] == mode]
    train_delta = train_subset['delta'].mean()
    test_delta = test_subset['delta'].mean()
    consistent = (train_delta > 0) == (test_delta > 0)
    print(f"  {mode}: train_delta={train_delta:.4f}, test_delta={test_delta:.4f}, consistent={consistent}")

# Save split
split_info = {
    'train_patients': sorted(train_patients.tolist()),
    'test_patients': sorted(test_patients.tolist()),
    'random_seed': 42,
    'n_train': n_train,
    'n_test': n_test
}

with open(f"{output_dir}/cv_split.json", "w") as f:
    json.dump(split_info, f, indent=2)

train_df.to_csv(f"{output_dir}/cv_train_results.csv", index=False)
test_df.to_csv(f"{output_dir}/cv_test_results.csv", index=False)

print(f"\nSplit saved to {output_dir}/cv_split.json")
print(f"Train results saved to {output_dir}/cv_train_results.csv")
print(f"Test results saved to {output_dir}/cv_test_results.csv")

# Full dataset comparison
print("\n" + "="*60)
print("FULL DATASET (n=62) - FOR COMPARISON")
print("="*60)
for mode in ['Standard DTI', 'Atlas DTI', 'DKI']:
    subset = df[df['mode'] == mode]
    aniso_mean = subset['dsc_aniso'].mean()
    iso_mean = subset['dsc_iso'].mean()
    delta_mean = subset['delta'].mean()
    win_pct = (subset['delta'] > 0).sum() / len(subset) * 100
    w, p = wilcoxon(subset['dsc_aniso'], subset['dsc_iso'])
    print(f"  {mode}: aniso={aniso_mean:.4f}, iso={iso_mean:.4f}, delta={delta_mean:.4f}, wins={win_pct:.0f}%, p={p:.6f}")