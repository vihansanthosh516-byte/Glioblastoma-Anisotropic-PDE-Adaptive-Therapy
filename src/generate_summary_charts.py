import os
import matplotlib.pyplot as plt
import pandas as pd

# Ensure output directory exists
os.makedirs("outputs", exist_ok=True)

# 1. Ingest real validation metrics data
data = {
    'PatientID': ['0041', '0045', '0022', '0003 (TP1-TP2)', '0013', '0037', '0003 (TP2-TP5)'],
    'Delta_Days': [91, 70, 52, 56, 98, 105, 140],
    'V0': [3985, 47734, 154238, 84540, 43633, 60937, 85818],
    'V1': [9666, 94403, 222866, 85818, 56472, 90226, 193095],
    'Rho': [0.0098, 0.0105, 0.0087, 0.0050, 0.0050, 0.0050, 0.0067],
    'D': [0.0049, 0.0054, 0.0040, 0.0010, 0.0010, 0.0010, 0.0024],
    'Status': ['Best Fit', 'Excellent', 'Good', 'Slow Growth', 'Slow Growth', 'Slow Growth', 'Good']
}

df = pd.DataFrame(data)

# Separate active progressing tumors from slow-growth boundary cases
df_active = df[~df['Status'].str.contains('Slow')].copy()

# 2. Generate side-by-side summary figures
fig, axes = plt.subplots(1, 2, figsize=(14, 6), constrained_layout=True)

# Plot A: Tumor Volume Progression (V0 vs V1)
axes[0].scatter(df['V0'], df['V1'], c=df['Delta_Days'], cmap='coolwarm', s=100, edgecolors='black', zorder=3)
max_val = max(df['V0'].max(), df['V1'].max()) * 1.05
axes[0].plot([0, max_val], [0, max_val], 'k--', alpha=0.5, label='No Growth Line ($V_0 = V_1$)')
axes[0].set_title(r'Longitudinal Volume Expansion ($V_0 \rightarrow V_1$)', fontsize=12, fontweight='bold')
axes[0].set_xlabel('Baseline Volume $V_0$ ($\text{mm}^3$)', fontsize=10)
axes[0].set_ylabel('Follow-up Volume $V_1$ ($\text{mm}^3$)', fontsize=10)
axes[0].grid(True, linestyle=':', alpha=0.6)
axes[0].legend(loc='upper left')

# Annotate points on volume plot
for _, row in df.iterrows():
    axes[0].annotate(f"Pt {row['PatientID']}", (row['V0'], row['V1']), textcoords="offset points", xytext=(5,5), ha='left', fontsize=8)

# Plot B: Estimated Biomechanical Parameter Space ($\rho$ vs $D$) - Note the 'r' prefix for raw string!
scatter = axes[1].scatter(df_active['Rho'], df_active['D'], c=df_active['Delta_Days'], cmap='viridis', s=120, edgecolors='black', zorder=3)
axes[1].set_title(r'Inferred Biomechanical Parameters ($\rho$ vs $D$)', fontsize=12, fontweight='bold')
axes[1].set_xlabel(r'Proliferation Rate $\rho$ ($\text{day}^{-1}$)', fontsize=10)
axes[1].set_ylabel(r'Diffusivity $D$ ($\text{mm}^2/\text{day}$)', fontsize=10)
axes[1].grid(True, linestyle=':', alpha=0.6)
cbar = fig.colorbar(scatter, ax=axes[1])
cbar.set_label(r'Time Interval ($\Delta$ Days)', fontsize=9)

# Annotate active parameter points
for _, row in df_active.iterrows():
    axes[1].annotate(f"Pt {row['PatientID']}", (row['Rho'], row['D']), textcoords="offset points", xytext=(6,6), ha='left', fontsize=9, fontweight='semibold')

# Save output chart
output_path = 'outputs/validation_summary_dashboard.png'
plt.savefig(output_path, dpi=300)
print(f"Summary validation plots successfully compiled and saved to: {output_path}")