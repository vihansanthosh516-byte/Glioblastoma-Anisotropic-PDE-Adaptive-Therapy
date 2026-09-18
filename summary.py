import pandas as pd
from scipy.stats import wilcoxon

for mode in ['standard', 'atlas', 'dki']:
    df = pd.read_csv(f'C:/Users/vihan/Downloads/output/results_{mode}.csv')
    w, p = wilcoxon(df['dsc_aniso'], df['dsc_iso'])
    print(f'{mode}: N={len(df)}, aniso={df["dsc_aniso"].mean():.4f}, iso={df["dsc_iso"].mean():.4f}, delta={df["delta"].mean():.4f}, aniso_wins={sum(df["delta"] > 0)}/{len(df)}, p={p:.6f}')