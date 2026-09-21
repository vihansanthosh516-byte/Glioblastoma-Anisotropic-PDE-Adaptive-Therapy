import pandas as pd

expr = pd.read_csv('data/tcga_gbm_expression.tsv', sep='\t').set_index('sample')
clin = pd.read_csv('data/tcga_gbm_clinical.csv')

targets = ['S100A6', 'S100A11', 'S100A8', 'CCL3L1']
expr_subset = expr.loc[targets]

expr_patients = expr_subset.T
expr_patients.index.name = 'sample'

merged = clin.merge(expr_patients, left_on='sample', right_index=True, how='inner')

print('Merged patients:', len(merged))
print('Columns:', list(merged.columns))
print()
print('Survival summary:')
print('  Events:', int(merged['overall_survival'].sum()))
print('  Median OS:', round(merged['overall_survival_time'].median(), 1), 'days')
print()
print('Expression ranges:')
for g in targets:
    print('  ' + g + ':', round(merged[g].min(), 2), 'to', round(merged[g].max(), 2), '(log2 TPM)')
print()
print('ZNF10* genes in expression file:')
znf_matches = expr.index[expr.index.str.contains('ZNF10', na=False)]
print(list(znf_matches)[:20])