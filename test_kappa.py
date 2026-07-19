import numpy as np
import sys
sys.path.insert(0, r'C:\Users\avneesh\science-Fair-20206-vs\src')
from dose_response_model import hill_equation, bliss_combination, KAPPA_SYNERGY

drug_a = {'ec50_tumor': 0.8, 'emax_tumor': 0.65, 'hill': 1.2, 'ec50_healthy': 4.0, 'emax_healthy': 0.20}
drug_b = {'ec50_tumor': 2.0, 'emax_tumor': 0.45, 'hill': 1.8, 'ec50_healthy': 20.0, 'emax_healthy': 0.08}
synergy = -0.0180

C1_range = np.linspace(0, 8, 10)
C2_range = np.linspace(0, 25, 10)

print('TI grid with kappa (should show variation across C1 and C2):')
print('C1\\C2', end='')
for C2 in C2_range:
    print(f'{C2:6.1f}', end='')
print()

for C1 in C1_range:
    print(f'{C1:5.1f}', end='')
    for C2 in C2_range:
        e1_t = hill_equation(C1, drug_a['ec50_tumor'], drug_a['emax_tumor'], drug_a['hill'])
        e2_t = hill_equation(C2, drug_b['ec50_tumor'], drug_b['emax_tumor'], drug_b['hill'])
        tumor_eff = bliss_combination(e1_t, e2_t, synergy, C1, C2, KAPPA_SYNERGY)
        
        e1_h = hill_equation(C1, drug_a['ec50_healthy'], drug_a['emax_healthy'], drug_a['hill'])
        e2_h = hill_equation(C2, drug_b['ec50_healthy'], drug_b['emax_healthy'], drug_b['hill'])
        healthy_eff = bliss_combination(e1_h, e2_h, synergy, C1, C2, KAPPA_SYNERGY)
        
        if tumor_eff <= 0 or healthy_eff <= 0:
            ti = -10.0
        else:
            ti = np.log2(tumor_eff / healthy_eff)
        print(f'{ti:6.2f}', end='')
    print()