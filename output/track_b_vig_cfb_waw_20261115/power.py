"""Power of the pooled CFB + WAW test against persistence only. No patient data is read.

Assumptions: effect 0.027 (output/track_b_code_audit_20261114/loco4_all_baselines_stdout.txt, pooled VIG minus persistence),
per-patient SD of the paired difference 0.154 (output/track_b_vig_lumiere_20261114/result.json). Normal approximation.
Level: this is the second outside run of VIG, so the two-sided level is 0.05 / 2 = 0.025 (97.5% CI).
"""
import numpy as np
from scipy.stats import norm

EFFECT, SD = 0.027, 0.154
for level in (0.05, 0.025):
    zc = norm.ppf(1 - level / 2)
    need = int(np.ceil(((zc + norm.ppf(0.80)) * SD / EFFECT) ** 2))
    row = ", ".join(f"n={n}: {norm.cdf(EFFECT * np.sqrt(n) / SD - zc):.2f}" for n in (100, 120, 158, 186))
    print(f"two-sided level {level}: critical z {zc:.3f} | patients for 80% power {need} | power at {row}")
