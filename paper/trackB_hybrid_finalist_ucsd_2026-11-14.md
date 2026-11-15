# Track B hybrid finalist UCSD run — locked before running (2026-11-14)

Dev results (all 6 failed the dev rule, 99.64% CI): moe 1.453, compartment 1.463, pkpd 1.473, stack 1.474, deconv 1.480, switch 1.590 (patient MAE; pop_rate 1.528, v8 1.449). Files: output/track_b_hybrid_<name>/result.json.

Per protocol the finalist is the lowest dev MAE: **moe** (gated mixture of experts). It runs ONCE on UCSD with the rules of `paper/trackB_v8_ucsd_external_protocol_2026-11-14.md` (fit on all MU, comparator pop_rate, 95% CI, same pairs). moe also uses v6 edge-shape and v7 brightness features; on UCSD they are computed with the same formulas from the input session's BraTS seg, T1post and FLAIR.

Honesty: no candidate passed dev, and moe did not beat v8 on dev. UCSD was used once already (v8). This is a second exposure; reported as such whatever the result.
