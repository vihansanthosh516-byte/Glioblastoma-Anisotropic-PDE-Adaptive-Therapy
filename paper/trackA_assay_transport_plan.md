# Track A assay control transport stress plan

Simulation only. Use the same paired feature draws, logistic regressions, replicate count, patient counts, artifact shifts and gate threshold from the [previous plan](trackA_assay_gate_plan.md), with fresh seed 20261129. Test planted spatial signal and reversing artifact. Cross control sample sizes 8 and 40 with control sensitivities 0, 0.5 and 1. Sensitivity multiplies the control artifact mean (0.80); the observed paired assay artifact remains unchanged. Control means are zero in spatial signal. Use identical draws for all arms within each setting. Report each setting and inclusion fraction. Keep tissue cohorts closed.

The unchanged gate passes a setting if it improves by at least 0.05 over expression plus composition with planted signal and stays within 0.03 of chance when spatial is an artifact. Overall pass requires every setting to pass. Zero sensitivity is a deliberate failure test for an unobserved artifact. A pass under full sensitivity alone is conditional evidence, not deployment validation.

Compare with unconditional spatial inclusion and unconditional omission using the same paired score. Prior art includes [negative-control normalization with replicate samples](https://pmc.ncbi.nlm.nih.gov/articles/PMC4404308/). The specific gate is a candidate application. Novelty remains unresolved.
