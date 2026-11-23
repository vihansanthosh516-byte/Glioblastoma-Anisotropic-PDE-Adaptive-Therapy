# Track A assay control gate result

**Simulation. The assay control gate passed every frozen check.** With planted spatial signal, it scored 0.817 versus 0.682 for composition plus expression. Under the reversing artifact, it rose from 0.216 for always including spatial to 0.496. In the mixed regime it scored 0.687 versus 0.688 for composition plus expression. The gate kept spatial in every pure spatial-signal replicate and excluded it in every pure artifact replicate. [Saved result](../output/trackA_assay_gate/result.json) | [Frozen plan](trackA_assay_gate_plan.md)

The control measurements were independent of biological samples, but the simulator gave them a direct readout of a large training artifact. The gate restored the baseline by removing the spatial feature. It did not recover useful spatial information in mixed cases. A stronger test should vary control sample size and simulate artifact directions and shifts the controls do not capture. No patient data were opened.
