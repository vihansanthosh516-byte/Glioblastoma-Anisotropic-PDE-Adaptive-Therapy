# Track A positive-control mismatch result

**Simulation. A passing positive control did not ensure coverage of the actual assay artifact.** The prior rule qualified all mismatched replicates and scored 0.542 versus 0.683 without spatial. The audit-dependent fallback matched 0.683 and passed the frozen checks. Its alignment flag came from the simulator. This is an upper-bound illustration of the value of independent qualification, not evidence that misalignment can be detected on real samples. [Saved results](../output/trackA_probe_mismatch/result.json) | [Frozen plan](trackA_probe_mismatch_plan.md)

The next data milestone is to establish whether paired technical acquisitions and independently annotated regions exist together in a fresh cohort. A practical audit must document which perturbations were challenged and whether each challenge affects the actual features used by the predictor. Keep the mismatch result visible in the paper's limitations or linked development record.
