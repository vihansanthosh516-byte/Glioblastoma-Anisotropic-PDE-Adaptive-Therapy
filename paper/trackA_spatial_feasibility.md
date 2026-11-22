# Track A spatial feasibility gate

Simulation with exact enumeration. Use a ring of eight malignant positions with a myeloid bridge between each consecutive pair. Four malignant positions carry state A, and four carry state B. Each malignant position has the same degree. Each map therefore has the same state counts and direct state-to-myeloid contact totals.

Assign the synthetic core pattern AAAABBBB and edge pattern ABABABAB. Count the fraction of myeloid bridges connecting malignant positions of opposite state. Enumerate all assignments of four A marks to eight positions. Compute the exact minimum fraction of changed malignant marks needed for the edge score to reach or fall below the core score, holding core fixed. Composition, direct contact totals and graph geometry serve as tie controls.

Pass feasibility if edge exceeds core, all controls tie, a positive exact disruption budget exists, and changing an unused depth annotation alone leaves the statistic unchanged. The annotation check is a software control, not a realistic model of depth-induced label errors. Because synthetic labels are assigned by design, this only verifies identifiability in a toy geometry. It cannot establish prediction, biological validity or novelty. No patient data. See [parallel research](trackA_parallel_research.md).
