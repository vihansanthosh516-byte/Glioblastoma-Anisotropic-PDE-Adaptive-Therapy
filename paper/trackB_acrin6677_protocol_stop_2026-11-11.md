# Track B ACRIN 6677 protocol stop (2026-11-11, America/Chicago)

The archive SHA-256 matches the [handoff plan](NEXT_STEPS_PLAN_2026-11-12.md). No patient measurement, outcome, or treatment arm was read.

The [handoff plan, Section 5](NEXT_STEPS_PLAN_2026-11-12.md#5-step-4-the-one-patient-run) says only falling histories decide the test because M equals B2 on non-declining histories. The [frozen protocol, Section 5](trackB_acrin6677_protocol_2026-11-12.md#5-endpoint-test-and-decision-rule) compares M with whichever baseline has the lowest cohort median error. If B1 or B3 is selected, M can differ from that baseline on rising histories. Restricting the test to falling histories could change the paired median, p-value, and verdict. The authorized amendment defines the model branches and measurement fallback, but does not resolve this conflict in the analysis population.

Before opening patient tables, clarify that all eligible histories enter the comparison against the selected baseline and that the falling-history statement applies only to differences between M and B2. No amendment was frozen, no code or synthetic tests were run, and no verdict exists.
