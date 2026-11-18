# Track A order reversal simulation plan

Simulation follow-up to the [block voting result](trackA_block_simulation_result_2026-11-17.md). Use the prior generator, methods, depths and sample sizes from the [frozen plan](trackA_block_simulation_plan_2026-11-17.md) with fresh seed 20261119. Keep all real cohorts closed.

Run the unequal-depth diffuse shared shift, unequal-depth null, and two unequal-depth reversal regimes. In the reversal regimes replace the shared edge shift with +1 to gene indices 0, 2, 4 and -1 to indices 1, 3, 5. Apply it to all edge blocks in diffuse reversal and one uniformly sampled edge block per patient in localized reversal. The changes create order-sensitive signal without forcing every patient to reverse.

Candidate: combine block logistic and soft block voting by averaging their continuous paired margins after dividing each by its root mean square margin on training patients (use scale 1 for a zero margin vector). This replaces the sign-only proposal before any draws, since disagreement between two nonzero signs would always force a tie. Score ties as in the prior plan. This prediction requires the paired specimens together. Compare against every prior arm on identical draws. Pass only if mean concordance exceeds every comparator by at least 0.01 in both reversal regimes, is within 0.02 of 0.50 under null, and trails the best arm by at most 0.02 under shared shifts. This is a fixed ensemble screening experiment. Novelty remains unestablished.

