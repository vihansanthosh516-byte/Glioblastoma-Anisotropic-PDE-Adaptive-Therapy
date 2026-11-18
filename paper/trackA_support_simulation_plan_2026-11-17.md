# Track A support threshold simulation plan

Status: post hoc proposal motivated by the [Artzi failure](trackA_artzi_result_2026-11-17.md). This simulation screens a design for future study. Artzi, Yu and HGCC remain spent.

Simulate paired gene order votes directly. Each independent replicate has 10 training patients and 100 fresh test patients, each with a core and edge vote on each of 30 synthetic pairs. Core votes are independent fair signs. For the first five pairs, draw a shared patient response indicator and, when it is on, set core to -1 and edge to +1. Otherwise draw core and edge independently as fair signs. Remaining pairs always use independent fair signs. Response probabilities are 0, 0.45 and 0.75. Use 500 replicates per regime and NumPy seed 20261117.

For each pair choose orientation by larger positive fraction of training patient contrasts, then larger mean signed contrast, then positive orientation. Keep up to five pairs with positive contrasts for at least 70% of training patients, sorted by support, signed mean, then pair index. Average their oriented votes. Empty selection ties. The simulation starts with pair votes and skips gene screening.

Compare with the original all-patient veto and a mean-contrast vote that uses the same cap and keeps pairs with strictly positive oriented mean, ranked by mean then index. All models use the same simulated train/test draws and mean patient paired concordance, scoring ties as 0.5.

Pass only if, in both signal regimes, the candidate averages at least 0.02 above each comparator, and in the null regime its mean is within 0.02 of 0.50. These are simulation screening thresholds. Threshold relaxation has related prior art in robust and stability based selection. No new method claim or real-cohort scoring follows from a simulation pass alone. See the [prior-art inventory](trackA_prior_art_full_2026-11-16.md).
