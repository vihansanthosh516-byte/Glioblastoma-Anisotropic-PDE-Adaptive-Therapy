# Track A spatial stress simulation

This fresh simulation extends the [exact toy gate](trackA_spatial_feasibility.md). It tests a fixed spatial score. No patient data. Use the portable JavaScript runner written for this plan.

Use seed 20261124 and a uint32 linear congruential generator (multiplier 1664525, increment 1013904223). Simulate 1000 independent pairs per regime. Each pair shares a graph of 24 malignant positions on a cycle, with a myeloid bridge on each cycle edge. Independently add a bridge from i to i+3 modulo 24 with probability 0.25 for each i. This varies neighborhood connectivity, without modeling physical coordinates. Rotate both latent maps by a shared uniform offset. Core marks are a contiguous half cycle. Edge marks alternate in signal regimes, and copy core in the null.

For each region, independently mark each position low depth with probability 0.25 for core and 0.75 for edge. Flip state independently with probability 0.20 at low depth and 0.02 otherwise. Run clean signal (no flips), noisy signal and noisy null in that order. The noisy null deliberately couples label noise to region, but has no latent biological difference.

For each observed map compute the fraction of bridges with opposite endpoint states minus 2*nA*nB/(24*23), the expected fraction under global state permutation. Larger values predict edge. Compare against unadjusted bridge fraction, observed fraction A, degree-weighted fraction A and geometry-only tie. Also evaluate a depth-matched score that uses only bridges joining two high-depth positions and subtracts the same permutation expectation calculated over the high-depth positions. With no such bridges, return zero. This is the candidate. All methods score the same patients by edge greater than core, ties half.

Pass if the candidate exceeds every baseline by at least 0.01 under noisy signal, scores at least 0.80 under clean signal, and is within 0.03 of chance under noisy null. Report the rate of maps with no high-depth bridges. This test may expose loss of coverage under unequal depth. It cannot establish real-cohort performance or novelty.
