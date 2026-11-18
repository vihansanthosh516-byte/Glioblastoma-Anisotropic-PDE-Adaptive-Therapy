# Track A tissue block voting result

**Simulation. Soft block voting failed the frozen rule.** In unequal-depth diffuse signal, soft votes scored 0.494 and block logistic regression scored 0.9735. In localized signal they scored 0.499 and 0.606, respectively. The null check passed. [Saved results](../output/trackA_block_simulation_20261117/result.json) | [Frozen plan](trackA_block_simulation_plan_2026-11-17.md)

The simulator shifted both genes within every affected pair by the same amount. This preserves their latent within-block order, even though the shift changes abundance relative to other genes. The simulation shows a limitation of the fixed pair representation under shared shifts. It cannot establish how well block voting works when real order reversals exist. The next simulation needs fresh draws with separate common-shift and order-reversal regimes, fixed before scoring. Keep paired logistic regression as a comparator. No patient data were used.
