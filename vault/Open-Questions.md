# Open Questions (updated 2026-10-04)

## Scientific
- [x] Is the high-rho adaptive failure a setpoint artifact? No: THRESHOLD_OFF=0.50 did not restore holidays. 0.80 kept.
- [x] Is the rho > 0.02 rule real? It holds inside the cohort (10/10 vs 0/51, permutation p 5e-5, leave-one-out 98.4%) but n=10 and post hoc. LUMIERE cannot test it (3-4 patients above 0.02).
- [x] Does the Track A to Track B link work? No. Score is 1.0 for all patients.
- [ ] What is the leak-free C-GAT accuracy? Waiting for script 88 (5 folds).
- [ ] How should Track A be reframed? Patient-level LR/RF/scVI-LR give 60-66%; the lookup gives 82.8% on the random split.
- [ ] Why is the tract alignment lower under aniso (0.45 vs 0.50)? Not checked.
- [ ] Do the 42 negative-rho patients belong in the paper as a finding?
- [ ] Should the paper claim progression non-inferiority or final-mass non-inferiority? Data supports progression only.
- [ ] Does the CGGA age effect (HR 1.009) conflict with TCGA (1.031)? Probably a cohort difference (IDH-mutant secondary GBM in CGGA). IDH-wild-type primary gives 1.018.
- [ ] Why is the pooled CGGA multivariate c-index 0.48? Not checked.
- [ ] Should RHUH-GBM (n=40, no MGMT) be analysed, or listed as future work?

## Engineering
- [ ] Script 45 hardcoded "8 patients" strings at lines about 1711 and 1775.
- [ ] Vault RAG index needs a re-run after note edits.
- [ ] `output/cgat/gat_edge_*.npy` are git-ignored and missing; script 88 rebuilds the graph.

## Writing
- [ ] Track B section order: resistance first, then dose sparing, then the stratified rule?
- [ ] Where to put the MCP/RAG tooling (methods or supplement)?
- [ ] Track A claims: reframe or remove after script 88.
