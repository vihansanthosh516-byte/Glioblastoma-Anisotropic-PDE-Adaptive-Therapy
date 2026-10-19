> **HISTORICAL BACKGROUND:** October 4 old-GBM questions, corrected 2026-10-19; current goal/status is START_HERE.md. Resolved items below cite later ledger entries. Remaining questions are proposals, not instructions to run patient analyses. Original progression-noninferiority suggestion is unsupported by the shorter adaptive TTP in master_cohort_summary; see corrected archived synthesis and ledger 2j.

# Open Questions (updated 2026-10-04)

## Scientific
- [x] Is the high-rho adaptive failure a setpoint artifact? In the assumed-treatment simulation, the 0.50 threshold did not support that explanation; this does not establish biological causation. 0.80 kept historically; source/provenance limits in Decisions.md and ledger 2j.
- [x] Is the rho > 0.02 rule real? It holds inside the cohort (10/10 vs 0/51, permutation p 5e-5, leave-one-out 98.4%) but n=10 and post hoc. LUMIERE cannot test it (3-4 patients above 0.02).
- [x] Does the Track A to Track B link work? No. Ledger A24 reports score 1.0 for all 61; exact linked-score artifact/key remains missing here.
- [x] Script 88 completed: clean PCA graph mean accuracy 63.8% (SD 13.3); cVAE no-label-edge 74.2% retains upstream feature leakage, published-edge 75.8% retains label leakage. Source: `output/cgat_leak_free.json:variants`, ledger 2o.
- [ ] How should Track A be reframed? Patient-level LR/RF/scVI-LR give 60-66%; the lookup gives 82.8% on the random split.
- [ ] Why is the tract alignment lower under aniso (0.45 vs 0.50)? Not checked.
- [ ] Do the 42 negative-rho patients belong in the paper as a finding?
- [ ] Should the paper claim progression non-inferiority or final-mass non-inferiority? Original interpretation withdrawn: no progression non-inferiority is established; adaptive TTP was shorter in the saved master summary.
- [ ] Does the CGGA age effect (HR 1.009) conflict with TCGA (1.031)? The cause is unproved; subgroup restriction changes the estimate (ledger 2k), not a demonstrated explanation. IDH-wild-type primary gives 1.018.
- [ ] Why is the pooled CGGA multivariate c-index 0.48? Not checked.
- [ ] Historical RHUH question superseded by H-2 scoring/H-4 training (ledger 2ax/2bf); H-2 now spent. MGMT-field availability needs a dated source dictionary.

## Engineering
- [ ] Script 45 still has hardcoded "8 patients" report strings (current lines about 1739/1803, checked 2026-10-19); source correction tracked under chunk 057, not resolved by this note.
- [ ] Separate vault RAG index freshness is unverified; the refreshed graphify index is not evidence that vault RAG was rebuilt.
- [ ] `output/cgat/gat_edge_*.npy` are git-ignored and missing; script 88 clean-PCA variant builds its own kNN graph; cVAE variants still reference those missing original files.

## Writing
- [ ] Track B section order: resistance first, then dose sparing, then the stratified rule?
- [ ] Where to put the MCP/RAG tooling (methods or supplement)?
- [ ] Track A claims: use completed patient-fold results (ledger 2o), preserving feature/label leakage and donor confounding caveats.
