# Open Questions

## Scientific
- [x] Is the high-rho adaptive failure a setpoint artifact? Tested at THRESHOLD_OFF=0.50: 43% drug, 5.8 holidays (vs 32%, 48.9 at 0.80). Holidays not restored, artifact not supported. 0.80 retained.
- [ ] Do the 42 negative-rho responders belong in the paper as a 
      positive finding (TMZ response rate ~41%)?
- [ ] Should the paper claim "progression non-inferiority" or "final-mass 
      non-inferiority"? Current data supports progression only.
- [ ] Do we extend the horizon past 500 days to see if MTD eventually progresses?

## Engineering  
- [ ] Scripts 46, 47, 48 need per-patient patches (same pattern as 43, 44)
- [ ] Script 45 hardcoded "8 patients" strings at lines ~1711, ~1775
- [ ] Vault RAG index needs a re-run after each note edit

## Writing
- [ ] Track B section structure — resistance first, then dose sparing, then 
      the stratified rule?
- [ ] Where to put the MCP/RAG tooling in the paper (methods? supplementary?)