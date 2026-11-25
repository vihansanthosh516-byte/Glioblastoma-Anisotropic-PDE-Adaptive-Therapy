# Track A region: one labelled section to the rest of the tumor (DRAFT plan, 2026-11-24)

Status: DRAFT. Not locked. No histopathology label values have been read; only whether each section has any non-zero label (coverage). Prior art: `trackA_few_patient_designs_prior_art_2026-11-24.md`.

## Coverage (presence only; 91 of 97 files complete, 6 AT15 files still downloading)

| Donor | Labelled sections | Sites (from section name `site_replicate`) | Spots in labelled sections |
|---|---:|---|---:|
| AT10 | 8 of 8 | 4 sites x 2 replicates | 26,182 |
| AT3 | 8 of 8 | 4 sites x 2 replicates | 14,589 |
| AT4 | 8 of 8 | 4 sites x 2 replicates | 20,638 |
| AT13 | 4 of 8 | sites 1 and 4 x 2 replicates | 13,950 |
| AT15 | 2 of 3 so far | different naming; 6 files pending | 20,300 |
| AT6 | 2 of 8 | site 4 only x 2 replicates | 5,390 |

Donors with at least two labelled sites (needed for cross-site transfer): AT10, AT3, AT4, AT13 (4), plus AT15 if its pending files add a second site. AT6 supports replicate transfer only.

Label source: GBM-Space README, "percentage of spot overlap with IvyGAP histopathology annotations". Eight categories: cellular tumor, infiltrating tumor (grey matter), leading edge (grey matter), microvascular proliferation, hyperplastic blood vessels, necrosis, perinecrotic zone, pseudopalisading cells. Who drew the annotations is not yet verified from the paper.

## Claim (within-patient only)

Given one pathologist-labelled section, a method maps core vs infiltrating/edge tissue on the same patient's sections from other sampling sites, from expression alone, with per-spot confidence that is honest on the new site.

## Proposed task definition (to fix before lock)

- Spot label: cellular tumor if CT overlap >= 50%; edge if IT + LE overlap >= 50%; other spots excluded from scoring.
- Source: one labelled section. Targets: every labelled section from a different site of the same donor.
- Unit of replication: donor. Every source-target pair is scored; donor score = mean over pairs.

## Candidate idea (to design and check on generated data before lock)

Per-patient classifier trained on the source section, with (a) within-section rank normalization to absorb section-level shifts, and (b) split-conformal abstention calibrated on the source site's replicate section, so the method reports which target spots it is confident about. The new part to test by ablation is the replicate-calibrated confidence under site shift.

## Comparators (same pairs)

1. Alignment-based label transfer (PASTE2).
2. Ivy GAP marker scoring (no training).
3. GBM-Space spatial niche abundances, best single niche.
4. Plain classifier on the source section (no normalization, no abstention).
5. GBM-Space cell-state abundances classifier.

## Open questions before lock

- Read the GBM-Space Methods for who annotated, which sections, and resolution.
- Primary metric: AUROC (ranking) vs selective accuracy at fixed coverage (tests the confidence part). One must be primary.
- Pass rule across only 4-5 donors: e.g., candidate beats every comparator in every donor (sign test cannot reach p < 0.05 with 4 donors; state this limit).
- Check edge-spot counts per section after the lock rule is fixed (counts only).

## Weak points known now

- 4 (maybe 5) donors with cross-site labels; results are a within-patient demonstration, not generalization to new patients.
- Ivy GAP categories may miss some tumor structure (Greenwald 2024).
- GBM-Space has already studied this atlas; the contribution must be the transfer method, not the biology.
