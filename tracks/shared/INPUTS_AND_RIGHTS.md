# Integrated-phase inputs, rights, and checksums

## Allowed input boundary

The integrated simulation will use synthetic signal settings inspired by the Track A question; it will not use observed A5 metrics to parameterize signal accuracy. A Track B aggregate may inform a history-only MRI scale only if the exact Duke agreement permits that use. The integrated simulation will not use or reopen raw expression matrices, patient MRI, masks, ADC maps, or patient-level files. It will not link patient identities across Tracks A and B. This restriction applies to the new integrated phase; the completed Duke exploratory analysis already read target masks.

## Before an input is used

For each permitted aggregate, create a local input-manifest row containing: stable input name, source URL/record, local path, content type, SHA-256, source/license record, written agreement or permission basis, allowed operations, access date, and whether the row may be committed. Hash only the exact file approved for this phase. Keep restricted paths, patient identifiers, and any non-shareable derived information outside Git.

Separate dataset terms from article terms. A public paper license does not automatically grant rights to the dataset, derived summaries, trained weights, or report figures. When the Duke agreement is silent or unclear about a summary, omit that summary and use the predeclared standard-normal MRI-like fallback. Do not infer permission from user access alone.

## Checksums

The integrated run manifest will hash the exact permitted aggregate files and all code/configuration inputs before evaluation. It will record the hash algorithm, relative input key, file size, and digest. The manifest must not contain private absolute paths or identifying metadata. Raw patient files are outside the integrated-phase input set and will not be hashed as part of this plan.

The current planning manifest is [`configs/trackABC_input_manifest_v2.csv`](../../configs/trackABC_input_manifest_v2.csv). It retains the local A5 summary's checksum for provenance only; those observed metrics are not integrated-model inputs. It does not include a Duke summary. The template for additional permitted inputs is [`configs/trackABC_input_manifest_template.csv`](../../configs/trackABC_input_manifest_template.csv).

## Source records

- Track A source and protocol: [GEO GSE117891](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117891) and [`paper/trackA_yu2020/`](../../paper/trackA_yu2020/).
- Track B source record: [Duke Zenodo record 20548528](https://zenodo.org/records/20548528). The accepted access agreement and current record terms govern every use.
- Simulation and aggregate literature sources: see the [integrated plan bibliography](../../paper/trackABC_integrated_research_plan_master_2026-11-07.md#bibliography).
