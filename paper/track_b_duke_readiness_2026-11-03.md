# Track B Duke readiness (2026-11-03)

## Current status

The user reports that Duke granted download access. The exact accepted license and written grant are not saved in the project. Keep access and derived-output rights open until those terms are recorded.

The local `C:\Users\vihan\DATASET_FINAL.zip` was 76,144,640 bytes at the latest check and is growing. The Zenodo record API snapshot in `C:\Users\vihan\Downloads\20548528.json` lists one `DATASET_FINAL.zip` file, size 17,043,569,844 bytes, MD5 `99580fb8ccb3daed186dc6afe5b5bf2f`. The download is incomplete. No archive member has been opened. No MRI, mask, landmark, or patient-label voxel data has been read. No patient scoring has occurred.

The user authorized metadata-only inventory: archive filenames, anonymized patient/session identifiers, scan counts and dates, treatment/IDH metadata. This scope excludes MRI, mask, and landmark-coordinate contents.

## Source facts to carry forward

The [Zenodo record](https://zenodo.org/records/20548528) describes 61 patients and 503 sessions, with 3–18 sessions per patient. It lists 1 mm atlas-standardized NIfTI images and automated tumor segmentations. It describes labels 1 (necrotic core), 2 (nonenhancing T2-FLAIR abnormality), 3 (enhancing tumor), and 4 (resection cavity). Its displayed 136.3 GB is the record's data volume; the local single ZIP is about 17.0 GB compressed. These totals do not establish eligible patient triplets.

The [author usage guide](https://github.com/deshanyang/GBM-DIR-QA) says sessions are not registered to each other, despite atlas alignment. It says dates were extracted automatically and may be off by several days. Any longitudinal forecast needs a predeclared per-patient registration and interval policy. Automated segmentation is a source label, not ground-truth tumor.

Zenodo lists CC BY-NC-ND 4.0 or a custom Duke license. Download permission alone does not establish permission to share derived masks, predictions, figures, or model weights. Record the exact accepted agreement and school approval before any such use or sharing.

## Prepared inventory

`src/fk_duke_archive_inventory.py` is ready for the completed archive. It:

1. Checks exact file size and the expected Zenodo MD5, then records local SHA-256.
2. Lists ZIP central-directory names and sizes without opening or decompressing members.
3. Can read the header row of one explicitly named, non-landmark CSV. It does not extract archive members to disk.
4. Writes a local inventory under `output/track_b_duke_inventory_20261103/` and refuses to overwrite an existing record.

Synthetic tests for this tool passed in a scratch copy before this note was written. They do not validate the real archive or its metadata. Keep inventories and patient-level metadata local; do not commit them.

After the archive reaches the expected size and checksum passes, run the inventory first. Review filenames and identify metadata tables. Read only approved patient/session/date/treatment/IDH metadata fields. Do not open landmark CSVs, NIfTI images, brain masks, or segmentation masks in this metadata stage.

## Gates after metadata inventory

Before any MRI or mask voxel read:

- Save the written Duke access grant and exact license terms.
- Save written school approval with date, approver, and conditions.
- Reconcile the archive's patient/session map, IDH status, treatments, surgeries, and scan intervals from metadata.
- Check patient-level overlap with cohorts already used in this project, within the permitted terms.
- Count independent patients with eligible three-visit sequences. Freeze the minimum before viewing target outcomes.
- Write and adopt a Duke-specific protocol amendment covering labels, endpoint, support, registration, interval, folds, baselines, exclusions, failures, multiplicity, precision, and runtime.
- Obtain a separate, explicit read scope for MRI-only registration/support QC. Freeze QC rules before applying them.
- Freeze the cohort, code, settings, and hashes before a locked patient evaluation.

If metadata cannot establish a sufficient independent cohort or the rights/approval gates fail, stop the patient-run branch. Do not substitute LUMIERE and describe it as independent confirmation.

## Command after download completes

From the project root, run:

```powershell
py -3.14 -m src.fk_duke_archive_inventory "$env:USERPROFILE\DATASET_FINAL.zip"
```

The command fails closed on a partial file, wrong size, or checksum mismatch. It does not score patients or read image/mask voxels.
