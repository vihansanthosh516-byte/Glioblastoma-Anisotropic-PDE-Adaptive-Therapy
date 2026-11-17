# Track A Artzi Zenodo metadata audit

Dylan Harwood directed the user to [Zenodo 22791462](https://zenodo.org/records/22791462). The record lists Final_meta.csv and a separate expression matrix. The audit read the entire metadata CSV, verified MD5 53057e9514bda0a070edc3f3202bb35b against the [record API](https://zenodo.org/api/records/22791462), and did not open expression.

Counting distinct nonblank Patient codes by Setting gives 28 P-CT and 12 P-IT patients. The intersection is 10 patient codes. Restricting both to Malignant = Malignant retains all 10. These counts come from Final_meta.csv in the linked release. No recurrent core or surgical margin records were included in this overlap.

This is a candidate paired denominator. The CSV supplies Patient, Setting, coreID, fov and Malignant but no explicit resection identifier. Verify from methods or author documentation that both primary setting codes refer to the same initial surgery and that the settings were assigned independently of RNA. Confirm IDH-wildtype status and donor independence from development datasets before assigning the test role.

The metadata field Subtype is NA in P-IT, including its malignant cells. P-CT Subtype values do not supply a shared region-independent malignant-state annotation. Do not apply the state-matched candidate based on these labels alone. Verify whether label encodes a defensible shared cell state or freeze assignment using a development-only reference and check coverage. The alternative is a separately prespecified malignant pseudobulk benchmark. No model was scored.
