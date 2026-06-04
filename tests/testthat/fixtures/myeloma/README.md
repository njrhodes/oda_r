# Myeloma fixture provenance

This fixture is included for oda CTA parity and regression testing.

Source:
- The underlying myeloma data are public data from the survminer R package.
- survminer documents `myeloma` as extracted from publicly available gene expression data, GEO ID GSE4581.
- The survminer data frame has 256 rows and 12 columns, including survival status, time, treatment, molecular group, chr1q21 status, and gene-expression columns.

Use in oda:
- The fixture is transformed into CTA/MegaODA-style input files for parity testing.
- It is used to validate weighted CTA behavior, MINDENOM behavior, LOO STABLE behavior, missing-code handling, endpoint counts, and no-tree behavior.
- Reference CTA.exe outputs are included as golden anchors for MINDENOM = 1, 30, and 56.

Privacy:
- This is not private institutional data.
- No PHI is included.
- No local path or user-specific data dependency is required.

References:
- survminer myeloma dataset documentation: https://rpkgs.datanovia.com/survminer/reference/myeloma.html
- GEO accession: GSE4581
