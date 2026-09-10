# Data

No patient-level or confidential research data are stored in this
repository.

For testing the workflow, use a publicly available single-cell
RNA-seq dataset.

Expected 10X structure:

```text
data/
└── raw_matrix/
    ├── barcodes.tsv.gz
    ├── features.tsv.gz
    └── matrix.mtx.gz
```

For multi-sample integration:

```text

data/
└── samples/
    ├── sample1/
    │   ├── barcodes.tsv.gz
    │   ├── features.tsv.gz
    │   └── matrix.mtx.gz
    │
    └── sample2/
        ├── barcodes.tsv.gz
        ├── features.tsv.gz
        └── matrix.mtx.gz
```

The datasets used for demonstration should be publicly available and
properly cited in the repository documentation.


---

# 14. `results/README.md`

```markdown
# Results

This directory contains generated analysis outputs.

Examples include:

- QC summaries
- Highly variable genes
- Seurat objects
- Cluster markers
- Differential expression tables

Large generated files and raw datasets should not be committed to
GitHub.

The repository is intended to contain selected lightweight results
and figures demonstrating the analysis workflow.
