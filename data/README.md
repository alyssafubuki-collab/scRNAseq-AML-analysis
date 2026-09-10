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

