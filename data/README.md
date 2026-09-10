Raw public data are downloaded automatically by:

```text
scripts/01_download_GSE145410.R
```

The workflow uses GEO accession **GSE145410**.

The downloaded files are deliberately excluded from GitHub by `.gitignore`.

Expected local structure after download:

```text
data/
└── GSE145410/
    └── extracted/
        ├── GSM4317809_.../
        ├── GSM4317810_.../
        └── ...
```

Do not place CRCM patient-level data in this repository.
