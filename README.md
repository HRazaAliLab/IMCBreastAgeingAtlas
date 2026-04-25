# Single-cell spatial atlas of the ageing human breast

Code to reproduce analyses and figures from:

Gupta P*, Lee E*, Masqué Soler N, et al. **Single-cell spatial atlas of the aging human breast.** _Nat Aging_. 2026;6(4):916-931. doi:[10.1038/s43587-026-01104-3](https://www.nature.com/articles/s43587-026-01104-3)

Using imaging mass cytometry (IMC), this study profiles >3 million cells from 527 normal human breast samples to characterize age-associated changes in cellular composition, tissue architecture, and spatial interactions.

---

## Overview

This repository contains the analysis and figure-generation code associated with the manuscript, together with documentation describing the released datasets.

- The **Zenodo record** provides the archival snapshot used for publication.
- This GitHub repository hosts the actively maintained codebase.

Detailed documentation:
- **DataGuide.pdf** – released datasets, formats, column definitions, and usage notes  
- **CodeGuide.pdf** – codebase structure and figure-to-script mapping

---

## Repository structure

```text
code/
├── header.R                 shared project header
├── generalUtilities.R       helper and plotting utilities
├── main/                    main figure code
├── ext/                     extended data figure code
├── supp/                    supplementary figure code
└── condaEnv.yml             reference software environment
```

Processed and raw data are described in **DataGuide.pdf**.

---

## Quick start

### Requirements
- Linux or macOS  
- R ≥ 4.2  

Required R packages can be installed from CRAN or Bioconductor as needed.  
The provided `condaEnv.yml` documents the reference environment used during development but is not a strict requirement.

Some scripts use fork-based parallelisation (`mclapply`) and may require modification on Windows.

### Load the project header

```r
source(here::here("code/header.R"))
```

This loads required packages, helper functions, and shared settings.

### Run a figure script

For example:

```r
source("code/main/Fig1/mkClusterHeatmaps.R")
```

Outputs are written to a figure-specific subdirectory under `scratch/`.

A complete figure-to-script mapping is provided in **CodeGuide.pdf**.

---

## Data availability

All data required to reproduce the figures in the manuscript are available via the associated Zenodo record.

See **DataGuide.pdf** for:
- Dataset descriptions  
- File formats  
- Column definitions  
- Intended use and limitations
