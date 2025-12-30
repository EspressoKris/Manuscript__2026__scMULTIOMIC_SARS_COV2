# Project information:

**Coded in/requires**:

![Static Badge](https://img.shields.io/badge/v4.3.3-blue?style=flat&logo=R)

## Overview
This repository contains the analysis code and relevant documentation for the single-cell multi-omics profiling of SARS-CoV-2 infection using scBEAM-BCR technology.

## Contact
**First Author:** Carl Graham (carl.1.graham@kcl.ac.uk)   
**Senior Author:** Katie J Doores (katie.doores@kcl.ac.uk) 

## General Description
Single-cell characterization of immune responses to SARS-CoV-2 infection with synchroneous profiling of:
- Gene expression (scRNA-seq)
- Target surface protein capture (CITE-seq)
- Target antigen capture (BEAM-seq)
- B-cell receptor repertoires (BCR-seq)

## Graphical Abstract 


## Citation 

> Graham, C.

--- 

## Data Availability

### Raw Data
- **GEO Accession**: [GSExxxxxx](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSExxxxxx)
- **SRA**: [SRPxxxxxx](https://www.ncbi.nlm.nih.gov/sra?term=SRPxxxxxx)

### Processed Data
- **Zenodo**: [DOI:10.5281/zenodo.xxxxxxx](https://doi.org/10.5281/zenodo.xxxxxxx)
- Seurat objects
- Count matrices
- Metadata files

## Directory Structure

```
.
├── README.md
├── data/
│   ├── raw/                    # Raw sequencing data (not tracked)
│   ├── processed/              # Processed count matrices
│   └── metadata/               # Sample and clinical metadata
├── scripts/
│   ├── 01_preprocessing/       # Quality control and filtering
│   ├── 02_integration/         # Data integration and batch correction
│   ├── 03_clustering/          # Cell type identification
│   ├── 04_differential/        # Differential expression analysis
│   ├── 05_bcr_analysis/        # BCR repertoire analysis
│   └── 06_visualization/       # Figure generation
├── results/
│   ├── figures/                # Publication-ready figures
│   ├── tables/                 # Supplementary tables
│   └── reports/                # Analysis reports
├── notebooks/                  # R Markdown/Jupyter notebooks
└── environment/                # Conda/Docker specifications
```

--- 

## Dependencies

### R Packages
- Seurat (v5.0+)

## Key Analyses

1. **Quality Control**: `scripts/01_preprocessing/`
2. **Data Integration**: `scripts/02_integration/`
3. **Cell Type Annotation**: `scripts/03_clustering/`
4. **Differential Expression**: `scripts/04_differential/`
5. **BCR Repertoire Analysis**: `scripts/05_bcr_analysis/`
6. **Figure Generation**: `scripts/06_visualization/`

--- 

## License
This project is licensed under the CC4BY License - see LICENSE file for details.
```