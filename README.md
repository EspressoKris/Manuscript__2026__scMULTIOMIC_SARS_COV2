# Beyond polyclonal antibody responses: dissecting susceptibility to newly emerging SARS-CoV-2 variants

![R](https://img.shields.io/badge/R-4.4.1-blue?style=flat&logo=R)

Analysis code for the single-cell multi-omic profiling of memory B-cell responses to SARS-CoV-2 vaccination and omicron breakthrough infection.

## Summary

This study investigated six donors sampled after three doses of the BNT162b2 mRNA vaccine and again following an omicron lineage breakthrough infection (BTI). Droplet-based scRNA-seq, BCR-seq and CITE-seq were combined to resolve the gene expression, variable region sequence and antigen binding of individual B-cells. Molecular barcoded antigen probes identified Spike specific B-cells, the sub-domain they targeted and their specificity for the ancestral or BTI strains. Gene expression classified the B-cells into naive (transitional and naive) and memory (activated, resting, IgM and atypical) subsets. Variable region sequence analysis yielded complementarity-determining region lengths and isotype, and inferred germline gene usage, somatic hypermutation and clonal relationships. Monoclonal antibodies representing the MBC compartment were expressed recombinantly and characterised in binding, neutralisation and effector function assays. Together, this work showed how repeated single-strain vaccination can shape MBC responses before and after infection with an antigenically drifted variant.

## Contributors

- **Carl Graham** (first author) · carl.1.graham@kcl.ac.uk
- **Kristian Gurashi** (computational analysis) · kristian.gurashi@ndcls.ox.ac.uk
- **Katie J. Doores** (senior author) · katie.doores@kcl.ac.uk

## Citation

> Graham, C., ... Doores, K. J. *Beyond polyclonal antibody responses: dissecting susceptibility to newly emerging SARS-CoV-2 variants.* (in preparation)

## Graphical Abstract

<!-- add figure -->

---

## Study design

Six donors were sampled twice: after three vaccine doses (**T5**) and after omicron breakthrough infection (**T8**). From each sample, general B cells (**GB**) and antigen-specific B cells (**ASB**) were sorted, labelled with hashtags, and pooled three donors per 10x Genomics run. For each cell we measured gene expression, surface proteins, antigen binding and the antibody (BCR) sequence.

| Run | Time point | Donors |
|---|---|---|
| L364A | T5 | VA013, VA022, VA035 |
| L364B | T8 | VA013, VA022, VA035 |
| L364C | T8 | VA025, VA038, VA059 |
| L364D | T5 | VA025, VA038, VA059 |

## Analyses

1. **Read the sequencing data** — counts per cell and antibody sequences (Cell Ranger).
2. **Assign each cell to its donor** — using hashtags and genetic differences (souporcell).
3. **Identify cell types** — quality control, integration and annotation (Seurat).
4. **Analyse antibody sequences** — genes, clones and mutations (Immcantation).
5. **Link antibodies to the virus** — which Spike region and which strain each B cell binds, plus figures.

## Repository

```
01_CellRanger/      # 1. sequencing data
02_Demultiplexing/  # 2. donor assignment
03_Downstream/      # 3–5. cell types, antibodies, figures (scripts numbered in run order)
```

## Data

- **GEO:** [GSE352908](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE352908)
  - Raw: one set of FASTQs per 10x lane and sequencing library (libraries are patient multiplexed).
  - Processed: one HDF5 file per donor, time point and cell population, plus antibody sequences per run.
  - Resources: resources used for data analysis are downloadable at [Zenodo](https://zenodo.org/records/23914018/files/Resources.tar.gz?download=1&preview=1)

## Requirements

R 4.4.1. Exact package versions are saved next to each script (`__Session_Info.txt`).

Running the pipeline also requires external resources (reference genomes, containers, annotation datasets and germline references). See [`Data/10x/Resources/README.md`](Data/10x/Resources/README.md) for what is needed and where to download it.

## License

CC BY 4.0 — see `LICENSE`.