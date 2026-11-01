# Resources

External references used by the pipeline. These files are not tracked in git; download them into `Data/10x/Resources/` with the layout below.

```
Resources/
├── Common_Variants/          # souporcell common SNPs (GRCh38)
├── Gene_lists/               # HGNC immunoglobulin genes
├── Genomes/                  # 10x Genomics references
├── Reference_10x_Datasets/   # cell-type annotation references
├── Singularity/              # souporcell and Immcantation containers
└── VDJ_References/           # Immcantation repository
```

## Genomes (10x Genomics)

| Folder | Used by | Download |
|---|---|---|
| `refdata-gex-GRCh38-2024-A` | Cell Ranger, souporcell | [10x Genomics reference downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads) |
| `refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0` | Cell Ranger (VDJ-B) | [10x Genomics reference downloads](https://www.10xgenomics.com/support/software/cell-ranger/downloads) |

```bash
cd Genomes
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
curl -O https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz
tar -xzf refdata-gex-GRCh38-2024-A.tar.gz
tar -xzf refdata-cellranger-vdj-GRCh38-alts-ensembl-7.1.0.tar.gz
```

## Containers (Singularity)

| File | Used by | Version (from the image) | Source |
|---|---|---|---|
| `souporcell_v2.0.sif` | `02_Demultiplexing/*/Souporcell.sh` | souporcell 2.0, commit `790b281` (Singularity Hub build, 2 Mar 2021); freebayes 1.3.1, vartrix 1.1.14, minimap2 2.7, bedtools 2.28.0 | Singularity Hub archive: `shub://wheaton5/souporcell:latest` |
| `immcantation_suite_4.7.0.sif` | `03_Downstream/06a__IgBLAST.sh`, `06b__AIRR_DB.sh` | Immcantation 4.7.0 (released 21 Jan 2026; image built 6 Feb 2026) | [Immcantation](https://immcantation.readthedocs.io) Docker image |

```bash
cd Singularity
singularity build immcantation_suite_4.7.0.sif docker://immcantation/suite:4.7.0
```

Key tools inside the Immcantation 4.7.0 image: IgBLAST 1.22.0, Change-O 1.3.4, pRESTO 0.7.7, alakazam 1.4.2, shazam 1.3.1, scoper 1.4.0, dowser 2.4.0.999, tigger 1.1.0. Check with `singularity exec immcantation_suite_4.7.0.sif versions report`.

The souporcell image was built on Singularity Hub on 2 March 2021 and contains souporcell commit [`790b281`](https://github.com/wheaton5/souporcell/commit/790b281dbe4de38cb8bde4e9c4a533e554190624) — 61 commits after the [2.0 release](https://github.com/wheaton5/souporcell/releases/tag/2.0) and before v2.5. Singularity Hub has since been archived; its final build of the same recipe (19 April 2021) was made from the same commit, so it contains identical souporcell code:

```bash
singularity pull souporcell_v2.0.sif shub://wheaton5/souporcell:latest
```

## Common variants

`common_variants_grch38.fixed.vcf`: common SNPs for souporcell genotyping on GRCh38. The original file is linked from the [souporcell README](https://github.com/wheaton5/souporcell).

## Gene lists

`HGNC_Immunoglobulin_genes.txt`: immunoglobulin genes, excluded from variable features in `01__GEX_Integration.Rmd`. Download from the [HGNC gene groups](https://www.genenames.org/data/genegroup/) (immunoglobulin gene groups).

## Annotation references

Used by `02__Annotation.R` (SingleR).

| File | Dataset | Source |
|---|---|---|
| `TrianaS_VeltenL_Healthy_Reference_ALL_dataset.rds` | Triana *et al.*, *Nat Immunol* 2021 — healthy bone marrow proteo-genomic reference | [figshare](https://figshare.com/projects/Single-cell_proteo-genomic_reference_maps_of_the_human_hematopoietic_system/94469) |
| `Stewart2021_BCells.rds` | Stewart *et al.*, *Front Immunol* 2021 — peripheral B-cell subsets | [doi:10.3389/fimmu.2021.602539](https://doi.org/10.3389/fimmu.2021.602539) (see Data Availability) |
| `Terekhova_PBMCs2023/b_cells_*` | Terekhova *et al.*, *Immunity* 2023 — healthy blood atlas, B-cell subset | [Synapse syn49637038](https://www.synapse.org/#!Synapse:syn49637038/) (registration required) |

## VDJ references

`VDJ_References/immcantation/`: clone of the [Immcantation repository](https://bitbucket.org/kleinstein/immcantation) at commit `40b86f8` (4.5.0-70, 23 Jan 2025). Its `fetch_imgtdb.sh` downloads the IMGT germline sequences into `scripts/human/vdj/`, which are used by `06b__AIRR_DB.sh` (MakeDb) and `06c__VDJ_Immcantation_Analysis.Rmd` (germline reconstruction).

```bash
cd VDJ_References
git clone https://bitbucket.org/kleinstein/immcantation.git
cd immcantation && git checkout 40b86f8
cd scripts && ./fetch_imgtdb.sh     # also run automatically by 06b__AIRR_DB.sh if germlines are missing
```

**Germline version used in this analysis:** IMGT/GENE-DB ([imgt.org/genedb](https://www.imgt.org/genedb)), downloaded **6 February 2026** (recorded in `scripts/IMGT.yaml`). `fetch_imgtdb.sh` always downloads the current IMGT release, so a fresh download may differ; `germlines_md5.txt` lists checksums of the FASTA files used.

## Archived copy

All resources used in this analysis were packaged into a single zip archive deposited in [Zenodo](https://doi.org/10.5281/zenodo.xxxxxxx) (DOI: 10.5281/zenodo.xxxxxxx).