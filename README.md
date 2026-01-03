# SOC Multiome Analysis

This repository contains code and shell scripts to reproduce the multi‑omic analysis of single‑nucleus RNA‑sequencing (snRNA‑seq) and ATAC‑sequencing (snATAC‑seq) data from the basolateral amygdala (BLA) and piriform cortex (PC) during a second‑order conditioning (SOC) paradigm.  The raw data are not included; instead, this repository focuses on processing and analysis pipelines designed to run on 10x Genomics Multiome data.

## Data overview

The experiments profiled two brain regions (BLA and PC) under several behavioural conditions: encoding (P_A and U_A) and retrieval (P_R and U_R).  Each library corresponds to one mouse and one brain region.  Counts are generated using Cell Ranger ARC and stored as Seurat objects for downstream analysis.

### Behaviour groups

| Abbreviation | Description            |
|--------------|------------------------|
| **P_A**      | Paired conditioning, acquisation phase |
| **U_A**      | Unpaired conditioning, acquisation phase |
| **P_R**      | Paired conditioning, retrieval phase |
| **U_R**      | Unpaired conditioning, retrieval phase |

These labels appear in metadata columns of the Seurat objects and in file and directory names throughout the repository.

## Repository structure


Throughout the analysis the code writes intermediate and final files into a `results/` directory created relative to the project root.  Subdirectories such as `results/VGLUT1_pipeline`, `results/ATAC_peak_annotation`, `results/RNA_CPM` and `results/exports` will be created automatically.

## R analysis pipelines

All R scripts are written to be sourced from an R session; they assume that the required data objects (e.g. `BLA` and `PC` Seurat objects) are already loaded in memory and that relative paths point to the appropriate input files.  The scripts make extensive use of the [`Seurat`](https://satijalab.org/seurat/) and [`Signac`](https://github.com/timoast/signac) frameworks for multi‑modal analysis.

### `SOC_VGLUT1_ATAC.R`

This script orchestrates a multi‑step pipeline for each brain region:

1. **Subset VGLUT1 neurons:** The function `subset_VGLUT1` selects cells with positive residual expression of the Slc17a7 gene (encoding VGLUT1) in a specified cluster:contentReference[oaicite:11]{index=11}.  It adds a logical `VGLUT1_flag` column and subsets the Seurat object accordingly.
2. **Identify Fos⁺ cells:** `subset_Fos` further subsets the VGLUT1 cells to those expressing an immediate early gene (default `IEG1`):contentReference[oaicite:12]{index=12}.
3. **Weighted Nearest Neighbour (WNN) integration:** `run_WNN` runs PCA on the RNA assay and LSI on the ATAC assay, then integrates modalities using weighted nearest neighbours and computes a UMAP embedding:contentReference[oaicite:13]{index=13}.
4. **Compute FRiP metrics:** `calc_FRiP_df` calculates the fraction of reads in peaks (FRiP) for each cell and also computes a transformed `FRiP2` metric:contentReference[oaicite:14]{index=14}.  `calc_pseudobulk_FRiP` then aggregates FRiP2 values by library and behaviour group to produce pseudo‑bulk summaries:contentReference[oaicite:15]{index=15}.
5. **Export results:** For each region, the script saves RDS objects, writes per‑cell and pseudo‑bulk FRiP data to Excel sheets, and finally collects results across regions into a combined workbook:contentReference[oaicite:16]{index=16}.

To run the pipeline for BLA and PC, load your Seurat objects (`BLA` and `PC`) into your R environment, ensure the required metadata columns (`seurat_clusters`, `behavior_group`, `IEG1`, `library_id`, `fragment_counts`) are present, and then:

```r
source("SOC_VGLUT1_ATAC.R")


