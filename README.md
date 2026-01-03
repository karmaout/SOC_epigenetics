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

## Analysis pipelines

All R scripts are written to be sourced from an R session; they assume that the required data objects (e.g. `BLA` and `PC` Seurat objects) are already loaded in memory and that relative paths point to the appropriate input files.  The scripts make extensive use of the [`Seurat`](https://satijalab.org/seurat/) and [`Signac`](https://github.com/timoast/signac) frameworks for multi‑modal analysis. The pipeline combines cell-type–specific subsetting, multi-modal integration, chromatin accessibility quantification, transcription factor motif analysis, and RNA–ATAC correlation to link neuronal activity, chromatin state, and gene expression.



## Analysis schematic
```
10x Multiome data (RNA + ATAC)
│
├── Brain regions
│ ├── BLA
│ └── PC
│
├── Cell selection
│ ├── Cluster-based filtering
│ ├── VGLUT1⁺ neurons (Slc17a7)
│ └── Fos⁺ vs Fos⁻ (IEG1)
│
├── Multi-modal integration
│ ├── RNA: PCA
│ ├── ATAC: LSI
│ └── WNN integration → joint UMAP
│
├── Chromatin accessibility (ATAC)
│ ├── Peak matrix construction
│ ├── FRiP calculation
│ │ ├── Cell-level FRiP
│ │ └── Pseudo-bulk FRiP (library × condition)
│ ├── Motif discovery (HOMER)
│ ├── Motif–peak–gene annotation
│ └── Gene-level ATAC CPM
│
├── Gene expression (RNA)
│ ├── CPM calculation
│ └── Gene-level pseudo-bulk (library × condition)
│
├── RNA–ATAC integration
│ ├── Merge RNA & ATAC gene-level CPM
│ ├── Correlation analysis
│ │ ├── Spearman
│ │ └── Pearson
│ └── Region- and condition-specific comparisons
│
└── Outputs
├── RDS objects (VGLUT1 / Fos⁺)
├── Excel workbooks
│ ├── FRiP-based accessibility metrics
│ ├── Motif annotations
│ └── RNA–ATAC correlation tables
└── Figures and tables for main results
```


**Note:** FRiP is used here as a *primary measure of chromatin accessibility* and contributes directly to main figures, rather than serving solely as a quality control metric.

---

### Pipeline steps

1. **Subset VGLUT1 neurons**  
   The function `subset_VGLUT1` identifies excitatory neurons based on positive residual expression of *Slc17a7* (VGLUT1) within a specified cluster.  
   A logical column (`VGLUT1_flag`) is added to the metadata, and the Seurat object is subset accordingly.

2. **Identify Fos⁺ cells**  
   The function `subset_Fos` further subsets VGLUT1 neurons based on expression of an immediate early gene (default: `IEG1`), enabling separation of Fos⁺ and Fos⁻ populations for activity-dependent analyses.

3. **Weighted Nearest Neighbour (WNN) integration**  
   The function `run_WNN` performs:
   - PCA on the RNA assay  
   - LSI on the ATAC assay  
   - Multi-modal integration using Seurat’s weighted nearest neighbour (WNN) framework  
   - UMAP embedding for joint RNA–ATAC visualization

4. **Motif discovery and annotation (ATAC)**  
   Chromatin accessibility peaks from VGLUT1 and VGLUT1 Fos⁺ neurons are used for transcription factor motif analysis.  
   Motif discovery is performed externally using **HOMER**, and motif hit files are integrated back into the analysis to:
   - Identify transcription factor binding motifs enriched in accessible chromatin  
   - Annotate peaks with motif presence and associated gene symbols  
   - Enable stratification of peaks and genes based on motif occupancy (e.g. AP-1 / Fos-family motifs)

5. **RNA–ATAC correlation analysis (pseudo-bulk)**  
   To assess coupling between transcriptional output and chromatin accessibility:
   - RNA expression is summarized as pseudo-bulk CPM per library and behaviour group  
   - ATAC accessibility is aggregated to gene-level CPM based on annotated peaks  
   - Spearman and Pearson correlations are computed across regions and behavioural conditions  

   This analysis identifies genes showing coordinated or divergent regulation between chromatin accessibility and transcription.

6. **FRiP-based chromatin accessibility analysis**  
   - `calc_FRiP_df` computes FRiP and a transformed metric (`FRiP2`) at the single-cell level  
   - `calc_pseudobulk_FRiP` aggregates FRiP2 at the library × condition level  

   FRiP is treated as a biologically interpretable measure of chromatin accessibility and is used directly in downstream analyses and figures.

7. **Export results**  
   For each brain region, the pipeline:
   - Saves processed Seurat objects (VGLUT1 and VGLUT1 Fos⁺) as RDS files  
   - Exports FRiP metrics, motif annotations, and RNA–ATAC correlation results to Excel workbooks  
   - Merges results across regions into combined output files for downstream statistics and figure generation
