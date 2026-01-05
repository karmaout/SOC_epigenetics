
## -----------------------------
## Subset VGLUT1
## -----------------------------
subset_VGLUT1 <- function(obj, cluster_id = 0, sct_assay = "SCT", gene = "Slc17a7") {
  DefaultAssay(obj) <- sct_assay
  slc_resid <- GetAssayData(obj, assay = sct_assay, slot = "data")[gene, , drop = TRUE]
  clust <- obj@meta.data$seurat_clusters
  flag <- (slc_resid > 0) & (clust %in% c(cluster_id, as.character(cluster_id)))
  obj$VGLUT1_flag <- flag
  subset(obj, cells = names(which(flag)))
}

## -----------------------------
## Subset Fos+
## -----------------------------
subset_Fos <- function(obj, fos_col = "IEG1") {
  stopifnot(fos_col %in% colnames(obj@meta.data))
  subset(obj, subset = VGLUT1_flag & obj@meta.data[[fos_col]] > 0)
}

## -----------------------------
## WNN pipeline (RNA + ATAC)
## -----------------------------
run_WNN <- function(obj,
                    rna_assay = "SCT",
                    atac_assay = "ATAC",
                    npcs = 50,
                    dims_rna = 1:30,
                    dims_atac = 2:30) {
  
  DefaultAssay(obj) <- rna_assay
  if (!"pca" %in% names(obj@reductions)) {
    obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  }
  
  DefaultAssay(obj) <- atac_assay
  if (!"lsi" %in% names(obj@reductions)) {
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj, min.cutoff = "q0")
    obj <- RunSVD(obj)
  }
  
  obj <- FindMultiModalNeighbors(
    obj,
    reduction.list = list("pca", "lsi"),
    dims.list = list(dims_rna, dims_atac)
  )
  
  RunUMAP(
    obj,
    nn.name = "weighted.nn",
    reduction.name = "wnn.umap",
    reduction.key = "wnnUMAP_"
  )
}

## -----------------------------
## FRiP (cell level)
## -----------------------------
calc_FRiP_df <- function(obj,
                         assay = "ATAC",
                         fragments = "fragment_counts",
                         group_col = "behavior_group") {
  obj <- FRiP(obj, assay = assay, total.fragments = fragments)
  obj@meta.data %>%
    filter(!is.na(.data[[group_col]])) %>%
    mutate(
      FRiP  = as.numeric(FRiP),
      FRiP2 = FRiP / 2
    ) %>%
    filter(is.finite(FRiP2))
}

## -----------------------------
## FRiP (pseudo-bulk: library level)
## -----------------------------
calc_pseudobulk_FRiP <- function(df,
                                 lib_col = "library_id",
                                 group_col = "behavior_group") {
  df %>%
    group_by(.data[[lib_col]], .data[[group_col]]) %>%
    summarise(
      n_cells = n(),
      FRiP2_mean   = mean(FRiP2, na.rm = TRUE),
      FRiP2_median = median(FRiP2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(
      library_id     = .data[[lib_col]],
      behavior_group = .data[[group_col]]
    )
}

## -----------------------------
## Excel writer
## -----------------------------
write_xlsx_multi <- function(path, sheets) {
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}



## =========================================================
## 1) Objects to process
## =========================================================
obj_list <- list(
  BLA = BLA,   # must already exist in memory
  PC  = PC
)

params <- list(
  vglut1_cluster = 0,
  vglut1_gene    = "Slc17a7",
  fos_col        = "IEG1",
  rna_assay      = "SCT",
  atac_assay     = "ATAC",
  frag_col       = "fragment_counts",
  group_col      = "behavior_group",
  lib_col        = "library_id"
)

res <- list()

for (region in names(obj_list)) {
  
  message("▶ Processing ", region)
  
  obj <- obj_list[[region]]
  
  ## 1) VGLUT1
  vglut1 <- subset_VGLUT1(
    obj,
    cluster_id = params$vglut1_cluster,
    gene       = params$vglut1_gene
  )
  
  ## 2) Fos+
  vglut1_fos <- subset_Fos(vglut1, fos_col = params$fos_col)
  
  ## 3) WNN
  vglut1     <- run_WNN(vglut1)
  vglut1_fos <- run_WNN(vglut1_fos)
  
  ## 4) Save objects
  saveRDS(vglut1,     file.path(OUT_DIR, paste0(region, "_VGLUT1.rds")))
  saveRDS(vglut1_fos, file.path(OUT_DIR, paste0(region, "_VGLUT1_Fos.rds")))
  
  ## 5) FRiP
  frip     <- calc_FRiP_df(vglut1,     fragments = params$frag_col)
  frip_fos <- calc_FRiP_df(vglut1_fos, fragments = params$frag_col)
  
  frip     <- tibble::rownames_to_column(frip, "cell_id")
  frip_fos <- tibble::rownames_to_column(frip_fos, "cell_id")
  
  ## 6) Pseudobulk
  pb     <- calc_pseudobulk_FRiP(frip)
  pb_fos <- calc_pseudobulk_FRiP(frip_fos)
  
  ## 7) Per-region Excel
  write_xlsx_multi(
    file.path(OUT_DIR, paste0(region, "_FRiP.xlsx")),
    list(
      VGLUT1_cells     = frip,
      VGLUT1_Fos_cells = frip_fos,
      VGLUT1_pb        = pb,
      VGLUT1_Fos_pb    = pb_fos
    )
  )
  
  res[[region]] <- list(
    VGLUT1     = vglut1,
    VGLUT1_Fos = vglut1_fos,
    FRiP      = frip,
    FRiP_Fos  = frip_fos,
    PB        = pb,
    PB_Fos    = pb_fos
  )
}


write_xlsx_multi(
  file.path(OUT_DIR, "ALL_regions_FRiP.xlsx"),
  list(
    BLA_cells     = res$BLA$FRiP,
    BLA_Fos_cells = res$BLA$FRiP_Fos,
    PC_cells      = res$PC$FRiP,
    PC_Fos_cells  = res$PC$FRiP_Fos,
    BLA_pb        = res$BLA$PB,
    PC_pb         = res$PC$PB
  )
)

message("✅ Pipeline finished")

