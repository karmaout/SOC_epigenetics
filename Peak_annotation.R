## =========================================================
## 0) Packages
## =========================================================
suppressPackageStartupMessages({
  library(Signac)
  library(Seurat)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(data.table)
})

OUT_DIR <- "results/ATAC_peak_annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1) Build TxDb 
## =========================================================
txdb <- makeTxDbFromGFF(GTF_FILE, format = "gtf")

## =========================================================
## 2) Helpers
## =========================================================
get_peaks_gr <- function(obj, assay = "ATAC", keep_std = TRUE) {
  DefaultAssay(obj) <- assay
  gr <- granges(obj)
  peak_ids <- rownames(obj[[assay]])
  mcols(gr)$peak_id <- peak_ids
  if (keep_std) gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr
}

harmonize_seqlevels <- function(gr, txdb) {
  seqlevelsStyle(gr) <- seqlevelsStyle(txdb)
  common <- intersect(seqlevels(gr), seqlevels(txdb))
  keepSeqlevels(gr, common, pruning.mode = "coarse")
}

annotate_peaks <- function(gr, txdb, tss = c(-2000, 2000), annoDb = "org.Rn.eg.db") {
  anno <- annotatePeak(
    peak      = gr,
    TxDb      = txdb,
    tssRegion = tss,
    annoDb    = annoDb
  )
  anno_df <- as.data.frame(anno) %>% relocate(peak_id, .before = 1)
  list(anno = anno, anno_df = anno_df)
}

export_chip_plots <- function(anno, prefix, out_dir = OUT_DIR) {
  pdf(file.path(out_dir, paste0(prefix, "_ChIPseeker_plots.pdf")), width = 7, height = 6)
  print(plotAnnoPie(anno))
  print(plotAnnoBar(anno))
  dev.off()
  
  # TIFF outputs
  tiff_opts <- list(width = 7, height = 6, units = "in", res = 300, compression = "lzw")
  
  do.call(tiff, c(list(filename = file.path(out_dir, paste0(prefix, "_anno_pie.tiff"))), tiff_opts))
  plotAnnoPie(anno)
  dev.off()
  
  do.call(tiff, c(list(filename = file.path(out_dir, paste0(prefix, "_anno_bar.tiff"))), tiff_opts))
  plotAnnoBar(anno)
  dev.off()
}

compute_peak_metrics_by_group <- function(obj,
                                          assay = "ATAC",
                                          group_col = "behavior_group",
                                          group_levels = c("U_A","P_A","U_R","P_R","U_E","P_E")) {
  DefaultAssay(obj) <- assay
  obj[[group_col]] <- factor(obj[[group_col]][,1], levels = group_levels)
  
  counts <- GetAssayData(obj, assay = assay, slot = "counts")  # dgCMatrix
  meta   <- obj@meta.data
  levs   <- levels(droplevels(obj[[group_col]][,1]))
  levs   <- levs[!is.na(levs)]
  
  grp_colsums <- function(mat) {
    lapply(levs, function(g) {
      cells_g <- rownames(meta)[meta[[group_col]] == g]
      if (length(cells_g) == 0) return(Matrix::rowSums(mat[, 0, drop = FALSE]))
      Matrix::rowSums(mat[, cells_g, drop = FALSE])
    }) |> stats::setNames(levs)
  }
  
  pseudo_counts <- grp_colsums(counts)
  lib_sizes     <- vapply(pseudo_counts, sum, numeric(1))
  
  cpm_list <- lapply(levs, function(g) {
    s <- pseudo_counts[[g]]
    lib <- lib_sizes[[g]]
    if (lib == 0) rep(0, length(s)) else (s / lib) * 1e6
  })
  names(cpm_list) <- levs
  
  frac_list <- lapply(levs, function(g) {
    cells_g <- rownames(meta)[meta[[group_col]] == g]
    if (length(cells_g) == 0) return(rep(NA_real_, nrow(counts)))
    mat_g <- counts[, cells_g, drop = FALSE]
    Matrix::rowSums(mat_g > 0) / ncol(mat_g)
  })
  names(frac_list) <- levs
  
  metrics_df <- do.call(cbind, c(
    setNames(cpm_list,  paste0("cpm_",  levs)),
    setNames(frac_list, paste0("frac_", levs))
  )) |> as.data.frame()
  
  metrics_df$peak_id <- rownames(obj[[assay]])
  metrics_df
}

add_log2fc <- function(df, a, b) {
  ca <- paste0("cpm_", a); cb <- paste0("cpm_", b)
  if (all(c(ca, cb) %in% colnames(df))) {
    df[[paste0("log2FC_", b, "_vs_", a)]] <- log2((df[[cb]] + 1) / (df[[ca]] + 1))
  }
  df
}

deg_filter_and_label <- function(df,
                                 deg_common,
                                 deg_pc,
                                 deg_bla,
                                 symbol_col = "SYMBOL") {
  deg_all <- unique(c(deg_common, deg_pc, deg_bla))
  df %>%
    mutate(
      SYMBOL = str_trim(as.character(.data[[symbol_col]])),
      DEG_category = case_when(
        SYMBOL %in% deg_common ~ "PC&BLA common",
        SYMBOL %in% deg_pc     ~ "PC-specific",
        SYMBOL %in% deg_bla    ~ "BLA-specific",
        TRUE                   ~ "Other"
      )
    ) %>%
    filter(SYMBOL %in% deg_all)
}

plot_deg_boxplots <- function(hits_df, prefix, out_dir = OUT_DIR) {
  cat_levels <- c("PC&BLA common","PC-specific","BLA-specific")
  
  df_long <- hits_df %>%
    filter(DEG_category %in% cat_levels) %>%
    mutate(DEG_category = factor(DEG_category, levels = cat_levels)) %>%
    pivot_longer(
      cols = c(log2FC_P_A_vs_U_A, log2FC_P_R_vs_U_R),
      names_to = "comparison",
      values_to = "log2FC"
    ) %>%
    mutate(
      comparison = recode(
        comparison,
        "log2FC_P_A_vs_U_A" = "Encoding: P_A vs U_A",
        "log2FC_P_R_vs_U_R" = "Retrieval: P_R vs U_R"
      )
    ) %>%
    drop_na(log2FC)
  
  p_grid <- ggplot(df_long, aes(x = 1, y = log2FC, fill = DEG_category)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.08, alpha = 0.4, size = 0.9) +
    facet_grid(rows = vars(DEG_category), cols = vars(comparison)) +
    labs(x = NULL, y = "log2 Fold Change") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(face = "bold"))
  
  df_long <- df_long %>% mutate(xgroup = interaction(DEG_category, comparison, sep = " • ", drop = TRUE))
  p_single <- ggplot(df_long, aes(x = xgroup, y = log2FC, fill = DEG_category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 0.9) +
    labs(x = NULL, y = "log2 Fold Change") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  
  ggsave(file.path(out_dir, paste0(prefix, "_DEG_peaks_boxplot_grid.pdf")), p_grid, width = 8, height = 9)
  ggsave(file.path(out_dir, paste0(prefix, "_DEG_peaks_boxplot_single.pdf")), p_single, width = 11, height = 6)
}

export_fragments_by_group <- function(obj,
                                      out_dir = file.path(OUT_DIR, "fragments_by_group"),
                                      assay = "ATAC",
                                      group_col = "behavior_group",
                                      group_levels = c("U_A","P_A","U_R","P_R")) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  DefaultAssay(obj) <- assay
  obj[[group_col]] <- factor(obj[[group_col]][,1], levels = group_levels)
  levs <- levels(droplevels(obj[[group_col]][,1]))
  
  frags_path <- Fragments(obj, assay = assay)[[1]]@path
  message("Using fragment file: ", frags_path)
  
  frag_dt <- fread(frags_path, header = FALSE, showProgress = FALSE)
  if (ncol(frag_dt) == 4) {
    setnames(frag_dt, c("chr","start","end","barcode"))
    frag_dt[, count := 1L]
  } else if (ncol(frag_dt) >= 5) {
    setnames(frag_dt, c("chr","start","end","barcode","count")[seq_len(5)])
  } else {
    stop("Unexpected fragments file format: ", ncol(frag_dt), " columns.")
  }
  
  for (g in levs) {
    cells_g <- colnames(obj)[obj@meta.data[[group_col]] == g]
    if (length(cells_g) == 0) next
    
    outfile <- file.path(out_dir, paste0("fragments_", g, ".tsv.gz"))
    message("Writing: ", outfile)
    
    dt_sub <- frag_dt[barcode %chin% cells_g]
    if (nrow(dt_sub) == 0L) next
    
    fwrite(dt_sub, file = outfile, sep = "\t", quote = FALSE, compress = "gzip")
  }
  
  invisible(TRUE)
}

batch_coverage_plots <- function(obj, regions_named,
                                 prefix,
                                 out_dir = file.path(OUT_DIR, "coverageplots"),
                                 group_col = "behavior_group",
                                 idents = c("U_A","P_A","U_R","P_R"),
                                 extend_up = 2000,
                                 extend_down = 200) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (nm in names(regions_named)) {
    region <- regions_named[[nm]]
    p <- CoveragePlot(
      object = obj,
      region = region,
      group.by = group_col,
      idents = idents,
      extend.upstream = extend_up,
      extend.downstream = extend_down,
      peaks = TRUE,
      annotation = FALSE
    )
    ggsave(file.path(out_dir, paste0(prefix, "_", nm, ".pdf")), p, width = 10, height = 7)
  }
}
