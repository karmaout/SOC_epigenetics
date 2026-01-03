suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(openxlsx)
})

## =========================================================
## 0) Config (repo-relative)
## =========================================================
CFG <- list(
  in_dir  = "data/motifs",          # where hit_*.txt lives (edit as needed)
  anno_dir = "results/ATAC_peak_annotation/02_metrics",  # where annotation csv lives
  out_dir = "results/motif_hits",
  out_xlsx = "motif_hits_BLA_PC_annotated.xlsx",
  
  # HOMER hit file is typically BED-like (0-based start). Convert start+1 to match your peak_id.
  homer_start_is_0based = TRUE,
  
  # If your annotation peak_id format differs, adjust in build_peak_id()
  peak_id_sep = "-"
)

dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

## Inputs (edit file names only; directories above handle paths)
hit_files <- list(
  BLA = file.path(CFG$in_dir, "BLA_peak_motif_hits_v2.txt"),
  PC  = file.path(CFG$in_dir, "PC_peak_motif_hits_v2.txt")
)

anno_files <- list(
  BLA = file.path(CFG$anno_dir, "BLA_VGLUT1_ATAC_peaks_ChIPseeker_annotation_wGroups.csv"),
  PC  = file.path(CFG$anno_dir, "PC_VGLUT1_ATAC_peaks_ChIPseeker_annotation_wGroups.csv")
)

## =========================================================
## 1) Helpers
## =========================================================
assert_exists <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  invisible(TRUE)
}

build_peak_id <- function(chr, peak_start, peak_end, start_is_0based = TRUE, sep = "-") {
  chr <- as.character(chr)
  start_1based <- if (start_is_0based) peak_start + 1 else peak_start
  paste(chr, start_1based, peak_end, sep = sep)
}

read_homer_hits <- function(hit_file) {
  # HOMER hit file: 13 columns, no header (your original assumption)
  hits <- readr::read_tsv(hit_file, col_names = FALSE, show_col_types = FALSE)
  
  if (ncol(hits) < 13) {
    stop("Hit file has < 13 columns: ", hit_file, " (ncol=", ncol(hits), ")", call. = FALSE)
  }
  
  hits <- hits[, 1:13]
  colnames(hits) <- c(
    "motif_rank", "motif_pos", "sequence", "motif_name",
    "strand", "score", "chr", "peak_start", "peak_end",
    "peak_center", "hit_start", "hit_end", "hit_center"
  )
  
  hits
}

read_peak_annotation <- function(anno_file) {
  anno <- readr::read_csv(anno_file, show_col_types = FALSE)
  
  needed <- c("peak_id", "SYMBOL", "GENENAME", "annotation", "distanceToTSS")
  missing <- setdiff(needed, colnames(anno))
  if (length(missing) > 0) {
    stop("Annotation file missing columns: ", paste(missing, collapse = ", "),
         "\nFile: ", anno_file, call. = FALSE)
  }
  
  anno %>%
    select(all_of(needed)) %>%
    mutate(peak_id = as.character(peak_id)) %>%
    distinct(peak_id, .keep_all = TRUE)   # prevent accidental 1-to-many join inflation
}

annotate_hits_with_genes <- function(hit_file, anno_file, cfg = CFG) {
  assert_exists(hit_file)
  assert_exists(anno_file)
  
  hits <- read_homer_hits(hit_file) %>%
    mutate(
      peak_id = build_peak_id(chr, peak_start, peak_end,
                              start_is_0based = cfg$homer_start_is_0based,
                              sep = cfg$peak_id_sep)
    )
  
  anno <- read_peak_annotation(anno_file)
  
  out <- hits %>% left_join(anno, by = "peak_id")
  
  na_rate <- mean(is.na(out$SYMBOL)) * 100
  message(
    basename(hit_file),
    " -> matched SYMBOL: ",
    sum(!is.na(out$SYMBOL)), "/", nrow(out),
    " (NA: ", round(na_rate, 2), "%)"
  )
  
  out
}

write_multi_sheet_xlsx <- function(path, sheets_named_list) {
  wb <- createWorkbook()
  for (nm in names(sheets_named_list)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets_named_list[[nm]])
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}

## =========================================================
## 2) Run (BLA + PC)
## =========================================================
hits_BLA_annot <- annotate_hits_with_genes(hit_files$BLA, anno_files$BLA)
hits_PC_annot  <- annotate_hits_with_genes(hit_files$PC,  anno_files$PC)

out_xlsx_path <- file.path(CFG$out_dir, CFG$out_xlsx)

write_multi_sheet_xlsx(
  out_xlsx_path,
  list(
    BLA_motif_hits = hits_BLA_annot,
    PC_motif_hits  = hits_PC_annot
  )
)

message("✅ Wrote Excel: ", out_xlsx_path)
