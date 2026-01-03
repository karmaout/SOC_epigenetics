## =========================================================
## run_BLA_PC_VGLUT1_Correlation_pipeline.R
##
## One-script integrated pipeline:
##  1) Gene sets (DEG / GO)
##  2) RNA pseudo-bulk CPM (BLA vs PC) by group
##  3) ATAC gene-level CPM from peak tables (BLA vs PC)
##  4) Merge RNA+ATAC, export Excel workbooks
##  5) Correlation summaries across Excel sheets
##
## Inputs expected in memory:
##   - BLA_VGLUT1 (Seurat object)
##   - PC_VGLUT1  (Seurat object)
##
## Inputs expected on disk:
##   - results/ATAC_peak_annotation/02_metrics/BLA_ATAC_peaks_annotated_with_metrics.csv
##   - results/ATAC_peak_annotation/02_metrics/PC_ATAC_peaks_annotated_with_metrics.csv
##
## Outputs:
##   - results/RNA_CPM/...
##   - results/ATAC_gene_CPM/...
##   - results/exports/...
## =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(readr)
  library(openxlsx)
})

## =========================================================
## 0) Global config (repo-relative)
## =========================================================
CFG <- list(
  # behavior groups to compute for RNA (edit freely)
  rna_groups = c("P_R", "U_R", "P_A"),
  
  # peak CSV inputs (from your ATAC peak annotation pipeline)
  bla_peak_csv = "results/ATAC_peak_annotation/02_metrics/BLA_ATAC_peaks_annotated_with_metrics.csv",
  pc_peak_csv  = "results/ATAC_peak_annotation/02_metrics/PC_ATAC_peaks_annotated_with_metrics.csv",
  
  # output folders
  out_rna  = "results/RNA_CPM",
  out_atac = "results/ATAC_gene_CPM",
  out_exp  = "results/exports"
)

dir.create(CFG$out_rna,  showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_atac, showWarnings = FALSE, recursive = TRUE)
dir.create(CFG$out_exp,  showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1) Gene sets (from your script) — keep centralized here
## =========================================================
## ---- PR gene sets ----
pr_deg_common <- c("Arc","Bdnf","Dusp6","Egr1","Egr4","Elmo1","Fos","Fosb",
                   "Fosl2","Lingo1","Nr4a1","Ntrk2","Pde10a","Pde7b","Scg2",
                   "Tiparp","Trib1")

pr_deg_bla_specific <- c(
  "Abhd15","Ackr1","Acot3","Arl9","Asb12","Atp13a4","Bdkrb2","Bmp2","Cd2","Ddr2",
  "Fam178b","Gpa33","Gprc5a","Htra4","Irgm","Kcnk3","Kcnmb1","Klhl30","Lxn","Mapk4",
  "Maml3","NEWGENE-1310011","Ncam1","Nptx2","Plcxd2","Poted","Ptx3","Rab5al1",
  "RGD1563072","RT1-DOb","Slc35e1","Snap25","Sorbs2","Stc2","Tmc8","Ttc29","Zfp773-ps1","Zp3r"
)

pr_deg_pc_specific <- c(
  "Egr2","Nr4a3","Nr4a2","Vgf","Gadd45g","Gadd45b","Hmga1","Homer1","Nfil3","Alkal2",
  "Mas1","R3hdm1","Idi1","Irs2","Abra","Arl4d","Sik1","Spry2","Cdkn1a","Pim3",
  "Lonrf3","Ccn1","Tfap2c","Ky","Gast","Hhipl2","Atf3","Gnrhr","Gch1","Prl4a1",
  "Rab20","Msh5","Osbpl8","Prx","Actg2","Cd200r1","Nap1l5","Slc6a17","Sorcs1","Amer3",
  "Rapgef5","Numb","Smad7","Atp1b1","Trerf1","Sik2","Ryr2","Rheb","Nrn1","Grik2",
  "Hspa8","Ralgapa2","Plk2","Pank1","Atp1a1","Kcnv1","Frmpd4","Hmgcr","Gpr158","Pam",
  "Pim1","Frmd4a","Ppp6r2","Dnajc1","Baz1a","Ier5","Ptgs2","Cplane1","Rps6ka2"
)

## ---- UR gene sets ----
ur_deg_common <- c("Fos","Tiparp","Nr4a1","Egr4","Egr1","Abra","Bdnf","Homer1","Dusp6","Gprc5a","Sik1","Ptgs2")

ur_deg_bla_specific <- c(
  "Egr2","Nr4a2","Adrb1","Adcyap1","Bach1","Fosb","Btg1","Rasl11a","Commd2","Thrb",
  "Syngr4","Serpinh1","Cox7a1","Itgam","Lipm","Hamp","As3mt","Samd3","Cep55","Ccdc172",
  "Gng5","Pdlim5","Trib3","Crispld1","Tnfrsf14","Cyp2j16","Gpr3","St18","Esrp1",
  "Atp1b4","Shisa8","Arid5a","Il18rap","Ptprv","Adamts1","Ephb4"
)

ur_deg_pc_specific <- c(
  "Npas4","Arc","Vgf","Zdbf2","Ntrk2","Etv5","Nr4a3","Col7a1","Pde10a","Slc6a17","Fosl2",
  "Hmgcr","R3hdm1","Rhod","Gast","Lamb3","Myoc","Itga7","Cd3e","Plet1"
)

## ---- GO-targeted sets ----
pr_GO_deg_common <- c("Arc","Bdnf","Egr1","Fos","Nr4a1","Ntrk2")
pr_GO_deg_pc_specific <- c("Adgtf5","Atf3","Atp1a1","Ccn1","Cdkn1a","Egr2","Frmpd4","Grik2",
                           "Hmga1","Hmgcr","Homer1","Hspa8","Irs2","Nr4a2","Nr4a3","Numb","Pim1",
                           "Plk2","Ptgs2","Rheb","Sik1","Sik2","Smad7","Spry2","Trerf1","Vgf")
pr_GO_deg_bla_specific <- c("Bdkrb2","Bmp2","Ncam1","Irgm","Nptx2","Rab5al1","Snap25")

ur_GO_common <- character(0)
ur_GO_bla_specific <- c("Adamts1","Adcyap1","Adrb1","As3mt","Bach1","Bdnf","Ccdc172","Crispld1",
                        "Cyp2j16","Dusp6","Egr1","Egr2","Egr4","Ephb4","Fos","Gpr3","Gprc5a",
                        "Hamp","Homer1","Itgam","Nr4a1","Nr4a2","Pdlim5","Ptgs2","Ptprv",
                        "Serpinh1","Shisa8","Sik1","St18","Trib3")
ur_GO_pc_specific <- c("Rhod","Itga7","Cd3e","Plet1","Lamb3","Myoc","Col7a1")

## ---- PE sets ----
pe_common <- c("Arc","Bdnf","Dusp6","Egr1","Egr2","Egr4","Fos","Fosl2","Homer1","Nr4a1","Nr4a2","Nr4a3",
               "Ntrk2","Pde10a","Ptgs2","Trib1","Vgf")

pe_bla_specific <- c("Alox15","Ankrd30a","Cavin3","Cd1d1","Cd300lg","Ecscr","Ednrb","Fgf21","Nid2","Nrarp","Olig2","Twist1")
pe_pc_specific <- c("Abhd12b","AC120310.1","Adgrf1","Adora3","Arhgap42","Cdh13","Cenpm","Cib4","Diras2","Epor",
                    "Fam110d","Gadd45b","Gpr39","Gpr62","Hectd2","Hspb3","Igf1r","Lonrf3","Lpar3","Mas1","Ms4a2",
                    "Nfil3","Npas4","Osbpl8","Pam","Phf21b","Prok2","R3hdm1","Rapgef5","Scimp","Sorcs1","Tll1",
                    "Tpbgl","Tpte2","Tusc1","Zdbf2","Tiparp")

## =========================================================
## 2) Helpers
## =========================================================
## ---- Seurat v5 layer-safe counts pull ----
get_counts_layer <- function(obj, assay = "RNA", layer_preference = c("counts","data")) {
  a <- obj[[assay]]
  lay <- layer_preference[layer_preference %in% names(a@layers)][1]
  if (is.na(lay)) stop("No requested layers found in assay@layers. Available: ", paste(names(a@layers), collapse=", "))
  m <- a@layers[[lay]]
  rownames(m) <- rownames(a)
  colnames(m) <- colnames(obj)
  m
}

get_group_cpm <- function(counts_mat, obj, group_label, group_col = "behavior_group") {
  stopifnot(group_col %in% colnames(obj@meta.data))
  cells <- colnames(obj)[obj@meta.data[[group_col]] == group_label]
  if (length(cells) == 0) return(setNames(rep(NA_real_, nrow(counts_mat)), rownames(counts_mat)))
  
  sub  <- counts_mat[, cells, drop = FALSE]
  sums <- Matrix::rowSums(sub)
  lib  <- sum(sums)
  if (lib == 0) return(setNames(rep(0, length(sums)), names(sums)))
  sums / lib * 1e6
}

make_rna_cpm_table <- function(bla_obj, pc_obj, group_label, assay = "RNA") {
  DefaultAssay(bla_obj) <- assay
  DefaultAssay(pc_obj)  <- assay
  
  bla_counts <- get_counts_layer(bla_obj, assay = assay)
  pc_counts  <- get_counts_layer(pc_obj,  assay = assay)
  
  common_genes <- intersect(rownames(bla_counts), rownames(pc_counts))
  bla_counts <- bla_counts[common_genes, , drop = FALSE]
  pc_counts  <- pc_counts[common_genes,  , drop = FALSE]
  
  cpm_bla <- get_group_cpm(bla_counts, bla_obj, group_label)
  cpm_pc  <- get_group_cpm(pc_counts,  pc_obj,  group_label)
  
  data.frame(
    gene        = common_genes,
    RNA_CPM_BLA = as.numeric(cpm_bla[common_genes]),
    RNA_CPM_PC  = as.numeric(cpm_pc[common_genes]),
    logCPM_BLA  = log1p(as.numeric(cpm_bla[common_genes])),
    logCPM_PC   = log1p(as.numeric(cpm_pc[common_genes])),
    group       = group_label,
    stringsAsFactors = FALSE
  )
}

read_peak_csv_to_gene_cpm <- function(path, region_label,
                                      groups = c("U_R","P_R","P_A","U_A","U_E","P_E")) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  stopifnot("SYMBOL" %in% colnames(df))
  
  want <- paste0("cpm_", groups)
  have <- intersect(want, colnames(df))
  
  df %>%
    filter(!is.na(SYMBOL) & SYMBOL != "") %>%
    group_by(SYMBOL) %>%
    summarise(across(all_of(have), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    rename(gene = SYMBOL) %>%
    rename_with(~ paste0("ATAC_", region_label, "_", .x), all_of(have))
}

annotate_deg_table <- function(df, deg_common, deg_bla, deg_pc) {
  all_genes <- unique(c(deg_common, deg_bla, deg_pc))
  df %>%
    filter(gene %in% all_genes) %>%
    mutate(
      category = case_when(
        gene %in% deg_common ~ "common",
        gene %in% deg_bla    ~ "BLA-specific",
        gene %in% deg_pc     ~ "PC-specific",
        TRUE ~ NA_character_
      ),
      category = factor(category, levels = c("common","BLA-specific","PC-specific"))
    ) %>%
    arrange(category, gene)
}

load_and_merge_group <- function(group_label, rna_df_all, atac_gene_df) {
  # RNA subset for group
  rna <- rna_df_all %>%
    filter(group == group_label) %>%
    select(gene, RNA_CPM_BLA, RNA_CPM_PC, logCPM_BLA, logCPM_PC)
  
  # ATAC columns (gene-level)
  col_bla <- paste0("ATAC_BLA_cpm_", group_label)
  col_pc  <- paste0("ATAC_PC_cpm_",  group_label)
  
  atac2 <- atac_gene_df %>%
    mutate(
      ATAC_CPM_BLA = if (col_bla %in% colnames(atac_gene_df)) .data[[col_bla]] else 0,
      ATAC_CPM_PC  = if (col_pc  %in% colnames(atac_gene_df)) .data[[col_pc]]  else 0
    ) %>%
    select(gene, ATAC_CPM_BLA, ATAC_CPM_PC)
  
  full_join(rna, atac2, by = "gene") %>%
    mutate(across(where(is.numeric), ~ replace_na(.x, 0)))
}

cor_pair <- function(x, y, method = c("spearman","pearson")) {
  method <- match.arg(method)
  ok <- complete.cases(x, y)
  if (sum(ok) < 3) return(list(estimate = NA_real_, p_value = NA_real_, n = sum(ok)))
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method, exact = FALSE))
  list(estimate = unname(ct$estimate), p_value = ct$p.value, n = sum(ok))
}

cor_summary_for_sheet <- function(df, sheet_name,
                                  RNA_BLA = "RNA_CPM_BLA", RNA_PC = "RNA_CPM_PC",
                                  ATAC_BLA = "ATAC_CPM_BLA", ATAC_PC = "ATAC_CPM_PC") {
  df <- df %>%
    mutate(
      RNA_CPM_BLA  = as.numeric(.data[[RNA_BLA]]),
      RNA_CPM_PC   = as.numeric(.data[[RNA_PC]]),
      ATAC_CPM_BLA = as.numeric(.data[[ATAC_BLA]]),
      ATAC_CPM_PC  = as.numeric(.data[[ATAC_PC]]),
      log_RNA_BLA  = log10(RNA_CPM_BLA  + 1),
      log_RNA_PC   = log10(RNA_CPM_PC   + 1),
      log_ATAC_BLA = log10(ATAC_CPM_BLA + 1),
      log_ATAC_PC  = log10(ATAC_CPM_PC  + 1)
    )
  
  combos <- list(
    list(assay="RNA",  scale="raw",   x=df$RNA_CPM_BLA,  y=df$RNA_CPM_PC),
    list(assay="RNA",  scale="log10", x=df$log_RNA_BLA,  y=df$log_RNA_PC),
    list(assay="ATAC", scale="raw",   x=df$ATAC_CPM_BLA, y=df$ATAC_CPM_PC),
    list(assay="ATAC", scale="log10", x=df$log_ATAC_BLA, y=df$log_ATAC_PC)
  )
  
  out <- list()
  for (cmb in combos) {
    for (m in c("spearman","pearson")) {
      r <- cor_pair(cmb$x, cmb$y, method = m)
      out[[length(out) + 1]] <- data.frame(
        sheet = sheet_name, assay = cmb$assay, scale = cmb$scale,
        method = m, estimate = r$estimate, p_value = r$p.value, n_genes = r$n
      )
    }
  }
  bind_rows(out)
}

run_xlsx_correlation_summary <- function(xlsx_path, out_csv) {
  sheets <- openxlsx::getSheetNames(xlsx_path)
  df_list <- lapply(sheets, function(s) openxlsx::read.xlsx(xlsx_path, sheet = s))
  names(df_list) <- sheets
  
  cor_all <- bind_rows(lapply(names(df_list), function(s) cor_summary_for_sheet(df_list[[s]], s)))
  readr::write_csv(cor_all, out_csv)
  cor_all
}

## =========================================================
## 3) RNA pseudo-bulk CPM (BLA vs PC) for configured groups
## =========================================================
message("Step 1/4: RNA CPM by group ...")

rna_all <- bind_rows(lapply(CFG$rna_groups, function(g) {
  make_rna_cpm_table(BLA_VGLUT1, PC_VGLUT1, g)
}))

rna_out_csv <- file.path(CFG$out_rna, "BLA_PC_VGLUT1_RNA_CPM_by_group.csv")
readr::write_csv(rna_all, rna_out_csv)

## =========================================================
## 4) ATAC gene-level CPM from peak CSVs
## =========================================================
message("Step 2/4: ATAC gene CPM from peak CSVs ...")

if (!file.exists(CFG$bla_peak_csv)) stop("Missing file: ", CFG$bla_peak_csv)
if (!file.exists(CFG$pc_peak_csv))  stop("Missing file: ", CFG$pc_peak_csv)

bla_gene <- read_peak_csv_to_gene_cpm(CFG$bla_peak_csv, "BLA")
pc_gene  <- read_peak_csv_to_gene_cpm(CFG$pc_peak_csv,  "PC")

atac_gene <- full_join(bla_gene, pc_gene, by = "gene") %>%
  mutate(across(starts_with("ATAC_"), ~ replace_na(.x, 0)))

atac_out_csv <- file.path(CFG$out_atac, "BLA_PC_VGLUT1_ATAC_gene_CPM.csv")
readr::write_csv(atac_gene, atac_out_csv)

## =========================================================
## 5) Merge RNA+ATAC per group and export Excel workbooks
## =========================================================
message("Step 3/4: Merge + export Excel ...")

df_PR <- load_and_merge_group("P_R", rna_all, atac_gene)
df_UR <- load_and_merge_group("U_R", rna_all, atac_gene)
df_PA <- load_and_merge_group("P_A", rna_all, atac_gene)

sheets_PR_UR <- list(
  PR_group_PR_GO_DEGs = annotate_deg_table(df_PR, pr_GO_deg_common, pr_GO_deg_bla_specific, pr_GO_deg_pc_specific),
  PR_group_PR_DEGs    = annotate_deg_table(df_PR, pr_deg_common,    pr_deg_bla_specific,    pr_deg_pc_specific),
  UR_group_PR_GO_DEGs = annotate_deg_table(df_UR, pr_GO_deg_common, pr_GO_deg_bla_specific, pr_GO_deg_pc_specific),
  UR_group_PR_DEGs    = annotate_deg_table(df_UR, pr_deg_common,    pr_deg_bla_specific,    pr_deg_pc_specific),
  
  PR_group_UR_GO_DEGs = annotate_deg_table(df_PR, ur_GO_common,     ur_GO_bla_specific,     ur_GO_pc_specific),
  PR_group_UR_DEGs    = annotate_deg_table(df_PR, ur_deg_common,    ur_deg_bla_specific,    ur_deg_pc_specific),
  UR_group_UR_GO_DEGs = annotate_deg_table(df_UR, ur_GO_common,     ur_GO_bla_specific,     ur_GO_pc_specific),
  UR_group_UR_DEGs    = annotate_deg_table(df_UR, ur_deg_common,    ur_deg_bla_specific,    ur_deg_pc_specific),
  
  UR_global = df_UR,
  PR_global = df_PR
)

xlsx1 <- file.path(CFG$out_exp, "BLA_PC_VGLUT1_PR_UR_DEG_GO_RNA_ATAC_CPM.xlsx")
openxlsx::write.xlsx(sheets_PR_UR, file = xlsx1, overwrite = TRUE)

sheets_PA <- list(
  PA_group_DEGs = annotate_deg_table(df_PA, pe_common, pe_bla_specific, pe_pc_specific),
  PA_global     = df_PA
)
xlsx2 <- file.path(CFG$out_exp, "BLA_PC_VGLUT1_PA_RNA_ATAC_CPM.xlsx")
openxlsx::write.xlsx(sheets_PA, file = xlsx2, overwrite = TRUE)

## =========================================================
## 6) Correlation summaries across Excel sheets
## =========================================================
message("Step 4/4: Correlation summaries ...")

out_cor1 <- file.path(CFG$out_exp, "BLA_PC_VGLUT1_PR_UR_DEG_GO_RNA_ATAC_correlation_summary.csv")
cor1 <- run_xlsx_correlation_summary(xlsx1, out_cor1)

out_cor2 <- file.path(CFG$out_exp, "BLA_PC_VGLUT1_PA_RNA_ATAC_correlation_summary.csv")
cor2 <- run_xlsx_correlation_summary(xlsx2, out_cor2)

message("✅ DONE
- RNA CPM:  ", rna_out_csv, "
- ATAC CPM: ", atac_out_csv, "
- Excel:    ", xlsx1, "
- Excel:    ", xlsx2, "
- Cor:      ", out_cor1, "
- Cor:      ", out_cor2)
