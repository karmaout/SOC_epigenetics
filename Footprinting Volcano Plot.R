df <- read_xlsx("./bindetect_results_PR_UR.xlsx")


label_col <- intersect(c("motif_name","name","TF","motif","motif_id","id"), names(df))[1]
if (is.na(label_col)) stop("No motif label column found.")

p_value_threshold <- 0.05
neglog_threshold  <- -log10(p_value_threshold)

# DATA PROCESSING ---

df <- df %>%
  mutate(
    U_R_P_R_pvalue = as.numeric(U_R_P_R_pvalue),
    neglog10_p = -log10(pmax(U_R_P_R_pvalue, 1e-300)),
    PR_minus_UR = -as.numeric(U_R_P_R_change), # Ensure direction is correct!
    
    is_AP1 = cluster == "C_FOSL1::JUN",
    is_highlighted_col = tolower(as.character(U_R_P_R_highlighted)) %in% c("true","t","1","yes","y"),
    is_significant_basic = is_highlighted_col & (neglog10_p >= neglog_threshold)
  )

# CLUSTER GROUPING LOGIC ---

sig_clusters <- df %>%
  filter(is_significant_basic, !is_AP1) %>%
  count(cluster) %>%
  arrange(desc(n))

N_top <- 12
top_clusters <- head(sig_clusters$cluster, N_top)

df <- df %>%
  mutate(
    cluster_clean = gsub("^C_", "", cluster),
    color_group = case_when(
      is_AP1 & is_highlighted_col ~ "AP-1 (Highlighted)",
      !is_significant_basic ~ "Not significant",
      cluster %in% top_clusters ~ cluster_clean,
      TRUE ~ "Other Significant"
    )
  )

# DYNAMIC COLOR PALETTE ---

all_groups <- unique(df$color_group)
cluster_groups <- setdiff(all_groups, c("AP-1 (Highlighted)", "Not significant", "Other Significant"))

my_colors <- c("Not significant" = "grey85")

if ("AP-1 (Highlighted)" %in% all_groups) {
  my_colors["AP-1 (Highlighted)"] <- "firebrick"
}

if (length(cluster_groups) > 0) {
  cluster_cols <- hue_pal()(length(cluster_groups))
  names(cluster_cols) <- cluster_groups
  my_colors <- c(my_colors, cluster_cols)
}

if ("Other Significant" %in% all_groups) {
  my_colors["Other Significant"] <- "grey40"
}

# PLOT ---

ggplot(df, aes(x = PR_minus_UR, y = neglog10_p, color = color_group)) +
  
  # 1. Background Dots (Non-significant) - INCREASED SIZE
  geom_point(data = df %>% filter(color_group == "Not significant"), 
             size = 3, alpha = 0.5) +
  
  # 2. Foreground Dots (Significant) - INCREASED SIZE to 5
  geom_point(data = df %>% filter(color_group != "Not significant"), 
             size = 5, alpha = 0.85) +
  
  # 3. Labels (AP-1 Only) - INCREASED SIZE to 7
  geom_text_repel(
    data = df %>% filter(color_group == "AP-1 (Highlighted)"),
    aes(label = .data[[label_col]]),
    size = 7,                 # Much bigger text
    fontface = "bold",        # Make it bold
    max.overlaps = 50,
    show.legend = FALSE,
    box.padding = 0.6,
    point.padding = 0.5
  ) +
  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_hline(yintercept = neglog_threshold, linetype = "dotted", color = "grey50", linewidth = 1) +
  
  scale_color_manual(values = my_colors) +
  
  theme_classic() +
  
  # 4. GLOBAL FONT SCALING
  theme(
    text = element_text(size = 20),                # Base font size for everything
    plot.title = element_text(size = 24, face = "bold"),
    plot.subtitle = element_text(size = 18, color = "grey30"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 18, color = "black"),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "right" # or "bottom" if it's too wide
  ) +
  
  labs(
    title = "Differential TF Binding by Cluster",
    subtitle = paste("Showing top", N_top, "significant clusters"),
    x = "Δ footprint score (P_R − U_R)",
    y = expression(-log[10](p-value)),
    color = "Motif Family"
  ) +
  
  # Make the dots in the legend bigger too
  guides(color = guide_legend(override.aes = list(size = 6)))


ggsave(
  filename = "volcano_footprint.tiff",
  plot = last_plot(),  # Explicitly saves the last plot
  width = 12,          # Width in inches
  height = 10,         # Height in inches
  dpi = 300,           # High resolution (300 is standard print, use 600 for extra high)
  device = "tiff",
  compression = "lzw"  # Lossless compression to keep file size down
)
