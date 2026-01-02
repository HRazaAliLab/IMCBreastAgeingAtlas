library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellsclin[, AgeGroup := fifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
cellsclin[, ClusterID := as.integer(ClusterID)]

avg_area <- cellsclin[, .(
  AvgArea = mean(CellArea, na.rm = TRUE)
), by = .(PatientID, ClusterID, isKi67Pos)]

summary_results <- rbindlist(lapply(unique(avg_area$ClusterID), function(cluster) {
  dtc <- avg_area[ClusterID == cluster, .(PatientID, isKi67Pos, AvgArea)]
  # wide
  wide <- dcast(dtc, PatientID ~ isKi67Pos, value.var = "AvgArea")
  if (!("FALSE" %in% names(wide)) || !("TRUE" %in% names(wide))) {
    return(data.table(
      Cluster = as.character(cluster),
      PValue = NA_real_, CliffsDelta = NA_real_,
      CID_Lower = NA_real_, CID_Upper = NA_real_,
      Mean_Positive = NA_real_, Mean_Negative = NA_real_,
      patients_negative = sum(!is.na(wide$`FALSE`)),
      patients_positive = sum(!is.na(wide$`TRUE`))
    ))
  }
  setnames(wide, c("FALSE", "TRUE"), c("negative", "positive"))
  wide <- wide[!is.na(negative) & !is.na(positive)]
  
  if (nrow(wide) <= 1) {
    return(data.table(
      Cluster = as.character(cluster),
      PValue = NA_real_, CliffsDelta = NA_real_,
      CID_Lower = NA_real_, CID_Upper = NA_real_,
      Mean_Positive = mean(wide$positive, na.rm = TRUE),
      Mean_Negative = mean(wide$negative, na.rm = TRUE),
      patients_negative = sum(!is.na(wide$negative)),
      patients_positive = sum(!is.na(wide$positive))
    ))
  }
  
  ttest <- t.test(wide$positive, wide$negative, paired = TRUE)
  cliff <- cliff.delta(wide$positive, wide$negative, paired = TRUE)
  
  data.table(
    Cluster = as.character(cluster),
    PValue = ttest$p.value,
    CliffsDelta = cliff$estimate,
    CID_Lower = cliff$conf.int[1],
    CID_Upper = cliff$conf.int[2],
    Mean_Positive = mean(wide$positive, na.rm = TRUE),
    Mean_Negative = mean(wide$negative, na.rm = TRUE),
    patients_negative = sum(!is.na(wide$negative)),
    patients_positive = sum(!is.na(wide$positive))
  )
}))

summary_results[, adjpValue := p.adjust(PValue, method = "BH")]
summary_results[, alpha := fifelse(!is.na(adjpValue) & adjpValue < 0.05, "opaque", "translucent")]
summary_results[, adjPValueLabel := sapply(format_custom_pval(adjpValue), mkEnumPower)]

summary_results[, point_size := 4]

summary_results[, Cluster := as.integer(Cluster)]

summary_results_typed <- merge(summary_results, annotations, by.x = "Cluster", by.y = "ClusterID", all.x = TRUE)
summary_results_typed <- as.data.table(summary_results_typed)
summary_results_typed[, Type := factor(Type, levels = c("Epithelial", "Immune", "Stromal"))]
summary_results_typed[, Cluster := factor(Cluster, levels = unique(Cluster[order(Type, -CliffsDelta)]))]

labels <- setNames(annotations$BackupLabel, annotations$ClusterID)

broadcolors <- getNormalBreastProjectColours()$BroadPhenotype
combined_colors <- c(broadcolors, "opaque" = "black", "translucent" = "darkgray")

p_effect <- ggplot(summary_results_typed, aes(x = Cluster, y = CliffsDelta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
  geom_errorbar(aes(ymin = CID_Lower, ymax = CID_Upper, alpha = alpha), width = 0.2) +
  geom_point(shape = 21, aes(color = alpha, fill = alpha, size = point_size)) +
  scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
  scale_color_manual(values = combined_colors) +
  coord_flip(clip = "off") +
  scale_y_continuous(limits = c(-0.3, 0.8),
                     breaks = c(-0.3, 0, 0.8),
                     labels = function(x) ifelse(x == 0, "0", format(x))) +
  labs(x = "enrichment of larger cells", y = "Cliff's delta effect size") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 12),
    plot.margin = unit(c(0,40,0,0), "pt"),
    legend.position = "bottom",
    strip.text.y = element_blank()
  ) +
  guides(color = "none", alpha = "none", fill = "none") +
  scale_size_identity() +
  geom_text(aes(label = TeX(adjPValueLabel, output = "character"), y = 0.8, color = alpha),
            hjust = 0, size = 3.5, parse = TRUE) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free") +
  scale_x_discrete(labels = sapply(labels[levels(summary_results_typed$Cluster)], function(x) TeX(mkTeX(x, ignoreChar="/"))))

phenotypecolors <- setNames(summary_results_typed$Colour, summary_results_typed$Cluster)
phenotypelabels <- setNames(summary_results_typed$BackupLabel, summary_results_typed$Cluster)

p_colorbar <- ggplot(summary_results_typed, aes(x = factor(1), y = Cluster, fill = Cluster)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_y_discrete(labels = sapply(phenotypelabels[levels(summary_results_typed$Cluster)], function(x) TeX(mkTeX(x, ignoreChar="/")))) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = phenotypecolors) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(hjust = 1, color = "black", size = 12),
    axis.ticks.length.y = unit(0, "pt"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free")

counts_pos <- summary_results_typed[, .(Cluster, Type, patients_positive)]
counts_pos[, Group := "patients_positive"]
setnames(counts_pos, "patients_positive", "Count")

p_counts <- ggplot(counts_pos, aes(x = Cluster, y = Count)) +
  geom_bar(stat = "identity", fill = "gray") +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 300)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12),
    plot.margin = unit(c(30, 20, 5, 15), "pt"),
    legend.position = "none",
    strip.text.y = element_blank()
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free")

combined <- p_colorbar + p_effect + p_counts + plot_layout(widths = c(4, 100, 20), ncol = 3)
yaxislabel <- ggplot() +
  labs(y = "enrichment of larger cells") +
  theme_classic() +
  theme(plot.margin = unit(c(0,0,0,0), "pt"),
        axis.title.y = element_text(color = "black", size = 15)) +
  guides(x = "none", y = "none")
combined <- yaxislabel + combined + plot_layout(widths = c(1, 1000))

ggsave(file.path(outdir, "cellAreaKi67PositivityPhenotypes.pdf"),
       plot = combined, width = 5, height = 6, units = "in")

#also write table of n by phenotype (Ki67+ vs Ki67- paired) + stats
# summary_results_typed[, Cluster_plot_order := as.integer(Cluster)]
# 
# n_table <- summary_results_typed[, .(
#   Cluster,
#   BackupLabel,
#   Type,
#   patients_negative,
#   patients_positive,
#   PValue,
#   adjpValue,
#   CliffsDelta,
#   CID_Lower,
#   CID_Upper,
#   Mean_Negative,
#   Mean_Positive,
#   Cluster_plot_order
# )]
# 
# setnames(
#   n_table,
#   old = c(
#     "Cluster", "BackupLabel", "Type",
#     "patients_negative", "patients_positive",
#     "PValue", "adjpValue",
#     "CliffsDelta", "CID_Lower", "CID_Upper",
#     "Mean_Negative", "Mean_Positive",
#     "Cluster_plot_order"
#   ),
#   new = c(
#     "ClusterID", "PhenotypeLabel_Fallback", "Compartment",
#     "n_Ki67_negative", "n_Ki67_positive",
#     "p_value", "p_adj_BH",
#     "cliffs_delta", "CI_lower", "CI_upper",
#     "mean_Ki67_negative", "mean_Ki67_positive",
#     "plot_order"
#   )
# )
# 
# fwrite(
#   n_table,
#   file = file.path(outdir, "ExtFig5dKi67PositiveNegativeCellSizePhenotypes_n_byPhenotype.csv")
# )
