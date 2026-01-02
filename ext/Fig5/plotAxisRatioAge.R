library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cells[, AxisRatio := MajorAxisLength / MinorAxisLength]

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellsclin[, AgeGroup := fifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
cellsclin[, ClusterID := as.integer(ClusterID)]

avg_axis <- cellsclin[, .(AvgAxisRatio = mean(AxisRatio, na.rm = TRUE)),
                      by = .(PatientID, ClusterID, AgeGroup)]

summary_results <- rbindlist(lapply(unique(avg_axis$ClusterID), function(cluster) {
  cd <- avg_axis[ClusterID == cluster]
  
  tt <- t.test(AvgAxisRatio ~ AgeGroup, data = cd, var.equal = FALSE)
  cliff <- cliff.delta(AvgAxisRatio ~ AgeGroup, data = cd,
                       pooled = FALSE, hedges.correction = FALSE)
  means <- tapply(cd$AvgAxisRatio, cd$AgeGroup, mean)
  
  data.table(
    Cluster = as.character(cluster),
    PValue = tt$p.value,
    CliffsDelta = cliff$estimate,
    CID_Lower = cliff$conf.int[1],
    CID_Upper = cliff$conf.int[2],
    Mean_Postmenopausal = means["Postmenopausal"],
    Mean_Premenopausal  = means["Premenopausal"],
    patients_pre  = nrow(cd[AgeGroup == "Premenopausal"]),
    patients_post = nrow(cd[AgeGroup == "Postmenopausal"])
  )
}))

summary_results[, adjPValue := p.adjust(PValue, method = "BH")]
summary_results[, alpha := fifelse(adjPValue < 0.05, "opaque", "translucent")]
summary_results[, adjPValueLabel := sapply(format_custom_pval(adjPValue), mkEnumPower)]

max_ci <- max(summary_results$CID_Upper - summary_results$CID_Lower, na.rm = TRUE)
summary_results[, point_size := {
  ciw <- CID_Upper - CID_Lower
  min_size <- 2; max_size <- 5
  pmin(pmax(min_size + (max_size - min_size) * (ciw / max_ci), min_size), max_size)
}]

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
  coord_flip(clip = "off", ylim = c(min(summary_results_typed$CID_Lower), max(summary_results_typed$CID_Upper))) +
  scale_y_continuous(breaks = c(-0.4, 0, 0.4), labels = function(x) ifelse(x == 0, "0", format(x))) +
  labs(x = "increasing axis ratio", y = "Cliff's delta effect size") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.title.x = element_text(color = "black", size = 12),
    plot.margin = unit(c(0, 80, 5, 0), "pt"),
    legend.position = "bottom",
    strip.text.y = element_blank()
  ) +
  guides(color = "none", alpha = "none", fill = "none") +
  scale_size_identity() +
  geom_text(aes(label = TeX(adjPValueLabel, output = "character"), y = 0.45, color = alpha),
            hjust = 0, size = 4, parse = TRUE) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free") +
  scale_x_discrete(labels = sapply(labels[levels(summary_results_typed$Cluster)], function(x) TeX(mkTeX(x, ignoreChar="/"))))

phenotypecolors <- setNames(summary_results_typed$Colour, summary_results_typed$Cluster)
phenotypelabels <- setNames(summary_results_typed$BackupLabel, summary_results_typed$Cluster)

p_colorbar <- ggplot(summary_results_typed, aes(x = factor(1), y = Cluster, fill = Cluster)) +
  geom_tile(color = NA) +
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

summary_long <- melt(
  summary_results_typed,
  id.vars = c("Cluster", "Type"),
  measure.vars = c("patients_pre", "patients_post"),
  variable.name = "Group", value.name = "Count"
)
summary_long[, Group := factor(Group, levels = c("patients_post", "patients_pre"))]
levels(summary_long$Group) <- c("Postmenopause", "Premenopause")

p_counts <- ggplot(summary_long, aes(x = Cluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Premenopause" = "gray", "Postmenopause" = "dimgray")) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 500)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 12),
    axis.title.x = element_text(colour = "black", size = 12),
    plot.margin = unit(c(30, 20, 5, 0), "pt"),
    legend.position = "none",
    strip.text.y = element_blank()
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free")

combined <- p_colorbar + p_effect + p_counts + plot_layout(widths = c(4, 100, 20), ncol = 3)
yaxislabel <- ggplot() +
  labs(y = "increasing major/minor axis ratio") +
  theme_classic() +
  theme(plot.margin = unit(c(0,0,0,0), "pt"), axis.title.y = element_text(color = "black", size = 15)) +
  guides(x = "none", y = "none")
combined <- yaxislabel + combined + plot_layout(widths = c(1, 1000))

ggsave(file.path(outdir, "axisRatioAge.pdf"),
       plot = combined, width = 5, height = 6, units = "in")

# Also write table of n by phenotype + stats

# summary_results_typed[, Cluster_plot_order := as.integer(Cluster)]
# n_table <- summary_results_typed[, .(
#   Cluster,
#   BackupLabel,
#   Type,
#   patients_pre,
#   patients_post,
#   PValue,
#   adjPValue,
#   CliffsDelta,
#   CID_Lower,
#   CID_Upper,
#   Mean_Premenopausal,
#   Mean_Postmenopausal,
#   Cluster_plot_order
# )]
# setnames(
#   n_table,
#   old = c(
#     "Cluster", "BackupLabel", "Type",
#     "patients_pre", "patients_post",
#     "PValue", "adjPValue",
#     "CliffsDelta", "CID_Lower", "CID_Upper",
#     "Mean_Premenopausal", "Mean_Postmenopausal",
#     "Cluster_plot_order"
#   ),
#   new = c(
#     "ClusterID", "PhenotypeLabel_Fallback", "Compartment",
#     "n_Premenopause", "n_Postmenopause",
#     "p_value", "p_adj_BH",
#     "cliffs_delta", "CI_lower", "CI_upper",
#     "mean_Premenopause", "mean_Postmenopause",
#     "plot_order"
#   )
# )
# data.table::fwrite(
#   n_table,
#   file = file.path(outdir, "ExtFig5gPrePostMenopausalAxisRatio_n_byPhenotype.csv")
# )
