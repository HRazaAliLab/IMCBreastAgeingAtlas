library(here)
source(here("code", "header.R"))

outdir <- here("scratch/main/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

annotations[, ClusterID := as.integer(as.character(ClusterID))]

cellsclin <- merge(cells, clinical, by = "ImageID")
cellsclin[, ClusterID := as.integer(as.character(ClusterID))]

positiveKi67Fractions <- getPositiveFractions(
  d = cellsclin,
  clusterColumn = "ClusterID",
  idCols = "PatientID"
)

avg_ki67 <- merge(
  positiveKi67Fractions,
  unique(clinical[, .(Age, PatientID)]),
  by = "PatientID",
  all.x = TRUE
)

avg_ki67[, AgeGroup := ifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
setnames(avg_ki67, "positiveFraction", "AvgKi67Proportion")

summary_results_ki67 <- avg_ki67[
  , {
    cluster_data <- na.omit(.SD)
    
    test_result <- wilcox.test(
      AvgKi67Proportion ~ AgeGroup,
      data = cluster_data,
      exact = FALSE
    )
    
    cliff_d_result <- cliff.delta(
      AvgKi67Proportion ~ AgeGroup,
      data = cluster_data
    )
    
    means <- tapply(cluster_data$AvgKi67Proportion, cluster_data$AgeGroup, mean)
    
    list(
      UStatistic = as.numeric(test_result$statistic),
      PValue = test_result$p.value,
      CliffsDelta = cliff_d_result$estimate,
      CID_Lower = cliff_d_result$conf.int[1],
      CID_Upper = cliff_d_result$conf.int[2],
      Mean_Postmenopausal = means["Postmenopausal"],
      Mean_Premenopausal  = means["Premenopausal"],
      patients_pre  = sum(cluster_data$AgeGroup == "Premenopausal"),
      patients_post = sum(cluster_data$AgeGroup == "Postmenopausal")
    )
  },
  by = ClusterID
]

summary_results_ki67[, Cluster := as.integer(ClusterID)]
summary_results_ki67[, adjPValue := p.adjust(PValue, method = "BH")]
summary_results_ki67[, alpha := ifelse(adjPValue < 0.05, "opaque", "translucent")]
summary_results_ki67[, point_size := 4]

summary_results_typed <- merge(
  summary_results_ki67,
  annotations,
  by.x = "Cluster",
  by.y = "ClusterID",
  all.x = TRUE
)

summary_results_typed[, Type := factor(Type, levels = c("Epithelial", "Immune", "Stromal"))]
summary_results_typed[, Cluster := factor(Cluster, levels = unique(Cluster[order(-CliffsDelta)]))]

broadcolors <- getNormalBreastProjectColours()$BroadPhenotype
summary_results_typed[, adjPValueLabel := sapply(format_custom_pval(adjPValue), mkEnumPower)]

combined_colors <- c(broadcolors, "opaque" = "black", "translucent" = "darkgray")

p_effect <- ggplot(summary_results_typed, aes(x = Cluster, y = CliffsDelta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
  geom_errorbar(aes(ymin = CID_Lower, ymax = CID_Upper, alpha = alpha), width = 0.2) +
  geom_point(shape = 21, aes(color = alpha, fill = alpha, size = point_size)) +
  scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
  scale_color_manual(values = combined_colors) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-0.65, 0.3),
    breaks = c(-0.6, -0.3, 0, 0.3),
    labels = c("-0.6", "-0.3", "0", "0.3")
  ) +
  labs(x = NULL, y = "Cliff's delta effect size (Ki67 proportion)") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    axis.text.x  = element_text(colour = "black", size = 14),
    axis.title.x = element_text(color = "black", size = 14),
    plot.margin  = unit(c(0, 45, 5, 0), "pt"),
    legend.position = "bottom",
    strip.text.y = element_blank()
  ) +
  guides(color = "none", alpha = "none", fill = "none") +
  scale_size_identity() +
  geom_text(
    aes(
      label = TeX(adjPValueLabel, output = "character"),
      color = alpha,
      y = 0.3,
      x = Cluster
    ),
    parse = TRUE, hjust = 0, size = 4.5
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free")

labels <- setNames(annotations$TeXClusterLabel, annotations$ClusterID)
p_effect <- p_effect + scale_x_discrete(labels = sapply(labels, TeX))

phenotypeclust_colorbar_plot <- ggplot(summary_results_typed, aes(x = factor(1), y = Cluster, fill = Cluster)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_y_discrete(labels = sapply(setNames(summary_results_typed$BackupTeXClusterLabel, summary_results_typed$Cluster), TeX)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = setNames(summary_results_typed$Colour, summary_results_typed$Cluster)) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(hjust = 1, color = "black", size = 14),
    axis.ticks.length.y = unit(0, "pt"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free")

summary_results_long <- melt(
  summary_results_typed,
  id.vars = c("Cluster", "Type"),
  measure.vars = c("patients_pre", "patients_post"),
  variable.name = "Group",
  value.name = "Count"
)

summary_results_long[, Group := factor(Group, levels = c("patients_post", "patients_pre"))]
levels(summary_results_long$Group) <- c("Postmenopause", "Premenopause")

bar_plot <- ggplot(summary_results_long, aes(x = Cluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Premenopause" = "gray", "Postmenopause" = "dimgray")) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 500)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 14),
    axis.title.x = element_text(colour = "black", size = 14),
    plot.margin = unit(c(30, 20, 5, 0), "pt"),
    legend.position = "none",
    strip.text.y = element_blank()
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free")

combinedplot <- phenotypeclust_colorbar_plot + p_effect + bar_plot +
  plot_layout(widths = c(4, 100, 20), ncol = 3)

proliferativefractiony <- ggplot() +
  labs(y = "higher proliferative fraction") +
  theme_classic() +
  theme(
    plot.margin = margin(0,0,0,0, unit = "cm"),
    axis.title.y = element_text(color = "black", size = 15)
  ) +
  guides(x = "none", y = "none")

combinedplot <- proliferativefractiony + combinedplot + plot_layout(widths = c(1, 1000))

ggsave(
  file.path(outdir, "Ki67ProportionEffectSizes.pdf"),
  plot = combinedplot,
  width = 5, height = 6, units = "in"
)

#Additionally print the table containing number of patients for each comparison
# n_table <- summary_results_typed[, .(
#   Cluster,
#   TeXClusterLabel,
#   BackupTeXClusterLabel,
#   Type,
#   patients_pre,
#   patients_post,
#   PValue,
#   adjPValue,
#   CliffsDelta,
#   CID_Lower,
#   CID_Upper
# )]
# 
# setnames(
#   n_table,
#   old = c("Cluster","TeXClusterLabel","BackupTeXClusterLabel","Type","patients_pre","patients_post",
#           "PValue","adjPValue","CliffsDelta","CID_Lower","CID_Upper"),
#   new = c("ClusterID","PhenotypeLabel_TeX","PhenotypeLabel_Fallback","Compartment",
#           "n_Premenopause","n_Postmenopause",
#           "p_value","p_adj_BH","cliffs_delta","CI_lower","CI_upper")
# )
# 
# write.csv(
#   n_table,
#   file = file.path(outdir, "Ki67Proliferation_n_byPhenotype.csv"),
#   row.names = FALSE
# )
