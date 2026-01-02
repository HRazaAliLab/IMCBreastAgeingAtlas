library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
panel <- getPanel()
clinical <- getClinical()
annotations <- getCellClusters()
epiannotations <- annotations[isEpithelial == TRUE]
epiannotations[, normalizedID := .GRP, by = ClusterID]

cellsclin <- merge(cells, clinical, by = "ImageID")

epithelial <- cellsclin[isEpithelial == TRUE]

epithelial[, isNoPop := (isFOXA1Pos == FALSE &
                           isGATA3Pos == FALSE &
                           isERPos == FALSE &
                           isPRPos == FALSE &
                           isARPos == FALSE &
                           isHER2Pos == FALSE)]

epithelial[, isPrimaryPop := (isFOXA1Pos == TRUE &
                                isGATA3Pos == TRUE &
                                isERPos == TRUE &
                                isPRPos == TRUE &
                                isARPos == TRUE &
                                isHER2Pos == FALSE)]

epithelial[, isSecondaryPop := (isFOXA1Pos == FALSE &
                                  isGATA3Pos == FALSE &
                                  isERPos == TRUE &
                                  isPRPos == FALSE &
                                  isARPos == FALSE &
                                  isHER2Pos == FALSE)]

unique_patients <- unique(epithelial[, .(PatientID)])
unique_clusters <- unique(epithelial[, .(ClusterID)])
full_grid <- CJ(PatientID = unique_patients$PatientID, ClusterID = unique_clusters$ClusterID)

total_primary <- epithelial[isPrimaryPop == TRUE, .(TotalPrimary = .N), by = PatientID]
total_secondary <- epithelial[isSecondaryPop == TRUE, .(TotalSecondary = .N), by = PatientID]

primary_counts <- epithelial[isPrimaryPop == TRUE, .(PrimaryCount = .N), by = .(PatientID, ClusterID)]
secondary_counts <- epithelial[isSecondaryPop == TRUE, .(SecondaryCount = .N), by = .(PatientID, ClusterID)]

proportions <- merge(full_grid, primary_counts, by = c("PatientID", "ClusterID"), all.x = TRUE)
proportions <- merge(proportions, total_primary, by = "PatientID", all.x = TRUE)
proportions <- merge(proportions, secondary_counts, by = c("PatientID", "ClusterID"), all.x = TRUE)
proportions <- merge(proportions, total_secondary, by = "PatientID", all.x = TRUE)

proportions[is.na(PrimaryCount), PrimaryCount := 0]
proportions[is.na(SecondaryCount), SecondaryCount := 0]
proportions[is.na(TotalPrimary), TotalPrimary := 0]
proportions[is.na(TotalSecondary), TotalSecondary := 0]

proportions[, PrimaryProportion := ifelse(TotalPrimary == 0, 0, PrimaryCount / TotalPrimary)]
proportions[, SecondaryProportion := ifelse(TotalSecondary == 0, 0, SecondaryCount / TotalSecondary)]

final_proportions <- proportions[, .(PatientID, ClusterID, PrimaryProportion, SecondaryProportion)]

# ------------------ effect size summary ------------------
summarizeEffectSize <- function(data) {
  summary_results <- data.table(
    Cluster = character(),
    CliffDelta = numeric(),
    LowerCI = numeric(),
    UpperCI = numeric(),
    pValue = numeric()
  )
  
  for (cluster in unique(data$ClusterID)) {
    cluster_data <- data[ClusterID == cluster]
    
    test_result <- wilcox.test(cluster_data$PrimaryProportion,
                               cluster_data$SecondaryProportion,
                               paired = TRUE)
    
    delta <- effsize::cliff.delta(cluster_data$SecondaryProportion, cluster_data$PrimaryProportion)
    
    summary_results <- rbind(
      summary_results,
      data.table(
        Cluster = as.character(cluster),
        CliffDelta = delta$estimate,
        LowerCI = delta$conf.int[1],
        UpperCI = delta$conf.int[2],
        pValue = test_result$p.value
      )
    )
  }
  
  summary_results[, adjpValue := p.adjust(pValue, method = "BH")]
  summary_results <- summary_results[order(-abs(CliffDelta))]
  summary_results[, alpha := ifelse(adjpValue < 0.05, "opaque", "translucent")]
  summary_results[]
}

plotEffectSize <- function(summaryData, clusterLabels) {
  summaryData[, point_size := 3]
  
  p <- ggplot(summaryData, aes(x = reorder(Cluster, -CliffDelta), y = CliffDelta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, alpha = alpha), width = 0.2) +
    geom_point(shape = 21, color = "black", aes(fill = alpha, size = point_size)) +
    scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = "Cliff's delta effect size") +
    theme_classic() +
    theme(
      axis.title.y = element_text(colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      axis.ticks.length.y = unit(0, "cm"),
      axis.text.x = element_text(colour = "black", size = 15),
      axis.title.x = element_text(size = 15),
      plot.margin = unit(c(40, 80, 5, 0), "pt"),
      legend.position = "bottom"
    ) +
    scale_size_identity() +
    guides(color = "none", alpha = "none", fill = "none", size = "none")
  
  p
}

createColorBarPlot <- function(summaryData, clusterLabels) {
  summaryData$Cluster <- as.integer(summaryData$Cluster)
  clusterLabels <- merge(clusterLabels, summaryData, by.x = "ClusterID", by.y = "Cluster")
  clusterLabels <- clusterLabels[order(CliffDelta)]
  clusterLabels[, TeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
  
  clusternames <- setNames(clusterLabels$TeXClusterLabel, clusterLabels$normalizedID)
  colors <- setNames(clusterLabels$Colour, clusterLabels$normalizedID)
  clusterLabels[, normalizedID := factor(normalizedID, levels = names(colors))]
  
  ggplot(clusterLabels, aes(x = factor(1), y = normalizedID, fill = normalizedID)) +
    geom_tile(color = "white", size = 0.3) +
    scale_y_discrete(limits = rev(as.factor(clusterLabels$normalizedID)),
                     labels = rev(sapply(clusternames, TeX))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 15),
      axis.ticks.length.y = unit(0, "pt"),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
}

# ------------------ build + save ------------------
summary_results <- summarizeEffectSize(final_proportions)
summary_results$adjpValueLabel <- sapply(format_custom_pval(summary_results$adjpValue), mkEnumPower)

effectSizePlot <- plotEffectSize(summary_results, epiannotations) +
  scale_y_continuous(breaks = c(-0.7, 0, 0.7),
                     labels = function(x) ifelse(x == 0, "0", format(x))) +
  geom_text(
    aes(y = 0.9,
        x = reorder(Cluster, -CliffDelta),
        label = TeX(adjpValueLabel, output = "character"),
        color = ifelse(alpha == "translucent", "darkgray", "black")),
    parse = TRUE, hjust = 0, size = 5
  ) +
  scale_color_identity()

epicolorbarplot <- createColorBarPlot(summary_results, epiannotations)
epiplot <- epicolorbarplot + effectSizePlot + plot_layout(widths = c(0.3, 7))

pdf(file.path(outdir, "popClusterEnrichment.pdf"), width = 4, height = 3.5)
print(epiplot)
dev.off()

# cluster enrichment boxplots (Primary vs Secondary) ------------------

# long data: Patient x Cluster x (HR+/ER+)
dt <- melt(
  final_proportions,
  id.vars = c("PatientID", "ClusterID"),
  measure.vars = c("PrimaryProportion", "SecondaryProportion"),
  variable.name = "Population",
  value.name = "Proportion"
)[, PopulationLabel := factor(
  Population,
  levels = c("PrimaryProportion", "SecondaryProportion"),
  labels = c("HR+", "ER+")
)]

labs_dt <- copy(epiannotations)[, .(
  ClusterID = as.character(ClusterID),
  TeXClusterLabel = mkTeX(BackupLabel, ignoreChar = "/")
)]
dt[, ClusterID := as.character(ClusterID)]
dt <- merge(dt, labs_dt, by = "ClusterID", all.x = TRUE)

# choose 4 clusters
clusters_to_plot <- c("1", "2", "5", "6")
dt <- dt[ClusterID %chin% clusters_to_plot]
dt[, ClusterID := factor(ClusterID, levels = clusters_to_plot)]

# paired Wilcoxon p per cluster + TeX label
wide <- dcast(dt, PatientID + ClusterID + TeXClusterLabel ~ PopulationLabel, value.var = "Proportion")
pvals <- wide[, .(p = wilcox.test(`HR+`, `ER+`, paired = TRUE)$p.value), by = .(ClusterID, TeXClusterLabel)]
pvals[, pLabelChar := TeX(paste0("p=", sapply(format_custom_pval(p), mkEnumPower)), output = "character")]

ymax <- 1.0

createClusterBoxplot_HRvsER <- function(data, ptab, cid, ymax) {
  d <- data[ClusterID == cid]
  ggplot(d, aes(PopulationLabel, Proportion, fill = PopulationLabel)) +
    geom_boxplot(outlier.size = 0.1, outlier.color = "black") +
    scale_fill_manual(values = c("HR+" = "gray", "ER+" = "#8B4513")) +
    scale_y_continuous(
      limits = c(0, 1.2), breaks = c(0, ymax),
      expand = expansion(mult = c(0, 0.1)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title =TeX(d[, unique(TeXClusterLabel)][1])) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(vjust = -3, hjust = 0.5, size = 16),
      plot.title = element_text(hjust = 0.5, size = 14, vjust = -0.2),
      plot.margin = unit(c(0, 0, 10, 0), "pt"),
      legend.position = "none"
    ) +
    geom_text(
      data = data.frame(x = 1.5, y = ymax, lab = ptab[ClusterID == cid, pLabelChar][1]),
      aes(x = x, y = y, label = lab),
      inherit.aes = FALSE, parse = TRUE, size = 4.5, hjust = 0.5, vjust = 0
    )
}

p <- (createClusterBoxplot_HRvsER(dt, pvals, clusters_to_plot[1], ymax) +
        createClusterBoxplot_HRvsER(dt, pvals, clusters_to_plot[2], ymax)) /
  (createClusterBoxplot_HRvsER(dt, pvals, clusters_to_plot[3], ymax) +
     createClusterBoxplot_HRvsER(dt, pvals, clusters_to_plot[4], ymax))

ylab <- ggplot() +
  labs(y = "Phenotype prop. in\nHR-population") +
  theme_classic() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "pt"),
        axis.title.y = element_text(color = "black", size = 15)) +
  guides(x = "none", y = "none")

ggsave(
  ylab + p + plot_layout(widths = c(1, 1000)),
  filename = file.path(outdir, "boxplotsPopClusterEnrichment.pdf"),
  width = 4.2, height = 3.6, units = "in"
)
