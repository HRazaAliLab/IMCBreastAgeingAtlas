library(here)
source(here("code", "header.R"))

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

epiannotations <- annotations[isEpithelial == TRUE]
tmeannotations <- annotations[isEpithelial == FALSE]
epiannotations[, normalizedID := .GRP, by = ClusterID]
tmeannotations[, normalizedID := .GRP, by = ClusterID]

tissueAreasPatients <- clinical[, .(
  totalTissueArea = sum(TissueArea, na.rm = TRUE)
), by = PatientID]

outdir <- here("scratch/main/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)
epithelial <- cellsclin[isEpithelial == TRUE]
microenvironment <- cellsclin[isEpithelial == FALSE]

prepareEffectSizeData <- function(cellData, tissueAreasPatients) {
  merged_data <- getCellCounts(cellData, clusterColumn = "ClusterID", idCols = "PatientID")
  
  agemapping <- unique(cellData[, .(PatientID, Age)])
  merged_data <- merge(merged_data, agemapping, by = "PatientID", all.x = TRUE)
  
  merged_data[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
  setnames(merged_data, "proportion", "Proportion")
  
  merged_data <- merge(merged_data, tissueAreasPatients, by = "PatientID", all.x = TRUE)
  merged_data[, Density := 1e6 * nCellsPerType / totalTissueArea]
  
  merged_data
}

summarizeEffectSize <- function(data) {
  summary_results <- data.frame(
    Cluster                = character(),
    CliffDelta             = numeric(),
    LowerCI                = numeric(),
    UpperCI                = numeric(),
    pValue                 = numeric(),
    AvgProportionBelow50   = numeric(),
    SEProportionBelow50    = numeric(),
    AvgProportionAbove50   = numeric(),
    SEProportionAbove50    = numeric(),
    patientsAbove50        = numeric(),
    patientsBelow50        = numeric()
  )
  
  for (cluster in unique(data$ClusterID)) {
    cluster_data <- data[ClusterID == cluster]
    
    below_50_prop <- cluster_data[AgeGroup == "Below 50",
                                  .(MeanProportion = mean(Proportion, na.rm = TRUE),
                                    SEProportion   = sd(Proportion, na.rm = TRUE) / sqrt(.N))]
    above_50_prop <- cluster_data[AgeGroup == "Above 50",
                                  .(MeanProportion = mean(Proportion, na.rm = TRUE),
                                    SEProportion   = sd(Proportion, na.rm = TRUE) / sqrt(.N))]
    
    test_result <- wilcox.test(Density ~ AgeGroup, data = cluster_data)
    delta <- cliff.delta(Density ~ AgeGroup, data = cluster_data)
    
    summary_row <- data.frame(
      Cluster               = cluster,
      CliffDelta            = delta$estimate,
      LowerCI               = delta$conf.int[1],
      UpperCI               = delta$conf.int[2],
      pValue                = test_result$p.value,
      AvgProportionBelow50  = below_50_prop$MeanProportion,
      SEProportionBelow50   = below_50_prop$SEProportion,
      AvgProportionAbove50  = above_50_prop$MeanProportion,
      SEProportionAbove50   = above_50_prop$SEProportion,
      patientsAbove50       = nrow(cluster_data[AgeGroup == "Above 50"]),
      patientsBelow50       = nrow(cluster_data[AgeGroup == "Below 50"])
    )
    
    summary_results <- rbind(summary_results, summary_row)
  }
  
  summary_results$adjpValue <- p.adjust(summary_results$pValue, method = "BH")
  summary_results <- summary_results[order(-abs(summary_results$CliffDelta)), ]
  summary_results$alpha <- ifelse(summary_results$adjpValue < 0.05, "opaque", "translucent")
  summary_results$Cluster <- as.character(summary_results$Cluster)
  
  summary_results
}


plotEffectSize <- function(summaryData, clusterLabels) {
  summaryData$point_size <- 5
  labels <- setNames(clusterLabels$TeXClusterLabel, clusterLabels$normalizedID)
  
  ggplot(summaryData, aes(x = reorder(Cluster, -CliffDelta), y = CliffDelta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, alpha = alpha), width = 0.2) +
    geom_point(shape = 21, color = "black", aes(fill = alpha, size = point_size)) +
    scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
    coord_flip(
      clip = "off",
      ylim = c(min(summaryData$LowerCI), max(summaryData$UpperCI)),
      xlim = c(1, nrow(summaryData))
    ) +
    labs(x = NULL, y = "Cliff's delta effect size (Density)") +
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
    scale_x_discrete(labels = sapply(labels, TeX)) +
    guides(color = "none", alpha = "none", fill = "none", size = "none")
}

createColorBarPlot <- function(summaryData, clusterLabels) {
  dt <- data.table::as.data.table(summaryData)
  cl <- data.table::as.data.table(clusterLabels)
  
  dt[, Cluster := as.integer(Cluster)]
  mergedData <- merge(dt, cl, by.x = "Cluster", by.y = "ClusterID", all.x = TRUE)
  
  mergedData[, TeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
  
  mergedData <- mergedData[order(-CliffDelta)]
  
  map <- unique(mergedData[, .(normalizedID, TeXClusterLabel, Colour)])
  map <- map[rev(seq_len(.N))]
  
  clusternames <- setNames(map$TeXClusterLabel, map$normalizedID)
  colors       <- setNames(map$Colour,       map$normalizedID)
  
  mergedData[, normalizedID := factor(normalizedID, levels = names(colors))]
  
  ggplot(mergedData, aes(x = factor(1), y = normalizedID, fill = normalizedID)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_vline(xintercept = 1.5, color = "black", linewidth = 1) +
    scale_y_discrete(
      limits = rev(levels(mergedData$normalizedID)),
      labels = rev(sapply(clusternames, TeX))
    ) +
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


createProportionBarChart <- function(summaryData, clusterLabels, top_break) {
  dt <- data.table::as.data.table(summaryData)
  cl <- data.table::as.data.table(clusterLabels)
  
  dt[, Cluster := as.integer(Cluster)]
  dt <- merge(dt, cl, by.x = "Cluster", by.y = "ClusterID", all.x = TRUE)

  dt <- dt[order(-CliffDelta)]
  dt[, Cluster := factor(Cluster, levels = unique(Cluster))]
  
  long_summary <- data.table::melt(
    dt,
    id.vars = c(
      "Cluster", "CliffDelta", "LowerCI", "UpperCI", "pValue",
      "AvgProportionBelow50", "SEProportionBelow50",
      "AvgProportionAbove50", "SEProportionAbove50",
      "patientsAbove50", "patientsBelow50", "adjpValue", "alpha",
      "Colour", "PrintOrder", "BackupTeXClusterLabel", "normalizedID", "BackupLabel"
    ),
    measure.vars = c("AvgProportionBelow50", "AvgProportionAbove50"),
    variable.name = "AgeGroup",
    value.name = "Proportion"
  )
  
  long_summary[, AgeGroup := factor(
    AgeGroup,
    levels = c("AvgProportionAbove50", "AvgProportionBelow50"),
    labels = c("Postmenopausal", "Premenopausal")
  )]
  
  ggplot(long_summary, aes(x = Cluster, y = Proportion, fill = AgeGroup)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 15),
      axis.ticks.length.y = unit(0, "cm"),
      axis.title.x = element_text(vjust = -1, size = 13),
      plot.margin = unit(c(0.2, 10, 0, 0), "pt"),
      legend.position = "none",
      legend.title = element_blank()
    ) +
    labs(x = NULL, y = NULL) +
    scale_fill_manual(values = c("Premenopausal" = "gray", "Postmenopausal" = "dimgray")) +
    coord_flip() +
    scale_y_continuous(
      breaks = c(0, top_break),
      labels = function(x) ifelse(x == 0, "0", format(x)),
      expand = c(0, 0)
    )
}

# ============================================================
# Epithelial
# ============================================================

epi_merged_data <- prepareEffectSizeData(epithelial, tissueAreasPatients)
epi_summary_results <- adt(summarizeEffectSize(epi_merged_data))
epi_summary_results$adjpValueLabel <- sapply(format_custom_pval(epi_summary_results$adjpValue), mkEnumPower)

epi_proportion_bar_chart <- createProportionBarChart(epi_summary_results, epiannotations, top_break = 0.2)

epiEffectSizePlot <- plotEffectSize(epi_summary_results, epiannotations) +
  scale_y_continuous(
    breaks = c(-0.5, 0, 0.3),
    labels = function(x) ifelse(x == 0, "0", format(x))
  ) +
  geom_text(
    aes(
      y = 0.3,
      x = reorder(Cluster, -CliffDelta),
      label = TeX(adjpValueLabel, output = "character"),
      color = ifelse(alpha == "translucent", "darkgray", "black")
    ),
    parse = TRUE, hjust = 0, size = 5
  ) +
  scale_color_identity()

epicolorbarplot <- createColorBarPlot(epi_summary_results, epiannotations)

epicellPhenotypeplot <- epicolorbarplot + epi_proportion_bar_chart + epiEffectSizePlot +
  plot_layout(widths = c(0.3, 2, 7))

pdf(file.path(outdir, "effectSizesEpi.pdf"), width = 5.5, height = 5)
print(epicellPhenotypeplot)
dev.off()

# ============================================================
# Microenvironment
# ============================================================

micro_merged_data <- prepareEffectSizeData(microenvironment, tissueAreasPatients)
micro_summary_results <- adt(summarizeEffectSize(micro_merged_data))
micro_summary_results$adjpValueLabel <- sapply(format_custom_pval(micro_summary_results$adjpValue), mkEnumPower)

tme_proportion_bar_chart <- createProportionBarChart(micro_summary_results, tmeannotations, top_break = 0.3)

tmeEffectSizePlot <- plotEffectSize(micro_summary_results, tmeannotations) +
  coord_flip(
    clip = "off",
    ylim = c(min(micro_summary_results$LowerCI), 0.5),
    xlim = c(1, nrow(micro_summary_results))
  ) +
  scale_y_continuous(
    breaks = c(-0.5, 0, 0.5),
    labels = function(x) ifelse(x == 0, "0", format(x))
  ) +
  geom_text(
    aes(
      y = 0.5,
      x = reorder(Cluster, -CliffDelta),
      label = TeX(adjpValueLabel, output = "character"),
      color = ifelse(alpha == "translucent", "darkgray", "black")
    ),
    parse = TRUE, hjust = 0, size = 5
  ) +
  scale_color_identity()

tmecolorbarplot <- createColorBarPlot(micro_summary_results, tmeannotations)
tmeplot <- tmecolorbarplot + tme_proportion_bar_chart + tmeEffectSizePlot +
  plot_layout(widths = c(0.3, 2, 7))

pdf(file.path(outdir, "effectSizesTME.pdf"), width = 4.5, height = 4.4)
print(tmeplot)
dev.off()
