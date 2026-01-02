library(here)
source(here("code", "header.R"))

cells    <- getCells()
clinical <- getClinical()
annotations <- getCellClusters()

epiannotations <- annotations[isEpithelial == TRUE]
tmeannotations <- annotations[isEpithelial == FALSE]

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
  merged_data[, Density := 1e6 * nCellsPerType / totalTissueArea]  # cells/mm^2
  
  merged_data
}

summarizeEffectSize <- function(data) {
  summary_results <- data.frame(
    Cluster    = character(),
    CliffDelta = numeric(),
    LowerCI    = numeric(),
    UpperCI    = numeric(),
    pValue     = numeric()
  )
  
  for (cluster in unique(data$ClusterID)) {
    cluster_data <- data[ClusterID == cluster]
    
    test_result <- wilcox.test(Density ~ AgeGroup, data = cluster_data)
    delta <- cliff.delta(Density ~ AgeGroup, data = cluster_data)
    
    summary_results <- rbind(summary_results, data.frame(
      Cluster    = cluster,
      CliffDelta = delta$estimate,
      LowerCI    = delta$conf.int[1],
      UpperCI    = delta$conf.int[2],
      pValue     = test_result$p.value
    ))
  }
  
  summary_results$adjpValue <- p.adjust(summary_results$pValue, method = "BH")
  summary_results <- summary_results[order(-abs(summary_results$CliffDelta)), ]
  summary_results
}

createClusterBoxplot <- function(mergedData, summaryData, clusterLabels, cluster_id, ymax) {
  cluster_specific_data <- mergedData[ClusterID == cluster_id]
  cluster_specific_data <- na.omit(cluster_specific_data, cols = "Age")
  
  cluster_specific_data[, Age_Group := ifelse(Age >= 50, "postmenopausal", "premenopausal")]
  cluster_specific_data[, Age_Group := factor(Age_Group, levels = c("premenopausal", "postmenopausal"))]
  
  clusterLabels[, TeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
  label <- clusterLabels[ClusterID == cluster_id, TeXClusterLabel]
  
  adjPVal <- summaryData[summaryData$Cluster == cluster_id, "adjpValue"]
  if (is.data.frame(adjPVal) || data.table::is.data.table(adjPVal)) adjPVal <- adjPVal[[1]]
  adjPVal <- as.numeric(adjPVal[1])
  
  formattedPVal <- if (!is.na(adjPVal) && adjPVal > 1e-4) {
    formatC(adjPVal, format = "f", digits = 2)
  } else if (!is.na(adjPVal)) {
    mkEnumPower(formatC(adjPVal, format = "e", digits = 0))
  } else {
    "NA"
  }
  
  ggplot(cluster_specific_data, aes(x = Age_Group, y = Density, fill = Age_Group)) +
    geom_boxplot(outlier.size = 0.1, outlier.color = "black") +
    scale_fill_manual(values = c("premenopausal" = "gray", "postmenopausal" = "dimgray")) +
    theme_classic() +
    scale_x_discrete(labels = c("pre", "post")) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = c(0, floor(ymax / 0.05) * 0.05),
      expand = expand_scale(mult = c(0, 0.1)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title = TeX(label)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(vjust = -3, hjust = 0.5, size = 16),
      plot.title = element_text(hjust = 0.5, size = 14, vjust = -0.2),
      plot.margin = unit(c(0, 0, 10, 0), "pt"),
      legend.position = "none"
    ) +
    annotate(
      "text", x = 1.5, y = ymax,
      label = TeX(paste0("p=", formattedPVal)),
      size = 5, color = "black", hjust = 0.5
    )
}

epi_merged_data <- prepareEffectSizeData(epithelial, tissueAreasPatients)
epi_summary_results <- summarizeEffectSize(epi_merged_data)

micro_merged_data <- prepareEffectSizeData(microenvironment, tissueAreasPatients)
micro_summary_results <- summarizeEffectSize(micro_merged_data)

epiboxplot1 <- createClusterBoxplot(epi_merged_data, epi_summary_results, epiannotations, 9,  ymax = 400)
epiboxplot2 <- createClusterBoxplot(epi_merged_data, epi_summary_results, epiannotations, 5,  ymax = 500)
epiboxplot3 <- createClusterBoxplot(epi_merged_data, epi_summary_results, epiannotations, 3,  ymax = 200)
epiboxplot4 <- createClusterBoxplot(epi_merged_data, epi_summary_results, epiannotations, 10, ymax = 400)
epicombined <- epiboxplot1 + epiboxplot2 + epiboxplot3 + epiboxplot4 + plot_layout(nrow = 1)

tmeboxplot1 <- createClusterBoxplot(micro_merged_data, micro_summary_results, tmeannotations, 14, ymax = 100)
tmeboxplot2 <- createClusterBoxplot(micro_merged_data, micro_summary_results, tmeannotations, 13, ymax = 100)
tmeboxplot3 <- createClusterBoxplot(micro_merged_data, micro_summary_results, tmeannotations, 20, ymax = 400)
tmeboxplot4 <- createClusterBoxplot(micro_merged_data, micro_summary_results, tmeannotations, 15, ymax = 100)
tmecombined <- tmeboxplot1 + tmeboxplot2 + tmeboxplot3 + tmeboxplot4 + plot_layout(nrow = 1)

fullboxplots <- epicombined / tmecombined

blanklabelploty <- ggplot() +
  labs(y = TeX("density (cells/mm$^2$)")) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    axis.title.y = element_text(color = "black", size = 14)
  ) +
  guides(x = "none", y = "none")

fullboxplots <- blanklabelploty + fullboxplots + plot_layout(widths = c(1, 1000))

pdf(file.path(outdir, "boxplotsPrePostMenopause.pdf"), width = 8, height = 5.5)
print(fullboxplots)
dev.off()
