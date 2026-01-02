library(here)
source(here("code", "header.R"))

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID")

outdir <- here("scratch/ext/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cellContext <- read_parquet_adt(here("data", "derived", "cellContext.parquet"))
cellContext <- cellContext[, .(ImageID, CellID, TissueStructure)]

cellsclin <- merge(
  cellsclin,
  cellContext,
  by = c("ImageID", "CellID"),
  all.x = TRUE
)

epithelial <- cellsclin[isEpithelial == TRUE]

epithelial[, isPrimaryPop := (isFOXA1Pos & isGATA3Pos & isERPos & isPRPos & isARPos & !isHER2Pos)]
epithelial[, isSecondaryPop := (!isFOXA1Pos & !isGATA3Pos & isERPos & !isPRPos & !isARPos & !isHER2Pos)]

prepareEffectSizeData <- function(cellData) {
  patients_both_types <- cellData[TissueStructure %in% c("duct", "lobule"),
                                  .(has_duct = "duct" %in% TissueStructure,
                                    has_lobule = "lobule" %in% TissueStructure),
                                  by = PatientID][has_duct & has_lobule]$PatientID
  filtered_data <- cellData[PatientID %in% patients_both_types]
  
  all_patients <- unique(filtered_data[, .(PatientID)])
  phenotypes <- c("duct", "lobule")
  populations <- c("isPrimaryPop", "isSecondaryPop")
  complete_grid <- CJ(PatientID = all_patients$PatientID,
                      Population = populations,
                      TissueStructure = phenotypes)
  
  melted_data <- melt(filtered_data,
                      id.vars = c("PatientID", "TissueStructure"),
                      measure.vars = populations,
                      variable.name = "Population",
                      value.name = "isInPopulation")
  
  melted_data <- melted_data[isInPopulation == 1]
  
  cluster_counts <- melted_data[, .(Count = .N), by = .(PatientID, Population, TissueStructure)]
  cluster_counts_complete <- merge(complete_grid, cluster_counts,
                                   by = c("PatientID", "Population", "TissueStructure"),
                                   all = TRUE)
  cluster_counts_complete[is.na(Count), Count := 0]
  
  phenotype_totals <- filtered_data[, .(TotalCount = .N), by = .(PatientID, TissueStructure)]
  data_with_totals <- merge(cluster_counts_complete, phenotype_totals,
                            by = c("PatientID", "TissueStructure"))
  
  data_with_totals[, Proportion := Count / TotalCount]
  data_with_totals <- unique(data_with_totals)
  
  return(data_with_totals)
}

summarizeEffectSize <- function(data) {
  summary_results <- data.frame(
    Population = character(),
    MedianProportionDuct = numeric(),
    MedianProportionLobule = numeric(),
    CliffDelta = numeric(),
    LowerCI = numeric(),
    UpperCI = numeric(),
    pValue = numeric(),
    patients_duct = numeric(),
    patients_lobule = numeric()
  )
  
  for (population in unique(data$Population)) {
    population_data <- data[Population == population]
    
    duct_data <- population_data[TissueStructure == "duct"]
    lobule_data <- population_data[TissueStructure == "lobule"]
    merged_data <- merge(duct_data, lobule_data, by = "PatientID", suffixes = c("_duct", "_lobule"))
    
    if (nrow(merged_data) > 0) {
      median_proportion_duct <- median(merged_data$Proportion_duct)
      median_proportion_lobule <- median(merged_data$Proportion_lobule)
      
      test_result <- wilcox.test(merged_data$Proportion_lobule, merged_data$Proportion_duct, paired = TRUE)
      delta_result <- cliff.delta(merged_data$Proportion_duct, merged_data$Proportion_lobule)
      
      summary_results <- rbind(summary_results, data.frame(
        Population = population,
        MedianProportionDuct = median_proportion_duct,
        MedianProportionLobule = median_proportion_lobule,
        CliffDelta = -delta_result$estimate,
        LowerCI = -delta_result$conf.int[2],
        UpperCI = -delta_result$conf.int[1],
        pValue = test_result$p.value,
        patients_duct = nrow(duct_data),
        patients_lobule = nrow(lobule_data)
      ))
    }
  }
  
  summary_results$adjpValue <- p.adjust(summary_results$pValue, method = "BH")
  summary_results <- summary_results[order(-abs(summary_results$CliffDelta)),]
  summary_results$alpha <- ifelse(summary_results$adjpValue < 0.05, "opaque", "translucent")
  
  return(summary_results)
}

epi_merged <- prepareEffectSizeData(epithelial)
summary_results <- summarizeEffectSize(epi_merged)

createClusterBoxplot <- function(mergedData, summaryData, clusterLabels, population, ymax) {
  # filter data for the specific population
  population_specific_data <- mergedData[Population == population]
  
  population_specific_data[, Population := factor(
    Population,
    levels = c("isPrimaryPop", "isSecondaryPop"),
    labels = c("Primary (High Hormone)", "Secondary (ER+/FOXA1-/GATA3-)")
  )]
  
  population_label <- ifelse(population == "isPrimaryPop",
                             "Primary (High Hormone)",
                             "Secondary (ER+/FOXA1-/GATA3-)")
  
  adjPVal <- summaryData[summaryData$Population == population, "adjpValue"]
  if (adjPVal > 1e-4) {
    formattedPVal <- formatC(adjPVal, format = "f", digits = 1)
  } else {
    formattedPVal <- mkEnumPower(formatC(adjPVal, format = "e", digits = 0))
  }
  
  plot <- ggplot(population_specific_data,
                 aes(x = TissueStructure, y = Proportion, fill = TissueStructure)) +
    geom_boxplot(outlier.size = 0.1, outlier.color = "black") +
    scale_fill_manual(values = c(
      "duct" = getNB13Colours()$DuctLobule[["duct"]],
      "lobule" = getNB13Colours()$DuctLobule[["lobule"]]
    )) +
    theme_classic() +
    scale_x_discrete(labels = c("Duct", "Lobule")) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = c(0, ymax),
      expand = expansion(mult = c(0, 0.13)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 15),
      axis.title.y = element_text(vjust = -3, hjust = 0.5, size = 16),
      plot.title = element_text(hjust = 0.5, size = 14, vjust = -0.2),
      plot.margin = unit(c(0, 5, 10, 0), units = "pt"),
      legend.position = "none"
    )
  
  plot <- plot + annotate(
    "text",
    x = 1.5, y = ymax,
    label = TeX(paste0("p=", formattedPVal)),
    size = 5.5, color = "black", angle = 0, hjust = 0.5
  )
  
  return(plot)
}

plot1 <- createClusterBoxplot(epi_merged, summary_results, clusterLabels = epiannotations,
                              population = "isPrimaryPop", ymax = 0.15)
plot2 <- createClusterBoxplot(epi_merged, summary_results, clusterLabels = epiannotations,
                              population = "isSecondaryPop", ymax = 0.15)

combinedplots <- plot1 + plot2 + plot_layout(nrow = 1)

ggsave(combinedplots,
       width = 4.5, height = 1.75, units = "in",
       filename = file.path(outdir, "ductLobulePlot.pdf"))