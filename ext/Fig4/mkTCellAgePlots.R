library(here)
source(here("code", "header.R"))

outDir <- here("scratch/ext/Fig4")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age, TissueArea)], by = "ImageID", all.x = TRUE)

tissueAreasPatients <- clinical[, .(totalTissueArea = sum(TissueArea, na.rm = TRUE)), by = PatientID]

prepareEffectSizeData <- function(cellData) {
  cd4_t_cells <- cellData[ClusterID == 12]
  cd8_t_cells <- cellData[ClusterID == 13]
  
  cd4_proportions <- cd4_t_cells[, .(
    PD1Proportion = sum(`isPD-1Pos`, na.rm = TRUE) / .N,
    TCF1Proportion = sum(isTCF1Pos, na.rm = TRUE) / .N,
    GZMBProportion = sum(isGZMBPos, na.rm = TRUE) / .N
  ), by = PatientID]
  cd4_proportions[, ClusterID := 12]
  cd4_proportions[, CellType := "CD4"]
  
  cd8_proportions <- cd8_t_cells[, .(
    PD1Proportion = sum(`isPD-1Pos`, na.rm = TRUE) / .N,
    TCF1Proportion = sum(isTCF1Pos, na.rm = TRUE) / .N,
    GZMBProportion = sum(isGZMBPos, na.rm = TRUE) / .N
  ), by = PatientID]
  cd8_proportions[, ClusterID := 13]
  cd8_proportions[, CellType := "CD8"]
  
  merged_data <- rbind(cd4_proportions, cd8_proportions)
  merged_data[is.na(PD1Proportion), PD1Proportion := 0]
  merged_data[is.na(TCF1Proportion), TCF1Proportion := 0]
  merged_data[is.na(GZMBProportion), GZMBProportion := 0]
  
  agemapping <- unique(cellData[, .(PatientID, Age)])
  merged_data <- merge(merged_data, agemapping, by = "PatientID", all.x = TRUE)
  merged_data[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
  
  merged_data
}

prepareDensityData <- function(cellData, tissueAreasPatients) {
  cd4_t_cells <- cellData[ClusterID == 12]
  cd8_t_cells <- cellData[ClusterID == 13]
  
  cd4_counts <- cd4_t_cells[, .(
    PD1Count = sum(`isPD-1Pos`, na.rm = TRUE),
    TCF1Count = sum(isTCF1Pos, na.rm = TRUE),
    GZMBCount = sum(isGZMBPos, na.rm = TRUE)
  ), by = PatientID]
  cd4_counts[, ClusterID := 12]
  cd4_counts[, CellType := "CD4"]
  
  cd8_counts <- cd8_t_cells[, .(
    PD1Count = sum(`isPD-1Pos`, na.rm = TRUE),
    TCF1Count = sum(isTCF1Pos, na.rm = TRUE),
    GZMBCount = sum(isGZMBPos, na.rm = TRUE)
  ), by = PatientID]
  cd8_counts[, ClusterID := 13]
  cd8_counts[, CellType := "CD8"]
  
  merged_counts <- rbind(cd4_counts, cd8_counts)
  merged_counts[is.na(PD1Count), PD1Count := 0]
  merged_counts[is.na(TCF1Count), TCF1Count := 0]
  merged_counts[is.na(GZMBCount), GZMBCount := 0]
  
  merged_counts <- merge(merged_counts, tissueAreasPatients, by = "PatientID", all.x = TRUE)
  
  merged_counts[, PD1Density := PD1Count / totalTissueArea]
  merged_counts[, TCF1Density := TCF1Count / totalTissueArea]
  merged_counts[, GZMBDensity := GZMBCount / totalTissueArea]
  
  agemapping <- unique(cellData[, .(PatientID, Age)])
  merged_counts <- merge(merged_counts, agemapping, by = "PatientID", all.x = TRUE)
  merged_counts[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
  
  merged_counts
}

summarizeEffectSize <- function(data) {
  summary_results <- data.frame(
    Cluster = character(),
    Marker = character(),
    CliffDelta = numeric(),
    LowerCI = numeric(),
    UpperCI = numeric(),
    pValue = numeric(),
    AvgProportionBelow50 = numeric(),
    SEBelow50 = numeric(),
    AvgProportionAbove50 = numeric(),
    SEAbove50 = numeric(),
    patientsAbove50 = numeric(),
    patientsBelow50 = numeric()
  )
  
  markers <- c("PD1Proportion", "TCF1Proportion", "GZMBProportion")
  
  for (cluster in unique(data$ClusterID)) {
    cluster_data <- data[ClusterID == cluster]
    
    for (marker in markers) {
      cluster_data[, MarkerValue := get(marker)]
      cluster_data_marker <- cluster_data[, .(AgeGroup, MarkerValue)]
      cluster_data_marker <- na.omit(cluster_data_marker)
      
      below_50 <- cluster_data_marker[AgeGroup == "Below 50", .(
        Mean = mean(MarkerValue, na.rm = TRUE),
        SE = sd(MarkerValue, na.rm = TRUE) / sqrt(.N)
      )]
      above_50 <- cluster_data_marker[AgeGroup == "Above 50", .(
        Mean = mean(MarkerValue, na.rm = TRUE),
        SE = sd(MarkerValue, na.rm = TRUE) / sqrt(.N)
      )]
      
      test_result <- wilcox.test(MarkerValue ~ AgeGroup, data = cluster_data_marker)
      delta <- cliff.delta(MarkerValue ~ AgeGroup, data = cluster_data_marker)
      
      summary_row <- data.frame(
        Cluster = if (cluster == 12) "CD4" else "CD8",
        Marker = marker,
        CliffDelta = delta$estimate,
        LowerCI = delta$conf.int[1],
        UpperCI = delta$conf.int[2],
        pValue = test_result$p.value,
        AvgProportionBelow50 = below_50$Mean,
        SEBelow50 = below_50$SE,
        AvgProportionAbove50 = above_50$Mean,
        SEAbove50 = above_50$SE,
        patientsAbove50 = nrow(cluster_data_marker[AgeGroup == "Above 50"]),
        patientsBelow50 = nrow(cluster_data_marker[AgeGroup == "Below 50"])
      )
      
      summary_results <- rbind(summary_results, summary_row)
    }
  }
  
  summary_results$adjpValue <- p.adjust(summary_results$pValue, method = "BH")
  summary_results <- summary_results[order(-abs(summary_results$CliffDelta)), ]
  summary_results$alpha <- ifelse(summary_results$adjpValue < 0.05, "opaque", "translucent")
  summary_results
}

plotEffectSize <- function(summaryData) {
  scale_point_size <- function(LowerCI, UpperCI) {
    ci_width <- UpperCI - LowerCI
    min_size <- 2
    max_size <- 5
    max_ci_width_allowed <- max(summaryData$UpperCI - summaryData$LowerCI, na.rm = TRUE)
    size <- min_size + (max_size - min_size) * (ci_width / max_ci_width_allowed)
    size <- pmin(pmax(size, min_size), max_size)
    size
  }
  
  summaryData$point_size <- mapply(scale_point_size, summaryData$LowerCI, summaryData$UpperCI)
  
  ggplot(summaryData, aes(x = reorder(Label, -CliffDelta), y = CliffDelta)) +
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
      axis.text.y = element_text(colour = "black", size = 15),
      axis.text.x = element_text(colour = "black", size = 15),
      axis.title.x = element_text(size = 15),
      plot.margin = unit(c(40, 80, 5, 0), "pt"),
      legend.position = "bottom"
    ) +
    scale_size_identity() +
    scale_x_discrete(labels = function(x) sapply(x, TeX)) +
    guides(color = "none", alpha = "none", fill = "none", size = "none")
}

cellsclin_analysis <- na.omit(cellsclin[, .(
  PatientID, Age, ClusterID, `isPD-1Pos`, isTCF1Pos, isGZMBPos
)])
Tcell_merged_data <- prepareEffectSizeData(cellsclin_analysis)

Tcell_summary_results <- summarizeEffectSize(Tcell_merged_data)

Tcell_summary_results$adjpValueLabel <- sapply(format_custom_pval(Tcell_summary_results$adjpValue), mkEnumPower)
Tcell_summary_results$Label <- with(Tcell_summary_results, paste0("%", Marker, " ", Cluster, "+T"))
Tcell_summary_results <- Tcell_summary_results[Tcell_summary_results$Label != "%GZMBProportion CD4+T", ]
Tcell_summary_results$Label <- gsub("PD1Proportion", "PD-1", Tcell_summary_results$Label)
Tcell_summary_results$Label <- gsub("TCF1Proportion", "TCF1", Tcell_summary_results$Label)
Tcell_summary_results$Label <- gsub("GZMBProportion", "GZMB", Tcell_summary_results$Label)
Tcell_summary_results$Label <- gsub("CD4\\+T", "CD4$^+$T", Tcell_summary_results$Label)
Tcell_summary_results$Label <- gsub("CD8\\+T", "CD8$^+$T", Tcell_summary_results$Label)

TcellEffectSizePlot <- plotEffectSize(Tcell_summary_results) +
  scale_y_continuous(breaks = c(-0.4, 0, 0.3), labels = function(x) ifelse(x == 0, "0", format(x))) +
  geom_text(aes(
    y = 0.45, x = Label,
    label = TeX(adjpValueLabel, output = "character"),
    color = ifelse(alpha == "translucent", "darkgray", "black")
  ), parse = TRUE, hjust = 0, size = 5) +
  scale_color_identity()

ggsave(TcellEffectSizePlot, filename = file.path(outDir, "TAge.pdf"), units = "in", width = 4, height = 3.2)

# ----- boxplots -----
Tcell_merged_data <- na.omit(Tcell_merged_data)
Tcell_summary_results <- summarizeEffectSize(Tcell_merged_data)
Tcell_summary_results$adjpValueLabel <- sapply(format_custom_pval(Tcell_summary_results$adjpValue), mkEnumPower)
Tcell_summary_results$Label <- with(Tcell_summary_results, paste0("%", Marker, " ", Cluster, "+T"))

long_data <- melt(
  as.data.table(Tcell_merged_data),
  id.vars = c("PatientID", "Age", "AgeGroup", "ClusterID", "CellType"),
  measure.vars = c("PD1Proportion", "TCF1Proportion", "GZMBProportion"),
  variable.name = "Marker",
  value.name = "Proportion"
)
long_data$AgeGroup <- factor(long_data$AgeGroup, levels = c("Below 50", "Above 50"))

colors <- c("Below 50" = "gray", "Above 50" = "dimgray")
labels <- c("Below 50" = "pre", "Above 50" = "post")

create_boxplot <- function(data, marker, cell_type, colors, ymax_limit) {
  pLabel <- Tcell_summary_results$adjpValueLabel[
    Tcell_summary_results$Marker == marker & Tcell_summary_results$Cluster == cell_type
  ]
  
  ggplot(data[data$Marker == marker & data$CellType == cell_type, ],
         aes(x = AgeGroup, y = Proportion, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.1) +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 17),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_blank(),
      plot.margin = unit(c(0, 5, 10, 0), units = "pt"),
      legend.position = "none"
    ) +
    scale_fill_manual(values = colors, labels = labels) +
    scale_x_discrete(labels = labels) +
    scale_y_continuous(
      breaks = c(0, ymax_limit),
      limits = c(0, ymax_limit),
      expand = expand_scale(mult = c(0, 0.2)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    annotate("text", x = 1.5, y = ymax_limit,
             label = TeX(paste0("p=", pLabel), output = "character"),
             parse = TRUE, size = 5.3)
}

p1 <- create_boxplot(long_data, "PD1Proportion", "CD4", colors, 0.6)
p2 <- create_boxplot(long_data, "TCF1Proportion", "CD4", colors, 0.7)
p3 <- create_boxplot(long_data, "GZMBProportion", "CD4", colors, 0.6)
p4 <- create_boxplot(long_data, "PD1Proportion", "CD8", colors, 0.6)
p5 <- create_boxplot(long_data, "TCF1Proportion", "CD8", colors, 0.6)
p6 <- create_boxplot(long_data, "GZMBProportion", "CD8", colors, 0.6)

combined_plot <- (p1 | p2 | p3) / (p4 | p5 | p6)

ggsave(combined_plot, filename = file.path(outDir, "boxplotsTAge.pdf"), units = "in", width = 5, height = 3)
