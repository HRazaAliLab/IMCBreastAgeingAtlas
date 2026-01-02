library(here)
source(here("code", "header.R"))

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

epiannotations <- annotations[isEpithelial == TRUE]
tmeannotations <- annotations[isEpithelial == FALSE]
epiannotations[, normalizedID := .GRP, by = ClusterID]
tmeannotations[, normalizedID := .GRP, by = ClusterID]

outdir <- here("scratch/main/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)
epithelial <- cellsclin[isEpithelial == TRUE]
microenvironment <- cellsclin[isEpithelial == FALSE]

# ---- Patient order from clustering on proportions ----
merged_counts <- getCellCounts(cellsclin, clusterColumn = "ClusterID", "PatientID")
setnames(merged_counts, "proportion", "ClusterProportion")
proportion_matrix <- dcast(merged_counts, PatientID ~ ClusterID, value.var = "ClusterProportion", fill = 0)

hc <- getHC(proportion_matrix[, -1])
patient_order <- proportion_matrix$PatientID[hc$order]

#to visualise the dendrogram:
#library(ggdendro)
# dendrogram_plot <- ggdendrogram(
#   as.dendrogram(hc),
#   rotate = FALSE,
#   labels = FALSE,
#   theme_dendro = FALSE
# ) +
#   theme_void() +
#   theme(plot.margin = unit(c(0, 2, 0, 0), "pt"))

age_data <- unique(cellsclin[, .(PatientID, Age = Age)])
age_data[, PatientID := factor(PatientID, levels = patient_order)]
age_data[, AgeGroup := ifelse(Age < 50, "Below 50", "Above 50")]
age_colors <- c("Below 50" = "gray", "Above 50" = "dimgray")

age_scatter_plot <- ggplot(age_data, aes(x = PatientID, y = Age, fill = AgeGroup, color = AgeGroup)) +
  geom_point(size = 0.3, shape = 21) +
  geom_smooth(
    aes(x = as.numeric(PatientID), group = 1),
    method = "loess", se = FALSE, color = "black",
    linewidth = 0.75, linetype = "dashed"
  ) +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(vjust = -1.3, size = 13),
    plot.margin = unit(c(5, 2, 0, 0), "pt"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "Age") +
  scale_y_continuous(breaks = c(20, 70), expand = c(0.02, 0.02))

makeCountsBarChart <- function(cellData, clusterLabels, patient_order) {
  cluster_counts <- getCellCounts(
    d = cellData,
    clusterColumn = "ClusterID",
    idCols = "PatientID",
    sepBy = "isEpithelial"
  )[, .(PatientID, ClusterID, nCellsPerType)]
  
  cluster_counts[, ClusterID := factor(
    ClusterID,
    levels = as.character(clusterLabels$ClusterID[order(clusterLabels$PrintOrder)])
  )]
  
  cluster_counts_long <- melt(
    cluster_counts,
    id.vars = c("PatientID", "ClusterID"),
    variable.name = "Cluster",
    value.name = "Count"
  )
  
  cluster_counts_long[, PatientID := factor(PatientID, levels = patient_order)]
  
  cluster_labels <- setNames(as.character(clusterLabels$BackupTeXClusterLabel), clusterLabels$normalizedID)
  
  ggplot(cluster_counts_long, aes(x = PatientID, y = Count, fill = ClusterID)) +
    geom_bar(stat = "identity", width = 1.01) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.line.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(vjust = -1.5, hjust = -0.2, size = 13),
      plot.margin = unit(c(0.2, 0, 0, 0), "cm"),
      legend.position = "none"
    ) +
    labs(x = NULL, y = expression(italic(n) ~ "cells")) +
    scale_fill_manual(values = clusterLabels$Colour, labels = sapply(cluster_labels, TeX), name = NULL) +
    guides(fill = guide_legend(ncol = 2)) +
    scale_y_continuous(
      expand = c(0, 0),
      breaks = c(floor(max(cellData[, .N, by = PatientID]$N) / 1000) * 1000),
      labels = c(scales::comma(floor(max(cellData[, .N, by = PatientID]$N) / 1000) * 1000))
    )
}

makeProportionBarChart <- function(cellData, clusterLabels, patient_order, titleY) {
  cluster_counts <- getCellCounts(
    cellData,
    clusterColumn = "ClusterID",
    idCols = "PatientID"
  )[, .(PatientID, ClusterID, proportion)]
  
  proportion_matrix <- dcast(cluster_counts, PatientID ~ ClusterID, value.var = "proportion", fill = 0)
  proportion_matrix[, PatientID := factor(PatientID, levels = patient_order)]
  
  long_format <- melt(proportion_matrix, id.vars = "PatientID", variable.name = "Cluster", value.name = "Proportion")
  long_format[, Cluster := factor(Cluster, levels = as.character(clusterLabels$ClusterID[order(clusterLabels$PrintOrder)]))]
  
  ggplot(long_format, aes(x = PatientID, y = Proportion, fill = Cluster)) +
    geom_bar(stat = "identity", width = 1.01) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(vjust = -10, hjust = 0.5, size = 13),
      plot.margin = unit(c(0, 0, 0.55, 0), "cm"),
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.45, "cm"),
      legend.position = "right",
      legend.justification = "left",
      legend.box.just = "left",
      legend.margin = unit(c(0, 5, 0, 5), "pt")
    ) +
    labs(x = NULL, y = titleY) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 1), labels = c("0", "1")) +
    scale_fill_manual(
      values = clusterLabels$Colour,
      labels = sapply(setNames(as.character(clusterLabels$BackupTeXClusterLabel), clusterLabels$ClusterID), TeX),
      name = NULL
    ) +
    guides(fill = guide_legend(ncol = 1))
}

# ---- assemble ----
epithelial_count_bar_chart <- makeCountsBarChart(epithelial, epiannotations, patient_order)
microenvironment_count_bar_chart <- makeCountsBarChart(microenvironment, tmeannotations, patient_order)

epithelial_proportion_bar_chart <- makeProportionBarChart(epithelial, epiannotations, patient_order, "Epithelial\nproportion")
microenvironment_proportion_bar_chart <- makeProportionBarChart(microenvironment, tmeannotations, patient_order, "Microenvironment\nproportion")

combined_plot <- age_scatter_plot /
  epithelial_count_bar_chart /
  epithelial_proportion_bar_chart /
  microenvironment_count_bar_chart /
  microenvironment_proportion_bar_chart +
  plot_layout(ncol = 1, heights = c(1, 0.9, 2.15, 0.9, 2.15)) +
  plot_annotation(theme = theme(plot.margin = unit(c(0, 0, 30, 0), "pt")))

ggsave(
  filename = file.path(outdir, "patientCompositionClustering.pdf"),
  plot = combined_plot, width = 8.5, height = 7, units = "in"
)
