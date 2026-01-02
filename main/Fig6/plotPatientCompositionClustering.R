library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
clinical <- getClinical()

cells[, SpatialClusterID := as.character(SpatialClusterID)]
clinical[, PatientID := as.character(PatientID)]

cellsclin <- merge(
  cells,
  clinical[, .(ImageID, PatientID, Age)],
  by = "ImageID",
  all.x = TRUE
)

comm0 <- fread(here("data", "derived", "spatialClusterAnnotation.csv"))
comm0[, SpatialClusterID := as.character(SpatialClusterID)]

comm_ord <- comm0[order(PrintOrder)]

spatial_levels <- comm_ord$SpatialClusterID
color_map      <- setNames(comm_ord$Colour, comm_ord$SpatialClusterID)
cluster_labels <- setNames(as.character(comm_ord$Label), comm_ord$SpatialClusterID)

cluster_counts <- cellsclin[, .(CellCount = .N), by = .(PatientID, SpatialClusterID)]
total_counts   <- cellsclin[, .(TotalCount = .N), by = PatientID]
merged_counts  <- merge(cluster_counts, total_counts, by = "PatientID")
merged_counts[, ClusterProportion := CellCount / TotalCount]

merged_counts[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]

proportion_matrix <- dcast(
  merged_counts,
  PatientID ~ SpatialClusterID,
  value.var = "ClusterProportion",
  fill = 0
)

mat <- as.matrix(proportion_matrix[, -1, with = FALSE])

hc <- if (exists("getHC")) getHC(mat) else hclust(dist(mat))
patient_order <- proportion_matrix$PatientID[hc$order]

age_data <- unique(cellsclin[, .(PatientID, Age)])
age_data[, PatientID := factor(PatientID, levels = patient_order)]
age_data[, AgeGroup := ifelse(Age < 50, "Below 50", "Above 50")]
age_colors <- c("Below 50" = "gray", "Above 50" = "dimgray")

age_scatter_plot <- ggplot(age_data, aes(x = PatientID, y = Age, fill = AgeGroup, color = AgeGroup)) +
  geom_point(size = 0.3, shape = 21) +
  scale_color_manual(values = age_colors) +
  scale_fill_manual(values = age_colors) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title.y = element_text(vjust = -1.2, size = 11),
    plot.margin = unit(c(5, 2, 0, 0), "pt"),
    legend.position = "none"
  ) +
  labs(x = NULL, y = "age") +
  scale_y_continuous(breaks = c(20, 70), expand = c(0.02, 0.02))

makeCountsBarChart <- function(cellData, patientOrder, spatial_levels, color_map, cluster_labels) {
  cluster_counts <- cellData[, .(CellCount = .N), by = .(PatientID, SpatialClusterID)]
  total_counts   <- cellData[, .(TotalCells = .N), by = PatientID]
  
  max_total_count <- max(total_counts$TotalCells, na.rm = TRUE)
  rounded_max_count <- floor(max_total_count / 1000) * 1000
  
  cluster_counts[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
  cluster_counts[, PatientID := factor(PatientID, levels = patientOrder)]
  
  ggplot(cluster_counts, aes(x = PatientID, y = CellCount, fill = SpatialClusterID)) +
    geom_bar(stat = "identity", width = 1.3) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(vjust = -1.2, hjust = 0.4, size = 11),
      plot.margin = unit(c(0.2, 0, 0, 0), "cm"),
      legend.position = "none"
    ) +
    labs(x = NULL, y = expression(italic(n) ~ "cells")) +
    scale_fill_manual(values = color_map, labels = sapply(cluster_labels, TeX), name = NULL) +
    scale_y_continuous(expand = c(0, 0), breaks = c(rounded_max_count),
                       labels = c(scales::comma(rounded_max_count)))
}

spatial_count_bar_chart <- makeCountsBarChart(
  cellData = cellsclin,
  patientOrder = patient_order,
  spatial_levels = spatial_levels,
  color_map = color_map,
  cluster_labels = cluster_labels
)

makeProportionBarChart <- function(cellData, patientOrder, spatial_levels, color_map, cluster_labels, titleY) {
  cluster_counts <- cellData[, .(CellCount = .N), by = .(PatientID, SpatialClusterID)]
  total_counts   <- cellData[, .(TotalCells = .N), by = PatientID]
  merged_counts  <- merge(cluster_counts, total_counts, by = "PatientID")
  merged_counts[, ClusterProportion := CellCount / TotalCells]
  
  merged_counts[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
  
  proportion_matrix <- dcast(
    merged_counts,
    PatientID ~ SpatialClusterID,
    value.var = "ClusterProportion",
    fill = 0
  )
  proportion_matrix[, PatientID := factor(PatientID, levels = patientOrder)]
  
  long_format <- melt(
    proportion_matrix,
    id.vars = "PatientID",
    variable.name = "SpatialClusterID",
    value.name = "Proportion"
  )
  long_format[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
  
  ggplot(long_format, aes(x = PatientID, y = Proportion, fill = SpatialClusterID)) +
    geom_bar(stat = "identity", width = 1.3) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 10),
      axis.title.y = element_text(vjust = -10, size = 11),
      plot.margin = unit(c(0, 0, 0.2, 0), "cm"),
      legend.text = element_text(size = 6.5),
      legend.key.size = unit(0.35, "cm"),
      legend.position = "right",
      legend.justification = "left",
      legend.margin = unit(c(0,5,0,0), "pt"),
      legend.title = element_blank()
    ) +
    labs(x = NULL, y = titleY) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 1), labels = c("0", "1")) +
    scale_fill_manual(values = color_map, labels = sapply(cluster_labels, TeX)) +
    guides(fill = guide_legend(ncol = 1))
}

spatial_proportion_bar_chart <- makeProportionBarChart(
  cellData = cellsclin,
  patientOrder = patient_order,
  spatial_levels = spatial_levels,
  color_map = color_map,
  cluster_labels = cluster_labels,
  titleY = "Spatial cluster\nproportion"
)

combined_plot <- age_scatter_plot / spatial_count_bar_chart / spatial_proportion_bar_chart +
  plot_layout(ncol = 1, heights = c(1, 1, 2.5))

ggsave(
  filename = file.path(outdir, "PatientCompositionSpatialClustering.pdf"),
  plot = combined_plot,
  width = 6.5, height = 3, units = "in"
)
