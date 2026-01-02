library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig9")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
annotations <- getCellClusters()
dist <- read_parquet_adt(here("data/derived/cellCompartmentDistances.parquet"))

clin_min <- unique(clinical[, .(ImageID, PatientID, Age)], by = "ImageID")
cellsclin <- merge(
  cells,
  clin_min,
  by = "ImageID",
  all.x = TRUE
)

dist_keep <- intersect(
  c("ImageID", "CellID", "distanceToAdipose", "distanceToVascular", "distanceToLymphatic"),
  names(dist)
)
cellsclin <- merge(
  cellsclin,
  dist[, ..dist_keep],
  by = c("ImageID", "CellID"),
  all.x = TRUE
)

collapse_epithelial <- function(dt, cluster_col = "ClusterID") {
  dt <- copy(dt)
  suppressWarnings(dt[, (cluster_col) := as.integer(get(cluster_col))])
  dt[get(cluster_col) %in% 1:11, (cluster_col) := 26L]
  dt[, (cluster_col) := factor(get(cluster_col))]
  dt[]
}

make_annotation_subset <- function(annotations) {
  ann <- copy(annotations)[, .(ClusterID, Colour, BackupTeXClusterLabel)]
  suppressWarnings(ann[, ClusterID := as.integer(ClusterID)])
  ann <- ann[!ClusterID %in% 1:11]
  
  new_row <- data.table(
    ClusterID = 26L,
    Colour = unname(getNormalBreastProjectColours()$BroadPhenotype["Epithelial"]),
    BackupTeXClusterLabel = "epithelial"
  )
  ann <- unique(rbind(ann, new_row, fill = TRUE), by = "ClusterID")
  
  ann[, ClusterID := as.character(ClusterID)]
  ann[]
}

annotation_subset <- make_annotation_subset(annotations)

plot_distance_combo <- function(
    dt_in,
    distance_col,
    out_pdf,
    y_limits,
    y_breaks,
    y_label = "",
    bar_min = 5000
) {
  dt <- copy(dt_in)
  
  keep <- c("ImageID", "CellID", "PatientID", "Age", "ClusterID", distance_col)
  keep <- keep[keep %in% names(dt)]
  dt <- dt[, ..keep]

  dt <- dt[!is.na(get(distance_col))]
  if (nrow(dt) == 0) stop("No non-NA distances found for: ", distance_col)
  
  dt <- collapse_epithelial(dt, "ClusterID")
  dt[, distance_val := get(distance_col)]
  cluster_medians <- dt[, .(median_distance = median(distance_val, na.rm = TRUE)), by = ClusterID]
  
  dt[, ClusterID := factor(ClusterID, levels = rev(cluster_medians[order(median_distance), ClusterID]))]
  
  color_mapping <- setNames(annotation_subset$Colour, annotation_subset$ClusterID)
  label_mapping <- setNames(annotation_subset$BackupTeXClusterLabel, annotation_subset$ClusterID)
  
  mainplot <- ggplot(dt, aes(x = as.factor(ClusterID), y = distance_val, fill = as.factor(ClusterID))) +
    geom_boxplot(outlier.size = 0.1, outlier.alpha = 0) +
    coord_flip() +
    labs(title = "", x = "", y = "") +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "black", size = 18),
      axis.text.y = element_text(color = "black", size = 15),
      legend.position = "none",
      plot.margin = unit(c(0,0,0,0), units = "pt")
    ) +
    scale_y_continuous(limits = y_limits, expand = c(0, 0), breaks = y_breaks) +
    scale_fill_manual(values = color_mapping) +
    scale_x_discrete(labels = sapply(label_mapping, TeX))
  
  cell_counts <- dt[, .N, by = ClusterID]
  barplot <- ggplot(cell_counts, aes(x = as.factor(ClusterID), y = N, fill = as.factor(ClusterID))) +
    geom_bar(stat = "identity", fill = "gray") +
    coord_flip(ylim = c(bar_min, NA)) +
    labs(title = "", x = "", y = "") +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "black", size = 18),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0,10,0,0), units = "pt")
    ) +
    scale_x_discrete(labels = sapply(label_mapping, TeX)) +
    scale_y_log10(expand = c(0, 0), breaks = c(5000, 50000, 1000000), labels = c("5K", "50K", "1M"))
  
  combined <- mainplot + barplot + plot_layout(widths = c(7,3))
  ggsave(combined, filename = file.path(outdir, out_pdf), width = 6, height = 4)
  cat("Wrote: ", file.path(outdir, out_pdf), "\n", sep = "")
}

plot_distance_combo(
  dt_in       = cellsclin,
  distance_col = "distanceToAdipose",
  out_pdf      = "adiposeTMEdistances.pdf",
  y_limits     = c(0, 1200),
  y_breaks     = c(0, 500, 1000)
)

plot_distance_combo(
  dt_in       = cellsclin,
  distance_col = "distanceToVascular",
  out_pdf      = "vesselTMEdistances.pdf",
  y_limits     = c(0, 700),
  y_breaks     = c(0, 300, 600)
)

plot_distance_combo(
  dt_in       = cellsclin,
  distance_col = "distanceToLymphatic",
  out_pdf      = "lymphaticTMEdistances.pdf",
  y_limits     = c(0, 500),
  y_breaks     = c(0, 200, 400)
)
