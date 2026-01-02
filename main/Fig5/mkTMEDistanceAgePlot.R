library(here)
source(here("code", "header.R"))

outdir <- here("scratch/main/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()[, .(ImageID, PatientID, Age)]
annotations <- getCellClusters()

annotations[, ClusterID := as.integer(ClusterID)]
tmeannotations <- annotations[isEpithelial == FALSE]
tmeannotations[, normalizedID := .GRP, by = ClusterID]

if (!("X" %in% names(cells))) setnames(cells, "CenterX", "X", skip_absent = TRUE)
if (!("Y" %in% names(cells))) setnames(cells, "CenterY", "Y", skip_absent = TRUE)

cellsclin <- merge(
  cells[, .(ImageID, CellID, ClusterID, isEpithelial, X, Y)],
  clinical[, .(ImageID, PatientID, Age)],
  by = "ImageID",
  all.x = TRUE
)

cellsclin <- cellsclin[!is.na(X) & !is.na(Y) & !is.na(isEpithelial)]

epi_cells <- cellsclin[isEpithelial == TRUE,  .(ImageID, CellID, X, Y)]
tme_cells <- cellsclin[isEpithelial == FALSE, .(ClusterID, PatientID, ImageID, CellID, X, Y)]

unique_images <- unique(tme_cells$ImageID)
cat("Computing nearest-epithelium distances for", length(unique_images), "images\n")

mkDataPerImg <- function(imgID) {
  dat <- cellsclin[ImageID == imgID, .(ImageID, CellID, isEpithelial, X, Y)]
  if (nrow(dat[isEpithelial == TRUE]) == 0L) return(NULL)
  
  dat[, EpiTME := fifelse(isEpithelial, "Epi", "TME")]
  setnames(dat, c("X","Y"), c("Location_Center_X","Location_Center_Y"))
  
  id_map <- dat[, .(ImageID, CellID, isEpithelial, EpiTME)]
  d <- as.data.table(DistanceToNearestNeighbour(dat, "EpiTME"))
  stopifnot(nrow(d) == nrow(id_map))
  
  d <- cbind(id_map, d)
  d[EpiTME == "TME", .(ImageID, CellID, shortest_distance = DistToEpi)]
}

n_cores <- max(1, parallel::detectCores() - 1)
combined <- rbindlist(parallel::mclapply(unique_images, mkDataPerImg, mc.cores = n_cores), fill = TRUE)

tme_cells <- merge(tme_cells, combined, by = c("ImageID", "CellID"), all.x = TRUE)
cat("Done. Missing distance rows:", sum(is.na(tme_cells$shortest_distance)), "\n")

req_cols <- c("ImageID", "ClusterID", "shortest_distance")
miss <- setdiff(req_cols, names(tme_cells))
if (length(miss) > 0) stop("tme_cells missing columns: ", paste(miss, collapse = ", "))

tme_cells_clin <- merge(
  tme_cells,
  clinical[, setdiff(names(clinical), "PatientID"), with = FALSE],
  by = "ImageID",
  all.x = TRUE
)

median_shortest_distance <- tme_cells_clin[
  !is.na(PatientID) & !is.na(Age) & !is.na(shortest_distance) & !is.na(ClusterID),
  .(MedianDistance = median(shortest_distance, na.rm = TRUE)),
  by = .(PatientID, ClusterID, Age)
]
median_shortest_distance[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]

patients_per_cluster <- median_shortest_distance[, .(PatientCount = uniqueN(PatientID)), by = ClusterID]
setnames(patients_per_cluster, "ClusterID", "Cluster")

summary_results <- data.table()

for (cluster in sort(unique(median_shortest_distance$ClusterID))) {
  cluster_data <- median_shortest_distance[ClusterID == cluster]
  cluster_data <- cluster_data[!is.na(MedianDistance) & !is.na(AgeGroup)]
  if (uniqueN(cluster_data$AgeGroup) < 2) next
  
  test_result <- t.test(MedianDistance ~ AgeGroup, data = cluster_data, var.equal = FALSE)
  cluster_data[, AgeGroup := factor(AgeGroup, levels = c("Above 50", "Below 50"))]
  cd <- cliff.delta(MedianDistance ~ AgeGroup, data = cluster_data)
  
  temp_results <- data.table(
    Cluster = as.integer(cluster),
    TStatistic = unname(test_result$statistic),
    PValue = test_result$p.value,
    CliffsDelta = unname(cd$estimate),
    CID_Lower = unname(cd$conf.int[1]),
    CID_Upper = unname(cd$conf.int[2]),
    Mean_Above50 = mean(cluster_data[AgeGroup == "Above 50"]$MedianDistance, na.rm = TRUE),
    Mean_Below50 = mean(cluster_data[AgeGroup == "Below 50"]$MedianDistance, na.rm = TRUE),
    patients_premenopause = sum(cluster_data$AgeGroup == "Below 50", na.rm = TRUE),
    patients_postmenopause = sum(cluster_data$AgeGroup == "Above 50", na.rm = TRUE)
  )
  
  summary_results <- rbind(summary_results, temp_results, fill = TRUE)
}

if (nrow(summary_results) == 0) stop("No clusters had both age groups after filtering.")

summary_results[, adjPValue := p.adjust(PValue, method = "BH")]
summary_results[, alpha := ifelse(adjPValue < 0.05, "opaque", "translucent")]

patients_per_cluster[, Cluster := as.integer(Cluster)]
summary_results <- merge(summary_results, patients_per_cluster, by = "Cluster", all.x = TRUE)

summary_results <- summary_results[order(CliffsDelta, decreasing = TRUE)]
cluster_order <- summary_results$Cluster  # integer vector in desired order
summary_results[, Cluster := factor(Cluster, levels = cluster_order)]

summary_results[, point_size := 5]
summary_results[, adjPValueLabel := sapply(format_custom_pval(adjPValue), mkEnumPower)]

tmeannotations[, BackupTeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
labels <- setNames(tmeannotations$BackupTeXClusterLabel, as.character(tmeannotations$ClusterID))

p <- ggplot(summary_results, aes(x = Cluster, y = CliffsDelta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
  geom_errorbar(aes(ymin = CID_Lower, ymax = CID_Upper, alpha = alpha), width = 0.2) +
  geom_point(shape = 21, aes(fill = alpha, size = point_size)) +
  scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
  coord_flip(clip = "off", ylim = c(min(summary_results$CID_Lower), max(summary_results$CID_Upper))) +
  labs(x = NULL, y = "Cliff's delta effect size", title = NULL) +
  theme_classic() +
  scale_y_continuous(breaks = c(-0.6, -0.15, 0, 0.3), labels = function(x) ifelse(x == 0, "0", format(x))) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.title.x = element_text(colour = "black", size = 11),
    plot.margin = unit(c(30, 40, 5, 0), "pt"),
    legend.position = "bottom"
  ) +
  guides(color = "none", alpha = "none", fill = "none") +
  scale_size_identity() +
  geom_text(
    aes(label = TeX(adjPValueLabel, output = "character"), color = alpha, y = 0.45, x = Cluster),
    parse = TRUE, hjust = 0, size = 3.75
  ) +
  scale_color_manual(values = c("opaque" = "black", "translucent" = "darkgray")) +
  scale_x_discrete(labels = sapply(labels, TeX))

createColorBarPlot <- function(summaryData, clusterLabels, cluster_order) {

  summaryData[, ClusterJoin := as.integer(as.character(Cluster))]
  
  clusterLabels[, ClusterID := as.integer(ClusterID)]
  
  clusterLabels <- merge(
    clusterLabels,
    summaryData[, .(ClusterJoin, CliffsDelta)],
    by.x = "ClusterID",
    by.y = "ClusterJoin",
    all.x = FALSE
  )
  
  clusterLabels[, ClusterID := factor(ClusterID, levels = cluster_order)]
  clusterLabels <- clusterLabels[order(ClusterID)]
  
  clusterLabels[, TeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
  clusternames <- setNames(clusterLabels$TeXClusterLabel, clusterLabels$normalizedID)
  colors <- setNames(clusterLabels$Colour, clusterLabels$normalizedID)
  clusterLabels[, normalizedID := factor(normalizedID, levels = names(colors))]
  
  ggplot(clusterLabels, aes(x = factor(1), y = normalizedID, fill = normalizedID)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_y_discrete(
      limits = levels(clusterLabels$normalizedID),
      labels = sapply(clusternames[levels(clusterLabels$normalizedID)], TeX)
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(y = "Distance to epithelium") +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.title.y = element_text(color = "black", size = 11),
      axis.text.y = element_text(color = "black", size = 11),
      axis.ticks.length.y = unit(0, "pt"),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), units = "cm")
    )
}

summary_results_long <- melt(
  summary_results,
  id.vars = "Cluster",
  measure.vars = c("patients_premenopause", "patients_postmenopause"),
  variable.name = "Group",
  value.name = "Count"
)
summary_results_long[, Cluster := factor(Cluster, levels = cluster_order)]
summary_results_long[, Group := factor(Group, levels = c("patients_postmenopause", "patients_premenopause"))]
levels(summary_results_long$Group) <- c("Postmenopause", "Premenopause")

bar_plot <- ggplot(summary_results_long, aes(x = Cluster, y = Count, fill = Group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Premenopause" = "gray", "Postmenopause" = "dimgray")) +
  coord_flip() +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 500)) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(colour = "black", size = 11),
    axis.title.x = element_text(colour = "black", size = 11),
    plot.margin = unit(c(30, 20, 5, 0), "pt"),
    legend.position = "none"
  )
summary_results
annotations
plotcombined <- createColorBarPlot(summary_results, tmeannotations, cluster_order) + p + bar_plot +
  plot_layout(widths = c(0.3, 6, 2))

pdf(file.path(outdir, "TMEDistanceAgePlot.pdf"), width = 4, height = 3.6)
print(plotcombined)
dev.off()

message("Wrote: ", file.path(outdir, "TMEDistanceAgePlot.pdf"))

# #Export result table
# summary_results_typed <- merge(
#   copy(summary_results)[, Cluster := as.character(Cluster)],
#   tmeannotations[, .(ClusterID, BackupLabel, Type)][, ClusterID := as.character(ClusterID)],
#   by.x = "Cluster",
#   by.y = "ClusterID",
#   all.x = TRUE
# )
# 
# n_table <- summary_results_typed[, .(
#   ClusterID              = as.integer(as.character(Cluster)),
#   PhenotypeLabel_Fallback = BackupLabel,
#   Compartment            = Type,
#   n_Premenopause         = patients_premenopause,
#   n_Postmenopause        = patients_postmenopause,
#   p_value                = PValue,
#   p_adj_BH               = adjPValue,
#   cliffs_delta           = CliffsDelta,
#   CI_lower               = CID_Lower,
#   CI_upper               = CID_Upper,
#   mean_premenopause      = Mean_Below50,
#   mean_postmenopause     = Mean_Above50
# )]
# 
# n_table <- n_table[order(match(ClusterID, cluster_order))]
# 
# fwrite(
#   n_table,
#   file = file.path(outdir, "Fig5c_TMEDistance_AgeEffectSizes_byPhenotype.csv")
# )
