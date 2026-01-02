library(here)
source(here("code", "header.R"))

outdir <- here("scratch/main/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()[, .(ImageID, PatientID, Age)]
annotations <- getCellClusters()

annotations[, ClusterID := as.integer(ClusterID)]
tmeannotations <- annotations[isEpithelial == FALSE]
tmeannotations[, TeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]

if (!("X" %in% names(cells))) setnames(cells, "CenterX", "X", skip_absent = TRUE)
if (!("Y" %in% names(cells))) setnames(cells, "CenterY", "Y", skip_absent = TRUE)

cellsclin <- merge(
  cells[, .(ImageID, CellID, ClusterID, isEpithelial, X, Y)],
  clinical[, .(ImageID, PatientID, Age)],
  by = "ImageID",
  all.x = TRUE
)
cellsclin <- cellsclin[!is.na(X) & !is.na(Y) & !is.na(isEpithelial)]
cellsclin[, ClusterID := as.integer(ClusterID)]

tme_cells <- cellsclin[isEpithelial == FALSE, .(ClusterID, PatientID, ImageID, CellID, X, Y)]

unique_images <- unique(tme_cells$ImageID)
cat("Computing nearest-epithelium distances for", length(unique_images), "images\n")

mkDataPerImg <- function(imgID) {
  dat <- cellsclin[ImageID == imgID, .(ImageID, CellID, isEpithelial, X, Y)]
  if (nrow(dat[isEpithelial == TRUE]) == 0L) return(NULL)
  
  dat[, EpiTME := fifelse(isEpithelial, "Epi", "TME")]
  setnames(dat, c("X", "Y"), c("Location_Center_X", "Location_Center_Y"))
  
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

tme_cells <- merge(
  tme_cells,
  clinical[, setdiff(names(clinical), "PatientID"), with = FALSE],
  by = "ImageID",
  all.x = TRUE
)

median_shortest_distance <- tme_cells[
  !is.na(PatientID) & !is.na(Age) & !is.na(shortest_distance) & !is.na(ClusterID),
  .(MedianDistance = median(shortest_distance, na.rm = TRUE)),
  by = .(PatientID, ClusterID, Age)
]
median_shortest_distance[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]

summary_results <- data.table()
for (cluster in sort(unique(median_shortest_distance$ClusterID))) {
  cluster_data <- median_shortest_distance[ClusterID == cluster]
  cluster_data <- cluster_data[!is.na(MedianDistance) & !is.na(AgeGroup)]
  if (uniqueN(cluster_data$AgeGroup) < 2) next
  
  test_result <- t.test(MedianDistance ~ AgeGroup, data = cluster_data, var.equal = FALSE)
  summary_results <- rbind(
    summary_results,
    data.table(ClusterID = as.integer(cluster), PValue = test_result$p.value),
    fill = TRUE
  )
}
summary_results[, adjPValue := p.adjust(PValue, method = "BH")]

createDistanceBoxplot <- function(tmeData, summaryData, cluster_id, ymax, tmeann) {
  cluster_id <- as.integer(cluster_id)
  
  cluster_specific_data <- tmeData[ClusterID == cluster_id]
  cluster_specific_data <- cluster_specific_data[!is.na(AgeGroup) & !is.na(MedianDistance)]
  cluster_specific_data[, AgeGroup := factor(AgeGroup, levels = c("Below 50", "Above 50"))]
  
  adjPVal <- summaryData[ClusterID == cluster_id, adjPValue]
  if (length(adjPVal) == 0) adjPVal <- NA_real_
  
  label <- tmeann[ClusterID == cluster_id, TeXClusterLabel]
  if (length(label) == 0) label <- mkTeX(as.character(cluster_id))
  
  ggplot(cluster_specific_data, aes(x = AgeGroup, y = MedianDistance, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("Below 50" = "gray", "Above 50" = "dimgray")) +
    theme_classic(base_family = "ArialMT") +
    scale_x_discrete(labels = c("Below 50" = "pre", "Above 50" = "post")) +
    scale_y_continuous(
      limits = c(0, ymax),
      expand = expand_scale(mult = c(0, 0.1)),
      breaks = c(0, (floor(ymax / 100) * 100) / 2, floor(ymax / 100) * 100),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title = TeX(label)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.title.y = element_text(vjust = -3, hjust = 0.5, size = 18),
      axis.text.y = element_text(size = 17),
      plot.title = element_text(size = 18, hjust = 0.5),
      plot.margin = unit(c(10, 0, 10, 0), units = "pt"),
      legend.position = "none"
    ) +
    annotate(
      "text", x = 1.5, y = ymax,
      label = TeX(paste0("p=", mkEnumPower(ifelse(
        is.na(adjPVal), "NA",
        ifelse(adjPVal < 0.001, formatC(adjPVal, format = "e", digits = 0), sprintf("%.3f", adjPVal))
      )))),
      size = 5.5, color = "black", hjust = 0.5
    )
}


distancePlot1 <- createDistanceBoxplot(median_shortest_distance, summary_results, 21, ymax = 650, tmeann = tmeannotations)
distancePlot2 <- createDistanceBoxplot(median_shortest_distance, summary_results, 19, ymax = 600, tmeann = tmeannotations)
distancePlot3 <- createDistanceBoxplot(median_shortest_distance, summary_results, 25, ymax = 610, tmeann = tmeannotations)
distancePlot4 <- createDistanceBoxplot(median_shortest_distance, summary_results, 14, ymax = 840, tmeann = tmeannotations)

distancecombined <- distancePlot1 + distancePlot2 + distancePlot3 + distancePlot4 + plot_layout(ncol = 2)

distanceblanklabelploty <- ggplot() +
  labs(y = "median distance to\nepithelium (µm)") +
  theme_classic() +
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "cm"),
    axis.title.y = element_text(color = "black", size = 17)
  ) +
  guides(x = "none", y = "none")

distancecombined <- distanceblanklabelploty + distancecombined + plot_layout(widths = c(1, 1000))

pdf(file.path(outdir, "TMEDistanceAgeBoxplots.pdf"), width = 4, height = 4.5)
print(distancecombined)
dev.off()

message("Wrote: ", file.path(outdir, "TMEDistanceAgeBoxplots.pdf"))
