#Requires the output of doBivariateRipleyKCross.R

library(here)
source(here("code/header.R"))

outdir <- here("scratch", "main", "Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
clinical <- getClinical()
tmeann   <- getCellClusters()

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

# distance per image (parallel)
mkDataPerImg <- function(imgID) {
  dat <- cells[ImageID == imgID, .(ImageID, CellID, isEpithelial, X, Y)]
  if (nrow(dat[isEpithelial == TRUE]) == 0) return(NULL)
  
  dat[, EpiTME := ifelse(isEpithelial, "Epi", "TME")]
  
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

ripley_fp <- file.path(outdir, "bivariateRipleyKAllResultsDT.rds")
if (!file.exists(ripley_fp)) {
  stop("Missing Ripley results: ", ripley_fp, "\nRun Script B first (Ripley Kcross).")
}
results <- readRDS(ripley_fp)
results <- as.data.table(results)

proportion_dt <- results[, .(
  ProportionOutOfRangeAbove = sum(OutOfRangeAbove == TRUE, na.rm = TRUE) / .N,
  ProportionOutOfRangeBelow = sum(OutOfRangeBelow == TRUE, na.rm = TRUE) / .N,
  images = sum(!is.na(OutOfRangeAbove) & !is.na(OutOfRangeBelow))
), by = ClusterID]

cluster_medians <- tme_cells[, .(median_distance = median(shortest_distance, na.rm = TRUE)), by = ClusterID]
setorder(cluster_medians, -median_distance)
cluster_medians[, ClusterID := factor(ClusterID, levels = unique(ClusterID))]

tmeann <- tmeann[isEpithelial == FALSE]
tmeann[, ClusterID := as.factor(ClusterID)]
tme_cells[, ClusterID := factor(ClusterID, levels = levels(cluster_medians$ClusterID))]
proportion_dt[, ClusterID := factor(ClusterID, levels = levels(tme_cells$ClusterID))]

if (!("Colour" %in% names(tmeann))) stop("Cluster annotations missing 'Colour'.")
color_mapping <- setNames(tmeann$Colour, tmeann$ClusterID)

if (!("BackupTeXClusterLabel" %in% names(tmeann))) stop("Cluster annotations missing 'BackupTeXClusterLabel'.")
labels <- setNames(tmeann$BackupTeXClusterLabel, tmeann$ClusterID)

sampled_tme_cells <- tme_cells[!is.na(shortest_distance)][
  , .SD[sample(.N, min(.N, 5000))], by = ClusterID
]

p_main <- ggplot(sampled_tme_cells, aes(x = ClusterID, y = shortest_distance, fill = ClusterID)) +
  geom_point(position = position_jitter(width = 0.45, height = 0), shape = 16, size = 0.1, alpha = 0.2) +
  geom_boxplot(alpha = 0.75, width = 0.6, outlier.shape = NA) +
  coord_flip(ylim = c(0, 800)) +
  labs(x = NULL, y = "nearest distance to epithelium (µm)") +
  scale_fill_manual(values = color_mapping) +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 400, 800)) +
  theme(
    axis.title.y = element_text(colour = "black"),
    axis.title.x = element_text(size = 13),
    axis.text.y = element_text(colour = "black", size = 8.5),
    axis.text.x = element_text(colour = "black", size = 8.5),
    plot.margin = unit(c(20, 20, 5, 20), "pt"),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = sapply(labels, TeX))

p_text_1 <- ggplot(proportion_dt, aes(x = ClusterID, y = 0, fill = ProportionOutOfRangeAbove)) +
  geom_tile(aes(height = 1, width = 1), color = "white") +
  scale_fill_gradient(low = "black", high = "white", limits = c(0, 1)) +
  geom_text(aes(label = sprintf("%.0f", ProportionOutOfRangeAbove * 100),
                color = ProportionOutOfRangeAbove > 0.5),
            vjust = 0.5, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("white", "black")) +
  theme_void() +
  theme(plot.margin = unit(c(0, 0, 0, 30), "pt"),
        legend.position = "none") +
  coord_flip(clip = "off", ylim = c(0, 1)) +
  scale_x_discrete(limits = levels(tme_cells$ClusterID))

p_text_2 <- ggplot(proportion_dt, aes(x = ClusterID, y = 0, fill = ProportionOutOfRangeBelow)) +
  geom_tile(aes(height = 1, width = 1), color = "white") +
  scale_fill_gradient(low = "black", high = "white", limits = c(0, 1)) +
  geom_text(aes(label = sprintf("%.0f", ProportionOutOfRangeBelow * 100),
                color = ProportionOutOfRangeBelow > 0.5),
            vjust = 0.5, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("white", "black")) +
  theme_void() +
  theme(plot.margin = unit(c(0, 5, 0, 0), "pt"),
        legend.position = "none") +
  coord_flip(clip = "off", ylim = c(0, 1)) +
  scale_x_discrete(limits = levels(tme_cells$ClusterID))

p_text_3 <- ggplot(proportion_dt, aes(x = ClusterID, y = images)) +
  geom_bar(stat = "identity", fill = "dimgray") +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 1300), labels = c("0", "1250")) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 8.5),
    plot.margin = unit(c(0, 30, 0, 10), "pt")
  )

cluster_abundances <- tme_cells[, .N, by = ClusterID]
cluster_abundances[, ClusterID := factor(ClusterID, levels = levels(tme_cells$ClusterID))]

p_text_0 <- ggplot(cluster_abundances, aes(x = ClusterID, y = N)) +
  geom_bar(stat = "identity", fill = "dimgray") +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(breaks = c(0, 500000), labels = c("0", "0.5M")) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 8.5),
    plot.margin = unit(c(0, 0, 0, 10), "pt")
  )

p_combined <- p_main + p_text_0 + p_text_1 + p_text_2 + p_text_3 +
  plot_layout(widths = c(11.5, 2.5, 1, 1, 2.5))

out_pdf <- file.path(outdir, "TMEDistancesFromEpithelium.pdf")
ggsave(out_pdf, plot = p_combined, width = 7, height = 3.2, units = "in")
message("Wrote: ", out_pdf)
