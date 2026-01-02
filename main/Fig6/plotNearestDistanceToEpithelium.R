library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()

communityannotations <- fread(here("data", "derived", "spatialClusterAnnotation.csv"))
communityannotations <- as.data.table(communityannotations)
communityannotations[, SpatialClusterID := as.character(SpatialClusterID)]

cells[, SpatialClusterID := as.character(SpatialClusterID)]
cells[, CellID := as.character(CellID)]

nearest_dist_to_epi_one_image <- function(imgID, cells_dt) {
  dat <- cells_dt[ImageID == imgID]
  if (nrow(dat) == 0) return(NULL)
  
  epi <- dat[isEpithelial == TRUE,  .(x = CenterX, y = CenterY)]
  tme <- dat[isEpithelial == FALSE, .(ImageID, CellID, SpatialClusterID, x = CenterX, y = CenterY)]
  
  if (nrow(epi) == 0 || nrow(tme) == 0) return(NULL)
  
  ex <- epi$x
  ey <- epi$y
  
  min_d2 <- vapply(seq_len(nrow(tme)), function(i) {
    dx <- ex - tme$x[i]
    dy <- ey - tme$y[i]
    min(dx*dx + dy*dy, na.rm = TRUE)
  }, numeric(1))
  
  tme[, shortest_distance := sqrt(min_d2)]
  tme[, .(ImageID, CellID, SpatialClusterID, shortest_distance)]
}

unique_images <- unique(cells$ImageID)
cat("Starting distance calc for", length(unique_images), "images.\n")

cores <- max(1, detectCores() - 1)
dist_list <- mclapply(unique_images, nearest_dist_to_epi_one_image, cells_dt = cells, mc.cores = cores)
distance_cells <- rbindlist(dist_list, use.names = TRUE, fill = TRUE)

cat("All distances calculated. Rows:", nrow(distance_cells), "\n")

cluster_medians <- distance_cells[, .(median_distance = median(shortest_distance, na.rm = TRUE)),
                                  by = .(SpatialClusterID)]
setorder(cluster_medians, -median_distance)

distance_cells[, SpatialClusterID := factor(SpatialClusterID, levels = cluster_medians$SpatialClusterID)]
communityannotations[, SpatialClusterID := factor(SpatialClusterID, levels = cluster_medians$SpatialClusterID)]

color_mapping <- setNames(as.character(communityannotations$Colour), as.character(communityannotations$SpatialClusterID))
labels_map    <- setNames(as.character(communityannotations$Label),  as.character(communityannotations$SpatialClusterID))

desired_clusters <- c("3", "2", "6", "8", "7")
distance_cells <- distance_cells[as.character(SpatialClusterID) %in% desired_clusters]
distance_cells[, SpatialClusterID := factor(as.character(SpatialClusterID), levels = desired_clusters)]

set.seed(1)
sampled_tme_cells <- distance_cells[, .SD[sample(.N, min(.N, 10000))], by = SpatialClusterID]

p_main <- ggplot(sampled_tme_cells, aes(x = SpatialClusterID, y = shortest_distance, fill = SpatialClusterID)) +
  geom_point(position = position_jitter(width = 0.45, height = 0),
             shape = 16, size = 0.1, alpha = 0.2) +
  geom_boxplot(alpha = 0.75, width = 0.6, outlier.shape = NA) +
  coord_flip(ylim = c(0, 800)) +
  labs(x = NULL, y = "shortest distance to epithelium (µm)") +
  scale_fill_manual(values = color_mapping) +
  scale_x_discrete(labels = sapply(labels_map, TeX)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 400, 800)) +
  theme_classic() +
  theme(
    axis.title.y = element_text(colour = "black", size = 17),
    axis.title.x = element_text(size = 14),
    axis.text.y  = element_text(colour = "black", size = 16),
    axis.text.x  = element_text(colour = "black", size = 16),
    plot.margin  = unit(c(20, 50, 5, 20), "pt"),
    legend.position = "none"
  )

ggsave(file.path(outdir, "nearestDistanceToEpithelium.pdf"),
       p_main, width = 5, height = 4)
