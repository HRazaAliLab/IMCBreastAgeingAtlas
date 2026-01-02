library(here)
source(here("code", "header.R"))

#Dependency required
library(spatstat)

outdir <- here("scratch/main/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
panel <- getPanel()
clinical <- getClinical()[, .(ImageID, PatientID, Age, ImageWidth, ImageHeight)]
annotations <- getCellClusters()
tmeannotations <- annotations[isEpithelial == FALSE]
tmeannotations[, normalizedID := .GRP, by = ClusterID]

if (!("X" %in% names(cells))) setnames(cells, "CenterX", "X", skip_absent = TRUE)
if (!("Y" %in% names(cells))) setnames(cells, "CenterY", "Y", skip_absent = TRUE)

imgs <- unique(clinical[, .(ImageID, ImageWidth, ImageHeight)], by = "ImageID")

epithelial_cells <- cells[isEpithelial == TRUE,  .(ImageID, X, Y)]
tme_cells        <- cells[isEpithelial == FALSE, .(ClusterID, ImageID, X, Y)]

plotKcrossEnvelope <- function(imageID, clusterID, simulationnumber = 199) {
  img_row <- imgs[ImageID == imageID]
  if (nrow(img_row) != 1) stop("ImageID not found (or duplicated) in clinical/imgs: ", imageID)
  
  current_epithelial_cells <- epithelial_cells[ImageID == imageID, .(x = X, y = Y, type = "epithelial")]
  current_tme_cells <- tme_cells[ImageID == imageID & ClusterID == clusterID, .(x = X, y = Y, type = "TME")]
  
  if (nrow(current_epithelial_cells) == 0 || nrow(current_tme_cells) == 0) {
    return(ggplot() + theme_void() + ggtitle(paste("Missing epi/TME:", imageID, clusterID)))
  }
  
  boundary <- owin(
    xrange = c(0, img_row$ImageWidth),
    yrange = c(0, img_row$ImageHeight)
  )
  
  combined_cells <- rbind(current_epithelial_cells, current_tme_cells)
  combined_ppp <- ppp(
    x = combined_cells$x,
    y = combined_cells$y,
    marks = factor(combined_cells$type),
    window = boundary
  )
  
  mark_chr <- as.character(marks(combined_ppp))
  if (!("epithelial" %in% mark_chr) || !("TME" %in% mark_chr)) {
    return(ggplot() + theme_void() + ggtitle(paste("Missing marks:", imageID, clusterID)))
  }
  
  env <- envelope(
    combined_ppp,
    fun = Kcross,
    i = "epithelial",
    j = "TME",
    nsim = simulationnumber,
    correction = "border",
    alternative = "two.sided",
    maxnerr = simulationnumber * 2,
    savefuns = TRUE
  )
  
  setnames(env, c("r", "obs", "theo", "lo", "hi"),
           c("distance", "observed", "theoretical", "low", "high"))
  env_df <- as.data.frame(env)
  
  max_value <- max(env_df$observed, env_df$theoretical, na.rm = TRUE)
  magnitude <- 10^floor(log10(max_value))
  rounded_max <- ceiling(max_value / magnitude) * magnitude
  custom_breaks <- c(0, rounded_max)
  
  clusterLabel <- tmeannotations[ClusterID == clusterID, BackupTeXClusterLabel]
  if (length(clusterLabel) == 0) clusterLabel <- as.character(clusterID)
  
  ggplot(env_df, aes(x = distance)) +
    geom_ribbon(aes(ymin = low, ymax = high), fill = "lightblue", alpha = 0.5) +
    geom_line(aes(y = observed), color = "black") +
    geom_line(aes(y = theoretical), color = "red", linetype = "dashed") +
    labs(title = TeX(clusterLabel), x = NULL, y = NULL) +
    scale_x_continuous(expand = c(0, 0), breaks = c(0, 150, 300)) +
    scale_y_continuous(breaks = custom_breaks, expand = c(0, 0)) +
    theme_classic() +
    theme(
      axis.text.x = element_text(color = "black", size = 17),
      axis.text.y = element_blank(),
      axis.title.x = element_text(color = "black", size = 18),
      axis.ticks.y = element_blank(),
      axis.title.y = element_text(color = "black", size = 18),
      plot.title = element_text(hjust = 0.5, size = 19),
      plot.margin = unit(c(0, 15, 0, 0), "pt")
    )
}

# Example panels
image_example <- "NormalBreastImgA_R2C10"
plot1 <- plotKcrossEnvelope(imageID = image_example, clusterID = 20)
plot2 <- plotKcrossEnvelope(imageID = image_example, clusterID = 13)
plot3 <- plotKcrossEnvelope(imageID = image_example, clusterID = 18)
plot4 <- plotKcrossEnvelope(imageID = image_example, clusterID = 23)

combinedplotfinal <- (plot1 + plot2 + plot3 + plot4 + plot_layout(ncol = 4))

blanklabelploty <- ggplot() +
  labs(y = TeX("$K_{cross}$")) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    axis.title.y = element_text(color = "black", size = 22)
  ) +
  guides(x = "none", y = "none")

blanklabelplotx <- ggplot() +
  labs(x = "distance from epithelium (µm)") +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    axis.title.x = element_text(color = "black", size = 22)
  ) +
  guides(x = "none", y = "none")

combinedplotfinal <- combinedplotfinal / blanklabelplotx + plot_layout(ncol = 1, heights = c(1000, 1))
combinedplotfinal <- blanklabelploty + combinedplotfinal + plot_layout(ncol = 2, widths = c(1, 1000))

ggsave(
  filename = file.path(outdir, "KCrossFunctionExamples.pdf"),
  plot = combinedplotfinal,
  width = 10.5, height = 3, units = "in"
)

plotTMEPhenotype <- function(imageID, clusterID) {
  dat <- cells[ImageID == imageID]
  cellSub <- dat[, .(CellID, isEpithelial, ClusterID)]
  
  imgPaths <- getImagePaths(imageID)
  contours <- getContour(imgPaths$Cell)
  
  contours <- merge(contours, cellSub, by = "CellID", all.x = TRUE)
  
  specified_color <- tmeannotations[ClusterID == clusterID, Colour]
  if (length(specified_color) == 0) specified_color <- "yellow"
    
  p <- mkRGBplot(imageID, noRed = TRUE)
  p$plot +
    geom_polygon(
      data = contours,
      mapping = aes(x = x, y = y, group = CellID),
      colour = "darkgray", linewidth = 0.1
    ) +
    geom_polygon(
      data = contours[isEpithelial == TRUE],
      mapping = aes(x = x, y = y, group = CellID),
      colour = "darkgray", fill = "darkgreen", linewidth = 0.1
    ) +
    geom_polygon(
      data = contours[ClusterID == clusterID],
      mapping = aes(x = x, y = y, group = CellID),
      colour = specified_color, fill = specified_color, linewidth = 2
    ) +
    geom_polygon(
      data = contours[ClusterID == clusterID],
      mapping = aes(x = x, y = y, group = CellID),
      colour = specified_color, fill = specified_color, linewidth = 1.6
    ) +
    geom_polygon(
      data = contours[ClusterID == clusterID],
      mapping = aes(x = x, y = y, group = CellID),
      colour = specified_color, fill = specified_color, linewidth = 1.2
    ) +
    geom_polygon(
      data = contours[ClusterID == clusterID],
      mapping = aes(x = x, y = y, group = CellID),
      colour = specified_color, fill = specified_color, linewidth = 0.8
    ) +
    geom_polygon(
      data = contours[ClusterID == clusterID],
      mapping = aes(x = x, y = y, group = CellID),
      colour = specified_color, fill = specified_color, linewidth = 0.4
    )
}

plotimage1 <- plotTMEPhenotype(imageID = image_example, clusterID = 20)
plotimage2 <- plotTMEPhenotype(imageID = image_example, clusterID = 13)
plotimage3 <- plotTMEPhenotype(imageID = image_example, clusterID = 18)
plotimage4 <- plotTMEPhenotype(imageID = image_example, clusterID = 23)

combinedimage <- plotimage1 + plotimage2 + plotimage3 + plotimage4 + plot_layout(ncol = 4)

ggsave(
  filename = file.path(outdir, "singleTMEPhenotypeExamples.pdf"),
  plot = combinedimage,
  width = 14, height = 4, units = "in"
)

message("Wrote: ", file.path(outdir, "KCrossFunctionExamples.pdf"))
message("Wrote: ", file.path(outdir, "singleTMEPhenotypeExamples.pdf"))
