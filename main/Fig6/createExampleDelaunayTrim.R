library(here)
source(here("code", "header.R"))

#additional dependency
library(deldir)

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
clinical <- getClinical()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)

# ============================================================
# RGB + contours + trimmed overlay
# ============================================================
imgID <- "NormalBreastImgK_R3C10"
exampleimage <- cellsclin[ImageID == imgID]

dt <- deldir(exampleimage$CenterX, exampleimage$CenterY)
edges <- data.frame(
  x1 = dt$delsgs$x1, y1 = dt$delsgs$y1,
  x2 = dt$delsgs$x2, y2 = dt$delsgs$y2
)
edges$distance <- sqrt((edges$x1 - edges$x2)^2 + (edges$y1 - edges$y2)^2)
cutoff <- quantile(edges$distance, 0.99, na.rm = TRUE)
trimmed_edges <- as.data.table(edges[edges$distance < cutoff, ])
edge_list <- trimmed_edges[, .(x1, y1, x2, y2)]

cellSub <- exampleimage[, .(CellID, SpatialClusterID, isEpithelial)]
imgPaths <- getImagePaths(imgID)
contours <- getContour(imgPaths$Cell)
contours <- merge(contours, cellSub, by = "CellID", all.x = TRUE)

p <- mkRGBplot(imgID, RGB = c("pH2AX", "FOXP3", "IDO"), noRed = TRUE)$plot +
  geom_polygon(data = contours,
               aes(x = x, y = y, group = CellID),
               colour = "darkgray", fill = NA, linewidth = 0.2) +
  geom_polygon(data = contours[isEpithelial == TRUE],
               aes(x = x, y = y, group = CellID),
               colour = "#6C7ACE", fill = NA, linewidth = 0.2) +
  geom_segment(data = edge_list,
               aes(x = x1, y = y1, xend = x2, yend = y2),
               colour = "white", linewidth = 0.05) +
  labs(title = "graph edges, trimmed at 99th percentile (pixel distance)")

ggsave(file.path(outdir, "exampleDelaunayTrim99.pdf"), p, width = 11, height = 7)

# ============================================================
# Delaunay graph edges only within same neighbourhood label
# ============================================================
plot_delaunay_graph <- function(image_id, cells_dt, comm_ann) {
  image_data <- cells_dt[ImageID == image_id]
  if (nrow(image_data) == 0) stop("No cells for ImageID: ", image_id)
  
  image_data[, CellID := as.character(CellID)]
  image_data[, SpatialClusterID := as.character(SpatialClusterID)]
  comm_ann[, SpatialClusterID := as.character(SpatialClusterID)]
  
  dt <- deldir(image_data$CenterX, image_data$CenterY)
  
  edges <- data.table(
    from = image_data$CellID[dt$delsgs$ind1],
    to   = image_data$CellID[dt$delsgs$ind2]
  )
  
  x_from <- image_data$CenterX[match(edges$from, image_data$CellID)]
  y_from <- image_data$CenterY[match(edges$from, image_data$CellID)]
  x_to   <- image_data$CenterX[match(edges$to,   image_data$CellID)]
  y_to   <- image_data$CenterY[match(edges$to,   image_data$CellID)]
  
  edges[, length := sqrt((x_from - x_to)^2 + (y_from - y_to)^2)]
  edges <- edges[!is.na(length)]
  
  threshold <- quantile(edges$length, 0.99, na.rm = TRUE)
  edges <- edges[length <= threshold]
  
  from_sp <- image_data$SpatialClusterID[match(edges$from, image_data$CellID)]
  to_sp   <- image_data$SpatialClusterID[match(edges$to,   image_data$CellID)]
  
  edges[, from_label := comm_ann$Label[match(from_sp, comm_ann$SpatialClusterID)]]
  edges[, to_label   := comm_ann$Label[match(to_sp,   comm_ann$SpatialClusterID)]]
  
  edges <- edges[!is.na(from_label) & !is.na(to_label)]
  edges <- edges[from_label == to_label]
  
  ggplot() +
    theme_void() +
    theme(panel.background = element_rect(fill = "black", color = NA),
          plot.background  = element_rect(fill = "black", color = NA)) +
    labs(title = paste("Neighbourhood-only Delaunay:", image_id)) +
    geom_segment(
      data = edges,
      aes(
        x    = image_data$CenterX[match(from, image_data$CellID)],
        y    = image_data$CenterY[match(from, image_data$CellID)],
        xend = image_data$CenterX[match(to,   image_data$CellID)],
        yend = image_data$CenterY[match(to,   image_data$CellID)],
        color = from_label
      ),
      linewidth = 0.2
    ) +
    scale_color_manual(values = setNames(comm_ann$Colour, comm_ann$Label)) +
    theme(legend.position = "none") +
    coord_fixed(ratio = 1)
}


communityexample <- plot_delaunay_graph(imgID, cellsclin, communityannotations)
ggsave(filename = file.path(outdir, "communityIdentificationExample.pdf"),
       plot = communityexample, width = 7, height = 7)
