library(here)
source(here("code", "header.R"))

#Required dependency below
library(deldir)

outdir <- here("scratch", "main", "Fig7")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
ctx   <- getCellContext()

cells <- merge(cells, ctx[, .(ImageID, CellID, TissueStructure)], by=c("ImageID","CellID"), all.x=TRUE)

plot_specific_image <- function(image_id, dt_cells) {
  d <- dt_cells[ImageID == image_id]
  d <- d[!is.na(TissueStructure)]
  if (nrow(d) < 3) stop("Not enough cells w/ TissueStructure in this image.")
  
  dt <- deldir(d$CenterX, d$CenterY)
  edges <- with(dt, data.frame(
    from = d$CellID[delsgs$ind1],
    to   = d$CellID[delsgs$ind2]
  ))
  
  edges$length <- sqrt(
    (d$CenterX[match(edges$from, d$CellID)] - d$CenterX[match(edges$to, d$CellID)])^2 +
      (d$CenterY[match(edges$from, d$CellID)] - d$CenterY[match(edges$to, d$CellID)])^2
  )
  thr <- quantile(edges$length, 0.99, na.rm=TRUE)
  edges <- edges[edges$length <= thr, , drop=FALSE]
  
  edges$from_label <- d$TissueStructure[match(edges$from, d$CellID)]
  edges$to_label   <- d$TissueStructure[match(edges$to,   d$CellID)]
  edges <- edges[edges$from_label == edges$to_label, , drop=FALSE]
  
  color_mapping <- c("duct"="#4169E1", "lobule"="#DC143C")
  
  ggplot() +
    theme_void() +
    theme(panel.background = element_rect(fill="black", color=NA),
          plot.background  = element_rect(fill="black", color=NA)) +
    geom_segment(
      data = edges,
      aes(
        x = d$CenterX[match(from, d$CellID)],
        y = d$CenterY[match(from, d$CellID)],
        xend = d$CenterX[match(to, d$CellID)],
        yend = d$CenterY[match(to, d$CellID)],
        color = from_label
      ),
      linewidth = 0.2
    ) +
    scale_color_manual(values = color_mapping) +
    theme(legend.position="none") +
    coord_fixed(ratio=1)
}

p <- plot_specific_image("NormalBreastImgA_R2C10", cells)
ggsave(filename = file.path(outdir, "ductLobuleStructureExample.pdf"),
       plot = p, width = 5, height = 5, units = "in")
