library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellsclin[, ClusterID := as.integer(ClusterID)]

cellsclin[, CellPhenotype_BroadCategory := NA_character_]
cellsclin[isEpithelial == TRUE, CellPhenotype_BroadCategory := "Epithelial"]
for (cid in unique(cellsclin$ClusterID)) {
  if (cid %in% annotations$ClusterID) {
    type <- annotations[ClusterID == cid, Type][1]
    cellsclin[isEpithelial == FALSE & ClusterID == cid, CellPhenotype_BroadCategory := type]
  }
}

categories <- c("Epithelial", "Immune", "Stromal")
categorynames <- c("epithelial", "immune", "stromal")

ymin_breaks <- c(40, 30, 30)
ymax_breaks <- c(80, 70, 60)

broad_cols <- getNormalBreastProjectColours()$BroadPhenotype

plots_list <- list()

for (i in seq_along(categories)) {
  cat_i <- categories[i]
  point_color <- broad_cols[cat_i]
  ylo <- ymin_breaks[i]
  yhi <- ymax_breaks[i]
  
  subset_dt <- cellsclin[CellPhenotype_BroadCategory == cat_i]
  # mean area per patient
  pts <- subset_dt[, .(MeanArea = mean(CellArea, na.rm = TRUE)), by = .(PatientID, Age)]
  
  ct <- cor.test(pts$Age, pts$MeanArea, method = "spearman")
  
  p <- ggplot(pts, aes(x = Age, y = MeanArea)) +
    geom_point(size = 0.1, alpha = 0.3, color = point_color) +
    geom_smooth(method = "loess", color = point_color, se = FALSE) +
    labs(subtitle = categorynames[i]) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(limits = c(ylo, yhi), expand = c(0,0), breaks = c(ylo, yhi)) +
    scale_x_continuous(breaks = c(20, 70), limits = c(16, 70), expand = c(0,0)) +
    annotate("text", x = 68, y = ((yhi - ylo) * 0.90) + ylo,
             label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(ct$p.value)))),
             size = 4.5, hjust = 1, color = "black") +
    annotate("text", x = 68, y = ((yhi - ylo) * 0.75) + ylo,
             label = TeX(paste0("$\\rho$ = ", format(ct$estimate, digits = 2))),
             size = 4.5, hjust = 1, color = "black") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.subtitle = element_text(color = "black", size = 15, hjust = 0.5),
      axis.text.x = element_text(color = "black", size = 15),
      axis.text.y = element_text(color = "black", size = 15)
    )
  
  if (i == 3) {
    p <- p + labs(subtitle = "stromal", x = "age") +
      theme(axis.title.x = element_text(size = 15))
  }
  plots_list[[cat_i]] <- p
}

plot_grid <- wrap_plots(plots_list, ncol = 1)

blanklabelploty <- ggplot() +
  labs(y = "average cell area (µm²)") +
  theme_classic() +
  theme(plot.margin = unit(c(0,5,0,0), "pt"),
        axis.title.y = element_text(color = "black", size = 16)) +
  guides(x = "none", y = "none")

plot_grid <- blanklabelploty + plot_grid + plot_layout(widths = c(1, 1000))

ggsave(file.path(outdir, "cellAreaBroadCategoryScatterplots.pdf"),
       plot_grid, width = 2.5, height = 7, units = "in")
