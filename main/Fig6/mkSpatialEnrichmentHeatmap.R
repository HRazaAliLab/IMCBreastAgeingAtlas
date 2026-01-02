library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()

ann0  <- getCellClusters()
comm0 <- fread(here("data", "derived", "spatialClusterAnnotation.csv"))

cells[, ClusterID := as.character(ClusterID)]
cells[, SpatialClusterID := as.character(SpatialClusterID)]
ann0[,  ClusterID := as.character(ClusterID)]
comm0[, SpatialClusterID := as.character(SpatialClusterID)]
comm0[, ClusterID := SpatialClusterID]

pgclust_totals <- cells[, .(pgclust_count = .N), by = ClusterID]
spatial_totals <- cells[, .(spatialclust_count = .N), by = SpatialClusterID]

cells2 <- merge(cells, pgclust_totals, by = "ClusterID")
cells2 <- merge(cells2, spatial_totals, by = "SpatialClusterID")

total_cells <- nrow(cells2)
cells2[, expected_count := (pgclust_count / total_cells) * spatialclust_count, by = SpatialClusterID]

observed_counts <- cells2[, .(observed_count = .N), by = .(ClusterID, SpatialClusterID)]
enrichment_dt <- merge(
  observed_counts,
  unique(cells2[, .(ClusterID, SpatialClusterID, expected_count)]),
  by = c("ClusterID", "SpatialClusterID")
)
enrichment_dt[, enrichment := observed_count / expected_count]
enrichment_dt <- unique(enrichment_dt[, .(ClusterID, SpatialClusterID, enrichment)])

ann_ord <- copy(ann0)[order(Type, PrintOrder)]
pheno_levels <- ann_ord$ClusterID

comm_ord <- copy(comm0)[rev(order(PrintOrder))]
spatial_levels <- comm_ord$ClusterID

enrichment_dt <- merge(enrichment_dt, ann0[, .(ClusterID, Type)], by = "ClusterID", all.x = TRUE)

enrichment_dt[, ClusterID := factor(ClusterID, levels = pheno_levels)]
enrichment_dt[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]

labels <- setNames(ann0$TeXClusterLabel, ann0$ClusterID)
communitylabels <- setNames(comm_ord$Label, comm_ord$ClusterID)

# enrichment heatmap
enrichmentplot <- ggplot(enrichment_dt, aes(x = ClusterID, y = SpatialClusterID, fill = enrichment)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white", high = "black",
    na.value = "white",
    limits = c(1, max(enrichment_dt$enrichment, na.rm = TRUE)),
    breaks = c(1, 8)
  ) +
  theme_classic() +
  scale_x_discrete(labels = sapply(labels, TeX)) +
  scale_y_discrete(labels = sapply(communitylabels, TeX)) +
  labs(x = NULL, y = NULL, fill = "enrichment") +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0,0,0,0), "pt"),
    legend.position = "none",
    panel.border = element_blank(),
    strip.text.x = element_blank()
  ) +
  facet_grid(cols = vars(Type), scales = "free_x", space = "free")

#spatial cluster colorbar (y)
colors_sp <- setNames(comm_ord$Colour, comm_ord$ClusterID)
comm_bar <- copy(comm_ord)
comm_bar[, ClusterID := factor(ClusterID, levels = spatial_levels)]

spatialclust_colorbar_plot <- ggplot(comm_bar, aes(x = factor(1), y = ClusterID, fill = ClusterID)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_y_discrete(limits = as.factor(spatial_levels), labels = sapply(communitylabels, TeX)) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = colors_sp) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 14),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  )

# Phenotype cluster colorbar (x)
ann_bar <- copy(ann0)
ann_bar[, BackupTeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
phenotypelabels <- setNames(ann_bar$BackupTeXClusterLabel, ann_bar$ClusterID)
phenotypecolors <- setNames(ann_bar$Colour, ann_bar$ClusterID)

ann_bar[, ClusterID := factor(ClusterID, levels = pheno_levels)]

phenotypeclust_colorbar_plot <- ggplot(ann_bar, aes(y = factor(1), x = ClusterID, fill = ClusterID)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_x_discrete(labels = sapply(phenotypelabels, TeX)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = phenotypecolors) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 14),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  facet_grid(~ Type, scales = "free_x", space = "free")

# Marker proportion tiles by neighbourhood
spatial_totals2 <- cells2[, .(total_count = .N), by = SpatialClusterID]

calculateProportions <- function(conditionCol) {
  positiveCounts <- cells2[get(conditionCol) == TRUE, .N, by = SpatialClusterID]
  proportions <- merge(positiveCounts, spatial_totals2, by = "SpatialClusterID", all.y = TRUE)
  proportions[is.na(N), N := 0]
  proportions[, Percent := (N / total_count) * 100]
  proportions[, .(SpatialClusterID, Percent)]
}

ki67Proportions <- calculateProportions("isKi67Pos")
ERProportions   <- calculateProportions("isERPos")
PRProportions   <- calculateProportions("isPRPos")
ARProportions   <- calculateProportions("isARPos")

ki67Proportions[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
ERProportions[,   SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
PRProportions[,   SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
ARProportions[,   SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]

tile_plot <- function(dt, subtitle, lim_hi, fmt = "%.0f", high_col = "darkblue") {
  ggplot(dt, aes(x = factor(1), y = SpatialClusterID, fill = Percent)) +
    geom_tile(color = "white", linewidth = 0.1) +
    geom_text(aes(label = sprintf(fmt, Percent)), vjust = 0.5, color = "white", size = 4.4) +
    scale_fill_gradient(low = "darkgrey", high = high_col, limits = c(0, lim_hi)) +
    labs(subtitle = subtitle) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.line = element_blank(),
      legend.position = "none",
      plot.subtitle = element_text(hjust = 0.5, vjust = -1.5, size = 12),
      plot.margin = unit(c(0,0,0,0), "pt")
    )
}

ki67_plot <- tile_plot(ki67Proportions, "Ki67", lim_hi = 3.8, fmt = "%.1f", high_col = "darkblue")
ER_plot   <- tile_plot(ERProportions,   "ER",   lim_hi = 40,  fmt = "%.0f", high_col = "darkred")
PR_plot   <- tile_plot(PRProportions,   "PR",   lim_hi = 40,  fmt = "%.0f", high_col = "darkred")
AR_plot   <- tile_plot(ARProportions,   "AR",   lim_hi = 40,  fmt = "%.0f", high_col = "darkred")

combinedplot <- spatialclust_colorbar_plot + enrichmentplot + ki67_plot + ER_plot + PR_plot + AR_plot +
  plot_spacer() + phenotypeclust_colorbar_plot +
  plot_spacer() + plot_spacer() + plot_spacer() + plot_spacer() +
  plot_layout(widths = c(0.1, 7, 0.52, 0.45, 0.45, 0.45),
              heights = c(5, 0.2), nrow = 2)

pdf(file.path(outdir, "spatialEnrichmentHeatmap.pdf"), width = 9.2, height = 4)
print(combinedplot)
dev.off()
