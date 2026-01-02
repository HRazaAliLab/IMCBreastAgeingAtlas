library(here)
source(here("code", "header.R"))

outdir <- here("scratch/main/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cells[, ClusterID := as.integer(as.character(ClusterID))]
annotations[, ClusterID := as.integer(as.character(ClusterID))]

cells[, CellPhenotype_BroadCategory := NA_character_]
cells[isEpithelial == TRUE, CellPhenotype_BroadCategory := "Epithelial"]

annot_map <- unique(annotations[, .(ClusterID, Type)])
cells[isEpithelial == FALSE,
      CellPhenotype_BroadCategory := annot_map[.SD, on = "ClusterID", x.Type]]

cells[!CellPhenotype_BroadCategory %in% c("Epithelial", "Immune", "Stromal"),
      CellPhenotype_BroadCategory := NA_character_]

cells[, CellPhenotype_BroadCategory := factor(
  CellPhenotype_BroadCategory,
  levels = c("Epithelial", "Immune", "Stromal")
)]

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

categories    <- c("Epithelial", "Immune", "Stromal")
categorynames <- c("epithelial", "immune", "stromal")

plots_list <- list()
specific_y_breaks <- c(0.3, 0.15, 0.05)

broadphenotypecolors <- getNormalBreastProjectColours()$BroadPhenotype

for (i in seq_along(categories)) {
  category <- categories[i]
  y_breaks <- specific_y_breaks[i]
  point_color <- broadphenotypecolors[i]
  
  subset_data <- cellsclin[CellPhenotype_BroadCategory == category]
  proportions <- subset_data[, .(Proportion = mean(isKi67Pos, na.rm = TRUE)), by = .(PatientID, Age)]
  proportions <- na.omit(proportions)
  
  correlation_test <- cor.test(proportions$Age, proportions$Proportion, method = "spearman")
  
  p <- ggplot(proportions, aes(x = Age, y = sqrt(Proportion))) +
    geom_point(size = 0.1, alpha = 0.3, color = point_color) +
    geom_smooth(method = "loess", color = point_color, se = FALSE) +
    labs(subtitle = categorynames[i]) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, y_breaks),
      expand = c(0,0),
      breaks = c(0, y_breaks),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 70), limits = c(15, 70), expand = c(0,0)) +
    annotate("text", x = 68, y = y_breaks * 0.90,
             label = TeX(paste0("pval = ", mkEnumPower(format_custom_pval(correlation_test$p.value)))),
             size = 4.5, hjust = 1, color = "black") +
    annotate("text", x = 68, y = y_breaks * 0.65,
             label = TeX(paste0("$\\rho$ = ", format(correlation_test$estimate, digits = 2))),
             size = 4.5, hjust = 1, color = "black") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.subtitle = element_text(color = "black", size = 14, hjust = 0.5),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14)
    )
  
  if (i == 3) {
    p <- p + labs(subtitle = "stromal", x = "age") +
      theme(axis.title.x = element_text(size = 12))
  }
  
  plots_list[[category]] <- p
}

plot_grid <- wrap_plots(plots_list, ncol = 1)

blanklabelploty <- ggplot() +
  labs(y = "proliferative proportion of compartment (sqrt)") +
  theme_classic() +
  guides(x = "none", y = "none")

plot_grid <- blanklabelploty + plot_grid + plot_layout(widths = c(1, 1000))

pdf(file.path(outdir, "Ki67CompartmentScatterplots.pdf"), width = 2.5, height = 5)
print(plot_grid)
dev.off()
