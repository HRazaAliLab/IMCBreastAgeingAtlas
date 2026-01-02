library(here)
source(here("code", "header.R"))

cells        <- getCells()
clinical     <- getClinical()
annotations  <- getCellClusters()

outdir <- here("scratch/main/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

tissueAreasPatients <- clinical[, .(
  totalTissueArea = sum(TissueArea, na.rm = TRUE)
), by = PatientID]

cells[, CellPhenotype_BroadCategory := NA_character_]
cells[isEpithelial == TRUE, CellPhenotype_BroadCategory := "Epithelial"]

for (i in unique(cells$ClusterID)) {
  if (i %in% annotations$ClusterID) {
    category <- annotations[ClusterID == i, Type]
    cells[isEpithelial == FALSE & ClusterID == i,
          CellPhenotype_BroadCategory := category]
  }
}

cells[CellPhenotype_BroadCategory == "Stromal", CellPhenotype_BroadCategory := "Stromal"]
cells[CellPhenotype_BroadCategory == "Immune",  CellPhenotype_BroadCategory := "Immune"]

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

categories  <- c("Epithelial", "Immune", "Stromal")
total_cells <- cellsclin[, .(TotalCells = .N), by = .(PatientID, Age)]
all_patients <- unique(cellsclin[, .(PatientID, Age)])

combined_proportions <- rbindlist(lapply(categories, function(category) {
  dt <- cellsclin[CellPhenotype_BroadCategory == category,
                  .(CellCount = .N), by = .(PatientID, Age)]
  dt <- merge(all_patients, dt, by = c("PatientID", "Age"), all.x = TRUE)
  dt[is.na(CellCount), CellCount := 0]
  dt[, Category := category]
  dt
}), fill = TRUE)

combined_proportions <- merge(combined_proportions, total_cells, by = c("PatientID", "Age"))
combined_proportions <- merge(combined_proportions, tissueAreasPatients, by = "PatientID", all.x = TRUE)

combined_proportions[, Proportion := CellCount / TotalCells]
combined_proportions[, Density := 1e6 * CellCount / totalTissueArea]
combined_proportions <- na.omit(combined_proportions)

broadphenotypecolors <- getNormalBreastProjectColours()$BroadPhenotype

ymaxes <- c(1300, 300, 1200)
ymins  <- c(0, 0, 0)

plots_list <- list()
i <- 1

for (category in categories) {
  dt <- combined_proportions[Category == category]
  correlation_test <- cor.test(dt$Age, dt$Density, method = "spearman")
  point_color <- broadphenotypecolors[i]
  
  p <- ggplot(dt, aes(x = Age, y = Density)) +
    geom_point(size = 0.1, alpha = 0.3, color = point_color) +
    geom_smooth(method = "lm", color = point_color, se = FALSE) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(ymins[i], ymaxes[i]), expand = c(0, 0),
      breaks = c(ymins[i], ymaxes[i]),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(expand = c(0, 0), limits = c(15, 70), breaks = c(20, 70)) +
    annotate(
      "text", x = 16, y = ((ymaxes[i] - ymins[i]) * 0.07) + ymins[i],
      label = TeX(paste0("$\\rho$ = ", format(correlation_test$estimate, digits = 2))),
      size = 6.5, hjust = 0, color = "black"
    ) +
    annotate(
      "text", x = 16, y = ymaxes[i] - ((ymaxes[i] - ymins[i]) * 0.07),
      label = TeX(paste0(
        "pval = ",
        mkEnumPower(ifelse(
          correlation_test$p.value < 0.001,
          sprintf("%.0e", signif(correlation_test$p.value, 1)),
          signif(correlation_test$p.value, 1)
        ))
      )),
      parse = TRUE, size = 6.5, hjust = 0, color = "black"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(color = "black", size = 18),
      axis.text.y = element_text(color = "black", size = 18),
      plot.margin = unit(c(10, 10, 0, 0), "pt")
    )
  
  if (category == "Stromal") {
    p <- p + labs(x = "Age") + theme(axis.title.x = element_text(size = 19))
  }
  
  plots_list[[category]] <- p
  i <- i + 1
}

plot_grid <- wrap_plots(plots_list, nrow = 3)

blanklabelploty <- ggplot() +
  labs(y = TeX("density (cells/mm$^2$)")) +
  theme_classic() +
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "pt"),
    axis.title.y = element_text(size = 19)
  ) +
  guides(x = "none", y = "none")

plot_grid <- blanklabelploty + plot_grid + plot_layout(widths = c(1, 1000))

pdf(file.path(outdir, "ageScatterplotsCompartments.pdf"), width = 3, height = 6)
print(plot_grid)
dev.off()
