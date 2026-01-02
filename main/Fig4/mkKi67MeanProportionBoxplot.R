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

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

mean_proportions_by_patient <- cellsclin[
  , .(Proportion = mean(isKi67Pos, na.rm = TRUE)),
  by = .(PatientID, CellPhenotype_BroadCategory)
]

broadphenotypecolors <- getNormalBreastProjectColours()$BroadPhenotype

box_plot <- ggplot(
  mean_proportions_by_patient,
  aes(x = CellPhenotype_BroadCategory, y = Proportion, colour = CellPhenotype_BroadCategory)
) +
  geom_boxplot(outlier.colour = NULL, outlier.size = 0.1) +
  scale_colour_manual(values = broadphenotypecolors) +
  theme_classic() +
  labs(x = NULL, y = NULL, fill = "Cell Phenotype") +
  theme(
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title = element_text(size = 14),
    legend.position = "none",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.23, 0.02)),
    breaks = c(0, 0.2),
    labels = function(x) ifelse(x == 0, "0", format(x))
  )

ggsave(
  plot = box_plot,
  filename = file.path(outdir, "Ki67MeanProportionBoxplot.pdf"),
  width = 2, height = 1
)