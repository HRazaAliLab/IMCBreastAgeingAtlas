library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID")

epithelial <- cellsclin[isEpithelial == TRUE]

epithelial[, isPrimaryPop := (isFOXA1Pos & isGATA3Pos & isERPos & isPRPos & isARPos & !isHER2Pos)]
epithelial[, isSecondaryPop := (!isFOXA1Pos & !isGATA3Pos & isERPos & !isPRPos & !isARPos & !isHER2Pos)]

aggregated_proportions <- epithelial[, .(
  TotalCells = .N,
  TotalPrimaryPop = sum(isPrimaryPop, na.rm = TRUE),
  TotalSecondaryPop = sum(isSecondaryPop, na.rm = TRUE)
), by = .(PatientID, Age)]

aggregated_proportions[, ProportionPrimary := TotalPrimaryPop / TotalCells]
aggregated_proportions[, ProportionSecondary := TotalSecondaryPop / TotalCells]
aggregated_proportions[is.na(ProportionPrimary), ProportionPrimary := 0]
aggregated_proportions[is.na(ProportionSecondary), ProportionSecondary := 0]

long_proportions <- melt(
  aggregated_proportions,
  id.vars = c("PatientID", "Age"),
  measure.vars = c("ProportionPrimary", "ProportionSecondary"),
  variable.name = "Population",
  value.name = "Proportion"
)

long_proportions[, Population := factor(
  Population,
  levels = c("ProportionPrimary", "ProportionSecondary"),
  labels = c("Primary (High Hormone)", "Secondary (ER+/FOXA1-/GATA3-)")
)]

primary_color <- "#4D4D4D"
secondary_color <- "#8B4513"
    
p <- ggplot(long_proportions, aes(x = Age, y = Proportion, color = Population)) +
  geom_point(alpha = 0.2, size = 0.1) +
  geom_smooth(method = "loess", se = FALSE) +
  scale_color_manual(values = c(
    "Primary (High Hormone)" = primary_color,
    "Secondary (ER+/FOXA1-/GATA3-)" = secondary_color
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", size = 18),
    axis.text.y = element_text(color = "black", size = 18),
    legend.position = "none",
    plot.margin = unit(c(10,10,10,10), "pt")
  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,0.2), breaks = c(0,0.2), labels = c("0", "0.1")) +
  scale_x_continuous(expand = c(0,0), limits = c(18, 70), breaks = c(20, 70))

ggsave(p, filename = file.path(outdir, "ageAssocPopulations.pdf"),
       units = "in", width = 2.5, height = 4)
