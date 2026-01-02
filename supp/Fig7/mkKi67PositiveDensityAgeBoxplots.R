##########################
# Ki67+ cell density across age, per compartment (boxplots)
##########################

library(here)
library(data.table)
library(ggplot2)
library(patchwork)

source(here("code", "header.R"))

outdir <- here("scratch/main/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID")

# ---- total tissue area per patient from clinical (NO tissueAreas file) ----
tissueAreasPatients <- clinical[, .(
  totalTissueArea = sum(TissueArea, na.rm = TRUE)
), by = PatientID]

# ---- patient grid ----
patients <- unique(cellsclin[!is.na(Age), .(PatientID, Age)])
comps <- unique(cellsclin$CellPhenotype_BroadCategory)

grid <- CJ(
  PatientID = patients$PatientID,
  CellPhenotype_BroadCategory = comps,
  unique = TRUE
)[, Age := patients[.SD, on = "PatientID", Age]]

# ---- Ki67+ counts per patient x compartment ----
ki67_counts <- cellsclin[
  !is.na(Age) & isKi67Pos == TRUE,
  .(Ki67_Count = .N),
  by = .(PatientID, CellPhenotype_BroadCategory)
]

ki67_density <- merge(
  grid, ki67_counts,
  by = c("PatientID", "CellPhenotype_BroadCategory"),
  all.x = TRUE
)
ki67_density[is.na(Ki67_Count), Ki67_Count := 0]

ki67_density <- merge(
  ki67_density, tissueAreasPatients,
  by = "PatientID", all.x = TRUE
)

ki67_density[, Density := 1e6 * Ki67_Count / totalTissueArea]
ki67_density[, AgeGroup := ifelse(Age < 50, "Below 50", "50 and above")]
ki67_density[, AgeGroup := factor(AgeGroup, levels = c("Below 50", "50 and above"))]

# ---- Wilcoxon per compartment (density) ----
test_results <- rbindlist(lapply(unique(ki67_density$CellPhenotype_BroadCategory), function(cat) {
  sub <- ki67_density[CellPhenotype_BroadCategory == cat]
  wt <- wilcox.test(Density ~ AgeGroup, data = sub, exact = FALSE, correct = TRUE)
  data.table(Compartment = cat, p_value = wt$p.value)
}))
test_results[, p_adj := p.adjust(p_value, method = "BH")]
print(test_results)

# ---- plotting ----
ymaxs <- c(Epithelial = 50, Immune = 5, Stromal = 5)

make_panel <- function(cat) {
  dat <- ki67_density[CellPhenotype_BroadCategory == cat]
  ymax <- ymaxs[[cat]]
  
  ggplot(dat, aes(AgeGroup, Density, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.1, width = 0.6) +
    scale_fill_manual(values = c("Below 50" = "gray80", "50 and above" = "gray40")) +
    scale_y_continuous(limits = c(0, ymax), breaks = c(0, ymax), expand = c(0, 0)) +
    theme_classic() +
    labs(
      x = NULL,
      y = if (cat == "Epithelial") expression("Ki67+ density (cells/"*mm^2*")") else NULL
    ) +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = 12),
      legend.position = "none"
    )
}

panels <- lapply(names(ymaxs), make_panel)
combined_plot <- panels[[1]] + panels[[2]] + panels[[3]] + plot_layout(ncol = 3, guides = "collect")

ggsave(
  filename = file.path(outdir, "Ki67_Density_By_AgeBoxplots.pdf"),
  plot = combined_plot,
  width = 4.5,
  height = 2
)
