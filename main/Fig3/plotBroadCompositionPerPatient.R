library(here)
source(here("code", "header.R"))

cells    <- getCells()
clinical <- getClinical()
annotations <- getCellClusters()

outdir <- here("scratch/main/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells[, CellPhenotype_BroadCategory := NA_character_]
cells[isEpithelial == TRUE, CellPhenotype_BroadCategory := "Epithelial"]

for (i in unique(cells$ClusterID)) {
  if (i %in% annotations$ClusterID) {
    category <- annotations[ClusterID == i, Type]
    cells[isEpithelial == FALSE & ClusterID == i,
          CellPhenotype_BroadCategory := category]
  }
}

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

broadphenotypecolors <- getNormalBreastProjectColours()$BroadPhenotype

cell_counts <- cellsclin[, .N, by = .(PatientID, CellPhenotype_BroadCategory)]
cell_counts[, Total := sum(N), by = PatientID]
cell_counts[, Proportion := fifelse(CellPhenotype_BroadCategory == "Stromal", N / Total, 0), by = PatientID]
proportions <- cell_counts[, .(StromalProportion = max(Proportion)), by = PatientID]
ordered_patient_ids <- proportions[order(StromalProportion)]$PatientID
cell_counts[, PatientID := factor(PatientID, levels = ordered_patient_ids)]

cell_counts[, Total := sum(N), by = PatientID]
cell_counts[, Proportion := N / Total]

proportions_per_patient <- cell_counts[, .(PatientID, CellPhenotype_BroadCategory, Proportion)]
average_proportions <- proportions_per_patient[, .(AverageProportion = mean(Proportion)), by = CellPhenotype_BroadCategory]

average_proportions[CellPhenotype_BroadCategory %in% c("Epithelial", "Immune", "Stromal"),
                    CellPhenotype_BroadCategory := factor(CellPhenotype_BroadCategory,
                                                          levels = c("Epithelial", "Immune", "Stromal"))]
average_proportions <- average_proportions[order(-CellPhenotype_BroadCategory)]

average_proportions[, pct := AverageProportion * 100]
average_proportions[, floor_pct := floor(pct)]
average_proportions[, frac := pct - floor_pct]
to_allocate <- 100 - sum(average_proportions$floor_pct)
average_proportions[order(-frac)[1:to_allocate], floor_pct := floor_pct + 1]
average_proportions[, RoundedPct := floor_pct]

average_proportions[, CumulativeProportion := cumsum(AverageProportion)]
average_proportions[, LabelPosition := CumulativeProportion - (AverageProportion / 2)]

broadplot <- ggplot(cell_counts, aes(x = PatientID, y = N, fill = CellPhenotype_BroadCategory)) +
  geom_bar(stat = "identity", position = "fill", width = 1.5) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0, 1)) +
  scale_fill_manual(values = broadphenotypecolors,
                    labels = c("Epithelial" = "epithelial", "Immune" = "immune", "Stromal" = "stromal")) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(y = "cell prop. per patient") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(0, 1, 0, 0), "pt"),
    legend.box = "horizontal",
    legend.text = element_text(size = 16)
  ) +
  guides(fill = guide_legend(nrow = 1))

average_plot <- ggplot(average_proportions, aes(x = 1, y = AverageProportion, fill = CellPhenotype_BroadCategory)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format(), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = broadphenotypecolors,
                    labels = c("Epithelial" = "epithelial", "Immune" = "immune", "Stromal" = "stromal")) +
  labs(title = "mean %") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.3, vjust = -4, size = 15),
    plot.margin = unit(c(0, 10, 10, 0), "pt"),
    plot.title.position = "plot"
  ) +
  geom_text(aes(label = RoundedPct, y = LabelPosition), color = "black", size = 5.4)

main <- broadplot + average_plot + plot_layout(widths = c(8, 1))
ggsave(file.path(outdir, "broadCompositionPerPatient.pdf"),
       plot = main, width = 4.5, height = 5, units = "in")

print(average_proportions)

library(here)
source(here("code", "header.R"))

cells    <- getCells()
clinical <- getClinical()
annotations <- getCellClusters()

outdir <- here("scratch/main/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- Broad phenotype labels ----
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

broadphenotypecolors <- getNormalBreastProjectColours()$BroadPhenotype

# ---- Counts + ordering by stromal proportion ----
cell_counts <- cellsclin[, .N, by = .(PatientID, CellPhenotype_BroadCategory)]
cell_counts[, Total := sum(N), by = PatientID]
cell_counts[, Proportion := fifelse(CellPhenotype_BroadCategory == "Stromal", N / Total, 0), by = PatientID]
proportions <- cell_counts[, .(StromalProportion = max(Proportion)), by = PatientID]
ordered_patient_ids <- proportions[order(StromalProportion)]$PatientID
cell_counts[, PatientID := factor(PatientID, levels = ordered_patient_ids)]

# recompute true proportions
cell_counts[, Total := sum(N), by = PatientID]
cell_counts[, Proportion := N / Total]

# ---- Mean composition (rounded to integer sum 100) ----
proportions_per_patient <- cell_counts[, .(PatientID, CellPhenotype_BroadCategory, Proportion)]
average_proportions <- proportions_per_patient[, .(AverageProportion = mean(Proportion)), by = CellPhenotype_BroadCategory]

average_proportions[CellPhenotype_BroadCategory %in% c("Epithelial", "Immune", "Stromal"),
                    CellPhenotype_BroadCategory := factor(CellPhenotype_BroadCategory,
                                                          levels = c("Epithelial", "Immune", "Stromal"))]
average_proportions <- average_proportions[order(-CellPhenotype_BroadCategory)]

average_proportions[, pct := AverageProportion * 100]
average_proportions[, floor_pct := floor(pct)]
average_proportions[, frac := pct - floor_pct]
to_allocate <- 100 - sum(average_proportions$floor_pct)
average_proportions[order(-frac)[1:to_allocate], floor_pct := floor_pct + 1]
average_proportions[, RoundedPct := floor_pct]

average_proportions[, CumulativeProportion := cumsum(AverageProportion)]
average_proportions[, LabelPosition := CumulativeProportion - (AverageProportion / 2)]

broadplot <- ggplot(cell_counts, aes(x = PatientID, y = N, fill = CellPhenotype_BroadCategory)) +
  geom_bar(stat = "identity", position = "fill", width = 1.5) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)), breaks = c(0, 1)) +
  scale_fill_manual(values = broadphenotypecolors,
                    labels = c("Epithelial" = "epithelial", "Immune" = "immune", "Stromal" = "stromal")) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(y = "cell prop. per patient") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 17),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = unit(c(0, 1, 0, 0), "pt"),
    legend.box = "horizontal",
    legend.text = element_text(size = 16)
  ) +
  guides(fill = guide_legend(nrow = 1))

average_plot <- ggplot(average_proportions, aes(x = 1, y = AverageProportion, fill = CellPhenotype_BroadCategory)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = broadphenotypecolors,
                    labels = c("Epithelial" = "epithelial", "Immune" = "immune", "Stromal" = "stromal")) +
  labs(title = "mean %") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.3, vjust = -4, size = 15),
    plot.margin = unit(c(0, 10, 10, 0), "pt"),
    plot.title.position = "plot"
  ) +
  geom_text(aes(label = RoundedPct, y = LabelPosition), color = "black", size = 5.4)

main <- broadplot + average_plot + plot_layout(widths = c(8, 1))
ggsave(file.path(outdir, "BroadCompositionPerPatient.pdf"),
       plot = main, width = 4.5, height = 5, units = "in")

print(average_proportions)
