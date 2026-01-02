library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cells[, AxisRatio := MajorAxisLength / MinorAxisLength]

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellsclin[, AgeGroup := fifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
cellsclin[, ClusterID := as.integer(ClusterID)]

cellsclin[, CellPhenotype_BroadCategory := NA_character_]
cellsclin[isEpithelial == TRUE, CellPhenotype_BroadCategory := "Epithelial"]
for (cid in unique(cellsclin$ClusterID)) {
  if (cid %in% annotations$ClusterID) {
    type <- annotations[ClusterID == cid, Type][1]
    cellsclin[isEpithelial == FALSE & ClusterID == cid, CellPhenotype_BroadCategory := type]
  }
}

avg_axis_broad <- cellsclin[, .(AvgAxisRatio = mean(AxisRatio, na.rm = TRUE)),
                            by = .(PatientID, CellPhenotype_BroadCategory, AgeGroup)]

avg_axis_broad[, CellPhenotype_BroadCategory :=
                 factor(CellPhenotype_BroadCategory, levels = c("Epithelial", "Immune", "Stromal"))
]

patient_counts_broad <- unique(cellsclin[, .(CellPhenotype_BroadCategory, PatientID)])[,
                                                                                       .N, by = CellPhenotype_BroadCategory
]
avg_axis_broad <- merge(avg_axis_broad, patient_counts_broad, by = "CellPhenotype_BroadCategory", all.x = TRUE)
avg_axis_broad[, annotation_position := (1.6 - 1.2) * 0.015 + 1.2]

broad_cols <- getNormalBreastProjectColours()$BroadPhenotype

p7a <- ggplot(avg_axis_broad, aes(x = AvgAxisRatio, y = factor(CellPhenotype_BroadCategory, levels = c("Stromal", "Immune", "Epithelial")),
                                  fill = CellPhenotype_BroadCategory)) +
  geom_boxplot(outlier.size = 0.1) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  scale_fill_manual(values = broad_cols) +
  scale_x_continuous(breaks = c(1.2, 1.4, 1.6), limits = c(1.2, 1.6),
                     expand = expand_scale(mult = c(0.1, 0)), position = "top") +
  theme(
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    axis.title.y = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin = unit(c(40, 20, 0, 5), "pt")
  ) +
  scale_y_discrete(labels = c("stromal", "immune", "epithelial")) +
  geom_text(aes(label = N, x = annotation_position),
            check_overlap = TRUE, size = 4, vjust = 0.5, hjust = 0.9)

ggsave(
  p7a,
  filename = file.path(outdir, "axisRatioBroadPhenotypeComparison.pdf"),
  units = "in", width = 4.5, height = 1.75
)

avg_axis <- cellsclin[, .(AvgAxisRatio = mean(AxisRatio, na.rm = TRUE)),
                      by = .(PatientID, ClusterID, AgeGroup)]
avg_axis_typed <- merge(avg_axis, annotations, by = "ClusterID", all.x = TRUE)
avg_axis_typed[, Type := factor(Type, levels = c("Epithelial", "Immune", "Stromal"))]
avg_axis_typed[, ClusterID_f := factor(ClusterID, levels = unique(ClusterID[order(Type, PrintOrder)]))]

labels <- setNames(annotations$BackupLabel, annotations$ClusterID)
colors <- setNames(annotations$Colour, annotations$ClusterID)

type_labels <- c("epithelial", "immune", "stromal")
label_mapping <- setNames(type_labels, c("Epithelial", "Immune", "Stromal"))

patient_counts <- unique(cellsclin[, .(ClusterID, PatientID)])[, .N, by = ClusterID]
avg_axis_typed <- merge(avg_axis_typed, patient_counts, by = "ClusterID", all.x = TRUE)
avg_axis_typed[, annotation_position := (2 - 1) * 0.015 + 1]

p7b <- ggplot(avg_axis_typed, aes(x = AvgAxisRatio, y = factor(ClusterID_f), fill = factor(ClusterID))) +
  geom_boxplot(outlier.size = 0.1) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  scale_x_continuous(breaks = c(1, 1.5, 2), limits = c(1, 2),
                     expand = expand_scale(mult = c(0.1, 0)), position = "top") +
  theme(
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 14),
    legend.position = "none",
    plot.margin = unit(c(40, 20, 0, 20), "pt"),
    strip.text.y = element_blank(),
    strip.background = element_blank()
  ) +
  facet_grid(rows = vars(Type), scales = "free_y", space = "free",
             labeller = labeller(Type = label_mapping)) +
  scale_y_discrete(labels = sapply(labels[levels(avg_axis_typed$ClusterID_f)], function(x) TeX(mkTeX(x, ignoreChar="/")))) +
  geom_text(aes(label = N, x = annotation_position),
            size = 4, vjust = 0.5, hjust = 1.05, check_overlap = TRUE)

ggsave(
  p7b,
  filename = file.path(outdir, "axisRatioPhenotypeComparison.pdf"),
  units = "in", width = 5.2, height = 6
)
