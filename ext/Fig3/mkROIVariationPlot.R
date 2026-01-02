# ---------------------- ROI VARIATION (VARIANCE OF PROPORTION) ----------------------
library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

annotations <- getCellClusters()
epiannotations <- annotations[isEpithelial == TRUE]
tmeannotations <- annotations[isEpithelial == FALSE]

cluster_variance <- function(dt, annot_dt, group_label) {
  
  counts <- dt[, .N, by = .(PatientID, ImageID, ClusterID)]
  totals <- counts[, .(Total = sum(N)), by = .(PatientID, ImageID)]
  counts <- merge(counts, totals, by = c("PatientID", "ImageID"))
  
  counts[, Proportion := N / Total]
  
  variation <- counts[, .(VarProp = var(Proportion, na.rm = TRUE)),
                      by = .(PatientID, ClusterID)]
  variation[, ClusterID := as.character(ClusterID)]
  
  annot <- copy(annot_dt)
  annot[, ClusterID := as.character(ClusterID)]
  
  variation <- merge(
    variation,
    annot[, .(ClusterID, Colour, BackupTeXClusterLabel)],
    by = "ClusterID",
    all.x = TRUE
  )
  
  variation[, Group := group_label]
  
  cols <- setNames(annot$Colour, annot$ClusterID)
  variation[, Fill := cols[ClusterID]]
  
  variation
}

epi_variation <- cluster_variance(
  cellsclin[isEpithelial == TRUE],
  epiannotations,
  "Epithelial"
)

tme_variation <- cluster_variance(
  cellsclin[isEpithelial == FALSE],
  tmeannotations,
  "Microenvironment"
)

combined_variation <- rbind(epi_variation, tme_variation)

p <- ggplot(combined_variation,
            aes(x = factor(ClusterID), y = VarProp, fill = Fill)) +
  geom_boxplot(outlier.size = 0.1, outlier.alpha = 0.2) +
  facet_wrap(~ Group, scales = "free_x") +
  labs(x = NULL, y = TeX("ROI Variation\n(Variance of Proportion)")) +
  theme_classic() +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(0, 0.3),
    breaks = c(0, 0.3),
    labels = c("0", "0.3")
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                               color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 15),
    legend.position = "none",
    plot.margin = unit(c(20, 20, 20, 20), "pt"),
    strip.text.x = element_text(size = 15, color = "black", hjust = 0.5),
    strip.background = element_blank()
  ) +
  scale_fill_identity() +
  scale_x_discrete(labels = function(x)
    sapply(x, function(clID)
      TeX(combined_variation[ClusterID == clID,
                             unique(BackupTeXClusterLabel)])))

ggsave(
  p,
  filename = file.path(outdir, "ROIVariationPlot.pdf"),
  units = "in",
  width = 7.75,
  height = 4
)
