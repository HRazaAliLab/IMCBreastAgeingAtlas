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

avg_axis <- cellsclin[, .(AvgAxisRatio = mean(AxisRatio, na.rm = TRUE)),
                      by = .(PatientID, ClusterID, AgeGroup)]

summary_p <- rbindlist(lapply(unique(avg_axis$ClusterID), function(cluster) {
  cd <- avg_axis[ClusterID == cluster]
  tt <- t.test(AvgAxisRatio ~ AgeGroup, data = cd, var.equal = FALSE)
  data.table(Cluster = as.character(cluster), PValue = tt$p.value)
}))
summary_p[, adjPValue := p.adjust(PValue, method = "BH")]

createAxisRatioBoxplot <- function(patientData, summaryData, cluster_id, ymin, ymax) {
  dt <- patientData[ClusterID == cluster_id]
  dt[, AgeGroup := factor(AgeGroup, levels = c("Premenopausal", "Postmenopausal"))]
  dt <- na.omit(dt, cols = "AgeGroup")
  
  adjPVal <- summaryData[Cluster == as.character(cluster_id), adjPValue]
  label <- annotations[ClusterID == cluster_id, BackupLabel]
  
  ggplot(dt, aes(x = AgeGroup, y = AvgAxisRatio, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.2) +
    scale_fill_manual(values = c("Premenopausal" = "gray", "Postmenopausal" = "dimgray")) +
    theme_classic() +
    scale_x_discrete(labels = c("pre", "post")) +
    scale_y_continuous(limits = c(ymin, ymax),
                       breaks = c(ymin, floor(ymax / 0.05) * 0.05),
                       expand = expand_scale(mult = c(0, 0.08))) +
    labs(x = NULL, y = NULL, title = TeX(mkTeX(label, ignoreChar="/"))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 17),
      plot.title = element_text(size = 18, hjust = 0.5),
      legend.position = "none",
      plot.margin = unit(c(10,0,20,0), "pt")
    ) +
    annotate(
      "text", x = 1.5, y = ymax,
      label = TeX(mkEnumPower(paste0("p=", format_custom_pval(adjPVal)))),
      size = 6, color = "black", hjust = 0.5
    )
}

axisratioplot1 <- createAxisRatioBoxplot(avg_axis, summary_p, 25, ymin = 1.3, ymax = 1.71)
axisratioplot2 <- createAxisRatioBoxplot(avg_axis, summary_p, 15, ymin = 1.1, ymax = 1.5)

combined <- axisratioplot1 / axisratioplot2 + plot_layout(ncol = 1)

ylab <- ggplot() +
  labs(y = "major/minor axis ratio") +
  theme_classic() +
  theme(plot.margin = unit(c(10,0,0,0), "pt"),
        axis.title.y = element_text(color = "black", size = 17)) +
  guides(x = "none", y = "none")

combined <- ylab + combined + plot_layout(widths = c(1, 1000))

ggsave(file.path(outdir, "axisRatioAgeBoxplots.pdf"),
       combined, width = 2.5, height = 6.5, units = "in")
