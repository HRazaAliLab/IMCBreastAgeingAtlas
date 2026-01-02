library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellsclin[, AgeGroup := fifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
cellsclin[, ClusterID := as.integer(ClusterID)]

avg_area <- cellsclin[, .(AvgArea = mean(CellArea, na.rm = TRUE)), by = .(PatientID, ClusterID, AgeGroup)]

summary_p <- rbindlist(lapply(unique(avg_area$ClusterID), function(cluster) {
  cd <- avg_area[ClusterID == cluster]
  w <- wilcox.test(AvgArea ~ AgeGroup, data = cd, exact = FALSE)
  data.table(Cluster = as.character(cluster), PValue = w$p.value)
}))
summary_p[, adjPValue := p.adjust(PValue, method = "BH")]

createAreaBoxplot <- function(patientData, summaryData, cluster_id, ymin, ymax) {
  dt <- patientData[ClusterID == cluster_id]
  dt[, AgeGroup := factor(AgeGroup, levels = c("Premenopausal", "Postmenopausal"))]
  dt <- na.omit(dt, cols = "AgeGroup")
  
  adjPVal <- summaryData[Cluster == as.character(cluster_id), adjPValue]
  label <- annotations[ClusterID == cluster_id, BackupLabel]
  
  ggplot(dt, aes(x = AgeGroup, y = AvgArea, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("Premenopausal" = "gray", "Postmenopausal" = "dimgray")) +
    theme_classic() +
    scale_x_discrete(labels = c("pre", "post")) +
    scale_y_continuous(limits = c(ymin, ymax),
                       breaks = c(ymin, floor(ymax/5)*5),
                       expand = expand_scale(mult = c(0, 0.08))) +
    labs(x = NULL, y = NULL, title = TeX(mkTeX(label, ignoreChar="/"))) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 14),
      plot.title = element_text(size = 14, hjust = 0.5),
      legend.position = "none",
      plot.margin = unit(c(10,0,10,0), "pt")
    ) +
    annotate("text", x = 1.5, y = ymax,
             label = TeX(mkEnumPower(paste0("p=", format_custom_pval(adjPVal)))),
             size = 4.5, color = "black", hjust = 0.5)
}

areaplot1 <- createAreaBoxplot(avg_area, summary_p, 25, ymin = 15, ymax = 70)
areaplot2 <- createAreaBoxplot(avg_area, summary_p, 8,  ymin = 20, ymax = 100)

combined <- areaplot1 / areaplot2 + plot_layout(ncol = 1)

ylab <- ggplot() + labs(y = "average area per patient (µm²)") +
  theme_classic() +
  theme(plot.margin = unit(c(0,0,0,0), "pt"),
        axis.title.y = element_text(color = "black", size = 16)) +
  guides(x = "none", y = "none")

combined <- ylab + combined + plot_layout(widths = c(1, 1000))

ggsave(file.path(outdir, "cellAreaAgeBoxplots.pdf"),
       combined, width = 2, height = 6.25, units = "in")
