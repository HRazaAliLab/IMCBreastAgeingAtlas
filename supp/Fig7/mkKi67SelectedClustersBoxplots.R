library(here)
source(here("code", "header.R"))

outdir <- here("scratch/supp/Fig7")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

annotations[, ClusterID := as.integer(as.character(ClusterID))]

cellsclin <- merge(cells, clinical, by = "ImageID")
cellsclin[, ClusterID := as.integer(as.character(ClusterID))]

positiveKi67Fractions <- getPositiveFractions(d = cellsclin, clusterColumn = "ClusterID", idCols = "PatientID")
avg_ki67 <- merge(positiveKi67Fractions, unique(clinical[, .(Age, PatientID)]), by = "PatientID", all.x = TRUE)
avg_ki67[, AgeGroup := ifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
setnames(avg_ki67, "positiveFraction", "AvgKi67Proportion")

summary_results_ki67 <- avg_ki67[
  , {
    cluster_data <- na.omit(.SD)
    test_result <- wilcox.test(AvgKi67Proportion ~ AgeGroup, data = cluster_data, exact = FALSE)
    cliff_d_result <- cliff.delta(AvgKi67Proportion ~ AgeGroup, data = cluster_data)
    list(
      PValue = test_result$p.value,
      CliffsDelta = cliff_d_result$estimate,
      CID_Lower = cliff_d_result$conf.int[1],
      CID_Upper = cliff_d_result$conf.int[2]
    )
  },
  by = ClusterID
]

summary_results_ki67[, Cluster := as.integer(ClusterID)]
summary_results_ki67[, adjPValue := p.adjust(PValue, method = "BH")]

summary_results_typed <- merge(summary_results_ki67, annotations, by.x = "Cluster", by.y = "ClusterID", all.x = TRUE)

get_adj_p <- function(summaryData, cluster_id) {
  x <- summaryData[summaryData$Cluster == cluster_id, "adjPValue"]
  if (is.data.frame(x) || data.table::is.data.table(x)) x <- x[[1]]
  as.numeric(x[1])
}

createKi67Boxplot <- function(patientData, summaryData, cluster_id, ymax) {
  cluster_specific_data <- patientData[ClusterID == cluster_id]
  cluster_specific_data[, AgeGroup := factor(AgeGroup, levels = c("Premenopausal", "Postmenopausal"))]
  cluster_specific_data <- na.omit(cluster_specific_data, cols = "AgeGroup")
  
  adjPVal <- get_adj_p(summaryData, cluster_id)
  label <- annotations[ClusterID == cluster_id, BackupTeXClusterLabel]
  
  p_label <- if (!is.na(adjPVal) && adjPVal < 0.001) {
    mkEnumPower(formatC(adjPVal, format = "e", digits = 0))
  } else if (!is.na(adjPVal)) {
    sprintf("%.3f", adjPVal)
  } else {
    "NA"
  }
  
  ggplot(cluster_specific_data, aes(x = AgeGroup, y = AvgKi67Proportion, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("Premenopausal" = "gray", "Postmenopausal" = "dimgray")) +
    theme_classic() +
    scale_x_discrete(labels = c("pre", "post")) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = c(0, floor(ymax / 0.05) * 0.05),
      expand = expand_scale(mult = c(0, 0.1)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title = TeX(label)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black", size = 14),
      axis.title.y = element_text(vjust = -3, hjust = 0.5),
      plot.title = element_text(size = 15, hjust = 0.5, vjust = -0.5),
      legend.position = "none"
    ) +
    annotate("text", x = 1.5, y = ymax,
             label = TeX(paste0("p=", p_label)),
             size = 5.5, color = "black", hjust = 0.5)
}

# ---- chosen clusters ----
ki67plot1 <- createKi67Boxplot(avg_ki67, summary_results_typed, 5,  ymax = 0.44)
ki67plot2 <- createKi67Boxplot(avg_ki67, summary_results_typed, 16, ymax = 0.20)
ki67plot3 <- createKi67Boxplot(avg_ki67, summary_results_typed, 21, ymax = 0.06)
ki67plot4 <- createKi67Boxplot(avg_ki67, summary_results_typed, 25, ymax = 0.05)

ki67combined <- (ki67plot1 + ki67plot2 + ki67plot3 + ki67plot4) + plot_layout(ncol = 4)

ki67blanklabelploty <- ggplot() +
  labs(y = "proliferative fraction") +
  theme_classic() +
  theme(
    plot.margin = margin(0,0,0,0, unit = "cm"),
    axis.title.y = element_text(color = "black", size = 15)
  ) +
  guides(x = "none", y = "none")

ki67combined <- ki67blanklabelploty + ki67combined + plot_layout(widths = c(1, 1000))

ggsave(
  file.path(outdir, "Ki67SelectedClustersBoxplots.pdf"),
  ki67combined,
  width = 7, height = 2.7, units = "in"
)
