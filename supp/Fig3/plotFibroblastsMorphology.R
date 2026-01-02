library(here)
source(here::here("code/header.R"))

outdir <- here("scratch/supp/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()

clinical <- getClinical()
cells <- merge(cells, clinical[, .(ImageID, PatientID)], by = "ImageID", all.x = TRUE)

tme <- cells[isEpithelial == FALSE]

# fibroblasts = cluster 25
tme[, ClusterType := ifelse(as.character(ClusterID) == "25", "fibroblasts", "other microenvironment")]

tme[, AxisRatio := MajorAxisLength / MinorAxisLength]

avg_ratio <- tme[, .(Avg_AxisRatio = mean(AxisRatio, na.rm = TRUE)),
                 by = .(PatientID, ClusterType)]

fibro <- avg_ratio[ClusterType == "fibroblasts"]
other <- avg_ratio[ClusterType == "other microenvironment"]
paired <- merge(fibro, other, by = "PatientID", suffixes = c(".fibro", ".other"))

test_ratio <- wilcox.test(
  paired$Avg_AxisRatio.fibro,
  paired$Avg_AxisRatio.other,
  paired = TRUE
)

avg_ratio[, ClusterType := factor(
  ClusterType,
  levels = c("fibroblasts", "other microenvironment"),
  labels = c("fibroblasts", "other TME")
)]

plot_ratio <- ggplot(avg_ratio, aes(x = ClusterType, y = Avg_AxisRatio, fill = ClusterType)) +
  geom_boxplot(outlier.size = 0.2) +
  scale_fill_manual(values = c("fibroblasts" = "#CDAABA", "other TME" = "gray")) +
  theme_classic() +

scale_y_continuous(limits = c(1.3, 1.7), breaks = c(1.3, 1.7)) +
  labs(x = NULL, y = "major / minor axis ratio") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.text = element_text(color = "black"),
    axis.title.y = element_text(vjust = -1),
    legend.position = "none"
  ) +
  annotate(
    "text",
    x = 1.5, y = 1.68,
    label = TeX(paste0("p = ", mkEnumPower(format_custom_pval(test_ratio$p.value)))),
    parse = TRUE,
    size = 5
  )

pdf(file.path(outdir, "fibroblastsMorphology.pdf"), width = 2.5, height = 5)
print(plot_ratio)
dev.off()
