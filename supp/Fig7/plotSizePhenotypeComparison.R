library(here)
source(here("code", "header.R"))

outdir <- here("scratch/supp/Fig7")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells       <- getCells()
clinical    <- getClinical()
annotations <- getCellClusters()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellsclin[, AgeGroup := fifelse(Age >= 50, "Postmenopausal", "Premenopausal")]
cellsclin[, ClusterID := as.integer(ClusterID)]

avg_area <- cellsclin[, .(
  AvgArea = mean(CellArea, na.rm = TRUE)
), by = .(PatientID, ClusterID, AgeGroup)]

avg_area_typed <- merge(avg_area, annotations, by = "ClusterID", all.x = TRUE)
avg_area_typed[, Type := factor(Type, levels = c("Epithelial", "Immune", "Stromal"))]

avg_area_typed[, ClusterID_f := factor(
  ClusterID,
  levels = unique(ClusterID[order(Type, PrintOrder)])
)]

labels <- setNames(annotations$BackupLabel, annotations$ClusterID)
colors <- setNames(annotations$Colour, annotations$ClusterID)

type_labels <- c("Epithelial", "Immune", "Stromal")
label_mapping <- setNames(type_labels, type_labels)

p1 <- ggplot(avg_area_typed, aes(x = ClusterID_f, y = AvgArea, fill = factor(ClusterID))) +
  geom_boxplot(outlier.size = 0.1) +
  labs(x = NULL, y = TeX("average area ($\\mu$m$^2$)")) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  scale_y_continuous(breaks = c(0, 60, 120), limits = c(0, 120), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = 13),
    axis.text.y = element_text(color = "black", size = 13),
    axis.title.y = element_text(color = "black", size = 13),
    legend.position = "none",
    plot.margin = unit(c(40, 0, 0, 5), "pt"),
    strip.text.x = element_text(size = 13, color = "black", hjust = 0.5),
    strip.background = element_blank()
  ) +
  facet_grid(cols = vars(Type), scales = "free_x", space = "free",
             labeller = labeller(Type = label_mapping)) +
  scale_x_discrete(labels = sapply(labels[levels(avg_area_typed$ClusterID_f)],
                                   function(x) TeX(mkTeX(x, ignoreChar="/"))))

ggsave(
  p1,
  filename = file.path(outdir, "sizePhenotypeComparison.pdf"),
  units = "in", width = 7.5, height = 3.2
)
