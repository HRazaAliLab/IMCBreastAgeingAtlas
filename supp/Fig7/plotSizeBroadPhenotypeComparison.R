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

cellsclin[, CellPhenotype_BroadCategory := NA_character_]
cellsclin[isEpithelial == TRUE, CellPhenotype_BroadCategory := "Epithelial"]

for (cid in unique(cellsclin$ClusterID)) {
  if (cid %in% annotations$ClusterID) {
    type <- annotations[ClusterID == cid, Type][1]
    cellsclin[isEpithelial == FALSE & ClusterID == cid, CellPhenotype_BroadCategory := type]
  }
}

avg_area_broad <- cellsclin[, .(
  AvgArea = mean(CellArea, na.rm = TRUE)
), by = .(PatientID, CellPhenotype_BroadCategory, AgeGroup, isKi67Pos)]

avg_area_broad[, CellPhenotype_BroadCategory :=
                 factor(CellPhenotype_BroadCategory, levels = c("Epithelial", "Immune", "Stromal"))
]

lighten_hex <- function(cols, amount = 0.7) {
  stopifnot(is.numeric(amount), amount >= 0, amount <= 1)
  out <- vapply(cols, function(cl) {
    rgb0 <- grDevices::col2rgb(cl) / 255
    rgb1 <- rgb0 + (1 - rgb0) * amount
    grDevices::rgb(rgb1[1], rgb1[2], rgb1[3])
  }, character(1))
}

broad_cols <- getNormalBreastProjectColours()$BroadPhenotype
lighter <- lighten_hex(broad_cols, amount = 0.7)

combined_colors <- c(
  Epithelial_TRUE  = lighter["Epithelial"],
  Epithelial_FALSE = broad_cols["Epithelial"],
  Immune_TRUE      = lighter["Immune"],
  Immune_FALSE     = broad_cols["Immune"],
  Stromal_TRUE     = lighter["Stromal"],
  Stromal_FALSE    = broad_cols["Stromal"]
)
names(combined_colors) <- sub("\\..*$", "", names(combined_colors))

avg_area_broad[, CombinedCategory := paste0(CellPhenotype_BroadCategory, "_", isKi67Pos)]

results <- avg_area_broad[, .(
  p_value = wilcox.test(AvgArea[isKi67Pos == TRUE], AvgArea[isKi67Pos == FALSE])$p.value
), by = CellPhenotype_BroadCategory]

label_data <- results[, .(
  CellPhenotype_BroadCategory,
  x = 1.5,
  y = 125,
  label = mkEnumPower(format_custom_pval(p_value))
)]

p1b <- ggplot(avg_area_broad, aes(x = factor(isKi67Pos), y = AvgArea, fill = CombinedCategory)) +
  geom_boxplot(outlier.size = 0.1) +
  labs(x = NULL, y = TeX("mean area ($\\mu$m$^2$)")) +
  theme_classic() +
  scale_fill_manual(values = combined_colors) +
  scale_y_continuous(breaks = c(0, 60, 120), limits = c(0, 130), expand = c(0, 0)) +
  theme(
    axis.text.x = element_text(angle = 0, color = "black", size = 18),
    strip.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 18),
    legend.position = "none",
    plot.margin = unit(c(40, 0, 0, 5), "pt")
  ) +
  scale_x_discrete(labels = c("TRUE" = "+", "FALSE" = "-")) +
  facet_wrap(~ CellPhenotype_BroadCategory, scales = "free_x") +
  geom_text(
    data = label_data,
    aes(x = x, y = y, label = TeX(paste0(label), output = "character")),
    parse = TRUE, size = 4.5, color = "black", hjust = 0.5,
    inherit.aes = FALSE
  )

ggsave(
  p1b,
  filename = file.path(outdir, "sizeBroadPhenotypeComparison.pdf"),
  units = "in", width = 3.5, height = 3.2
)
