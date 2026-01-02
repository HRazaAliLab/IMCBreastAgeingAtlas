library(here)
source(here::here("code/header.R"))

outdir <- here("scratch/supp/Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
cells <- merge(cells, clinical[, .(ImageID, PatientID)], by = "ImageID", all.x = TRUE)

epithelial <- cells[isEpithelial == TRUE]

markers <- c("FOXA1", "GATA3", "ER", "PR", "AR", "HER2")

need_is <- c(
  "isFOXA1Pos","isGATA3Pos","isERPos","isPRPos","isARPos","isHER2Pos"
)

epithelial_numeric <- epithelial[
  ,
  c(list(PatientID = PatientID), lapply(.SD, as.integer)),
  .SDcols = need_is,
  by = PatientID
]

setnames(
  epithelial_numeric,
  old = c("isFOXA1Pos", "isGATA3Pos", "isERPos", "isPRPos", "isARPos", "isHER2Pos"),
  new = markers
)

epithelial_numeric[, combination := apply(.SD, 1, function(row) {
  paste(names(row)[row == 1], collapse = ",")
}), .SDcols = markers]

all_combinations <- do.call(
  CJ,
  setNames(rep(list(c(0L, 1L)), length(markers)), markers)
)

all_combinations <- all_combinations[rowSums(all_combinations) > 0]

all_combinations[, combination := apply(.SD, 1, function(x) {
  present <- markers[x == 1]
  if (length(present) > 0) paste(present, collapse = ",") else ""
}), .SDcols = markers]

all_combinations <- all_combinations[combination != ""]

counts_combinations <- epithelial_numeric[, .N, by = combination][order(-N)]
counts_combinations[, combination := factor(combination, levels = combination)]

counts_combinations_all <- merge(all_combinations, counts_combinations, by = "combination", all.x = TRUE)
counts_combinations_all[is.na(N), N := 0L]

counts_combinations_all[, num_markers := rowSums(.SD), .SDcols = markers]
counts_combinations_all[, sort_key := sapply(strsplit(combination, ",", fixed = TRUE), function(x) {
  paste(match(x, markers), collapse = ",")
})]
setorder(counts_combinations_all, num_markers, sort_key)

counts_combinations_all[, combo_id := factor(combination, levels = combination)]

heatmap_data <- melt(
  counts_combinations_all,
  id.vars = c("combination", "N", "num_markers", "combo_id"),
  measure.vars = markers,
  variable.name = "marker",
  value.name = "present"
)

heatmap_data[, marker := factor(marker, levels = rev(markers))]

heatmap_plot <- ggplot(heatmap_data, aes(x = combo_id, y = marker)) +
  geom_tile(aes(fill = factor(present)), color = "white") +
  scale_fill_manual(values = c("0" = "white", "1" = "dimgray"), name = "Presence") +
  theme_classic() +
  labs(x = NULL, y = NULL) +
  theme(
    axis.ticks.x = element_blank(),
    axis.text.y  = element_text(size = 10, color = "black"),
    axis.text.x  = element_blank(),
    axis.title   = element_text(color = "black", size = 14),
    panel.grid   = element_blank(),
    legend.position = "none"
  ) +
  scale_y_discrete(expand = c(0, 0))

counts_df <- counts_combinations_all[, .(combination, N)]
counts_df[, combination := factor(combination, levels = levels(counts_combinations_all$combo_id))]

bar_chart <- ggplot(counts_df, aes(x = combination, y = N)) +
  geom_col(width = 0.6, fill = "dimgray") +
  theme_classic() +
  labs(x = NULL, y = "Counts") +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 40000), labels = c("0", "40,000")) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0, "cm"),
    axis.text.y = element_text(color = "black", size = 10),
    axis.title  = element_text(color = "black", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

combined_plot <- bar_chart / heatmap_plot + plot_layout(heights = c(4, 4))

ggsave(
  filename = file.path(outdir, "fullCombinationPlot.pdf"),
  plot = combined_plot,
  width = 13, height = 2.9, units = "in"
)

message("Saved: ", file.path(outdir, "fullCombinationPlot.pdf"))
