library(here)
source(here("code", "header.R"))

outdir <- here("scratch/supp/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
cells <- merge(cells, clinical[, .(ImageID, PatientID)], by = "ImageID", all.x = TRUE)

epithelial <- cells[isEpithelial == TRUE & !is.na(PatientID)]

# ER+ epithelial cells with Ki67 status
erpos_cells <- epithelial[isERPos == TRUE]
erpos_cells[, Ki67Status := ifelse(isKi67Pos, "Ki67+", "Ki67-")]

p_cell <- ggplot(erpos_cells,
                 aes(x = Ki67Status, y = ER, fill = Ki67Status)) +
  geom_boxplot(outlier.size = 0.2, outlier.stroke = 0.2) +
  theme_classic() +
  scale_fill_manual(values = c("Ki67+" = "white", "Ki67-" = "white")) +
  labs(x = NULL, y = "Normalized ER intensity") +
  scale_y_continuous(
    breaks = c(0, 6),
    limits = c(0, 7),
    expand = c(0, 0)
  ) +
  scale_x_discrete(labels = c("", "")) +
  theme(
    axis.title.y = element_text(size = 12, color = "black", vjust = 2),
    axis.text.x  = element_text(size = 11, color = "black", vjust = 2),
    axis.text.y  = element_text(size = 11, color = "black"),
    legend.position = "none"
  )

ggsave(
  plot = p_cell,
  filename = file.path(outdir, "ERKi67CellBoxplot.pdf"),
  units = "in",
  width = 1.5,
  height = 2.5
)
