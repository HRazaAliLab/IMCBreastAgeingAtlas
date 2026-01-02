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

# Patient-level mean ER per Ki67 status
patient_er_means <- erpos_cells[
  , .(meanER = mean(ER, na.rm = TRUE)),
  by = .(PatientID, Ki67Status)
]

# Keep only paired patients
validPatients <- patient_er_means[
  , .(nStatus = uniqueN(Ki67Status)),
  by = PatientID
][nStatus == 2, PatientID]

patient_er_means <- patient_er_means[PatientID %in% validPatients]

# Paired Wilcoxon
wide <- dcast(patient_er_means, PatientID ~ Ki67Status, value.var = "meanER")
pval <- wilcox.test(wide[["Ki67+"]], wide[["Ki67-"]], paired = TRUE)$p.value

formattedPVal <- if (pval > 1e-4) {
  formatC(pval, format = "f", digits = 2)
} else {
  mkEnumPower(formatC(pval, format = "e", digits = 0))
}

p_patient <- ggplot(
  patient_er_means,
  aes(x = Ki67Status, y = meanER, fill = Ki67Status)
) +
  geom_boxplot(outlier.size = 0.1) +
  theme_classic() +
  scale_fill_manual(values = c("Ki67+" = "white", "Ki67-" = "white")) +
  labs(
    x = NULL,
    y = "Mean ER per\npatient",
    subtitle = TeX(paste0("Paired Wilcoxon p=", formattedPVal))
  ) +
  scale_y_continuous(
    breaks = c(0, 6),
    limits = c(0, 7),
    expand = c(0, 0)
  ) +
  scale_x_discrete(labels = c("", "")) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11, color = "black"),
    axis.title.y = element_text(size = 12, color = "black", vjust = 2),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

ggsave(
  plot = p_patient,
  filename = file.path(outdir, "ERKi67PatientBoxplot.pdf"),
  units = "in",
  width = 1.7,
  height = 1.4
)
