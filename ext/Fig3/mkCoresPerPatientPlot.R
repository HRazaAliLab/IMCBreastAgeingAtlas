# ---------------------- CORES PER PATIENT ----------------------
library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()

cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

cellsclin[, candidateBase := gsub("_split", "", ImageID)]
cellsclin[, baseImageID := ifelse(any(ImageID == candidateBase), candidateBase, ImageID),
          by = .(PatientID, candidateBase)]

roi_counts <- cellsclin[, .(nROIs = uniqueN(baseImageID)), by = PatientID]
roi_summary <- roi_counts[, .(Count = .N), by = nROIs]

p <- ggplot(roi_summary, aes(x = factor(nROIs), y = Count)) +
  geom_bar(stat = "identity", fill = "darkgrey", width = 0.95) +
  labs(x = "Cores per patient", y = "Number of patients") +
  theme_classic() +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200), breaks = c(0, 200)) +
  theme(
    axis.text = element_text(color = "black", size = 15),
    plot.margin = unit(c(10, 10, 10, 10), "pt")
  )

ggsave(p, filename = file.path(outdir, "CoresPerPatientPlot.pdf"),
       units = "in", width = 3, height = 3)
