# ---------------------- IMAGE-TO-IMAGE PEARSON CORRELATION ----------------------
library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig3")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

avg_roi_corr <- function(dt) {
  counts <- dt[, .N, by = .(PatientID, ImageID, ClusterID)]
  wide <- dcast(counts, PatientID + ImageID ~ ClusterID, value.var = "N", fill = 0)
  
  pheno_cols <- setdiff(names(wide), c("PatientID", "ImageID"))
  wide[, (pheno_cols) := lapply(.SD, as.numeric), .SDcols = pheno_cols]
  wide[, Total := rowSums(.SD), .SDcols = pheno_cols]
  wide[Total > 0, (pheno_cols) := lapply(.SD, function(x) x / Total), .SDcols = pheno_cols]
  
  patient_ids <- unique(wide$PatientID)
  
  rbindlist(lapply(patient_ids, function(pid) {
    patient_data <- wide[PatientID == pid]
    
    if (nrow(patient_data) < 2) {
      return(data.table(PatientID = pid, Similarity = NA_real_))
    }
    
    phenotype_cols <- setdiff(names(patient_data), c("PatientID", "ImageID", "Total"))
    prop_mat <- as.matrix(patient_data[, ..phenotype_cols])
    
    corr_mat <- cor(t(prop_mat), method = "pearson")
    lower_tri <- corr_mat[lower.tri(corr_mat)]
    
    data.table(PatientID = pid, Similarity = mean(lower_tri, na.rm = TRUE))
  }))
}

# Epithelial
pearson_epi <- avg_roi_corr(cellsclin[isEpithelial == TRUE])
pearson_epi[, `:=`(CellGroup = "Epithelial", Metric = "Pearson")]

# Microenvironment
pearson_tme <- avg_roi_corr(cellsclin[isEpithelial == FALSE])
pearson_tme[, `:=`(CellGroup = "Microenvironment", Metric = "Pearson")]

pearson_dt <- rbind(pearson_epi, pearson_tme)

# optional: add core counts (not used in plot)
cores_dt <- cellsclin[, .(TotalCores = uniqueN(ImageID)), by = PatientID]
pearson_dt <- merge(pearson_dt, cores_dt, by = "PatientID", all.x = TRUE)

pearson_dt[, Group := factor(CellGroup, levels = c("Epithelial", "Microenvironment"))]

p <- ggplot(pearson_dt, aes(x = Group, y = Similarity, fill = Group)) +
  geom_boxplot(outlier.size = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray", size = 0.5) +
  labs(x = "", y = "Image-to-image Pearson correlation") +
  theme_classic() +
  scale_y_continuous(limits = c(-1, 1), breaks = c(-1, 0, 1), expand = c(0.01, 0)) +
  scale_fill_manual(values = c("Epithelial" = "#6C7ACE",
                               "Microenvironment" = "#E8DDA6")) +
  theme(
    legend.position = "none",
    axis.text = element_text(color = "black", size = 20),
    plot.margin = unit(c(10, 10, 10, 10), "pt")
  )

ggsave(p, filename = file.path(outdir, "ImageCorrelationPlot.pdf"),
       units = "in", width = 5, height = 4)
