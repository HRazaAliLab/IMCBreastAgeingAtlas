library(here)
source(here("code", "header.R"))

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID")

outdir <- here("scratch/ext/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

epithelial <- cellsclin[isEpithelial == TRUE]

epithelial[, isPrimaryPop := (isFOXA1Pos & isGATA3Pos & isERPos & isPRPos & isARPos & !isHER2Pos)]
epithelial[, isSecondaryPop := (!isFOXA1Pos & !isGATA3Pos & isERPos & !isPRPos & !isARPos & !isHER2Pos)]

primary_color <- "#4D4D4D"
  library(here)
source(here("code", "header.R"))

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID")

outdir <- here("scratch/ext/Fig4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

epithelial <- cellsclin[isEpithelial == TRUE]

epithelial[, isPrimaryPop := (isFOXA1Pos & isGATA3Pos & isERPos & isPRPos & isARPos & !isHER2Pos)]
epithelial[, isSecondaryPop := (!isFOXA1Pos & !isGATA3Pos & isERPos & !isPRPos & !isARPos & !isHER2Pos)]

primary_color <- "#4D4D4D"
secondary_color <- "#8B4513"
    
  # ---------------- proliferative fractions ----------------
aggregated_data <- epithelial[, .(
  TotalPrimaryPop = sum(isPrimaryPop, na.rm = TRUE),
  TotalPrimaryKi67Pos = sum(isPrimaryPop & isKi67Pos, na.rm = TRUE),
  TotalSecondaryPop = sum(isSecondaryPop, na.rm = TRUE),
  TotalSecondaryKi67Pos = sum(isSecondaryPop & isKi67Pos, na.rm = TRUE)
), by = .(PatientID, Age)]

aggregated_data[, PercentKi67Primary := (TotalPrimaryKi67Pos / TotalPrimaryPop)]
aggregated_data[, PercentKi67Secondary := (TotalSecondaryKi67Pos / TotalSecondaryPop)]
aggregated_data[is.na(PercentKi67Primary), PercentKi67Primary := 0]
aggregated_data[is.na(PercentKi67Secondary), PercentKi67Secondary := 0]

paired <- aggregated_data[, .(PercentKi67Primary, PercentKi67Secondary), by = PatientID]
wilcox_test <- wilcox.test(paired$PercentKi67Primary, paired$PercentKi67Secondary, paired = TRUE)
pval <- wilcox_test$p.value
formattedPVal <- if (pval > 1e-4) formatC(pval, format = "f", digits = 2) else mkEnumPower(formatC(pval, format = "e", digits = 0))

long_agg <- melt(
  aggregated_data,
  id.vars = c("PatientID", "Age"),
  measure.vars = c("PercentKi67Primary", "PercentKi67Secondary"),
  variable.name = "Population",
  value.name = "PercentKi67"
)

long_agg[, Population := factor(
  Population,
  levels = c("PercentKi67Primary", "PercentKi67Secondary"),
  labels = c("Primary (High Hormone)", "Secondary (ER+/FOXA1-/GATA3-)")
)]

p_prolif <- ggplot(long_agg, aes(x = Population, y = PercentKi67, fill = Population)) +
  geom_boxplot(outlier.size = 0.1, outlier.color = "black") +
  scale_fill_manual(values = c(
    "Primary (High Hormone)" = primary_color,
    "Secondary (ER+/FOXA1-/GATA3-)" = secondary_color
  )) +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 14),
    legend.position = "none"
  ) +
  scale_y_continuous(limits = c(0, 0.3), breaks = c(0, 0.3), labels = c("0", "0.3"),
                     expand = expansion(mult = c(0, 0.1))) +
  annotate("text", x = 1.5, y = 0.3,
           label = TeX(paste0("p=", formattedPVal)),
           size = 3.5, hjust = 0.5)

ggsave(p_prolif, filename = file.path(outdir, "proliferativeFractionsPopulations.pdf"),
       units = "in", width = 1.5, height = 1.5)


# ---------------- area + axis ratio ----------------

calculate_and_plot <- function(data, measure_col, y_limits, y_breaks, out_fn) {
  
  if (!(measure_col %in% names(data))) stop(paste0(measure_col, " not found in data"))

  y_label <- switch(
    measure_col,
    "CellArea"   = "Cell area",
    "AxisRatio"  = "Axis ratio (major / minor)",
    measure_col  # fallback
  )
  
  aggregated_data <- data[, .(
    TotalPrimaryPop       = sum(isPrimaryPop, na.rm = TRUE),
    TotalPrimaryMeasure   = sum(get(measure_col)[isPrimaryPop], na.rm = TRUE),
    TotalSecondaryPop     = sum(isSecondaryPop, na.rm = TRUE),
    TotalSecondaryMeasure = sum(get(measure_col)[isSecondaryPop], na.rm = TRUE)
  ), by = .(PatientID, Age)]
  
  aggregated_data[, MeanPrimary   := TotalPrimaryMeasure / TotalPrimaryPop]
  aggregated_data[, MeanSecondary := TotalSecondaryMeasure / TotalSecondaryPop]
  
  aggregated_data[is.na(MeanPrimary),   MeanPrimary := 0]
  aggregated_data[is.na(MeanSecondary), MeanSecondary := 0]
  
  paired_data <- aggregated_data[, .(MeanPrimary, MeanSecondary), by = PatientID]
  paired_data <- na.omit(paired_data)
  
  w <- wilcox.test(paired_data$MeanPrimary, paired_data$MeanSecondary, paired = TRUE)
  pval <- w$p.value
  formattedPVal <- if (pval > 1e-4) formatC(pval, format="f", digits=1) else mkEnumPower(formatC(pval, format="e", digits=0))
  
  long_aggregated_data <- melt(
    aggregated_data,
    id.vars = c("PatientID", "Age"),
    measure.vars = c("MeanPrimary", "MeanSecondary"),
    variable.name = "Population",
    value.name = "MeanMeasure"
  )
  
  long_aggregated_data[, Population := factor(
    Population,
    levels = c("MeanPrimary", "MeanSecondary"),
    labels = c("Primary (High Hormone)", "Secondary (ER+/FOXA1-/GATA3-)")
  )]
  
  p <- ggplot(long_aggregated_data, aes(x = Population, y = MeanMeasure, fill = Population)) +
    geom_boxplot(outlier.size = 0.1, outlier.color = "black") +
    scale_fill_manual(values = c(
      "Primary (High Hormone)" = primary_color,
      "Secondary (ER+/FOXA1-/GATA3-)" = secondary_color
    )) +
    theme_classic() +
    labs(y = y_label) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 12),
      axis.title.y = element_text(size = 12),
      legend.position = "none"
    ) +
    scale_y_continuous(limits = y_limits, breaks = y_breaks, expand = expansion(mult = c(0, 0.1))) +
    annotate("text", x = 1.5, y = y_limits[2],
             label = TeX(paste0("p=", formattedPVal)),
             size = 4, hjust = 0.5)
  
  ggsave(p, filename = file.path(outdir, out_fn), units = "in", width = 2, height = 2)
}


# Plot area (CellArea)
calculate_and_plot(
  epithelial,
  measure_col = "CellArea",
  y_limits = c(30, 140),
  y_breaks = c(30, 140),
  out_fn = "areaPlot.pdf"
)

# Plot axis ratio (MajorAxisLength / MinorAxisLength)
epithelial[, AxisRatio := MajorAxisLength / MinorAxisLength]

calculate_and_plot(
  epithelial,
  measure_col = "AxisRatio",
  y_limits = c(1.2, 1.7),
  y_breaks = c(1.2, 1.7),
  out_fn = "axisRatioPlot.pdf"
)
