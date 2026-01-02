library(here)
source(here("code", "header.R"))

outDir <- here("scratch/ext/Fig4")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
clinical <- read_parquet_adt(here('data/derived/clinical.parquet'))

cells <- getCells()
clinical <- getClinical()

cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age, TissueArea)], by = "ImageID", all.x = TRUE)

microenvironment <- cellsclin[isEpithelial == FALSE]

# ----- NK cell investigation -----
cd45_counts <- microenvironment[, .(count = .N, cd45_positive = sum(isCD45Pos)), by = isCD45Pos]
total_cells <- nrow(microenvironment)
cd45_positive_count <- cd45_counts[isCD45Pos == TRUE, count]
cd45_positive_percentage <- (cd45_positive_count / total_cells) * 100

cat("Total cells:", total_cells, "\n")
cat("CD45+ cells count:", cd45_positive_count, "\n")
cat("CD45+ cells percentage:", cd45_positive_percentage, "%\n\n")

cd45_positive_population <- microenvironment[isCD45Pos == TRUE]
cd56_in_cd45_positive_count <- cd45_positive_population[, sum(isCD56Pos)]
cd56_in_cd45_positive_percentage <- (cd56_in_cd45_positive_count / cd45_positive_count) * 100

cat("CD56+ in CD45+ population count:", cd56_in_cd45_positive_count, "\n")
cat("CD56+ in CD45+ population percentage:", cd56_in_cd45_positive_percentage, "%\n\n")

cd45_negative_population <- microenvironment[isCD45Pos == FALSE]
cd56_in_cd45_negative_count <- cd45_negative_population[, sum(isCD56Pos)]
cd45_negative_count <- nrow(cd45_negative_population)
cd56_in_cd45_negative_percentage <- (cd56_in_cd45_negative_count / cd45_negative_count) * 100

cat("CD45- cells count:", cd45_negative_count, "\n")
cat("CD56+ in CD45- population count:", cd56_in_cd45_negative_count, "\n")
cat("CD56+ in CD45- population percentage:", cd56_in_cd45_negative_percentage, "%\n")

plot_data <- data.frame(
  Category = c("%CD45+ of all TME cells", "%CD56+ of CD45+ cells", "%CD56+ of CD45- cells"),
  Percentage = c(cd45_positive_percentage, cd56_in_cd45_positive_percentage, cd56_in_cd45_negative_percentage),
  Count = c(cd45_positive_count, cd56_in_cd45_positive_count, cd56_in_cd45_negative_count)
)
plot_data$Category <- factor(plot_data$Category, levels = rev(unique(plot_data$Category)))

NKcounts <- ggplot(plot_data, aes(x = Percentage, y = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = scales::comma(Count)), hjust = -0.1, size = 5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Cell Population Percentages and Counts", x = "", y = "") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  theme(
    axis.title.x = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    axis.text.y = element_blank(),
    plot.title = element_blank(),
    plot.margin = unit(c(0, 40, 0, 0), units = "pt")
  )

ggsave(NKcounts, filename = file.path(outDir, "NKCounts.pdf"), units = "in", width = 2, height = 2)

nk_cells <- microenvironment[isCD45Pos == TRUE & isCD56Pos == TRUE, .N, by = PatientID]
setnames(nk_cells, "N", "NK_cell_counts")

tissueAreasPatients <- clinical[, .(totalTissueArea = sum(TissueArea, na.rm = TRUE)), by = PatientID]

clinical_pat <- unique(clinical[, .(PatientID, Age)])

merged_data <- merge(clinical_pat, nk_cells, by = "PatientID", all.x = TRUE)
merged_data <- merge(merged_data, tissueAreasPatients, by = "PatientID", all.x = TRUE)

merged_data[is.na(NK_cell_counts), NK_cell_counts := 0]

merged_data[, NK_cell_density := NK_cell_counts * 1000 * 1000 / totalTissueArea]

result <- merged_data[, .(PatientID, Age, totalTissueArea, NK_cell_counts, NK_cell_density)]
result <- unique(result)
result <- na.omit(result)

correlation_test <- cor.test(result$Age, result$NK_cell_density, method = "spearman")
correlation_coefficient <- correlation_test$estimate
p_value <- correlation_test$p.value
n <- nrow(result)

NKage <- ggplot(result, aes(x = Age, y = NK_cell_density)) +
  geom_point(size = 0.1, alpha = 0.3) +
  geom_smooth(method = "loess", se = FALSE, color = "dimgray") +
  labs(title = "", x = "age", y = TeX("density (NK cells/mm$^2$)")) +
  theme_classic() +
  theme(
    axis.title.x = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    axis.text.y = element_text(color = "black", size = 15)
  ) +
  scale_x_continuous(expand = c(0, 0), limits = c(15, 86), breaks = c(20, 50, 80)) +
  scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 5), breaks = c(0, 5)) +
  annotate("text", x = 17, y = 5 * 0.925,
           label = TeX(paste0("$\\rho$ = ", round(correlation_coefficient, 2))),
           size = 5.5, hjust = 0) +
  annotate("text", x = 17, y = 5 * 0.775,
           label = TeX(paste0("p-value = ", format_custom_pval(p_value))),
           size = 5.5, hjust = 0) +
  annotate("text", x = 17, y = 5 * 0.625,
           label = TeX(paste0("$\\textit{n} = ", n)),
           size = 5.5, hjust = 0)

ggsave(NKage, filename = file.path(outDir, "NKAge.pdf"), units = "in", width = 2.3, height = 3)
