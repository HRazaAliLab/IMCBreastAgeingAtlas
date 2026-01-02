library(here)
source(here::here("code/header.R"))

outdir <- here("scratch/main/Fig2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()

clinical <- getClinical()
cells <- merge(cells, clinical[, .(ImageID, PatientID)], by = "ImageID", all.x = TRUE)

epithelial <- cells[isEpithelial == TRUE]

markers <- c("FOXA1", "GATA3", "ER", "PR", "AR", "HER2", "CK8/18", "CK5/14")
numcols <- 15

need_is <- c("isFOXA1Pos","isGATA3Pos","isERPos","isPRPos","isARPos","isHER2Pos","isCK8/18Pos","isCK5/14Pos","isKi67Pos")

epithelial_numeric <- epithelial[, c(list(PatientID = PatientID), lapply(.SD, as.integer)),
                                 .SDcols = need_is, by = PatientID]

setnames(epithelial_numeric,
         old = c("isFOXA1Pos", "isGATA3Pos", "isERPos", "isPRPos", "isARPos", "isHER2Pos", "isCK8/18Pos", "isCK5/14Pos"),
         new = markers)

epithelial_numeric[, combination := apply(.SD, 1, function(row) {
  paste(names(row)[row == 1], collapse = ",")
}), .SDcols = markers]

counts_combinations <- epithelial_numeric[, .N, by = combination][order(-N)]
counts_combinations[, combination := factor(combination, levels = combination)]

counts_combinations <- counts_combinations[-1][1:numcols]

bar_chart <- ggplot(counts_combinations, aes(x = combination, y = N)) +
  geom_col(width = 0.6, fill = 'darkgray') +
  theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    axis.title.y = element_text(color = "black", size = 16),
    plot.margin = unit(c(0,0,0,0), units = "cm"),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA)
  ) +
  coord_cartesian() +
  scale_x_discrete(expand = c(0.05,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,350000), labels = c("0", "350,000")) +
  labs(x = element_blank(), y = "intersection size\n(cell counts)")

counts_of_subjects <- melt(epithelial_numeric, id.vars = "combination",
                           measure.vars = markers,
                           variable.name = "subject", value.name = "present")[,
                                                                              .(counts = sum(present)), by = subject][order(-counts)]

points_data <- counts_combinations[, .(combination,
                                       subjects = unlist(strsplit(as.character(combination), ","))),
                                   by = combination]
points_data <- points_data[, .(combination, subjects)]
points_data[, subjects := factor(subjects, levels = rev(markers))]

subject_labels <- setNames(paste0(levels(points_data$subjects), "$^+$"), levels(points_data$subjects))
point_chart <- ggplot(points_data, aes(x = combination, y = subjects)) +
  geom_line(aes(group = combination), col = 'black') +
  geom_point(size = 3, col = 'black') +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    axis.ticks.x = element_blank(),
    panel.grid.major = element_line(colour = "darkgray", linetype = "dotted"),
    panel.grid.minor = element_line(colour = "darkgray", linetype = "dotted"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = unit(c(0,0,0,0), units = "cm"),
    axis.line.x = element_blank(),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA)
  ) +
  labs(x = element_blank(), y = element_blank()) +
  scale_x_discrete(drop = TRUE, expand = c(0.05,0.02)) +
  scale_y_discrete(expand = c(0.1,0.1), labels = sapply(subject_labels, TeX))

counts_of_subjects[, subject := factor(subject, levels = rev(markers))]
subject_bars <- ggplot(counts_of_subjects, aes(x = counts, y = subject)) +
  geom_col(width = 0.6, fill = 'darkgray') +
  theme_classic()+
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 15),
    axis.text.x = element_text(color = "black", size = 15),
    axis.title.x = element_text(color = "black", size = 15),
    plot.margin = unit(c(0,0,0,0), units = "cm"),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA)
  ) +
  coord_cartesian(expand = FALSE) +
  labs(x = "counts", y = NULL) +
  scale_y_discrete(labels = sapply(subject_labels, TeX)) +
  scale_x_continuous(breaks = c(0,500000), labels = c("0", "500,000"))

patient_counts <- epithelial_numeric[, .N, by = .(PatientID, combination)]
patient_counts <- patient_counts[N > 10, .N, by = combination][order(-N)]
patient_counts <- patient_counts[combination %in% counts_combinations$combination]
patient_counts[, combination := factor(combination, levels = counts_combinations$combination)]
patient_bar_chart <- ggplot(patient_counts, aes(x = combination, y = N)) +
  geom_col(width = 0.6, fill = 'darkgray') +
  theme_classic() +
  labs(x = "", title = TeX("$\\textit{n}$ patients (with > 10 cells)"), y = "") +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_text(color = "black", size = 15, vjust = -1),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 15),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5),
    plot.margin = unit(c(0,0,0,0), units = "cm"),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA)
  ) +
  coord_cartesian() +
  scale_x_discrete(expand = c(0.05,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 400))

ki67_percent <- epithelial_numeric[, .(percent_positive = mean(isKi67Pos) * 100), by = combination]
ki67_percent <- ki67_percent[combination %in% counts_combinations$combination]
ki67_percent[, combination := factor(combination, levels = counts_combinations$combination)]
ki67_bar_chart <- ggplot(ki67_percent, aes(x = combination, y = percent_positive)) +
  geom_col(width = 0.6, fill = 'darkgray') +
  theme_classic() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    plot.title = element_text(color = "black", size = 14, hjust = 0.5),
    axis.title.y = element_text(color = "black", size = 16, vjust = -1),
    axis.text.y = element_text(color = "black", size = 15),
    plot.margin = unit(c(0,0,0,0), units = "cm"),
    panel.background = element_rect(fill = NA, colour = NA),
    plot.background = element_rect(fill = NA, colour = NA)
  ) +
  coord_cartesian() +
  scale_x_discrete(expand = c(0.05,0)) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 5), labels = c("0", "5")) +
  labs(x = element_blank(), title = TeX("%Ki67$^+$"), y = "")

upset_plot <- plot_spacer() + bar_chart +
  subject_bars +
  point_chart + plot_spacer() + ki67_bar_chart + plot_spacer() + patient_bar_chart +
  plot_layout(ncol = 2, widths = c(0.5, 1), heights = c(2,1.5,0.4,0.4))

pdf(file.path(outdir, "UpSetPlot.pdf"),
    width = 6, height = 7, bg = "transparent")
print(upset_plot)
dev.off()
