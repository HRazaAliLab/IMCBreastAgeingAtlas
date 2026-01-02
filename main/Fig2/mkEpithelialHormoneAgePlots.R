library(here)
source(here("code", "header.R"))

outdir <- here("scratch/main/Fig2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical, by = "ImageID", all.x = TRUE)

epithelialclin <- cellsclin[isEpithelial == TRUE & !is.na(PatientID) & !is.na(Age)]

markers_is <- c("isFOXA1Pos", "isGATA3Pos", "isERPos", "isPRPos", "isARPos", "isHER2Pos")

# Per-patient proportion positive
proportions_list <- lapply(markers_is, function(marker) {
  proportions <- epithelialclin[, .(Proportion = mean(get(marker), na.rm = TRUE)), by = PatientID]
  unique_ages  <- unique(epithelialclin[, .(PatientID, Age)])
  merge(proportions, unique_ages, by = "PatientID")
})

print(nrow(proportions_list[[1]])) # number of patients used

format_p_value <- function(p) {
  p_rounded <- signif(p, digits = 1)
  if (p_rounded < 0.001) formatC(p_rounded, format = "e", digits = 1) else formatC(p_rounded, format = "f", digits = 2)
}

plots <- list()
ymaxes <- c(0.35, 0.3, 0.3, 0.3, 0.3, 0.1)

for (i in seq_along(markers_is)) {
  marker_name <- gsub("is(.*)Pos", "\\1", markers_is[i])
  marker_name <- paste0(marker_name, "$^+$")
  
  correlation_test <- suppressWarnings(
    cor.test(proportions_list[[i]]$Age, proportions_list[[i]]$Proportion, method = "spearman")
  )
  
  p <- ggplot(proportions_list[[i]], aes(x = Age, y = Proportion)) +
    geom_point(size = 0.1, alpha = 0.2, color = "black") +
    geom_smooth(method = "loess", color = "dimgray", se = FALSE) +
    labs(subtitle = TeX(marker_name)) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, ymaxes[[i]]),
      breaks = c(0, ymaxes[[i]]),
      expand = c(0, 0),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 70), expand = c(0, 0), limits = c(NA, 70)) +
    annotate(
      "text", x = 17, y = ymaxes[[i]] * 0.95,
      label = TeX(paste0("$\\rho$=", format(correlation_test$estimate, digits = 2))),
      size = 4.5, hjust = 0, color = "black"
    ) +
    annotate(
      "text", x = 17, y = ymaxes[[i]] * 0.80,
      label = TeX(paste0("p-value=", mkEnumPower(format_p_value(correlation_test$p.value)))),
      size = 4.5, hjust = 0, color = "black"
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.subtitle = element_text(
        hjust = 0.5, color = "black", size = 14,
        margin = margin(t = -20, unit = "pt"), vjust = 1
      ),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      plot.margin = margin(t = 25, r = 10, b = 0, l = 5, unit = "pt")
    )
  
  plots[[i]] <- p
}

plot_grid <- wrap_plots(plots, ncol = 3)

blanklabelploty <- ggplot() +
  labs(y = "epithelial proportion") +
  theme_classic() +
  guides(x = "none", y = "none") +
  theme(
    axis.title.y = element_text(color = "black", size = 14),
    plot.margin = unit(c(0, 0, 0, 0), units = "pt")
  )

blanklabelplotx <- ggplot() +
  labs(x = "age") +
  theme_classic() +
  guides(x = "none", y = "none") +
  theme(
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    axis.title.x = element_text(color = "black", size = 14)
  )

plot_grid <- blanklabelploty + plot_grid + plot_layout(widths = c(1, 1000))
plot_grid <- plot_grid / blanklabelplotx + plot_layout(heights = c(1000, 1))

pdf(file.path(outdir, "epithelialHormoneAgePlots.pdf"), width = 6.5, height = 4.5)
print(plot_grid)
dev.off()