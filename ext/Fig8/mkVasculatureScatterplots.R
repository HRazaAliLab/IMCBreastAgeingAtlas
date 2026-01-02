library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig8")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

inDir <- here("data/derived")
compartmentSummaryPath <- file.path(inDir, "compartmentSummary.parquet")
dt <- read_parquet_adt(compartmentSummaryPath)

make_three_vasculature <- function(dt) {
  
  all_results_dt <- copy(dt)
  
  # 1) proportion vs age (sqrt)
  all_results_dt <- na.omit(all_results_dt, cols = c("Age", "vascularAreaFraction"))
  correlation_test <- cor.test(all_results_dt$Age, all_results_dt$vascularAreaFraction, method = "spearman")
  
  maxy1 <- 0.15
  p <- ggplot(all_results_dt, aes(x = Age, y = sqrt(vascularAreaFraction))) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE, color = "#8B004B") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, maxy1), expand = c(0,0), breaks = c(0, maxy1),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 50, 80), limits = c(15, 80), expand = c(0,0)) +
    annotate("text", x = 17, y = maxy1 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(correlation_test$estimate, digits = 2))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy1 * 0.8,
             label = TeX(paste0("pval=", mkEnumPower(format_custom_pval(correlation_test$p.value)))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy1 * 0.1,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(all_results_dt))),
             size = 5, hjust = 0, color = "black") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      plot.margin = unit(c(10,10,10,10), "pt")
    ) +
    xlab("age") +
    ylab("vascular area / total tissue area\nper patient (sqrt)")
  
  # 2) mean area vs age
  all_results_dt <- na.omit(all_results_dt, cols = c("Age", "vascularMeanArea"))
  correlation_test2 <- cor.test(all_results_dt$Age, all_results_dt$vascularMeanArea, method = "spearman")
  
  maxy2 <- 1200
  q <- ggplot(all_results_dt, aes(x = Age, y = vascularMeanArea)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", color = "#8B004B", se = FALSE) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, maxy2), expand = c(0,0), breaks = c(0, maxy2),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 50, 80), limits = c(15, 80), expand = c(0,0)) +
    annotate("text", x = 17, y = maxy2 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(correlation_test2$estimate, digits = 2))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy2 * 0.8,
             label = TeX(paste0("pval=", mkEnumPower(format_custom_pval(correlation_test2$p.value)))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy2 * 0.1,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(all_results_dt))),
             size = 5, hjust = 0, color = "black") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      plot.margin = unit(c(10,10,10,10), "pt")
    ) +
    xlab("age") +
    ylab(TeX("average vascular area ($\\mu m^2$)"))
  
  # 3) epithelial distance vs age
  all_results_dt <- na.omit(all_results_dt, cols = c("Age", "vascularMeanDistanceToEpithelium"))
  correlation_test3 <- cor.test(all_results_dt$Age, all_results_dt$vascularMeanDistanceToEpithelium, method = "spearman")
  
  maxy3 <- 400
  r <- ggplot(all_results_dt, aes(x = Age, y = vascularMeanDistanceToEpithelium)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", color = "#8B004B", se = FALSE) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, maxy3), expand = c(0,0), breaks = c(0, maxy3),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 50, 80), limits = c(15, 80), expand = c(0,0)) +
    annotate("text", x = 17, y = maxy3 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(correlation_test3$estimate, digits = 2))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy3 * 0.8,
             label = TeX(paste0("pval=", mkEnumPower(format_custom_pval(correlation_test3$p.value)))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy3 * 0.1,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(all_results_dt))),
             size = 5, hjust = 0, color = "black") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      plot.margin = unit(c(10,10,10,10), "pt")
    ) +
    xlab("age") +
    ylab("mean nearest vascular distance\nfrom epithelial cells (μm)")
  
  p + q + r + plot_layout(nrow = 1)
}

combined_vasculature <- make_three_vasculature(dt)

pdf(file.path(outdir, "vasculatureScatterplots.pdf"), width = 6, height = 2)
print(combined_vasculature)
dev.off()

cat("Wrote: ", file.path(outdir, "vasculatureScatterplots.pdf"), "\n", sep = "")
