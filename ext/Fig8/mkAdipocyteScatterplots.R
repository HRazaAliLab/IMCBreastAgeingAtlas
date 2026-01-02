library(here)
source(here("code", "header.R"))

outdir <- here("scratch/ext/Fig8")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

inDir <- here("data/derived")
compartmentSummaryPath <- file.path(inDir, "compartmentSummary.parquet")
dt <- read_parquet_adt(compartmentSummaryPath)

make_three_adipose <- function(dt) {
  
  # 1) proportion vs age (sqrt)
  dt1 <- na.omit(copy(dt), cols = c("Age","adiposeAreaFraction"))
  cor1 <- cor.test(dt1$Age, dt1$adiposeAreaFraction, method = "spearman")
  maxy1 <- 0.2
  
  p <- ggplot(dt1, aes(x = Age, y = sqrt(adiposeAreaFraction))) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE, color = "#B8860B") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, maxy1), expand = c(0,0), breaks = c(0, maxy1),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 50, 80), limits = c(15, 80), expand = c(0,0)) +
    annotate("text", x = 17, y = maxy1 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(cor1$estimate, digits = 2))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy1 * 0.8,
             label = TeX(paste0("pval=", mkEnumPower(format_custom_pval(cor1$p.value)))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy1 * 0.1,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(dt1))),
             size = 5, hjust = 0, color = "black") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      plot.margin = unit(c(10,10,10,10), "pt")
    ) +
    xlab("age") +
    ylab("fat area / total tissue area\nper patient (sqrt)")
  
  # 2) mean area vs age
  dt2 <- na.omit(copy(dt), cols = c("Age","adiposeMeanArea"))
  cor2 <- cor.test(dt2$Age, dt2$adiposeMeanArea, method = "spearman")
  maxy2 <- 2500
  
  q <- ggplot(dt2, aes(x = Age, y = adiposeMeanArea)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE, color = "#B8860B") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, maxy2), expand = c(0,0), breaks = c(0, maxy2),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 50, 80), limits = c(15, 80), expand = c(0,0)) +
    annotate("text", x = 17, y = maxy2 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(cor2$estimate, digits = 2))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy2 * 0.8,
             label = TeX(paste0("pval=", mkEnumPower(format_custom_pval(cor2$p.value)))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy2 * 0.1,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(dt2))),
             size = 5, hjust = 0, color = "black") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      plot.margin = unit(c(10,10,10,10), "pt")
    ) +
    xlab("age") +
    ylab(TeX("average adipocyte area ($\\mu m^2$)"))
  
  # 3) epithelial distance vs age
  dt3 <- na.omit(copy(dt), cols = c("Age","adiposeMeanDistanceToEpithelium"))
  cor3 <- cor.test(dt3$Age, dt3$adiposeMeanDistanceToEpithelium, method = "spearman")
  maxy3 <- 900
  
  r <- ggplot(dt3, aes(x = Age, y = adiposeMeanDistanceToEpithelium)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE, color = "#B8860B") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, maxy3), expand = c(0,0), breaks = c(0, maxy3),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(breaks = c(20, 50, 80), limits = c(15, 80), expand = c(0,0)) +
    annotate("text", x = 17, y = maxy3 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(cor3$estimate, digits = 2))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy3 * 0.8,
             label = TeX(paste0("pval=", mkEnumPower(format_custom_pval(cor3$p.value)))),
             size = 5, hjust = 0, color = "black") +
    annotate("text", x = 17, y = maxy3 * 0.1,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(dt3))),
             size = 5, hjust = 0, color = "black") +
    theme(
      axis.title.x = element_blank(), axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 10),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(color = "black", size = 14),
      plot.margin = unit(c(10,10,10,10), "pt")
    ) +
    xlab("age") +
    ylab("mean nearest adipocyte distance\nfrom epithelial cells (μm)")
  
  p + q + r + plot_layout(nrow = 1)
}

combined <- make_three_adipose(dt)

pdf(file.path(outdir, "adipocyteScatterplots.pdf"), width = 6, height = 2)
print(combined)
dev.off()

cat("Wrote: ", file.path(outdir, "adipocyteScatterplots.pdf"), "\n", sep = "")