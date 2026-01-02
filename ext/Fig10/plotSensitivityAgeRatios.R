#Requires output from prepareSensitivityAgeRatios.py.

library(here)
source(here("code/header.R"))

outdir <- here("scratch", "ext", "Fig10")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ratios_fp <- file.path(outdir, "ratios.csv")
stats_fp  <- file.path(outdir, "panel_stats.csv")
series_fp <- file.path(outdir, "panel_series.csv")
roi_fp    <- file.path(outdir, "roi_order.csv")

ratios <- fread(ratios_fp)
stats  <- fread(stats_fp)
series <- fread(series_fp)
roi_order <- fread(roi_fp)

roi_order[, ImageID := as.character(ImageID)]
ratios[, ImageID := as.character(ImageID)]

stats[,  `:=`(Size = as.integer(Size),
              NImages = as.integer(NImages),
              p_value = as.numeric(p_value))]

series[, `:=`(Size = as.integer(Size),
              NImages = as.integer(NImages),
              Age = as.numeric(Age))]

setkey(roi_order, ImageID)
ratios[, ImageID := factor(ImageID, levels = roi_order$ImageID, ordered = TRUE)]

sizes <- c(400, 600, 800, 1000, 1200)
image_nums <- c(50, 100, 400, 700, 1400)

fmt_p <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  if (length(p) == 0 || is.na(p)) return("= NA")
  if (p < 0.005) return("< 0.005")
  sprintf("= %.2f", p)
}

plot_panel <- function(size_val, nimgs_val) {
  rois_n <- head(roi_order$ImageID, min(nimgs_val, nrow(roi_order)))
  
  pts <- ratios[as.character(ImageID) %in% rois_n & Size == size_val]
  
  ln <- series[Size == size_val & NImages == nimgs_val]
  setorder(ln, Age)
  
  pv <- stats[Size == size_val & NImages == nimgs_val, p_value]
  pv_str <- fmt_p(pv[1])
  
  max_x <- if (nrow(ln) > 0) max(ln$Age) * 0.99 else 80 * 0.99
  
  ggplot() +
    geom_point(data = pts, aes(x = Age, y = Ratio),
               size = 0.2, alpha = 0.2, color = "black") +
    geom_line(data = ln, aes(x = Age, y = linear_y),
              color = "darkgray", linewidth = 0.8) +
    geom_line(data = ln, aes(x = Age, y = spline_y),
              color = "dimgray", linewidth = 0.8) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(color = "black", size = 12),
      plot.margin = unit(c(10, 10, 10, 10), units = "pt")
    ) +
    scale_x_continuous(expand = c(0, 0), breaks = c(20, 70), limits = c(20, 80)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(0, 1), limits = c(0, 1)) +
    coord_cartesian(clip = "off") +
    annotate("text", x = 21, y = 0.1,
             label = paste0("p-value ", pv_str),
             hjust = 0, size = 4.5, color = "black")
}

plot_list <- list()
for (s in sizes) {
  for (n in image_nums) {
    plot_list[[paste0("size_", s, "_imgs_", n)]] <- plot_panel(s, n)
  }
}

combined_plot <- wrap_plots(plot_list, ncol = length(image_nums))

ggsave(
  filename = file.path(outdir, "sensitivityAgeRatios.pdf"),
  plot = combined_plot,
  width = 14, height = 9, units = "in"
)
