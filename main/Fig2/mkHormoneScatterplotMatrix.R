library(here)
source(here::here("code/header.R"))

outdir <- here("scratch/main/Fig2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
panel <- getPanel()

cells[, (panel_cols) := lapply(.SD, function(x) as.numeric(scale(clip99(x)))), .SDcols = panel_cols]

epithelial <- cells[isEpithelial == TRUE]

markers <- c("FOXA1", "GATA3", "ER", "PR", "AR", "HER2")
colors <- viridisLite::viridis(256)
threshold <- 1
max_counts <- 20000
colors[1:threshold] <- "white"
  
plot_list <- vector("list", length(markers)^2)

for (i in 1:length(markers)) {
  for (j in 1:length(markers)) {
    
    plot_index <- (i - 1) * length(markers) + j
    
    if (j < i) {
      
      # make a tiny table for this panel (keeps memory low)
      dt <- epithelial[, .(x = get(markers[j]), y = get(markers[i]))]
      
      spearman_cor <- suppressWarnings(cor.test(dt$y, dt$x, method = "spearman"))
      cor_value <- spearman_cor$estimate
      
      x_min <- min(dt$x, na.rm = TRUE)
      x_max <- max(dt$x, na.rm = TRUE)
      y_min <- min(dt$y, na.rm = TRUE)
      y_max <- max(dt$y, na.rm = TRUE)
      
      sampled_data <- dt[sample.int(nrow(dt), min(1000, nrow(dt))), ]
      loess_fit <- loess(y ~ x, data = sampled_data)
      dt[, loess_smooth := predict(loess_fit, newdata = data.frame(x = x))]
      
      p <- ggplot(dt, aes(x = x, y = y)) +
        geom_hex() +
        geom_line(aes(y = loess_smooth), color = "darkgray", linetype = "dashed", size = 1.5) +
        scale_fill_gradientn(
          colours = colors,
          limits = c(1, max_counts),
          na.value = "white",
          oob = scales::oob_squish
        ) +
        theme_classic() +
        labs(x = markers[j], y = markers[i]) +
        scale_y_continuous(expand = c(0, 0), limits = c(NA, y_max - 0.2)) +
        scale_x_continuous(expand = c(0, 0), limits = c(NA, x_max - 0.2)) +
        guides(fill = "none") +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()
        )
      
      if (i == 6) {
        p <- p + theme(axis.title.x = element_text(size = 20, color = "black",
                                                   margin = margin(t = 5, r = 0, b = 0, l = 0)))
      }
      if (j == 1) {
        p <- p + theme(axis.title.y = element_text(size = 20, color = "black", hjust = 1,
                                                   margin = margin(t = 0, r = 5, b = 0, l = 0)))
      }
      if (i == 6 & j == 1) {
        p <- p + theme(
          axis.title.x = element_text(size = 20, color = "black",
                                      margin = margin(t = 5, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(size = 20, color = "black", hjust = 1,
                                      margin = margin(t = 0, r = 5, b = 0, l = 0))
        )
      }
      
      if (i == 5 | i == 6) {
        p <- p + annotate("text",
                          x = x_max - 0.4,
                          y = ((y_max - y_min) * 0.07) + y_min,
                          label = TeX(paste0("$\\rho$ = ", round(cor_value, 2))),
                          color = "black", hjust = 1, size = 7.5)
      } else if (i == 2 | i == 3 | i == 4) {
        p <- p + annotate("text",
                          x = x_max - 0.4,
                          y = ((y_max - y_min) * 0.07) + y_min,
                          label = TeX(paste0("$\\rho$ = ", round(cor_value, 2))),
                          color = "white", hjust = 1, size = 7.5)
      }
      
      plot_list[[plot_index]] <- p
      
    } else {
      plot_list[[plot_index]] <- ggplot() + theme_void()
    }
  }
}

plot_grid <- wrap_plots(plot_list, ncol = 6)
blanklabelplotx <- ggplot() + labs(x = "rescaled intensities") + theme_classic() +
  guides(x = "none", y = "none") + theme(axis.title.x = element_text(size = 25))
plot_grid <- plot_grid / blanklabelplotx +
  plot_layout(heights = c(1000, 1)) +
  plot_annotation(theme = theme(plot.margin = margin(l = 30, b = 30)))

ggsave(file.path(outdir, "hormoneScatterplotMatrix.pdf"),
       plot_grid, width = 14, height = 14)