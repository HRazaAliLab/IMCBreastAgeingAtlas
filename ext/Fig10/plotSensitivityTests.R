# Requires output from prepareSensitivityCS.py, prepareSensitivityKL.py, and prepareSensitivityCLF.py for plotting.

library(here)
source(here("code", "header.R"))

base_dir <- here("scratch/ext/Fig10")
outdir <- base_dir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

kl_master  <- fread(file.path(base_dir, "kl_master.csv"))
kl_roc     <- fread(file.path(base_dir, "kl_roc.csv"))
cs_master  <- fread(file.path(base_dir, "cs_master.csv"))
cs_roc     <- fread(file.path(base_dir, "cs_roc.csv"))
clf_master <- fread(file.path(base_dir, "clf_master.csv"))
clf_roc    <- fread(file.path(base_dir, "clf_roc.csv"))

kl_master[, variable := as.numeric(variable)]
cs_master[, variable := as.numeric(variable)]
clf_master[, variable := as.numeric(variable)]
kl_roc[, V1 := as.numeric(V1)]
cs_roc[, V1 := as.numeric(V1)]
clf_roc[, V1 := as.numeric(V1)]

kl_master_epi <- kl_master[which == "epi"]
kl_master_imm <- kl_master[which == "imm"]
kl_master_str <- kl_master[which == "str"]

createScatterplot <- function(data, title, linecolor) {
  plot <- ggplot(data, aes(x = variable, y = value)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE, color = linecolor) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(
      limits = c(0, 1), expand = c(0, 0),
      breaks = c(0, 1),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(
      expand = c(0.02, 0),
      limits = c(400, 1200),
      breaks = c(400, 1200)
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      plot.margin = unit(c(10, 15, 10, 15), "pt")
    ) +
    labs(subtitle = title)
  return(plot)
}

createConnectedScatterplot <- function(data, ymax, y_column, linecolor) {
  plot <- ggplot(data, aes_string(x = "V1", y = y_column)) +
    geom_line(color = linecolor, size = 0.8) +
    geom_vline(xintercept = 600, linetype = "dashed", color = "dimgray", size = 0.5) +
    geom_point(size = 1, alpha = 0.8) +
    theme_classic() +
    coord_cartesian(clip = "off") +
    scale_y_continuous(
      expand = c(0, 0),
      limits = c(0, NA),
      breaks = c(0, ymax),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    scale_x_continuous(
      expand = c(0.02, 0),
      limits = c(400, 1200),
      breaks = c(400, 1200)
    ) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 12),
      axis.text.y = element_text(color = "black", size = 12),
      plot.margin = unit(c(0, 15, 0, 15), "pt")
    )
  
  return(plot)
}

plot_epi <- createScatterplot(kl_master_epi, title = "Epithelial cells", linecolor = "#6C7ACE")
plot_imm <- createScatterplot(kl_master_imm, title = "Immune cells",      linecolor = "#7CA876")
plot_str <- createScatterplot(kl_master_str, title = "Stromal cells",     linecolor = "#EFAA9C")

plot_epi_roc <- createConnectedScatterplot(kl_roc, ymax = 0.01, y_column = "epi", linecolor = "#6C7ACE")
plot_imm_roc <- createConnectedScatterplot(kl_roc, ymax = 0.03, y_column = "imm", linecolor = "#7CA876")
plot_str_roc <- createConnectedScatterplot(kl_roc, ymax = 0.02, y_column = "str", linecolor = "#EFAA9C")

plot_cs <- ggplot(cs_master, aes(x = variable, y = value)) +
  geom_point(size = 0.1, alpha = 0.2) +
  geom_smooth(method = "loess", se = FALSE, color = "dimgray") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.5, 1), expand = c(0, 0),
    breaks = c(0.5, 1),
    labels = function(x) ifelse(x == 1, "1", format(x))
  ) +
  scale_x_continuous(
    expand = c(0.02, 0),
    limits = c(400, 1200),
    breaks = c(400, 1200)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    plot.margin = unit(c(10, 15, 10, 15), "pt")
  )

plot_cs_roc <- ggplot(cs_roc, aes(x = V1, y = cosine_sim)) +
  geom_line(color = "dimgray", size = 0.8) +
  geom_vline(xintercept = 600, linetype = "dashed", color = "dimgray", size = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(-0.04, -0.01),
    breaks = c(-0.04, -0.01),
    labels = function(x) ifelse(x == 0, "0", format(x))
  ) +
  scale_x_continuous(
    expand = c(0.02, 0),
    limits = c(400, 1200),
    breaks = c(400, 1200)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    plot.margin = unit(c(0, 15, 0, 15), "pt")
  )

clf_master_label         <- clf_master[which == "label"]
clf_master_spatial_label <- clf_master[which == "spatial_label"]

setorder(clf_master_label, variable)
setorder(clf_master_spatial_label, variable)

plot_clf_label <- ggplot(clf_master_label, aes(x = variable, y = value)) +
  geom_line(color = "dimgray", size = 0.8) +
  geom_point(size = 1, alpha = 0.8) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0, 1), expand = c(0, 0),
    breaks = c(0, 1),
    labels = function(x) ifelse(x == 1, "1", format(x))
  ) +
  scale_x_continuous(
    expand = c(0.02, 0),
    limits = c(0, 1400),
    breaks = c(0, 1400)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    plot.margin = unit(c(10, 15, 10, 15), "pt")
  )

plot_clf_spatial_label <- ggplot(clf_master_spatial_label, aes(x = variable, y = value)) +
  geom_line(color = "dimgray", size = 0.8) +
  geom_point(size = 1, alpha = 0.8) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.6, 0.67), expand = c(0, 0),
    breaks = c(0.6, 0.67),
    labels = function(x) format(x)
  ) +
  scale_x_continuous(
    expand = c(0.02, 0),
    limits = c(0, 1400),
    breaks = c(0, 1400)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    plot.margin = unit(c(10, 15, 10, 15), "pt")
  )

plot_clf_label_roc <- ggplot(clf_roc, aes(x = V1, y = label)) +
  geom_line(color = "dimgray", size = 0.8) +
  geom_vline(xintercept = 120, linetype = "dashed", color = "dimgray", size = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(-0.3, NA),
    breaks = c(-0.3, 0),
    labels = function(x) ifelse(x == 0, "0", format(x))
  ) +
  scale_x_continuous(
    expand = c(0.02, 0),
    limits = c(0, 1500),
    breaks = c(0, 1500)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    plot.margin = unit(c(0, 15, 0, 15), "pt")
  )

plot_clf_spatial_label_roc <- ggplot(clf_roc, aes(x = V1, y = spatial_label)) +
  geom_line(color = "dimgray", size = 0.8) +
  geom_vline(xintercept = 120, linetype = "dashed", color = "dimgray", size = 0.5) +
  geom_point(size = 1, alpha = 0.8) +
  theme_classic() +
  coord_cartesian(clip = "off") +
  scale_y_continuous(
    expand = c(0, 0),
    breaks = c(-0.03, 0),
    labels = function(x) ifelse(x == 0, "0", format(x))
  ) +
  scale_x_continuous(
    expand = c(0.02, 0),
    limits = c(0, 1500),
    breaks = c(0, 1500)
  ) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.text.y = element_text(color = "black", size = 12),
    plot.margin = unit(c(0, 15, 0, 15), "pt")
  )

combinedplot <- plot_epi + plot_imm + plot_str + plot_cs + plot_clf_label + plot_clf_spatial_label +
  plot_epi_roc + plot_imm_roc + plot_str_roc + plot_cs_roc + plot_clf_label_roc + plot_clf_spatial_label_roc +
  plot_layout(ncol = 6)

ggsave(
  combinedplot,
  width = 14, height = 3, units = "in",
  filename = file.path(outdir, "sensitivityTests.pdf")
)
