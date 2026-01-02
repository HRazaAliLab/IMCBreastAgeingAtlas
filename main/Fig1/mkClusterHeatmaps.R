library(here)
source(here::here("code/header.R"))

outdir <- here("scratch/main/Fig1")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

cells <- getCells()
annotations <- getCellClusters()

epiannotations <- annotations[isEpithelial == TRUE]
tmeannotations <- annotations[isEpithelial == FALSE]

panel <- getPanel()

select_markers <- function(x, pattern, invert = FALSE) {
  grep(pattern, x, value = TRUE, invert = invert)
}

epiMarkers <- select_markers(
  panel$panel,
  "CK5/14|CK8/18|CK7|panCK|HLA-ABC|GATA|ER|PR|AR|TCF|CD15|LDHA|SMA|Calp|FOXA1"
)

tmeMarkers <- select_markers(
  panel$panel,
  "CK|Ki67|DNA|H3|ER|PR|AR|FOXA1|pH2AX|GATA3|Cyclin D1|IDO|FOXP3",
  invert = TRUE
)

marker_to_cells_col <- function(x) {
  x <- gsub("CK5/14", "CK5_14", x, fixed = TRUE)
  x <- gsub("CK8/18", "CK8_18", x, fixed = TRUE)
  x <- gsub("PD-1",  "PD1",  x, fixed = TRUE)
  x <- gsub("PD-L1", "PDL1", x, fixed = TRUE)
  x <- gsub("/", "_", x, fixed = TRUE)
  x
}

mkTilePlot <- function(props, subtitle, limits, label_fmt = "%.1f", hi = NULL) {
  p <- ggplot(props, aes(x = factor(1), y = ClusterID, fill = Percent)) +
    geom_tile(color = "white", linewidth = 0.1) +
    geom_text(aes(label = sprintf(label_fmt, Percent)),
              vjust = 0.5, color = "white", size = 4.4) +
    scale_x_discrete(expand = c(0, 0)) +
    labs(subtitle = subtitle, x = NULL, y = NULL) +
    theme_classic() +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line  = element_blank(),
      legend.position = "none",
      plot.subtitle = element_text(color = "black", size = 13, hjust = 0.5, vjust = -2),
      plot.margin = unit(c(0, 0, 0, 0), units = "cm")
    )
  
  if (is.null(hi)) {
    p <- p + scale_fill_gradient(low = "darkgrey", high = "darkblue", limits = limits)
  } else {
    p <- p + scale_fill_gradient(low = "darkgrey", high = hi, limits = limits)
  }
  
  p
}

mkggNormalClusterHeatMap <- function(dat, clusterlabels, clusterColumns, ByVar = "ClusterID", epithelial = FALSE) {
  
  stopifnot(ByVar %in% names(dat))
  stopifnot(all(c("ClusterID", "BackupLabel", "Colour") %in% names(clusterlabels)))
  
  # Align marker names between panel and cells
  clusterColumns <- marker_to_cells_col(clusterColumns)
  clusterColumns <- intersect(clusterColumns, names(dat))
  if (length(clusterColumns) == 0) stop("No marker columns found in `dat` after mapping panel names to cells colnames.")
  
  # Median marker intensity by cluster
  medians <- dat[, lapply(.SD, median, na.rm = TRUE), by = ByVar, .SDcols = clusterColumns]
  
  rescale01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (is.infinite(rng[1]) || diff(rng) == 0) return(rep(0, length(x)))
    (x - rng[1]) / diff(rng)
  }
  
  # Variable order uses *cells* column names
  variableorder_tme <- c(
    "CD45","CD3","CD4","CD8","TCF1","GZMB","PD1","CD20","CD79a","CD56","CD11c",
    "HLA-ABC","HLA-DR","CD15","MPO","CD68","CD163","PDL1","PDGFRbeta","CD31",
    "Caveolin-1","Vimentin","PDPN","SMA","Calponin","LDHA","CA9"
  )
  
  variableorder_epi <- c(
    "panCK","FOXA1","GATA3","ER","PR","AR","HER2","LDHA","CK8_18","HLA-ABC",
    "CK7","CD15","TCF1","CK5_14","SMA","Calponin"
  )
  
  variableorder <- if (epithelial) variableorder_epi else variableorder_tme
  
  # Rescale medians across clusters (0-1 within each marker)
  medians[, (clusterColumns) := lapply(.SD, rescale01), .SDcols = clusterColumns]
  
  medians_melted <- data.table::melt(medians, id.vars = ByVar)
  data.table::setnames(medians_melted, ByVar, "ClusterID")
  medians_melted$variable <- factor(medians_melted$variable, levels = variableorder)
  medians_melted$ClusterID <- factor(medians_melted$ClusterID, levels = rev(clusterlabels$ClusterID))
  
  # Counts per cluster
  count_data <- dat[, .N, by = ByVar]
  count_data_melted <- data.table::melt(count_data, id.vars = ByVar)
  data.table::setnames(count_data_melted, ByVar, "ClusterID")
  count_data_melted$ClusterID <- factor(count_data_melted$ClusterID, levels = rev(clusterlabels$ClusterID))
  
  # Axis labels: keys must be cells colnames; values can be pretty/TeX
  variablelabels <- if (epithelial) {
    setNames(
      c("panCK","FOXA1","GATA3","ER","PR","AR","HER2","LDHA","CK8/18","HLA-ABC",
        "CK7","CD15","TCF1","CK5/14","SMA","Calponin"),
      variableorder
    )
  } else {
    setNames(
      c("CD45","CD3","CD4","CD8","TCF1","GZMB","PD-1","CD20","CD79a","CD56","CD11c",
        "HLA-ABC","HLA-DR","CD15","MPO","CD68","CD163","PD-L1","$PDGFR\\beta$","CD31",
        "Caveolin-1","Vimentin","PDPN","SMA","Calponin","LDHA","CA9"),
      variableorder
    )
  }
  
  # Heatmap (marker medians)
  heatmap_plot <- ggplot(medians_melted, aes(x = variable, y = ClusterID, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "white", high = "black", mid = "white",
                         limits = c(0, 1), breaks = c(0, 1), name = NULL) +
    scale_x_discrete(expand = c(0, 0), labels = sapply(variablelabels, TeX)) +
    labs(x = NULL, y = NULL, title = ifelse(epithelial, "epithelial", "microenvironment")) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, color = "black", size = 14, hjust = 1, vjust = 0.5,
                                 margin = margin(t = -3, unit = "pt")),
      axis.line.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks  = element_blank(),
      axis.line.x = element_blank(),
      plot.title  = element_text(color = "black", hjust = 0.5, vjust = -4, size = 18),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.text = element_text(size = 16, color = "black"),
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  
  # Cluster label colorbar
  clusterlabels <- data.table::copy(clusterlabels)
  clusterlabels[, normalizedID := .GRP, by = ClusterID]
  clusterlabels[, BackupTeXClusterLabel := mkTeX(BackupLabel, ignoreChar = "/")]
  
  clusternames <- setNames(clusterlabels$BackupTeXClusterLabel, clusterlabels$normalizedID)
  
  colorbar_plot <- ggplot(clusterlabels, aes(x = factor(1), y = rev(normalizedID), fill = I(Colour))) +
    geom_tile(color = NA) +
    scale_y_discrete(limits = rev(as.factor(clusterlabels$normalizedID)),
                     labels = rev(sapply(clusternames, TeX))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_identity() +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 14),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0, 0, 0, 0), "cm")
    )
  
  # Counts bar
  if (epithelial) {
    countbreaks <- c(0, 200000)
    countlabels <- c(0, "0.2M")
    count_margin <- 0.6
  } else {
    countbreaks <- c(0, 500000)
    countlabels <- c(0, "0.5M")
    count_margin <- 0.2
  }
  
  count_plot <- ggplot(count_data_melted, aes(x = value, y = ClusterID)) +
    geom_bar(stat = "identity", fill = "dimgray") +
    scale_x_continuous(breaks = countbreaks, labels = countlabels) +
    theme_classic() +
    theme(
      axis.text.y = element_blank(),
      axis.title  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.line.x  = element_blank(),
      axis.text.x  = element_text(color = "black", size = 15),
      plot.margin  = unit(c(0, count_margin, 0, 0), units = "cm")
    )
  
  # % positive tiles (booleans already exist in `cells`)
  totalCounts <- dat[!is.na(get(ByVar)), .N, by = ByVar]
  
  calculateProportions <- function(conditionCol) {
    stopifnot(conditionCol %in% names(dat))
    pos <- dat[!is.na(get(ByVar)) & get(conditionCol), .N, by = ByVar]
    props <- merge(pos, totalCounts, by = ByVar, all.y = TRUE, suffixes = c(".pos", ".total"))
    props$N.pos[is.na(props$N.pos)] <- 0
    props[, Percent := (N.pos / N.total) * 100]
    props[, ClusterID := factor(get(ByVar), levels = rev(clusterlabels$ClusterID))]
    props
  }
  
  ki67_props <- calculateProportions("isKi67Pos")
  ki67_margin <- if (epithelial) 0.1 else 0.3
  
  ki67_plot <- mkTilePlot(
    props = ki67_props,
    subtitle = TeX(ifelse(epithelial, "Ki67", "%Ki67$^+$")),
    limits = c(0, 11),
    label_fmt = "%.1f"
  ) + theme(plot.margin = unit(c(0, ki67_margin, 0, 0.1), units = "cm"))
  
  # Layout (patchwork) — use wrap_plots to avoid ggplot_add pitfalls
  if (epithelial) {
    ER_plot <- mkTilePlot(calculateProportions("isERPos"),  "ER", limits = c(0, 100), label_fmt = "%.0f", hi = "darkred")
    PR_plot <- mkTilePlot(calculateProportions("isPRPos"),  "PR", limits = c(0, 100), label_fmt = "%.0f", hi = "darkred")
    AR_plot <- mkTilePlot(calculateProportions("isARPos"),  "AR", limits = c(0, 100), label_fmt = "%.0f", hi = "darkred")
    
    return(
      patchwork::wrap_plots(colorbar_plot, heatmap_plot, count_plot, ki67_plot, ER_plot, PR_plot, AR_plot) +
        patchwork::plot_layout(widths = c(0.5, 20, 3, 2.4, 2, 2, 2))
    )
  }
  
  patchwork::wrap_plots(colorbar_plot, heatmap_plot, count_plot, ki67_plot) +
    patchwork::plot_layout(widths = c(0.3, 20, 3, 1.3))
}

epiheatmap <- mkggNormalClusterHeatMap(
  dat = cells[isEpithelial == TRUE],
  clusterlabels = epiannotations,
  clusterColumns = epiMarkers,
  ByVar = "ClusterID",
  epithelial = TRUE
)

ggsave(
  filename = file.path(outdir, "epithelialHeatmap.pdf"),
  plot = epiheatmap,
  units = "in",
  width = 6.6,
  height = 5
)

tmeheatmap <- mkggNormalClusterHeatMap(
  dat = cells[isEpithelial == FALSE],
  clusterlabels = tmeannotations,
  clusterColumns = tmeMarkers,
  ByVar = "ClusterID",
  epithelial = FALSE
)

ggsave(
  filename = file.path(outdir, "microenvironmentHeatmap.pdf"),
  plot = tmeheatmap,
  units = "in",
  width = 8.1,
  height = 5.5
)