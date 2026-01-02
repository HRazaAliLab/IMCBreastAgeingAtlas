library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "main", "Fig7")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
ctx      <- getCellContext()
clinical <- getClinical()
ann      <- getCellClusters()

cells <- merge(
  cells,
  ctx[, .(ImageID, CellID, TissueStructure)],
  by = c("ImageID", "CellID"),
  all.x = TRUE
)
cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)

epi <- cellsclin[isEpithelial == TRUE]
epi <- epi[!is.na(TissueStructure)]
epi <- epi[TissueStructure %in% c("duct", "lobule")]

epiann <- ann[isEpithelial == TRUE]

prepareEffectSizeData_epi <- function(epi_dt) {
  #patients with both duct + lobule for comparison
  pats <- epi_dt[, .(
    has_duct   = any(TissueStructure == "duct"),
    has_lobule = any(TissueStructure == "lobule")
  ), by = PatientID][has_duct & has_lobule, PatientID]
  
  d <- epi_dt[PatientID %in% pats]
  
  complete_grid <- CJ(
    PatientID = unique(d$PatientID),
    ClusterID = unique(d$ClusterID),
    TissueStructure = c("duct", "lobule"),
    sorted = FALSE
  )
  
  counts <- d[, .(Count = .N), by = .(PatientID, ClusterID, TissueStructure)]
  counts <- merge(complete_grid, counts, by = c("PatientID", "ClusterID", "TissueStructure"), all.x = TRUE)
  counts[is.na(Count), Count := 0L]
  
  totals <- d[, .(TotalCount = .N), by = .(PatientID, TissueStructure)]
  counts <- merge(counts, totals, by = c("PatientID", "TissueStructure"), all.x = TRUE)
  counts[, Proportion := fifelse(TotalCount > 0, Count / TotalCount, 0)]
  
  counts
}

summarizeEffectSize_epi <- function(dt) {
  out <- data.table(
    Cluster = character(),
    MedianProportionDuct = numeric(),
    MedianProportionLobule = numeric(),
    CliffDelta = numeric(),
    LowerCI = numeric(),
    UpperCI = numeric(),
    pValue = numeric()
  )
  
  for (cl in unique(dt$ClusterID)) {
    x <- dt[ClusterID == cl]
    duct   <- x[TissueStructure == "duct",   .(PatientID, Proportion)]
    lobule <- x[TissueStructure == "lobule", .(PatientID, Proportion)]
    paired <- merge(duct, lobule, by = "PatientID", suffixes = c("_duct", "_lobule"))
    if (nrow(paired) == 0) next
    
    wt <- wilcox.test(paired$Proportion_lobule, paired$Proportion_duct, paired = TRUE)
    
    cd <- cliff.delta(paired$Proportion_duct, paired$Proportion_lobule)
    
    out <- rbind(out, data.table(
      Cluster = as.character(cl),
      MedianProportionDuct   = median(paired$Proportion_duct),
      MedianProportionLobule = median(paired$Proportion_lobule),
      CliffDelta = -as.numeric(cd$estimate),
      LowerCI = -cd$conf.int[2],
      UpperCI = -cd$conf.int[1],
      pValue = wt$p.value
    ))
  }
  
  out[, adjpValue := p.adjust(pValue, method = "BH")]
  out[, alpha := ifelse(adjpValue < 0.05, "opaque", "translucent")]
  out[order(-abs(CliffDelta))]
}

plotEffectSize <- function(summary_dt, labels_map) {
  summary_dt <- as.data.table(summary_dt)
  summary_dt[, Cluster := as.character(Cluster)]
  
  ord <- summary_dt[order(-CliffDelta), Cluster]
  summary_dt[, ClusterOrd := factor(Cluster, levels = ord)]
  
  summary_dt[, point_size := {
    ci_width <- UpperCI - LowerCI
    min_size <- 2
    max_size <- 5
    max_ci_width_allowed <- max(ci_width, na.rm = TRUE)
    if (is.finite(max_ci_width_allowed) && max_ci_width_allowed > 0) {
      sz <- min_size + (max_size - min_size) * (ci_width / max_ci_width_allowed)
      pmin(pmax(sz, min_size), max_size)
    } else {
      rep(max_size, .N)
    }
  }]
  
  ggplot(summary_dt, aes(x = ClusterOrd, y = CliffDelta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, alpha = alpha), width = 0.2) +
    geom_point(shape = 21, color = "black", aes(fill = alpha, size = point_size)) +
    scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
    scale_size_identity() +
    coord_flip(
      clip = "off",
      ylim = c(min(summary_dt$LowerCI, na.rm = TRUE), max(summary_dt$UpperCI, na.rm = TRUE))
    ) +
    labs(x = NULL, y = "Cliff's delta effect size") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.ticks.length.y = unit(0, "cm"),
      axis.text.x  = element_text(colour = "black", size = 15),
      axis.title.x = element_text(size = 15),
      plot.margin  = unit(c(40, 100, 5, 0), "pt"),
      legend.position = "none"
    ) +
    scale_x_discrete(labels = sapply(labels_map, TeX)) +
    guides(color = "none", alpha = "none", fill = "none", size = "none")
}

createColorBarPlot <- function(summary_dt, ann_dt) {
  summary_dt <- as.data.table(summary_dt)
  ann_dt     <- as.data.table(ann_dt)
  
  summary_dt[, Cluster := as.character(Cluster)]
  ann_dt[, ClusterID := as.character(ClusterID)]
  
  tmp <- merge(
    ann_dt[, .(ClusterID, BackupLabel, Colour)],
    summary_dt[, .(Cluster, CliffDelta)],
    by.x = "ClusterID", by.y = "Cluster",
    all.x = FALSE, all.y = TRUE
  )
  setorder(tmp, CliffDelta)
  
  tmp[, TeXLabel := mkTeX(BackupLabel, ignoreChar = "/")]
  tmp[, ClusterID := factor(ClusterID, levels = tmp$ClusterID)]
  
  clusternames <- setNames(tmp$TeXLabel, tmp$ClusterID)
  colors       <- setNames(tmp$Colour,   tmp$ClusterID)
  
  ggplot(tmp, aes(x = factor(1), y = ClusterID, fill = ClusterID)) +
    geom_tile(color = NA) +
    scale_y_discrete(
      limits = rev(levels(tmp$ClusterID)),
      labels = rev(sapply(clusternames, TeX))
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(
      axis.title  = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 15),
      axis.ticks  = element_blank(),
      panel.grid  = element_blank(),
      legend.position = "none",
      plot.margin = unit(c(0,0,0,0), "cm")
    )
}

epi_merged <- prepareEffectSizeData_epi(epi)
epi_summary <- summarizeEffectSize_epi(epi_merged)

epi_summary[, adjpValueLabel := if (exists("format_custom_pval") && exists("mkEnumPower")) {
  sapply(format_custom_pval(adjpValue), mkEnumPower)
} else sprintf("%.3g", adjpValue)]

epiann[, ClusterID := as.character(ClusterID)]
labels_map <- setNames(epiann$BackupTeXClusterLabel, epiann$ClusterID)

p_eff <- plotEffectSize(epi_summary, labels_map) +
  scale_y_continuous(
    breaks = c(-0.7, 0, 0.5),
    labels = function(x) ifelse(x == 0, "0", as.character(x))
  ) +
  geom_text(
    aes(
      y = 0.65, x = ClusterOrd,
      label = TeX(adjpValueLabel, output = "character"),
      color = ifelse(alpha == "translucent", "darkgray", "black")
    ),
    parse = TRUE, hjust = 0, size = 4.5
  ) +
  scale_color_identity()

p_cb  <- createColorBarPlot(epi_summary, epiann)

final <- p_cb + p_eff + plot_layout(widths = c(0.3, 7))

pdf(file.path(outdir, "ductLobuleEpithelialComposition.pdf"), width = 4.5, height = 4)
print(final)
dev.off()
