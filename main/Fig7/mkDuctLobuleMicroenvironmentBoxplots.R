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
  ctx[, .(ImageID, CellID, TissueStructure, isNearDuct, isNearLobule)],
  by = c("ImageID","CellID"),
  all.x = TRUE
)
cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)

epi <- cellsclin[isEpithelial == TRUE & TissueStructure %in% c("duct","lobule")]
pats <- epi[, .(
  has_duct   = any(TissueStructure == "duct"),
  has_lobule = any(TissueStructure == "lobule")
), by = PatientID][has_duct & has_lobule, PatientID]

tme    <- cellsclin[isEpithelial == FALSE & PatientID %in% pats]
tmeann <- ann[isEpithelial == FALSE]

prepareEffectSizeDataMicroenvironment <- function(tme_dt) {
  tme_dt <- as.data.table(copy(tme_dt))

  both <- tme_dt[isNearDuct %in% TRUE & isNearLobule %in% TRUE]
  both2 <- copy(both); both2[, CellPhenotype := "NearLobule"]
  
  tme_dt[, CellPhenotype := NA_character_]
  tme_dt[isNearDuct %in% TRUE & !(isNearLobule %in% TRUE), CellPhenotype := "NearDuct"]
  tme_dt[isNearLobule %in% TRUE & !(isNearDuct %in% TRUE), CellPhenotype := "NearLobule"]
  tme_dt[!(isNearDuct %in% TRUE) & !(isNearLobule %in% TRUE), CellPhenotype := "ConnectiveTissue"]
  
  tme_dt[isNearDuct %in% TRUE & isNearLobule %in% TRUE, CellPhenotype := "NearDuct"]
  
  combined <- rbindlist(list(tme_dt, both2), use.names = TRUE, fill = TRUE)
  combined[is.na(CellPhenotype), CellPhenotype := "ConnectiveTissue"]
  
  grid <- CJ(
    PatientID     = unique(combined$PatientID),
    ClusterID     = unique(combined$ClusterID),
    CellPhenotype = c("NearDuct","NearLobule","ConnectiveTissue"),
    sorted = FALSE
  )
  
  counts <- combined[, .(Count = .N), by = .(PatientID, ClusterID, CellPhenotype)]
  counts <- merge(grid, counts, by = c("PatientID","ClusterID","CellPhenotype"), all.x = TRUE)
  counts[is.na(Count), Count := 0L]
  
  totals <- combined[, .(TotalCount = .N), by = .(PatientID, CellPhenotype)]
  counts <- merge(counts, totals, by = c("PatientID","CellPhenotype"), all.x = TRUE)
  
  counts[, Proportion := fifelse(TotalCount > 0, Count / TotalCount, 0)]
  counts
}

summarizePvalsTME <- function(dt) {
  dt <- as.data.table(dt)
  
  out <- data.table(Cluster = character(), pValue = numeric())
  
  for (cl in unique(dt$ClusterID)) {
    x <- dt[ClusterID == cl]
    
    d1 <- x[CellPhenotype == "NearDuct",   .(PatientID, Proportion)]
    d2 <- x[CellPhenotype == "NearLobule", .(PatientID, Proportion)]
    paired <- merge(d1, d2, by = "PatientID", suffixes = c("_duct","_lobule"))
    
    if (nrow(paired) == 0) next
    paired <- paired[is.finite(Proportion_duct) & is.finite(Proportion_lobule)]
    if (nrow(paired) == 0) next
    
    wt <- wilcox.test(paired$Proportion_duct, paired$Proportion_lobule, paired = TRUE)
    out <- rbind(out, data.table(Cluster = as.character(cl), pValue = wt$p.value))
  }
  
  out[, adjpValue := p.adjust(pValue, method = "BH")]
  out
}

tme_merged <- prepareEffectSizeDataMicroenvironment(tme)
tme_pvals  <- summarizePvalsTME(tme_merged)

createPairedTMEClusterPlot <- function(mergedData, pvals, clusterLabels, cluster_id, ymax) {
  cluster_id <- as.character(cluster_id)
  
  d <- mergedData[ClusterID == cluster_id & CellPhenotype %in% c("NearDuct","NearLobule")]
  label <- clusterLabels[ClusterID == cluster_id, BackupTeXClusterLabel][1]
  adjP  <- pvals[Cluster == cluster_id, adjpValue][1]
  
  wide <- dcast(d, PatientID ~ CellPhenotype, value.var = "Proportion")
  wide <- na.omit(wide)
  
  ptxt <- if (is.na(adjP)) "p=NA" else paste0("p=", mkEnumPower(format_custom_pval(adjP)))
  
  ggplot(melt(wide, id.vars="PatientID"), aes(x = variable, y = value, fill = variable)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("NearDuct" = "#4169E1", "NearLobule" = "#DC143C")) +
    theme_classic() +
    scale_x_discrete(labels = c("duct","lobule")) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = c(0, floor(ymax/0.05)*0.05),
      expand = expand_scale(mult = c(0, 0.1)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title = TeX(label)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black", size = 20),
      plot.title = element_text(size = 20, hjust = 0.5),
      plot.margin = unit(c(0,0,20,0), "pt"),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none"
    ) +
    annotate("text", x = 1.5, y = ymax*0.95, label = TeX(ptxt), parse = TRUE, size = 6.5, hjust = 0.5)
}

p1 <- createPairedTMEClusterPlot(tme_merged, tme_pvals, tmeann, "21", ymax=0.43)
p2 <- createPairedTMEClusterPlot(tme_merged, tme_pvals, tmeann, "16", ymax=0.20)
p3 <- createPairedTMEClusterPlot(tme_merged, tme_pvals, tmeann, "18", ymax=0.06)
p4 <- createPairedTMEClusterPlot(tme_merged, tme_pvals, tmeann, "24", ymax=0.21)
p5 <- createPairedTMEClusterPlot(tme_merged, tme_pvals, tmeann, "20", ymax=0.74)
p6 <- createPairedTMEClusterPlot(tme_merged, tme_pvals, tmeann, "14", ymax=0.185)

combined <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 3)

ylab <- ggplot() +
  labs(y = "microenvironment proportion") +
  theme_classic() +
  theme(
    plot.margin = margin(0,0,0,0, unit="cm"),
    axis.title.y = element_text(color="black", size=20)
  ) +
  guides(x = "none", y = "none")

combined <- ylab + combined + plot_layout(widths = c(1, 1000))

ggsave(
  file.path(outdir, "ductLobuleMicroenvironmentBoxplots.pdf"),
  combined, width = 7, height = 5, units = "in"
)

