library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "supp", "Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
ctx      <- getCellContext()
clinical <- getClinical()
ann      <- getCellClusters()

cells <- merge(cells, ctx[, .(ImageID, CellID, TissueStructure)], by=c("ImageID","CellID"), all.x=TRUE)
cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by="ImageID", all.x=TRUE)

epi <- cellsclin[isEpithelial==TRUE & TissueStructure %in% c("duct","lobule")]
epi[, ClusterID := as.character(ClusterID)]

epiann <- ann[isEpithelial==TRUE]
epiann[, ClusterID := as.character(ClusterID)]

prepareEffectSizeData_epi <- function(epi_dt) {
  epi_dt <- as.data.table(copy(epi_dt))
  
  pats <- epi_dt[, .(
    has_duct   = any(TissueStructure=="duct"),
    has_lobule = any(TissueStructure=="lobule")
  ), by=PatientID][has_duct & has_lobule, PatientID]
  
  d <- epi_dt[PatientID %in% pats]
  
  grid <- CJ(
    PatientID       = unique(d$PatientID),
    ClusterID       = unique(d$ClusterID),
    TissueStructure = c("duct","lobule"),
    sorted = FALSE
  )
  
  counts <- d[, .(Count=.N), by=.(PatientID,ClusterID,TissueStructure)]
  counts <- merge(grid, counts, by=c("PatientID","ClusterID","TissueStructure"), all.x=TRUE)
  counts[is.na(Count), Count:=0L]
  
  totals <- d[, .(TotalCount=.N), by=.(PatientID,TissueStructure)]
  counts <- merge(counts, totals, by=c("PatientID","TissueStructure"), all.x=TRUE)
  
  counts[, Proportion := fifelse(TotalCount>0, Count/TotalCount, 0)]
  counts
}

summarizePvals_epi <- function(dt) {
  dt <- as.data.table(dt)
  
  out <- data.table(Cluster=character(), pValue=numeric())
  
  for (cl in unique(dt$ClusterID)) {
    x <- dt[ClusterID==cl]
    
    duct <- x[TissueStructure=="duct",   .(PatientID, Proportion)]
    lob  <- x[TissueStructure=="lobule", .(PatientID, Proportion)]
    paired <- merge(duct, lob, by="PatientID", suffixes=c("_duct","_lobule"))
    
    if (nrow(paired) == 0) next
    paired <- paired[is.finite(Proportion_duct) & is.finite(Proportion_lobule)]
    if (nrow(paired) == 0) next
    
    wt <- wilcox.test(paired$Proportion_lobule, paired$Proportion_duct, paired=TRUE)
    out <- rbind(out, data.table(Cluster=as.character(cl), pValue=wt$p.value))
  }
  
  out[, adjpValue := p.adjust(pValue, method="BH")]
  out
}

epi_merged <- prepareEffectSizeData_epi(epi)
epi_pvals  <- summarizePvals_epi(epi_merged)

createPairedClusterPlot <- function(mergedData, pvals, annotations, cluster_id, ymax) {
  cluster_id <- as.character(cluster_id)
  
  d <- mergedData[ClusterID==cluster_id]
  label <- annotations[ClusterID==cluster_id, BackupTeXClusterLabel][1]
  adjP  <- pvals[Cluster==cluster_id, adjpValue][1]
  
  d[, TissueStructure := factor(TissueStructure, levels=c("duct","lobule"))]
  
  wide <- dcast(d, PatientID ~ TissueStructure, value.var="Proportion")
  wide <- na.omit(wide)
  setorder(wide, PatientID)
  
  ptxt <- if (is.na(adjP)) "p=NA" else paste0("p=", mkEnumPower(format_custom_pval(adjP)))
  
  ggplot(melt(wide, id.vars="PatientID"), aes(x=variable, y=value, fill=variable)) +
    geom_boxplot(outlier.size=0.1) +
    scale_fill_manual(values=c("duct"="#4169E1","lobule"="#DC143C")) +
    theme_classic() +
    scale_x_discrete(labels=c("duct","lobule")) +
    scale_y_continuous(
      limits=c(0,ymax),
      breaks=c(0,floor(ymax/0.05)*0.05),
      expand=expand_scale(mult=c(0,0.1)),
      labels=function(x) ifelse(x==0,"0",format(x))
    ) +
    labs(x=NULL, y=NULL, title=TeX(label)) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text=element_text(color="black", size=20),
      plot.title=element_text(size=20, hjust=0.5),
      plot.margin=unit(c(0,0,20,0),"pt"),
      legend.position="none"
    ) +
    annotate("text", x=1.5, y=ymax*0.95, label=TeX(ptxt), parse=TRUE, size=6.5, hjust=0.5)
}

p10 <- createPairedClusterPlot(epi_merged, epi_pvals, epiann, "10", ymax=0.66)
p11 <- createPairedClusterPlot(epi_merged, epi_pvals, epiann, "11", ymax=0.20)

combined <- p10 + p11

ylab <- ggplot() +
  labs(y="epithelial\nproportion") +
  theme_classic() +
  theme(
    plot.margin=margin(0,0,0,0, unit="cm"),
    axis.title.y=element_text(color="black", size=22)
  ) +
  guides(x="none", y="none")

combined <- ylab + combined + plot_layout(widths=c(1,1000))

ggsave(
  file.path(outdir, "ductLobuleBasalMyoepithelialBoxplots.pdf"),
  combined, width=5, height=2.3, units="in"
)
