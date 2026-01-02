library(here)
source(here::here("code", "header.R"))

outDir <- here("scratch/ext/Fig4")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()
annotations <- getCellClusters()

cellContext <- read_parquet_adt(here("data", "derived", "cellContext.parquet"))
cellContext <- data.table::setDT(cellContext)[, .(ImageID, CellID, TissueStructure)]

cells <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cells <- merge(cells, cellContext, by = c("ImageID", "CellID"), all.x = TRUE)

cells[, ClusterID := as.integer(as.character(ClusterID))]
annotations[, ClusterID := as.integer(ClusterID)]
cells <- merge(cells, annotations[, .(ClusterID, TeXClusterLabel)], by = "ClusterID", all.x = TRUE)

cells <- cells[!is.na(TissueStructure)]

cells[, totalCells := .N, by = ImageID]
cells[, nCellsDuctal := .N, by = .(TissueStructure, ImageID)]
ductProp <- cells[TissueStructure == "duct"]
ductProp[, proportionDuctal := nCellsDuctal / totalCells]
ductProp <- ductProp[, .SD[1], by = ImageID][, .(ImageID, proportionDuctal)]

cellCounts <- getCellCounts(cells, clusterColumn = "ClusterID", idCols = "ImageID")

cellCounts <- merge(cellCounts, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
cellCounts[, overFifty := Age >= 50]
cellCounts <- merge(cellCounts, annotations[, .(ClusterID, TeXClusterLabel, isEpithelial)], by = "ClusterID", all.x = TRUE)
cellCounts <- merge(cellCounts, ductProp, by = "ImageID", all = FALSE)
cellCounts <- cellCounts[!is.na(Age)]
cellCounts[, adjustBy := sqrt(proportionDuctal)]

# ----------------------------
# Model runner
# ----------------------------
doClusterComparison <- function(countDat,
                                imgID = "ImageID",
                                ptID = "PatientID",
                                clusterCol = "TeXClusterLabel",
                                compareCol,
                                countColumn = "nCellsPerType",
                                totalColumn = "totalCells",
                                addPredictors = NULL) {
  
  stopifnot(data.table::is.data.table(countDat))
  d <- data.table::copy(countDat)
  
  stopifnot(compareCol %in% names(d))
  stopifnot(all(c(imgID, ptID, clusterCol) %in% names(d)))
  stopifnot(all(c(countColumn, totalColumn) %in% names(d)))
  
  d <- d[!is.na(get(clusterCol))][!is.na(get(compareCol))]
  
  predictors <- if (is.null(addPredictors)) compareCol else paste(compareCol, addPredictors, sep = " + ")
  predictors <- paste0(predictors, " + ")
  
  doComparisonBylme <- function(cluster, d, imgID, ptID) {
    d[, indicator := get(clusterCol) == cluster]
    fmla <- paste0(countColumn, "/", totalColumn, " ~ ", predictors, "(1|", ptID, ") + (1|", imgID, ")")
    
    lmeOut <- lme4::glmer(
      as.formula(fmla),
      family = binomial,
      weights = get(totalColumn),
      data = d[(indicator)]
    )
    
    mklinfctMatrix <- function(fitted) {
      m <- rep(0, nrow(coef(summary(fitted))))
      m[2] <- 1
      matrix(m, 1)
    }
    
    ht <- multcomp::glht(lmeOut, linfct = mklinfctMatrix(lmeOut))
    result <- summary(ht)
    
    effect <- log2(exp(coef(result)))
    pV <- result$test$pvalues[1]
    
    ci <- confint(ht)$confint
    ci <- c(ci[, "lwr"], ci[, "upr"])
    ci <- sapply(ci, function(x) log2(exp(x)))
    
    out <- data.table::data.table(
      ClusterID = cluster,
      nImages = d[(indicator), data.table::uniqueN(get(imgID))],
      nPatients = d[(indicator), data.table::uniqueN(get(ptID))],
      log2OR = effect,
      lci = ci[1],
      uci = ci[2],
      p = pV
    )
    
    d[, indicator := NULL]
    out
  }
  
  clusters_u <- d[, unique(get(clusterCol))]
  res <- data.table::rbindlist(lapply(clusters_u, doComparisonBylme, d = d, imgID = imgID, ptID = ptID))
  res[, adjP := p.adjust(p, method = "BH")]
  res
}

unadjusted <- doClusterComparison(cellCounts, clusterCol = "TeXClusterLabel", compareCol = "overFifty")
adjusted   <- doClusterComparison(cellCounts, clusterCol = "TeXClusterLabel", compareCol = "overFifty", addPredictors = "adjustBy")

estimateCols <- grep("ClusterID", names(adjusted), invert = TRUE, value = TRUE)
data.table::setnames(adjusted, estimateCols, paste0("Adj", estimateCols))

result <- merge(unadjusted, adjusted, by = "ClusterID")
data.table::setnames(result, "ClusterID", "TeXClusterLabel")

result <- merge(
  result,
  annotations[, .(TeXClusterLabel, isEpithelial)],
  by = "TeXClusterLabel",
  all.x = TRUE
)

ann <- data.table::copy(annotations)
ann[, isEpithelial := as.logical(isEpithelial)]
ann[, TeXClusterLabel := as.character(TeXClusterLabel)]

if (!"normalizedID" %in% names(ann)) {
  ann[, normalizedID := .I]
}
ann[, normalizedID := as.character(normalizedID)]

epiannotations <- ann[isEpithelial == TRUE, .(TeXClusterLabel, normalizedID, Colour, BackupLabel)]
tmeannotations <- ann[isEpithelial == FALSE, .(TeXClusterLabel, normalizedID, Colour, BackupLabel)]

prepareCellPhenotypesAdjusted <- function(data) {
  data <- data.table::as.data.table(data)
  data[, alpha := ifelse(adjP < 0.05, "opaque", "translucent")]
  data[, Adjalpha := ifelse(AdjadjP < 0.05, "opaqueadjusted", "translucentadjusted")]
  data <- data[order(-(log2OR))]
  return(data)
}

plotLog2OR <- function(summaryData, clusterLabels, offset = 0.15) {
  summaryData <- data.table::as.data.table(summaryData)
  summaryData$point_size <- 3
  summaryData$TeXClusterLabel <- as.numeric(factor(summaryData$TeXClusterLabel, levels = unique(summaryData$TeXClusterLabel)))
  labels <- setNames(as.character(clusterLabels$TeXClusterLabel), as.character(clusterLabels$normalizedID))
  p <- ggplot(summaryData) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_segment(aes(x = TeXClusterLabel - offset, xend = TeXClusterLabel - offset, y = lci, yend = uci, alpha = alpha), color = "black") +
    geom_segment(aes(x = TeXClusterLabel + offset, xend = TeXClusterLabel + offset, y = Adjlci, yend = Adjuci, alpha = Adjalpha), color = "black") +
    geom_point(aes(x = TeXClusterLabel - offset, y = log2OR, fill = alpha, size = point_size, color = alpha), shape = 21) +
    geom_point(aes(x = TeXClusterLabel + offset, y = Adjlog2OR, fill = Adjalpha, size = point_size, color = Adjalpha), shape = 21) +
    scale_color_manual(values = c("opaque" = "black", "translucent" = "gray", "opaqueadjusted" = "black", "translucentadjusted" = "gray")) +
    scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.3, "opaqueadjusted" = 1, "translucentadjusted" = 0.3)) +
    scale_fill_manual(values = c("opaque" = "#3D3D3D", "translucent" = "#3D3D3D", "opaqueadjusted" = "lightblue", "translucentadjusted" = "lightblue")) +
    coord_flip(clip = "off", ylim = c(min(c(summaryData$lci, summaryData$Adjlci)), max(c(summaryData$uci, summaryData$Adjuci))), xlim = c(1, nrow(summaryData))) +
    labs(x = NULL, y = TeX("log$_{2}$(OR)")) +
    theme_classic() +
    theme(axis.line.y = element_line(color = "black", size = 1),
          axis.title.y = element_text(colour = "black"),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.ticks.length.y = unit(0, "cm"),
          axis.text.x = element_text(colour = "black", size = 15),
          axis.title.x = element_text(size = 15),
          plot.margin = unit(c(40, 160, 5, 0), "pt"),
          legend.position = "bottom") +
    scale_size_identity() +
    scale_x_discrete(labels = sapply(labels, TeX)) +
    guides(color = "none", alpha = "none", fill = "none", size = "none")
  return(p)
}

createColorBarPlot <- function(summaryData, clusterLabels) {
  mergedData <- merge(clusterLabels, summaryData, by = "TeXClusterLabel")
  mergedData <- mergedData[order(log2OR)]
  mergedData[, TeXClusterLabel := mkTeX(BackupLabel, ignoreChar = '/')]
  clusternames <- setNames(mergedData$TeXClusterLabel, mergedData$normalizedID)
  colors <- setNames(mergedData$Colour, mergedData$normalizedID)
  mergedData[, normalizedID := factor(normalizedID, levels = names(colors))]
  colorbar_plot <- ggplot(mergedData, aes(x = factor(1), y = normalizedID, fill = normalizedID)) +
    geom_tile(color = "white", size = 0.3) +
    geom_vline(xintercept = 1.5, color = "black", size = 1) +
    scale_y_discrete(limits = rev(as.factor(mergedData$normalizedID)), labels = rev(sapply(clusternames, TeX))) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 15),
          axis.ticks.length.y = unit(0, "pt"),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          panel.border = element_blank(),
          legend.position = "none",
          plot.margin = unit(c(0,0,0,0), "cm"))
  return(colorbar_plot)
}

cellPhenotypesAdjusted <- data.table::as.data.table(result)

epiannotations$TeXClusterLabel <- as.character(epiannotations$TeXClusterLabel)
epiannotations$normalizedID <- as.character(epiannotations$normalizedID)

epicellPhenotypesAdjusted <- cellPhenotypesAdjusted[cellPhenotypesAdjusted$isEpithelial == TRUE, ]
epicellPhenotypesAdjusted <- prepareCellPhenotypesAdjusted(epicellPhenotypesAdjusted)
epicellPhenotypesAdjusted$adjpValueLabel <- mkEnumPower(format_custom_pval(epicellPhenotypesAdjusted$adjP))
epicellPhenotypesAdjusted$AdjadjpValueLabel <- mkEnumPower(format_custom_pval(epicellPhenotypesAdjusted$AdjadjP))

epicellPhenotypesEffectSizePlot <- plotLog2OR(epicellPhenotypesAdjusted, epiannotations)
epicellPhenotypesEffectSizePlot <- epicellPhenotypesEffectSizePlot +
  scale_y_continuous(breaks = c(-2, 0, 2), labels = function(x) ifelse(x == 0, "0", format(x))) +
  geom_text(aes(y = 2.4, x = TeXClusterLabel, label = TeX(adjpValueLabel, output = "character"), color = alpha),
            parse = TRUE, hjust = 0, size = 4.5) +
  geom_text(aes(y = 5.1, x = TeXClusterLabel, label = TeX(AdjadjpValueLabel, output = "character"), color = Adjalpha),
            parse = TRUE, hjust = 0, size = 4.5)

epicolorbar_plot <- createColorBarPlot(epicellPhenotypesAdjusted, epiannotations)
epicellPhenotypeplot <- epicolorbar_plot + epicellPhenotypesEffectSizePlot + plot_layout(widths = c(0.3, 7))

pdf(file.path(outDir, "EpiDuctLobuleAdjustedPhenotypes.pdf"), width = 5.5, height = 4.5)
print(epicellPhenotypeplot)
dev.off()

tmeannotations$TeXClusterLabel <- as.character(tmeannotations$TeXClusterLabel)
tmeannotations$normalizedID <- as.character(tmeannotations$normalizedID)

tmecellPhenotypesAdjusted <- cellPhenotypesAdjusted[cellPhenotypesAdjusted$isEpithelial == FALSE, ]
tmecellPhenotypesAdjusted <- prepareCellPhenotypesAdjusted(tmecellPhenotypesAdjusted)
tmecellPhenotypesAdjusted$adjpValueLabel <- mkEnumPower(format_custom_pval(tmecellPhenotypesAdjusted$adjP))
tmecellPhenotypesAdjusted$AdjadjpValueLabel <- mkEnumPower(format_custom_pval(tmecellPhenotypesAdjusted$AdjadjP))

tmecellPhenotypesEffectSizePlot <- plotLog2OR(tmecellPhenotypesAdjusted, tmeannotations)
tmecellPhenotypesEffectSizePlot <- tmecellPhenotypesEffectSizePlot +
  coord_flip(clip = "off",
             ylim = c(min(c(tmecellPhenotypesAdjusted$lci, tmecellPhenotypesAdjusted$Adjlci)), 1),
             xlim = c(1, nrow(tmecellPhenotypesAdjusted))) +
  scale_y_continuous(breaks = c(-2, 0, 1), labels = function(x) ifelse(x == 0, "0", format(x))) +
  geom_text(aes(y = 1.4, x = TeXClusterLabel, label = TeX(adjpValueLabel, output = "character"), color = alpha),
            parse = TRUE, hjust = 0, size = 4.5) +
  geom_text(aes(y = 3.6, x = TeXClusterLabel, label = TeX(AdjadjpValueLabel, output = "character"), color = Adjalpha),
            parse = TRUE, hjust = 0, size = 4.5)

tmecolorbar_plot <- createColorBarPlot(tmecellPhenotypesAdjusted, tmeannotations)
tmecellPhenotypeplot <- tmecolorbar_plot + tmecellPhenotypesEffectSizePlot + plot_layout(widths = c(0.3, 7))

pdf(file.path(outDir, "TMEDuctLobuleAdjustedPhenotypes.pdf"), width = 5.5, height = 4.5)
print(tmecellPhenotypeplot)
dev.off()

# ============================
# Proportion boxplots
# ============================

createProportionBoxplot <- function(clusterLabels, ymin, ymax, cluster_id, pvals) {
  cellcounts <- getCellCounts(cells, clusterColumn = "ClusterID", idCols = "PatientID")
  fibroblasts <- cellcounts[ClusterID == cluster_id]
  
  agemapping <- unique(cells[, .(PatientID, Age)])
  merged_data <- merge(fibroblasts, agemapping, by = "PatientID", all.x = TRUE)
  
  merged_data[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
  
  data.table::setnames(merged_data, "proportion", "Proportion")
  proportions <- merged_data
  proportions[, AgeGroup := factor(AgeGroup, levels = c("Below 50", "Above 50"))]
  proportions <- na.omit(proportions)
  
  print(nrow(proportions[AgeGroup == "Below 50"]))
  print(nrow(proportions[AgeGroup == "Above 50"]))
  
  # --- FIX: map ClusterID -> TeXClusterLabel using `annotations` (has ClusterID) ---
  label <- annotations[ClusterID == cluster_id, unique(TeXClusterLabel)]
  label <- as.character(label)
  
  # pull adjP from model results using the TeX label
  adjPVal <- pvals[TeXClusterLabel == label, adjP][1]
  adjPVal <- mkEnumPower(format_custom_pval(adjPVal))
  
  plot <- ggplot(proportions, aes(x = AgeGroup, y = Proportion, fill = AgeGroup)) +
    geom_boxplot(outlier.size = 0.1, outlier.color = "black") +
    scale_fill_manual(values = c("Below 50" = "gray", "Above 50" = "dimgray")) +
    theme_classic() +
    scale_x_discrete(labels = c("pre", "post")) +
    scale_y_continuous(
      limits = c(ymin, ymax),
      breaks = c(ymin, ymax),
      expand = expand_scale(mult = c(0, 0.1)),
      labels = c(ymin, ymax)
    ) +
    labs(x = NULL, y = NULL, title = TeX(label)) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 16),
      axis.title.y = element_text(vjust = -3, hjust = 0.5, size = 16),
      plot.title = element_text(hjust = 0.5, size = 16),
      plot.margin = unit(c(0, 5, 10, 0), units = "pt"),
      legend.position = "none"
    )
  
  plot <- plot + annotate(
    "text",
    x = 1.5,
    y = ymax,
    label = TeX(paste0("p=", adjPVal)),
    size = 5.5,
    color = "black",
    angle = 0,
    hjust = 0.5
  )
  
  return(plot)
}

plot1 <- createProportionBoxplot(tmeannotations, 0.05, 0.7, 25, cellPhenotypesAdjusted)
plot2 <- createProportionBoxplot(tmeannotations, 0.02, 0.15, 19, cellPhenotypesAdjusted)

combinedplots <- plot1 + plot2 + plot_layout(nrow = 1)

ggsave(
  combinedplots,
  filename = file.path(outDir, "boxplotsDuctLobuleAdjustedPhenotypes.pdf"),
  width = 3.5,
  height = 2.5
)

