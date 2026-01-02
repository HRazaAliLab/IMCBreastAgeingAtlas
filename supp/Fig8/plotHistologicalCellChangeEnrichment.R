library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "supp", "Fig8")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
ctx   <- getCellContext()
ann   <- getCellClusters()

epi_flags <- c("isSecretoryChange", "isApocrineMetaplasia", "isFibrocysticChange",
               "isColumnarChange", "isHyperplasia")
tme_flag  <- "isNearNormalEpithelium"

epi_flags <- epi_flags[epi_flags %in% names(ctx)]
has_tme   <- tme_flag %in% names(ctx)

keep_cols <- c("ImageID", "CellID", epi_flags)
if (has_tme) keep_cols <- c(keep_cols, tme_flag)

dt <- merge(
  cells[, .(ImageID, CellID, ClusterID, isEpithelial)],
  ctx[, ..keep_cols],
  by = c("ImageID", "CellID"),
  all.x = TRUE
)

to01 <- function(x) {
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x != 0))
  if (is.character(x)) return(as.integer(tolower(x) %in% c("true","t","1","yes","y")))
  as.integer(NA)
}

epi <- dt[isEpithelial == TRUE]
epi <- epi[complete.cases(epi[, ..epi_flags])]

if (nrow(epi) > 0 && length(epi_flags) > 0) {
  epi[, (epi_flags) := lapply(.SD, to01), .SDcols = epi_flags]
  
  epi[, NormalHistology := as.integer(rowSums(.SD, na.rm = TRUE) == 0), .SDcols = epi_flags]
  
  epi_flagset <- c(epi_flags, "NormalHistology")
  
  observed_epi <- epi[, lapply(.SD, sum, na.rm = TRUE), by = ClusterID, .SDcols = epi_flagset]
  
  total_cells_epi <- nrow(epi)
  prevalence_epi  <- epi[, lapply(.SD, sum, na.rm = TRUE), .SDcols = epi_flagset] / total_cells_epi
  cluster_totals  <- epi[, .(total_cells_in_ClusterID = .N), by = ClusterID]
  
  expected_epi <- copy(cluster_totals)
  for (nm in epi_flagset) expected_epi[, (nm) := total_cells_in_ClusterID * as.numeric(prevalence_epi[[nm]])]
  
  enrich_epi <- merge(observed_epi, expected_epi, by = "ClusterID", suffixes = c("_obs", "_exp"))
  for (nm in epi_flagset) {
    enrich_epi[, (paste0(nm, "_enrichment")) :=
                 get(paste0(nm, "_obs")) / pmax(get(paste0(nm, "_exp")), 1e-12)]
  }
  
  enrich_epi <- enrich_epi[, c("ClusterID", paste0(epi_flagset, "_enrichment")), with = FALSE]
  enrich_epi[, Source := "Epithelial"]
} else {
  enrich_epi <- data.table()
}

tme <- dt[isEpithelial == FALSE]

if (has_tme) {
  tme_flagset <- c(epi_flags, tme_flag)
  tme_flagset <- tme_flagset[tme_flagset %in% names(tme)]
  if (length(tme_flagset) == 0) stop("No TME histology proximity flags found.")

  tme <- tme[complete.cases(tme[, ..tme_flagset])]
  tme[, (tme_flagset) := lapply(.SD, to01), .SDcols = tme_flagset]

  tme[, NormalHistology := get(tme_flag)]

  tme_enrich_flags <- c(epi_flags, "NormalHistology")
  tme_enrich_flags <- tme_enrich_flags[tme_enrich_flags %in% names(tme)]
  
  observed_tme <- tme[, lapply(.SD, sum, na.rm = TRUE), by = ClusterID, .SDcols = tme_enrich_flags]

  total_cells_tme <- nrow(tme)
  prevalence_tme  <- tme[, lapply(.SD, sum, na.rm = TRUE), .SDcols = tme_enrich_flags] / total_cells_tme
  
  cluster_totals_tme <- tme[, .(total_cells_in_ClusterID = .N), by = ClusterID]
  expected_tme <- copy(cluster_totals_tme)
  for (nm in tme_enrich_flags) expected_tme[, (nm) := total_cells_in_ClusterID * as.numeric(prevalence_tme[[nm]])]
  
  enrich_tme <- merge(observed_tme, expected_tme, by = "ClusterID", suffixes = c("_obs","_exp"))
  for (nm in tme_enrich_flags) {
    enrich_tme[, (paste0(nm, "_enrichment")) := get(paste0(nm, "_obs")) / pmax(get(paste0(nm, "_exp")), 1e-12)]
  }
  
  enrich_tme <- enrich_tme[, c("ClusterID", paste0(tme_enrich_flags, "_enrichment")), with = FALSE]
  enrich_tme[, Source := "Microenvironment"]
  
} else {
  enrich_tme <- data.table()
}
sum(dt[isEpithelial == FALSE]$isApocrineMetaplasia)

combined <- rbindlist(list(enrich_epi, enrich_tme), fill = TRUE)
if (nrow(combined) == 0) stop("No rows to plot after filtering / missing calls.")

annotations <- copy(ann)

annotations[, ClusterID := as.factor(ClusterID)]
combined[,   ClusterID := as.factor(ClusterID)]
combined[,   ClusterID := factor(ClusterID, levels = annotations$ClusterID)]

combined <- merge(
  combined,
  annotations[, .(ClusterID, Type)],
  by = "ClusterID",
  all.x = TRUE
)

labels          <- setNames(annotations$BackupTeXClusterLabel, annotations$ClusterID)
phenotypecolors <- setNames(annotations$Colour,              annotations$ClusterID)
phenotypelabels <- setNames(annotations$BackupTeXClusterLabel, annotations$ClusterID)

long <- melt(
  combined,
  id.vars      = c("ClusterID", "Type"),
  measure.vars = patterns("_enrichment$"),
  variable.name = "EnrichmentType",
  value.name    = "EnrichmentValue"
)
long[, CellChange := gsub("_enrichment$", "", EnrichmentType)]

change_map <- c(
  isHyperplasia        = "Hyperplasia",
  isColumnarChange     = "ColumnarChange",
  isFibrocysticChange  = "FibrocysticChange",
  isApocrineMetaplasia = "ApocrineMetaplasia",
  isSecretoryChange    = "SecretoryChange",
  NormalHistology      = "NormalHistology"
)
long[, CellChange := change_map[as.character(CellChange)]]

ordered_levels <- c("Hyperplasia","ColumnarChange","FibrocysticChange","ApocrineMetaplasia","SecretoryChange","NormalHistology")
ordered_levels <- ordered_levels[ordered_levels %in% unique(long$CellChange)]
long[, CellChange := factor(CellChange, levels = ordered_levels)]

enrichmentplot <- ggplot(long, aes(x = ClusterID, y = factor(CellChange), fill = EnrichmentValue)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white", high = "black", na.value = "white",
    limits = c(1, max(long$EnrichmentValue, na.rm = TRUE)),
    breaks = c(1, 15),
    labels = c(TeX("$\\leq$1"), "15")
  ) +
  theme_classic() +
  scale_x_discrete(labels = sapply(labels, TeX)) +
  scale_y_discrete(labels = c(
    "Hyperplasia"        = "hyperplasia",
    "ColumnarChange"     = "columnar cell change",
    "FibrocysticChange"  = "fibrocystic change",
    "ApocrineMetaplasia" = "apocrine metaplasia",
    "SecretoryChange"    = "secretory cell change",
    "NormalHistology"    = "normal"
  )) +
  labs(x = "phenotype", y = NULL, fill = "enrichment") +
  guides(color = "none", size = guide_legend(nrow = 3, title.position = "top", direction = "vertical")) +
  theme(
    axis.text.x        = element_blank(),
    axis.text.y        = element_text(color = "black", size = 12),
    axis.ticks.length.x = unit(0.6, units = "pt"),
    axis.title.x       = element_blank(),
    axis.ticks.x       = element_blank(),
    strip.text.x       = element_blank(),
    legend.position    = "right",
    legend.title       = element_text(size = 12),
    plot.margin        = unit(c(0,0,0,0), units = "cm"),
    legend.text        = element_text(color = "black", size = 10),
    legend.ticks       = element_blank()
  ) +
  facet_grid(cols = vars(Type), scales = "free_x", space = "free")

annotations[, ClusterID := factor(ClusterID, levels = levels(long$ClusterID))]
annotations[, ClusterID := factor(ClusterID, levels = unique(as.character(ClusterID))), by = Type]

phenotypeclust_colorbar_plot <- ggplot(annotations, aes(y = factor(1), x = ClusterID, fill = ClusterID)) +
  geom_tile(color = "white", size = 0.3) +
  scale_x_discrete(labels = sapply(phenotypelabels, TeX)) +
  scale_y_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = phenotypecolors) +
  theme_minimal() +
  theme(
    axis.title         = element_blank(),
    axis.text.y        = element_blank(),
    axis.text.x        = element_text(angle = 90, vjust = 0.5, hjust = 1, color = "black", size = 11),
    axis.ticks.length.x = unit(0, "pt"),
    axis.ticks         = element_blank(),
    panel.grid         = element_blank(),
    panel.border       = element_blank(),
    strip.text.x       = element_blank(),
    legend.position    = "none",
    plot.margin        = unit(c(0,0,0,0), units = "cm")
  ) +
  facet_grid(~ Type, scales = "free_x", space = "free")

final <- enrichmentplot / phenotypeclust_colorbar_plot + plot_layout(heights = c(5, 0.2))

ggsave(
  filename = file.path(outdir, "histologicalCellChangeEnrichment.pdf"),
  plot = final, units = "in", width = 9, height = 3
)