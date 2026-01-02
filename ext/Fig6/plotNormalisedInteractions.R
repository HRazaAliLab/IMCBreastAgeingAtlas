library(here)
source(here("code", "header.R"))

#Requires the following additional dependencies
library(SpatialExperiment)
library(imcRtools)
library(BiocParallel)

outdir <- here("scratch/ext/Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- getCells()
clinical <- getClinical()

if (!("X" %in% names(cells))) setnames(cells, old = "CenterX", new = "X", skip_absent = TRUE)
if (!("Y" %in% names(cells))) setnames(cells, old = "CenterY", new = "Y", skip_absent = TRUE)

cellsclin <- merge(
  cells[, .(ImageID, CellID, X, Y, isEpithelial)],
  clinical[, .(ImageID, PatientID, Age)],
  by = "ImageID",
  all.x = TRUE
)

cellsclin[, CellNb := seq_len(.N)]

# -----------------------------
#build SpatialExperiment + Delaunay graph
# -----------------------------
spe <- SpatialExperiment(
  assays = list(counts = matrix(0, nrow = 0, ncol = nrow(cellsclin))),  # empty assay placeholder
  spatialCoords = as.matrix(cellsclin[, .(X, Y)])
)
colData(spe) <- S4Vectors::DataFrame(cellsclin[, .(ImageID, CellID, CellNb, isEpithelial, PatientID, Age)])
colnames(spe) <- as.character(cellsclin$CellNb)

n_workers <- max(1, BiocParallel::bpparam()$workers - 1)
bp <- BiocParallel::MulticoreParam(workers = n_workers, RNGseed = 13, progressbar = TRUE)

spe <- buildSpatialGraph(
  spe,
  img_id = "ImageID",
  type = "delaunay",
  max_dist = 20,
  coords = c("X", "Y"),
  BPPARAM = bp
)

to_list   <- spe@int_colData@listData$colPairs@listData$delaunay_interaction_graph@hits@to
from_list <- spe@int_colData@listData$colPairs@listData$delaunay_interaction_graph@hits@from

edges <- data.table(source = from_list, target = to_list)

cell_meta <- as.data.table(colData(spe))
cell_meta[, CellNb := as.integer(as.character(CellNb))]

edges <- merge(edges, cell_meta[, .(CellNb, ImageID, isEpithelial)], by.x = "source", by.y = "CellNb", all.x = TRUE)
setnames(edges, "isEpithelial", "IsEpithelial_source")
edges <- merge(edges, cell_meta[, .(CellNb, isEpithelial)], by.x = "target", by.y = "CellNb", all.x = TRUE)
setnames(edges, "isEpithelial", "IsEpithelial_target")

edges <- edges[!is.na(ImageID)]

simplify_edges <- function(dt) {
  # undirected: sort endpoints, then unique
  dt2 <- copy(dt)
  dt2[, a := pmin(source, target)]
  dt2[, b := pmax(source, target)]
  unique(dt2[, .(a, b)])
}

epi_epi <- edges[IsEpithelial_source == TRUE & IsEpithelial_target == TRUE, .(source, target, ImageID)]
epi_epi_s <- epi_epi[, simplify_edges(.SD), by = ImageID]
epi_epi_counts <- epi_epi_s[, .(interactions = .N), by = ImageID]

epi_counts <- cellsclin[isEpithelial == TRUE, .(counts = .N), by = .(ImageID, PatientID, Age)]
epi_epi_final <- merge(epi_counts, epi_epi_counts, by = "ImageID", all.x = TRUE)
epi_epi_final[is.na(interactions), interactions := 0]
epi_epi_agg <- epi_epi_final[, .(
  aggregated_counts = sum(counts),
  aggregated_interactions = sum(interactions)
), by = .(PatientID, Age)]
epi_epi_agg[, normalizedinteractions := aggregated_interactions / aggregated_counts]

ct1 <- cor.test(epi_epi_agg$Age, epi_epi_agg$normalizedinteractions, method = "spearman")

miny1 <- 2
maxy1 <- 3
p1 <- ggplot(epi_epi_agg, aes(x = Age, y = normalizedinteractions)) +
  geom_point(size = 0.2, alpha = 0.2, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "dimgray") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(miny1, maxy1), expand = c(0, 0), breaks = c(miny1, maxy1)) +
  scale_x_continuous(breaks = c(20, 80), limits = c(15, 80), expand = c(0, 0)) +
  annotate("text", x = 16, y = ((maxy1 - miny1) * 0.9) + miny1,
           label = TeX(paste0("$\\rho$ = ", format(ct1$estimate, digits = 1))),
           size = 5, hjust = 0, color = "black") +
  annotate("text", x = 16, y = ((maxy1 - miny1) * 0.7) + miny1,
           label = TeX(paste0("p-value = ", format(ct1$p.value, digits = 1))),
           size = 5, hjust = 0, color = "black") +
  annotate("text", x = 16, y = ((maxy1 - miny1) * 0.1) + miny1,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(epi_epi_agg))),
           size = 5, hjust = 0, color = "black") +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    plot.margin = unit(c(5, 5, 5, 5), "pt")
  )

# Epi–TME interactions (heterotypic)
epi_tme <- edges[
  (IsEpithelial_source == TRUE & IsEpithelial_target == FALSE) |
    (IsEpithelial_source == FALSE & IsEpithelial_target == TRUE),
  .(source, target, ImageID)
]
epi_tme_s <- epi_tme[, simplify_edges(.SD), by = ImageID]
epi_tme_counts <- epi_tme_s[, .(interactions = .N), by = ImageID]

all_counts <- cellsclin[, .(counts = .N), by = .(ImageID, PatientID, Age)]
epi_tme_final <- merge(all_counts, epi_tme_counts, by = "ImageID", all.x = TRUE)
epi_tme_final[is.na(interactions), interactions := 0]
epi_tme_agg <- epi_tme_final[, .(
  aggregated_counts = sum(counts),
  aggregated_interactions = sum(interactions)
), by = .(PatientID, Age)]
epi_tme_agg[, normalizedinteractions := aggregated_interactions / aggregated_counts]

ct2 <- cor.test(epi_tme_agg$Age, epi_tme_agg$normalizedinteractions, method = "spearman")

miny2 <- 0
maxy2 <- 0.6
p2 <- ggplot(epi_tme_agg, aes(x = Age, y = normalizedinteractions)) +
  geom_point(size = 0.2, alpha = 0.2, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "dimgray") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(miny2, maxy2), expand = c(0, 0), breaks = c(miny2, maxy2)) +
  scale_x_continuous(breaks = c(20, 80), limits = c(15, 80), expand = c(0, 0)) +
  annotate("text", x = 16, y = ((maxy2 - miny2) * 0.9) + miny2,
           label = TeX(paste0("$\\rho$ = ", format(ct2$estimate, digits = 1))),
           size = 5, hjust = 0, color = "black") +
  annotate("text", x = 16, y = ((maxy2 - miny2) * 0.7) + miny2,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(ct2$p.value)))),
           size = 5, hjust = 0, color = "black") +
  annotate("text", x = 16, y = ((maxy2 - miny2) * 0.1) + miny2,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(epi_tme_agg))),
           size = 5, hjust = 0, color = "black") +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    plot.margin = unit(c(5, 5, 5, 5), "pt")
  )

# TME–TME interactions
tme_tme <- edges[IsEpithelial_source == FALSE & IsEpithelial_target == FALSE, .(source, target, ImageID)]
tme_tme_s <- tme_tme[, simplify_edges(.SD), by = ImageID]
tme_tme_counts <- tme_tme_s[, .(interactions = .N), by = ImageID]

tme_counts <- cellsclin[isEpithelial == FALSE, .(counts = .N), by = .(ImageID, PatientID, Age)]
tme_tme_final <- merge(tme_counts, tme_tme_counts, by = "ImageID", all.x = TRUE)
tme_tme_final[is.na(interactions), interactions := 0]
tme_tme_agg <- tme_tme_final[, .(
  aggregated_counts = sum(counts),
  aggregated_interactions = sum(interactions)
), by = .(PatientID, Age)]
tme_tme_agg[, normalizedinteractions := aggregated_interactions / aggregated_counts]

ct3 <- cor.test(tme_tme_agg$Age, tme_tme_agg$normalizedinteractions, method = "spearman")

miny3 <- 0.8
maxy3 <- 1.6
p3 <- ggplot(tme_tme_agg, aes(x = Age, y = normalizedinteractions)) +
  geom_point(size = 0.2, alpha = 0.2, color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "dimgray") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(miny3, maxy3), expand = c(0, 0), breaks = c(miny3, maxy3)) +
  scale_x_continuous(breaks = c(20, 80), limits = c(15, 80), expand = c(0, 0)) +
  annotate("text", x = 16, y = ((maxy3 - miny3) * 0.9) + miny3,
           label = TeX(paste0("$\\rho$ = ", format(ct3$estimate, digits = 1))),
           size = 5, hjust = 0, color = "black") +
  annotate("text", x = 16, y = ((maxy3 - miny3) * 0.7) + miny3,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(ct3$p.value)))),
           size = 5, hjust = 0, color = "black") +
  annotate("text", x = 16, y = ((maxy3 - miny3) * 0.1) + miny3,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(tme_tme_agg))),
           size = 5, hjust = 0, color = "black") +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(color = "black", size = 14),
    plot.margin = unit(c(5, 5, 5, 5), "pt")
  )

interactionplottotal <- p1 + p2 + p3 + plot_layout(ncol = 1)

out_fp <- file.path(outdir, "normalisedInteractions.pdf")
ggsave(interactionplottotal, filename = out_fp, width = 2.2, height = 4)

message("Wrote: ", out_fp)
