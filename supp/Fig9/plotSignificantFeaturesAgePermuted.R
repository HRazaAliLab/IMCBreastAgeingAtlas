# NOTE: This script is memory- and CPU-intensive (large per-cell tables + image mask operations).
# Recommended: run on a machine/HPC node with ample RAM and multiple cores.

library(here)
source(here("code/header.R"))

# Dependencies:
#   DEswan (https://lehallib.github.io/DEswan/) — Lehallier et al., Nature Medicine (2019)
#   Used for age-dependent non-linear association testing.
library(DEswan)

library(igraph)
library(spatstat)
library(EBImage)

DERIVED_DIR <- here("data", "derived")
MASK_DIR    <- here("data/raw/images")
OUT_DIR     <- here("scratch", "supp", "Fig9")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FP_CELLS       <- file.path(DERIVED_DIR, "cells.parquet")
FP_CLINICAL    <- file.path(DERIVED_DIR, "clinical.parquet")
FP_CTX         <- file.path(DERIVED_DIR, "cellContext.parquet")
FP_ADJ         <- file.path(DERIVED_DIR, "cellAdjacency.parquet")
FP_COMP_SUM    <- file.path(DERIVED_DIR, "compartmentSummary.parquet")
FP_CLUSTER_ANN <- file.path(DERIVED_DIR, "cellClusterAnnotation.csv")
FP_SPATIAL_ANN <- file.path(DERIVED_DIR, "spatialClusterAnnotation.csv")

stop_if_missing <- function(fp) if (!file.exists(fp)) stop("Missing: ", fp)
for (fp in c(FP_CELLS, FP_CLINICAL, FP_CTX, FP_ADJ, FP_COMP_SUM, FP_CLUSTER_ANN, FP_SPATIAL_ANN)) stop_if_missing(fp)

cells    <- read_parquet_adt(FP_CELLS)
clinical <- read_parquet_adt(FP_CLINICAL)
ctx      <- read_parquet_adt(FP_CTX)
adj      <- read_parquet_adt(FP_ADJ)

comp_sum <- read_parquet_adt(FP_COMP_SUM)
clu_ann  <- fread(FP_CLUSTER_ANN)
sp_ann   <- fread(FP_SPATIAL_ANN)

panel <- getPanel()

tissue_by_patient <- clinical[, .(totalTissueArea = sum(TissueArea)), by = PatientID]
tissue_by_patient <- na.omit(tissue_by_patient)

clin_min <- unique(clinical[, .(ImageID, PatientID, Age)], by = "ImageID")

cluster_label_map <- setNames(clu_ann$BackupLabel, as.character(clu_ann$ClusterID))
spatial_label_map <- setNames(sp_ann$Label, as.character(sp_ann$SpatialClusterID))

label_for_cluster <- function(cid) {
  nm <- cluster_label_map[[as.character(cid)]]
  if (is.null(nm) || is.na(nm) || nm == "") nm <- paste0("Cluster ", cid)
  nm
}
label_for_spatial <- function(sid) {
  nm <- spatial_label_map[[as.character(sid)]]
  if (is.null(nm) || is.na(nm) || nm == "") nm <- paste0("Spatial ", sid)
  nm
}

# ----------------------------
# 1) Mean marker expression by patient
# ----------------------------
panel_markers <- intersect(panel$panel, names(cells))
if (length(panel_markers) == 0) stop("No panel markers found in cells.parquet that match panel$panel")

cells[, (panel_markers) := lapply(.SD, function(x) as.numeric(scale(clip99(x)))), .SDcols = panel_markers]

cellsclin <- merge(cells, clin_min, by = "ImageID", all.x = TRUE)
cellsclin <- cellsclin[!is.na(PatientID) & !is.na(Age)]

expr_by_patient <- cellsclin[
  , c(.(
    Age = Age[1]
  ), lapply(.SD, mean, na.rm = TRUE)),
  by = PatientID,
  .SDcols = panel_markers
]
expr_by_patient <- na.omit(expr_by_patient)

columns_to_remove <- c("Cyclin D1", "IDO", "FOXP3", "pH2AX", "DNA1", "DNA2", "CK19")
drop_present <- intersect(columns_to_remove, names(expr_by_patient))
if (length(drop_present) > 0) expr_by_patient[, (drop_present) := NULL]

expr_marker_cols <- setdiff(names(expr_by_patient), c("PatientID","Age"))
setnames(expr_by_patient, expr_marker_cols, paste("mean", expr_marker_cols, "intensity"))

# ----------------------------
# 2) Cluster proportions & densities (patient-level averages)
# ----------------------------
cellsclin_clean <- cellsclin[!is.na(PatientID) & !is.na(Age) & !is.na(ClusterID)]
cell_counts <- cellsclin_clean[, .N, by = .(PatientID, ClusterID, Age)]
cell_counts[, TotalCells := sum(N), by = PatientID]
cell_counts[, Proportion := fifelse(TotalCells > 0, N / TotalCells, 0)]

cell_counts <- merge(cell_counts, tissue_by_patient, by = "PatientID", all.x = TRUE)
cell_counts[, Density := N / totalTissueArea]

final_dt_avg <- cell_counts[
  , .(
    Average_Proportion = mean(Proportion, na.rm = TRUE),
    Average_Density    = mean(Density, na.rm = TRUE)
  ),
  by = .(PatientID, ClusterID, Age)
]

prop_wide <- dcast(final_dt_avg, PatientID + Age ~ ClusterID, value.var = "Average_Proportion", fill = 0)
dens_wide <- dcast(final_dt_avg, PatientID + Age ~ ClusterID, value.var = "Average_Density", fill = 0)

prop_cols <- setdiff(names(prop_wide), c("PatientID","Age"))
dens_cols <- setdiff(names(dens_wide), c("PatientID","Age"))

prop_new <- sapply(as.integer(prop_cols), function(cid) paste(label_for_cluster(cid), "proportion"))
dens_new <- sapply(as.integer(dens_cols), function(cid) paste(label_for_cluster(cid), "density"))

setnames(prop_wide, prop_cols, prop_new)
setnames(dens_wide, dens_cols, dens_new)

propdens_by_patient <- merge(prop_wide, dens_wide, by = c("PatientID","Age"), all = TRUE)

# ----------------------------
# 3) Morphology per patient+cluster
# ----------------------------
cellsclin[, AxisRatio := MajorAxisLength / MinorAxisLength]
morpho_summary <- cellsclin[
  , .(
    AvgCellArea  = mean(CellArea, na.rm = TRUE),
    AvgAxisRatio = mean(AxisRatio, na.rm = TRUE)
  ),
  by = .(PatientID, Age, ClusterID)
]

area_wide <- dcast(morpho_summary, PatientID + Age ~ ClusterID, value.var = "AvgCellArea", fill = NA_real_)
axis_wide <- dcast(morpho_summary, PatientID + Age ~ ClusterID, value.var = "AvgAxisRatio", fill = NA_real_)

area_cols <- setdiff(names(area_wide), c("PatientID","Age"))
axis_cols <- setdiff(names(axis_wide), c("PatientID","Age"))

area_new <- sapply(as.integer(area_cols), function(cid) paste("mean", label_for_cluster(cid), "cell area"))
axis_new <- sapply(as.integer(axis_cols), function(cid) paste("mean", label_for_cluster(cid), "cell axis ratio"))

setnames(area_wide, area_cols, area_new)
setnames(axis_wide, axis_cols, axis_new)

morpho_wide <- merge(area_wide, axis_wide, by = c("PatientID","Age"), all = TRUE)

# ----------------------------
# 4) Cluster-specific Ki67 positive fraction (isKi67Pos) per patient
# ----------------------------
if (!("isKi67Pos" %in% names(cellsclin))) {
  warning("isKi67Pos not found in cells.parquet; skipping cluster Ki67 fractions.")
  ki67_cluster_wide <- data.table(PatientID = unique(cellsclin$PatientID))
  ki67_cluster_wide[, Age := cellsclin[match(PatientID, PatientID), Age]]
} else {
  ki67_cluster <- cellsclin[
    , .(Ki67PosFrac = mean(as.logical(isKi67Pos), na.rm = TRUE)),
    by = .(PatientID, Age, ClusterID)
  ]
  ki67_cluster_wide <- dcast(ki67_cluster, PatientID + Age ~ ClusterID, value.var = "Ki67PosFrac", fill = 0)
  
  ki67_cols <- setdiff(names(ki67_cluster_wide), c("PatientID","Age"))
  ki67_new <- sapply(as.integer(ki67_cols), function(cid) paste(label_for_cluster(cid), "proliferative fraction"))
  setnames(ki67_cluster_wide, ki67_cols, ki67_new)
}

# ----------------------------
# 5) Epithelial-only curated call fractions
# ----------------------------
epi_markers <- c("isKi67Pos","isFOXA1Pos","isGATA3Pos","isERPos","isPRPos","isARPos","isHER2Pos")
epi_markers <- intersect(epi_markers, names(cellsclin))

epi_frac <- NULL
if (length(epi_markers) > 0) {
  epi_cells <- cellsclin[isEpithelial == TRUE]
  epi_frac <- epi_cells[
    , lapply(.SD, function(x) mean(as.logical(x), na.rm = TRUE)),
    by = .(PatientID, Age),
    .SDcols = epi_markers
  ]
  
  new_names <- c(
    "epithelium $Ki67^{+}$ fraction",
    "epithelium $FOXA1^{+}$ fraction",
    "epithelium $GATA3^{+}$ fraction",
    "epithelium $ER^{+}$ fraction",
    "epithelium $PR^{+}$ fraction",
    "epithelium $AR^{+}$ fraction",
    "epithelium $HER2^{+}$ fraction"
  )
  map <- setNames(new_names, c("isKi67Pos","isFOXA1Pos","isGATA3Pos","isERPos","isPRPos","isARPos","isHER2Pos"))
  setnames(epi_frac, epi_markers, unname(map[epi_markers]))
} else {
  warning("No epithelial curated calls found; skipping epithelial fractions.")
}

# ----------------------------
# 6) Interaction normalisation using cellAdjacency
# ----------------------------
cell_min <- cellsclin[, .(ImageID, CellID, isEpithelial, PatientID, Age)]
setkey(cell_min, ImageID, CellID)

edges <- copy(adj)
setkey(edges, ImageID, FromCellID)
edges <- edges[cell_min, on = .(ImageID = ImageID, FromCellID = CellID)]
setnames(edges, "isEpithelial", "IsEpithelial_source")

setkey(edges, ImageID, ToCellID)
edges <- edges[cell_min, on = .(ImageID = ImageID, ToCellID = CellID)]
setnames(edges, "isEpithelial", "IsEpithelial_target")

edges <- edges[!is.na(IsEpithelial_source) & !is.na(IsEpithelial_target)]
edges <- unique(edges, by = c("ImageID","FromCellID","ToCellID"))

epi_epi_int <- edges[IsEpithelial_source & IsEpithelial_target, .(interactions = .N), by = ImageID]
epi_tme_int <- edges[(IsEpithelial_source & !IsEpithelial_target) | (!IsEpithelial_source & IsEpithelial_target),
                     .(interactions = .N), by = ImageID]
tme_tme_int <- edges[!IsEpithelial_source & !IsEpithelial_target, .(interactions = .N), by = ImageID]

epi_counts_img <- cellsclin[isEpithelial == TRUE, .(counts = .N), by = .(ImageID, PatientID, Age)]
all_counts_img <- cellsclin[, .(counts = .N), by = .(ImageID, PatientID, Age)]
tme_counts_img <- cellsclin[isEpithelial == FALSE, .(counts = .N), by = .(ImageID, PatientID, Age)]

merge_int_norm <- function(count_dt, int_dt, out_col) {
  x <- merge(count_dt, int_dt, by = "ImageID", all.x = TRUE)
  x[is.na(interactions), interactions := 0L]
  
  agg <- x[, .(
    aggregated_counts = sum(counts),
    aggregated_interactions = sum(interactions)
  ), by = .(PatientID, Age)]
  
  agg[, (out_col) := fifelse(aggregated_counts > 0, aggregated_interactions / aggregated_counts, 0)]
  agg[, c("PatientID", "Age", out_col), with = FALSE]
}

epi_epi_norm <- merge_int_norm(epi_counts_img, epi_epi_int, "EpiEpi_norm")
epi_tme_norm <- merge_int_norm(all_counts_img, epi_tme_int, "EpiTME_norm")
tme_tme_norm <- merge_int_norm(tme_counts_img, tme_tme_int, "TMETME_norm")

interaction_cols <- Reduce(function(a,b) merge(a,b, by = c("PatientID","Age"), all = TRUE),
                           list(epi_epi_norm, epi_tme_norm, tme_tme_norm))

setnames(interaction_cols,
         c("EpiEpi_norm","EpiTME_norm","TMETME_norm"),
         c("normalised Epi-Epi interactions",
           "normalised Epi-TME interactions",
           "normalised TME-TME interactions"))

# ----------------------------
# 7) Mean distance from TME clusters to nearest epithelium (per image)
# ----------------------------
nn_available <- requireNamespace("FNN", quietly = TRUE)

dist_img_one <- function(imgID) {
  dat <- cellsclin[ImageID == imgID & !is.na(CenterX) & !is.na(CenterY),
                   .(CellID, ClusterID, isEpithelial, PatientID, Age, CenterX, CenterY)]
  if (nrow(dat) == 0) return(NULL)
  epi <- dat[isEpithelial == TRUE]
  tme <- dat[isEpithelial == FALSE]
  if (nrow(epi) == 0 || nrow(tme) == 0) return(NULL)
  
  epi_xy <- as.matrix(epi[, .(CenterX, CenterY)])
  tme_xy <- as.matrix(tme[, .(CenterX, CenterY)])
  
  if (nn_available) {
    nn <- FNN::get.knnx(epi_xy, tme_xy, k = 1)
    d <- as.numeric(nn$nn.dist[,1])
  } else {
    chunk <- 2000L
    d <- numeric(nrow(tme_xy))
    for (i in seq(1, nrow(tme_xy), by = chunk)) {
      j <- min(i + chunk - 1, nrow(tme_xy))
      block <- tme_xy[i:j, , drop = FALSE]
      d[i:j] <- apply(block, 1, function(p) {
        dx <- epi_xy[,1] - p[1]
        dy <- epi_xy[,2] - p[2]
        sqrt(min(dx*dx + dy*dy))
      })
    }
  }
  
  out <- tme[, .(PatientID, Age, ClusterID)]
  out[, shortest_distance := d]
  out
}

img_ids  <- unique(cellsclin$ImageID)
dist_tme <- rbindlist(lapply(img_ids, dist_img_one), fill = TRUE)

dist_summary <- dist_tme[, .(MeanDistToEpi = mean(shortest_distance, na.rm = TRUE)),
                         by = .(PatientID, Age, ClusterID)]
dist_wide <- dcast(dist_summary, PatientID + Age ~ ClusterID, value.var = "MeanDistToEpi", fill = NA_real_)

dist_cols <- setdiff(names(dist_wide), c("PatientID","Age"))
dist_new  <- sapply(as.integer(dist_cols), function(cid) paste("mean", label_for_cluster(cid), "distance to epithelium"))
setnames(dist_wide, dist_cols, dist_new)

# ----------------------------
# 8) Community stats via connected components within each SpatialClusterID
# ----------------------------
cell_sp <- cellsclin[, .(ImageID, CellID, SpatialClusterID, PatientID, Age)]
cell_sp <- cell_sp[!is.na(SpatialClusterID)]

edges_u <- unique(edges[, .(ImageID, a = FromCellID, b = ToCellID)])
community_rows <- list()

for (img in unique(cell_sp$ImageID)) {
  dat_img <- cell_sp[ImageID == img]
  if (nrow(dat_img) == 0) next
  
  by_sp <- split(dat_img$CellID, dat_img$SpatialClusterID)
  e_img <- edges_u[ImageID == img]
  if (nrow(e_img) == 0) next
  
  for (sid in names(by_sp)) {
    vids <- by_sp[[sid]]
    if (length(vids) < 2) next
    
    e_sub <- e_img[a %in% vids & b %in% vids, .(a,b)]
    if (nrow(e_sub) == 0) next
    
    g <- igraph::graph_from_data_frame(e_sub, directed = FALSE, vertices = data.frame(name = vids))
    comps <- igraph::components(g)$membership
    cc_dt <- data.table(CellID = as.integer(names(comps)), Component = as.integer(comps))
    cc_sizes <- cc_dt[, .(CommunitySize = .N), by = Component]
    
    community_rows[[length(community_rows) + 1L]] <- data.table(
      ImageID = img,
      SpatialClusterID = as.integer(sid),
      Component = cc_sizes$Component,
      CommunitySize = cc_sizes$CommunitySize
    )
  }
}

community_summary <- rbindlist(community_rows, fill = TRUE)
if (nrow(community_summary) == 0) {
  warning("No communities constructed from adjacency + SpatialClusterID; community features will be absent.")
  comm_feat <- data.table(PatientID = unique(cellsclin$PatientID), Age = unique(cellsclin$Age)[1])
} else {
  community_summary <- merge(community_summary, clin_min[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)
  
  community_counts <- community_summary[, .N, by = .(PatientID, Age, SpatialClusterID)]
  community_counts <- merge(community_counts, tissue_by_patient, by = "PatientID", all.x = TRUE)
  community_counts[, Density := fifelse(!is.na(totalTissueArea) & totalTissueArea > 0, N / totalTissueArea, NA_real_)]
  
  total_communities <- community_counts[, .(TotalCommunities = sum(N)), by = .(PatientID, Age)]
  community_counts <- merge(community_counts, total_communities, by = c("PatientID","Age"), all.x = TRUE)
  community_counts[, Proportion := fifelse(TotalCommunities > 0, N / TotalCommunities, NA_real_)]
  
  size_wide <- dcast(community_counts, PatientID + Age ~ SpatialClusterID, value.var = "N", fill = 0)
  dens_wide <- dcast(community_counts, PatientID + Age ~ SpatialClusterID, value.var = "Density", fill = NA_real_)
  prop_wide <- dcast(community_counts, PatientID + Age ~ SpatialClusterID, value.var = "Proportion", fill = NA_real_)
  
  scols <- setdiff(names(size_wide), c("PatientID","Age"))
  dcols <- setdiff(names(dens_wide), c("PatientID","Age"))
  pcols <- setdiff(names(prop_wide), c("PatientID","Age"))
  
  setnames(size_wide, scols, sapply(as.integer(scols), function(sid) paste0(label_for_spatial(sid), " community size")))
  setnames(dens_wide, dcols, sapply(as.integer(dcols), function(sid) paste0(label_for_spatial(sid), " communities per $mm^{2}$")))
  setnames(prop_wide, pcols, sapply(as.integer(pcols), function(sid) paste0(label_for_spatial(sid), " community proportion")))
  
  comm_feat <- Reduce(function(a,b) merge(a,b, by = c("PatientID","Age"), all = TRUE),
                      list(size_wide, dens_wide, prop_wide))
}

# ----------------------------
# 9) Duct/Lobule “structures” from epithelial TissueStructure + adjacency components
# ----------------------------
ductlob_rows <- list()
if ("TissueStructure" %in% names(ctx)) {
  epi_ctx <- merge(cellsclin[isEpithelial == TRUE, .(ImageID, CellID, PatientID, Age)],
                   ctx[, .(ImageID, CellID, TissueStructure)],
                   by = c("ImageID","CellID"), all.x = TRUE)
  epi_ctx <- epi_ctx[!is.na(TissueStructure) & TissueStructure %in% c("duct","lobule")]
  
  for (img in unique(epi_ctx$ImageID)) {
    dat_img <- epi_ctx[ImageID == img]
    e_img <- edges_u[ImageID == img]
    if (nrow(dat_img) == 0 || nrow(e_img) == 0) next
    
    for (lab in c("duct","lobule")) {
      vids <- dat_img[TissueStructure == lab, CellID]
      if (length(vids) < 2) next
      e_sub <- e_img[a %in% vids & b %in% vids, .(a,b)]
      if (nrow(e_sub) == 0) next
      
      g <- igraph::graph_from_data_frame(e_sub, directed = FALSE, vertices = data.frame(name = vids))
      comps <- igraph::components(g)$membership
      comp_sizes <- as.data.table(table(comps))
      setnames(comp_sizes, c("comps","Ncells"))
      
      ductlob_rows[[length(ductlob_rows) + 1L]] <- data.table(
        ImageID = img,
        PatientID = dat_img$PatientID[1],
        Age = dat_img$Age[1],
        TissueStructure = lab,
        CommunitySize = as.integer(comp_sizes$Ncells)
      )
    }
  }
}

ductlob_dt <- rbindlist(ductlob_rows, fill = TRUE)
if (nrow(ductlob_dt) > 0) {
  ductlob_counts <- ductlob_dt[, .N, by = .(PatientID, Age, TissueStructure)]
  ductlob_counts <- merge(ductlob_counts, tissue_by_patient, by = "PatientID", all.x = TRUE)
  ductlob_counts[, Density := N * 1e6 / totalTissueArea]
  
  mean_size_dt <- ductlob_dt[, .(MeanCommunitySize = mean(CommunitySize, na.rm = TRUE)),
                             by = .(PatientID, Age, TissueStructure)]
  ductlob_counts <- merge(ductlob_counts, mean_size_dt, by = c("PatientID","Age","TissueStructure"), all.x = TRUE)
  
  duct_dt <- ductlob_counts[TissueStructure == "duct",
                            .(PatientID, Age,
                              `duct structures per $mm^{2}$` = Density,
                              `cells per duct` = MeanCommunitySize)]
  lob_dt  <- ductlob_counts[TissueStructure == "lobule",
                            .(PatientID, Age,
                              `lobule structures per $mm^{2}$` = Density,
                              `cells per lobule` = MeanCommunitySize)]
  ductlob_feat <- merge(duct_dt, lob_dt, by = c("PatientID","Age"), all = TRUE)
} else {
  warning("Could not compute duct/lobule structures; skipping these features.")
  ductlob_feat <- data.table(PatientID = unique(cellsclin$PatientID), Age = unique(cellsclin$Age)[1])
}

# ----------------------------
# 10) Compartment summary (patient-level)
# ----------------------------
comp_feat <- as.data.table(comp_sum)
comp_feat[, PatientID := as.character(PatientID)]
comp_feat[, Age := as.numeric(Age)]

comp_keep <- intersect(
  c("PatientID","Age",
    "adiposeAreaFraction","adiposeMeanArea","adiposeMeanDistanceToEpithelium",
    "vascularAreaFraction","vascularMeanArea","vascularMeanDistanceToEpithelium",
    "lymphaticAreaFraction","lymphaticMeanArea","lymphaticMeanDistanceToEpithelium"),
  names(comp_feat)
)
comp_feat <- comp_feat[, ..comp_keep]

rename_comp <- c(
  adiposeAreaFraction = "adipose tissue area fraction",
  adiposeMeanArea = "mean adipocyte area",
  adiposeMeanDistanceToEpithelium = "mean distance adipocyte $\\rightarrow$ epithelium",
  vascularAreaFraction = "vascular tissue area fraction",
  vascularMeanArea = "mean blood vessel area",
  vascularMeanDistanceToEpithelium = "mean distance blood vessel $\\rightarrow$ epithelium",
  lymphaticAreaFraction = "lymphatic tissue area fraction",
  lymphaticMeanArea = "mean lymphatic vessel area",
  lymphaticMeanDistanceToEpithelium = "mean distance lymphatic vessel $\\rightarrow$ epithelium"
)
rename_comp <- rename_comp[names(rename_comp) %in% names(comp_feat)]
setnames(comp_feat, names(rename_comp), unname(rename_comp))

# ----------------------------
# Combine features into combined_dt (PatientID+Age)
# ----------------------------
combined_dt <- Reduce(function(a,b) merge(a,b, by = c("PatientID","Age"), all = TRUE),
                      list(
                        expr_by_patient,
                        propdens_by_patient,
                        morpho_wide,
                        ki67_cluster_wide,
                        if (!is.null(epi_frac)) epi_frac else data.table(PatientID=expr_by_patient$PatientID, Age=expr_by_patient$Age),
                        interaction_cols,
                        dist_wide,
                        comm_feat,
                        ductlob_feat,
                        comp_feat
                      ))

drop_present <- intersect(columns_to_remove, names(combined_dt))
if (length(drop_present) > 0) combined_dt[, (drop_present) := NULL]

# ----------------------------
# Myoepi thickness (mask-based)
# ----------------------------
compute_myoepi_patient_dt <- function(
    cells_dt,
    ctx_dt,
    clinical_dt,
    MASKDIR,
    cores = max(1L, parallel::detectCores() - 1L)
) {
  cells_dt    <- as.data.table(cells_dt)
  ctx_dt      <- as.data.table(ctx_dt)
  clinical_dt <- as.data.table(clinical_dt)
  
  epi <- merge(
    cells_dt,
    ctx_dt[, .(ImageID, CellID, TissueStructure)],
    by = c("ImageID", "CellID"),
    all.x = TRUE
  )
  epi <- merge(
    epi,
    clinical_dt[, .(ImageID, PatientID, Age, ImageHeight)],
    by = "ImageID",
    all.x = TRUE
  )
  
  epi <- epi[
    isEpithelial == TRUE &
      !is.na(TissueStructure) &
      !is.na(CenterX) & !is.na(CenterY) &
      !is.na(ImageHeight),
    .(ImageID, CellID, CenterX, CenterY, TissueStructure, ImageHeight, PatientID, Age)
  ]
  
  read_myo_mask <- function(image_id) {
    f <- file.path(MASKDIR, paste0(image_id, "MyoepiMask.tiff"))
    if (!file.exists(f)) return(NULL)
    m <- EBImage::readImage(f)
    m <- EBImage::normalize(m)
    m <- m > 0
    m <- EBImage::flip(m)
    m
  }
  
  process_one_image <- function(image_id) {
    mask <- read_myo_mask(image_id)
    if (is.null(mask)) return(NULL)
    
    dist_map <- EBImage::distmap(mask) * 2
    
    thickness_df <- as.data.frame(as.table(dist_map))
    colnames(thickness_df) <- c("x", "y", "thickness")
    thickness_df <- thickness_df[thickness_df$thickness > 0, , drop = FALSE]
    if (nrow(thickness_df) == 0) return(NULL)
    
    thickness_dt <- as.data.table(thickness_df)
    thickness_dt[, `:=`(
      x = as.numeric(x),
      y = as.numeric(y),
      thickness = as.numeric(thickness)
    )]
    
    epi_sub <- epi[ImageID == image_id]
    if (nrow(epi_sub) == 0) return(NULL)
    
    epi_sub[, CenterY_legacy := ImageHeight - CenterY]
    
    myo_ppp <- ppp(thickness_dt$x, thickness_dt$y,
                   c(0, max(thickness_dt$x)),
                   c(0, max(thickness_dt$y)))
    epi_ppp <- ppp(epi_sub$CenterX, epi_sub$CenterY_legacy,
                   c(0, max(thickness_dt$x)),
                   c(0, max(thickness_dt$y)))
    
    nn <- nncross(myo_ppp, epi_ppp)
    thickness_dt[, TissueStructure := epi_sub$TissueStructure[nn$which]]
    
    duct_dt   <- thickness_dt[TissueStructure == "duct"]
    lobule_dt <- thickness_dt[TissueStructure == "lobule"]
    
    data.table(
      ImageID = image_id,
      avg_thickness_duct   = if (nrow(duct_dt)   > 0) mean(duct_dt$thickness,   na.rm = TRUE) else NA_real_,
      area_duct            = nrow(duct_dt),
      avg_thickness_lobule = if (nrow(lobule_dt) > 0) mean(lobule_dt$thickness, na.rm = TRUE) else NA_real_,
      area_lobule          = nrow(lobule_dt)
    )
  }
  
  image_ids <- unique(epi$ImageID)
  
  img_res <- if (.Platform$OS.type == "unix") {
    parallel::mclapply(image_ids,
                       function(id) tryCatch(process_one_image(id), error = function(e) NULL),
                       mc.cores = cores
    )
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { library(data.table); library(EBImage); library(spatstat) })
    parallel::clusterExport(cl, varlist = c("epi","MASKDIR","read_myo_mask","process_one_image"),
                            envir = environment())
    parallel::parLapply(cl, image_ids, function(id) tryCatch(process_one_image(id), error = function(e) NULL))
  }
  
  image_dt <- rbindlist(img_res, fill = TRUE)
  if (nrow(image_dt) == 0) {
    return(data.table(
      PatientID = character(),
      Age = numeric(),
      `mean duct myoepi layer thickness` = numeric(),
      `mean lobule myoepi layer thickness` = numeric()
    ))
  }
  
  image_dt <- merge(
    image_dt,
    unique(epi[, .(ImageID, PatientID, Age)], by = "ImageID"),
    by = "ImageID",
    all.x = TRUE
  )
  
  patient_dt <- image_dt[, {
    valid_d <- !is.na(avg_thickness_duct) & !is.na(area_duct)
    d_th <- avg_thickness_duct[valid_d]
    d_ar <- area_duct[valid_d]
    duct_wmean <- if (length(d_th) > 0 && sum(d_ar, na.rm = TRUE) > 0) {
      sum(d_th * d_ar, na.rm = TRUE) / sum(d_ar, na.rm = TRUE)
    } else NA_real_
    
    valid_l <- !is.na(avg_thickness_lobule) & !is.na(area_lobule)
    l_th <- avg_thickness_lobule[valid_l]
    l_ar <- area_lobule[valid_l]
    lobule_wmean <- if (length(l_th) > 0 && sum(l_ar, na.rm = TRUE) > 0) {
      sum(l_th * l_ar, na.rm = TRUE) / sum(l_ar, na.rm = TRUE)
    } else NA_real_
    
    .(
      Age = Age[1],
      `mean duct myoepi layer thickness`   = duct_wmean,
      `mean lobule myoepi layer thickness` = lobule_wmean
    )
  }, by = PatientID]
  
  patient_dt[]
}

myo_patient <- compute_myoepi_patient_dt(
  cells_dt    = cells,
  ctx_dt      = ctx,
  clinical_dt = clinical,
  MASKDIR     = MASK_DIR
)

combined_dt <- merge(
  combined_dt,
  myo_patient[, .(PatientID, Age,
                  `mean duct myoepi layer thickness`,
                  `mean lobule myoepi layer thickness`)],
  by = c("PatientID","Age"),
  all.x = TRUE
)

# ----------------------------
# Build DEswan inputs
# ----------------------------
feature_cols <- setdiff(names(combined_dt), c("PatientID","Age"))
data.df      <- combined_dt[, ..feature_cols]
qt           <- combined_dt$Age


set.seed(13)
qt_perm <- sample(qt)

res_perm <- DEswan(
  data.df       = as.data.frame(data.df),
  qt            = qt_perm,
  window.center = seq(20, 70, 1),
  buckets.size  = 10
)
p_perm_wide <- reshape.DEswan(res_perm, parameter = 1, factor = "qt")
q_perm_wide <- q.DEswan(p_perm_wide, method = "BH")
setDT(q_perm_wide)

signif_counts_perm <- q_perm_wide[, lapply(.SD, function(x) sum(x < alpha, na.rm = TRUE)), .SDcols = -1]
vperm <- as.numeric(signif_counts_perm[1, ])
names(vperm) <- gsub("^X", "", names(signif_counts_perm))

plot_data_perm <- data.table(
  Age_Window_Center    = as.numeric(names(vperm)),
  Significant_Features = vperm
)

p_perm <- ggplot(plot_data_perm, aes(x = Age_Window_Center, y = Significant_Features)) +
  geom_point(size = 2, color = "black") +
  geom_line(size = 1, color = "black") +
  theme_classic() +
  labs(
    subtitle = "Permuted patient identities\n(bucket size = 10, q < 0.05)",
    x = "Age (years)",
    y = "No. of significant features"
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 70), limits = c(0, 70)) +
  theme(
    axis.title = element_text(color = "black", size = 13),
    axis.text  = element_text(color = "black", size = 13),
    plot.subtitle = element_text(color = "black", size = 13),
    plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm")
  )

out_perm <- file.path(OUT_DIR, "significantFeaturesAgePermuted.pdf")
ggsave(p_perm, filename = out_perm, units = "in", width = 4.5, height = 3)
cat("Wrote: ", out_perm, "\n", sep = "")
