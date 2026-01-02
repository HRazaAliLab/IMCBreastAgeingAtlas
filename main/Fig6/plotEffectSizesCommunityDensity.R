library(here)
source(here("code", "header.R"))

#Dependencies below needed if community_summary.csv not already created.
# library(deldir)
# library(igraph)

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
clinical <- getClinical()

clinical[, PatientID := as.character(PatientID)]
clinical[, Age := as.numeric(Age)]

tissueAreasPatients <- clinical[, .(totalTissueArea = sum(TissueArea, na.rm = TRUE)), by = PatientID]

communityannotations <- fread(here("data", "derived", "spatialClusterAnnotation.csv"))
if (!("SpatialClusterID" %in% names(communityannotations)) && ("ClusterID" %in% names(communityannotations))) {
  setnames(communityannotations, "ClusterID", "SpatialClusterID")
}
communityannotations <- as.data.table(communityannotations)
communityannotations[, SpatialClusterID := as.character(SpatialClusterID)]

if (!("PrintOrder" %in% names(communityannotations))) {
  communityannotations[, PrintOrder := suppressWarnings(as.integer(SpatialClusterID))]
}
comm_ord <- communityannotations[order(PrintOrder)]

spatial_levels <- as.character(comm_ord$SpatialClusterID)
clusternames   <- setNames(as.character(comm_ord$Label),  as.character(comm_ord$SpatialClusterID))
colors         <- setNames(as.character(comm_ord$Colour), as.character(comm_ord$SpatialClusterID))

community_summary_path <- file.path(outdir, "community_summary.parquet")
community_summary_csv  <- file.path(outdir, "community_summary.csv")

analyze_clusters_parallel <- function(cells_dt) {
  cells_dt <- as.data.table(cells_dt)
  
  req <- c("ImageID", "CellID", "CenterX", "CenterY", "SpatialClusterID", "ClusterID")
  missing <- setdiff(req, names(cells_dt))
  if (length(missing) > 0) stop("cells missing required columns: ", paste(missing, collapse = ", "))
  
  cells_dt[, CellID := as.character(CellID)]
  cells_dt[, SpatialClusterID := as.character(SpatialClusterID)]
  
  image_ids <- unique(cells_dt$ImageID)
  
  process_image <- function(image_id) {
    image_data <- cells_dt[ImageID == image_id]
    if (nrow(image_data) < 3) return(NULL)
    
    dt <- deldir(image_data$CenterX, image_data$CenterY)
    edges <- with(dt, data.frame(
      from = image_data$CellID[delsgs$ind1],
      to   = image_data$CellID[delsgs$ind2]
    ))
    
    g <- graph_from_data_frame(edges, directed = FALSE)
    
    image_results <- list()
    
    for (cluster in unique(image_data$SpatialClusterID)) {
      vids <- which(V(g)$name %in% as.character(image_data[SpatialClusterID == cluster, CellID]))
      if (length(vids) == 0) next
      subg <- induced_subgraph(g, vids = vids)
      
      comps <- components(subg)
      
      for (i in unique(comps$membership)) {
        members <- V(subg)[comps$membership == i]$name
        if (length(members) == 0) next
        
        pgclust_values <- image_data[CellID %in% members, ClusterID]
        
        image_results[[length(image_results) + 1]] <- data.table(
          ImageID          = image_id,
          SpatialClusterID = cluster,
          CommunitySize    = length(members),
          ShannonDiversity = calculate_shannon_diversity(pgclust_values),
          AvgVertexDegree  = mean(degree(subg, v = V(subg)[comps$membership == i]))
        )
      }
    }
    
    if (length(image_results) == 0) return(NULL)
    rbindlist(image_results, use.names = TRUE, fill = TRUE)
  }
  
  cores <- max(1, detectCores() - 1)
  results <- mclapply(image_ids, process_image, mc.cores = cores)
  rbindlist(results, use.names = TRUE, fill = TRUE)
}

if (file.exists(community_summary_path)) {
  community_summary <- as.data.table(read_parquet_adt(community_summary_path))
} else if (file.exists(community_summary_csv)) {
  community_summary <- fread(community_summary_csv)
} else {
  message("community_summary not found; computing now (this can be heavy).")
  community_summary <- analyze_clusters_parallel(cells)
  if (exists("write_parquet_adt")) write_parquet_adt(community_summary, community_summary_path)
  fwrite(community_summary, community_summary_csv)
}

community_summary <- as.data.table(community_summary)
community_summary <- community_summary[CommunitySize >= 10]
community_summary[, SpatialClusterID := as.character(SpatialClusterID)]

prepareEffectSizeData <- function(community_summary, clinical, tissueAreasPatients, spatial_levels) {
  community_summary <- merge(
    community_summary,
    clinical[, .(ImageID, PatientID, Age)],
    by = "ImageID",
    all.x = TRUE
  )
  
  community_summary[, SpatialClusterID := factor(SpatialClusterID, levels = spatial_levels)]
  
  all_combinations <- CJ(
    PatientID        = unique(community_summary$PatientID),
    SpatialClusterID  = factor(spatial_levels, levels = spatial_levels),
    sorted = FALSE
  )
  
  community_counts <- community_summary[, .N, by = .(PatientID, SpatialClusterID)]
  
  community_counts <- merge(
    all_combinations, community_counts,
    by = c("PatientID", "SpatialClusterID"),
    all.x = TRUE
  )
  community_counts[is.na(N), N := 0]
  
  community_counts <- merge(community_counts, tissueAreasPatients, by = "PatientID", all.x = TRUE)
  community_counts[, Density := N * 1e6 / totalTissueArea]
  
  agemapping <- unique(community_summary[, .(PatientID, Age)])
  merged_data <- merge(community_counts, agemapping, by = "PatientID", all.x = TRUE)
  merged_data[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
  
  merged_data <- na.omit(merged_data)
  setnames(merged_data, "N", "NumberOfCommunities")
  merged_data
}

summarizeEffectSize <- function(dt) {
  clusters <- levels(dt$SpatialClusterID)
  
  out <- data.table(
    Cluster    = character(),
    CliffDelta = numeric(),
    LowerCI    = numeric(),
    UpperCI    = numeric(),
    pValue     = numeric(),
    AvgBelow50 = numeric(),
    SEBelow50  = numeric(),
    AvgAbove50 = numeric(),
    SEAbove50  = numeric(),
    nBelow50   = integer(),
    nAbove50   = integer()
  )
  
  for (cluster in clusters) {
    cluster_dt <- dt[SpatialClusterID == cluster]
    if (nrow(cluster_dt) == 0) next
    
    below <- cluster_dt[AgeGroup == "Below 50",
                        .(Mean = mean(Density), SE = sd(Density) / sqrt(.N), n = .N)]
    above <- cluster_dt[AgeGroup == "Above 50",
                        .(Mean = mean(Density), SE = sd(Density) / sqrt(.N), n = .N)]
    
    test_result <- wilcox.test(Density ~ AgeGroup, data = cluster_dt)
    delta <- cliff.delta(Density ~ AgeGroup, data = cluster_dt)
    
    out <- rbind(out, data.table(
      Cluster    = as.character(cluster),
      CliffDelta = as.numeric(delta$estimate),
      LowerCI    = delta$conf.int[1],
      UpperCI    = delta$conf.int[2],
      pValue     = test_result$p.value,
      AvgBelow50 = below$Mean, SEBelow50 = below$SE,
      AvgAbove50 = above$Mean, SEAbove50 = above$SE,
      nBelow50   = below$n,    nAbove50  = above$n
    ))
  }
  
  out[, adjpValue := p.adjust(pValue, method = "BH")]
  out[, alpha := ifelse(adjpValue < 0.05, "opaque", "translucent")]
  out
}

plotEffectSize <- function(summary_dt, clusternames, spatial_levels) {
  summary_dt[, Cluster := factor(Cluster, levels = spatial_levels)]
  
  ggplot(summary_dt, aes(x = Cluster, y = CliffDelta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
    geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, alpha = alpha), width = 0.2) +
    geom_point(shape = 21, color = "black", aes(alpha = alpha), fill = "lightblue", size = 5) +
    coord_flip(clip = "off") +
    labs(x = NULL, y = "Cliff's delta effect size") +
    theme_classic() +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x  = element_text(colour = "black", size = 13),
      axis.title.x = element_text(size = 13),
      plot.margin  = unit(c(80, 100, 5, 0), "pt"),
      legend.position = "none"
    ) +
    scale_x_discrete(labels = sapply(clusternames, TeX))
}

spatial_merged_data <- prepareEffectSizeData(
  community_summary, clinical, tissueAreasPatients, spatial_levels
)

spatial_summary_results <- summarizeEffectSize(spatial_merged_data)

if (exists("format_custom_pval") && exists("mkEnumPower")) {
  spatial_summary_results[, adjpValueLabel := sapply(format_custom_pval(adjpValue), mkEnumPower)]
} else {
  spatial_summary_results[, adjpValueLabel := sprintf("%.3g", adjpValue)]
}

spatial_summary_results <- spatial_summary_results[order(-abs(CliffDelta))]

plot_levels <- spatial_summary_results[order(-CliffDelta), Cluster]
plot_levels <- as.character(plot_levels)

spatial_summary_results[, Cluster := factor(Cluster, levels = plot_levels)]

labels_map <- clusternames  # names are SpatialClusterID

p_eff <- ggplot(spatial_summary_results, aes(x = Cluster, y = CliffDelta)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  scale_alpha_manual(values = c("opaque" = 1, "translucent" = 0.5)) +
  geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI, alpha = alpha), width = 0.2) +
  geom_point(shape = 21, color = "black", aes(fill = alpha), size = 5) +
  scale_fill_manual(values = c("opaque" = "lightblue", "translucent" = "#3D3D3D")) +
  coord_flip(
    clip = "off",
    ylim = c(min(spatial_summary_results$LowerCI, na.rm = TRUE),
             max(spatial_summary_results$UpperCI, na.rm = TRUE))
  ) +
  labs(x = NULL, y = "Cliff's delta effect size") +
  theme_classic() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm"),
    axis.text.x  = element_text(colour = "black", size = 13),
    axis.title.x = element_text(size = 13),
    plot.margin  = unit(c(80, 100, 5, 0), "pt"),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = sapply(labels_map, TeX)) +
  geom_text(
    aes(
      y = 0.3,
      label = TeX(adjpValueLabel, output = "character"),
      color = ifelse(alpha == "translucent", "darkgray", "black")
    ),
    parse = TRUE, hjust = 0, size = 4.5
  ) +
  scale_color_identity() +
  scale_y_continuous(breaks = c(-0.6, 0, 0.2))

plot_for_colorbar <- merge(
  spatial_summary_results[, .(Cluster = as.character(Cluster), CliffDelta)],
  communityannotations[, .(SpatialClusterID, Label, Colour)],
  by.x = "Cluster", by.y = "SpatialClusterID",
  all.x = TRUE
)
plot_for_colorbar <- as.data.table(plot_for_colorbar)
plot_for_colorbar <- plot_for_colorbar[order(CliffDelta)]
plot_for_colorbar[, Cluster := factor(Cluster, levels = as.character(Cluster))]

clusternames_cb <- setNames(plot_for_colorbar$Label, plot_for_colorbar$Cluster)
colors_cb       <- setNames(plot_for_colorbar$Colour, plot_for_colorbar$Cluster)

colorbar_plot <- ggplot(plot_for_colorbar, aes(x = factor(1), y = Cluster, fill = Cluster)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_y_discrete(
    limits = rev(levels(plot_for_colorbar$Cluster)),
    labels = rev(sapply(clusternames_cb, TeX))
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_fill_manual(values = colors_cb) +
  theme_minimal() +
  theme(
    axis.title  = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0,0,0,0), "cm")
  )

prepostcombinedplot <- colorbar_plot + p_eff + plot_layout(widths = c(0.3, 7))

pdf(file.path(outdir, "effectSizesCommunityDensity.pdf"), width = 4.5, height = 4)
print(prepostcombinedplot)
dev.off()