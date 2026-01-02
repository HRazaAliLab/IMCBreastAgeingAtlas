library(here)
source(here("code", "header.R"))

#Dependencies below needed if community_summary.csv not already created.
# library(deldir)
# library(igraph)

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

community_summary_csv  <- file.path(outdir, "community_summary.csv")
community_summary_path <- file.path(outdir, "community_summary.parquet")

cells    <- getCells()
clinical <- getClinical()

clinical[, PatientID := as.character(PatientID)]
clinical[, Age := as.numeric(Age)]

tissueAreasPatients <- clinical[, .(totalTissueArea = sum(TissueArea, na.rm = TRUE)), by = PatientID]

communityannotations <- fread(here("data", "derived", "spatialClusterAnnotation.csv"))
communityannotations <- as.data.table(communityannotations)
communityannotations[, SpatialClusterID := as.character(SpatialClusterID)]

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

#load if exists; otherwise compute + cache to CSV (requires a lot of memory to generate)
if (file.exists(community_summary_csv)) {
  community_summary <- fread(community_summary_csv)
} else if (file.exists(community_summary_path)) {
  community_summary <- as.data.table(read_parquet_adt(community_summary_path))
} else {
  message("community_summary not found; computing now (this can be heavy).")
  community_summary <- analyze_clusters_parallel(cells)
  fwrite(community_summary, community_summary_csv)
  if (exists("write_parquet_adt")) write_parquet_adt(community_summary, community_summary_path)
}

community_summary <- as.data.table(community_summary)
community_summary <- community_summary[CommunitySize >= 10]
community_summary[, SpatialClusterID := as.character(SpatialClusterID)]

spatial_merged_data <- merge(
  community_summary,
  clinical[, .(ImageID, PatientID, Age)],
  by = "ImageID",
  all.x = TRUE
)

all_combinations <- CJ(
  PatientID        = unique(spatial_merged_data$PatientID),
  SpatialClusterID = unique(spatial_merged_data$SpatialClusterID),
  sorted = FALSE
)

community_counts <- spatial_merged_data[, .N, by = .(PatientID, SpatialClusterID)]
community_counts <- merge(all_combinations, community_counts,
                          by = c("PatientID", "SpatialClusterID"),
                          all.x = TRUE)
community_counts[is.na(N), N := 0]

community_counts <- merge(community_counts, tissueAreasPatients, by = "PatientID", all.x = TRUE)
community_counts[, Density := N * 1e6 / totalTissueArea]

agemapping <- unique(spatial_merged_data[, .(PatientID, Age)])
spatial_merged_data2 <- merge(community_counts, agemapping, by = "PatientID", all.x = TRUE)
spatial_merged_data2[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
spatial_merged_data2 <- na.omit(spatial_merged_data2)

setnames(spatial_merged_data2, "N", "NumberOfCommunities")

summarize_pvals_only <- function(dt) {
  clusters <- unique(dt$SpatialClusterID)
  out <- data.table(Cluster = character(), pValue = numeric())
  
  for (cluster in clusters) {
    cluster_dt <- dt[SpatialClusterID == cluster]
    if (nrow(cluster_dt) == 0) next
    test_result <- wilcox.test(Density ~ AgeGroup, data = cluster_dt)
    out <- rbind(out, data.table(Cluster = as.character(cluster), pValue = test_result$p.value))
  }
  
  out[, adjpValue := p.adjust(pValue, method = "BH")]
  out
}

spatial_summary_results <- summarize_pvals_only(spatial_merged_data2)

createClusterBoxplot <- function(cluster_id, ymax) {
  cluster_id <- as.character(cluster_id)
  
  cluster_specific_data <- spatial_merged_data2[SpatialClusterID == cluster_id]
  cluster_specific_data[, Age_Group := ifelse(Age >= 50, "postmenopausal", "premenopausal")]
  cluster_specific_data[, Age_Group := factor(Age_Group, levels = c("premenopausal", "postmenopausal"))]
  cluster_specific_data <- na.omit(cluster_specific_data, cols = "Age_Group")
  
  label <- communityannotations[SpatialClusterID == cluster_id, Label][1]
  adjPVal <- spatial_summary_results[Cluster == cluster_id, adjpValue][1]
  
  if (is.na(label)) label <- cluster_id
  if (is.na(adjPVal)) adjPVal <- NA_real_
  
  ptxt <- if (is.na(adjPVal)) {
    TeX("p=NA")
  } else {
    TeX(paste0("p=", mkEnumPower(ifelse(adjPVal < 0.001,
                                        formatC(adjPVal, format = "e", digits = 0),
                                        sprintf("%.3f", adjPVal)))))
  }
  
  ggplot(cluster_specific_data, aes(x = Age_Group, y = Density, fill = Age_Group)) +
    geom_boxplot(outlier.size = 0.1) +
    scale_fill_manual(values = c("premenopausal" = "gray", "postmenopausal" = "dimgray")) +
    theme_classic() +
    scale_x_discrete(labels = c("pre", "post")) +
    scale_y_continuous(
      limits = c(0, ymax),
      breaks = c(0, floor(ymax / 0.05) * 0.05),
      expand = expand_scale(mult = c(0, 0.1)),
      labels = function(x) ifelse(x == 0, "0", format(x))
    ) +
    labs(x = NULL, y = NULL, title = TeX(label)) +
    theme(
      axis.text.x = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.y = element_text(size = 23),
      axis.ticks.x = element_blank(),
      plot.title = element_text(size = 25, hjust = 0.5),
      legend.position = "none",
      plot.margin = unit(c(10,40,20,0), "pt")
    ) +
    annotate(
      "text",
      x = 1.5, y = ymax * 0.96,
      label = ptxt,
      parse = TRUE, size = 8, color = "black", hjust = 0.5
    )
}

spatialboxplot1 <- createClusterBoxplot("1", ymax = 20)
spatialboxplot2 <- createClusterBoxplot("0", ymax = 10)

spatialcombined <- spatialboxplot1 / spatialboxplot2 + plot_layout(nrow = 2)

spatialblanklabelploty <- ggplot() +
  labs(y = TeX("density (communities/mm$^2$)")) +
  theme_classic() +
  theme(
    plot.margin = margin(0,10,0,0, unit = "pt"),
    axis.title.y = element_text(size = 23)
  ) +
  guides(x = "none", y = "none")

spatialcombined <- spatialblanklabelploty + spatialcombined + plot_layout(widths = c(1, 1000))

pdf(file.path(outdir, "SpatialPrePostMenopauseBoxplots.pdf"), width = 3.5, height = 5)
print(spatialcombined)
dev.off()
