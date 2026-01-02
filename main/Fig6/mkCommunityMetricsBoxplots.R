library(here)
source(here("code", "header.R"))

#Additional dependencies required for this analysis.
library(deldir)
library(igraph)

outdir <- here("scratch", "main", "Fig6")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells <- as.data.table(getCells())

communityannotations <- fread(here("data", "derived", "spatialClusterAnnotation.csv"))

communityannotations <- as.data.table(communityannotations)

cells[, SpatialClusterID := as.character(SpatialClusterID)]
communityannotations[, SpatialClusterID := as.character(SpatialClusterID)]

cells[, SpatialClusterID_num := suppressWarnings(as.integer(SpatialClusterID))]
communityannotations[, SpatialClusterID_num := suppressWarnings(as.integer(SpatialClusterID))]

if (anyNA(cells$SpatialClusterID)) stop("cells$SpatialClusterID has NA after coercion.")
if (anyNA(communityannotations$SpatialClusterID)) stop("communityannotations$SpatialClusterID has NA after coercion.")

community_summary_path <- file.path(outdir, "community_summary.parquet")
community_summary_csv  <- file.path(outdir, "community_summary.csv")

calculate_shannon_diversity <- function(data) {
  if (length(data) > 1) {
    freq <- table(data) / length(data)
    freq <- freq[freq > 0]
    return(-sum(freq * log(freq)))
  } else {
    return(0)
  }
}

#community summary computation
analyze_clusters_parallel <- function(cells_dt) {
  cells_dt <- as.data.table(cells_dt)
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
      subg <- induced_subgraph(g, vids = vids)
      
      comps <- components(subg)
      
      for (i in unique(comps$membership)) {
        members <- V(subg)[comps$membership == i]$name
        
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
    
    rbindlist(image_results, use.names = TRUE, fill = TRUE)
  }
  
  cores <- max(1, detectCores() - 1)
  results <- mclapply(image_ids, process_image, mc.cores = cores)
  rbindlist(results, use.names = TRUE, fill = TRUE)
}

#load cached community_summary if exists; otherwise compute - requires a lot of memory
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
community_summary[, SpatialClusterID := as.character(SpatialClusterID)]
community_summary <- community_summary[CommunitySize >= 10]

avg_vertex_degree <- community_summary[
  , .(AvgVertexDegree = median(AvgVertexDegree, na.rm = TRUE)),
  by = SpatialClusterID
]
avg_vertex_degree[, SpatialClusterID_num := suppressWarnings(as.integer(SpatialClusterID))]
setorder(avg_vertex_degree, AvgVertexDegree, SpatialClusterID_num)

order_ids <- avg_vertex_degree$SpatialClusterID
plot_levels <- rev(order_ids)

community_summary[, SpatialClusterID := factor(SpatialClusterID, levels = plot_levels)]

communityannotations_ordered <- communityannotations[match(order_ids, communityannotations$SpatialClusterID)]

communitylabels <- setNames(as.character(communityannotations_ordered$Label),
                            as.character(communityannotations_ordered$SpatialClusterID))
color_mapping   <- setNames(as.character(communityannotations_ordered$Colour),
                            as.character(communityannotations_ordered$SpatialClusterID))

label_fun <- function(x) sapply(as.character(x), function(id) TeX(communitylabels[[id]]))

cellscounts <- ggplot(community_summary, aes(x = SpatialClusterID, y = CommunitySize, fill = SpatialClusterID)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = NULL, subtitle = "cells per\ncommunity") +
  scale_x_discrete(labels = label_fun) +
  scale_y_log10(expand = c(0,0), limits = c(10, 1000), breaks = c(10, 100, 1000)) +
  scale_fill_manual(values = color_mapping, breaks = names(color_mapping)) +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = "black", size = 13),
    axis.text.x = element_text(color = "black", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    plot.margin = unit(c(20,15,0,0), "pt"),
    legend.position = "none"
  )

diversity <- ggplot(community_summary, aes(x = SpatialClusterID, y = ShannonDiversity, fill = SpatialClusterID)) +
  geom_boxplot(outlier.size = 0.1) +
  coord_flip(clip = "off") +
  labs(x = NULL, y = NULL, subtitle = "Shannon\ndiversity") +
  scale_x_discrete(labels = label_fun) +
  scale_y_continuous(expand = c(0,0), breaks = c(0,1,2,3), limits = c(0, 3)) +
  scale_fill_manual(values = color_mapping, breaks = names(color_mapping)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 13),
    plot.margin = unit(c(20,5,0,10), "pt"),
    legend.position = "none"
  )

degree <- ggplot(community_summary, aes(x = SpatialClusterID, y = AvgVertexDegree, fill = SpatialClusterID)) +
  geom_boxplot(outlier.size = 0.01) +
  coord_flip() +
  labs(x = NULL, y = NULL, subtitle = "average\nvertex degree") +
  scale_x_discrete(labels = label_fun) +
  scale_y_continuous(expand = c(0,0), breaks = c(2, 4, 6), limits = c(1.667, 6.333)) +
  scale_fill_manual(values = color_mapping, breaks = names(color_mapping)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 13),
    plot.margin = unit(c(20,5,0,5), "pt"),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "none"
  )

community_counts <- community_summary[, .(NumberOfCommunities = .N), by = .(SpatialClusterID)]

bar_levels <- rev(levels(community_summary$SpatialClusterID))
community_counts[, SpatialClusterID := factor(as.character(SpatialClusterID), levels = bar_levels)]

communitycounts <- ggplot(community_counts, aes(x = SpatialClusterID, y = NumberOfCommunities)) +
  geom_bar(stat = "identity", fill = "lightgray") +
  coord_flip() +
  labs(x = NULL, y = NULL, subtitle = "# of\ncommunities") +
  theme_classic() +
  scale_x_discrete(labels = label_fun) +
  scale_y_continuous(expand = c(0,0), breaks = c(0, 5000)) +
  theme(
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(color = "black", size = 13),
    plot.margin = unit(c(20,5,0,5), "pt"),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    legend.position = "none"
  )

combined <- cellscounts + diversity + degree + communitycounts +
  plot_layout(nrow = 1, widths = c(3,3,3,2.5)) +
  plot_annotation(theme = theme(plot.margin = unit(c(0,20,0,0), "pt")))

ggsave(
  combined,
  filename = file.path(outdir, "communityMetricsBoxplots.R.pdf"),
  width = 6.5, height = 3, units = "in"
)
