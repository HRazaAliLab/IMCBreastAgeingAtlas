library(here)
source(here("code", "header.R"))

#Dependencies required to generate community_summary_structures.csv.
library(deldir)
library(igraph)

outdir <- here("scratch", "main", "Fig7")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cache_csv <- file.path(outdir, "community_summary_structures.csv")

cells    <- getCells()
ctx      <- getCellContext()
clinical <- getClinical()

tissueAreasPatients <- clinical[, .(totalTissueArea = sum(TissueArea, na.rm = TRUE)), by = PatientID]

cells <- merge(cells,
               ctx[, .(ImageID, CellID, TissueStructure)],
               by = c("ImageID", "CellID"),
               all.x = TRUE
)

cellsclin <- merge(cells,
                   clinical[, .(ImageID, PatientID, Age)],
                   by = "ImageID",
                   all.x = TRUE
)

cellsclin[, CellID := as.character(CellID)]
cellsclin <- cellsclin[!is.na(CenterX) & !is.na(CenterY)]

analyze_clusters_parallel <- function(cellsclin_dt) {
  cellsclin_dt <- as.data.table(cellsclin_dt)
  cellsclin_dt[, CellID := as.character(CellID)]
  image_ids <- unique(cellsclin_dt$ImageID)
  
  process_image <- function(image_id) {
    image_data <- cellsclin_dt[ImageID == image_id]
    if (nrow(image_data) < 3) return(NULL)
    
    dt <- deldir(image_data$CenterX, image_data$CenterY)
    
    edges <- with(dt, data.frame(
      from = image_data$CellID[delsgs$ind1],
      to   = image_data$CellID[delsgs$ind2],
      stringsAsFactors = FALSE
    ))
    if (nrow(edges) == 0) return(NULL)
    
    g <- graph_from_data_frame(edges, directed = FALSE)
    
    image_results <- list()
    for (label in unique(image_data$TissueStructure)) {
      if (is.na(label)) next
      
      vertices_in_label <- which(V(g)$name %in% image_data[TissueStructure == label, CellID])
      if (length(vertices_in_label) == 0) next
      
      subg <- induced_subgraph(g, vids = vertices_in_label)
      comps <- components(subg)
      
      for (k in unique(comps$membership)) {
        members <- V(subg)[comps$membership == k]$name
        if (length(members) == 0) next
        
        image_results[[length(image_results) + 1]] <- data.table(
          ImageID = image_id,
          TissueStructure = label,
          CommunitySize = length(members)
        )
      }
    }
    
    if (length(image_results) == 0) return(NULL)
    rbindlist(image_results, use.names = TRUE)
  }
  
  cores <- max(1, parallel::detectCores() - 1)
  if (.Platform$OS.type == "unix") {
    rbindlist(parallel::mclapply(image_ids, process_image, mc.cores = cores), fill = TRUE)
  } else {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { library(deldir); library(igraph); library(data.table) })
    parallel::clusterExport(cl, varlist = c("cellsclin_dt", "process_image"), envir = environment())
    res <- parallel::parLapply(cl, image_ids, process_image)
    rbindlist(res, fill = TRUE)
  }
}

#build community summary (requires a lot of memory) or load from memory
if (file.exists(cache_csv)) {
  community_summary <- fread(cache_csv)
} else {
  message("Computing community_summary (heavy).")
  community_summary <- analyze_clusters_parallel(cellsclin)
  fwrite(community_summary, cache_csv)
}
community_summary <- as.data.table(community_summary)

community_summary <- merge(
  community_summary,
  unique(clinical[, .(ImageID, PatientID, Age)]),
  by = "ImageID",
  all.x = TRUE
)
community_summary <- community_summary[TissueStructure != "microenvironment"]

prepareEffectSizeData <- function(community_summary, tissueAreasPatients_dt) {
  
  all_combinations <- CJ(
    PatientID = unique(community_summary$PatientID),
    TissueStructure = unique(community_summary$TissueStructure)
  )
  
  community_counts <- community_summary[, .N, by = .(PatientID, TissueStructure)]
  community_counts <- merge(
    all_combinations, community_counts,
    by = c("PatientID", "TissueStructure"), all.x = TRUE
  )
  community_counts[is.na(N), N := 0L]
  
  community_counts <- merge(
    community_counts,
    tissueAreasPatients_dt[, .(PatientID, totalTissueArea)],
    by = "PatientID",
    all.x = TRUE
  )
  community_counts[, Density := N * 1e6 / totalTissueArea]
  
  mean_community_size <- community_summary[, .(
    MeanCommunitySize = mean(CommunitySize, na.rm = TRUE)
  ), by = .(PatientID, TissueStructure)]
  
  community_counts <- merge(
    community_counts, mean_community_size,
    by = c("PatientID", "TissueStructure"), all.x = TRUE
  )
  community_counts[is.na(MeanCommunitySize), MeanCommunitySize := 0]
  
  agemapping <- unique(community_summary[, .(PatientID, Age)])
  merged_data <- merge(community_counts, agemapping, by = "PatientID", all.x = TRUE)
  
  merged_data[, AgeGroup := ifelse(Age >= 50, "Above 50", "Below 50")]
  merged_data <- na.omit(merged_data)
  
  setnames(merged_data, "N", "NumberOfCommunities")
  merged_data[]
}

community_merged_data <- prepareEffectSizeData(community_summary, tissueAreasPatients)

cor_results_structures <- list()
for (structure in unique(community_merged_data$TissueStructure)) {
  structure_data <- community_merged_data[TissueStructure == structure]
  cor_test <- cor.test(structure_data$Age, structure_data$Density, method = "spearman")
  cor_results_structures[[structure]] <- list(rho = cor_test$estimate, p.value = cor_test$p.value)
}

ductdata <- community_merged_data[TissueStructure == "duct"]
ductagecommunityplot <- ggplot(ductdata, aes(x = Age, y = Density, color = TissueStructure)) +
  geom_point(size = 0.1, alpha = 0.45) +
  geom_smooth(aes(group = TissueStructure), method = "loess", se = FALSE) +
  scale_color_manual(values = c("duct" = "#4169E1", "lobule" = "#DC143C")) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 10), expand = c(0.001, 0.01), breaks = c(0, 10), labels = c(0, 10)) +
  scale_x_continuous(expand = c(0.0001, 0.01), limits = c(15, 86), breaks = c(20, 50, 80)) +
  labs(x = NULL, y = TeX("density (structures/mm$^2$)")) +
  theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16),
    plot.margin = unit(c(10, 5, 0, 0), units = "pt"),
    legend.position = "none"
  ) +
  annotate("text", x = 16, y = 9.5,
           label = TeX(paste0("$\\rho$ = ", format(cor_results_structures[["duct"]]$rho, digits = 2))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 8.75,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(cor_results_structures[["duct"]]$p.value)))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 8,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(ductdata))),
           size = 5.5, color = "black", hjust = 0)

avg_density_over_50  <- mean(ductdata[Age >= 50]$Density, na.rm = TRUE)
avg_density_under_50 <- mean(ductdata[Age < 50]$Density,  na.rm = TRUE)
fold_change <- avg_density_over_50 / avg_density_under_50
cat("DUCT Average fold change:", fold_change, "\n")
cat("DUCT Log2 average fold change:", log2(fold_change), "\n")

lobuledata <- community_merged_data[TissueStructure == "lobule"]
lobuleagecommunityplot <- ggplot(lobuledata, aes(x = Age, y = Density, color = TissueStructure)) +
  geom_point(size = 0.1, alpha = 0.45) +
  geom_smooth(aes(group = TissueStructure), method = "loess", se = FALSE, fullrange = TRUE) +
  scale_color_manual(values = c("duct" = "#4169E1", "lobule" = "#DC143C")) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 10), expand = c(0.001, 0.01), breaks = c(0, 10), labels = c(0, 10)) +
  scale_x_continuous(expand = c(0.0001, 0.01), limits = c(15, 86), breaks = c(20, 50, 80)) +
  labs(x = NULL, y = TeX("density (communities/mm$^2$)")) +
  theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = unit(c(10, 0, 0, 0), units = "pt"),
    legend.position = "none"
  ) +
  annotate("text", x = 16, y = 9.5,
           label = TeX(paste0("$\\rho$ = ", format(cor_results_structures[["lobule"]]$rho, digits = 2))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 8.75,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(cor_results_structures[["lobule"]]$p.value)))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 8,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(lobuledata))),
           size = 5.5, color = "black", hjust = 0)

avg_density_over_50  <- mean(lobuledata[Age >= 50]$Density, na.rm = TRUE)
avg_density_under_50 <- mean(lobuledata[Age < 50]$Density,  na.rm = TRUE)
fold_change <- avg_density_over_50 / avg_density_under_50
cat("LOBULE Average fold change:", fold_change, "\n")
cat("LOBULE Log2 average fold change:", log2(fold_change), "\n")

ductlobulestructurescombined <- ductagecommunityplot + lobuleagecommunityplot
ductlobulex <- ggplot() + labs(x = "age") + theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        axis.title.x = element_text(color = "black", size = 16)) +
  guides(x = "none", y = "none")
ductlobulestructurescombined <- ductlobulestructurescombined / ductlobulex + plot_layout(heights = c(1000, 1))

cor_results <- list()

ductdatafiltered <- ductdata[MeanCommunitySize != 0]
cor_test <- cor.test(ductdatafiltered$Age, ductdatafiltered$MeanCommunitySize, method = "spearman")
cor_results[["duct"]] <- list(rho = cor_test$estimate, p.value = cor_test$p.value)

ductagecellsplot <- ggplot(ductdatafiltered, aes(x = Age, y = MeanCommunitySize, color = TissueStructure)) +
  geom_point(size = 0.1, alpha = 0.45) +
  geom_smooth(aes(group = TissueStructure), method = "loess", se = FALSE) +
  scale_color_manual(values = c("duct" = "#4169E1", "lobule" = "#DC143C")) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 200), expand = c(0.001, 0.01), breaks = c(0, 200), labels = c(0, 200)) +
  scale_x_continuous(expand = c(0.0001, 0.01), limits = c(15, 86), breaks = c(20, 50, 80)) +
  labs(x = NULL, y = TeX("mean cells per structure")) +
  theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_text(color = "black", size = 16),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_text(color = "black", size = 16),
    plot.margin = unit(c(10, 5, 0, 0), units = "pt"),
    legend.position = "none"
  ) +
  annotate("text", x = 16, y = 190,
           label = TeX(paste0("$\\rho$ = ", format(cor_results[["duct"]]$rho, digits = 2))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 175,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(cor_results[["duct"]]$p.value)))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 160,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(ductdatafiltered))),
           size = 5.5, color = "black", hjust = 0)

lobuledatafiltered <- lobuledata[MeanCommunitySize != 0]
cor_test <- cor.test(lobuledatafiltered$Age, lobuledatafiltered$MeanCommunitySize, method = "spearman")
cor_results[["lobule"]] <- list(rho = cor_test$estimate, p.value = cor_test$p.value)

lobuleagecellsplot <- ggplot(lobuledatafiltered, aes(x = Age, y = MeanCommunitySize, color = TissueStructure)) +
  geom_point(size = 0.1, alpha = 0.45) +
  geom_smooth(aes(group = TissueStructure), method = "loess", se = FALSE) +
  scale_color_manual(values = c("duct" = "#4169E1", "lobule" = "#DC143C")) +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 200), expand = c(0.001, 0.01), breaks = c(0, 200), labels = c(0, 200)) +
  scale_x_continuous(expand = c(0.0001, 0.01), limits = c(15, 86), breaks = c(20, 50, 80)) +
  labs(x = NULL, y = TeX("mean cells per structure")) +
  theme(
    axis.text.x = element_text(color = "black", size = 16),
    axis.text.y = element_blank(),
    axis.title.x = element_text(color = "black", size = 16),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    plot.margin  = unit(c(10, 0, 0, 0), units = "pt"),
    legend.position = "none"
  ) +
  annotate("text", x = 16, y = 190,
           label = TeX(paste0("$\\rho$ = ", format(cor_results[["lobule"]]$rho, digits = 2))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 175,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(cor_results[["lobule"]]$p.value)))),
           size = 5.5, color = "black", hjust = 0) +
  annotate("text", x = 16, y = 160,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(lobuledatafiltered))),
           size = 5.5, color = "black", hjust = 0)

ductlobulecellscombined <- ductagecellsplot + lobuleagecellsplot
ductlobulecellsx <- ggplot() + labs(x = "age") + theme_classic() +
  theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"),
        axis.title.x = element_text(color = "black", size = 16)) +
  guides(x = "none", y = "none")
ductlobulecellscombined <- ductlobulecellscombined / ductlobulecellsx + plot_layout(heights = c(1000, 1))

ggsave(file.path(outdir, "structuresDensityAge.pdf"),
       ductlobulestructurescombined, width = 5, height = 4.5, units = "in")
ggsave(file.path(outdir, "cellsDensityAge.pdf"),
       ductlobulecellscombined, width = 5, height = 4.5, units = "in")
