library(here)
source(here("code", "header.R"))

# Additional dependency required below
library(spatstat)

outdir <- here("scratch", "main", "Fig5")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
clinical <- getClinical()
clinical
imgs <- clinical[, .(ImageID, ImageWidth, ImageHeight)]
imgs <- unique(imgs, by = "ImageID")

if (!("X" %in% names(cells))) setnames(cells, "CenterX", "X", skip_absent = TRUE)
if (!("Y" %in% names(cells))) setnames(cells, "CenterY", "Y", skip_absent = TRUE)

epi <- cells[isEpithelial == TRUE,  .(ImageID, x = X, y = Y, type = "epithelial")]
tme <- cells[isEpithelial == FALSE, .(ImageID, ClusterID, x = X, y = Y, type = "TME")]

simulationnumber <- 199

process_one_image <- function(image_id) {
  img_row <- imgs[ImageID == image_id]
  if (nrow(img_row) != 1) return(NULL)
  
  boundary <- owin(
    xrange = c(0, img_row$ImageWidth),
    yrange = c(0, img_row$ImageHeight)
  )
  
  epi_i <- epi[ImageID == image_id]
  tme_i <- tme[ImageID == image_id]
  
  if (nrow(epi_i) == 0 || nrow(tme_i) == 0) return(NULL)
  
  clusters <- sort(unique(tme_i$ClusterID))
  
  res_list <- lapply(clusters, function(cl) {
    tme_cl <- tme_i[ClusterID == cl]
    
    combined <- rbind(
      epi_i[, .(x, y, type)],
      tme_cl[, .(x, y, type)]
    )
    
    combined_ppp <- ppp(
      x = combined$x, y = combined$y,
      window = boundary,
      marks = factor(combined$type)
    )
    
    mark_vals <- marks(combined_ppp)
    
    if (!("epithelial" %in% mark_vals) || !("TME" %in% mark_vals)) {
      return(data.table(
        ImageID = image_id,
        ClusterID = cl,
        OutOfRangeAbove = NA,
        OutOfRangeBelow = NA,
        ErrorMessage = NA
      ))
    }
    
    env <- tryCatch(
      envelope(
        combined_ppp,
        fun = Kcross,
        i = "epithelial",
        j = "TME",
        nsim = simulationnumber,
        correction = "border",
        alternative = "two.sided",
        maxnerr = simulationnumber * 2
      ),
      error = function(e) e
    )
    
    if (inherits(env, "error")) {
      return(data.table(
        ImageID = image_id, ClusterID = cl,
        OutOfRangeAbove = NA, OutOfRangeBelow = NA,
        ErrorMessage = as.character(env$message)
      ))
    }
    
    if (any(is.na(env$obs), is.na(env$hi), is.na(env$lo))) {
      return(data.table(ImageID = image_id, ClusterID = cl, OutOfRangeAbove = NA, OutOfRangeBelow = NA, ErrorMessage = "NA values in envelope"))
    }
    
    out_above <- any(env$obs > env$hi)
    out_below <- any(env$obs < env$lo)
    
    data.table(ImageID = image_id, ClusterID = cl, OutOfRangeAbove = out_above, OutOfRangeBelow = out_below, ErrorMessage = NA)
  })
  
  rbindlist(res_list, fill = TRUE)
}

# -------------------------
# run mode: local multicore
# -------------------------
image_ids <- unique(imgs$ImageID)
cat("Running bivariate Ripley Kcross for", length(image_ids), "images\n")

n_cores <- max(1, parallel::detectCores() - 1)
cat("Using", n_cores, "cores\n")

results_list <- parallel::mclapply(image_ids, process_one_image, mc.cores = n_cores)
all_results <- rbindlist(results_list, fill = TRUE)

saveRDS(all_results, file.path(outdir, "bivariateRipleyKAllResultsDT.rds"))
cat("Saved:", file.path(outdir, "bivariateRipleyKAllResultsDT.rds"), "\n")
