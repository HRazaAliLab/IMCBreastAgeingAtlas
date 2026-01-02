library(here)
source(here("code", "header.R"))

# Additional dependencies required by the script
library(EBImage)
library(spatstat)

maskdir <- here("data/raw/images")
plotdir <- here("scratch/ext/Fig8")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
ctx      <- getCellContext()
clinical <- getClinical()

cells <- merge(
  cells,
  ctx[, .(ImageID, CellID, TissueStructure)],
  by = c("ImageID", "CellID"),
  all.x = TRUE
)
cellsclin <- merge(
  cells,
  clinical[, .(ImageID, PatientID, Age, ImageHeight)],
  by = "ImageID",
  all.x = TRUE
)


epi <- cellsclin[
  isEpithelial == TRUE &
  !is.na(TissueStructure) &
  !is.na(CenterX) & !is.na(CenterY) &
  !is.na(ImageHeight),
  .(ImageID, CellID, CenterX, CenterY, TissueStructure, ImageHeight)
]

read_myo_mask <- function(image_id) {
  f <- file.path(maskdir, paste0(image_id, "MyoepiMask.tiff"))
  if (!file.exists(f)) return(NULL)
  m <- EBImage::readImage(f)
  m <- EBImage::normalize(m)
  m <- m > 0
  m <- EBImage::flip(m)
  m
}

process_image <- function(image_id) {
  mask <- read_myo_mask(image_id)
  if (is.null(mask)) return(NULL)
  
  # thickness map (pixels) -> multiply by 2 as distmap reports half-thickness
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
  
  # nearest epithelial cell for each myoepithelial pixel
  myo_ppp <- ppp(
    thickness_dt$x, thickness_dt$y,
    c(0, max(thickness_dt$x)),
    c(0, max(thickness_dt$y))
  )
  
  epi_ppp <- ppp(
    epi_sub$CenterX, epi_sub$CenterY_legacy,
    c(0, max(thickness_dt$x)),
    c(0, max(thickness_dt$y))
  )
  
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

image_ids <- unique(clinical$ImageID)
cores <- max(1, parallel::detectCores() - 1)

img_res <- if (.Platform$OS.type == "unix") {
  parallel::mclapply(
    image_ids,
    function(id) tryCatch(process_image(id), error = function(e) NULL),
    mc.cores = cores
  )
} else {
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    library(data.table)
    library(EBImage)
    library(spatstat)
  })
  parallel::clusterExport(
    cl,
    varlist = c("epi","maskdir","read_myo_mask","process_image"),
    envir = environment()
  )
  parallel::parLapply(cl, image_ids, function(id) tryCatch(process_image(id), error = function(e) NULL))
}

image_dt <- rbindlist(img_res, fill = TRUE)
image_dt <- merge(image_dt, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)

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
  
  .(Age = Age[1], avg_thickness_duct = duct_wmean, avg_thickness_lobule = lobule_wmean)
}, by = PatientID]

both <- na.omit(patient_dt[, .(PatientID, avg_thickness_duct, avg_thickness_lobule)])
wil <- wilcox.test(both$avg_thickness_duct, both$avg_thickness_lobule, paired = TRUE)

melted <- melt(
  both,
  id.vars = "PatientID",
  variable.name = "Structure",
  value.name = "Avg_Thickness"
)

p_box <- ggplot(melted, aes(x = Structure, y = Avg_Thickness, fill = Structure)) +
  geom_boxplot(outlier.size = 0.1) +
  scale_fill_manual(values = c("avg_thickness_duct" = "#4169E1", "avg_thickness_lobule" = "#DC143C")) +
  labs(
    title = NULL, x = NULL, y = NULL,
    subtitle = TeX(paste0(
      "p=", mkEnumPower(ifelse(wil$p.value < 0.001,
                               formatC(wil$p.value, format = "e", digits = 0),
                               sprintf("%.3f", wil$p.value)))
    ))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    plot.subtitle = element_text(size = 13, hjust = 0.5, vjust = -3, color = "black"),
    legend.position = "none"
  ) +
  scale_y_continuous(breaks = c(3, 12), limits = c(2.5, 11.5), expand = expand_scale(mult = c(0, 0.1))) +
  scale_x_discrete(labels = c("duct", "lobule"))

out_box <- file.path(plotdir, "myoepiThicknessDuctvsLobule.pdf")
ggsave(plot = p_box, filename = out_box, width = 2, height = 1.5)

duct <- na.omit(patient_dt[, .(PatientID, Age, avg_thickness_duct)])
cor_duct <- cor.test(duct$Age, duct$avg_thickness_duct, method = "spearman")

miny1 <- 3
maxy1 <- 12

p_duct <- ggplot(duct, aes(x = Age, y = avg_thickness_duct)) +
  geom_point(size = 0.1, alpha = 0.3, color = "#4169E1") +
  geom_smooth(method = "loess", se = FALSE, color = "#4169E1") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(miny1, maxy1), expand = c(0, 0), breaks = c(miny1, maxy1)) +
  scale_x_continuous(breaks = c(20, 80), limits = c(15, 80), expand = c(0, 0)) +
  annotate("text", x = 17, y = ((maxy1 - miny1) * 0.925) + miny1,
           label = TeX(paste0("$\\rho$ = ", format(cor_duct$estimate, digits = 2))),
           size = 6, hjust = 0, color = "black") +
  annotate("text", x = 17, y = ((maxy1 - miny1) * 0.8) + miny1,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(cor_duct$p.value)))),
           size = 6, hjust = 0, color = "black") +
  annotate("text", x = 17, y = ((maxy1 - miny1) * 0.1) + miny1,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(duct))),
           size = 6, hjust = 0, color = "black") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(color = "black", size = 17),
    axis.text.y  = element_text(color = "black", size = 17),
    plot.margin  = unit(c(10, 10, 10, 10), "pt")
  )

out_duct <- file.path(plotdir, "myoepiThicknessAgeDuct.pdf")
ggsave(plot = p_duct, filename = out_duct, width = 2.6, height = 2)

lob <- na.omit(patient_dt[, .(PatientID, Age, avg_thickness_lobule)])
cor_lob <- cor.test(lob$Age, lob$avg_thickness_lobule, method = "spearman")

miny2 <- 2
maxy2 <- 9

p_lob <- ggplot(lob, aes(x = Age, y = avg_thickness_lobule)) +
  geom_point(size = 0.1, alpha = 0.3, color = "#DC143C") +
  geom_smooth(method = "loess", se = FALSE, color = "#DC143C") +
  coord_cartesian(clip = "off") +
  theme_classic() +
  scale_y_continuous(limits = c(miny2, maxy2), expand = c(0, 0), breaks = c(miny2, maxy2)) +
  scale_x_continuous(breaks = c(20, 80), limits = c(15, 80), expand = c(0, 0)) +
  annotate("text", x = 17, y = ((maxy2 - miny2) * 0.925) + miny2,
           label = TeX(paste0("$\\rho$ = ", format(cor_lob$estimate, digits = 2))),
           size = 6, hjust = 0, color = "black") +
  annotate("text", x = 17, y = ((maxy2 - miny2) * 0.8) + miny2,
           label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(cor_lob$p.value)))),
           size = 6, hjust = 0, color = "black") +
  annotate("text", x = 17, y = ((maxy2 - miny2) * 0.1) + miny2,
           label = TeX(paste0("$\\textit{n}$ = ", nrow(lob))),
           size = 6, hjust = 0, color = "black") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(color = "black", size = 17),
    axis.text.y  = element_text(color = "black", size = 17),
    plot.margin  = unit(c(10, 10, 10, 10), "pt")
  )

out_lob <- file.path(plotdir, "myoepiThicknessAgeLobule.pdf")
ggsave(plot = p_lob, filename = out_lob, width = 2.6, height = 2)

cat("Done.\nWrote:\n",
    out_box, "\n",
    out_duct, "\n",
    out_lob, "\n")
