library(here)
source(here("code", "header.R"))

outdir <- here("scratch", "supp", "Fig8")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cells    <- getCells()
ctx      <- getCellContext()
clinical <- getClinical()

phenotypes <- c(
  "isSecretoryChange",
  "isApocrineMetaplasia",
  "isFibrocysticChange",
  "isColumnarChange",
  "isHyperplasia"
)
phenotypes <- phenotypes[phenotypes %in% names(ctx)]
if (length(phenotypes) == 0) stop("No epithelial histology flags found in cellContext.")

dt <- merge(
  cells[, .(ImageID, CellID, isEpithelial)],
  ctx[, c("ImageID", "CellID", phenotypes), with = FALSE],
  by = c("ImageID","CellID"),
  all.x = TRUE
)
dt <- merge(dt, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)

epi <- dt[isEpithelial == TRUE]

epi_copy <- copy(epi)
epi_copy <- na.omit(epi_copy)
if (nrow(epi_copy) == 0) stop("No epithelial rows remaining after na.omit().")

epi_copy[, AgeGroup := ifelse(Age < 50, "<50", ">=50")]

phenotypes_numeric <- lapply(epi_copy[, ..phenotypes], function(x) as.numeric(as.factor(x)) - 1)
epi_copy[, (phenotypes) := phenotypes_numeric]

patient_presence <- epi_copy[, lapply(.SD, max), by = .(PatientID, AgeGroup), .SDcols = phenotypes]

long_pat <- melt(
  patient_presence,
  id.vars = c("PatientID","AgeGroup"),
  measure.vars = phenotypes,
  variable.name = "Phenotype",
  value.name = "Presence"
)

props  <- long_pat[, .(Proportion = mean(Presence, na.rm = TRUE)), by = .(AgeGroup, Phenotype)]
counts <- long_pat[, .(Count      = sum(Presence, na.rm = TRUE)), by = .(AgeGroup, Phenotype)]
props  <- merge(props, counts, by = c("AgeGroup","Phenotype"), all = TRUE)

age_group_counts <- epi_copy[, .(nPatients = uniqueN(PatientID)), by = AgeGroup]
age_group_counts <- age_group_counts[order(factor(AgeGroup, levels = c("<50", ">=50")))]

props[, Phenotype := factor(Phenotype, levels = c(
  "isSecretoryChange",
  "isApocrineMetaplasia",
  "isFibrocysticChange",
  "isColumnarChange",
  "isHyperplasia"
))]

label_map <- c(
  "isHyperplasia"        = "hyperplasia",
  "isColumnarChange"     = "columnar\ncell change",
  "isFibrocysticChange"  = "fibrocystic\nchange",
  "isApocrineMetaplasia" = "apocrine\nmetaplasia",
  "isSecretoryChange"    = "secretory\ncell change"
)

fill_cols <- c("<50" = "#CBCBCB", ">=50" = "#737373")

p <- ggplot(props, aes(x = Phenotype, y = Proportion, fill = AgeGroup)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = Count, group = AgeGroup),
    position = position_dodge(width = 0.7),
    vjust = -0.25, color = "black", size = 5
  ) +
  scale_fill_manual(
    values = fill_cols,
    labels = c(
      "<50"  = paste0("Age < 50\nn = ",  age_group_counts[AgeGroup == "<50",  nPatients]),
      ">=50" = paste0("Age \u2265 50\nn = ", age_group_counts[AgeGroup == ">=50", nPatients])
    )
  ) +
  theme_classic() +
  labs(x = NULL, y = "proportion of patients\nwithin age group") +
  scale_x_discrete(labels = label_map) +
  scale_y_continuous(expand = c(0, 0), breaks = c(0, 0.10), limits = c(0, 0.13)) +
  theme(
    axis.text.x   = element_text(angle = 0, hjust = 0.5, colour = "black", size = 12),
    axis.text.y   = element_text(color = "black", size = 12),
    axis.title.y  = element_text(color = "black", size = 14),
    legend.title  = element_blank(),
    legend.text   = element_blank(),
    legend.key    = element_blank()
  ) + 
  theme(legend.position = "none")

ggsave(
  filename = file.path(outdir, "histologicalChangeAgeBarplots.pdf"),
  plot = p, units = "in", width = 7.5, height = 3
)

