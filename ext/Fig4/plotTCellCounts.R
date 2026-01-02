library(here)
source(here("code", "header.R"))

outDir <- here("scratch/ext/Fig4")
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)

# --- Public setup ---
cells <- getCells()
clinical <- getClinical()
cellsclin <- merge(cells, clinical[, .(ImageID, PatientID, Age)], by = "ImageID", all.x = TRUE)

cd4_t_cells <- cellsclin[ClusterID == 12]
cd8_t_cells <- cellsclin[ClusterID == 13]

cd4_counts <- data.table(
  Marker = c("PD-1+", "TCF1+", "GZMB+"),
  Count = c(sum(cd4_t_cells$`isPD-1Pos`, na.rm = TRUE),
            sum(cd4_t_cells$isTCF1Pos, na.rm = TRUE),
            sum(cd4_t_cells$isGZMBPos, na.rm = TRUE)),
  Total = nrow(cd4_t_cells)
)
cd4_counts[, Percentage := (Count / Total) * 100]
cd4_counts[, Type := "CD4+T"]

cd8_counts <- data.table(
  Marker = c("PD-1+", "TCF1+", "GZMB+"),
  Count = c(sum(cd8_t_cells$`isPD-1Pos`, na.rm = TRUE),
            sum(cd8_t_cells$isTCF1Pos, na.rm = TRUE),
            sum(cd8_t_cells$isGZMBPos, na.rm = TRUE)),
  Total = nrow(cd8_t_cells)
)
cd8_counts[, Percentage := (Count / Total) * 100]
cd8_counts[, Type := "CD8+T"]

combined_counts <- rbind(cd4_counts, cd8_counts)
combined_counts[, Type := factor(Type, levels = rev(unique(Type)))]

custom_colors <- c("CD4+T" = "#2EF068", "CD8+T" = "#E0FCCF")

Tcounts <- ggplot(combined_counts, aes(x = Marker, y = Percentage, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", size = 0.2) +
  geom_text(aes(label = scales::comma(Count)),
            position = position_dodge(width = 1), hjust = -0.1, size = 5) +
  labs(title = "", x = "", y = "", fill = "") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(color = "black", size = 15)) +
  scale_y_continuous(limits = c(0, 30),
                     expand = expansion(mult = c(0, 0.1)),
                     breaks = c(0, 15, 30)) +
  coord_flip(clip = "off") +
  scale_fill_manual(values = custom_colors)

ggsave(Tcounts, filename = file.path(outDir, "TCellCounts.pdf"), units = "in", width = 2.5, height = 3)
