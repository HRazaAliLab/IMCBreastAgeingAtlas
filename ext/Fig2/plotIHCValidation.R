library(here)
source(here::here("code/header.R"))

outdir <- here("scratch/ext/Fig2")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

file_paths <- c("ck5.csv", "er.csv", "foxa1.csv", "her2.csv", "ki67.csv", "pr.csv")
base_dir <- "scratch/FromEricLee/YVR Normal Breast"

# Load files into data.tables
ck5   <- fread(file.path(base_dir, file_paths[1]))
er    <- fread(file.path(base_dir, file_paths[2]))
foxa1 <- fread(file.path(base_dir, file_paths[3]))
her2  <- fread(file.path(base_dir, file_paths[4]))
ki67  <- fread(file.path(base_dir, file_paths[5]))
pr    <- fread(file.path(base_dir, file_paths[6]))

tables <- list(ck5, er, foxa1, her2, ki67, pr)
for (dt in tables) {
  setnames(dt, old = names(dt)[1:2], new = c("IHCpercentage", "IMCintensity"))
}

createScatterplot <- function(data, maxx1, maxy1, title) {
  correlation_test <- cor.test(data$IHCpercentage, data$IMCintensity, method = "spearman")
  
  plot <- ggplot(data, aes(x = IHCpercentage, y = IMCintensity)) +
    geom_point(size = 0.1, alpha = 0.2) +
    geom_smooth(method = "loess", se = FALSE, color = "dimgray") +
    coord_cartesian(clip = "off") +
    theme_classic() +
    scale_y_continuous(limits = c(0, maxy1), expand = c(0,0), breaks = c(0, maxy1), labels = function(x) ifelse(x == 0, "0", format(x))) +
    scale_x_continuous(expand = c(0.02,0), limits = c(0,maxx1), breaks = c(0, maxx1)) +
    annotate("text", x = min(data$IHCpercentage, na.rm = TRUE) * 1.05, y = maxy1 * 0.95,
             label = TeX(paste0("$\\rho$ = ", format(correlation_test$estimate, digits = 2))),
             size = 4.25, hjust = 0, color = "black") +
    annotate("text", x = min(data$IHCpercentage, na.rm = TRUE) * 1.05, y = maxy1 * 0.79,
             label = TeX(paste0("p-value = ", mkEnumPower(format_custom_pval(correlation_test$p.value)))),
             size = 4.25, hjust = 0, color = "black") +
    annotate("text", x = min(data$IHCpercentage, na.rm = TRUE) * 1.05, y = maxy1 * 0.61,
             label = TeX(paste0("$\\textit{n}$ = ", nrow(data))),
             size = 4.25, hjust = 0, color = "black") +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, color = "black", size = 12),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          plot.margin = unit(c(10,10,10,10), "pt")) +
    labs(subtitle = title)
  
  return(plot)
}

ck5plot <- createScatterplot(ck5, maxx1 = 100, maxy1 = 2.5, title = "CK5")
ki67plot <- createScatterplot(ki67, maxx1 = 10, maxy1 = 0.5, title = "Ki67")
foxa1plot <- createScatterplot(foxa1, maxx1 = 100, maxy1 = 0.7, title = "FOXA1")
erplot <- createScatterplot(er, maxx1 = 100, maxy1 = 1, title = "ER")
prplot <- createScatterplot(pr, maxx1 = 100, maxy1 = 0.7, title = "PR")
her2plot <- createScatterplot(her2, maxx1 = 10, maxy1 = 1, title = "HER2")

combinedplot <- ck5plot + ki67plot + foxa1plot + erplot + prplot + her2plot + plot_layout(ncol = 6)
ggsave(combinedplot, width = 14, height = 2, units = "in", filename = here(file.path(base_dir, "ScatterplotsIHCvsIMC.pdf")))

