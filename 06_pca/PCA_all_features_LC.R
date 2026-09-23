############################################################
## PCA of LC-HRMS nontarget chemical features (all detected
## features, not restricted to confirmed/annotated compounds),
## INQUIRE study.
## Produces: PCA scatter plots colored by country and by
## indoor/outdoor deployment, the variance-explained table, and
## the per-sample PCA coordinate table (source data).
############################################################

# ---- Load required packages ----
if(!require(data.table)) install.packages("data.table")
if(!require(FactoMineR)) install.packages("FactoMineR")
if(!require(factoextra)) install.packages("factoextra")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(ggplot2)) install.packages("ggplot2")

library(data.table)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(ggplot2)

# ---- Set your working/data directory here ----
# Expected input file: "PCA-like dataset_FINAL.csv", with one row per sample
# and columns: SampleID, Country, Deployment, followed by one column per
# nontarget chemical feature (concentration/intensity values).
data_dir <- "path/to/data"
setwd(data_dir)
file_path <- file.path(data_dir, "PCA-like dataset_FINAL.csv")

# ---- Load dataset ----
cat("Loading data... this may take a few seconds.\n")
df <- fread(file = file_path, nThread = parallel::detectCores(), showProgress = TRUE)
cat("File loaded. Dimensions:", dim(df), "\n")

# ---- Separate metadata (sample info) from the chemical feature matrix ----
if(!all(c("SampleID", "Country", "Deployment") %in% names(df))){
  stop("Columns 'SampleID', 'Country', and 'Deployment' must exist in the dataset.")
}

metadata <- df[, .(SampleID, Country, Deployment)]
pca_data <- df[, !c("SampleID", "Country", "Deployment"), with = FALSE]
rownames(pca_data) <- metadata$SampleID

# ---- Log1p-transform to reduce the influence of a small number of very
# high-intensity features/samples, while safely handling zeros ----
cat("Applying log1p transformation...\n")
pca_data_log <- log1p(as.matrix(pca_data))

# ---- Run PCA ----
cat("Running PCA...\n")
res_pca <- PCA(pca_data_log, graph = FALSE)

# ---- Variance explained per principal component ----
var_explained <- res_pca$eig[, 2]
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)
cat("PC1 variance explained:", pc1_var, "%\n")
cat("PC2 variance explained:", pc2_var, "%\n")

# ---- Export: variance explained per PC (source data) ----
variance_table_lc_all <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance_explained_pct = round(var_explained, 2),
  Cumulative_pct = round(cumsum(var_explained), 2)
)
write.csv(variance_table_lc_all, "PCA_VarianceExplained_LC_AllFeatures.csv", row.names = FALSE)
cat("Saved variance-explained table -> PCA_VarianceExplained_LC_AllFeatures.csv\n")

# ---- Per-sample PCA coordinates (PC1/PC2), with sample metadata attached ----
pca_coords <- as.data.frame(res_pca$ind$coord)
pca_coords$SampleID <- rownames(pca_coords)
pca_coords$Country <- metadata$Country
pca_coords$Deployment <- metadata$Deployment

# ---- Colors used consistently across all figures in the project ----
country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C"
)
inout_colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")

# ---- Plot: PCA colored by country ----
p_country <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Country)) +
  geom_point(size = 3, alpha = 0.9) +
  ggtitle("PCA of LC-HRMS nontarget features - Colored by Country") +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = country_colors) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

# ---- Plot: PCA colored by indoor/outdoor deployment, with 95% ellipses ----
p_deploy <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Deployment)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Deployment), type = "norm", level = 0.95, linetype = 2) +
  ggtitle("PCA of LC-HRMS nontarget features - Colored by Deployment") +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = inout_colors) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )

print(p_country)
print(p_deploy)

# ---- Save high-resolution PNGs ----
ggsave("PCA_Country_all_features.png",       p_country, width = 8, height = 6, dpi = 300, bg = "transparent")
ggsave("PCA_IndoorOutdoor_all_features.png", p_deploy,  width = 8, height = 6, dpi = 300, bg = "transparent")

cat("PCA complete. Plots saved as PNGs in the working directory.\n")

# ---- Export: per-sample PCA coordinates (source data for the figure) ----
write.csv(pca_coords[, c("SampleID", "Country", "Deployment", "Dim.1", "Dim.2")],
          "PCA_coords_LC_AllFeatures.csv", row.names = FALSE)
cat("Saved per-sample PCA coordinates -> PCA_coords_LC_AllFeatures.csv\n")
