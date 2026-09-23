############################################################
## PCA of LC-HRMS target and confidently annotated (Level 1-2)
## chemical features, INQUIRE study.
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
# Expected input file: "Targets and annotations_FINAL.csv", with one row per
# sample and columns: SampleID, Country, Deployment, followed by one column
# per chemical feature (concentration/intensity values).
data_dir <- "path/to/data"
setwd(data_dir)
file_path <- file.path(data_dir, "Targets and annotations_FINAL.csv")

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

# ---- Keep only numeric (chemical feature) columns ----
pca_data <- pca_data[, sapply(pca_data, is.numeric), with = FALSE]
cat("Numeric feature matrix dimensions:", dim(pca_data), "\n")
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
variance_table_lc_targets <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance_explained_pct = round(var_explained, 2),
  Cumulative_pct = round(cumsum(var_explained), 2)
)
write.csv(variance_table_lc_targets, "PCA_VarianceExplained_LC_TargetsAnnotations.csv", row.names = FALSE)
cat("Saved variance-explained table -> PCA_VarianceExplained_LC_TargetsAnnotations.csv\n")

# ---- Per-sample PCA coordinates (PC1/PC2), with sample metadata attached ----
pca_coords <- as.data.frame(res_pca$ind$coord)
pca_coords$SampleID  <- rownames(pca_coords)
pca_coords$Country    <- metadata$Country
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
  ggtitle("PCA of LC-HRMS target and annotated compounds - Colored by Country") +
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
  ggtitle("PCA of LC-HRMS target and annotated compounds - Colored by Deployment") +
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
ggsave("PCA_Country_targets_and_annotations.png",       p_country, width = 8, height = 6, dpi = 300, bg = "transparent")
ggsave("PCA_IndoorOutdoor_targets_and_annotations.png", p_deploy,  width = 8, height = 6, dpi = 300, bg = "transparent")

cat("PCA complete. Plots saved as PNGs in the working directory.\n")

# ---- Diagnostic: flag any sample that falls inside the OTHER deployment
# group's 95% confidence ellipse (Mahalanobis distance to the opposite
# group's mean/covariance below the chi-square 95% threshold for 2
# dimensions). This does not change the plotted data — it is a QC check
# for samples that look "misplaced" relative to their assigned group. ----
pca_coords$SampleID   <- as.character(metadata$SampleID)
pca_coords$Deployment <- as.character(metadata$Deployment)

indoor_points  <- subset(pca_coords, Deployment == "Indoor")[,  c("Dim.1", "Dim.2")]
outdoor_points <- subset(pca_coords, Deployment == "Outdoor")[, c("Dim.1", "Dim.2")]

indoor_mean  <- colMeans(indoor_points)
indoor_cov   <- cov(indoor_points)
outdoor_mean <- colMeans(outdoor_points)
outdoor_cov  <- cov(outdoor_points)

threshold <- qchisq(0.95, df = 2)  # 2D PCA

pca_coords$Inside_Other_Ellipse <- FALSE
pca_coords$Inside_Other_Ellipse[pca_coords$Deployment == "Indoor"] <-
  mahalanobis(indoor_points, center = outdoor_mean, cov = outdoor_cov) < threshold
pca_coords$Inside_Other_Ellipse[pca_coords$Deployment == "Outdoor"] <-
  mahalanobis(outdoor_points, center = indoor_mean, cov = indoor_cov) < threshold

misplaced_samples <- subset(
  pca_coords[, c("SampleID", "Deployment", "Dim.1", "Dim.2", "Inside_Other_Ellipse")],
  Inside_Other_Ellipse == TRUE
)
cat("Samples inside the opposite group's 95% confidence ellipse:\n")
print(misplaced_samples)

# ---- Export: per-sample PCA coordinates (source data for the figure) ----
write.csv(pca_coords[, c("SampleID", "Country", "Deployment", "Dim.1", "Dim.2", "Inside_Other_Ellipse")],
          "PCA_coords_LC_TargetsAnnotations.csv", row.names = FALSE)
cat("Saved per-sample PCA coordinates -> PCA_coords_LC_TargetsAnnotations.csv\n")
