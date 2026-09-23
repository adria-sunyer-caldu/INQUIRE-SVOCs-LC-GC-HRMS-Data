############################################################
## PCA of GC-HRMS chemical features, INQUIRE study.
## Run twice in this script: once restricted to confirmed compounds
## (ID level 1-2), and once on all detected GC-HRMS features
## (nontarget, unfiltered, ID level 1-5), for direct comparison
## against the equivalent LC-HRMS PCA scripts.
## Restricted to the main indoor/outdoor samplers only (the extra indoor
## samplers, codes with a _2/_3 suffix, are different indoor sampling
## locations, not replicates — excluded here for paired indoor/outdoor
## sample logic). Sample names are pseudonymized: {Country}_{Indoor|Outdoor}_{code}.
## Run separately from the LC-HRMS PCA: different instrument,
## different intensity scale, minimal chemical overlap between the
## two platforms, so the two PCAs are not merged into one.
############################################################

if(!require(data.table)) install.packages("data.table")
if(!require(FactoMineR)) install.packages("FactoMineR")
if(!require(factoextra)) install.packages("factoextra")
if(!require(RColorBrewer)) install.packages("RColorBrewer")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(readxl)) install.packages("readxl")
if(!require(dplyr)) install.packages("dplyr")

library(data.table)
library(FactoMineR)
library(factoextra)
library(RColorBrewer)
library(ggplot2)
library(readxl)
library(dplyr)

# ---- Toggle: restrict to confirmed compounds (ID level 1-2), or use ALL
# GC-HRMS features. The first half of this script (below) always runs with
# this toggle; the second half (further down) always runs on all features,
# regardless of this setting, for the all-features PCA. ----
FILTER_CONFIRMED_ONLY <- TRUE

# ---- Set your working/data directory here ----
# Expects two input files in this directory:
#   "GC_Level1_2_polished_FINAL.xlsx"    (confirmed compounds, ID level 1-2)
#   "GC_Level1_5_polished_FULL_final.xlsx" (all features, ID level 1-5)
out_dir <- "path/to/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

gc_file <- "GC_Level1_2_polished_FINAL.xlsx"

# ---- Load ----
cat("Loading GC data...\n")
gc_raw <- read_excel(gc_file)
cat("Raw dimensions:", dim(gc_raw), "\n")

meta_cols <- c("RT","RT2","Feature_ID","Database match","Formula","InChiKey","Smiles","CAS",
               "IUPAC name","Method","m/z","Ion","MS/MS","MW","Exact Mass","ID level",
               "Componant info","RT/RI LC-MS","RT/RI UoA","RI","LoD","LoQ","Uncertainty",
               "Quantification Method","Notes on Batch Effects","Concentration Type")
meta_cols <- intersect(meta_cols, names(gc_raw))

if (FILTER_CONFIRMED_ONLY) {
  gc_use <- gc_raw %>% filter(`ID level` %in% c(1, 2))
  cat("Restricted to confirmed compounds (ID level 1-2):", nrow(gc_use), "of", nrow(gc_raw), "\n")
  suffix <- "targets_and_annotations"
} else {
  gc_use <- gc_raw
  cat("Using ALL GC features (nontarget):", nrow(gc_use), "\n")
  suffix <- "all_features"
}

# ---- Sample columns: main indoor/outdoor samplers only ----
sample_cols_all <- setdiff(names(gc_raw), meta_cols)
sample_cols <- sample_cols_all[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", sample_cols_all)]
cat("Indoor/outdoor sample columns:", length(sample_cols), "\n")

# ---- Exclude NL_Indoor_AR: extreme PCA outlier ----
# Driven entirely by ~9 rare compounds (e.g. n-hexyl salicylate, heptadecane,
# chloroxylenol, several fragrance/personal-care esters) that are detected
# ONLY in this one sample and are zero in all other ~408 indoor/outdoor samples.
# Under FactoMineR::PCA's default unit-variance scaling, near-zero-variance
# features like these get massively inflated, which is why this single
# sample separated so far from every other point on PC2. Excluded here
# pending confirmation this isn't a lab/batch artifact; noted in the
# figure caption.
outlier_sample <- "NL_Indoor_AR"
if (outlier_sample %in% sample_cols) {
  cat("Excluding known PCA outlier sample:", outlier_sample, "\n")
  sample_cols <- setdiff(sample_cols, outlier_sample)
  cat("Indoor/outdoor sample columns after exclusion:", length(sample_cols), "\n")
}

gc_use[sample_cols] <- lapply(gc_use[sample_cols], as.numeric)

# ---- Build sample x feature matrix (transpose: samples as rows for PCA) ----
feature_mat <- as.matrix(gc_use[, sample_cols])
feature_mat[is.na(feature_mat)] <- 0
pca_data <- t(feature_mat)   # rows = samples, columns = features
rownames(pca_data) <- sample_cols

# ---- Metadata per sample ----
Country <- sub("_.*", "", sample_cols)
Deployment <- ifelse(grepl("_Indoor_", sample_cols), "Indoor", "Outdoor")

# ---- Log1p-transform to reduce the influence of a small number of very
# high-intensity features/samples, while safely handling zeros ----
cat("Applying log1p transformation...\n")
pca_data_log <- log1p(pca_data)

# ---- Drop zero-variance features (no signal in this subset, can occur after
# filtering to confirmed compounds only) ----
zero_var <- apply(pca_data_log, 2, function(x) sd(x) == 0)
if (any(zero_var)) {
  cat("Dropping", sum(zero_var), "zero-variance features before PCA\n")
  pca_data_log <- pca_data_log[, !zero_var]
}

# ---- Run PCA ----
cat("Running PCA...\n")
res_pca <- PCA(pca_data_log, graph = FALSE)

var_explained <- res_pca$eig[, 2]
pc1_var <- round(var_explained[1], 1)
pc2_var <- round(var_explained[2], 1)
cat("PC1 variance explained:", pc1_var, "%\n")
cat("PC2 variance explained:", pc2_var, "%\n")

# ---- Export: variance explained per PC (source data) ----
variance_table_gc <- data.frame(
  PC = paste0("PC", seq_along(var_explained)),
  Variance_explained_pct = round(var_explained, 2),
  Cumulative_pct = round(cumsum(var_explained), 2)
)
write.csv(variance_table_gc, paste0("PCA_VarianceExplained_GC_", suffix, ".csv"), row.names = FALSE)
cat("Saved variance-explained table -> PCA_VarianceExplained_GC_", suffix, ".csv\n", sep = "")

# ---- Per-sample PCA coordinates (PC1/PC2), with sample metadata attached ----
pca_coords <- as.data.frame(res_pca$ind$coord)
pca_coords$SampleID <- rownames(pca_coords)
pca_coords$Country <- Country
pca_coords$Deployment <- Deployment

# ---- Colors used consistently across all figures in the project
# (SL and SI both map to Slovenia's color: SL was an old country code used
# in some raw files, replaced by SI everywhere else in the project) ----
country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C", "SL" = "#FCA00C"
)
inout_colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")

# ---- Plot: PCA colored by country ----
p_country <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Country)) +
  geom_point(size = 3, alpha = 0.9) +
  ggtitle("PCA of GC-HRMS features - Colored by Country") +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = country_colors) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

# ---- Plot: PCA colored by indoor/outdoor deployment, with 95% ellipses ----
p_deploy <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Deployment)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Deployment), type = "norm", level = 0.95, linetype = 2) +
  ggtitle("PCA of GC-HRMS features - Colored by Deployment") +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = inout_colors) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

print(p_country)
print(p_deploy)

# ---- Save ----
ggsave(paste0("PCA_Country_GC_", suffix, ".png"), p_country, width = 8, height = 6, dpi = 300, bg = "white")
ggsave(paste0("PCA_IndoorOutdoor_GC_", suffix, ".png"), p_deploy, width = 8, height = 6, dpi = 300, bg = "white")

# ---- Diagnostic: flag any sample that falls inside the OTHER deployment
# group's 95% confidence ellipse (Mahalanobis distance to the opposite
# group's mean/covariance below the chi-square 95% threshold for 2
# dimensions). This does not change the plotted data — it is a QC check
# for samples that look "misplaced" relative to their assigned group. ----
indoor_points  <- subset(pca_coords, Deployment == "Indoor")[,  c("Dim.1", "Dim.2")]
outdoor_points <- subset(pca_coords, Deployment == "Outdoor")[, c("Dim.1", "Dim.2")]
indoor_mean <- colMeans(indoor_points); indoor_cov <- cov(indoor_points)
outdoor_mean <- colMeans(outdoor_points); outdoor_cov <- cov(outdoor_points)
threshold <- qchisq(0.95, df = 2)

pca_coords$Inside_Other_Ellipse <- FALSE
pca_coords$Inside_Other_Ellipse[pca_coords$Deployment == "Indoor"] <-
  mahalanobis(indoor_points, center = outdoor_mean, cov = outdoor_cov) < threshold
pca_coords$Inside_Other_Ellipse[pca_coords$Deployment == "Outdoor"] <-
  mahalanobis(outdoor_points, center = indoor_mean, cov = indoor_cov) < threshold

misplaced <- subset(pca_coords, Inside_Other_Ellipse == TRUE)
cat("\nSamples inside the opposite group's 95% confidence ellipse:", nrow(misplaced), "of", nrow(pca_coords), "\n")

# ---- Export: per-sample PCA coordinates (source data for the figure) ----
write.csv(pca_coords[, c("SampleID","Country","Deployment","Dim.1","Dim.2","Inside_Other_Ellipse")],
          paste0("PCA_coords_GC_", suffix, ".csv"), row.names = FALSE)

cat("\nDone. PCA plots and coordinate table saved in", out_dir, "\n")

############################################################
## GC-HRMS PCA — ALL FEATURES (nontarget, unfiltered, ID level 1-5),
## regardless of the FILTER_CONFIRMED_ONLY toggle set above. This is
## the true nontarget-vs-nontarget comparison against the LC-HRMS
## all-features PCA.
############################################################

setwd(out_dir)

cat("\nRunning GC PCA on ALL features (nontarget, unfiltered, ID level 1-5)...\n")
gc_full_file <- "GC_Level1_5_polished_FULL_final.xlsx"
gc_all <- read_excel(gc_full_file)
cat("Using ALL GC features (from", gc_full_file, "):", nrow(gc_all), "\n")

# Recompute sample columns against this file's own header — do not assume it
# matches the Level 1-2 file's column names exactly.
sample_cols_all_full <- setdiff(names(gc_all), meta_cols)
sample_cols_full <- sample_cols_all_full[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", sample_cols_all_full)]
cat("Indoor/outdoor sample columns (full file):", length(sample_cols_full), "\n")
if (!setequal(sample_cols_full, sample_cols)) {
  cat("WARNING: sample columns in the full Level 1-5 file differ from the Level 1-2 file's sample columns — check both files' headers before trusting this PCA.\n")
}
sample_cols <- sample_cols_full

# Recompute Country/Deployment to match this block's own sample_cols order
Country <- sub("_.*", "", sample_cols)
Deployment <- ifelse(grepl("_Indoor_", sample_cols), "Indoor", "Outdoor")

gc_all[sample_cols] <- lapply(gc_all[sample_cols], as.numeric)

all_feature_mat <- as.matrix(gc_all[, sample_cols])
all_feature_mat[is.na(all_feature_mat)] <- 0
all_pca_data <- t(all_feature_mat)
rownames(all_pca_data) <- sample_cols

all_pca_data_log <- log1p(all_pca_data)
all_zero_var <- apply(all_pca_data_log, 2, function(x) sd(x) == 0)
if (any(all_zero_var)) {
  cat("Dropping", sum(all_zero_var), "zero-variance features before PCA\n")
  all_pca_data_log <- all_pca_data_log[, !all_zero_var]
}

res_pca_all <- PCA(all_pca_data_log, graph = FALSE)
var_explained_all <- res_pca_all$eig[, 2]
pc1_var_all <- round(var_explained_all[1], 1)
pc2_var_all <- round(var_explained_all[2], 1)
cat("All-features PC1 variance explained:", pc1_var_all, "%\n")
cat("All-features PC2 variance explained:", pc2_var_all, "%\n")

# ---- Export: variance explained per PC (source data) ----
variance_table_gc_all <- data.frame(
  PC = paste0("PC", seq_along(var_explained_all)),
  Variance_explained_pct = round(var_explained_all, 2),
  Cumulative_pct = round(cumsum(var_explained_all), 2)
)
write.csv(variance_table_gc_all, "PCA_VarianceExplained_GC_AllFeatures.csv", row.names = FALSE)
cat("Saved variance-explained table -> PCA_VarianceExplained_GC_AllFeatures.csv\n")

pca_coords_all <- as.data.frame(res_pca_all$ind$coord)
pca_coords_all$SampleID <- rownames(pca_coords_all)
pca_coords_all$Country <- Country
pca_coords_all$Deployment <- Deployment

p_country_all <- ggplot(pca_coords_all, aes(x = Dim.1, y = Dim.2, color = Country)) +
  geom_point(size = 3, alpha = 0.9) +
  ggtitle("PCA of GC-HRMS nontarget features - Colored by Country") +
  xlab(paste0("PC1 (", pc1_var_all, "%)")) +
  ylab(paste0("PC2 (", pc2_var_all, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = country_colors) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

p_deploy_all <- ggplot(pca_coords_all, aes(x = Dim.1, y = Dim.2, color = Deployment)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Deployment), type = "norm", level = 0.95, linetype = 2) +
  ggtitle("PCA of GC-HRMS nontarget features - Colored by Deployment") +
  xlab(paste0("PC1 (", pc1_var_all, "%)")) +
  ylab(paste0("PC2 (", pc2_var_all, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = inout_colors) +
  theme(legend.position = "right", plot.title = element_text(face = "bold"))

print(p_country_all)
print(p_deploy_all)

ggsave("PCA_Country_GC_all_features.png", p_country_all, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("PCA_IndoorOutdoor_GC_all_features.png", p_deploy_all, width = 8, height = 6, dpi = 300, bg = "white")

write.csv(pca_coords_all[, c("SampleID","Country","Deployment","Dim.1","Dim.2")],
          "PCA_coords_GC_all_features.csv", row.names = FALSE)

cat("\nDone. All-features GC PCA plots and coordinates saved in", out_dir, "\n")
