############################################################
## PCA of INQUIRE GC-HRMS SVOCs – Colored by Country & Deployment
## Adapted from the LC-HRMS PCA scripts for direct visual comparison.
## Sample codes are pseudonymized: {Country}_{Indoor|Outdoor}_{code},
## e.g. PT_Indoor_AK. Restricted to the main indoor/outdoor samplers
## (extra indoor samplers with a _2/_3 suffix are excluded, since they
## have no outdoor partner for paired-sample logic).
## NOT merged with the LC PCA: different instrument, different scale,
## minimal chemical overlap between platforms.
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

# ---- Toggle: restrict to confirmed compounds (ID level 1-2) or use ALL GC features ----
FILTER_CONFIRMED_ONLY <- FALSE   # set to TRUE to run on the confirmed/annotated subset instead

# ---- Paths ----
out_dir <- "path/to/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

gc_file <- "GC_Level1_5_polished_FULL_zenodo.xlsx"   # pseudonymized, deposit-ready full GC nontarget feature set

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
# ONLY in this one sample and are zero in all other ~408 indoor/outdoor
# samples. Under FactoMineR::PCA's default unit-variance scaling,
# near-zero-variance features like these get massively inflated, which is
# why this single sample separates so far from every other point. Excluded
# here; noted in the figure caption.
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
Country    <- sub("_.*", "", sample_cols)
Deployment <- ifelse(grepl("_Indoor_", sample_cols), "Indoor", "Outdoor")

# ---- Log-transform (handle zeros) ----
cat("Applying log1p transformation...\n")
pca_data_log <- log1p(pca_data)

# Drop any zero-variance columns (features with no signal in this subset — can occur after filtering)
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

# ---- Extract coordinates ----
pca_coords <- as.data.frame(res_pca$ind$coord)
pca_coords$SampleID <- rownames(pca_coords)
pca_coords$Country <- Country
pca_coords$Deployment <- Deployment

# ---- Custom colors (same palette as Figure 1 map) ----
country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C"
)
inout_colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")

# ---- PCA Plot: Country ----
p_country <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Country)) +
  geom_point(size = 3, alpha = 0.9) +
  ggtitle(paste0("GC-HRMS PCA by Country (", suffix, ")")) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = country_colors) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 13)) +
  coord_cartesian(clip = "off")

# ---- PCA Plot: Deployment ----
p_deploy <- ggplot(pca_coords, aes(x = Dim.1, y = Dim.2, color = Deployment)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = Deployment), type = "norm", level = 0.95, linetype = 2) +
  ggtitle(paste0("GC-HRMS PCA by Deployment (", suffix, ")")) +
  xlab(paste0("PC1 (", pc1_var, "%)")) +
  ylab(paste0("PC2 (", pc2_var, "%)")) +
  theme_minimal(base_size = 13) +
  scale_color_manual(values = inout_colors) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 13)) +
  coord_cartesian(clip = "off")

print(p_country)
print(p_deploy)

# ---- Save ----
ggsave(paste0("PCA_Country_GC_", suffix, ".png"), p_country, width = 10, height = 6, dpi = 600, bg = "white")
ggsave(paste0("PCA_IndoorOutdoor_GC_", suffix, ".png"), p_deploy, width = 10, height = 6, dpi = 600, bg = "white")

# ---- Mahalanobis diagnostic (same as the LC targets/annotations script) ----
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
write.csv(pca_coords[, c("SampleID","Country","Deployment","Dim.1","Dim.2","Inside_Other_Ellipse")],
          paste0("PCA_coords_GC_", suffix, ".csv"), row.names = FALSE)

cat("\nDONE. PCA plots and coordinate table saved in", out_dir, "\n")

