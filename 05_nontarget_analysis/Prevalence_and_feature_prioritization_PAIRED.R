# ==========================================================
# Full CEC Feature Prioritization Script
# UPDATED: every figure now uses the PAIRED indoor/outdoor ratio
# (per-household indoor/(indoor+outdoor), averaged across pairs),
# computed ONCE, early, with index-verified pairing — not the old
# intensity-summed ratio (sum(indoor)/(sum(indoor)+sum(outdoor))
# pooled across all samples), which has been removed everywhere.
# ==========================================================

library(data.table)
library(readxl)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(progress)
library(grid)

# ---------------------------
# Set your working/data directory here
# ---------------------------
data_dir <- "path/to/data"
setwd(data_dir)
data_file <- "LC_All_Features_FULL_final_zenodo.xlsx"   # Zenodo LC deposit
cat("📁 Working directory and file defined:\n", getwd(), "\n", data_file, "\n\n")

# ---------------------------
# 1) Load file
# ---------------------------
cat("⏳ Loading file (no header read) ...\n")
raw_data <- as.data.table(read_excel(data_file, sheet = "Full_dataset_INQUIRE_nontarget_",
                                     col_names = FALSE, .name_repair = "minimal"))
cat("✅ Raw data loaded. dim(raw_data):", paste(dim(raw_data), collapse = " x "), "\n\n")
flush.console()

cat("----- RAW DATA PREVIEW (rows 1..6, cols 1..40) -----\n")
print(as.data.frame(raw_data[1:6, 1:min(40, ncol(raw_data)), with = FALSE]))
cat("---------------------------------------------------\n\n")
flush.console()

# ---------------------------
# 2) meta preview and first sample column fixed
# ---------------------------
meta_preview <- raw_data[1:5, ]
cat("Meta_preview dims:", paste(dim(meta_preview), collapse = " x "), "\n")
cat("Head of meta_preview row1 (Class) cols 1..50:\n")
print(as.character(unlist(meta_preview[1, 1:min(50, ncol(meta_preview)), with = FALSE])))
cat("\n")
flush.console()

first_sample_col <- 33
cat("🧭 first_sample_col set to", first_sample_col, "(AG). ncol(raw_data) =", ncol(raw_data), "\n\n")
flush.console()

# ---------------------------
# 3) Extract class_info and sample_names
# ---------------------------
cat("Extracting class_info and sample_names...\n")
class_info <- as.character(unlist(meta_preview[1, first_sample_col:ncol(meta_preview), with = FALSE]))
class_info <- trimws(class_info)
sample_names <- as.character(unlist(meta_preview[4, first_sample_col:ncol(meta_preview), with = FALSE]))
sample_names <- trimws(sample_names)

cat("Lengths: length(class_info) =", length(class_info), "  length(sample_names) =", length(sample_names), "\n")
cat("Example class_info [1:40]:\n"); print(class_info[1:min(40, length(class_info))]); cat("\n")
cat("Example sample_names [1:40]:\n"); print(sample_names[1:min(40, length(sample_names))]); cat("\n\n")
flush.console()

cat("table(class_info):\n"); print(table(class_info, useNA = "ifany")); cat("\n")
flush.console()

# ---------------------------
# 4) Extract feature metadata (A:AF) and sample_matrix
# ---------------------------
cat("Extracting feature_meta (cols 1..(first_sample_col-1)) and sample_matrix...\n")
feature_meta <- raw_data[5:nrow(raw_data), 1:(first_sample_col-1)]
feature_meta_colnames <- as.character(unlist(meta_preview[4, 1:(first_sample_col-1), with = FALSE]))
feature_meta_colnames <- make.names(feature_meta_colnames, unique = TRUE)
colnames(feature_meta) <- feature_meta_colnames

sample_matrix <- raw_data[5:nrow(raw_data), first_sample_col:ncol(raw_data)]
colnames(sample_matrix) <- sample_names

cat("feature_meta dim:", paste(dim(feature_meta), collapse = " x "), "\n")
cat("sample_matrix dim:", paste(dim(sample_matrix), collapse = " x "), "\n\n")
flush.console()

cat("First few entries of sample_matrix[1, ] (as characters):\n")
print(as.character(sample_matrix[1, 1:min(20, ncol(sample_matrix)), with = FALSE]))
cat("\nConverting sample_matrix to numeric now (only replacing commas) ...\n")
flush.console()

sample_matrix[] <- lapply(sample_matrix, function(x) as.numeric(gsub(",", "", as.character(x))))
cat("Conversion done. Check NAs or non-numeric columns (sum of NA per column):\n")
print(colSums(is.na(sample_matrix))[1:min(20, ncol(sample_matrix))])
cat("\n")
flush.console()

# ---------------------------
# 5) Sanity checks: alignment and shapes
# ---------------------------
cat("=== SANITY CHECKS ===\n")
cat("ncol(sample_matrix):", ncol(sample_matrix), "length(class_info):", length(class_info), "length(sample_names):", length(sample_names), "\n")
ok_align <- length(class_info) == ncol(sample_matrix)
cat("class_info length == ncol(sample_matrix) ? ", ok_align, "\n")
if(!ok_align){
  cat(">>> MISMATCH: Please check sample columns and first_sample_col. Printing first 60 colnames(sample_matrix):\n")
  print(colnames(sample_matrix)[1:min(60, ncol(sample_matrix))])
  cat("Printing class_info[1:60]:\n"); print(class_info[1:min(60, length(class_info))])
  stop("Alignment mismatch between class_info and sample_matrix columns.")
}
cat("table of class_info (again):\n"); print(table(class_info))
cat("======================\n\n")
flush.console()

# ---------------------------
# 6) Remove features with zero across all samples (unchanged)
# ---------------------------
cat("Filtering zero-total features (no change in logic)...\n")
keep_features <- rowSums(sample_matrix > 0, na.rm = TRUE) > 0
cat("zero-filter removed:", sum(!keep_features), "features\n")
sample_matrix <- sample_matrix[keep_features, ]
feature_meta <- feature_meta[keep_features, ]
cat("Remaining features:", nrow(sample_matrix), "\n\n")
flush.console()

# ---------------------------
# 7) Prevalence filter >=5% (unchanged)
# ---------------------------
prevalence_threshold <- 0.05
prevalence <- rowSums(sample_matrix > 0, na.rm = TRUE) / ncol(sample_matrix)
keep_prev <- prevalence >= prevalence_threshold
cat("Prevalence filter (>=5%): kept:", sum(keep_prev), "features\n")
sample_matrix <- sample_matrix[keep_prev, ]
feature_meta <- feature_meta[keep_prev, ]
cat("After prevalence filter: nrow(sample_matrix) =", nrow(sample_matrix), "\n\n")
flush.console()

# Assign feature IDs now, before anything else touches row order —
# everything downstream (paired ratio, results_numeric, nap_table,
# Universal section) uses THIS row order and THESE IDs.
if(is.null(rownames(sample_matrix))){
  rownames(sample_matrix) <- paste0("Feature_", seq_len(nrow(sample_matrix)))
}
rownames(feature_meta) <- rownames(sample_matrix)
alignment_ids <- as.character(feature_meta[[1]])
cat("Feature IDs assigned. n =", length(alignment_ids), "\n\n")
flush.console()

# ==========================================================
# NEW / MOVED UP: PAIRED INDOOR/OUTDOOR RATIO
# This is THE ratio used everywhere below (Score, Volcano,
# Heatmap top30, Violin/Strip plots, Universal section, NAP
# export). Computed once, index-verified (match() on sample
# names, not relying on column order), instead of the old
# pooled-intensity ratio.
# ==========================================================
cat("\n⏳ Computing PAIRED indoor/outdoor ratio (per household, index-verified)...\n")

extract_country <- function(s) {
  s <- as.character(s)
  # pseudonymized sample names: {Country}_{Indoor|Outdoor}_{code}
  if (!grepl("^[A-Z]{2}_(Indoor|Outdoor)_", s)) return("Unknown")
  sub("_.*", "", s)
}
country_info <- sapply(sample_names, extract_country)

# Parse household code / country from the pseudonymized sample names
# ({Country}_{Indoor|Outdoor}_{code}; indoor and outdoor samples of the same
# household share the same code). rep_num is "1" for the Indoor-coded and "2"
# for the Outdoor-coded sample; Indoor/Outdoor itself still comes from
# class_info, and the loop below checks that the two agree. Extra indoor
# samplers (_2/_3 suffix) have no outdoor partner and are skipped.
parse_sample_v2 <- function(s) {
  m <- regmatches(s, regexec("^([A-Z]{2})_(Indoor|Outdoor)_([A-Z]{2}(_[23])?)$", s))[[1]]

  house   <- if(length(m) > 0) m[4] else NA
  country <- if(length(m) > 0) m[2] else NA
  rep_num <- if(length(m) > 0) ifelse(m[3] == "Indoor", "1", "2") else NA

  data.frame(house = house, country = country, rep_num = rep_num, stringsAsFactors = FALSE)
}

sample_meta <- do.call(rbind, lapply(sample_names, parse_sample_v2))
sample_meta$sample_name <- sample_names
sample_meta$deploy      <- class_info

cat("Sample metadata parsed. Head:\n")
print(head(sample_meta, 10))
cat("\nDeploy counts per country:\n")
print(table(sample_meta$country, sample_meta$deploy, useNA = "ifany"))
flush.console()

# Find valid pairs: same house+country, both replicates present,
# and opposite environments (one Indoor, one Outdoor)
valid_pairs <- list()
unique_combos <- unique(sample_meta[, c("house", "country")])
unique_combos <- unique_combos[!is.na(unique_combos$house) & !is.na(unique_combos$country), ]

for(i in seq_len(nrow(unique_combos))){
  h  <- unique_combos$house[i]
  co <- unique_combos$country[i]
  subset_hc <- sample_meta[sample_meta$house == h & sample_meta$country == co, ]
  rep1 <- subset_hc[subset_hc$rep_num == "1", ]
  rep2 <- subset_hc[subset_hc$rep_num == "2", ]

  if(nrow(rep1) == 0 || nrow(rep2) == 0) {
    cat("⚠️ Missing replicate for house:", h, "country:", co, "\n")
    next
  }
  if(rep1$deploy == rep2$deploy) {
    cat("⚠️ Same environment for both replicates — house:", h, "country:", co,
        "| rep1:", rep1$deploy, "| rep2:", rep2$deploy, "\n")
    next
  }

  if(rep1$deploy == "Indoor") {
    indoor_name  <- rep1$sample_name
    outdoor_name <- rep2$sample_name
  } else {
    indoor_name  <- rep2$sample_name
    outdoor_name <- rep1$sample_name
  }

  valid_pairs[[length(valid_pairs) + 1]] <- data.frame(
    house = h, country = co,
    indoor_name = indoor_name, outdoor_name = outdoor_name,
    stringsAsFactors = FALSE
  )
}

valid_pairs_df <- do.call(rbind, valid_pairs)
cat("✅ Valid pairs found:", nrow(valid_pairs_df), "\n\n")
print(head(valid_pairs_df, 10))
flush.console()

# Index-verified pairing (match on sample_names, NOT on column order)
n_feats <- nrow(sample_matrix)
n_pairs <- nrow(valid_pairs_df)

sm_mat <- as.matrix(sample_matrix)
colnames(sm_mat) <- sample_names

indoor_idx  <- match(valid_pairs_df$indoor_name,  sample_names)
outdoor_idx <- match(valid_pairs_df$outdoor_name, sample_names)
cat("Indoor index NA count:", sum(is.na(indoor_idx)), "\n")
cat("Outdoor index NA count:", sum(is.na(outdoor_idx)), "\n")
stopifnot(sum(is.na(indoor_idx)) == 0, sum(is.na(outdoor_idx)) == 0)

ratio_mat <- matrix(NA_real_, nrow = n_pairs, ncol = n_feats)
pb3 <- progress_bar$new(total = n_pairs, format = "  [:bar] :percent :eta")
for (p in seq_len(n_pairs)) {
  pb3$tick()
  ind_vals <- sm_mat[, indoor_idx[p]]
  out_vals <- sm_mat[, outdoor_idx[p]]
  total    <- ind_vals + out_vals
  ratio_mat[p, ] <- ifelse(total > 0, ind_vals / total, NA_real_)
}

paired_ratio <- colMeans(ratio_mat, na.rm = TRUE)   # length == n_feats, same row order as sample_matrix/feature_meta
cat("\nPaired ratio summary (across", n_pairs, "household pairs):\n")
print(summary(paired_ratio))
flush.console()

# ---------------------------
# 8) DEBUG: Print first feature vector and paired ratio (explicit)
# ---------------------------
cat("----- DEBUG CHECKPOINT: First feature vector paired ratio -----\n")
cat("Feature 1 paired ratio (mean across pairs):", paired_ratio[1], "\n")
cat("N pairs with data for feature 1:", sum(!is.na(ratio_mat[, 1])), "\n")
cat("----- END CHECKPOINT -----\n\n")
flush.console()

# ---------------------------
# 9) Main loop: prevalence (unchanged concept) + paired ratio (NEW) + p-values
# ---------------------------
cat("⏳ Computing prevalence and p-values for all features...\n")
n_features <- nrow(sample_matrix)
indoor_ratio <- paired_ratio          # <-- PAIRED ratio now, not intensity-summed
prevalence_indoor <- numeric(n_features)
prevalence_outdoor <- numeric(n_features)
pvals <- numeric(n_features)
pvals_paired <- numeric(n_features)

pb <- progress_bar$new(total = n_features, format = "  [:bar] :percent :eta")
N_DEBUG <- min(20, n_features)

for(i in seq_len(n_features)){
  pb$tick()
  vals <- as.numeric(sample_matrix[i, ])
  indoor_vals <- vals[class_info == "Indoor"]
  outdoor_vals <- vals[class_info == "Outdoor"]

  prevalence_indoor[i] <- sum(indoor_vals > 0, na.rm = TRUE) / length(indoor_vals)
  prevalence_outdoor[i] <- sum(outdoor_vals > 0, na.rm = TRUE) / length(outdoor_vals)

  # Pvalue: UNPAIRED Wilcoxon rank-sum test on raw intensities (Indoor
  # samples vs Outdoor samples, pooled). Kept for backward compatibility.
  if(length(unique(c(indoor_vals, outdoor_vals))) == 1){
    pvals[i] <- 1
  } else {
    pvals[i] <- tryCatch({
      wilcox.test(indoor_vals, outdoor_vals, paired = FALSE)$p.value
    }, error = function(e){
      warning(sprintf("Wilcox error at feature %d: %s", i, e$message))
      1
    })
  }

  # Pvalue_paired: paired Wilcoxon signed-rank test on the per-household
  # ratio (ratio_mat[, i]) against the null of 0.5 (indoor == outdoor).
  # This is the test consistent with the paired Indoor_Ratio effect-size
  # axis used in the volcano plot below.
  pair_vals <- ratio_mat[, i]
  pair_vals <- pair_vals[!is.na(pair_vals)]
  if(length(pair_vals) < 2 || length(unique(pair_vals)) == 1){
    pvals_paired[i] <- 1
  } else {
    pvals_paired[i] <- tryCatch({
      wilcox.test(pair_vals, mu = 0.5, paired = FALSE)$p.value
    }, error = function(e){
      warning(sprintf("Paired wilcox error at feature %d: %s", i, e$message))
      1
    })
  }

  if(i <= N_DEBUG){
    cat("=== Debug Feature", i, "===\n")
    cat("paired indoor_ratio:", indoor_ratio[i], "\n")
    cat("prev_in:", prevalence_indoor[i], " prev_out:", prevalence_outdoor[i], "\n")
    cat("wilcox p:", pvals[i], "\n")
    cat("=== end debug ===\n\n")
    flush.console()
  }
}

cat("\n✅ Completed loop. Summary of computed vectors:\n")
cat("NA counts: indoor_ratio (paired) NA =", sum(is.na(indoor_ratio)), " pvals NA =", sum(is.na(pvals)), "\n")
cat("Paired indoor ratio distribution (summary):\n"); print(summary(indoor_ratio)); cat("\n")
cat("Prevalence indoor summary:\n"); print(summary(prevalence_indoor)); cat("\n")
cat("Prevalence outdoor summary:\n"); print(summary(prevalence_outdoor)); cat("\n")
cat("Pval summary:\n"); print(summary(pvals)); cat("\n")
flush.console()

# ---------------------------
# NEW SECTION: Country × Indoor/Outdoor % Contribution (unchanged — mean
# intensity contribution per country×class group, a different concept
# from the ratio, not affected by this update)
# ---------------------------
cat("\n⏳ Computing %Contribution per Country × Indoor/Outdoor subset...\n")

group_labels <- paste(class_info, country_info, sep = "_")
unique_groups <- sort(unique(group_labels))
cat("Detected Country × Class groups:\n")
print(unique_groups)
cat("\n")

n_groups <- length(unique_groups)
n_features_current <- nrow(sample_matrix)

percent_contrib_matrix <- matrix(NA, nrow = n_features_current, ncol = n_groups)
colnames(percent_contrib_matrix) <- paste0("Pct_", unique_groups)
rownames(percent_contrib_matrix) <- rownames(sample_matrix)

for(i in seq_len(n_features_current)){
  vals <- as.numeric(sample_matrix[i, ])
  group_means <- numeric(n_groups)
  for(g in seq_along(unique_groups)){
    group_vals <- vals[group_labels == unique_groups[g]]
    group_means[g] <- mean(group_vals, na.rm = TRUE)
  }
  total_mean_sum <- sum(group_means, na.rm = TRUE)
  if(total_mean_sum > 0){
    percent_contrib_matrix[i, ] <- (group_means / total_mean_sum) * 100
  } else {
    percent_contrib_matrix[i, ] <- NA
  }
}
cat("Finished computing % contributions.\n\n")
flush.console()

# ---------------------------
# 10) results_numeric / results
# ---------------------------
results_numeric <- data.frame(
  Feature_ID = rownames(sample_matrix),
  Indoor_Ratio = indoor_ratio,     # <-- now the PAIRED ratio
  Prevalence_Indoor = prevalence_indoor,
  Prevalence_Outdoor = prevalence_outdoor,
  Pvalue = pvals,
  Pvalue_paired = pvals_paired,
  Score = abs(indoor_ratio - 0.5) * -log10(pvals),
  percent_contrib_matrix
)
cat("Check results_numeric structure:\n"); str(results_numeric); cat("\n")
flush.console()

results <- cbind(feature_meta, results_numeric)
cat("After cbind: dim(results) =", paste(dim(results), collapse = " x "), "\n\n")
flush.console()

# ---------------------------
# 11) Sort and show top 30
# ---------------------------
cat("Sorting by Score and showing top 30 (detailed)...\n")
results_sorted <- results[order(-results$Score), ]
show_cols <- intersect(c("Feature_ID", "Indoor_Ratio", "Prevalence_Indoor", "Prevalence_Outdoor", "Pvalue", "Score"), colnames(results_sorted))
cat("Columns to display:", paste(show_cols, collapse = ", "), "\n\n")
print(results_sorted[1:30, ..show_cols])
flush.console()

cat("\nShowing raw intensity rows for top 6 features (first 10 sample columns):\n")
top_ids <- results_sorted$Feature_ID[1:6]
print(as.data.frame(sample_matrix)[top_ids, 1:min(10, ncol(sample_matrix)), drop = FALSE])
flush.console()

# ==========================================================
# STEP: Volcano Plot — Paired Ratio (-1 to 1), Fixed Axis (-2 to 2)
# ==========================================================
library(ggplot2)

log_p_min <- 10
p_thresh <- 10^-log_p_min

indoor_ratio_scaled <- 2 * (results_numeric$Indoor_Ratio - 0.5)   # now paired

volcano_data <- data.frame(
  IndoorRatio = indoor_ratio_scaled,
  NegLogP = -log10(results_numeric$Pvalue_paired),
  Feature = results_numeric$Feature_ID
)
volcano_data <- subset(volcano_data, NegLogP >= log_p_min)
volcano_data$Color <- ifelse(volcano_data$IndoorRatio > 0, "Indoor", "Outdoor")

sig_indoor <- sum(indoor_ratio_scaled > 0 & -log10(results_numeric$Pvalue_paired) >= log_p_min)
sig_outdoor <- sum(indoor_ratio_scaled < 0 & -log10(results_numeric$Pvalue_paired) >= log_p_min)
cat("Significant features (p <", format(p_thresh, scientific = TRUE), "):\n")
cat("Indoor-abundant (paired ratio):", sig_indoor, "\n")
cat("Outdoor-abundant (paired ratio):", sig_outdoor, "\n")

if(nrow(volcano_data) == 0){
  warning("No points meet volcano plotting threshold. Volcano plot will be empty.")
  p_volcano <- ggplot() + ggtitle("Volcano Plot: no points")
} else {
  p_volcano <- ggplot(volcano_data, aes(x = IndoorRatio, y = NegLogP)) +
    geom_point(aes(color = Color), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")) +
    theme_minimal() +
    xlab("Paired Indoor Ratio (scaled -1 to 1)") +
    ylab("-log10(P-value)") +
    ggtitle(paste0("Nontarget feature-level indoor/outdoor significance Wilcoxon (paired ratio)\n",
                   "p < ", format(p_thresh, scientific = TRUE),
                   " | indoor n = ", sig_indoor, ", outdoor n = ", sig_outdoor))
  y_max <- max(volcano_data$NegLogP, na.rm = TRUE)
  if(is.finite(y_max) && y_max > log_p_min) p_volcano <- p_volcano + ylim(log_p_min, y_max)
  p_volcano <- p_volcano + xlim(-1, 1)
}
print(p_volcano)
ggsave("volcano_plot_thresholded_PAIRED.png", p_volcano, width = 9, height = 10, dpi = 600)
cat("✅ Saved: volcano_plot_thresholded_PAIRED.png\n")

# ==========================================================
# STEP: Heatmaps (Top30 Indoor/Outdoor by paired ratio)
# ==========================================================
library(ComplexHeatmap)
library(circlize)
library(grid)

cat("\n=== Heatmaps (paired ratio) ===\n")

valid_rows <- !is.na(results_sorted$Indoor_Ratio)
top_indoor_features <- results_sorted[valid_rows & results_sorted$Indoor_Ratio > 0.5, ][1:30, ]
top_outdoor_features <- results_sorted[valid_rows & results_sorted$Indoor_Ratio < 0.5, ][1:30, ]
top_indoor_features <- top_indoor_features[!is.na(top_indoor_features$Feature_ID), , drop = FALSE]
top_outdoor_features <- top_outdoor_features[!is.na(top_outdoor_features$Feature_ID), , drop = FALSE]

cat("Top Indoor features (paired ratio):", nrow(top_indoor_features), "\n")
cat("Top Outdoor features (paired ratio):", nrow(top_outdoor_features), "\n\n")

ordered_samples <- c(sample_names[class_info == "Indoor"], sample_names[class_info == "Outdoor"])
ordered_class <- class_info[match(ordered_samples, sample_names)]
ordered_country <- country_info[match(ordered_samples, sample_names)]
unique_countries <- sort(unique(ordered_country))
cat("Detected countries:", paste(unique_countries, collapse = ", "), "\n\n")

class_colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")
set.seed(42)
if(length(unique_countries) > 0){
  country_colors <- structure(circlize::rand_color(length(unique_countries)), names = unique_countries)
} else {
  country_colors <- NULL
}

ha <- HeatmapAnnotation(
  Class = ordered_class,
  Country = ordered_country,
  col = list(Class = class_colors, Country = country_colors),
  annotation_name_side = "left"
)

scale_rows_safe <- function(mat) {
  t(apply(mat, 1, function(x) {
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) rep(0, length(x))
    else (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
}

sm_df <- as.data.frame(sample_matrix)
rownames(sm_df) <- rownames(sample_matrix)

if(nrow(top_indoor_features) > 0){
  mat_indoor <- as.matrix(sm_df[top_indoor_features$Feature_ID, ordered_samples, drop = FALSE])
} else {
  mat_indoor <- matrix(nrow = 0, ncol = length(ordered_samples))
}
if(nrow(top_outdoor_features) > 0){
  mat_outdoor <- as.matrix(sm_df[top_outdoor_features$Feature_ID, ordered_samples, drop = FALSE])
} else {
  mat_outdoor <- matrix(nrow = 0, ncol = length(ordered_samples))
}

mat_indoor_scaled <- if(nrow(mat_indoor) >= 1 && ncol(mat_indoor) >= 1) scale_rows_safe(mat_indoor) else mat_indoor
mat_outdoor_scaled <- if(nrow(mat_outdoor) >= 1 && ncol(mat_outdoor) >= 1) scale_rows_safe(mat_outdoor) else mat_outdoor

if(nrow(mat_indoor_scaled) >= 1 && ncol(mat_indoor_scaled) >= 1){
  ht_indoor <- Heatmap(
    mat_indoor_scaled, name = "Scaled Intensity", top_annotation = ha,
    show_row_names = TRUE, show_column_names = FALSE,
    cluster_rows = TRUE, cluster_columns = FALSE, column_split = ordered_class,
    column_title = "Top 30 Indoor Features (paired ratio, z-scored)",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
} else { ht_indoor <- NULL; warning("No indoor heatmap data to plot.") }

if(nrow(mat_outdoor_scaled) >= 1 && ncol(mat_outdoor_scaled) >= 1){
  ht_outdoor <- Heatmap(
    mat_outdoor_scaled, name = "Scaled Intensity", top_annotation = ha,
    show_row_names = TRUE, show_column_names = FALSE,
    cluster_rows = TRUE, cluster_columns = FALSE, column_split = ordered_class,
    column_title = "Top 30 Outdoor Features (paired ratio, z-scored)",
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
} else { ht_outdoor <- NULL; warning("No outdoor heatmap data to plot.") }

if(dev.cur() == 1) dev.new()
if(!is.null(ht_indoor)){ grid.newpage(); draw(ht_indoor) }
if(!is.null(ht_outdoor)){ grid.newpage(); draw(ht_outdoor) }

if(!is.null(ht_indoor)){
  png("Heatmap_Top30_Indoor_PAIRED.png", width = 2000, height = 1500, res = 200)
  draw(ht_indoor); dev.off()
  cat("✅ Saved: Heatmap_Top30_Indoor_PAIRED.png\n")
}
if(!is.null(ht_outdoor)){
  png("Heatmap_Top30_Outdoor_PAIRED.png", width = 2000, height = 1500, res = 200)
  draw(ht_outdoor); dev.off()
  cat("✅ Saved: Heatmap_Top30_Outdoor_PAIRED.png\n")
}

if(nrow(mat_indoor_scaled) >= 1) write.csv(mat_indoor_scaled, "heatmap_matrix_top30_indoor_scaled_PAIRED.csv", row.names = TRUE)
if(nrow(mat_outdoor_scaled) >= 1) write.csv(mat_outdoor_scaled, "heatmap_matrix_top30_outdoor_scaled_PAIRED.csv", row.names = TRUE)
write.csv(data.frame(Sample = ordered_samples, Class = ordered_class, Country = ordered_country),
          "heatmap_sample_annotations.csv", row.names = FALSE)

# ---- Save key datasets ----
write.csv(results_numeric, "results_numeric_all_features_PAIRED.csv", row.names = FALSE)
write.csv(results, "results_all_features_with_metadata_PAIRED.csv", row.names = FALSE)
write.csv(top_indoor_features, "top30_indoor_features_PAIRED.csv", row.names = FALSE)
write.csv(top_outdoor_features, "top30_outdoor_features_PAIRED.csv", row.names = FALSE)

volcano_data_full <- data.frame(
  Feature = results_numeric$Feature_ID,
  IndoorRatio = 2*(results_numeric$Indoor_Ratio - 0.5),
  NegLogP = -log10(results_numeric$Pvalue),
  Color = ifelse(results_numeric$Indoor_Ratio > 0.5, "Red", "Blue")
)
write.csv(volcano_data_full, "volcano_data_full_PAIRED.csv", row.names = FALSE)

summary_counts <- data.frame(
  Threshold_Log10P = log_p_min, Pvalue_Threshold = p_thresh,
  Significant_Indoor = sig_indoor, Significant_Outdoor = sig_outdoor
)
write.csv(summary_counts, "volcano_significant_counts_PAIRED.csv", row.names = FALSE)
cat("✅ Heatmaps & key datasets saved (paired ratio).\n")

# ==========================================================
# STEP: Strip Plot — Paired Indoor Ratio, all features retained after
# the >=5% detection-frequency filter (this is the version used in the
# paper; three earlier variants of this plot — a prevalence-based
# violin, a paired-ratio strip restricted to the -log10(p) significant
# subset, and a prevalence-based strip for all >=5% DF features — were
# exploratory and are not used in any published figure, so have been
# removed here).
# ==========================================================
ratio_data <- data.frame(
  IndoorRatio = 2 * (results_numeric$Indoor_Ratio - 0.5),
  NegLogP     = -log10(results_numeric$Pvalue)
)
n_indoor_ratio  <- sum(ratio_data$IndoorRatio >= 0, na.rm = TRUE)
n_outdoor_ratio <- sum(ratio_data$IndoorRatio <  0, na.rm = TRUE)
cat("Indoor-enriched (paired ratio, all 5%DF):", n_indoor_ratio, "| Outdoor-enriched:", n_outdoor_ratio, "\n")

dens_r   <- density(ratio_data$IndoorRatio, adjust = 1, from = -1, to = 1, n = 2048, na.rm = TRUE)
y_vals_r <- dens_r$x
x_vals_r <- dens_r$y / max(dens_r$y) * 0.4
ratio_data$Color <- ifelse(ratio_data$IndoorRatio >= 0, "Indoor", "Outdoor")
ratio_data$half_width <- approx(y_vals_r, x_vals_r, xout = ratio_data$IndoorRatio, rule = 2)$y
set.seed(42)
ratio_data$jitter_x <- 1 + runif(nrow(ratio_data), -ratio_data$half_width * 0.85, ratio_data$half_width * 0.85)

p_strip_ratio_all <- ggplot() +
  geom_point(data = ratio_data, aes(x = jitter_x, y = IndoorRatio, color = Color), size = 1.2, alpha = 0.5) +
  scale_color_manual(values = c("Indoor" = "#8B0000", "Outdoor" = "#003366")) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "dashed") +
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.5, 1.6)) +
  labs(title = "Feature Distribution by Paired Indoor Ratio",
       subtitle = paste0("n = ", nrow(ratio_data), " features (≥5% DF), mean across ", n_pairs, " household pairs"),
       x = NULL, y = "Paired Indoor Ratio (scaled: -1 = fully outdoor, +1 = fully indoor)") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"), axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12), plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
        axis.line = element_line(color = "black"), legend.position = "none")

print(p_strip_ratio_all)
ggsave("strip_paired_indoor_ratio_all5pct.png", p_strip_ratio_all, dpi = 600, width = 5, height = 9)
cat("✅ Saved: strip_paired_indoor_ratio_all5pct.png\n")

# ---- Save paired ratio results (feature-level, id-matched) ----
paired_results <- data.frame(
  Feature_ID        = rownames(sample_matrix),
  Mean_Paired_Ratio = paired_ratio,
  Scaled_Ratio      = 2 * (paired_ratio - 0.5),
  N_pairs_detected  = colSums(!is.na(ratio_mat))
)
write.csv(paired_results, "paired_indoor_ratio_per_feature.csv", row.names = FALSE)
write.csv(valid_pairs_df, "valid_pairs_used.csv", row.names = FALSE)
cat("✅ Saved: paired_indoor_ratio_per_feature.csv\n")
cat("✅ Saved: valid_pairs_used.csv\n")

# ==========================================================
# STEP: NAP-ready paired ratio table (Alignment ID, Cluster Index, Mode)
# ==========================================================
cat("\n⏳ Generating NAP-ready output files...\n")
stopifnot(length(alignment_ids) == nrow(sample_matrix))

parse_alignment_id <- function(ids) {
  num  <- as.integer(sub("_(NEG|POS|ESIpos|ESIPOS).*", "", ids, ignore.case = TRUE))
  mode <- ifelse(grepl("NEG", ids, ignore.case = TRUE), "NEG", "POS")
  data.frame(Cluster_Index = num, Mode = mode, stringsAsFactors = FALSE)
}
parsed <- parse_alignment_id(alignment_ids)

nap_table <- data.frame(
  Alignment_ID       = alignment_ids,
  Cluster_Index      = parsed$Cluster_Index,
  Mode               = parsed$Mode,
  MZ                 = as.numeric(feature_meta[[3]]),
  RT                 = as.numeric(feature_meta[[2]]),
  Mean_Paired_Ratio  = paired_ratio,          # same value as results_numeric$Indoor_Ratio, row-aligned
  Scaled_Ratio       = 2 * (paired_ratio - 0.5),
  N_pairs_detected   = colSums(!is.na(ratio_mat)),
  stringsAsFactors   = FALSE
)

nap_table$Category <- cut(
  nap_table$Mean_Paired_Ratio,
  breaks = c(-Inf, 0.15, 0.35, 0.65, 0.85, Inf),
  labels = c("Strong_Outdoor", "Outdoor", "Neutral", "Indoor", "Strong_Indoor"),
  right  = TRUE
)

cat("\nCategory distribution:\n")
print(table(nap_table$Category, useNA = "ifany"))

write.csv(nap_table, "NAP_paired_ratio_all_features.csv", row.names = FALSE)
cat("✅ Saved: NAP_paired_ratio_all_features.csv\n")

strong_indoor  <- nap_table[!is.na(nap_table$Category) & nap_table$Category == "Strong_Indoor", ]
strong_outdoor <- nap_table[!is.na(nap_table$Category) & nap_table$Category == "Strong_Outdoor", ]
write.csv(strong_indoor,  "NAP_strong_indoor_features.csv",  row.names = FALSE)
write.csv(strong_outdoor, "NAP_strong_outdoor_features.csv", row.names = FALSE)
cat(paste("✅ Saved: NAP_strong_indoor_features.csv  (n =", nrow(strong_indoor), ")\n"))
cat(paste("✅ Saved: NAP_strong_outdoor_features.csv (n =", nrow(strong_outdoor), ")\n"))

cat("\nFirst 10 rows of NAP table:\n")
print(head(nap_table[order(nap_table$Mean_Paired_Ratio, decreasing = TRUE), ], 10))

# ============================================================
# UNIVERSAL INDOOR/OUTDOOR ENRICHMENT (paired ratio + per-country DF)
# Since paired_ratio is now computed once at the top and carried
# through in row order (results_numeric, results, nap_table all
# share the same row order — no separate merge needed).
# ============================================================
cat("\n⏳ Computing universal indoor/outdoor enrichment (paired ratio)...\n")

countries_unique <- sort(unique(country_info))
cat("Countries detected:", paste(countries_unique, collapse = ", "), "\n")

n_feat <- nrow(sample_matrix)
df_indoor_by_country  <- matrix(NA_real_, nrow = n_feat, ncol = length(countries_unique))
df_outdoor_by_country <- matrix(NA_real_, nrow = n_feat, ncol = length(countries_unique))
colnames(df_indoor_by_country)  <- paste0("DF_Indoor_",  countries_unique)
colnames(df_outdoor_by_country) <- paste0("DF_Outdoor_", countries_unique)

for (ci in seq_along(countries_unique)) {
  ctry <- countries_unique[ci]
  idx_in  <- which(class_info == "Indoor"  & country_info == ctry)
  idx_out <- which(class_info == "Outdoor" & country_info == ctry)
  if (length(idx_in)  > 0)
    df_indoor_by_country[, ci]  <- rowSums(sm_mat[, idx_in,  drop = FALSE] > 0, na.rm = TRUE) / length(idx_in)
  if (length(idx_out) > 0)
    df_outdoor_by_country[, ci] <- rowSums(sm_mat[, idx_out, drop = FALSE] > 0, na.rm = TRUE) / length(idx_out)
}
cat("Per-country DF matrices computed.\n")

DF_THRESHOLD    <- 0.50
RATIO_THRESHOLD <- 0.85   # PAIRED ratio threshold now

min_df_indoor  <- apply(df_indoor_by_country,  1, min, na.rm = TRUE)
min_df_outdoor <- apply(df_outdoor_by_country, 1, min, na.rm = TRUE)

is_universal_indoor <- (
  !is.na(results_numeric$Indoor_Ratio) &
    results_numeric$Indoor_Ratio > RATIO_THRESHOLD &
    min_df_indoor > DF_THRESHOLD
)
is_universal_outdoor <- (
  !is.na(results_numeric$Indoor_Ratio) &
    results_numeric$Indoor_Ratio < (1 - RATIO_THRESHOLD) &
    min_df_outdoor > DF_THRESHOLD
)

cat(paste("Universal indoor features  (paired ratio >", RATIO_THRESHOLD,
          "& DF >", DF_THRESHOLD, "in ALL countries):", sum(is_universal_indoor), "\n"))
cat(paste("Universal outdoor features (paired ratio <", 1 - RATIO_THRESHOLD,
          "& DF >", DF_THRESHOLD, "in ALL countries):", sum(is_universal_outdoor), "\n\n"))

df_indoor_dt  <- as.data.frame(df_indoor_by_country)
df_outdoor_dt <- as.data.frame(df_outdoor_by_country)

universal_base <- cbind(
  results[, c("Alignment.ID","Average.Mz","Average.Rt.min.",
              "Metabolite.name","Adduct.type","Formula",
              "Ontology","INCHIKEY","SMILES")],
  results_numeric[, c("Indoor_Ratio","Prevalence_Indoor","Prevalence_Outdoor","Pvalue","Score")],
  Min_DF_Indoor  = min_df_indoor,
  Min_DF_Outdoor = min_df_outdoor,
  df_indoor_dt,
  df_outdoor_dt
)

universal_indoor_features  <- universal_base[is_universal_indoor,  ]
universal_outdoor_features <- universal_base[is_universal_outdoor, ]
universal_indoor_features  <- universal_indoor_features[order(universal_indoor_features$Indoor_Ratio, decreasing = TRUE), ]
universal_outdoor_features <- universal_outdoor_features[order(universal_outdoor_features$Indoor_Ratio, decreasing = FALSE), ]

write.csv(universal_indoor_features,  "universal_indoor_features_PAIRED.csv",  row.names = FALSE)
write.csv(universal_outdoor_features, "universal_outdoor_features_PAIRED.csv", row.names = FALSE)
cat("Saved: universal_indoor_features_PAIRED.csv\n")
cat("Saved: universal_outdoor_features_PAIRED.csv\n\n")

# ── SENSITIVITY ANALYSIS (paired ratio) ───────────────────────
thresholds <- c(0.30, 0.40, 0.50, 0.60, 0.70, 0.80)
sensitivity <- data.frame(
  DF_Threshold = thresholds,
  N_Universal_Indoor  = sapply(thresholds, function(t)
    sum(!is.na(results_numeric$Indoor_Ratio) & results_numeric$Indoor_Ratio > RATIO_THRESHOLD & min_df_indoor > t, na.rm = TRUE)),
  N_Universal_Outdoor = sapply(thresholds, function(t)
    sum(!is.na(results_numeric$Indoor_Ratio) & results_numeric$Indoor_Ratio < (1 - RATIO_THRESHOLD) & min_df_outdoor > t, na.rm = TRUE))
)
cat("Sensitivity analysis (paired ratio):\n")
print(sensitivity)
write.csv(sensitivity, "universal_sensitivity_analysis_PAIRED.csv", row.names = FALSE)

library(tidyr)
sens_long <- pivot_longer(sensitivity, cols = c("N_Universal_Indoor","N_Universal_Outdoor"), names_to = "Type", values_to = "N")
sens_long$Type <- ifelse(sens_long$Type == "N_Universal_Indoor", "Indoor", "Outdoor")

p_sens <- ggplot(sens_long, aes(x = factor(DF_Threshold * 100), y = N, fill = Type)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = N), position = position_dodge(0.65), vjust = -0.4, size = 3.2, color = "grey20") +
  scale_fill_manual(values = c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")) +
  labs(title = "Number of Features: Sensitivity to Detection Frequency Threshold",
       subtitle = paste0("Paired indoor ratio > 0.7 | Paired outdoor ratio < -0.7"),
       x = "Minimum detection frequency per country (%)", y = "Number of features", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))
print(p_sens)
ggsave("Fig_Universal_Sensitivity_PAIRED.png", p_sens, width = 8, height = 5, dpi = 300)
cat("Saved: Fig_Universal_Sensitivity_PAIRED.png\n")

# NOTE: an earlier version of this script also plotted RT-vs-m/z chemical
# space for the "Universal" indoor/outdoor feature set defined above (paired
# ratio threshold + detection in >=50% of samples in every country). That
# plot is not used in any published figure — the chemical-space figure in
# the paper (Figure 7b) uses the "Strong Indoor/Outdoor" definition instead
# (paired ratio only, no cross-country requirement; see below) — so it has
# been removed here. The Universal feature lists and sensitivity analysis
# above are still exported and used (Supplementary Figure 5).
library(patchwork)

cat("\n=== UNIVERSAL ENRICHMENT COMPLETE (PAIRED RATIO) ===\n")
cat("Universal indoor features: ",  nrow(universal_indoor_features),  "\n")
cat("Universal outdoor features:", nrow(universal_outdoor_features), "\n")
cat("All figures/tables in this run use the PAIRED indoor/outdoor ratio.\n")


# ============================================================
# NEW SECTION: Chemical Space — PAIRED RATIO ONLY (no DF /
# cross-country filter, no "Universal" gate). Strong Indoor
# (paired ratio > 0.85) vs Strong Outdoor (paired ratio < 0.15),
# drawn from the full feature set (results_numeric$Indoor_Ratio),
# same as the KMD Strong Indoor/Outdoor definition.
# ============================================================
cat("\n⏳ Chemical Space — paired ratio only (no cross-country DF requirement)...\n")

RATIO_THRESHOLD_CS <- 0.85   # same threshold as KMD/Universal, just no DF gate

is_strong_indoor_cs  <- !is.na(results_numeric$Indoor_Ratio) & results_numeric$Indoor_Ratio > RATIO_THRESHOLD_CS
is_strong_outdoor_cs <- !is.na(results_numeric$Indoor_Ratio) & results_numeric$Indoor_Ratio < (1 - RATIO_THRESHOLD_CS)

cat("Strong Indoor  (paired ratio >", RATIO_THRESHOLD_CS, "):", sum(is_strong_indoor_cs), "\n")
cat("Strong Outdoor (paired ratio <", 1 - RATIO_THRESHOLD_CS, "):", sum(is_strong_outdoor_cs), "\n\n")

strong_indoor_cs  <- results[is_strong_indoor_cs,  c("Alignment.ID","Average.Mz","Average.Rt.min.")]
strong_outdoor_cs <- results[is_strong_outdoor_cs, c("Alignment.ID","Average.Mz","Average.Rt.min.")]
strong_indoor_cs$Indoor_Ratio  <- results_numeric$Indoor_Ratio[is_strong_indoor_cs]
strong_outdoor_cs$Indoor_Ratio <- results_numeric$Indoor_Ratio[is_strong_outdoor_cs]

strong_combined_cs <- rbind(
  data.frame(Alignment.ID = strong_indoor_cs$Alignment.ID,
             MZ = as.numeric(strong_indoor_cs$Average.Mz),
             RT = as.numeric(strong_indoor_cs$Average.Rt.min.),
             Paired_Ratio = strong_indoor_cs$Indoor_Ratio,
             Type = "Strong Indoor"),
  data.frame(Alignment.ID = strong_outdoor_cs$Alignment.ID,
             MZ = as.numeric(strong_outdoor_cs$Average.Mz),
             RT = as.numeric(strong_outdoor_cs$Average.Rt.min.),
             Paired_Ratio = strong_outdoor_cs$Indoor_Ratio,
             Type = "Strong Outdoor")
)
strong_combined_cs <- strong_combined_cs[!is.na(strong_combined_cs$MZ) & !is.na(strong_combined_cs$RT), ]
strong_combined_cs$Type <- factor(strong_combined_cs$Type, levels = c("Strong Outdoor","Strong Indoor"))

# ── EXPORT: Figure 7b source data — written immediately, before any
# plotting, so it is saved even if the downstream patchwork/ggsave
# calls fail (plotting is not required for this export) ──
write.csv(strong_combined_cs, "chemspace_PAIRED_StrongOnly.csv", row.names = FALSE)
cat(paste("Saved: chemspace_PAIRED_StrongOnly.csv (n =", nrow(strong_combined_cs), ")\n"))

col_vals_cs <- c("Strong Indoor" = "#B03A2E", "Strong Outdoor" = "#1F618D")

p_space_cs <- ggplot(strong_combined_cs, aes(x = RT, y = MZ, color = Type)) +
  geom_point(size = 2, stroke = 0, shape = 16, alpha = 0.5) +
  scale_color_manual(values = col_vals_cs, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "Chemical Space: Strong Indoor vs Strong Outdoor Features (Paired Ratio Only)",
       subtitle = paste0("Strong Indoor n=", sum(is_strong_indoor_cs),
                         " | Strong Outdoor n=", sum(is_strong_outdoor_cs),
                         " | No detection-frequency / cross-country filter applied"),
       x = "Retention time (min)", y = expression(italic(m/z)),
       caption = paste0("Strong Indoor: paired ratio > ", RATIO_THRESHOLD_CS,
                        " | Strong Outdoor: paired ratio < ", 1 - RATIO_THRESHOLD_CS)) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(size = 8.5, color = "grey40", hjust = 0),
        plot.caption = element_text(size = 8, color = "grey50"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3))

p_top_cs <- ggplot(strong_combined_cs, aes(x = RT, fill = Type, color = Type)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  scale_fill_manual(values = col_vals_cs) + scale_color_manual(values = col_vals_cs) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  theme_void() + theme(legend.position = "none")

p_right_cs <- ggplot(strong_combined_cs, aes(x = MZ, fill = Type, color = Type)) +
  geom_density(alpha = 0.3, linewidth = 0.5) +
  scale_fill_manual(values = col_vals_cs) + scale_color_manual(values = col_vals_cs) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  coord_flip() + theme_void() + theme(legend.position = "none")

p_space_cs_final <- (p_top_cs + plot_spacer() + p_space_cs + p_right_cs) +
  plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))

print(p_space_cs_final)

print(p_space_cs)
print(p_top_cs)
print(p_right_cs)

ggsave("ONLY ChemicalSpace_PAIRED.png", p_space_cs, width = 12, height = 6, dpi = 600)
ggsave("ONLY Density TOP.png", p_top_cs, width = 12, height = 1, dpi = 600)
ggsave("ONLY Density Right.png", p_right_cs, width = 1, height = 6, dpi = 600)

ggsave("Fig_ChemicalSpace_PAIRED_StrongOnly.png", p_space_cs_final, width = 10, height = 6, dpi = 300)
ggsave("Fig_ChemicalSpace_PAIRED_StrongOnly_main.png", p_space_cs, width = 10, height = 6, dpi = 300)

write.csv(strong_combined_cs, "chemspace_PAIRED_StrongOnly.csv", row.names = FALSE)

cat("\n=== CHEMICAL SPACE (PAIRED RATIO ONLY) COMPLETE ===\n")
cat("Strong Indoor: ",  sum(is_strong_indoor_cs),  "\n")
cat("Strong Outdoor:", sum(is_strong_outdoor_cs), "\n")
cat("No cross-country DF requirement applied — matches KMD's Strong Indoor/Outdoor definition.\n")
