library(readxl)
library(ggplot2)

# ==========================================================
# GC equivalent of LC Step 20 / "Image 7" paired indoor ratio plot
# Full GC nontarget feature set (Zenodo GC deposit, ~1358 features)
# NOT the confirmed-only 106/159 subset used elsewhere in GC work
# ==========================================================

cat("⏳ Loading the GC nontarget dataset...\n")
# ---- Set your working/data directory here ----
data_dir <- "path/to/data"
setwd(data_dir)

gc_path <- "GC_Level1_5_polished_FULL_zenodo.xlsx"   # Zenodo GC deposit: 432 cols, 419 indoor/outdoor sample cols
raw <- read_excel(gc_path, col_names = TRUE)

feature_ids <- raw[["Feature_ID"]]
n_features_total <- length(feature_ids)
cat("Total feature rows loaded:", n_features_total, "\n")

# ---- Identify sample columns: {Country}_{Indoor|Outdoor}_{code} (main samplers only) ----
all_cols <- colnames(raw)
sample_cols <- all_cols[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", all_cols)]
cat("Sample columns matched (main indoor/outdoor samplers):", length(sample_cols), "\n")

# ---- Parse country / house / type from column names ----
parse_gc_sample <- function(s) {
  m <- regmatches(s, regexec("^([A-Z]{2})_(Indoor|Outdoor)_([A-Z]{2})$", s))[[1]]
  country <- m[2]
  house   <- m[4]   # household code, shared by the indoor and outdoor sample
  type    <- m[3]   # "Indoor" / "Outdoor"
  data.frame(sample_name = s, country = country, house = house, type = type,
             stringsAsFactors = FALSE)
}

sample_meta <- do.call(rbind, lapply(sample_cols, parse_gc_sample))

# ---- Find households with BOTH an indoor and an outdoor sample (extra indoor samplers excluded) ----
unique_combos <- unique(sample_meta[, c("house", "country")])

valid_pairs <- list()
for (i in seq_len(nrow(unique_combos))) {
  h  <- unique_combos$house[i]
  co <- unique_combos$country[i]
  subset_hc <- sample_meta[sample_meta$house == h & sample_meta$country == co, ]

  is1 <- subset_hc[subset_hc$type == "Indoor", ]
  os1 <- subset_hc[subset_hc$type == "Outdoor", ]

  if (nrow(is1) == 0 || nrow(os1) == 0) next  # missing pair, skip (matches known 1 missing house)

  valid_pairs[[length(valid_pairs) + 1]] <- data.frame(
    house = h, country = co,
    indoor_name  = is1$sample_name[1],
    outdoor_name = os1$sample_name[1],
    stringsAsFactors = FALSE
  )
}
valid_pairs_df <- do.call(rbind, valid_pairs)
cat("✅ Valid indoor/outdoor pairs found:", nrow(valid_pairs_df), "\n\n")

# ---- Build feature x sample matrix ----
sm_df <- as.data.frame(raw[, sample_cols])
sm_df <- as.data.frame(lapply(sm_df, function(x) suppressWarnings(as.numeric(x))))
rownames(sm_df) <- feature_ids

# ---- Compute paired ratio per feature across all valid pairs ----
# Same convention as LC: ratio = indoor / (indoor + outdoor), NA if total == 0
n_pairs <- nrow(valid_pairs_df)
paired_ratio_matrix <- matrix(NA_real_, nrow = n_pairs, ncol = n_features_total)
colnames(paired_ratio_matrix) <- feature_ids
rownames(paired_ratio_matrix) <- paste0(valid_pairs_df$house, "_", valid_pairs_df$country)

for (p in seq_len(n_pairs)) {
  ind <- as.numeric(sm_df[[valid_pairs_df$indoor_name[p]]])
  out <- as.numeric(sm_df[[valid_pairs_df$outdoor_name[p]]])

  total <- ind + out
  ratio <- ifelse(!is.na(total) & total > 0, ind / total, NA)

  paired_ratio_matrix[p, ] <- ratio
}

cat("✅ Paired ratio matrix computed. dim:", paste(dim(paired_ratio_matrix), collapse = " x "), "\n\n")

# ---- Mean paired ratio per feature ----
mean_paired_ratio <- colMeans(paired_ratio_matrix, na.rm = TRUE)

# ---- Detection frequency across pairs (non-NA fraction) = prevalence filter basis ----
detect_freq <- colSums(!is.na(paired_ratio_matrix)) / n_pairs

# ---- Apply >=5% prevalence filter (mirrors LC's 5% DF filtering step) ----
keep <- detect_freq >= 0.05 & !is.nan(mean_paired_ratio)
cat("Features passing >=5% DF filter:", sum(keep), "of", n_features_total, "\n\n")

# ---- Scale to -1 to +1 ----
paired_ratio_scaled <- 2 * (mean_paired_ratio - 0.5)

# ---- Build plot data ----
paired_data <- data.frame(
  PairedRatio = paired_ratio_scaled[keep],
  Feature_ID  = feature_ids[keep]
)
paired_data <- paired_data[!is.na(paired_data$PairedRatio), ]

n_indoor_paired  <- sum(paired_data$PairedRatio >= 0)
n_outdoor_paired <- sum(paired_data$PairedRatio <  0)
cat("Indoor-enriched (paired ratio):", n_indoor_paired, "\n")
cat("Outdoor-enriched (paired ratio):", n_outdoor_paired, "\n\n")

# ---- Density-constrained jitter ----
dens_p   <- density(paired_data$PairedRatio, adjust = 1.2, from = -1, to = 1, n = 2048)
y_vals_p <- dens_p$x
x_vals_p <- dens_p$y / max(dens_p$y) * 0.35

paired_data$Color      <- ifelse(paired_data$PairedRatio >= 0, "Indoor", "Outdoor")
paired_data$half_width <- approx(y_vals_p, x_vals_p, xout = paired_data$PairedRatio, rule = 2)$y

set.seed(42)
paired_data$jitter_x <- 1 + runif(nrow(paired_data),
                                  -paired_data$half_width * 0.85,
                                  paired_data$half_width * 0.85)

paired_data$group_dummy <- 1  # single-violin numeric x position (must match continuous x scale)

# ---- Build split-color violin polygon (red for y>=0 / indoor half, blue for y<0 / outdoor half) ----
# Reuses the same density object (dens_p/y_vals_p/x_vals_p) as the jitter width above,
# so the polygon outline and the dot jitter widths are the exact same shape.
pos_idx <- y_vals_p >= 0
neg_idx <- y_vals_p <= 0

violin_top <- data.frame(
  x = c(1 + x_vals_p[pos_idx], rev(1 - x_vals_p[pos_idx])),
  y = c(y_vals_p[pos_idx],     rev(y_vals_p[pos_idx]))
)
violin_bottom <- data.frame(
  x = c(1 + x_vals_p[neg_idx], rev(1 - x_vals_p[neg_idx])),
  y = c(y_vals_p[neg_idx],     rev(y_vals_p[neg_idx]))
)

# ---- Plot: split-color violin outline/fill behind jittered dots ----
p_strip_paired_gc <- ggplot() +

  geom_polygon(data = violin_top, aes(x = x, y = y),
               fill = "#8B0000", color = "#8B0000", alpha = 0.18, linewidth = 0.6) +
  geom_polygon(data = violin_bottom, aes(x = x, y = y),
               fill = "#003366", color = "#003366", alpha = 0.18, linewidth = 0.6) +

  geom_point(data = paired_data,
             aes(x = jitter_x, y = PairedRatio, color = Color),
             size = 1.1, alpha = 0.55) +

  scale_color_manual(values = c("Indoor" = "#8B0000", "Outdoor" = "#003366")) +

  geom_hline(yintercept = 0, color = "black", linewidth = 0.7, linetype = "dashed") +

  # Headroom above +-1 so points/violin tips at the extremes aren't cut off
  scale_y_continuous(limits = c(-1.08, 1.08), breaks = seq(-1, 1, 0.25), expand = c(0, 0)) +
  scale_x_continuous(limits = c(0.5, 1.5)) +

  labs(
    title    = "Feature Distribution by Mean Paired Indoor Ratio",
    subtitle = paste0("n = ", nrow(paired_data), " features (≥5% DF, mean across ",
                      n_pairs, " household pairs)"),
    x        = NULL,
    y        = "Mean Paired Indoor Ratio (−1 = fully outdoor, +1 = fully indoor)"
  ) +

  theme_minimal() +
  theme(
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    axis.ticks        = element_line(color = "black"),
    axis.text.x       = element_blank(),
    axis.text.y       = element_text(size = 12),
    plot.title        = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle     = element_text(size = 11, hjust = 0.5, color = "grey40"),
    axis.line         = element_line(color = "black"),
    legend.position   = "none"
  )

print(p_strip_paired_gc)
ggsave("GC_strip_paired_indoor_ratio.png", p_strip_paired_gc, dpi = 600, width = 5, height = 9)
cat("✅ Saved: GC_strip_paired_indoor_ratio.png\n")

# ---- Save results ----
gc_paired_results <- data.frame(
  Feature_ID        = feature_ids,
  Detection_Freq    = detect_freq,
  Mean_Paired_Ratio = mean_paired_ratio,
  Scaled_Ratio      = paired_ratio_scaled,
  Passed_5pct_DF    = keep
)
write.csv(gc_paired_results, "GC_paired_indoor_ratio_per_feature.csv", row.names = FALSE)
write.csv(valid_pairs_df, "GC_valid_pairs_used.csv", row.names = FALSE)
cat("✅ Saved: GC_paired_indoor_ratio_per_feature.csv\n")
cat("✅ Saved: GC_valid_pairs_used.csv\n")

# ==========================================================
# GC equivalent of LC volcano + sensitivity-threshold plots
# Reuses paired_ratio_matrix, mean_paired_ratio, detect_freq,
# valid_pairs_df (with $country) already computed above.
# ==========================================================

library(tidyr)

# ---------------------------------------------------------
# PART A: Paired-ratio volcano (feature-level significance)
# Paired Wilcoxon signed-rank test per feature, ratio vs null 0.5
# (matches the LC script's paired-test fix, not the earlier
# unpaired version)
# ---------------------------------------------------------

cat("\n⏳ Computing paired Wilcoxon signed-rank p-values per GC feature...\n")

n_feats_gc <- ncol(paired_ratio_matrix)
pvals_paired_gc <- numeric(n_feats_gc)

for (i in seq_len(n_feats_gc)) {
  pair_vals <- paired_ratio_matrix[, i]
  pair_vals <- pair_vals[!is.na(pair_vals)]
  if (length(pair_vals) < 2 || length(unique(pair_vals)) == 1) {
    pvals_paired_gc[i] <- 1
  } else {
    pvals_paired_gc[i] <- tryCatch({
      wilcox.test(pair_vals, mu = 0.5, paired = FALSE)$p.value
    }, error = function(e) 1)
  }
}
names(pvals_paired_gc) <- feature_ids

# ---- Export: add paired Wilcoxon p-values to the per-feature results table
# (gc_paired_results was already written above without this column — this
# adds it now that pvals_paired_gc has been computed) ----
gc_paired_results$Pvalue_paired <- pvals_paired_gc[gc_paired_results$Feature_ID]
write.csv(gc_paired_results, "GC_paired_indoor_ratio_per_feature_with_pvalues.csv", row.names = FALSE)
cat("✅ Saved: GC_paired_indoor_ratio_per_feature_with_pvalues.csv\n")

log_p_min_gc <- 2  # GC has far fewer features/pairs than LC (204 pairs, ~1358
                    # features vs LC's 205 pairs, 155,448 features), so p-values
                    # cannot reach LC's 1e-10 depth. Threshold set permissively;
                    # inspect the distribution and adjust before finalizing.
p_thresh_gc <- 10^-log_p_min_gc

volcano_data_gc <- data.frame(
  IndoorRatio = paired_ratio_scaled,
  NegLogP     = -log10(pvals_paired_gc),
  Feature     = feature_ids
)
volcano_data_gc <- subset(volcano_data_gc, NegLogP >= log_p_min_gc & keep)
volcano_data_gc$Color <- ifelse(volcano_data_gc$IndoorRatio > 0, "Indoor", "Outdoor")

sig_indoor_gc  <- sum(volcano_data_gc$IndoorRatio > 0)
sig_outdoor_gc <- sum(volcano_data_gc$IndoorRatio < 0)
cat("GC significant features (p <", format(p_thresh_gc, scientific = TRUE), "):\n")
cat("Indoor-abundant (paired ratio):", sig_indoor_gc, "\n")
cat("Outdoor-abundant (paired ratio):", sig_outdoor_gc, "\n")

if (nrow(volcano_data_gc) == 0) {
  warning("No GC points meet volcano plotting threshold (log_p_min_gc = ",
          log_p_min_gc, "). Lower the threshold and rerun.")
  p_volcano_gc <- ggplot() + ggtitle("GC Volcano Plot: no points at this threshold")
} else {
  p_volcano_gc <- ggplot(volcano_data_gc, aes(x = IndoorRatio, y = NegLogP)) +
    geom_point(aes(color = Color), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")) +
    theme_minimal() +
    xlab("Paired Indoor Ratio (scaled -1 to 1)") +
    ylab("-log10(P-value)") +
    ggtitle(paste0("Nontarget feature-level indoor/outdoor significance (GC, paired ratio)\n",
                    "p < ", format(p_thresh_gc, scientific = TRUE),
                    " | indoor n = ", sig_indoor_gc, ", outdoor n = ", sig_outdoor_gc)) +
    xlim(-1, 1)
}

print(p_volcano_gc)
ggsave("GC_volcano_plot_thresholded_PAIRED.png", p_volcano_gc, width = 8, height = 8, dpi = 300)
cat("✅ Saved: GC_volcano_plot_thresholded_PAIRED.png\n")

# ---------------------------------------------------------
# PART B: Sensitivity to detection-frequency threshold
# Universal Indoor/Outdoor = paired ratio > 0.7 / < -0.7,
# AND detected at >= threshold in EVERY country (mirrors LC's
# "applied in ALL 8 countries" logic)
# ---------------------------------------------------------

cat("\n⏳ Computing per-country detection frequency for GC sensitivity plot...\n")

countries_gc <- sort(unique(valid_pairs_df$country))
n_countries_gc <- length(countries_gc)

# Per-feature, per-country detection frequency (fraction of that country's
# pairs where the feature was detected, i.e. non-NA in paired_ratio_matrix)
detect_freq_by_country <- matrix(NA_real_, nrow = n_feats_gc, ncol = n_countries_gc,
                                  dimnames = list(feature_ids, countries_gc))
pair_country <- valid_pairs_df$country
names(pair_country) <- paste0(valid_pairs_df$house, "_", valid_pairs_df$country)
pair_country <- pair_country[rownames(paired_ratio_matrix)]

for (co in countries_gc) {
  rows_co <- which(pair_country == co)
  sub_mat <- paired_ratio_matrix[rows_co, , drop = FALSE]
  detect_freq_by_country[, co] <- colSums(!is.na(sub_mat)) / length(rows_co)
}

ratio_scaled_named <- paired_ratio_scaled
names(ratio_scaled_named) <- feature_ids

df_thresholds_gc <- seq(0.30, 0.80, by = 0.10)
sensitivity_gc <- data.frame(DF_Threshold = df_thresholds_gc,
                              N_Universal_Indoor  = NA_integer_,
                              N_Universal_Outdoor = NA_integer_)

for (t in seq_along(df_thresholds_gc)) {
  thresh <- df_thresholds_gc[t]
  meets_thresh_all_countries <- apply(detect_freq_by_country >= thresh, 1, all)

  is_universal_indoor  <- ratio_scaled_named > 0.7  & meets_thresh_all_countries
  is_universal_outdoor <- ratio_scaled_named < -0.7 & meets_thresh_all_countries

  sensitivity_gc$N_Universal_Indoor[t]  <- sum(is_universal_indoor, na.rm = TRUE)
  sensitivity_gc$N_Universal_Outdoor[t] <- sum(is_universal_outdoor, na.rm = TRUE)
}

cat("GC sensitivity table:\n")
print(sensitivity_gc)
write.csv(sensitivity_gc, "GC_universal_sensitivity_table.csv", row.names = FALSE)

sens_long_gc <- pivot_longer(sensitivity_gc,
                              cols = c("N_Universal_Indoor", "N_Universal_Outdoor"),
                              names_to = "Type", values_to = "N")
sens_long_gc$Type <- ifelse(sens_long_gc$Type == "N_Universal_Indoor", "Indoor", "Outdoor")

p_sens_gc <- ggplot(sens_long_gc, aes(x = factor(DF_Threshold * 100), y = N, fill = Type)) +
  geom_col(position = "dodge", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = N), position = position_dodge(0.65), vjust = -0.4, size = 3.2, color = "grey20") +
  scale_fill_manual(values = c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")) +
  labs(title = "Number of Features: Sensitivity to Detection Frequency Threshold (GC)",
       subtitle = "Paired indoor ratio > 0.7 | Paired outdoor ratio < -0.7",
       x = "Minimum detection frequency per country (%)", y = "Number of features", fill = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 12, hjust = 0),
        plot.subtitle = element_text(size = 9, color = "grey40"),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3))

print(p_sens_gc)
ggsave("GC_Fig_Universal_Sensitivity_PAIRED.png", p_sens_gc, width = 8, height = 5, dpi = 300)
cat("✅ Saved: GC_Fig_Universal_Sensitivity_PAIRED.png\n")
