# =============================================================================
# RT_Stability_Analysis.R
#
# Retention time stability across the full injection sequence, per confirmed
# LC-HRMS compound (internal standards + all confirmed targets/annotations).
# Feeds Supplementary Table S11 (RT Stability).
#
# For each compound: mean RT, SD, RSD%, min/max RT, and absolute drift range
# across every sample + field blank injection where that compound was
# detected (RT value > 0).
#
# Input: RT_targets_annotations_merged.txt (output of
# Rename_sample_columns_RT_Blanks.R) — one row per compound, one column per
# sample/field blank, values = retention time (min) in that injection.
#
# NOTE: this script identifies sample/field-blank columns by their renamed
# pattern (COUNTRY_Hn_io / FieldBlank_*) and treats every other column as
# metadata (kept as-is in the output, whatever it's actually called) —
# deliberately NOT hardcoding a specific metadata column name like "Name",
# since I have not verified the exact header text in your merged file.
# Check the metadata columns in the output against what you expect before
# trusting this for publication.
# =============================================================================

library(dplyr)

# ---- Set your working/data directory here ----
base_dir <- "path/to/data"
setwd(base_dir)

rt_file <- "RT_targets_annotations_POS_and_NEG_FINAL.txt"

# ---- Load merged RT data ----
rt_raw <- read.delim(rt_file, check.names = FALSE, stringsAsFactors = FALSE)

# Sample/field-blank RT columns, identified by the renamed pattern
# (COUNTRY_Hn_io, pseudonymized {Country}_{Indoor|Outdoor}_{code}, or\n# FieldBlank_COUNTRY_nn) produced by
# Rename_sample_columns_RT_Blanks.R. Everything else is treated as metadata.
sample_cols <- names(rt_raw)[grepl("^[A-Z]{2}_H[0-9]+_[12]$|^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?$", names(rt_raw)) |
                              grepl("^FieldBlank_", names(rt_raw))]
meta_cols <- setdiff(names(rt_raw), sample_cols)
cat("RT sample/field-blank columns identified:", length(sample_cols), "\n")
cat("Metadata columns kept as-is:", paste(meta_cols, collapse = ", "), "\n")

rt_raw[sample_cols] <- lapply(rt_raw[sample_cols], function(x) suppressWarnings(as.numeric(x)))

# ---- Per-compound RT stability across all injections where RT > 0 ----
compute_rt_stats <- function(row_vals) {
  vals <- row_vals[row_vals > 0 & !is.na(row_vals)]
  if (length(vals) < 2) {
    return(c(N_injections = length(vals), Mean_RT_min = NA, SD_RT_min = NA,
             RSD_pct = NA, Min_RT = NA, Max_RT = NA, Drift_range_min = NA))
  }
  m <- mean(vals); s <- sd(vals)
  c(N_injections = length(vals),
    Mean_RT_min = round(m, 4),
    SD_RT_min = round(s, 4),
    RSD_pct = if (m > 0) round(s / m * 100, 2) else NA,
    Min_RT = round(min(vals), 4),
    Max_RT = round(max(vals), 4),
    Drift_range_min = round(max(vals) - min(vals), 4))
}

rt_matrix <- as.matrix(rt_raw[, sample_cols])
stats_mat <- t(apply(rt_matrix, 1, compute_rt_stats))

rt_stability <- cbind(rt_raw[, meta_cols, drop = FALSE], as.data.frame(stats_mat))
rt_stability <- rt_stability[rt_stability$N_injections >= 2, ]

cat("Table S11 (RT stability): ", nrow(rt_stability), "compounds\n")
cat("RT RSD% distribution: median =", round(median(rt_stability$RSD_pct, na.rm = TRUE), 3),
    "| max =", round(max(rt_stability$RSD_pct, na.rm = TRUE), 3),
    "| min =", round(min(rt_stability$RSD_pct, na.rm = TRUE), 3), "\n")

write.csv(rt_stability, "RT_Stability_Table_S11.csv", row.names = FALSE)
cat("Saved: RT_Stability_Table_S11.csv\n")
