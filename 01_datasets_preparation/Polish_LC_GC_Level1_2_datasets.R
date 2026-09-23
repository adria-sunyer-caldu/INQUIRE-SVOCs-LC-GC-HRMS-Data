# =============================================================================
# Polish_LC_GC_Level1_2_datasets.R
#
# Purpose: produce a single cleaned, consistent Level 1-2 dataset per platform
# (LC-HRMS, GC-HRMS) to be used everywhere downstream.
#
# For each platform:
#   1. Restrict to Confidence/ID level 1 and 2 (2a counts as "2" for LC).
#   2. Remove compounds with ZERO detection across ALL samples (value >0 in
#      no sample at all) - these are annotated but never actually observed
#      (e.g. the two GC phthalates found on 2026-08-19: Bis(2-butoxyethyl)
#      phthalate, Di(2-methoxyethyl) phthalate).
#   3. Flag duplicated InChIKeys (same InChIKey appearing under >1 Feature ID)
#      so the user can decide manually whether to merge, keep both, or drop.
#      Does NOT automatically resolve duplicates - flag only.
#
# Inputs (exact filenames as provided):
#   LC: "Targets and annotations concentrations.xlsx" (sheet 1)
#       - sample columns matched by regex ^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?\s*$
#       - Confidence level column values: 1, "2a" (treated as level 2)
#   GC: "GC_Level1_5_polished_FULL_zenodo.xlsx" (Zenodo GC deposit)
#       - sample columns matched by regex ^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$
#       - "ID level" column values: 1,2,3,5
#
# Outputs (written to same directory as this script by default):
#   LC_Level1_2_polished.xlsx   - EXACT SAME COLUMN FORMAT as the original
#                                 input file (same columns, same order, same
#                                 names). Only ROWS are filtered: level 1/2
#                                 only, zero-detection-across-all-samples rows
#                                 removed. No columns added or changed.
#   GC_Level1_2_polished.xlsx   - same, for the GC file.
#   LC_removed_zero_detection.csv / GC_removed_zero_detection.csv
#         - the rows that were removed for zero detection, kept for the record
#           (original columns, minus the sample columns to keep it compact).
#   LC_duplicate_inchikeys.csv / GC_duplicate_inchikeys.csv
#         - ONLY the rows whose InChIKey is duplicated within the polished
#           file, with a Duplicate_Group_Size column added, for manual review.
#           This flag is NOT written into the main polished file, so the
#           polished file's format/columns are untouched.
# =============================================================================

library(readxl)
library(writexl)
library(dplyr)
library(stringr)

# ---- USER-EDITABLE PATHS ---------------------------------------------------
setwd("path/to/data")

lc_path  <- "Targets and annotations concentrations.xlsx"   # derived LC table, not deposited
gc_path  <- "GC_Level1_5_polished_FULL_zenodo.xlsx"
out_dir  <- "."
# -----------------------------------------------------------------------------

polish_platform <- function(df, sample_cols, inchikey_col, level_col,
                             levels_keep, name_col, platform_label) {

  cat(sprintf("\n=== %s ===\n", platform_label))
  cat(sprintf("Rows before level filter: %d\n", nrow(df)))
  original_cols <- names(df)  # preserve exact original column set/order

  # 1) Level filter (character comparison handles "2a" vs numeric 2/1 cleanly)
  df <- df %>%
    filter(as.character(.data[[level_col]]) %in% levels_keep)
  cat(sprintf("Rows after restricting to levels %s: %d\n",
              paste(levels_keep, collapse = "/"), nrow(df)))

  # 2) Zero-detection-across-all-samples (computed but NOT stored as a column
  #    in the returned data - used only to split rows in/out)
  sample_mat <- df %>% select(all_of(sample_cols)) %>% as.matrix()
  sample_mat[is.na(sample_mat)] <- 0
  row_max <- apply(sample_mat, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  row_max[!is.finite(row_max)] <- 0
  is_zero <- row_max <= 0

  n_zero <- sum(is_zero)
  cat(sprintf("Compounds with zero detection in ALL samples (removed): %d\n", n_zero))
  if (n_zero > 0) {
    cat("  ->", paste(df[[name_col]][is_zero], collapse = "; "), "\n")
  }

  removed_zero <- df[is_zero, original_cols, drop = FALSE]
  polished     <- df[!is_zero, original_cols, drop = FALSE]  # exact original columns/order

  # 3) Duplicate InChIKey check, computed on the FINAL polished set, reported
  #    separately - never written into the polished file itself.
  ik_clean <- trimws(polished[[inchikey_col]])
  dup_counts <- table(ik_clean[!is.na(ik_clean) & ik_clean != ""])
  dup_group_size <- ifelse(
    is.na(ik_clean) | ik_clean == "", 1L,
    as.integer(dup_counts[ik_clean])
  )
  is_dup <- dup_group_size > 1

  n_dup_groups <- sum(dup_counts > 1)
  n_dup_rows   <- sum(is_dup)
  cat(sprintf("Duplicate InChIKey groups in polished set: %d (%d total rows involved)\n",
              n_dup_groups, n_dup_rows))

  duplicates <- polished[is_dup, original_cols, drop = FALSE]
  if (nrow(duplicates) > 0) {
    duplicates$Duplicate_Group_Size <- dup_group_size[is_dup]
    duplicates <- duplicates[order(duplicates[[inchikey_col]]), ]
    for (k in names(dup_counts[dup_counts > 1])) {
      names_here <- polished[[name_col]][ik_clean == k]
      cat(sprintf("  InChIKey %s: %s\n", k, paste(names_here, collapse = " | ")))
    }
  }

  list(
    polished     = polished,      # exact original format, rows filtered only
    removed_zero = removed_zero,  # exact original format
    duplicates   = duplicates     # original format + 1 extra flag column, for review only
  )
}

# ---- LC ----------------------------------------------------------------
lc_raw <- read_excel(lc_path, sheet = 1)
lc_sample_cols <- grep("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?\\s*$", names(lc_raw), value = TRUE)
stopifnot(length(lc_sample_cols) > 0)

lc_res <- polish_platform(
  df            = lc_raw,
  sample_cols   = lc_sample_cols,
  inchikey_col  = "InChiKey",
  level_col     = "Confidence level",
  levels_keep   = c("1", "2", "2a"),
  name_col      = "Name",
  platform_label = "LC-HRMS"
)

# Main polished file: EXACT original column format, rows filtered only
write_xlsx(lc_res$polished, file.path(out_dir, "LC_Level1_2_polished.xlsx"))

# Side files for the record / manual review - not used downstream as "the" dataset
write.csv(lc_res$removed_zero %>% select(-all_of(lc_sample_cols)),
          file.path(out_dir, "LC_removed_zero_detection.csv"), row.names = FALSE)
if (nrow(lc_res$duplicates) > 0) {
  write.csv(lc_res$duplicates %>% select(-all_of(lc_sample_cols)),
            file.path(out_dir, "LC_duplicate_inchikeys.csv"), row.names = FALSE)
}

# ---- GC ----------------------------------------------------------------
gc_raw_all <- read_excel(gc_path)
gc_raw <- gc_raw_all

gc_sample_cols <- grep("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", names(gc_raw), value = TRUE)
stopifnot(length(gc_sample_cols) > 0)

# "IUPAC name" already exists in this file - used only for console reporting,
# never added/modified as a column.
gc_res <- polish_platform(
  df            = gc_raw,
  sample_cols   = gc_sample_cols,
  inchikey_col  = "InChiKey",
  level_col     = "ID level",
  levels_keep   = c("1", "2"),
  name_col      = "IUPAC name",
  platform_label = "GC-HRMS"
)

# Main polished file: EXACT original column format, rows filtered only
write_xlsx(gc_res$polished, file.path(out_dir, "GC_Level1_2_polished.xlsx"))

# Side files for the record / manual review - not used downstream as "the" dataset
write.csv(gc_res$removed_zero %>% select(-all_of(gc_sample_cols)),
          file.path(out_dir, "GC_removed_zero_detection.csv"), row.names = FALSE)
if (nrow(gc_res$duplicates) > 0) {
  write.csv(gc_res$duplicates %>% select(-all_of(gc_sample_cols)),
            file.path(out_dir, "GC_duplicate_inchikeys.csv"), row.names = FALSE)
}

cat("\nDone. Polished .xlsx files have the EXACT original column format (rows filtered only).\n")
cat("Duplicate InChIKeys are reported in separate *_duplicate_inchikeys.csv files -\n")
cat("nothing was merged or dropped for duplication in the polished files; review manually.\n")
cat("Written to:", normalizePath(out_dir), "\n")

# =============================================================================
# SECTION 2: FULL DATASET, LEVEL 1-5, NO LEVEL RESTRICTION
#
# Same cleaning as above (remove zero-detection-across-all-samples rows),
# but WITHOUT the level 1-2 restriction - keeps every annotation level
# present in the raw files (LC: 1/2/2a/3/4/5; GC: 1/2/3/4/5).
#
# Duplicate InChIKeys are still only FLAGGED (not removed) for every group,
# EXCEPT one specific, manually-confirmed case: the row named "Triglyme"
# shares its InChIKey with one other row in the dataset that is NOT an
# isomer (unlike every other duplicate group here, which genuinely are
# isomers and are correctly left alone). That one non-isomer duplicate
# partner is removed by InChIKey match (not by name, so this is robust to
# whatever that other compound happens to be called) - Triglyme's own row
# is always kept.
# =============================================================================

cat("\n\n============================================================\n")
cat("SECTION 2: FULL DATASET (LEVEL 1-5), ZERO-DETECTION REMOVED\n")
cat("============================================================\n")

drop_named_duplicate_of <- function(df, inchikey_col, name_col, keep_name,
                                     platform_label) {
  ik_target <- df[[inchikey_col]][
    !is.na(df[[name_col]]) & trimws(tolower(df[[name_col]])) == tolower(keep_name)
  ]
  ik_target <- unique(trimws(na.omit(ik_target)))

  if (length(ik_target) == 0) {
    cat(sprintf("  [%s] '%s' not found - skipping the special duplicate removal.\n",
                platform_label, keep_name))
    return(df)
  }
  if (length(ik_target) > 1) {
    cat(sprintf("  [%s] WARNING: multiple InChIKeys found under the name '%s' (%s) - skipping to avoid a wrong removal, check manually.\n",
                platform_label, keep_name, paste(ik_target, collapse = ", ")))
    return(df)
  }

  is_match   <- trimws(df[[inchikey_col]]) == ik_target & !is.na(df[[inchikey_col]])
  is_keep_row <- trimws(tolower(df[[name_col]])) == tolower(keep_name) & !is.na(df[[name_col]])
  to_drop <- is_match & !is_keep_row

  n_drop <- sum(to_drop)
  if (n_drop == 0) {
    cat(sprintf("  [%s] No other row shares %s's InChIKey (%s) - nothing to remove.\n",
                platform_label, keep_name, ik_target))
  } else {
    cat(sprintf("  [%s] Removing %d row(s) sharing %s's InChIKey (%s), keeping %s itself:\n",
                platform_label, n_drop, keep_name, ik_target, keep_name))
    cat("    ->", paste(df[[name_col]][to_drop], collapse = "; "), "\n")
  }

  df[!to_drop, , drop = FALSE]
}

polish_full_platform <- function(df, sample_cols, inchikey_col, name_col,
                                  keep_name, platform_label) {

  cat(sprintf("\n=== %s (full, level 1-5) ===\n", platform_label))
  original_cols <- names(df)
  cat(sprintf("Rows before any filtering: %d\n", nrow(df)))

  # 1) Special case: drop the non-isomer duplicate of `keep_name`, by InChIKey
  df <- drop_named_duplicate_of(df, inchikey_col, name_col, keep_name, platform_label)

  # 2) Zero-detection-across-all-samples removal (same logic as Section 1)
  sample_mat <- df %>% select(all_of(sample_cols)) %>% as.matrix()
  sample_mat[is.na(sample_mat)] <- 0
  row_max <- apply(sample_mat, 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  row_max[!is.finite(row_max)] <- 0
  is_zero <- row_max <= 0

  n_zero <- sum(is_zero)
  cat(sprintf("Compounds with zero detection in ALL samples (removed): %d\n", n_zero))
  if (n_zero > 0) {
    cat("  ->", paste(df[[name_col]][is_zero], collapse = "; "), "\n")
  }

  removed_zero <- df[is_zero, original_cols, drop = FALSE]
  polished     <- df[!is_zero, original_cols, drop = FALSE]

  cat(sprintf("Rows in final full (level 1-5) dataset: %d\n", nrow(polished)))

  # 3) Remaining duplicate InChIKeys - flagged only, exactly as in Section 1
  #    (these are the isomer groups the user wants left alone).
  ik_clean <- trimws(polished[[inchikey_col]])
  dup_counts <- table(ik_clean[!is.na(ik_clean) & ik_clean != ""])
  dup_group_size <- ifelse(
    is.na(ik_clean) | ik_clean == "", 1L,
    as.integer(dup_counts[ik_clean])
  )
  is_dup <- dup_group_size > 1
  n_dup_groups <- sum(dup_counts > 1)
  cat(sprintf("Remaining duplicate InChIKey groups (flagged only, left in place): %d\n",
              n_dup_groups))
  if (n_dup_groups > 0) {
    for (k in names(dup_counts[dup_counts > 1])) {
      names_here <- polished[[name_col]][ik_clean == k]
      cat(sprintf("  InChIKey %s: %s\n", k, paste(names_here, collapse = " | ")))
    }
  }
  duplicates <- polished[is_dup, original_cols, drop = FALSE]
  if (nrow(duplicates) > 0) duplicates$Duplicate_Group_Size <- dup_group_size[is_dup]

  list(polished = polished, removed_zero = removed_zero, duplicates = duplicates)
}

# ---- LC: full level 1-5 -----------------------------------------------
lc_full_res <- polish_full_platform(
  df           = lc_raw,                 # original raw read, before Section 1's level filter
  sample_cols  = lc_sample_cols,
  inchikey_col = "InChiKey",
  name_col     = "Name",
  keep_name    = "Triglyme",
  platform_label = "LC-HRMS"
)

write_xlsx(lc_full_res$polished, file.path(out_dir, "LC_Level1_5_polished_FULL.xlsx"))
write.csv(lc_full_res$removed_zero %>% select(-all_of(lc_sample_cols)),
          file.path(out_dir, "LC_Level1_5_removed_zero_detection.csv"), row.names = FALSE)
if (nrow(lc_full_res$duplicates) > 0) {
  write.csv(lc_full_res$duplicates %>% select(-all_of(lc_sample_cols)),
            file.path(out_dir, "LC_Level1_5_duplicate_inchikeys.csv"), row.names = FALSE)
}

# ---- GC: full level 1-5 -------------------------------------------------
gc_full_res <- polish_full_platform(
  df           = gc_raw,                 # original raw read, before Section 1's level filter
  sample_cols  = gc_sample_cols,
  inchikey_col = "InChiKey",
  name_col     = "IUPAC name",
  keep_name    = "Triglyme",
  platform_label = "GC-HRMS"
)

write_xlsx(gc_full_res$polished, file.path(out_dir, "GC_Level1_5_polished_FULL.xlsx"))
write.csv(gc_full_res$removed_zero %>% select(-all_of(gc_sample_cols)),
          file.path(out_dir, "GC_Level1_5_removed_zero_detection.csv"), row.names = FALSE)
if (nrow(gc_full_res$duplicates) > 0) {
  write.csv(gc_full_res$duplicates %>% select(-all_of(gc_sample_cols)),
            file.path(out_dir, "GC_Level1_5_duplicate_inchikeys.csv"), row.names = FALSE)
}

cat("\nDone. Full (Level 1-5) datasets written - zero-detection rows removed,\n")
cat("Triglyme's non-isomer InChIKey duplicate removed, all other duplicate\n")
cat("groups (isomers) left untouched and only flagged for the record.\n")
cat("Written to:", normalizePath(out_dir), "\n")

