# =============================================================================
# Field_Blank_Levels_Analysis.R
#
# Field blank levels per confirmed LC-HRMS compound.
#   - TARGET compounds (CC-quantified): field blank signal was already
#     subtracted from sample peak area upstream, prior to calibration-curve
#     quantification. Reported here for the record only.
#   - ANNOTATION compounds (semiquantified, not calibration-curve based):
#     blank/sample levels reported the same way, for the record. No
#     pass/fail decision or exclusion is made in this script — it reports
#     the numbers only.
# Feeds Supplementary Table S12 (Field Blank Levels).
#
# Compound type (target vs annotation) is read from the "Quantification"
# column of LC_Level1_2_polished_FINAL.xlsx ("Target" vs "Semiquantification"),
# joined by Feature_ID — the same column used throughout this project's other
# LC scripts to distinguish the two.
#
# Input: Blanks_targets_annotations_POS_and_NEG_FINAL.txt — one row per
# compound; columns for real samples (COUNTRY_Hn_io, or pseudonymized\n# {Country}_{Indoor|Outdoor}_{code} in deposited files), field blanks
# (FieldBlank_COUNTRY_nn), AND calibration-curve points (Cal*_ESIPOS/NEG_*)
# mixed in with the metadata columns — Cal* columns are excluded below
# (not treated as sample, blank, or metadata; not part of this table).
# =============================================================================

library(readxl)
library(dplyr)

# ---- Set your working/data directory here ----
base_dir <- "path/to/data"
setwd(base_dir)

blanks_file   <- "Blanks_targets_annotations_POS_and_NEG_FINAL.txt"
polished_file <- "LC_Level1_2_polished_FINAL.xlsx"

# ---- Load compound type (target vs annotation) from the polished dataset ----
polished <- read_excel(polished_file)
names(polished) <- trimws(names(polished))
stopifnot("Quantification" %in% names(polished), "Feature ID" %in% names(polished))

compound_type <- polished %>%
  transmute(Feature_ID = `Feature ID`,
            Compound_type = ifelse(Quantification == "Target", "target", "annotation"))
cat("Compound type reference loaded:", nrow(compound_type), "compounds",
    "(", sum(compound_type$Compound_type == "target"), "targets,",
    sum(compound_type$Compound_type == "annotation"), "annotations )\n")

# ---- Load blanks data ----
blanks_raw <- read.delim(blanks_file, check.names = FALSE, stringsAsFactors = FALSE)

sample_cols <- names(blanks_raw)[grepl("^[A-Z]{2}_H[0-9]+_[12]$|^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?$", names(blanks_raw))]
blank_cols  <- names(blanks_raw)[grepl("^FieldBlank_", names(blanks_raw))]
# Calibration-curve columns (Cal0_ESIPOS_111_CZ, Cal0-005_ESINEG_4_EE, etc.) —
# not samples, not field blanks, not real metadata either. Excluded entirely,
# same as Rename_sample_columns_RT_Blanks.R drops Cal_/QC_/NORMSample_ columns.
cal_cols <- names(blanks_raw)[grepl("^Cal[0-9-]*_ESI(POS|NEG)_", names(blanks_raw))]
meta_cols <- setdiff(names(blanks_raw), c(sample_cols, blank_cols, cal_cols))

cat("Sample columns:", length(sample_cols), "| Field blank columns:", length(blank_cols),
    "| Calibration columns excluded:", length(cal_cols), "\n")
cat("Metadata columns kept as-is:", paste(meta_cols, collapse = ", "), "\n")

blanks_raw[c(sample_cols, blank_cols)] <- lapply(blanks_raw[c(sample_cols, blank_cols)],
                                                   function(x) suppressWarnings(as.numeric(x)))

# ---- Join compound type by Feature_ID ----
stopifnot("Feature_ID" %in% names(blanks_raw))
blanks_raw <- blanks_raw %>% left_join(compound_type, by = "Feature_ID")
n_unmatched <- sum(is.na(blanks_raw$Compound_type))
if (n_unmatched > 0) {
  cat("WARNING:", n_unmatched, "compounds in the blanks file did not match a Feature_ID",
      "in the polished dataset — check these before trusting the table.\n")
}

# ---- Per-compound blank/sample levels (numbers only, no pass/fail decision) ----
blank_levels <- blanks_raw %>%
  rowwise() %>%
  mutate(
    N_field_blanks   = sum(!is.na(c_across(all_of(blank_cols)))),
    Mean_blank_area  = { v <- c_across(all_of(blank_cols)); if (all(is.na(v))) 0 else round(mean(v, na.rm = TRUE), 2) },
    Max_sample_area  = round(max(c_across(all_of(sample_cols)), na.rm = TRUE), 2),
    Mean_sample_area = round(mean(c_across(all_of(sample_cols)), na.rm = TRUE), 2)
  ) %>%
  ungroup() %>%
  mutate(
    Fold_change_max_sample_vs_mean_blank = ifelse(Mean_blank_area > 0,
                                                    round(Max_sample_area / Mean_blank_area, 2), NA),
    Blank_pct_of_mean_sample = ifelse(Mean_sample_area > 0,
                                        round(Mean_blank_area / Mean_sample_area * 100, 2), NA)
  ) %>%
  select(all_of(meta_cols), Compound_type,
         N_field_blanks, Mean_blank_area, Max_sample_area, Mean_sample_area,
         Fold_change_max_sample_vs_mean_blank, Blank_pct_of_mean_sample)

cat("\nTable S12 (field blank levels): ", nrow(blank_levels), "compounds\n")
cat("Targets:", sum(blank_levels$Compound_type == "target"),
    "| Annotations:", sum(blank_levels$Compound_type == "annotation"), "\n")

write.csv(blank_levels, "Field_Blank_Levels_Table_S12.csv", row.names = FALSE)
cat("Saved: Field_Blank_Levels_Table_S12.csv\n")
