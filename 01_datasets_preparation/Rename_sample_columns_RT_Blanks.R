# =============================================================================
# Rename_sample_columns_RT_Blanks.R
#
# For the filtered RT and Blanks target/annotation files:
#   - Real household samples: Sample_H{n}_{1|2}_ESI(POS|NEG)_{injnum}_{COUNTRY}
#       -> renamed to {COUNTRY}_H{n}_{1|2}, matching the polished LC dataset exactly
#   - Field blanks: FB_{n}_ESI(POS|NEG)_{injnum}_{COUNTRY}
#       -> renamed to FieldBlank_{COUNTRY}_{n, zero-padded}, keeping the same
#          n as already used (which follows injection order)
#   - Calibration curve (Cal_*), QC (QC_*), and normalized-sample duplicate
#     (NORMSample_*) columns are DROPPED entirely (not needed for this table)
#   - All other (metadata) columns are kept, untouched
#
# POS and NEG are merged into ONE output file each (one for RT, one for
# Blanks), with an added "Ionization" column (POS/NEG) and "Feature_ID"
# column (Alignment ID + polarity suffix, e.g. "67328_POS_SU") so rows stay
# uniquely traceable back to the polished dataset once combined.
#
# NOTE ON THE MERGE: the original multi-row header block (Class, Sample
# Type, Injection Order, Batch ID) is DROPPED in the merged output. Those
# rows describe each SAMPLE COLUMN, but a given household sample (e.g.
# CZ_H1_1) has a DIFFERENT injection order and batch number in the POS run
# vs the NEG run - once POS and NEG rows share the same column name, there
# is no single value left to put in that header cell. If you need
# batch/injection-order preserved, that requires a long-format table
# instead (one row per compound x sample, with batch/injection order as a
# per-row attribute) - let me know if you want that version instead.
# =============================================================================

base_dir <- "path/to/data"

classify_col <- function(col) {
  # Returns list(keep = TRUE/FALSE, new_name = character)
  m_sample <- regmatches(col, regexec(
    "^Sample_H([0-9]+)_([12])_ESI(POS|NEG)_[0-9]+_([A-Z]{2})$", col))[[1]]
  if (length(m_sample) == 5) {
    hh <- m_sample[2]; io <- m_sample[3]; country <- m_sample[5]
    return(list(keep = TRUE, new_name = sprintf("%s_H%s_%s", country, hh, io)))
  }

  m_fb <- regmatches(col, regexec(
    "^FB_([0-9]+)_ESI(POS|NEG)_[0-9]+_([A-Z]{2})$", col))[[1]]
  if (length(m_fb) == 4) {
    n <- as.integer(m_fb[2]); country <- m_fb[4]
    return(list(keep = TRUE, new_name = sprintf("FieldBlank_%s_%02d", country, n)))
  }

  if (grepl("^(Cal_|QC_|NORMSample_)", col)) {
    return(list(keep = FALSE, new_name = col))
  }

  list(keep = TRUE, new_name = col)
}

# Reads one filtered POS or NEG file, returns a data.frame with:
#   Feature_ID (with polarity suffix), Ionization, all metadata columns,
#   and renamed sample/field-blank columns. Cal/QC/NORMSample columns
#   dropped, multi-row sample-annotation headers dropped.
read_and_prepare <- function(fname, ionization) {
  in_path <- file.path(base_dir, fname)
  if (!file.exists(in_path)) {
    cat(sprintf("FILE NOT FOUND, skipping: %s\n", in_path)); return(NULL)
  }
  lines <- readLines(in_path, warn = FALSE)

  header_row_idx <- which(sapply(lines, function(l) {
    identical(trimws(strsplit(l, "\t", fixed = TRUE)[[1]][1]), "Alignment ID")
  }))[1]

  cols <- strsplit(lines[header_row_idx], "\t", fixed = TRUE)[[1]]
  classified <- lapply(cols, classify_col)
  keep_idx  <- vapply(classified, function(x) x$keep, logical(1))
  new_names <- vapply(classified, function(x) x$new_name, character(1))
  n_total   <- length(cols)

  data_lines <- lines[(header_row_idx + 1):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]

  rows <- lapply(data_lines, function(l) {
    fields <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(fields) < n_total) fields <- c(fields, rep("", n_total - length(fields)))
    fields[seq_len(n_total)][keep_idx]
  })
  df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(df) <- new_names[keep_idx]

  df$Feature_ID <- paste0(df[["Alignment ID"]], "_", ionization, "_SU")
  df$Ionization <- ionization
  df <- df[, c("Feature_ID", "Ionization", setdiff(names(df), c("Feature_ID","Ionization")))]

  n_dropped <- sum(!keep_idx)
  cat(sprintf("%s (%s): %d rows, %d cols kept (%d dropped as Cal/QC/NORMSample)\n",
              fname, ionization, nrow(df), sum(keep_idx), n_dropped))
  df
}

merge_and_write <- function(pos_file, neg_file, out_name) {
  pos_df <- read_and_prepare(pos_file, "POS")
  neg_df <- read_and_prepare(neg_file, "NEG")

  all_cols <- union(names(pos_df), names(neg_df))
  for (c in setdiff(all_cols, names(pos_df))) pos_df[[c]] <- NA
  for (c in setdiff(all_cols, names(neg_df))) neg_df[[c]] <- NA
  pos_df <- pos_df[, all_cols]
  neg_df <- neg_df[, all_cols]

  merged <- rbind(pos_df, neg_df)
  out_path <- file.path(base_dir, out_name)
  write.table(merged, out_path, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat(sprintf("Merged POS+NEG -> %s (%d rows, %d cols)\n", out_name, nrow(merged), ncol(merged)))
}

merge_and_write("RT_targets_annotations_POS.txt", "RT_targets_annotations_NEG.txt",
                 "RT_targets_annotations_merged.txt")
merge_and_write("Blanks_targets_annotations_POS.txt", "Blanks_targets_annotations_NEG.txt",
                 "Blanks_targets_annotations_merged.txt")

cat("\nDone.\n")

