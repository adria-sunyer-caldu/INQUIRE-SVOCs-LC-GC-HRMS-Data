# =============================================================================
# Internal standard normalization QC — RSD heatmap + before/after normalization
#
# Part 1: RSD% per internal standard (row) x batch (column), using
# QC_NORM (unspiked pool) replicates only - the cleanest true reproducibility
# signal, excluding QC_SPK_HIGH/LOW (deliberate spike-level differences,
# not replicate noise). All 42 ISs included, no exclusions.
#
# Part 2: both color variants of the "Internal standards: Raw vs Closest
# RT normalized" figure - one colored by sample Class (CC/FB/Indoor/Outdoor/
# QC), one colored by Batch ID.
#
# NOTE ON WHAT PART 2 SHOWS: for internal standards, "after normalization"
# collapses to ~1 per IS (an IS normalized by its own closest-RT IS is
# trivially 1) - confirmed directly against the real data: 99.3% of values
# are exactly 1, 0.7% are exactly 0 (where the reference IS wasn't detected
# in that injection). Stacked across all 16 (POS) / ~26 (NEG) internal
# standards, the expected flat sum is ~16 (POS) - this holds for 781/831
# injections; the remaining ~50 dip lower because at least one IS was 0
# (undetected) in that specific injection. That is the correct, expected
# behavior of this plot type - it validates the normalization procedure
# itself, not a "before vs after" comparison in the target-compound sense.
#
# This plot type is NOT meaningful for target/native compounds - those are
# expected to vary genuinely across samples (real exposure differences), so
# a flat "after" line would not be the right signal to look for there.
#
# Both parts share the same batch-lookup source (IS Areas POS/NEG.txt),
# read once here rather than twice.
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)

# ---- Set your working/data directory here ----
base_dir <- "path/to/data"
out_dir  <- base_dir

pos_path <- file.path(base_dir, "IS Areas POS.txt")
neg_path <- file.path(base_dir, "IS Areas NEG.txt")

before_pos_path <- file.path(base_dir, "IS_Areas_Before_Norm_POS.txt")
before_neg_path <- file.path(base_dir, "IS_Areas_Before_Norm_NEG.txt")
after_pos_path  <- file.path(base_dir, "IS_after_norm_POS.txt")
after_neg_path  <- file.path(base_dir, "IS_after_norm_NEG.txt")

# ---- Shared: batch lookup, built once from IS Areas POS/NEG.txt ----
# Batch ID source: IS Areas POS/NEG.txt DO have a real "Batch ID" header row;
# IS_Areas_Before_Norm_*.txt / IS_after_norm_*.txt do NOT (only 3 header
# rows: Class, Sample Type, Injection Order - confirmed directly against the
# files). The lookup is keyed by Injection Order (a plain shared number,
# verified to match 1:1 across files - unlike column names, which have
# inconsistent Cal-point formatting between the two file types and are NOT
# safe to join on).
build_batch_lookup <- function(path) {
  lines <- readLines(path, warn = FALSE)
  inj_row   <- strsplit(lines[3], "\t", fixed = TRUE)[[1]]
  batch_row <- strsplit(lines[4], "\t", fixed = TRUE)[[1]]
  ok <- grepl("^[0-9]+$", trimws(inj_row))
  setNames(trimws(batch_row[ok]), trimws(inj_row[ok]))
}
batch_lookup_pos <- build_batch_lookup(pos_path)
batch_lookup_neg <- build_batch_lookup(neg_path)

# ============================================================================
# PART 1: RSD% heatmap (QC_NORM replicates only), by batch
# ============================================================================

read_qc_norm_long <- function(path, ionization) {
  lines <- readLines(path, warn = FALSE)
  class_row <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  type_row  <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  batch_row <- strsplit(lines[4], "\t", fixed = TRUE)[[1]]
  header_idx <- which(sapply(lines, function(l) {
    identical(trimws(strsplit(l, "\t", fixed = TRUE)[[1]][1]), "Alignment ID")
  }))[1]
  cols <- strsplit(lines[header_idx], "\t", fixed = TRUE)[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]

  sample_idx <- which(trimws(type_row) == "QC" & grepl("QC_NORM", cols))

  out <- vector("list", length(data_lines))
  for (j in seq_along(data_lines)) {
    fields <- strsplit(data_lines[j], "\t", fixed = TRUE)[[1]]
    name <- fields[4]
    vals <- suppressWarnings(as.numeric(fields[sample_idx]))
    keep <- !is.na(vals)
    if (!any(keep)) next
    out[[j]] <- data.frame(
      Name = name, Ionization = ionization,
      Batch_ID = trimws(batch_row[sample_idx][keep]),
      Value = vals[keep]
    )
  }
  bind_rows(out)
}

pos_long <- read_qc_norm_long(pos_path, "POS")
neg_long <- read_qc_norm_long(neg_path, "NEG")
all_long <- bind_rows(pos_long, neg_long)

rsd_by_batch <- all_long %>%
  group_by(Name, Ionization, Batch_ID) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE), n = n(), .groups = "drop") %>%
  mutate(RSD_pct = ifelse(n >= 2 & Mean > 0, SD / Mean * 100, NA))

# Several IS names are identical between POS and NEG (e.g. Atrazine-2-hydroxy-d5,
# Oxybenzone-d3, Daidzein-d6). Sorting/leveling by Name alone would blend their
# rank across both polarities even though the data itself stays correctly
# split by the Ionization facet - fixed here by ranking within each
# Ionization independently, using a composite label so row identity is
# unambiguous either way.
rsd_by_batch <- rsd_by_batch %>%
  mutate(Name_label = paste0(Name, " (", Ionization, ")"))

name_order <- rsd_by_batch %>%
  group_by(Name_label, Ionization) %>% summarise(mean_rsd = mean(RSD_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(Ionization, desc(mean_rsd)) %>% pull(Name_label)
rsd_by_batch$Name_label <- factor(rsd_by_batch$Name_label, levels = rev(name_order))
rsd_by_batch$Batch_ID <- factor(rsd_by_batch$Batch_ID,
                                 levels = as.character(sort(as.numeric(unique(rsd_by_batch$Batch_ID)))))

p_heatmap <- ggplot(rsd_by_batch, aes(x = Batch_ID, y = Name_label, fill = RSD_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(RSD_pct, 0)), size = 4, color = "black") +
  facet_wrap(~Ionization, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey90", name = "RSD %") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 10), panel.grid = element_blank()) +
  labs(x = "Batch ID", y = NULL,
       title = "RSD% per internal standard (QC_NORM replicates only), by batch")
print(p_heatmap)
ggsave(file.path(out_dir, "Fig_IS_heatmap_QCnorm.png"), p_heatmap, width = 15, height = 10, dpi = 600)
cat("Saved: Fig_IS_heatmap_QCnorm.png\n")

write.csv(rsd_by_batch, file.path(out_dir, "IS_QCnorm_RSD_by_batch.csv"), row.names = FALSE)
cat("Saved underlying data: IS_QCnorm_RSD_by_batch.csv\n\n")

# ============================================================================
# PART 2: Before/after normalization, by Class and by Batch
# ============================================================================

read_long <- function(path, ionization, batch_lookup) {
  lines <- readLines(path, warn = FALSE)
  class_row <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  type_row  <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  inj_row   <- strsplit(lines[3], "\t", fixed = TRUE)[[1]]
  header_idx <- which(sapply(lines, function(l) {
    identical(trimws(strsplit(l, "\t", fixed = TRUE)[[1]][1]), "Alignment ID")
  }))[1]
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[nchar(trimws(data_lines)) > 0]
  sample_idx <- which(trimws(type_row) %in% c("Sample", "QC", "Blank", "Standard"))

  out <- vector("list", length(data_lines))
  for (j in seq_along(data_lines)) {
    fields <- strsplit(data_lines[j], "\t", fixed = TRUE)[[1]]
    name <- fields[4]
    vals <- suppressWarnings(as.numeric(fields[sample_idx]))
    keep <- !is.na(vals)
    if (!any(keep)) next
    inj_vals <- trimws(inj_row[sample_idx][keep])
    out[[j]] <- data.frame(
      Name = name, Ionization = ionization,
      InjectionOrder = suppressWarnings(as.numeric(inj_vals)),
      Class = trimws(class_row[sample_idx][keep]),
      Batch_ID = unname(batch_lookup[inj_vals]),   # looked up by injection order, not read directly
      Value = vals[keep]
    )
  }
  bind_rows(out)
}

before_long <- bind_rows(read_long(before_pos_path, "POS", batch_lookup_pos),
                          read_long(before_neg_path, "NEG", batch_lookup_neg))
after_long  <- bind_rows(read_long(after_pos_path, "POS", batch_lookup_pos),
                          read_long(after_neg_path, "NEG", batch_lookup_neg))

n_missing_batch <- sum(is.na(before_long$Batch_ID)) + sum(is.na(after_long$Batch_ID))
if (n_missing_batch > 0) {
  cat(sprintf("WARNING: %d rows could not be matched to a Batch ID by injection order - check that\n", n_missing_batch))
  cat("IS Areas POS.txt / IS Areas NEG.txt are the correct, unmodified batch-lookup source files.\n\n")
}

cat(sprintf("Before-norm: %d compounds, %d points\n", length(unique(before_long$Name)), nrow(before_long)))
cat(sprintf("After-norm:  %d compounds, %d points\n", length(unique(after_long$Name)), nrow(after_long)))

frac_zero_or_one <- mean(after_long$Value == 0 | after_long$Value == 1, na.rm = TRUE)
cat(sprintf("%.1f%% of after-norm values are exactly 0 or 1 (expected for IS-on-IS normalization)\n\n",
            frac_zero_or_one * 100))

# Trend line data: total stacked value per injection order (summed across all
# IS), computed SEPARATELY per Ionization - POS and NEG injection order
# numbers overlap completely (both run 1-831, confirmed against the real
# files), so combining them without this split would sum unrelated POS/NEG
# injections together at the same x-position, in both the bars and the trend line.
trend_before <- before_long %>% group_by(Ionization, InjectionOrder) %>% summarise(Total = sum(Value), .groups = "drop")
trend_after  <- after_long  %>% group_by(Ionization, InjectionOrder) %>% summarise(Total = sum(Value), .groups = "drop")

make_panel <- function(data, color_var, title_str, show_legend, trend_data, palette = NULL) {
  p <- ggplot(data, aes(x = InjectionOrder * 2, y = Value, fill = .data[[color_var]])) +
    geom_col(width = 2, show.legend = show_legend) +
    geom_smooth(data = trend_data, aes(x = InjectionOrder * 2, y = Total), inherit.aes = FALSE,
                color = "black", se = FALSE, method = "loess", linewidth = 0.6) +
    facet_wrap(~Ionization, ncol = 1, scales = "free_y") +
    scale_x_continuous(labels = function(x) x / 2) +
    theme_minimal(base_size = 11) +
    labs(title = title_str, x = "Injection Order", y = if (title_str == "Raw") "Value" else NULL,
         fill = color_var)
  if (!is.null(palette)) p <- p + scale_fill_manual(values = palette)
  p
}

# Explicit, high-contrast palette for Class - varies BOTH hue and luminance,
# not just hue, so thin interleaved bars of similar-brightness colors (e.g.
# the original olive Indoor vs sky-blue Outdoor) don't visually blend
# together. Loosely Okabe-Ito-based (colorblind-safe, high separation).
class_palette <- c(
  "CC"      = "#E76BF3",  # orchid/magenta
  "FB"      = "#F8766D",  # red
  "Indoor"  = "#A3A500",  # olive
  "Outdoor" = "#00B0F6",  # sky blue
  "QC"      = "#0000CC"   # deep blue
)

# ---- Variant 1: colored by Class ----
p1a <- make_panel(before_long, "Class", "Raw", TRUE, trend_before, palette = class_palette)
p1b <- make_panel(after_long,  "Class", "Closest RT normalized", FALSE, trend_after, palette = class_palette)
combined1 <- p1a | p1b
print(combined1)
ggsave(file.path(out_dir, "Fig_Before_After_Normalization_byClass.png"), combined1, width = 20, height = 8, dpi = 600)
cat("Saved: Fig_Before_After_Normalization_byClass.png\n")

# ---- Variant 2: colored by Batch ID ----
p2a <- make_panel(before_long, "Batch_ID", "Raw", TRUE, trend_before)
p2b <- make_panel(after_long,  "Batch_ID", "Closest RT normalized", FALSE, trend_after)
combined2 <- p2a | p2b
print(combined2)
ggsave(file.path(out_dir, "Fig_Before_After_Normalization_byBatch.png"), combined2, width = 20, height = 8, dpi = 600)
cat("Saved: Fig_Before_After_Normalization_byBatch.png\n")

# ---- Variant 3: Class (top) stacked over Batch (bottom) ----
combined_stacked <- combined1 / combined2
print(combined_stacked)
ggsave(file.path(out_dir, "Fig_Before_After_Normalization_stacked.png"), combined_stacked, width = 20, height = 16, dpi = 600)
cat("Saved: Fig_Before_After_Normalization_stacked.png\n")
