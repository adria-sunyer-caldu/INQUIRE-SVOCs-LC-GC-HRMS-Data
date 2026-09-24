# =============================================================================
# Fig_GC_Before_After_Normalization.R
#
# GC-HRMS internal standard signal before and after closest-eluting-IS
# normalization. Household samples and QC injections are plotted together
# in injection order.
#
# BATCH ASSIGNMENT: read from the NILU injection-batch export
# (dr21__nts_svoc_pdms_nilu__04_INQUIRE_samples_ID.csv, columns
# lab_sample_id, batch_inj), 16 true instrument batches, not the 8
# country-inferred batches used in Fig_GC_IS_heatmap.R.
#
# AFTER NORMALIZATION: no normalized-IS export exists for GC. Each IS is
# divided by itself (raw area / raw area), which is what closest-eluting-IS
# normalization does to an IS, giving exactly 1 per IS per injection. This
# is a computed panel, not an exported one; replace with NILU's normalized
# IS export if it becomes available.
#
# GAP-FILLED VALUES: the per-IS repeated minimum constant in IS_Area is an
# MZmine gap-fill value, not a measured area. It is plotted as exported in
# the Raw panel, where it is orders of magnitude below the measured areas
# and therefore invisible.
#
# 8 internal standards: 13C12 PBDE-183 was removed from the input file
# (gap-filled in 50% of injections, and never the closest-eluting IS for
# any feature).
#
# Classes available in IS_Area: household samples (Indoor/Outdoor) and QC.
# GC IS_Area contains no calibration (CC) or field blank (FB) injections.
#
# Source data output matches the "Suppl. Figure 42" sheet format:
# Name, Ionization, InjectionOrder, Class, Batch_ID, Value_raw,
# Value_normalized.
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

base_dir       <- "path/to/data"
out_dir        <- base_dir
xlsx_path      <- file.path(base_dir, "20260823_INQUIRE_NILU_IS-RT-FB_data_8IS.xlsx")
batch_csv_path <- file.path(base_dir, "dr21__nts_svoc_pdms_nilu__04_INQUIRE_samples_ID.csv")

# ---- Read IS areas ----
is_raw   <- read_excel(xlsx_path, sheet = "IS_Area", .name_repair = "minimal")
all_cols <- names(is_raw)[-(1:4)]
is_hh    <- grepl("^[A-Z]{2}_HH_[0-9]+_(IS[0-9]|OS1)$", all_cols)
is_qc    <- grepl("^4CPS_PDMS_[0-9]+$", all_cols)

# ---- Batch assignment from CSV ----
batch_map <- read.csv(batch_csv_path, stringsAsFactors = FALSE)

all_batch <- batch_map %>%
  select(lab_sample_id, batch_inj) %>%
  mutate(Batch_ID = as.integer(gsub("batch", "", batch_inj)))

seq_df <- data.frame(
  Column         = all_cols[is_hh | is_qc],
  InjectionOrder = which(is_hh | is_qc),
  stringsAsFactors = FALSE
) %>%
  left_join(all_batch %>% rename(Column = lab_sample_id), by = "Column") %>%
  mutate(
    Class    = ifelse(grepl("^4CPS_PDMS", Column), "QC",
               ifelse(grepl("_IS[0-9]$",  Column), "Indoor", "Outdoor")),
    Batch_ID = as.character(Batch_ID)
  )

cat("Batch assignment summary:\n")
print(table(seq_df$Batch_ID, seq_df$Class, useNA = "ifany"))

# ---- Long format ----
before_long <- lapply(seq_len(nrow(is_raw)), function(j) {
  v <- as.numeric(unlist(is_raw[j, seq_df$Column]))
  data.frame(Name           = is_raw$Name[j],
             Ionization     = "EI",
             InjectionOrder = seq_df$InjectionOrder,
             Class          = seq_df$Class,
             Batch_ID       = seq_df$Batch_ID,
             Value_raw      = v,
             stringsAsFactors = FALSE)
}) %>% bind_rows() %>% filter(!is.na(Value_raw))

after_long <- before_long %>%
  mutate(Value_normalized = Value_raw / Value_raw)   # IS normalized against itself = 1

cat(sprintf("Before-norm: %d IS, %d points\n", length(unique(before_long$Name)), nrow(before_long)))
cat(sprintf("After-norm:  %d IS, %d points\n", length(unique(after_long$Name)), nrow(after_long)))

trend_before <- before_long %>%
  group_by(Ionization, InjectionOrder) %>%
  summarise(Total = sum(Value_raw), .groups = "drop")
trend_after <- after_long %>%
  group_by(Ionization, InjectionOrder) %>%
  summarise(Total = sum(Value_normalized), .groups = "drop")

before_long$Batch_ID <- factor(before_long$Batch_ID, levels = as.character(1:16))
after_long$Batch_ID  <- factor(after_long$Batch_ID,  levels = as.character(1:16))

make_panel <- function(data, y_var, color_var, title_str, show_legend,
                       trend_data, y_trend, palette = NULL) {
  p <- ggplot(data, aes(x = InjectionOrder * 2,
                        y = .data[[y_var]],
                        fill = .data[[color_var]])) +
    geom_col(width = 2, show.legend = show_legend) +
    geom_smooth(data = trend_data,
                aes(x = InjectionOrder * 2, y = .data[[y_trend]]),
                inherit.aes = FALSE,
                color = "black", se = FALSE, method = "loess", linewidth = 0.6) +
    facet_wrap(~Ionization, ncol = 1, scales = "free_y") +
    scale_x_continuous(labels = function(x) x / 2) +
    theme_minimal(base_size = 11) +
    labs(title = title_str,
         x     = "Injection Order",
         y     = if (title_str == "Raw") "IS peak area (raw)" else "Normalized value",
         fill  = color_var)
  if (!is.null(palette)) p <- p + scale_fill_manual(values = palette)
  p
}

class_palette <- c("Indoor" = "#A3A500", "Outdoor" = "#00B0F6", "QC" = "#0000CC")

# ---- Variant 1: by Class ----
p1a <- make_panel(before_long, "Value_raw",        "Class", "Raw",
                  TRUE,  trend_before, "Total", palette = class_palette)
p1b <- make_panel(after_long,  "Value_normalized", "Class", "Closest RT normalized",
                  FALSE, trend_after,  "Total", palette = class_palette)
combined1 <- p1a | p1b
print(combined1)
ggsave(file.path(out_dir, "Fig_GC_Before_After_Normalization_byClass.png"),
       combined1, width = 20, height = 5, dpi = 600)
cat("Saved: Fig_GC_Before_After_Normalization_byClass.png\n")

# ---- Variant 2: by Batch ----
p2a <- make_panel(before_long, "Value_raw",        "Batch_ID", "Raw",
                  TRUE,  trend_before, "Total")
p2b <- make_panel(after_long,  "Value_normalized", "Batch_ID", "Closest RT normalized",
                  FALSE, trend_after,  "Total")
combined2 <- p2a | p2b
print(combined2)
ggsave(file.path(out_dir, "Fig_GC_Before_After_Normalization_byBatch.png"),
       combined2, width = 20, height = 5, dpi = 600)
cat("Saved: Fig_GC_Before_After_Normalization_byBatch.png\n")

# ---- Source data — matches Suppl. Figure 42 sheet format exactly ----
source_data_42 <- before_long %>%
  left_join(after_long %>% select(Name, InjectionOrder, Value_normalized),
            by = c("Name", "InjectionOrder")) %>%
  select(Name, Ionization, InjectionOrder, Class, Batch_ID,
         Value_raw, Value_normalized)

write.csv(source_data_42,
          file.path(out_dir, "GC_Before_After_Normalization_SourceData.csv"),
          row.names = FALSE)
cat("Saved: GC_Before_After_Normalization_SourceData.csv\n")
cat("  -> Replace sheet 'Suppl. Figure 42' in Source_Data_INQUIRE_FINAL.xlsx with this file.\n")
