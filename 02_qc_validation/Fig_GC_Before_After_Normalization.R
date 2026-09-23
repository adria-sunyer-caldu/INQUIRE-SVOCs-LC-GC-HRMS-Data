# =============================================================================
# Fig_GC_Before_After_Normalization.R
#
# GC-HRMS counterpart of Fig_Before_After_Normalization.R (LC). Same three
# outputs as LC: colored by Class (Indoor/Outdoor/QC, LC class_palette),
# colored by Batch ID (default ggplot hue), and both stacked. Household
# samples and QC injections are plotted together in injection order.
# Same titles, loess trend line and axis labels as the LC script.
#
# INJECTION ORDER: taken from the column order of IS_Area, which follows the
# acquisition sequence (QC injections 4CPS_PDMS_01-23 are numbered in
# injection order and sit inside the country blocks). Batch IDs 1-8 are
# numbered in run order; see Fig_GC_IS_heatmap.R for the assignment
# rule (QC -> batch of the preceding household sample).
#
# AFTER NORMALIZATION: no normalized-IS export exists for GC (LC had
# IS_after_norm_*.txt). Here each IS is divided by itself (exported area /
# exported area), which is what closest-eluting-IS normalization does to an
# IS, giving exactly 1 per IS per injection (stacked total = 8). This is a
# computed panel, not an exported one; replace with NILU's normalized IS
# export if it becomes available.
#
# GAP-FILLED VALUES: the per-IS repeated minimum constant in IS_Area is an
# MZmine gap-fill value (confirmed by NILU), not a measured area. It is
# plotted as exported in the Raw panel, where it is orders of magnitude
# below the measured areas and therefore invisible.
#
# 8 internal standards: 13C12 PBDE-183 was removed from the input file
# (gap-filled in 50% of injections, and never the closest-eluting IS for
# any feature).
#
# Classes available in IS_Area: household samples (Indoor/Outdoor) and QC.
# GC IS_Area contains no calibration (CC) or field blank (FB) injections.
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(patchwork)

base_dir <- "path/to/data"
out_dir  <- base_dir
xlsx_path <- file.path(base_dir, "20260823_INQUIRE_NILU_IS-RT-FB_data_8IS.xlsx")

boundary_override <- c()   # keep identical to Fig_GC_IS_heatmap.R

is_raw <- read_excel(xlsx_path, sheet = "IS_Area", .name_repair = "minimal")
all_cols <- names(is_raw)[-(1:4)]
is_hh <- grepl("^[A-Z]{2}_HH_[0-9]+_(IS[0-9]|OS1)$", all_cols)
is_qc <- grepl("^4CPS_PDMS_[0-9]+$", all_cols)

seq_df <- data.frame(Column = all_cols, InjectionOrder = seq_along(all_cols),
                     Country = ifelse(is_hh, substr(all_cols, 1, 2), NA),
                     stringsAsFactors = FALSE) %>%
  filter(is_hh | is_qc)
country_run_order <- unique(na.omit(seq_df$Country))
seq_df$Batch_country <- seq_df$Country
for (i in seq_len(nrow(seq_df))) {
  if (is.na(seq_df$Batch_country[i])) {
    prev <- seq_df$Batch_country[seq_len(i - 1)]
    prev <- prev[!is.na(prev)]
    seq_df$Batch_country[i] <- if (length(prev)) tail(prev, 1) else seq_df$Country[!is.na(seq_df$Country)][1]
  }
}
seq_df$Batch_ID <- as.character(match(seq_df$Batch_country, country_run_order))
for (qc in names(boundary_override)) seq_df$Batch_ID[seq_df$Column == qc] <- boundary_override[[qc]]
seq_df$Class <- ifelse(grepl("^4CPS_PDMS", seq_df$Column), "QC",
                ifelse(grepl("_IS[0-9]$", seq_df$Column), "Indoor", "Outdoor"))

before_long <- lapply(seq_len(nrow(is_raw)), function(j) {
  v_exported <- as.numeric(unlist(is_raw[j, seq_df$Column]))
  v <- v_exported
  data.frame(Name = is_raw$Name[j], Ionization = "EI",
             InjectionOrder = seq_df$InjectionOrder, Class = seq_df$Class,
             Batch_ID = seq_df$Batch_ID, Value = v, Value_exported = v_exported,
             stringsAsFactors = FALSE)
}) %>% bind_rows() %>% filter(!is.na(Value))

# Closest-eluting-IS normalization of an IS against itself: exported area /
# same exported area = 1 for every injection (no exported value is 0).
after_long <- before_long %>% mutate(Value = Value_exported / Value_exported)
stopifnot(all(after_long$Value == 1))

cat(sprintf("Before-norm: %d compounds, %d points\n", length(unique(before_long$Name)), nrow(before_long)))
cat(sprintf("After-norm:  %d compounds, %d points\n", length(unique(after_long$Name)), nrow(after_long)))
cat(sprintf("Internal standards: %d\n\n", length(unique(before_long$Name))))

trend_before <- before_long %>% group_by(Ionization, InjectionOrder) %>% summarise(Total = sum(Value), .groups = "drop")
trend_after  <- after_long  %>% group_by(Ionization, InjectionOrder) %>% summarise(Total = sum(Value), .groups = "drop")

before_long$Batch_ID <- factor(before_long$Batch_ID, levels = as.character(1:8))
after_long$Batch_ID  <- factor(after_long$Batch_ID,  levels = as.character(1:8))

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

# Same palette as LC class_palette, restricted to the classes present in GC
# IS_Area (no CC or FB injections exist in this file).
class_palette <- c(
  "Indoor"  = "#A3A500",  # olive
  "Outdoor" = "#00B0F6",  # sky blue
  "QC"      = "#0000CC"   # deep blue
)

# ---- Variant 1: colored by Class ----
p1a <- make_panel(before_long, "Class", "Raw", TRUE, trend_before, palette = class_palette)
p1b <- make_panel(after_long,  "Class", "Closest RT normalized", FALSE, trend_after, palette = class_palette)
combined1 <- p1a | p1b
print(combined1)
ggsave(file.path(out_dir, "Fig_GC_Before_After_Normalization_byClass.png"), combined1, width = 20, height = 5, dpi = 600)
cat("Saved: Fig_GC_Before_After_Normalization_byClass.png\n")

# ---- Variant 2: colored by Batch ID ----
p2a <- make_panel(before_long, "Batch_ID", "Raw", TRUE, trend_before)
p2b <- make_panel(after_long,  "Batch_ID", "Closest RT normalized", FALSE, trend_after)
combined2 <- p2a | p2b
print(combined2)
ggsave(file.path(out_dir, "Fig_GC_Before_After_Normalization_byBatch.png"), combined2, width = 20, height = 5, dpi = 600)
cat("Saved: Fig_GC_Before_After_Normalization_byBatch.png\n")

# ---- Variant 3: Class (top) stacked over Batch (bottom) ----
combined_stacked <- combined1 / combined2
print(combined_stacked)
ggsave(file.path(out_dir, "Fig_GC_Before_After_Normalization_stacked.png"), combined_stacked, width = 20, height = 10, dpi = 600)
cat("Saved: Fig_GC_Before_After_Normalization_stacked.png\n")

# Source Data
write.csv(before_long %>% rename(Value_raw = Value) %>% select(-Value_exported) %>%
            left_join(after_long %>% select(Name, InjectionOrder, Value_normalized = Value),
                      by = c("Name", "InjectionOrder")),
          file.path(out_dir, "GC_Before_After_Normalization_SourceData.csv"), row.names = FALSE)
cat("Saved underlying data: GC_Before_After_Normalization_SourceData.csv\n")
