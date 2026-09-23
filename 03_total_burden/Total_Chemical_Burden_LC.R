# ============================================================
# TOTAL CHEMICAL BURDEN ANALYSIS
# INQUIRE LC-HRMS Nontarget Dataset (ESI Positive)
#
# Approach: Sum of blank-subtracted, IS-normalized peak areas
#           per sample as a proxy for total chemical burden.
#
# Design: Paired indoor/outdoor (same household, same month).
#         Statistics: Wilcoxon signed-rank test (paired).
#         Multiple testing: Benjamini-Hochberg correction.
#
# Input: MS-DIAL output, Zenodo LC deposit (Excel)
#   - Rows 1-3: metadata (Class, Sample Type, Injection Order)
#   - Row 4:    column headers
#   - Row 5+:   feature data
#   - Cols 1-AF: feature metadata
#   - Cols AG+:  sample peak areas (already IS-normalized, blank-subtracted)
#   - Sample names are pseudonymized: {Country}_{Indoor|Outdoor}_{code};
#     indoor and outdoor samples of the same household share the same code,
#     extra indoor samplers carry a _2/_3 suffix and remain unpaired.
#
# Output: 4 publication-quality figures + 3 CSV result tables
# ============================================================

library(data.table)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(patchwork)

# ============================================================
# CONFIGURATION — EDIT THESE
# ============================================================

data_dir <- "path/to/data"
file_path <- file.path(data_dir, "LC_All_Features_FULL_final_zenodo.xlsx")   # Zenodo LC deposit

output_dir <- file.path(data_dir, "Total_Chemical_Burden_Output")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ============================================================
# STEP 1: READ FILE STRUCTURE AND EXTRACT CLASS LABELS
# ============================================================

cat("Step 1: Reading file structure...\n")

# Read only first 4 rows to get class labels and column names
# This is fast even for a large file because we stop after 4 rows
header_raw <- as.data.frame(read_excel(
  file_path,
  sheet        = "Full_dataset_INQUIRE_nontarget_",
  col_names    = FALSE,
  n_max        = 4,
  .name_repair = "minimal"
))

# Row 4 contains the actual column headers (Alignment ID, Average Rt, ..., Sample names)
col_names_all <- trimws(as.character(unlist(header_raw[4, ])))

# Identify sample columns by their pseudonymized code pattern
sample_col_indices <- which(grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?$", col_names_all))
cat(paste("  Total columns in file:         ", ncol(header_raw), "\n"))
cat(paste("  Sample columns found:          ", length(sample_col_indices), "\n"))
cat(paste("  First sample column index:     ", min(sample_col_indices), "\n"))

sample_names_all <- col_names_all[sample_col_indices]

# Extract class labels (Indoor/Outdoor) from Row 1, same column positions
class_labels_all <- trimws(as.character(unlist(header_raw[1, sample_col_indices])))
cat(paste("  Indoor samples:                ", sum(class_labels_all == "Indoor"), "\n"))
cat(paste("  Outdoor samples:               ", sum(class_labels_all == "Outdoor"), "\n"))
cat(paste("  Other (blanks/QCs/excluded):   ",
          sum(!class_labels_all %in% c("Indoor", "Outdoor")), "\n"))

# Keep only Indoor and Outdoor
keep_idx     <- which(class_labels_all %in% c("Indoor", "Outdoor"))
sample_names <- sample_names_all[keep_idx]
class_labels <- class_labels_all[keep_idx]

cat(paste("  Samples retained for analysis: ", length(sample_names), "\n\n"))

# ============================================================
# STEP 2: READ FEATURE AREA DATA (SAMPLE COLUMNS ONLY)
# ============================================================

cat("Step 2: Reading feature data (this may take 1-2 minutes for a large file)...\n")

# Read only column 1 (Alignment ID) + all sample columns
# skip=3 skips the 3 metadata rows; header=TRUE reads row 4 as header
cols_to_read <- c(1, sample_col_indices)

dt <- as.data.table(read_excel(
  file_path,
  sheet        = "Full_dataset_INQUIRE_nontarget_",
  skip         = 3,
  col_names    = TRUE,
  na           = c("", "NA", "null", "NULL"),
  .name_repair = "minimal"
))
dt <- dt[, cols_to_read, with = FALSE]
setnames(dt, c("Feature_ID", col_names_all[sample_col_indices]))
cat(paste("  Features loaded:", nrow(dt), "\n\n"))

# ============================================================
# STEP 3: KEEP ONLY RETAINED SAMPLES, REPLACE NA WITH 0
# ============================================================

cat("Step 3: Cleaning data...\n")

# Keep only the Indoor/Outdoor samples (drop blanks/QCs)
keep_cols <- c("Feature_ID", sample_names)
dt <- dt[, ..keep_cols]

# Replace NA with 0 (not detected = 0 area after blank subtraction)
for (col in sample_names) {
  set(dt, which(is.na(dt[[col]])), col, 0)
}

# Coerce to numeric (safety check)
for (col in sample_names) {
  set(dt, j = col, value = as.numeric(dt[[col]]))
}

cat("  Data cleaned.\n\n")

# ============================================================
# STEP 4: COMPUTE TOTAL CHEMICAL BURDEN PER SAMPLE
# ============================================================

cat("Step 4: Computing total chemical burden per sample...\n")

# Sum all feature areas per sample as a proxy for total chemical load
burden_vec <- dt[, lapply(.SD, function(x) sum(x, na.rm = TRUE)),
                 .SDcols = sample_names]

# Build metadata table
burden_dt <- data.table(
  Sample  = sample_names,
  Class   = class_labels,
  Burden  = as.numeric(burden_vec[1, ])
)

# Extract country (prefix of the sample name)
burden_dt[, Country := sub("_.*", "", Sample)]

# Extract household code (shared by the indoor and outdoor sample)
burden_dt[, House := sub("^[A-Z]{2}_(Indoor|Outdoor)_", "", Sample)]

# Verify extraction
cat("  Burden computed. Summary per country and class:\n")
print(burden_dt[, .(
  N        = .N,
  Median   = round(median(Burden), 0),
  Min      = round(min(Burden), 0),
  Max      = round(max(Burden), 0)
), by = .(Country, Class)][order(Country, Class)])

cat("\n")

# ============================================================
# STEP 5: BUILD PAIRED DATASET (one row per household)
# ============================================================

cat("Step 5: Building paired dataset...\n")

# Reshape to wide: one row per Country + House pair
burden_wide <- dcast(
  burden_dt,
  Country + House ~ Class,
  value.var = "Burden",
  fun.aggregate = mean
)

# Remove incomplete pairs (houses missing either indoor or outdoor)
n_before <- nrow(burden_wide)
burden_wide <- na.omit(burden_wide)
n_after  <- nrow(burden_wide)
cat(paste("  Pairs before NA removal:", n_before, "\n"))
cat(paste("  Complete pairs retained:", n_after, "\n\n"))

# Derived metrics per pair
burden_wide[, Log2FC       := log2(Indoor / Outdoor)]
burden_wide[, IndoorRatio  := Indoor / (Indoor + Outdoor)]
burden_wide[, Indoor_Higher:= Indoor > Outdoor]

# ============================================================
# STEP 6: STATISTICAL ANALYSIS — OVERALL
# ============================================================

cat("==============================================================\n")
cat("STATISTICAL RESULTS — OVERALL (all countries combined)\n")
cat("==============================================================\n")

wtest_overall <- wilcox.test(
  burden_wide$Indoor,
  burden_wide$Outdoor,
  paired      = TRUE,
  alternative = "two.sided",
  exact       = FALSE
)

n_total   <- nrow(burden_wide)
# Rank-biserial correlation as effect size for paired Wilcoxon
r_rb_overall <- (2 * wtest_overall$statistic) / (n_total * (n_total + 1) / 2) - 1

cat(paste("Test:                      Paired Wilcoxon signed-rank\n"))
cat(paste("N pairs:                  ", n_total, "\n"))
cat(paste("V statistic:              ", wtest_overall$statistic, "\n"))
cat(paste("p-value:                  ", format(wtest_overall$p.value, scientific = TRUE, digits = 3), "\n"))
cat(paste("Effect size (r_rb):       ", round(r_rb_overall, 3), "\n"))
cat(paste("Median indoor burden:     ", format(round(median(burden_wide$Indoor)), big.mark = ","), "\n"))
cat(paste("Median outdoor burden:    ", format(round(median(burden_wide$Outdoor)), big.mark = ","), "\n"))
cat(paste("Median log2 FC:           ", round(median(burden_wide$Log2FC), 3), "\n"))
cat(paste("Median fold change (x):   ", round(2^median(burden_wide$Log2FC), 2), "\n"))
cat(paste("% pairs Indoor > Outdoor: ", round(mean(burden_wide$Indoor_Higher) * 100, 1), "%\n"))
cat("\n")

# ============================================================
# STEP 7: STATISTICAL ANALYSIS — BY COUNTRY
# ============================================================

cat("==============================================================\n")
cat("STATISTICAL RESULTS — BY COUNTRY\n")
cat("==============================================================\n")

countries <- sort(unique(burden_wide$Country))

country_stats <- rbindlist(lapply(countries, function(ctry) {
  sub <- burden_wide[Country == ctry]
  n   <- nrow(sub)
  
  if (n < 5) {
    warning(paste("Skipping", ctry, "— only", n, "pairs"))
    return(NULL)
  }
  
  wt <- wilcox.test(sub$Indoor, sub$Outdoor,
                    paired = TRUE, exact = FALSE, alternative = "two.sided")
  r_rb <- (2 * wt$statistic) / (n * (n + 1) / 2) - 1
  
  data.table(
    Country           = ctry,
    N_pairs           = n,
    Median_Indoor     = median(sub$Indoor),
    Median_Outdoor    = median(sub$Outdoor),
    Median_Log2FC     = round(median(sub$Log2FC), 3),
    Median_FC         = round(2^median(sub$Log2FC), 2),
    Pct_Indoor_Higher = round(mean(sub$Indoor_Higher) * 100, 1),
    W_statistic       = wt$statistic,
    P_value           = wt$p.value,
    Effect_size_r     = round(r_rb, 3)
  )
}))

# BH multiple testing correction
country_stats[, P_adj := p.adjust(P_value, method = "BH")]
country_stats[, Sig   := fcase(
  P_adj < 0.001, "***",
  P_adj < 0.01,  "**",
  P_adj < 0.05,  "*",
  default        = "ns"
)]

print(country_stats[, .(Country, N_pairs, Median_FC, Pct_Indoor_Higher,
                        P_value = formatC(P_value, format = "e", digits = 2),
                        P_adj   = formatC(P_adj, format = "e", digits = 2),
                        Sig, Effect_size_r)])
cat("\n")

# ============================================================
# STEP 8: SAVE STATISTICAL RESULTS
# ============================================================

fwrite(country_stats, "TotalBurden_Statistics_ByCountry.csv")
fwrite(burden_wide,   "TotalBurden_PairedBurden_AllHouseholds.csv")

overall_out <- data.table(
  Test               = "Paired Wilcoxon signed-rank (all countries)",
  N_pairs            = n_total,
  V_statistic        = wtest_overall$statistic,
  P_value            = wtest_overall$p.value,
  Effect_size_r      = round(r_rb_overall, 3),
  Median_Log2FC      = round(median(burden_wide$Log2FC), 3),
  Median_FC          = round(2^median(burden_wide$Log2FC), 2),
  Pct_Indoor_Higher  = round(mean(burden_wide$Indoor_Higher) * 100, 1)
)
fwrite(overall_out, "TotalBurden_Statistics_Overall.csv")

# ============================================================
# STEP 9: PLOT 1 — OVERALL Indoor vs Outdoor
# ============================================================

cat("Step 9: Generating plots...\n")

col_in  <- "#B03A2E"
col_out <- "#1F618D"

# Long format for plotting
plot_long <- burden_dt[Class %in% c("Indoor", "Outdoor")]
plot_long[, Class        := factor(Class, levels = c("Outdoor", "Indoor"))]
plot_long[, Log10_Burden := log10(Burden + 1)]

# P-value label
p_val_str <- ifelse(
  wtest_overall$p.value < 2.2e-16, "p < 2.2\u00D710\u207B\u00B9\u2076",
  paste0("p = ", formatC(wtest_overall$p.value, format = "e", digits = 2))
)
fc_str    <- paste0(round(2^median(burden_wide$Log2FC), 1), "\u00D7 higher indoors")
subtitle_overall <- paste0(
  "n = ", n_total, " paired households | ",
  fc_str, " (median) | ", p_val_str
)

p1 <- ggplot(plot_long, aes(x = Class, y = Log10_Burden, fill = Class, color = Class)) +
  geom_violin(alpha = 0.25, trim = FALSE, linewidth = 0.7, adjust = 1.2) +
  geom_boxplot(width = 0.12, alpha = 0.9, outlier.shape = NA,
               linewidth = 0.6, color = "grey15") +
  geom_beeswarm(size = 1.4, alpha = 0.55, cex = 1.6, priority = "random") +
  scale_fill_manual(values  = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_color_manual(values = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_y_continuous(
    labels = function(x) format(10^x, scientific = FALSE, big.mark = ","),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title    = "Total SVOC Chemical Burden",
    subtitle = subtitle_overall,
    x        = NULL,
    y        = "Summed peak area (IS-normalized, blank-subtracted)",
    caption  = "Each point = one air sample. Box = IQR; whiskers = 1.5\u00D7IQR.\nData: LC-HRMS ESI(+), IS-normalized, blank-subtracted."
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position     = "none",
    plot.title          = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle       = element_text(size = 10, color = "grey30", hjust = 0),
    plot.caption        = element_text(size = 8, color = "grey50"),
    axis.text.x         = element_text(size = 13, face = "bold"),
    panel.grid.major.y  = element_line(color = "grey90", linewidth = 0.4)
  )

print(p1)
cat("  Saved: Fig_Chemical_Burden_Overall.png\n")

# ============================================================
# STEP 10: PLOT 2 — BY COUNTRY (with paired lines)
# ============================================================
# Long format for connecting lines
pairs_long <- melt(
  burden_wide[, .(Country, House,
                  Outdoor = log10(Outdoor + 1),
                  Indoor  = log10(Indoor  + 1))],
  id.vars       = c("Country", "House"),
  variable.name = "Class",
  value.name    = "Log10_Burden"
)
pairs_long[, Class := factor(Class, levels = c("Outdoor", "Indoor"))]

# Add fold change and significance label per country (for facet label)
label_dt <- country_stats[, .(
  Country,
  facet_label = paste0(Country, "  (", Median_FC, "\u00D7  ", Sig, ")")
)]

plot_long_c <- merge(burden_dt[Class %in% c("Indoor", "Outdoor")],
                     label_dt, by = "Country")
plot_long_c[, Class        := factor(Class, levels = c("Outdoor", "Indoor"))]
plot_long_c[, Log10_Burden := log10(Burden + 1)]

pairs_long <- merge(pairs_long, label_dt, by = "Country")

# Order countries by fold change (descending)
fc_order <- country_stats[order(-Median_FC), Country]
plot_long_c[, facet_label := factor(
  facet_label,
  levels = label_dt[match(fc_order, Country), facet_label]
)]
pairs_long[, facet_label := factor(
  facet_label,
  levels = label_dt[match(fc_order, Country), facet_label]
)]

p2 <- ggplot(plot_long_c, aes(x = Class, y = Log10_Burden,
                              fill = Class, color = Class)) +
  geom_line(
    data        = pairs_long,
    aes(x       = Class, y = Log10_Burden, group = House),
    color       = "grey65", alpha = 0.30, linewidth = 0.35,
    inherit.aes = FALSE
  ) +
  geom_violin(alpha = 0.22, trim = FALSE, linewidth = 0.55, adjust = 1.2) +
  geom_boxplot(width = 0.12, alpha = 0.88, outlier.shape = NA,
               linewidth = 0.55, color = "grey15") +
  geom_beeswarm(size = 1.0, alpha = 0.60, cex = 1.3, priority = "random") +
  scale_fill_manual(values  = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_color_manual(values = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_y_continuous(
    labels = function(x) scales::scientific(10^x),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  facet_wrap(~ facet_label, nrow = 2, scales = "fixed") +
  labs(
    title    = "Total SVOC Chemical Burden by Country",
    subtitle = "LC-HRMS ESI(+) | Ordered by indoor/outdoor fold change (FC)\nGrey lines connect paired samples from the same household",
    x        = NULL,
    y        = expression(Summed~peak~area~(log[10]~scale)),
    caption  = "FC = median fold change Indoor/Outdoor. Significance: * p<0.05, ** p<0.01, *** p<0.001 (BH-corrected paired Wilcoxon)."
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position     = "none",
    plot.title          = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle       = element_text(size = 9, color = "grey30", hjust = 0),
    plot.caption        = element_text(size = 7.5, color = "grey50"),
    strip.background    = element_rect(fill = "grey96", color = "grey70"),
    strip.text          = element_text(face = "bold", size = 9),
    panel.grid.major.y  = element_line(color = "grey90", linewidth = 0.35),
    panel.spacing       = unit(0.9, "lines"),
    axis.text.x         = element_text(size = 10, face = "bold")
  )
print(p2)
cat("  Saved: Fig_Chemical_Burden_ByCountry.png\n")

# ============================================================
# STEP 11: PLOT 3 — FOLD CHANGE SUMMARY BY COUNTRY
# ============================================================

# Per-household Log2FC boxplot by country (ordered by median)
fc_plot_data <- burden_wide[, .(Country, House, Log2FC)]
fc_plot_data[, Country := factor(Country,
                                   levels = country_stats[order(Median_Log2FC), Country])]

p3 <- ggplot(fc_plot_data, aes(x = Country, y = Log2FC)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey40", linewidth = 0.6) +
  geom_violin(fill = col_in, color = col_in, alpha = 0.22,
              trim = FALSE, linewidth = 0.55, adjust = 1.2) +
  geom_boxplot(width = 0.22, fill = col_in, color = "grey15",
               alpha = 0.85, outlier.shape = 21,
               outlier.size = 1.2, outlier.alpha = 0.5,
               linewidth = 0.55) +
  geom_beeswarm(color = col_in, alpha = 0.55, size = 1.2,
                cex = 1.3, priority = "random") +
  geom_text(
    data = country_stats[order(Median_Log2FC)],
    aes(x = Country, y = max(burden_wide$Log2FC) * 1.05,
        label = paste0(Median_FC, "\u00D7\n", Sig)),
    size = 3.0, fontface = "bold", color = "grey20",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    name     = expression(log[2]~fold~change~(Indoor/Outdoor)),
    sec.axis = sec_axis(
      ~ 2^.,
      name   = "Fold change (Indoor/Outdoor)",
      breaks = c(0.5, 1, 2, 4, 8, 16),
      labels = c("0.5\u00D7", "1\u00D7", "2\u00D7", "4\u00D7", "8\u00D7", "16\u00D7")
    )
  ) +
  labs(
    title    = "Indoor/Outdoor Fold Change in Total Chemical Burden per Household",
    subtitle = "Paired Wilcoxon signed-rank | BH-corrected | * p<0.05, ** p<0.01, *** p<0.001",
    x        = NULL,
    caption  = "Each violin/boxplot summarises log\u2082(Indoor/Outdoor) across all household pairs per country.\nDashed line: indoor = outdoor (log\u2082FC = 0)."
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle      = element_text(size = 9, color = "grey30"),
    plot.caption       = element_text(size = 8, color = "grey50"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    axis.text.x        = element_text(face = "bold", size = 11)
  )
print(p3)
cat("  Saved: Fig_ChemBurden_FoldChange_ByCountry.png\n")

# ============================================================
# STEP 12: LC-ONLY — HORIZONTAL, 8 COUNTRIES AS ROWS,
# INDOOR + OUTDOOR TOGETHER IN THE SAME ROW
# ============================================================
# One row per country, x-axis = burden. Indoor and Outdoor both
# plotted within that same row (small vertical offset), not as
# separate panels or facets.

cat("Step 12: Building horizontal by-country plot (LC only, Indoor+Outdoor combined per row)...\n")

lc_row_data <- burden_wide[, .(Country, House,
                                Outdoor = log10(Outdoor + 1),
                                Indoor  = log10(Indoor  + 1))]
lc_row_long <- melt(lc_row_data, id.vars = c("Country", "House"),
                     variable.name = "Class", value.name = "Log10_Burden")
lc_row_long[, Class := factor(Class, levels = c("Outdoor", "Indoor"))]

# Country order: by median fold change (reuse country_stats already computed)
country_order_h <- country_stats[order(-Median_FC), Country]
lc_row_long[, Country := factor(Country, levels = rev(country_order_h))]

lc_row_long[, y_num := as.numeric(Country)]
lc_row_long[, y_off := ifelse(Class == "Indoor", y_num + 0.18, y_num - 0.18)]

build_density_polygon <- function(vals, baseline, scale = 1, half = FALSE, side = 1) {
  d <- density(vals, n = 512)
  h <- d$y / max(d$y) * scale
  if (half) {
    data.frame(x = c(d$x, rev(d$x)), y = c(baseline + h * side, rep(baseline, length(d$x))))
  } else {
    data.frame(x = d$x, y = baseline + h)
  }
}

half_polys <- rbindlist(lapply(unique(paste(lc_row_long$Country, lc_row_long$Class)), function(g) {
  sub <- lc_row_long[paste(Country, Class) == g]
  if (nrow(sub) < 2) return(NULL)
  poly <- as.data.table(build_density_polygon(sub$Log10_Burden, baseline = sub$y_off[1],
                                                scale = 0.16, half = TRUE, side = 1))
  poly[, Country := sub$Country[1]]
  poly[, Class   := sub$Class[1]]
  poly
}))

box_stats_h <- lc_row_long[, .(
  ymin  = quantile(Log10_Burden, 0.25) - 1.5 * IQR(Log10_Burden),
  lower = quantile(Log10_Burden, 0.25),
  med   = median(Log10_Burden),
  upper = quantile(Log10_Burden, 0.75),
  ymax  = quantile(Log10_Burden, 0.75) + 1.5 * IQR(Log10_Burden)
), by = .(Country, Class, y_off)]

p_horizontal <- ggplot() +
  geom_polygon(data = half_polys, aes(x = x, y = y, group = interaction(Country, Class), fill = Class),
               alpha = 0.4, color = NA) +
  geom_boxplot(data = box_stats_h, aes(xmin = ymin, xlower = lower, xmiddle = med, xupper = upper, xmax = ymax,
                                         y = y_off - 0.05, group = interaction(Country, Class), fill = Class),
               stat = "identity", width = 0.06, alpha = 0.85, orientation = "y", color = "grey15") +
  geom_jitter(data = lc_row_long, aes(x = Log10_Burden, y = y_off - 0.11, color = Class),
              width = 0, height = 0.03, size = 0.8, alpha = 0.5) +
  scale_y_continuous(breaks = seq_along(levels(lc_row_long$Country)),
                      labels = levels(lc_row_long$Country)) +
  scale_x_continuous(labels = function(x) scales::scientific(10^x)) +
  scale_fill_manual(values  = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_color_manual(values = c("Indoor" = col_in, "Outdoor" = col_out)) +
  labs(
    title    = "Total SVOC Chemical Burden by Country (LC-HRMS)",
    subtitle = "Indoor and Outdoor shown together per country, ordered by fold change",
    x        = expression(Summed~peak~area~(log[10]~scale)),
    y        = NULL,
    fill     = NULL, color = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "top",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle      = element_text(size = 9, color = "grey30", hjust = 0),
    axis.text.y        = element_text(face = "bold", size = 11),
    panel.grid.major.y = element_blank()
  )

print(p_horizontal)
cat("  Saved: Fig_Chemical_Burden_ByCountry_Horizontal.png\n\n")

cat("\n==============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================\n")
cat(paste("Output directory:", output_dir, "\n\n"))
cat("Figures:\n")
cat("  Fig_Chemical_Burden_Overall.png              — Overall Indoor vs Outdoor\n")
cat("  Fig_Chemical_Burden_ByCountry.png            — By country with paired lines\n")
cat("  Fig_ChemBurden_FoldChange_ByCountry.png      — Fold change distribution per country\n")
cat("  Fig_Chemical_Burden_ByCountry_Horizontal.png — Horizontal, indoor+outdoor combined per row\n\n")
cat("Tables:\n")
cat("  TotalBurden_Statistics_Overall.csv          — Overall test result\n")
cat("  TotalBurden_Statistics_ByCountry.csv        — Per-country statistics\n")
cat("  TotalBurden_PairedBurden_AllHouseholds.csv  — Raw paired burden values\n")
