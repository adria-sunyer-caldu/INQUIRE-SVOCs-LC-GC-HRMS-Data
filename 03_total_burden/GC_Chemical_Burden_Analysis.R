# ============================================================
# TOTAL CHEMICAL BURDEN ANALYSIS — GC-HRMS
# INQUIRE Project | GC-HRMS Nontarget Data (NILU)
#
# Parallel script to the LC-HRMS total chemical burden analysis.
# Adapted for GC-HRMS file structure:
#   - Excel input (no MS-DIAL header rows)
#   - Sample naming (pseudonymized): {Country}_Indoor_{code} / {Country}_Outdoor_{code}
#     (indoor and outdoor samples of the same household share the same code;
#     extra indoor samplers with a _2/_3 suffix are not used)
#   - 4CPS_PDMS columns = QC samples (excluded)
#   - Italy indoor/outdoor labelling already corrected in this file
#
# Approach: Sum of all feature areas per sample, as a proxy for
# total chemical burden (same approach as the LC-HRMS script).
# ============================================================

library(readxl)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)
library(scales)

# ============================================================
# CONFIGURATION
# ============================================================

data_dir <- "path/to/data"
file_path  <- file.path(data_dir, "GC_Level1_5_polished_FULL_zenodo.xlsx")   # Zenodo GC deposit
output_dir <- file.path(data_dir, "Total_Chemical_Burden_Output_GC")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

# ============================================================
# STEP 1: READ DATA
# ============================================================

cat("Step 1: Reading GC-HRMS Excel file...\n")

raw <- as.data.table(read_excel(file_path))

# Drop rows with no RT (empty/header rows)
raw <- raw[!is.na(as.numeric(RT))]

cat(paste("  Features loaded:", nrow(raw), "\n"))

# Identify sample columns
all_cols <- names(raw)

# Indoor / Outdoor main samplers only; QC columns (4CPS_PDMS) not matched
is_cols  <- all_cols[grepl("^[A-Z]{2}_Indoor_[A-Z]{2}$", all_cols)]
os_cols  <- all_cols[grepl("^[A-Z]{2}_Outdoor_[A-Z]{2}$", all_cols)]
qc_cols  <- all_cols[grepl("^4CPS_PDMS", all_cols)]

cat(paste("  Indoor samples:        ", length(is_cols), "\n"))
cat(paste("  Outdoor samples:       ", length(os_cols), "\n"))
cat(paste("  QC samples excluded:   ", length(qc_cols), "\n\n"))

# ============================================================
# STEP 2: COMPUTE TOTAL CHEMICAL BURDEN PER SAMPLE
# ============================================================

cat("Step 2: Computing total chemical burden per sample...\n")

# Convert area columns to numeric, replace NA with 0
for (col in c(is_cols, os_cols)) {
  raw[, (col) := as.numeric(get(col))]
  raw[is.na(get(col)), (col) := 0]
}

# Sum all feature areas per sample
burden_indoor <- raw[, lapply(.SD, sum, na.rm = TRUE), .SDcols = is_cols]
burden_outdoor <- raw[, lapply(.SD, sum, na.rm = TRUE), .SDcols = os_cols]

# Long format
burden_long <- rbindlist(list(
  data.table(
    Sample  = is_cols,
    Burden  = as.numeric(burden_indoor[1, ]),
    Class   = "Indoor"
  ),
  data.table(
    Sample  = os_cols,
    Burden  = as.numeric(burden_outdoor[1, ]),
    Class   = "Outdoor"
  )
))

# Extract Country and household code from naming: {Country}_{Indoor|Outdoor}_{code}
burden_long[, Country := sub("_.*", "", Sample)]
burden_long[, House   := sub("^[A-Z]{2}_(Indoor|Outdoor)_", "", Sample)]

# Handle Slovenia abbreviation (SL in GC vs SI in LC)
burden_long[Country == "SL", Country := "SI"]

cat("  Burden computed. Summary:\n")
print(burden_long[, .(N = .N, Median = round(median(Burden))),
                  by = .(Country, Class)][order(Country, Class)])
cat("\n")

# ============================================================
# STEP 3: BUILD PAIRED DATASET
# ============================================================

cat("Step 3: Building paired dataset...\n")

burden_wide <- dcast(burden_long, Country + House ~ Class, value.var = "Burden",
                     fun.aggregate = mean)

n_before <- nrow(burden_wide)
burden_wide <- na.omit(burden_wide)
cat(paste("  Complete pairs:", nrow(burden_wide), "(removed", n_before - nrow(burden_wide), "incomplete)\n\n"))

burden_wide[, Log2FC        := log2(Indoor / Outdoor)]
burden_wide[, IndoorRatio   := Indoor / (Indoor + Outdoor)]
burden_wide[, Indoor_Higher := Indoor > Outdoor]

# ============================================================
# STEP 4: STATISTICAL ANALYSIS — OVERALL
# ============================================================

cat("==============================================================\n")
cat("STATISTICAL RESULTS — OVERALL (all countries combined)\n")
cat("==============================================================\n")

wtest <- wilcox.test(burden_wide$Indoor, burden_wide$Outdoor,
                     paired = TRUE, exact = FALSE, alternative = "two.sided")

n_total  <- nrow(burden_wide)
r_rb     <- (2 * wtest$statistic) / (n_total * (n_total + 1) / 2) - 1
med_fc   <- round(2^median(burden_wide$Log2FC), 2)
pct_high <- round(mean(burden_wide$Indoor_Higher) * 100, 1)

cat(paste("Test:                      Paired Wilcoxon signed-rank\n"))
cat(paste("N pairs:                  ", n_total, "\n"))
cat(paste("V statistic:              ", wtest$statistic, "\n"))
cat(paste("p-value:                  ", format(wtest$p.value, scientific = TRUE, digits = 3), "\n"))
cat(paste("Effect size (r_rb):       ", round(r_rb, 3), "\n"))
cat(paste("Median indoor burden:     ", format(round(median(burden_wide$Indoor)), big.mark = ","), "\n"))
cat(paste("Median outdoor burden:    ", format(round(median(burden_wide$Outdoor)), big.mark = ","), "\n"))
cat(paste("Median fold change:       ", med_fc, "x\n"))
cat(paste("Median log2 FC:           ", round(median(burden_wide$Log2FC), 3), "\n"))
cat(paste("% pairs Indoor > Outdoor: ", pct_high, "%\n\n"))

# ============================================================
# STEP 5: STATISTICAL ANALYSIS — BY COUNTRY
# ============================================================

cat("==============================================================\n")
cat("STATISTICAL RESULTS — BY COUNTRY\n")
cat("==============================================================\n")

countries <- sort(unique(burden_wide$Country))

country_stats <- rbindlist(lapply(countries, function(ctry) {
  sub <- burden_wide[Country == ctry]
  if (nrow(sub) < 5) return(NULL)
  wt  <- wilcox.test(sub$Indoor, sub$Outdoor, paired = TRUE,
                     exact = FALSE, alternative = "two.sided")
  r   <- (2 * wt$statistic) / (nrow(sub) * (nrow(sub) + 1) / 2) - 1
  data.table(
    Country           = ctry,
    N_pairs           = nrow(sub),
    Median_FC         = round(2^median(sub$Log2FC), 2),
    Pct_Indoor_Higher = round(mean(sub$Indoor_Higher) * 100, 1),
    W_statistic       = wt$statistic,
    P_value           = wt$p.value,
    Effect_size_r     = round(r, 3)
  )
}))

country_stats[, P_adj := p.adjust(P_value, method = "BH")]
country_stats[, Sig   := fcase(
  P_adj < 0.001, "***",
  P_adj < 0.01,  "**",
  P_adj < 0.05,  "*",
  default        = "ns"
)]

print(country_stats[, .(Country, N_pairs, Median_FC, Pct_Indoor_Higher,
                        P_value = formatC(P_value, format = "e", digits = 2),
                        P_adj   = formatC(P_adj,   format = "e", digits = 2),
                        Sig, Effect_size_r)])

# ============================================================
# STEP 6: SAVE RESULTS
# ============================================================

fwrite(country_stats,  "GC_TotalBurden_Statistics_ByCountry.csv")
fwrite(burden_wide,    "GC_TotalBurden_PairedBurden_AllHouseholds.csv")
fwrite(data.table(
  Test               = "Paired Wilcoxon signed-rank (all countries)",
  N_pairs            = n_total,
  V_statistic        = wtest$statistic,
  P_value            = wtest$p.value,
  Effect_size_r      = round(r_rb, 3),
  Median_Log2FC      = round(median(burden_wide$Log2FC), 3),
  Median_FC          = med_fc,
  Pct_Indoor_Higher  = pct_high
), "GC_TotalBurden_Statistics_Overall.csv")

# ============================================================
# STEP 7: PLOT 1 — OVERALL Indoor vs Outdoor
# ============================================================

cat("\nStep 7: Generating plots...\n")

col_in  <- "#B03A2E"
col_out <- "#1F618D"

plot_long <- burden_long[Class %in% c("Indoor", "Outdoor")]
plot_long[, Class        := factor(Class, levels = c("Outdoor", "Indoor"))]
plot_long[, Log10_Burden := log10(Burden + 1)]

p_val_str <- ifelse(
  wtest$p.value < 2.2e-16, "p < 2.2\u00D710\u207B\u00B9\u2076",
  paste0("p = ", formatC(wtest$p.value, format = "e", digits = 2))
)
subtitle_overall <- paste0(
  "n = ", n_total, " paired households | ",
  med_fc, "\u00D7 higher indoors (median) | ", p_val_str
)

p1 <- ggplot(plot_long, aes(x = Class, y = Log10_Burden, fill = Class, color = Class)) +
  geom_violin(alpha = 0.25, trim = FALSE, linewidth = 0.7, adjust = 1.2) +
  geom_boxplot(width = 0.12, alpha = 0.9, outlier.shape = NA,
               linewidth = 0.6, color = "grey15") +
  geom_beeswarm(size = 1.4, alpha = 0.55, cex = 1.6, priority = "random") +
  scale_fill_manual(values  = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_color_manual(values = c("Indoor" = col_in, "Outdoor" = col_out)) +
  scale_y_continuous(
    labels = function(x) scales::scientific(10^x),
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    title    = "Total SVOC Chemical Burden (GC-HRMS)",
    subtitle = subtitle_overall,
    x        = NULL,
    y        = "Summed peak area (IS-normalized, blank-subtracted)",
    caption  = "Each point = one air sample. Box = IQR; whiskers = 1.5\u00D7IQR.\nData: GC-HRMS EI, normalized, blank-subtracted. Source: NILU."
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 14, hjust = 0),
    plot.subtitle      = element_text(size = 10, color = "grey30", hjust = 0),
    plot.caption       = element_text(size = 8,  color = "grey50"),
    axis.text.x        = element_text(size = 13, face = "bold"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4)
  )

print(p1)
ggsave("GC_Fig_Chemical_Burden_Overall.png", p1, width = 4.5, height = 6, dpi = 600)
cat("  Saved: GC_Fig_Chemical_Burden_Overall.png\n")

# ============================================================
# STEP 8: PLOT 2 — BY COUNTRY (with paired lines)
# ============================================================

# Facet labels with FC and significance
label_dt <- country_stats[, .(
  Country,
  facet_label = paste0(Country, "  (", Median_FC, "\u00D7  ", Sig, ")")
)]

plot_long_c <- merge(plot_long, label_dt, by = "Country")
plot_long_c[, Class := factor(Class, levels = c("Outdoor", "Indoor"))]

# Paired lines
pairs_long <- melt(
  burden_wide[, .(Country, House,
                  Outdoor = log10(Outdoor + 1),
                  Indoor  = log10(Indoor  + 1))],
  id.vars       = c("Country", "House"),
  variable.name = "Class",
  value.name    = "Log10_Burden"
)
pairs_long[, Class := factor(Class, levels = c("Outdoor", "Indoor"))]
pairs_long <- merge(pairs_long, label_dt, by = "Country")

# Order by fold change descending
fc_order <- country_stats[order(-Median_FC), Country]
plot_long_c[, facet_label := factor(facet_label,
                                    levels = label_dt[match(fc_order, Country), facet_label])]
pairs_long[, facet_label := factor(facet_label,
                                   levels = label_dt[match(fc_order, Country), facet_label])]

p2 <- ggplot(plot_long_c, aes(x = Class, y = Log10_Burden,
                              fill = Class, color = Class)) +
  geom_line(data = pairs_long,
            aes(x = Class, y = Log10_Burden, group = House),
            color = "grey65", alpha = 0.30, linewidth = 0.35,
            inherit.aes = FALSE) +
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
    title    = "Total SVOC Chemical Burden by Country (GC-HRMS)",
    subtitle = "GC-HRMS EI | Ordered by indoor/outdoor fold change (FC)\nGrey lines connect paired samples from the same household",
    x        = NULL,
    y        = "Summed peak area (log\u2081\u2080 scale)",
    caption  = "FC = median fold change Indoor/Outdoor. * p<0.05, ** p<0.01, *** p<0.001 (BH-corrected paired Wilcoxon).\nSource: NILU GC-HRMS. Italy indoor/outdoor labelling corrected."
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position    = "none",
    plot.title         = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle      = element_text(size = 9, color = "grey30", hjust = 0),
    plot.caption       = element_text(size = 7.5, color = "grey50"),
    strip.background   = element_rect(fill = "grey96", color = "grey70"),
    strip.text         = element_text(face = "bold", size = 9),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.35),
    panel.spacing      = unit(0.9, "lines"),
    axis.text.x        = element_text(size = 10, face = "bold")
  )
print(p2)
ggsave("GC_Fig_Chemical_Burden_ByCountry.png", p2, width = 13, height = 7.5, dpi = 600)
cat("  Saved: GC_Fig_Chemical_Burden_ByCountry.png\n")

# ============================================================
# STEP 9: PLOT 3 — FOLD CHANGE BY COUNTRY
# ============================================================

fc_plot <- burden_wide[, .(Country, Log2FC)]
fc_plot <- merge(fc_plot, country_stats[, .(Country, Sig, Median_FC)], by = "Country")
fc_plot[, Country := factor(Country, levels = country_stats[order(Median_Log2FC = log2(Median_FC)), Country])]

p3 <- ggplot(fc_plot, aes(x = Country, y = Log2FC)) +
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
    data = country_stats[order(log2(Median_FC))],
    aes(x = Country,
        y = max(burden_wide$Log2FC, na.rm = TRUE) * 1.05,
        label = paste0(Median_FC, "\u00D7\n", Sig)),
    size = 3.0, fontface = "bold", color = "grey20",
    inherit.aes = FALSE
  ) +
  scale_y_continuous(
    name     = expression(log[2]~fold~change~(Indoor/Outdoor)),
    sec.axis = sec_axis(~ 2^.,
                        name   = "Fold change (Indoor/Outdoor)",
                        breaks = c(0.5, 1, 2, 4, 8, 16, 32),
                        labels = c("0.5\u00D7","1\u00D7","2\u00D7","4\u00D7","8\u00D7","16\u00D7","32\u00D7"))
  ) +
  labs(
    title    = "Indoor/Outdoor Fold Change in Total Chemical Burden (GC-HRMS)",
    subtitle = "Paired Wilcoxon signed-rank | BH-corrected | * p<0.05, ** p<0.01, *** p<0.001",
    x        = NULL,
    caption  = "Each violin/boxplot: log\u2082(Indoor/Outdoor) per household pair.\nDashed line = equal burden. Italy indoor/outdoor labelling corrected."
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
ggsave("GC_Fig_ChemBurden_FoldChange_ByCountry.png", p3, width = 10, height = 6, dpi = 600)
cat("  Saved: GC_Fig_ChemBurden_FoldChange_ByCountry.png\n")

# ============================================================
# DONE
# ============================================================

cat("\n==============================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================\n")
cat(paste("Output directory:", output_dir, "\n\n"))
cat("Figures:\n")
cat("  GC_Fig_Chemical_Burden_Overall.png\n")
cat("  GC_Fig_Chemical_Burden_ByCountry.png\n")
cat("  GC_Fig_ChemBurden_FoldChange_ByCountry.png\n\n")
cat("Tables:\n")
cat("  GC_TotalBurden_Statistics_Overall.csv\n")
cat("  GC_TotalBurden_Statistics_ByCountry.csv\n")
cat("  GC_TotalBurden_PairedBurden_AllHouseholds.csv\n")
