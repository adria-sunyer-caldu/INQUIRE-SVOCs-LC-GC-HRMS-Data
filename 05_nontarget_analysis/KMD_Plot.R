# ============================================================
# KENDRICK MASS DEFECT PLOT — LC-HRMS Nontarget Features
# Colored by Mean Paired Indoor/Outdoor Ratio
#
# Input:  NAP_paired_ratio_all_features.csv
# Output: Fig_KMD_All_IndoorRatio.png
#         Fig_KMD_StrongIndoor_vs_Outdoor.png
# ============================================================

library(data.table)
library(ggplot2)
library(scales)

# ── CONFIGURATION ─────────────────────────────────────────────
data_dir   <- "path/to/data"
input_file <- file.path(data_dir, "NAP_paired_ratio_all_features.csv")
out_dir    <- data_dir

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ── LOAD AND COMPUTE KMD ──────────────────────────────────────
cat("Loading data...\n")
df <- fread(input_file)
df <- df[!is.na(MZ) & !is.na(Mean_Paired_Ratio) & MZ > 50]
cat(paste("Features:", nrow(df), "\n"))

# Kendrick Mass Defect (CH2 base)
CH2_exact   <- 14.01565
CH2_nominal <- 14.0

df[, Kendrick_Mass := MZ * (CH2_nominal / CH2_exact)]
df[, Nominal_KM    := round(Kendrick_Mass)]
df[, KMD           := Nominal_KM - Kendrick_Mass]
df[, Category      := factor(Category,
    levels = c("Strong_Indoor","Indoor","Neutral","Outdoor","Strong_Outdoor"))]

cat("KMD range:", round(min(df$KMD), 3), "to", round(max(df$KMD), 3), "\n")

# ── COLOR SCALE ───────────────────────────────────────────────
# Blue (outdoor) → white (neutral) → red (indoor)
kmd_colors <- c(
  "0"   = "#1F618D",   # outdoor
  "0.5" = "#F5F5F5",   # neutral
  "1"   = "#B03A2E"    # indoor
)

# ── PLOT 1: ALL FEATURES ──────────────────────────────────────
cat("Plotting all features (this may take 30-60 seconds)...\n")

p1 <- ggplot(df, aes(x = MZ, y = KMD, color = Mean_Paired_Ratio)) +
  geom_point(size = 0.4, alpha = 0.35, stroke = 0, shape = 16) +
  scale_color_gradientn(
    colors = c("#1F618D", "#AED6F1", "#F5F5F5", "#F1948A", "#B03A2E"),
    values = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    name   = "Mean Paired Indoor Ratio\n(0 = outdoor  |  1 = indoor)",
    guide  = guide_colorbar(barwidth = 12, barheight = 0.8, title.position = "top")
  ) +
  scale_x_continuous(limits = c(50, quantile(df$MZ, 0.995)),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-0.55, 0.55),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(
    title    = "Kendrick Mass Defect Plot — LC-HRMS Nontarget Features (n = 155,448)",
    subtitle = "Colored by mean paired indoor/outdoor ratio across 205 households, 8 European countries",
    x        = expression(italic(m/z)),
    y        = expression("Kendrick Mass Defect (CH"[2]*")"),
    caption  = "Homologous series appear as horizontal rows of points at the same KMD value.\nRed = strongly indoor-dominant | Blue = strongly outdoor-dominant"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position   = "bottom",
    legend.title      = element_text(size = 9, hjust = 0.5),
    legend.text       = element_text(size = 8),
    plot.title        = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle     = element_text(size = 9, color = "grey40", hjust = 0),
    plot.caption      = element_text(size = 8, color = "grey50"),
    panel.grid.major  = element_line(color = "grey92", linewidth = 0.3),
    axis.text         = element_text(size = 10)
  )
print(p1)
ggsave(file.path(out_dir, "Fig_KMD_All_IndoorRatio.png"),
       p1, width = 12, height = 7, dpi = 300)
cat("Saved: Fig_KMD_All_IndoorRatio\n")

# ── PLOT 2: STRONG INDOOR vs STRONG OUTDOOR ───────────────────
cat("Plotting strong indoor vs outdoor...\n")

strong <- df[Category %in% c("Strong_Indoor", "Strong_Outdoor")]
strong[, Label := paste0(
  ifelse(Category == "Strong_Indoor",
         "Strong Indoor (ratio > 0.85)", "Strong Outdoor (ratio < 0.15)"),
  "  (n = ", format(.N, big.mark = ","), ")"),
  by = Category]
strong[, Category := factor(Category,
    levels = c("Strong_Indoor","Strong_Outdoor"))]

col_in  <- "#B03A2E"
col_out <- "#1F618D"

p2 <- ggplot(strong, aes(x = MZ, y = KMD,
                          color = Category)) +
  geom_point(size = 1.2, alpha = 0.5, stroke = 0, shape = 16) +
  scale_color_manual(
    values = c("Strong_Indoor"   = col_in,
               "Strong_Outdoor"  = col_out),
    labels = c("Strong_Indoor"   = paste0("Strong Indoor (ratio > 0.85)\nn = ",
                                    format(nrow(strong[Category=="Strong_Indoor"]), big.mark=",")),
               "Strong_Outdoor"  = paste0("Strong Outdoor (ratio < 0.15)\nn = ",
                                    format(nrow(strong[Category=="Strong_Outdoor"]), big.mark=",")))
  ) +
  facet_wrap(~ Category, nrow = 1, labeller = labeller(Category = c(
    Strong_Indoor  = paste0("Strong Indoor (ratio > 0.85)\nn = ",
                            format(nrow(strong[Category=="Strong_Indoor"]), big.mark=",")),
    Strong_Outdoor = paste0("Strong Outdoor (ratio < 0.15)\nn = ",
                            format(nrow(strong[Category=="Strong_Outdoor"]), big.mark=","))
  ))) +
  scale_x_continuous(limits = c(50, quantile(df$MZ, 0.995)),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(limits = c(-0.55, 0.55),
                     expand = expansion(mult = c(0.01, 0.01))) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4, linetype = "dashed") +
  labs(
    title   = "Kendrick Mass Defect: Strong Indoor vs Strong Outdoor Features",
    subtitle= "Homologous series appear as horizontal clusters at the same KMD value",
    x       = expression(italic(m/z)),
    y       = expression("Kendrick Mass Defect (CH"[2]*")"),
    caption = "Each point = one LC-HRMS nontarget feature. Horizontal rows indicate homologous series."
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(face = "bold", size = 13, hjust = 0),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    plot.caption     = element_text(size = 8, color = "grey50"),
    strip.background = element_rect(fill = "grey96", color = "grey70"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    axis.text        = element_text(size = 10)
  )
print(p2)
ggsave(file.path(out_dir, "Fig_KMD_StrongIndoor_vs_Outdoor.png"),
       p2, width = 12, height = 6, dpi = 600)
cat("Saved: Fig_KMD_StrongIndoor_vs_Outdoor\n")

# ── EXPORT: Figure 7a source data (Strong Indoor/Outdoor, KMD vs m/z) ──
fig7a_export <- strong[, .(Alignment_ID, MZ, Kendrick_Mass, Nominal_KM, KMD, Category)]
fwrite(fig7a_export, file.path(out_dir, "Fig7a_KMD_StrongIndoorOutdoor_SourceData.csv"))
cat(paste("Saved: Fig7a_KMD_StrongIndoorOutdoor_SourceData.csv (n =", nrow(fig7a_export), ")\n"))

cat("\nDone. Both plots saved to:", out_dir, "\n")
