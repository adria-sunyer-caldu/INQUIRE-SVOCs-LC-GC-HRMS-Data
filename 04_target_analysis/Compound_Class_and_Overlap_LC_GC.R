# ============================================================
# COMPOUND CLASS INDOOR/OUTDOOR ANALYSIS + VENN DIAGRAM
# INQUIRE Project | LC-HRMS + GC-HRMS
#
# Analyses:
#   1. Assign compound classes to LC targets+annotations and GC L1-3
#   2. Compute per-compound and per-class indoor/outdoor metrics
#   3. Statistical test per class (one-sample Wilcoxon, H0: median ratio = 0.5)
#   4. Figures: compound class barplot (LC), compound class barplot (GC),
#               combined dot/lollipop plot, Venn diagram
# ============================================================

library(readxl)
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(scales)
library(patchwork)

# ============================================================
# CONFIGURATION
# ============================================================

dir_path  <- "path/to/data"
lc_file   <- file.path(dir_path, "LC_Level1_2_polished_FINAL.xlsx")
gc_file   <- file.path(dir_path, "GC_Level1_2_polished_FINAL.xlsx")
out_dir   <- file.path(dir_path, "CompoundClass_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

col_in  <- "#B03A2E"
col_out <- "#1F618D"

# ============================================================
# HELPER: COMPOUND CLASS ASSIGNMENT — LC
# ============================================================

assign_class_lc <- function(name) {
  n <- tolower(name)
  
  # OPFRs
  if (grepl("phosphate|phosphoric|tris\\(|triethylphosph|diphenyl phosph|bis\\(2-ethylhexyl\\) phos|bis\\(1,3-dichloro|tris\\(chloropropyl\\)|tris\\(2-but|tris\\(2-chloro|tris\\(2-ethyl|tri-n-butyl phos|tri-o-cresyl|2-ethylhexyl diphenyl", n) &&
      !grepl("pfas|perfluoro|glycero|phosphocholine|phospho-|phospha|phosph.*acid.*fatty|inositol", n))
    return("OPFRs")
  
  # Phthalates
  if (grepl("phthalate|phthalic", n))
    return("Phthalates")
  
  # PFAS
  if (grepl("perfluoro|pfas|pfoa|pfos|pfbs|pfhx|fosaa|fbsa|trifluoroacet|trifluoromethane", n))
    return("PFAS")
  
  # PEGs / Glycols / Ethers
  if (grepl("ethylene glycol|propylene glycol|tetraglyme|triglyme|diglyme|monoglyme|methoxyethoxy|ethoxyethoxy|butoxyethoxy|glycol dim|glycol mono|polyethylene|oxatetrac|3,6,9,12-tetra", n))
    return("PEGs / Glycols")
  
  # Pesticides & biocides
  if (grepl("atrazin|terbutryn|terbutylazin|metolachlor|cyanazin|prometon|fenuron|pendimethalin|metalaxyl|carbofuran|carbendazim|malathion|imidacloprid|thiamethoxam|chlorantraniliprole|chlorothalonil|fluopyram|pyrimethanil|prosulfocarb|fipronil|bioallethrin|bentazon|2,4-dichlorophenoxy|3,5,6-trichloro-2-pyridinol|3-phenoxybenzoic|piperonyl|deet|icaridin|triclocarban|dcoit|pentachlorophenol|2-methyl-4-isothiazol", n))
    return("Pesticides / Biocides")
  
  # UV filters & cosmetics
  if (grepl("oxybenzone|octocrylene|benzophenone|methylbenzylidene camphor|ethyl 4-\\(dimethylamino\\)benz|crotamiton|paraben|methyl 4-hydroxy|tert-butyl 4-hydroxy|propylparaben", n))
    return("UV filters / Personal care")
  
  # Pharmaceuticals & drugs
  if (grepl("atenolol|carbamazepine|cotinine|hydroxycotinine|diclofenac|fluconazole|gabapentin|gemfibrozil|levetiracetam|metformin|o-desmethylvenlafaxine|paracetamol|pseudoephedrine|salbutamol|sulfadiazine|sulfamethoxazole|topiramate|trimethoprim|benzoylecgonine|nornicotine|myosmine|nicotine|caffeine|theobromine|anisomycin", n))
    return("Pharmaceuticals")
  
  # Fragrances / Musks
  if (grepl("galaxolide|linalool|nootkatone|sclareolide|myrtenal|cashmeran|camphor|ambrox|galaxolide|muskon", n))
    return("Fragrances / Musks")
  
  # Industrial chemicals / rubber
  if (grepl("diphenylguanidine|di-o-tolylguanidine|diphenyl sulfone|diphenylsulfoxide|hexamethylphosphoramide|isophorone|n,n-dimethyl-n'-p-tolyl|diphenylguanidine|diphenyl sulf", n))
    return("Industrial chemicals")
  
  # Surfactants / detergents
  if (grepl("lauryl|dodecylbenzene|diethanolamine|polidocanol|myristyl sulfate|lauroylethanolamine|diethanolamide|dodecanoyl|trimethyltetradecyl|cetrimonium", n))
    return("Surfactants / Detergents")
  
  # Lipids / fatty acids
  if (grepl("phosphocholine|glycerophosph|phospho-.*inositol|palmitoyl|hexadecanoyl|dilinoleoyl|galcer|cholesterol|palmitelaidic|hydroxydodecanoate|hydroxydecanoic|hydroxyoctadec|hydroxyhexa|brassylic|decanedioic|tg\\(|dg\\(|glycan", n))
    return("Lipids / Fatty acids")
  
  # Biogenic / natural
  if (grepl("adenine|adenosine|uridine|tocopherol|daidzein|hesperidin|kaempferol|genkwanin|sakuranetin|diosmin|myrcene|usnic|terrein|rupestonic|bisabolol|dihydroactinidiolide|loliolide|dimethylfraxetin|dehydrocostus|petasol|4-coumaric|p-coumaric|cinnamic|gentisyl|proline|tyrosine|pipecolic|nicotinamide|pyridoxine|histamine|trehalose|melezitose|urocanic|sucralose|acesulfam|sorbic acid|melamine|hydroquinone|pyrogallol|adenosin|theobromin|caffeine|cotinin", n))
    return("Biogenic / Natural")
  
  return("Other / Unclassified")
}

# ============================================================
# HELPER: COMPOUND CLASS ASSIGNMENT — GC
# ============================================================

assign_class_gc <- function(name) {
  n <- tolower(name)
  
  if (grepl("phthalate|phthalic", n))
    return("Phthalates")
  
  if (grepl("fluoranthene|pyrene|dibenzofuran|dibenzothiophen|fluorene|acenaphth|tetramethyl-biphenyl|diisopropylnaphthalene", n))
    return("PAHs (4+ rings)")
  
  if (grepl("naphthalene|methylnaphthalene|biphenyl|phenanthrene|anthracene", n))
    return("PAHs (2-3 rings)")
  
  if (grepl("toluene|xylene|trimethylbenzene|diisopropylnaphthalene|benzene, 1|benzene, \\(1-|benzene, \\(2-|ethylbenzene", n))
    return("BTEX / Alkylbenzenes")
  
  if (grepl("hexachlorobenzene|pentachlorophenol|chloroxylenol|bromophenyl|2-bromo|difluoro|trichloromethyl|tetramethyl-biphenyl|4-bromophenyl", n))
    return("Halogenated compounds")
  
  if (grepl("myrcene|pinene|limonene|terpinene|camphor|borneol|menthol|terpineol|carvone|geranyl|bornyl|citronellol|cedrene|ionone|cashmeran|nerolidol|ambrox|guaiol|liguloxide|amberonne|terpinyl|hexylcinnamal|sclareolide|linalool|bisabolol|nootkatone", n))
    return("Terpenes / Fragrance")
  
  if (grepl("salicylate|benzyl benzoate|phenoxyethanol|phenoxypropan|diphenyl ether|lilial|cyclohexyl sal|hexyl sal|2-ethylhexyl sal|n-hexyl sal|benzyl sal", n))
    return("UV filters / Fragrance esters")
  
  if (grepl("diisopropyl adipate|hexanoic acid, hexyl|adipic acid|bornyl acetate|geranyl acetate|isoamyl|dimethyl glutarate|tetrahydrofurfuryl|isopropyl myristate|isopropyl palmitate|methyl dehydroabietate|benzyl benzoate|tributyl|acetylcitrate|succinic acid|glutaric acid|oxalic acid|isophthalic|terephthalic|benzoic acid|butanedioic|dodecanoic acid, methyl|hexanoic acid|gamma-decalactone|gamma-undecalactone|delta-decalacton|dibutyl succ|dibutyl adipate|ethyl.*hexanoate|cyclopenta.*methyl ester|benzoate", n))
    return("Esters / Plasticizers")
  
  if (grepl("^nonane$|^undecane$|tridecane|tetradecane|hexadecane|heptadecane|eicosane|docosane|nonadecane|heneicosane|pentane, 2,2,4|octane, 2,4,6|butane, 2,2-|hexane, 3,3|dimethyl octane|heptamethyl|undecane, |octane, 3,3|nonane, 2-methyl|benzene, \\(1-butyl|benzene, \\(1-ethyl|benzene, \\(1-methyl", n))
    return("Alkanes / Waxes")
  
  if (grepl("aldehyde|decanal|dodecanal|benzaldehyde|nonanal|pelargonaldehyde|decylaldehyde|octanal|hexanal", n))
    return("Aldehydes")
  
  if (grepl("2-ethylhexanol|3-methoxy-3-methylbutanol|ethyl tert-butyl ether|propanol|butanol|borneol|terpineol|menthol|dodecenol|glucitol|methyl.*diol|hydroxymethyl|2-phenoxyethanol|phenoxypropan|1,4,7-trimethyl.*diol", n))
    return("Alcohols / Ethers")
  
  return("Other / Unclassified")
}

# ============================================================
# STEP 1: READ AND PROCESS LC DATA
# ============================================================

cat("Reading LC data...\n")
lc_raw <- as.data.table(read_excel(lc_file))

# Sample columns (pseudonymized codes): {Country}_{Indoor|Outdoor}_{code}, e.g. NL_Indoor_AA / NL_Outdoor_AA
lc_sample_cols <- names(lc_raw)[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?\\s*$", names(lc_raw))]
lc_in_cols  <- lc_sample_cols[grepl("_Indoor_", lc_sample_cols)]
lc_out_cols <- lc_sample_cols[grepl("_Outdoor_", lc_sample_cols)]

cat(paste("  LC features:", nrow(lc_raw), "| Indoor samples:", length(lc_in_cols),
          "| Outdoor samples:", length(lc_out_cols), "\n"))

# Convert to numeric
for (col in c(lc_in_cols, lc_out_cols)) {
  lc_raw[, (col) := as.numeric(get(col))]
  lc_raw[is.na(get(col)), (col) := 0]
}

# Compute per-compound means
lc_raw[, mean_indoor  := rowMeans(.SD, na.rm = TRUE), .SDcols = lc_in_cols]
lc_raw[, mean_outdoor := rowMeans(.SD, na.rm = TRUE), .SDcols = lc_out_cols]
lc_raw[, indoor_ratio := mean_indoor / (mean_indoor + mean_outdoor + 1e-9)]
lc_raw[, log2FC       := log2((mean_indoor + 1) / (mean_outdoor + 1))]
lc_raw[, Class        := sapply(Name, assign_class_lc)]
lc_raw[, Platform     := "LC-HRMS (ESI+/-)"]

lc_class_data <- lc_raw[, .(Name, Class, Platform, indoor_ratio, log2FC,
                            mean_indoor, mean_outdoor,
                            Confidence = `Confidence level`)]

cat("LC class distribution:\n")
print(lc_class_data[, .N, by = Class][order(-N)])

# ============================================================
# STEP 2: READ AND PROCESS GC DATA
# ============================================================

cat("\nReading GC data...\n")
gc_raw <- as.data.table(read_excel(gc_file))
gc_raw <- gc_raw[!is.na(as.numeric(RT))]

# Keep only annotated features (L1-3)
gc_ann <- gc_raw[`ID level` %in% c(1, 2, 3)]

gc_in_cols  <- names(gc_ann)[grepl("^[A-Z]{2}_Indoor_[A-Z]{2}$", names(gc_ann))]
gc_out_cols <- names(gc_ann)[grepl("^[A-Z]{2}_Outdoor_[A-Z]{2}$", names(gc_ann))]

# Exclude 4CPS_PDMS (QC)
cat(paste("  GC annotated features (L1-3):", nrow(gc_ann),
          "| Indoor:", length(gc_in_cols),
          "| Outdoor:", length(gc_out_cols), "\n"))

for (col in c(gc_in_cols, gc_out_cols)) {
  gc_ann[, (col) := as.numeric(get(col))]
  gc_ann[is.na(get(col)), (col) := 0]
}

gc_ann[, mean_indoor  := rowMeans(.SD, na.rm = TRUE), .SDcols = gc_in_cols]
gc_ann[, mean_outdoor := rowMeans(.SD, na.rm = TRUE), .SDcols = gc_out_cols]
gc_ann[, indoor_ratio := mean_indoor / (mean_indoor + mean_outdoor + 1e-9)]
gc_ann[, log2FC       := log2((mean_indoor + 1) / (mean_outdoor + 1))]
gc_ann[, Class        := sapply(`Database match`, assign_class_gc)]
gc_ann[, Platform     := "GC-HRMS (EI)"]
gc_ann[, Name         := `Database match`]
gc_ann[, Confidence   := as.character(`ID level`)]

gc_class_data <- gc_ann[, .(Name, Class, Platform, indoor_ratio, log2FC,
                            mean_indoor, mean_outdoor, Confidence)]

cat("GC class distribution:\n")
print(gc_class_data[, .N, by = Class][order(-N)])

# ============================================================
# STEP 3: STATISTICAL ANALYSIS PER CLASS
# ============================================================

cat("\n=== STATISTICAL ANALYSIS PER CLASS ===\n")

run_class_stats <- function(data, platform_name) {
  classes <- unique(data$Class)
  results <- rbindlist(lapply(classes, function(cls) {
    sub <- data[Class == cls]
    n   <- nrow(sub)
    if (n < 3) return(NULL)
    
    # One-sample Wilcoxon: is median indoor ratio != 0.5?
    wt <- tryCatch(
      wilcox.test(sub$indoor_ratio, mu = 0.5, alternative = "two.sided", exact = FALSE),
      error = function(e) list(statistic = NA, p.value = NA)
    )
    
    data.table(
      Platform         = platform_name,
      Class            = cls,
      N                = n,
      Median_IndoorRatio = round(median(sub$indoor_ratio), 3),
      Mean_IndoorRatio = round(mean(sub$indoor_ratio), 3),
      Median_Log2FC    = round(median(sub$log2FC), 2),
      Median_FC        = round(2^median(sub$log2FC), 2),
      Pct_IndoorDom    = round(mean(sub$indoor_ratio > 0.5) * 100, 1),
      W_statistic      = wt$statistic,
      P_value          = wt$p.value
    )
  }))
  
  results[, P_adj := p.adjust(P_value, method = "BH")]
  results[, Sig   := fcase(
    P_adj < 0.001, "***",
    P_adj < 0.01,  "**",
    P_adj < 0.05,  "*",
    default        = "ns"
  )]
  results[order(-Median_IndoorRatio)]
}

lc_stats <- run_class_stats(lc_class_data, "LC-HRMS (ESI+/-)")
gc_stats <- run_class_stats(gc_class_data, "GC-HRMS (EI)")

cat("\nLC compound class statistics:\n")
print(lc_stats[, .(Class, N, Median_FC, Pct_IndoorDom,
                   P_value = formatC(P_value, format = "e", digits = 2),
                   P_adj = formatC(P_adj, format = "e", digits = 2), Sig)])

cat("\nGC compound class statistics:\n")
print(gc_stats[, .(Class, N, Median_FC, Pct_IndoorDom,
                   P_value = formatC(P_value, format = "e", digits = 2),
                   P_adj = formatC(P_adj, format = "e", digits = 2), Sig)])

# Save stats
fwrite(lc_stats, file.path(out_dir, "LC_CompoundClass_Statistics.csv"))
fwrite(gc_stats, file.path(out_dir, "GC_CompoundClass_Statistics.csv"))

# ============================================================
# STEP 4: FIGURE — LC COMPOUND CLASS BARPLOT
# ============================================================

cat("\nGenerating figures...\n")

make_class_plot <- function(stats_dt, title_str) {
  stats_dt <- stats_dt[N >= 3 & Class != "Other / Unclassified"]

  # Convert 0-1 indoor ratio (0.5 = equal) to the paper's standard -1/+1
  # scale (-1 = fully outdoor, 0 = equal, +1 = fully indoor).
  stats_dt[, Ratio_pm1 := 2 * (Mean_IndoorRatio - 0.5)]
  stats_dt[, Class := factor(Class, levels = Class[order(Ratio_pm1)])]
  stats_dt[, point_color := ifelse(Ratio_pm1 > 0.2, col_in,
                                   ifelse(Ratio_pm1 > 0, "#7D6608", col_out))]

  offset    <- 0.03   # small, fixed gap between the dot and its label — same for every row
  edge_pad  <- 0.12   # extra scale room beyond +/-1, reserved for label text only

  ggplot(stats_dt, aes(x = Class, y = Ratio_pm1)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "grey50", linewidth = 0.7) +
    geom_segment(aes(x = Class, xend = Class,
                     y = 0, yend = Ratio_pm1,
                     color = point_color),
                 linewidth = 3, alpha = 0.7) +
    geom_point(aes(color = point_color), size = 7, alpha = 0.9) +
    geom_text(aes(y = Ratio_pm1 + ifelse(Ratio_pm1 >= 0, offset, -offset),
                  hjust = ifelse(Ratio_pm1 >= 0, 0, 1),
                  label = paste0(Sig, "  n=", N)),
              size = 4, color = "grey20") +
    scale_color_identity() +
    scale_y_continuous(
      limits = c(-1 - edge_pad, 1 + edge_pad),
      breaks = seq(-1, 1, by = 0.5),
      labels = scales::number_format(accuracy = 0.1),
      expand = expansion(mult = c(0, 0))
    ) +
    coord_flip(clip = "off") +
    labs(
      title   = title_str,
      x       = NULL,
      y       = "Mean Indoor Ratio (\u22121 = fully outdoor | +1 = fully indoor)"
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position    = "none",
      plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
      panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
      axis.text.y        = element_text(size = 9, color = "#1F4E5F"),
      axis.text.x        = element_text(size = 9),
      plot.margin        = margin(t = 5.5, r = 5.5, b = 5.5, l = 5.5)
    )
}

p_lc <- make_class_plot(lc_stats, "LC-HRMS")

p_gc <- make_class_plot(gc_stats, "GC-HRMS")

print (p_lc)
print(p_gc)
ggsave(file.path(out_dir, "Fig_LC_CompoundClass.png"), p_lc, width = 12, height = 9, dpi = 600)
ggsave(file.path(out_dir, "Fig_GC_CompoundClass.png"), p_gc, width = 12, height = 8, dpi = 600)
cat("  Saved: Fig_LC_CompoundClass and Fig_GC_CompoundClass\n")



# ============================================================
# STEP 6: VENN DIAGRAM — manual ggplot2
# ============================================================

lc_keys <- unique(na.omit(toupper(trimws(lc_raw$InChiKey))))
gc_keys <- unique(na.omit(toupper(trimws(gc_ann$InChiKey))))
overlap  <- intersect(lc_keys, gc_keys)
shared_names <- lc_raw[toupper(trimws(InChiKey)) %in% overlap, Name]

n_lc   <- length(lc_keys)
n_gc   <- length(gc_keys)
n_both <- length(overlap)

shared_label <- paste0("Shared (n=", n_both, "):\n",
                       paste(shared_names, collapse = "\n"))

p_venn <- ggplot() +
  # LC circle
  annotate("path",
           x = 0.38 + 0.32 * cos(seq(0, 2*pi, length.out = 300)),
           y = 0.5  + 0.38 * sin(seq(0, 2*pi, length.out = 300)),
           color = "#2471A3", linewidth = 1.8, alpha = 0.9) +
  annotate("polygon",
           x = 0.38 + 0.32 * cos(seq(0, 2*pi, length.out = 300)),
           y = 0.5  + 0.38 * sin(seq(0, 2*pi, length.out = 300)),
           fill = "#2471A3", alpha = 0.15) +
  # GC circle
  annotate("path",
           x = 0.62 + 0.22 * cos(seq(0, 2*pi, length.out = 300)),
           y = 0.5  + 0.26 * sin(seq(0, 2*pi, length.out = 300)),
           color = "#BA4A00", linewidth = 1.8, alpha = 0.9) +
  annotate("polygon",
           x = 0.62 + 0.22 * cos(seq(0, 2*pi, length.out = 300)),
           y = 0.5  + 0.26 * sin(seq(0, 2*pi, length.out = 300)),
           fill = "#BA4A00", alpha = 0.15) +
  # LC label inside
  annotate("text", x = 0.28, y = 0.5,
           label = paste0("LC-HRMS\n(ESI+/-)\nn = ", n_lc - n_both),
           size = 4.2, fontface = "bold", color = "#2471A3", hjust = 0.5) +
  # GC label inside
  annotate("text", x = 0.74, y = 0.5,
           label = paste0("GC-HRMS\n(EI)\nn = ", n_gc - n_both),
           size = 3.8, fontface = "bold", color = "#BA4A00", hjust = 0.5) +
  # Overlap label
  annotate("text", x = 0.535, y = 0.5,
           label = as.character(n_both),
           size = 4.5, fontface = "bold", color = "grey20", hjust = 0.5) +
  # Shared compound names below
  annotate("text", x = 0.535, y = 0.08,
           label = shared_label,
           size = 2.8, color = "grey30", fontface = "italic", hjust = 0.5,
           lineheight = 1.3) +
  # Platform labels outside
  annotate("text", x = 0.08, y = 0.92,
           label = "LC-HRMS\n(ESI+/-)",
           size = 4.0, fontface = "bold", color = "#2471A3", hjust = 0) +
  annotate("text", x = 0.78, y = 0.92,
           label = "GC-HRMS\n(EI)",
           size = 4.0, fontface = "bold", color = "#BA4A00", hjust = 0) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title    = "Annotated Compound Coverage: LC-HRMS vs GC-HRMS",
    subtitle = paste0("LC: ", nrow(lc_raw), " targets + annotations  |  ",
                      "GC: ", nrow(gc_ann), " level 1\u20133 features  |  ",
                      "Overlap: ", n_both, " compounds (<0.5%)")
  ) +
  theme_void() +
  theme(
    plot.title    = element_text(face = "bold", size = 12,
                                 hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 9, color = "grey30",
                                 hjust = 0.5, margin = margin(b = 10)),
    plot.margin   = margin(10, 10, 10, 10)
  )
print(p_venn)
ggsave(file.path(out_dir, "Fig_Venn_LC_GC.png"), p_venn,
       width = 7, height = 6, dpi = 300)
cat("  Saved: Fig_Venn_LC_GC.png\n")

# ============================================================
# DONE
# ============================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat(paste("Output:", out_dir, "\n"))
cat("Files:\n")
cat("  LC_CompoundClass_Statistics.csv\n")
cat("  GC_CompoundClass_Statistics.csv\n")
cat("  Fig_LC_CompoundClass.png\n")
cat("  Fig_GC_CompoundClass.png\n")
cat("  Fig_Venn_LC_GC.png\n")