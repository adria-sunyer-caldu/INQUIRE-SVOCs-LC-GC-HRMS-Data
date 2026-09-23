library(readxl)
library(data.table)
library(ggplot2)

# ============================================================
# Single-panel widespread-compound figure: all
# compounds on one continuous x-axis, grouped by class (dashed
# dividers + class labels above), sorted within each class,
# Indoor (red) / Outdoor (blue) paired bars.
#
# Two metrics per platform, per user request:
#   - Detection frequency (%) - robust to outliers, prevalence story
#   - LC: mean concentration | GC: mean summed peak area (GC has no
#     true concentration units, so "area" is used instead, matching
#     the "summed peak area" terminology already used throughout
#     the GC burden figures elsewhere in this project)
#
# OPFR classifier bug fixed: the old exclusion term "phospha" is a
# literal substring of "phosphate" itself, so it was silently
# disqualifying every true OPFR compound whose name contains the
# word "phosphate" (verified against the real target list - all 13
# true OPFRs were being misclassified). Replaced with the actually-
# intended "phosphatidyl" term, verified against all real OPFR and
# lipid-phosphate compound names in this dataset.
# ============================================================

dir_path <- "path/to/data"
lc_file  <- file.path(dir_path, "LC_Level1_2_polished_FINAL.xlsx")
gc_file  <- file.path(dir_path, "GC_Level1_2_polished_FINAL.xlsx")
out_dir  <- dir_path

col_in  <- "#B03A2E"
col_out <- "#1F618D"

# ---- Explicit rename map for specific ALL-CAPS display names that are
# NOT acronyms (DEET, PFAS-type acronyms, etc. are deliberately left
# untouched - only exact-matched below) ----
name_fixes <- c(
  "OMEGA-HYDROXYDODECANOATE" = "Omega-hydroxydodecanoate",
  "BETA-MYRCENE" = "Beta-myrcene",
  "2-ETHYLHEXANOL" = "2-ethylhexanol",
  "GAMMA-TERPINENE" = "Gamma-terpinene",
  "PELARGONALDEHYDE" = "Pelargonaldehyde",
  "CAMPHOR" = "Camphor",
  "MENTHOL" = "Menthol",
  "2-PHENOXYETHANOL" = "2-phenoxyethanol",
  "CARVONE" = "Carvone",
  "2-METHYLNAPHTHALENE" = "2-methylnaphthalene",
  "1,1-DIMETHYL-2-PHENYLETHYL ACETATE" = "1,1-dimethyl-2-phenylethyl acetate",
  "1-METHYLNAPHTHALENE" = "1-methylnaphthalene",
  "BIPHENYL" = "Biphenyl",
  "ALPHA-CEDRENE" = "Alpha-cedrene",
  "ALPHA-HEXYLCINNAMALDEHYDE" = "Alpha-hexylcinnamaldehyde",
  "BENZYL BENZOATE" = "Benzyl benzoate",
  "DOCOSANE" = "Docosane",
  "ADIPIC ACID DI-2-ETHYLHEXYL ESTER" = "Adipic acid di-2-ethylhexyl ester"
)
apply_name_fixes <- function(x) {
  hit <- match(x, names(name_fixes))
  ifelse(!is.na(hit), name_fixes[hit], x)
}

# ---- Classifiers ----
assign_class_lc <- function(name) {
  n <- tolower(name)
  if (grepl("phosphate|phosphoric|tris\\(2-|tris\\(1,3|triethylphosph|diphenyl phosph|bis\\(2-ethylhexyl\\) phos|bis\\(1,3-dichloro|tris\\(chloropropyl\\)|tri-n-butyl phos|tri-o-cresyl|2-ethylhexyl diphenyl", n) &&
      !grepl("pfas|perfluoro|glycero|phosphocholine|phosphatidyl|phospho-.*inositol|phosph.*acid.*fatty|inositol|triazine", n))
    return("OPFRs")
  if (grepl("phthalate|phthalic", n)) return("Phthalates")
  if (grepl("perfluoro|pfas|pfoa|pfos|pfbs|pfhx|fosaa|fbsa|trifluoroacet|trifluoromethane", n)) return("PFAS")
  if (grepl("ethylene glycol|propylene glycol|tetraglyme|triglyme|diglyme|monoglyme|methoxyethoxy|ethoxyethoxy|butoxyethoxy|glycol dim|glycol mono|polyethylene|oxatetrac|3,6,9,12-tetra", n)) return("PEGs / Glycols")
  if (grepl("atrazin|terbutryn|terbutylazin|metolachlor|cyanazin|prometon|fenuron|pendimethalin|metalaxyl|carbofuran|carbendazim|malathion|imidacloprid|thiamethoxam|chlorantraniliprole|chlorothalonil|fluopyram|pyrimethanil|prosulfocarb|fipronil|bioallethrin|bentazon|2,4-dichlorophenoxy|3,5,6-trichloro-2-pyridinol|3-phenoxybenzoic|piperonyl|deet|icaridin|triclocarban|dcoit|pentachlorophenol|2-methyl-4-isothiazol", n)) return("Pesticides / Biocides")
  if (grepl("oxybenzone|octocrylene|benzophenone|methylbenzylidene camphor|ethyl 4-\\(dimethylamino\\)benz|crotamiton|paraben|methyl 4-hydroxy|tert-butyl 4-hydroxy|propylparaben", n)) return("UV filters / Personal care")
  if (grepl("atenolol|carbamazepine|cotinine|hydroxycotinine|diclofenac|fluconazole|gabapentin|gemfibrozil|levetiracetam|metformin|o-desmethylvenlafaxine|paracetamol|pseudoephedrine|salbutamol|sulfadiazine|sulfamethoxazole|topiramate|trimethoprim|benzoylecgonine|nornicotine|myosmine|nicotine|caffeine|theobromine|anisomycin", n)) return("Stimulants")
  if (grepl("galaxolide|linalool|nootkatone|sclareolide|myrtenal|cashmeran|camphor|ambrox|galaxolide|muskon", n)) return("Fragrances / Musks")
  if (grepl("diphenylguanidine|di-o-tolylguanidine|benzothiazole|mercaptobenzothiazole|6ppd|hexamethoxymethylmelamine|cyclohexyl-benzothiazole|n-cyclohexyl-2-benzothiazol|dicyclohexylamine", n)) return("Tires / Rubber")
  if (grepl("diphenyl sulfone|diphenylsulfoxide|hexamethylphosphoramide|isophorone|n,n-dimethyl-n'-p-tolyl", n)) return("Industrial chemicals")
  if (grepl("lauryl|dodecylbenzene|diethanolamine|polidocanol|myristyl sulfate|lauroylethanolamine|diethanolamide|dodecanoyl|trimethyltetradecyl|cetrimonium", n)) return("Surfactants / Detergents")
  if (grepl("phosphocholine|glycerophosph|phospho-.*inositol|palmitoyl|hexadecanoyl|dilinoleoyl|galcer|cholesterol|palmitelaidic|hydroxydodecanoate|hydroxydecanoic|hydroxyoctadec|hydroxyhexa|brassylic|decanedioic|tg\\(|dg\\(|glycan", n)) return("Lipids / Fatty acids")
  if (grepl("adenine|adenosine|uridine|tocopherol|daidzein|hesperidin|kaempferol|genkwanin|sakuranetin|diosmin|myrcene|usnic|terrein|rupestonic|bisabolol|dihydroactinidiolide|loliolide|dimethylfraxetin|dehydrocostus|petasol|4-coumaric|p-coumaric|cinnamic|gentisyl|proline|tyrosine|pipecolic|nicotinamide|pyridoxine|histamine|trehalose|melezitose|urocanic|sucralose|acesulfam|sorbic acid|melamine|hydroquinone|pyrogallol|adenosin|theobromin|caffeine|cotinin", n)) return("Biogenic / Natural")
  return("Other / Unclassified")
}

assign_class_gc <- function(name) {
  n <- tolower(name)
  if (grepl("phthalate|phthalic", n)) return("Phthalates")
  if (grepl("fluoranthene|pyrene|dibenzofuran|dibenzothiophen|fluorene|acenaphth|tetramethyl-biphenyl|diisopropylnaphthalene", n)) return("PAHs (4+ rings)")
  if (grepl("naphthalene|methylnaphthalene|biphenyl|phenanthrene|anthracene", n)) return("PAHs (2-3 rings)")
  if (grepl("toluene|xylene|trimethylbenzene|diisopropylnaphthalene|benzene, 1|benzene, \\(1-|benzene, \\(2-|ethylbenzene", n)) return("BTEX / Alkylbenzenes")
  if (grepl("hexachlorobenzene|pentachlorophenol|chloroxylenol|bromophenyl|2-bromo|difluoro|trichloromethyl|tetramethyl-biphenyl|4-bromophenyl", n)) return("OCPs / Chlorophenols")
  if (grepl("myrcene|pinene|limonene|terpinene|camphor|borneol|menthol|terpineol|carvone|geranyl|bornyl|citronellol|cedrene|ionone|cashmeran|nerolidol|ambrox|guaiol|liguloxide|amberonne|terpinyl|hexylcinnamal|sclareolide|linalool|bisabolol|nootkatone", n)) return("Terpenes / Fragrance")
  if (grepl("salicylate|benzyl benzoate|phenoxyethanol|phenoxypropan|diphenyl ether|lilial|cyclohexyl sal|hexyl sal|2-ethylhexyl sal|n-hexyl sal|benzyl sal", n)) return("UV filters / Fragrance esters")
  if (grepl("diisopropyl adipate|hexanoic acid, hexyl|adipic acid|bornyl acetate|geranyl acetate|isoamyl|dimethyl glutarate|tetrahydrofurfuryl|isopropyl myristate|isopropyl palmitate|methyl dehydroabietate|benzyl benzoate|tributyl|acetylcitrate|succinic acid|glutaric acid|oxalic acid|isophthalic|terephthalic|benzoic acid|butanedioic|dodecanoic acid, methyl|hexanoic acid|gamma-decalactone|gamma-undecalactone|delta-decalacton|dibutyl succ|dibutyl adipate|ethyl.*hexanoate|cyclopenta.*methyl ester|benzoate", n)) return("Esters / Plasticizers")
  if (grepl("^nonane$|^undecane$|tridecane|tetradecane|hexadecane|heptadecane|eicosane|docosane|nonadecane|heneicosane|pentane, 2,2,4|octane, 2,4,6|butane, 2,2-|hexane, 3,3|dimethyl octane|heptamethyl|undecane, |octane, 3,3|nonane, 2-methyl|benzene, \\(1-butyl|benzene, \\(1-ethyl|benzene, \\(1-methyl", n)) return("Alkanes / Waxes")
  if (grepl("aldehyde|decanal|dodecanal|benzaldehyde|nonanal|pelargonaldehyde|decylaldehyde|octanal|hexanal", n)) return("Aldehydes")
  if (grepl("2-ethylhexanol|3-methoxy-3-methylbutanol|ethyl tert-butyl ether|propanol|butanol|borneol|terpineol|menthol|dodecenol|glucitol|methyl.*diol|hydroxymethyl|2-phenoxyethanol|phenoxypropan|1,4,7-trimethyl.*diol", n)) return("Alcohols / Ethers")
  return("Other / Unclassified")
}

# ---- Shared plotting function: single merged panel, class-grouped ----
# ---- Explicit class order: exogenous/anthropogenic first, endogenous/
# biogenic last (right side of x-axis), rather than alphabetical ----
lc_class_order <- c("OPFRs", "Phthalates", "PFAS", "PEGs / Glycols", "Tires / Rubber",
                     "Pesticides / Biocides", "UV filters / Personal care", "Industrial chemicals",
                     "Surfactants / Detergents", "Stimulants", "Fragrances / Musks",
                     "Lipids / Fatty acids", "Biogenic / Natural")
gc_class_order <- c("Phthalates", "OCPs / Chlorophenols", "Esters / Plasticizers",
                     "BTEX / Alkylbenzenes", "PAHs (2-3 rings)", "PAHs (4+ rings)",
                     "Alkanes / Waxes", "Aldehydes", "Alcohols / Ethers",
                     "UV filters / Fragrance esters", "Terpenes / Fragrance")

# ---- Prioritize top-N compounds per class (by Total), PLUS force-include
# the most outdoor-leaning compound(s) per class even if they wouldn't
# otherwise make the cut by total magnitude - so classes don't read as
# "everything is indoor" when real outdoor-dominant compounds exist
# (e.g. tire-wear markers) ----
top_n_per_class <- function(long_dt, n_default = 5, overrides = list(), n_outdoor_force = 1) {
  wide <- dcast(long_dt, Name + Class + Total ~ Deploy, value.var = "Value")
  wide[, outdoor_ratio := Outdoor / (Indoor + Outdoor + 1e-9)]

  setorder(wide, Class, -Total)
  wide[, rank_total := seq_len(.N), by = Class]
  keep_n <- sapply(wide$Class, function(cl) if (!is.null(overrides[[cl]])) overrides[[cl]] else n_default)

  setorder(wide, Class, -outdoor_ratio)
  wide[, rank_outdoor := seq_len(.N), by = Class]

  keep <- wide[rank_total <= keep_n | rank_outdoor <= n_outdoor_force]
  long_dt[Name %in% keep$Name]
}

# Classes to show in full (or a higher cap) rather than top-5, since
# they're small to begin with or central to the paper's story
class_overrides <- list(
  "OPFRs"          = 13,   # small class, all individually interesting
  "PEGs / Glycols" = 10,   # central to the R3 "ubiquitous unregulated" story
  "PFAS"           = 10
)

build_merged_panel <- function(long_dt, value_label, plot_title, plot_subtitle, log_scale = FALSE, class_order = NULL) {
  # Order: by explicit class_order if given (exogenous -> endogenous),
  # else alphabetical; then by Total descending within class
  if (!is.null(class_order)) {
    present <- unique(as.character(long_dt$Class))
    ordered_classes <- c(class_order[class_order %in% present], setdiff(present, class_order))
    long_dt[, Class := factor(Class, levels = ordered_classes)]
    setorder(long_dt, Class, -Total)
  } else {
    setorder(long_dt, Class, -Total)
  }
  ordered_names <- unique(long_dt$Name)
  long_dt[, Name := factor(Name, levels = rev(ordered_names))]  # reversed so first class reads top-to-bottom
  long_dt[, x_num := as.numeric(Name)]

  # Class boundaries for dashed dividers + label positions (single
  # column - no faceting, which was misaligning boundaries because
  # each free_y facet panel independently rescaled its own row range)
  class_bounds <- long_dt[, .(ymin = min(x_num), ymax = max(x_num)), by = Class]
  class_bounds[, ymid := (ymin + ymax) / 2]
  setorder(class_bounds, ymax)
  divider_y <- class_bounds$ymax[-nrow(class_bounds)] + 0.5

  x_label <- if (log_scale) max(long_dt$Value[long_dt$Value > 0], na.rm = TRUE) * 1.3 else max(long_dt$Value, na.rm = TRUE) * 1.02

  p <- ggplot(long_dt, aes(x = Value, y = Name, fill = Deploy))

  if (!log_scale) {
    # Darker reference gridlines at 25/50/75/100 for the 0-100% DF plots
    p <- p + geom_vline(xintercept = c(25, 50, 75, 100), color = "grey55", linewidth = 0.4)
  }

  p <- p +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = divider_y, linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_text(data = class_bounds, aes(x = x_label, y = ymid, label = Class),
              inherit.aes = FALSE, fontface = "bold", size = 3, hjust = 0) +
    scale_fill_manual(values = c("Indoor" = col_in, "Outdoor" = col_out),
                       guide = guide_legend(reverse = TRUE)) +
    coord_cartesian(clip = "off") +
    labs(title = plot_title, subtitle = plot_subtitle, y = NULL, x = value_label, fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_text(size = 7.5),
          plot.margin = margin(10, 120, 10, 10),
          legend.position = "top",
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())

  if (log_scale) {
    max_val <- max(long_dt$Value, na.rm = TRUE)
    # Powers of ten only - the previous mixed 1/5/10/50/100/500... breaks
    # were spaced too densely under pseudo-log compression and collided
    # into unreadable merged labels (e.g. "50" + "100" -> "50100")
    brks <- 10^(0:8)
    brks <- c(0, brks[brks <= max_val * 2])
    p <- p + scale_x_continuous(trans = scales::pseudo_log_trans(base = 10),
                                  breaks = brks,
                                  labels = scales::label_number(big.mark = ",", accuracy = 1))
  }
  p
}

# ---- Export the exact data behind a DF panel, in the exact order the
# compounds appear top-to-bottom in the plot (mirrors the ordering logic
# inside build_merged_panel: class_order first, then -Total within class;
# ordered_names as computed there reads top-to-bottom in the figure) ----
export_plot_data <- function(long_dt, class_order, out_csv) {
  present <- unique(as.character(long_dt$Class))
  ordered_classes <- c(class_order[class_order %in% present], setdiff(present, class_order))
  long_dt[, Class := factor(Class, levels = ordered_classes)]
  setorder(long_dt, Class, -Total)
  ordered_names <- unique(long_dt$Name)

  wide <- dcast(long_dt, Name + Class ~ Deploy, value.var = "Value")
  setnames(wide, c("Indoor", "Outdoor"), c("Indoor_detection_freq", "Outdoor_detection_freq"))
  wide[, Plot_order := match(Name, ordered_names)]
  setorder(wide, Plot_order)

  fwrite(wide[, .(Plot_order, Name, Class, Indoor_detection_freq, Outdoor_detection_freq)], out_csv)
  cat("  Exported:", out_csv, "(", nrow(wide), "compounds )\n")
}

# ---- Split a long_dt into two halves by class (for two separate panels
# the user will merge manually), roughly balanced by row count ----
split_by_class_half <- function(long_dt, class_order) {
  present <- unique(as.character(long_dt$Class))
  ordered_classes <- c(class_order[class_order %in% present], setdiff(present, class_order))
  sizes <- long_dt[, .N, by = Class][, .(Class, n_rows = N / 2)]
  sizes[, Class := factor(Class, levels = ordered_classes)]
  setorder(sizes, Class)
  sizes[, cum_n := cumsum(n_rows)]
  half_point <- sum(sizes$n_rows) / 2
  split_class <- as.character(sizes[cum_n >= half_point, Class][1])
  split_idx <- which(ordered_classes == split_class)
  list(
    part1 = long_dt[Class %in% ordered_classes[1:split_idx]],
    part2 = long_dt[Class %in% ordered_classes[(split_idx + 1):length(ordered_classes)]]
  )
}

# ============================================================
# LC: read, classify, compute DF% and mean concentration
# ============================================================
cat("Reading LC target+annotation concentrations...\n")
lc_raw <- as.data.table(read_excel(lc_file))

lc_sample_cols <- names(lc_raw)[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?\\s*$", names(lc_raw))]
lc_in_cols  <- lc_sample_cols[grepl("_Indoor_", lc_sample_cols)]
lc_out_cols <- lc_sample_cols[grepl("_Outdoor_", lc_sample_cols)]
cat("  LC compounds:", nrow(lc_raw), "| Indoor cols:", length(lc_in_cols), "| Outdoor cols:", length(lc_out_cols), "\n")

for (col in c(lc_in_cols, lc_out_cols)) {
  lc_raw[, (col) := as.numeric(get(col))]
  lc_raw[is.na(get(col)), (col) := 0]
}
lc_raw[, mean_indoor  := rowMeans(.SD, na.rm = TRUE), .SDcols = lc_in_cols]
lc_raw[, mean_outdoor := rowMeans(.SD, na.rm = TRUE), .SDcols = lc_out_cols]
lc_raw[, DF_indoor  := rowSums(.SD > 0) / length(lc_in_cols) * 100, .SDcols = lc_in_cols]
lc_raw[, DF_outdoor := rowSums(.SD > 0) / length(lc_out_cols) * 100, .SDcols = lc_out_cols]
lc_raw[, Name := apply_name_fixes(Name)]
lc_raw[, Class := sapply(Name, assign_class_lc)]
lc_raw <- lc_raw[Class != "Other / Unclassified"]

# NOTE: excluded compound names are applied AFTER top_n_per_class selection
# below (not here), specifically to avoid backfill: removing a name from
# the pool before selection just lets top_n_per_class pick a replacement
# from further down the ranking, silently introducing compounds that were
# never part of the original selection. Excluding after selection instead
# just shrinks that class - nothing new ever enters.
lc_exclude_names <- c(
  "Dicyclohexylamine",
  "Bis(2-ethylhexyl) phosphate", "Tris(2-ethylhexyl) phosphate",
  "Bis(1,3-dichloro-2-propyl) phosphate", "Tri-o-cresyl phosphate",
  "Perfluoro-tetradecanoic acid", "Perfluoro-n-hexadecanoic acid",
  "Perfluoro-tridecanoic acid", "Sodium perfluoro-1-dodecanesulfonate",
  "Perfluoro-n-hexanoic acid", "N-EtFOSAA",
  "Pendimethalin", "Diphenyl sulfone", "Hexamethylphosphoramide", "Atenolol"
)

lc_long_df <- rbind(
  lc_raw[, .(Name, Class, Total = DF_indoor + DF_outdoor, Deploy = "Indoor",  Value = DF_indoor)],
  lc_raw[, .(Name, Class, Total = DF_indoor + DF_outdoor, Deploy = "Outdoor", Value = DF_outdoor)]
)
lc_long_conc <- rbind(
  lc_raw[, .(Name, Class, Total = mean_indoor + mean_outdoor, Deploy = "Indoor",  Value = mean_indoor)],
  lc_raw[, .(Name, Class, Total = mean_indoor + mean_outdoor, Deploy = "Outdoor", Value = mean_outdoor)]
)
lc_long_df[, Deploy := factor(Deploy, levels = c("Outdoor", "Indoor"))]
lc_long_conc[, Deploy := factor(Deploy, levels = c("Outdoor", "Indoor"))]

lc_df_filtered   <- top_n_per_class(lc_long_df, overrides = class_overrides)
lc_conc_filtered <- top_n_per_class(lc_long_conc, overrides = class_overrides)

# Apply requested exclusions AFTER selection - shrinks the class, doesn't backfill
lc_df_filtered   <- lc_df_filtered[!(Name %in% lc_exclude_names)]
lc_conc_filtered <- lc_conc_filtered[!(Name %in% lc_exclude_names)]

# ---- Export Figure 4a source data: exact compounds + order as plotted ----
export_plot_data(copy(lc_df_filtered), lc_class_order, file.path(out_dir, "Fig4a_LC_DetectionFrequency_SourceData.csv"))

lc_df_split   <- split_by_class_half(lc_df_filtered, lc_class_order)
lc_conc_split <- split_by_class_half(lc_conc_filtered, lc_class_order)

p_lc_df_1 <- build_merged_panel(lc_df_split$part1, "Detection frequency (%)",
             "LC-HRMS Target & Annotated Compounds (1/2)", "Detection frequency, grouped by class",
             class_order = lc_class_order)
p_lc_df_2 <- build_merged_panel(lc_df_split$part2, "Detection frequency (%)",
             "LC-HRMS Target & Annotated Compounds (2/2)", "Detection frequency, grouped by class",
             class_order = lc_class_order)
p_lc_conc_1 <- build_merged_panel(lc_conc_split$part1, "Mean concentration (log scale)",
             "LC-HRMS Target & Annotated Compounds (1/2)", "Mean concentration, grouped by class",
             log_scale = TRUE, class_order = lc_class_order)
p_lc_conc_2 <- build_merged_panel(lc_conc_split$part2, "Mean concentration (log scale)",
             "LC-HRMS Target & Annotated Compounds (2/2)", "Mean concentration, grouped by class",
             log_scale = TRUE, class_order = lc_class_order)

ggsave(file.path(out_dir, "Fig_LC_MergedClassPanel_DF_part1.png"), p_lc_df_1, width = 8, height = 10, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_LC_MergedClassPanel_DF_part2.png"), p_lc_df_2, width = 8, height = 10, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_LC_MergedClassPanel_Concentration_part1.png"), p_lc_conc_1, width = 8, height = 10, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_LC_MergedClassPanel_Concentration_part2.png"), p_lc_conc_2, width = 8, height = 10, dpi = 600, limitsize = FALSE)
cat("Saved LC merged panels (DF + concentration), 2 parts each\n")

# ---- Also save unsplit, single full-panel versions ----
p_lc_df_full <- build_merged_panel(lc_df_filtered, "Detection frequency (%)",
             "LC-HRMS Target & Annotated Compounds", "Detection frequency, grouped by class",
             class_order = lc_class_order)
p_lc_conc_full <- build_merged_panel(lc_conc_filtered, "Mean concentration (log scale)",
             "LC-HRMS Target & Annotated Compounds", "Mean concentration, grouped by class",
             log_scale = TRUE, class_order = lc_class_order)

ggsave(file.path(out_dir, "Fig_LC_MergedClassPanel_DF_full.png"), p_lc_df_full, width = 8, height = 20, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_LC_MergedClassPanel_Concentration_full.png"), p_lc_conc_full, width = 8, height = 20, dpi = 600, limitsize = FALSE)
cat("Saved LC merged panels (DF + concentration), unsplit full versions\n")

# ============================================================
# GC: read, filter ID level 1-2, classify, compute DF% and mean area
# ============================================================
cat("\nReading GC data...\n")
gc_raw <- as.data.table(read_excel(gc_file))
gc_raw <- gc_raw[!is.na(as.numeric(RT))]
gc_ann <- gc_raw[`ID level` %in% c(1, 2)]
cat("  GC annotated features (L1-2 only):", nrow(gc_ann), "\n")

gc_in_cols  <- names(gc_ann)[grepl("^[A-Z]{2}_Indoor_[A-Z]{2}$", names(gc_ann))]
gc_out_cols <- names(gc_ann)[grepl("^[A-Z]{2}_Outdoor_[A-Z]{2}$", names(gc_ann))]

for (col in c(gc_in_cols, gc_out_cols)) {
  gc_ann[, (col) := as.numeric(get(col))]
  gc_ann[is.na(get(col)), (col) := 0]
}
gc_ann[, mean_indoor  := rowMeans(.SD, na.rm = TRUE), .SDcols = gc_in_cols]
gc_ann[, mean_outdoor := rowMeans(.SD, na.rm = TRUE), .SDcols = gc_out_cols]
gc_ann[, DF_indoor  := rowSums(.SD > 0) / length(gc_in_cols) * 100, .SDcols = gc_in_cols]
gc_ann[, DF_outdoor := rowSums(.SD > 0) / length(gc_out_cols) * 100, .SDcols = gc_out_cols]
gc_ann[, Class := sapply(`Database match`, assign_class_gc)]
gc_ann[, Name  := apply_name_fixes(`Database match`)]
gc_ann <- gc_ann[Class != "Other / Unclassified"]
gc_ann <- gc_ann[Class != "Phthalates"]

gc_long_df <- rbind(
  gc_ann[, .(Name, Class, Total = DF_indoor + DF_outdoor, Deploy = "Indoor",  Value = DF_indoor)],
  gc_ann[, .(Name, Class, Total = DF_indoor + DF_outdoor, Deploy = "Outdoor", Value = DF_outdoor)]
)
gc_long_area <- rbind(
  gc_ann[, .(Name, Class, Total = mean_indoor + mean_outdoor, Deploy = "Indoor",  Value = mean_indoor)],
  gc_ann[, .(Name, Class, Total = mean_indoor + mean_outdoor, Deploy = "Outdoor", Value = mean_outdoor)]
)
gc_long_df[, Deploy := factor(Deploy, levels = c("Outdoor", "Indoor"))]
gc_long_area[, Deploy := factor(Deploy, levels = c("Outdoor", "Indoor"))]

gc_df_filtered   <- top_n_per_class(gc_long_df)
gc_area_filtered <- top_n_per_class(gc_long_area)

# ---- Export Figure 4b source data: exact compounds + order as plotted ----
export_plot_data(copy(gc_df_filtered), gc_class_order, file.path(out_dir, "Fig4b_GC_DetectionFrequency_SourceData.csv"))

gc_df_split   <- split_by_class_half(gc_df_filtered, gc_class_order)
gc_area_split <- split_by_class_half(gc_area_filtered, gc_class_order)

p_gc_df_1 <- build_merged_panel(gc_df_split$part1, "Detection frequency (%)",
             "GC-HRMS Annotated Compounds (ID level 1-2) (1/2)", "Detection frequency, grouped by class",
             class_order = gc_class_order)
p_gc_df_2 <- build_merged_panel(gc_df_split$part2, "Detection frequency (%)",
             "GC-HRMS Annotated Compounds (ID level 1-2) (2/2)", "Detection frequency, grouped by class",
             class_order = gc_class_order)
p_gc_area_1 <- build_merged_panel(gc_area_split$part1, "Mean summed peak area (log scale)",
             "GC-HRMS Annotated Compounds (ID level 1-2) (1/2)", "Mean summed peak area, grouped by class",
             log_scale = TRUE, class_order = gc_class_order)
p_gc_area_2 <- build_merged_panel(gc_area_split$part2, "Mean summed peak area (log scale)",
             "GC-HRMS Annotated Compounds (ID level 1-2) (2/2)", "Mean summed peak area, grouped by class",
             log_scale = TRUE, class_order = gc_class_order)

ggsave(file.path(out_dir, "Fig_GC_MergedClassPanel_DF_part1.png"), p_gc_df_1, width = 7, height = 9, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_GC_MergedClassPanel_DF_part2.png"), p_gc_df_2, width = 7, height = 9, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_GC_MergedClassPanel_Area_part1.png"), p_gc_area_1, width = 7, height = 9, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_GC_MergedClassPanel_Area_part2.png"), p_gc_area_2, width = 7, height = 9, dpi = 600, limitsize = FALSE)
cat("Saved GC merged panels (DF + area), 2 parts each\n")

# ---- Also save unsplit, single full-panel versions ----
p_gc_df_full <- build_merged_panel(gc_df_filtered, "Detection frequency (%)",
             "GC-HRMS Annotated Compounds (ID level 1-2)", "Detection frequency, grouped by class",
             class_order = gc_class_order)
p_gc_area_full <- build_merged_panel(gc_area_filtered, "Mean summed peak area (log scale)",
             "GC-HRMS Annotated Compounds (ID level 1-2)", "Mean summed peak area, grouped by class",
             log_scale = TRUE, class_order = gc_class_order)

ggsave(file.path(out_dir, "Fig_GC_MergedClassPanel_DF_full.png"), p_gc_df_full, width = 8.5, height = 16, dpi = 600, limitsize = FALSE)
ggsave(file.path(out_dir, "Fig_GC_MergedClassPanel_Area_full.png"), p_gc_area_full, width = 10, height = 16, dpi = 600, limitsize = FALSE)
cat("Saved GC merged panels (DF + area), unsplit full versions\n")
