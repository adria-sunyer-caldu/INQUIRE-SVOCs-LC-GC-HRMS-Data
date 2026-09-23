library(data.table)
library(readxl)
library(ggplot2)
library(dplyr)

# ==========================================================
# Total SVOC chemical burden - real per-household ratios, both
# platforms, feeding a horizontal raincloud plot (replacement for
# the two 8-facet violin grids, panels c/d).
#
# LC-HRMS: deployment (Indoor/Outdoor) comes from the "Class" row
#   (row 1 of the raw CSV), matched by COLUMN POSITION to the
#   sample columns - NOT from the "_1_"/"_2_" suffix in the sample
#   name, which is not a reliable indoor/outdoor flag (confirmed
#   against the 10-row format example and the original LC script's
#   own parse_sample_v2() logic).
#   Country and household code are read from the pseudonymized sample
#   name: {Country}_{Indoor|Outdoor}_{code} (indoor and outdoor samples of
#   the same household share the same code).
#
# GC-HRMS: deployment comes directly from the sample name suffix
#   ({Country}_Indoor_{code} / {Country}_Outdoor_{code}), country from the
#   name prefix. Extra indoor samplers (_2/_3 suffix) are not used.
# ==========================================================

# ---------- Paths (edit these) ----------
data_dir <- "path/to/data"
lc_path <- file.path(data_dir, "LC_All_Features_FULL_final_zenodo.xlsx")   # Zenodo LC deposit
gc_path <- file.path(data_dir, "GC_Level1_5_polished_FULL_zenodo.xlsx")   # Zenodo GC deposit
out_dir <- data_dir
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==========================================================
# STEP 1: LC-HRMS - real per-household total burden ratios
# ==========================================================
cat("Reading LC header rows...\n")
lc_header <- as.data.frame(read_excel(lc_path, sheet = "Full_dataset_INQUIRE_nontarget_", col_names = FALSE,
                                      n_max = 4, .name_repair = "minimal"))
class_row   <- as.character(unlist(lc_header[1, ]))
colname_row <- trimws(as.character(unlist(lc_header[4, ])))
sample_col_idx <- which(grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}(_[23])?$", colname_row))
cat("LC sample columns found:", length(sample_col_idx), "\n")

sample_names <- colname_row[sample_col_idx]
deploy       <- class_row[sample_col_idx]   # Indoor / Outdoor, position-matched

# Parse house + country from the sample name (these two ARE reliable from the name)
parse_lc_name <- function(s) {
  m <- regmatches(s, regexec("^([A-Z]{2})_(Indoor|Outdoor)_([A-Z]{2}(_[23])?)$", s))[[1]]
  if (length(m) == 0) return(c(NA, NA))
  c(m[4], m[2])   # house (household code), country
}
lc_meta <- as.data.frame(t(sapply(sample_names, parse_lc_name)), stringsAsFactors = FALSE)
colnames(lc_meta) <- c("House", "Country")
lc_meta$SampleName <- sample_names
lc_meta$Deploy     <- deploy
lc_meta$ColIndex   <- sample_col_idx

cat("Reading LC data (sample columns only)...\n")
# Read the feature rows (4 metadata rows skipped), then keep the sample
# columns by position
lc_dt <- as.data.table(read_excel(lc_path, sheet = "Full_dataset_INQUIRE_nontarget_", skip = 4,
                                  col_names = FALSE, .name_repair = "minimal"))
lc_dt <- lc_dt[, sample_col_idx, with = FALSE]
setnames(lc_dt, as.character(sample_col_idx))
lc_dt[, (names(lc_dt)) := lapply(.SD, as.numeric)]

cat("LC feature rows loaded:", nrow(lc_dt), "\n")

# Total burden per sample = column sum across all features
lc_burden <- colSums(lc_dt, na.rm = TRUE)
lc_meta$Burden <- lc_burden[as.character(lc_meta$ColIndex)]

# Pair Indoor/Outdoor by House + Country, compute log ratio
lc_pairs <- lc_meta %>%
  group_by(House, Country) %>%
  summarise(
    Indoor  = Burden[Deploy == "Indoor"][1],
    Outdoor = Burden[Deploy == "Outdoor"][1],
    .groups = "drop"
  ) %>%
  filter(!is.na(Indoor), !is.na(Outdoor), Indoor > 0, Outdoor > 0) %>%
  mutate(LogRatio = log(Indoor / Outdoor), Platform = "LC-HRMS")

cat("LC valid pairs:", nrow(lc_pairs), "\n\n")

# ==========================================================
# STEP 2: GC-HRMS - real per-household total burden ratios
# ==========================================================
cat("Reading GC file...\n")
gc_raw <- read_excel(gc_path, col_names = TRUE)

gc_cols <- colnames(gc_raw)
gc_sample_cols <- gc_cols[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", gc_cols)]
cat("GC sample columns found:", length(gc_sample_cols), "\n")

parse_gc_name <- function(s) {
  m <- regmatches(s, regexec("^([A-Z]{2})_(Indoor|Outdoor)_([A-Z]{2})$", s))[[1]]
  c(house = m[4], country = m[2], type = m[3])
}
gc_meta <- as.data.frame(t(sapply(gc_sample_cols, parse_gc_name)), stringsAsFactors = FALSE)
colnames(gc_meta) <- c("House", "Country", "Type")
gc_meta$SampleName <- gc_sample_cols

gc_num <- as.data.frame(lapply(gc_raw[, gc_sample_cols], function(x) suppressWarnings(as.numeric(x))))
gc_burden <- colSums(gc_num, na.rm = TRUE)
gc_meta$Burden <- gc_burden[gc_meta$SampleName]

gc_pairs <- gc_meta %>%
  group_by(House, Country) %>%
  summarise(
    Indoor  = Burden[Type == "Indoor"][1],
    Outdoor = Burden[Type == "Outdoor"][1],
    .groups = "drop"
  ) %>%
  filter(!is.na(Indoor), !is.na(Outdoor), Indoor > 0, Outdoor > 0) %>%
  mutate(LogRatio = log(Indoor / Outdoor), Platform = "GC-HRMS")

cat("GC valid pairs:", nrow(gc_pairs), "\n\n")

# ==========================================================
# STEP 3: Combine and plot (real data, horizontal raincloud)
# ==========================================================
combined <- bind_rows(
  lc_pairs %>% select(Country, Platform, LogRatio),
  gc_pairs %>% select(Country, Platform, LogRatio)
)

country_order <- c("EE","PT","UK","CZ","SI","SE","IT","NL")
combined$Country[combined$Country == "SL"] <- "SI"  # unify display label (GC's SL was only for internal matching)
combined$Country <- factor(combined$Country, levels = rev(country_order))

col_lc <- "#5A2D75"  # dark purple - platform color, distinct from Indoor/Outdoor red/blue elsewhere
col_gc <- "#A85A1F"  # dark orange

combined$y_num <- as.numeric(combined$Country)
combined$y_off <- ifelse(combined$Platform == "LC-HRMS", combined$y_num + 0.22, combined$y_num - 0.22)

build_density_polygon <- function(vals, baseline, scale = 1, half = FALSE, side = 1) {
  d <- density(vals, n = 512)
  h <- d$y / max(d$y) * scale
  if (half) {
    data.frame(x = c(d$x, rev(d$x)), y = c(baseline + h * side, rep(baseline, length(d$x))))
  } else {
    data.frame(x = d$x, y = baseline + h)
  }
}

half_polys <- do.call(rbind, lapply(unique(paste(combined$Country, combined$Platform)), function(g) {
  sub <- combined[paste(combined$Country, combined$Platform) == g, ]
  if (nrow(sub) < 2) return(NULL)
  poly <- build_density_polygon(sub$LogRatio, baseline = sub$y_off[1], scale = 0.16, half = TRUE, side = 1)
  poly$Country <- sub$Country[1]; poly$Platform <- sub$Platform[1]
  poly
}))

box_stats <- combined %>% group_by(Country, Platform, y_off) %>%
  summarise(ymin = quantile(LogRatio, 0.25) - 1.5*IQR(LogRatio),
            lower = quantile(LogRatio, 0.25), med = median(LogRatio),
            upper = quantile(LogRatio, 0.75),
            ymax = quantile(LogRatio, 0.75) + 1.5*IQR(LogRatio), .groups = "drop")

p_raincloud <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_polygon(data = half_polys, aes(x = x, y = y, group = interaction(Country, Platform), fill = Platform, color = Platform),
               alpha = 0.55, linewidth = 0.5) +
  geom_boxplot(data = box_stats, aes(xmin = ymin, xlower = lower, xmiddle = med, xupper = upper, xmax = ymax,
                                       y = y_off - 0.05, group = interaction(Country, Platform), fill = Platform),
               stat = "identity", width = 0.06, alpha = 0.9, orientation = "y",
               color = "grey15", linewidth = 0.9) +
  geom_jitter(data = combined, aes(x = LogRatio, y = y_off - 0.11, color = Platform),
              width = 0, height = 0.03, size = 1.6, alpha = 0.6) +
  scale_y_continuous(breaks = seq_along(levels(combined$Country)), labels = levels(combined$Country)) +
  scale_fill_manual(values = c("LC-HRMS" = col_lc, "GC-HRMS" = col_gc)) +
  scale_color_manual(values = c("LC-HRMS" = col_lc, "GC-HRMS" = col_gc)) +
  labs(title = "Total SVOC Chemical Burden - Indoor/Outdoor Ratio by Country",
       subtitle = "Real data, both platforms",
       x = "log(Indoor/Outdoor total burden ratio)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

print(p_raincloud)
ggsave(file.path(out_dir, "Fig_Raincloud_LC_GC_ByCountry.png"), p_raincloud, dpi = 600, width = 8, height = 12)
write.csv(combined, file.path(out_dir, "Combined_LC_GC_ByCountry_Ratios.csv"), row.names = FALSE)
cat("Saved outputs to:", out_dir, "\n")
