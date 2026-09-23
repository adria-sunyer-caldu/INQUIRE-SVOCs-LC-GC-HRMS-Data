# =============================================================================
# LC_Violin_Plots.R
#
# Violin + jitter plot of one LC-HRMS compound's concentration, indoor vs
# outdoor, faceted by country (16 groups: Indoor/Outdoor x 8 countries).
# Used for the individual named-compound supplementary figures (SF18).
#
# Reads from the polished LC dataset (LC_Level1_2_polished_FINAL.xlsx),
# NOT the raw "Targets_and_annotations_with_adducts_and_spectra.xlsx" file
# an earlier version of this script used — matching every other LC script
# in this project, which reads the polished, zero-detection-filtered,
# duplicate-corrected dataset. Columns are identified by name (Name,
# Confidence level, etc.), not by hardcoded column position — the raw
# file's column 22/26+ positions do not apply to the polished file's
# narrower structure.
#
# Two y-axis scale modes, since some compounds have too wide a dynamic
# range for a fixed linear-ish break scale to work well:
#   SCALE_MODE <- "linear"     - fixed breaks, no transform
#   SCALE_MODE <- "pseudoscale" - x^(1/1.3) power transform, wider breaks
# Pick whichever looks right for the compound at hand.
#
# To plot a different compound: change CHEMICAL_TO_PLOT below and re-run.
# =============================================================================

library(dplyr)
library(ggplot2)
library(readxl)
library(tidyr)
library(scales)

# ---- Set your working/data directory here ----
base_dir <- "path/to/data"
setwd(base_dir)

lc_file <- "LC_Level1_2_polished_FINAL.xlsx"

# ---- Config: change these two for each compound ----
CHEMICAL_TO_PLOT <- "Tris(chloropropyl) phosphate"  # must match the "Name" column exactly
SCALE_MODE <- "linear"   # "linear" or "pseudoscale"

# ---- Load polished LC data ----
lc_raw <- read_excel(lc_file)
names(lc_raw) <- trimws(names(lc_raw))

meta_cols <- c("Feature ID","MZ","RT","Name","Molecular formula","InChiKey","SMILES",
               "CAS","Ion type","Confidence level","MS/MS","LOD","LOQ","Quantification")
sample_cols <- setdiff(names(lc_raw), meta_cols)
lc_raw[sample_cols] <- lapply(lc_raw[sample_cols], as.numeric)

if (!(CHEMICAL_TO_PLOT %in% lc_raw$Name)) {
  stop("'", CHEMICAL_TO_PLOT, "' not found in the Name column — check spelling against LC_Level1_2_polished_FINAL.xlsx.")
}

# ---- Reshape to long format, indoor/outdoor + country from sample column name ----
long_data <- lc_raw %>%
  filter(Name == CHEMICAL_TO_PLOT) %>%
  select(Name, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Concentration") %>%
  mutate(
    IO           = ifelse(grepl("_Indoor_", Sample), "Indoor", "Outdoor"),
    Country      = sub("_.*", "", Sample),
    Sample_Group = paste0(IO, "_", Country)
  )

group_order <- c(
  "Indoor_CZ", "Outdoor_CZ", "Indoor_EE", "Outdoor_EE", "Indoor_IT", "Outdoor_IT",
  "Indoor_NL", "Outdoor_NL", "Indoor_PT", "Outdoor_PT", "Indoor_SE", "Outdoor_SE",
  "Indoor_SI", "Outdoor_SI", "Indoor_UK", "Outdoor_UK"
)
long_data$Sample_Group <- factor(long_data$Sample_Group, levels = group_order)

group_colors <- setNames(
  rep(c("#B03A2E", "#1F618D"), 8),
  group_order
)

# ---- Build the plot ----
p <- ggplot(long_data, aes(x = Sample_Group, y = Concentration, fill = Sample_Group)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", alpha = 0.6,
              draw_quantiles = c(0.5), adjust = 1) +
  geom_jitter(width = 0.1, color = "black", size = 2, alpha = 0.5) +
  labs(title = CHEMICAL_TO_PLOT, x = "Sample Type", y = "Concentration (ng/mL)") +
  scale_fill_manual(values = group_colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.line = element_line(color = "black"),
    legend.position = "none"
  )

if (SCALE_MODE == "pseudoscale") {
  p <- p + scale_y_continuous(
    trans = scales::trans_new("sixth_root",
                              transform = function(x) x^(1/1.3),
                              inverse   = function(x) x^1.3),
    breaks = c(0, 50, 250, 500, 1000, 2500, 5000, 10000, 20000, 30000, 40000, 70000, 150000, 1500000),
    limits = c(0, NA),
    expand = c(0, 0)
  )
} else {
  p <- p + scale_y_continuous(
    breaks = c(0, 50, 500, 1000, 2500, 5000, 10000, 20000, 30000, 40000, 75000),
    limits = c(0, NA),
    expand = c(0, 0)
  )
}

print(p)

out_dir <- "For publication supplementary"
dir.create(out_dir, showWarnings = FALSE)
ggsave(file.path(out_dir, paste0(gsub("[^A-Za-z0-9]", "_", CHEMICAL_TO_PLOT), ".png")),
       p, dpi = 600, width = 12, height = 8)
cat("Saved:", file.path(out_dir, paste0(gsub("[^A-Za-z0-9]", "_", CHEMICAL_TO_PLOT), ".png")), "\n")
