# ============================================================================
# INQUIRE — Score-5 compound detection frequency by country (LC + GC combined)
# Builds the heatmap for main-text Figure 6 (12 compounds x 8 countries).
# "Detected" = concentration/area > 0 in that sample (matches the convention
# used elsewhere in the toxicity pipeline). Indoor + Outdoor samples pooled
# per country. If your pipeline instead defines "detected" via LOD/LOQ
# thresholds rather than simple >0, adjust the detected_lc()/detected_gc()
# functions below accordingly before trusting the output for publication.
# ============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

# ---- Set your working/data directory here ----
setwd("path/to/data")

# ---------------------------
# 1) LC: load target concentrations
# ---------------------------
lc_file <- "LC_Level1_2_polished_FINAL.xlsx"
lc_raw <- read_excel(lc_file)
names(lc_raw) <- trimws(names(lc_raw))

meta_cols <- c("Feature ID","MZ","RT","Name","Molecular formula","InChiKey","SMILES",
               "CAS","Ion type","Confidence level","MS/MS","LOD","LOQ","Quantification")
lc_sample_cols <- setdiff(names(lc_raw), meta_cols)

lc_targets <- c(
  "Triphenyl phosphate", "Diisobutyl phthalate", "Tris(2-chloroethyl) phosphate",
  "Dihexyl phthalate", "Diethylene glycol dimethyl ether", "Benzyl butyl phthalate",
  "Dipentyl phthalate", "Di(2-methoxyethyl) phthalate"
)

lc_sub <- lc_raw %>% filter(Name %in% lc_targets)
stopifnot(nrow(lc_sub) == length(lc_targets))  # fails loudly if a name doesn't match

lc_long <- lc_sub %>%
  select(Name, all_of(lc_sample_cols)) %>%
  pivot_longer(cols = -Name, names_to = "sample", values_to = "value") %>%
  mutate(
    value = as.numeric(value),
    country = str_match(sample, "^([A-Z]{2})_(Indoor|Outdoor)_[A-Z]{2}(_[23])?$")[, 2],
    deployment = ifelse(str_detect(sample, "_Indoor_"), "Indoor", "Outdoor"),
    detected = ifelse(!is.na(value) & value > 0, 1, 0)
  ) %>%
  filter(!is.na(country))

lc_df_by_country <- lc_long %>%
  group_by(Name, country) %>%
  summarise(DF = 100 * sum(detected) / n(), .groups = "drop")

# ---------------------------
# 2) GC: load the polished GC Level 1-2 dataset
# ---------------------------
gc_file <- "GC_Level1_2_polished_FINAL.xlsx"
gc_raw <- read_excel(gc_file)

gc_targets <- c("Lilial", "Anthracene", "Fluoranthene", "Pyrene", "Dibutyl phthalate")

gc_sample_cols <- names(gc_raw)[str_detect(names(gc_raw), "^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$")]

gc_sub <- gc_raw %>% filter(`Database match` %in% gc_targets)
stopifnot(nrow(gc_sub) == length(gc_targets))  # fails loudly if a name doesn't match

gc_long <- gc_sub %>%
  select(Name = `Database match`, all_of(gc_sample_cols)) %>%
  pivot_longer(cols = -Name, names_to = "sample", values_to = "value") %>%
  mutate(
    value = as.numeric(value),
    country_raw = str_match(sample, "^([A-Z]{2})_(Indoor|Outdoor)_[A-Z]{2}$")[, 2],
    country = ifelse(country_raw == "SL", "SI", country_raw),  # raw data uses SL for Slovenia
    deployment = ifelse(str_detect(sample, "_Indoor_"), "Indoor", "Outdoor"),
    detected = ifelse(!is.na(value) & value > 0, 1, 0)
  ) %>%
  filter(!is.na(country))

gc_df_by_country <- gc_long %>%
  group_by(Name, country) %>%
  summarise(DF = 100 * sum(detected) / n(), .groups = "drop")

# ---------------------------
# 3) Combine LC + GC, order rows/columns, plot
# ---------------------------
combined <- bind_rows(lc_df_by_country, gc_df_by_country)

compound_order <- c(lc_targets, gc_targets)
country_order  <- c("PT","UK","NL","SE","CZ","EE","IT","SI")

combined <- combined %>%
  mutate(
    Name = factor(Name, levels = rev(compound_order)),  # rev so first compound plots on top
    country = factor(country, levels = country_order)
  )

# Color palette used in the published figure
palette_stops <- c("#F7EAE6","#F4D4C7","#EFBBA7","#EAA186","#DC7D5A","#C95B36","#A34829","#863A20","#703019")

p_heatmap <- ggplot(combined, aes(x = country, y = Name, fill = DF)) +
  geom_tile(color = "black") +
  geom_text(aes(label = paste0(round(DF), "%"),
                color = DF > 55), size = 5) +
  scale_fill_gradientn(colours = palette_stops,
                       values = rescale(seq(0, 100, length.out = length(palette_stops)), to = c(0,1)),
                       limits = c(0, 100),
                       name = "Detection\nfrequency (%)") +
  scale_color_manual(values = c("TRUE" = "white", "FALSE" = "black"), guide = "none") +
  labs(title = "Score-5 compound detection frequency by country",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

print(p_heatmap)
ggsave("Score5_DF_heatmap.png", p_heatmap, width = 12, height = 7, dpi = 600)
cat("Saved: Score5_DF_heatmap.png\n")
