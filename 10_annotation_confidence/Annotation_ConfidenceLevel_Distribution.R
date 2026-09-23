############################################################
## Annotation Confidence Level Distribution — LC-HRMS & GC-HRMS
## Bar charts (pie charts made Level 1/2 unreadably small and
## mislabeled given how tiny those slices are relative to the rest).
##
## LC feature universe = the 155,448-feature FILTERED set (zero-total
## + >=5% prevalence filters already applied upstream), taken from
## NAP_xref_all_features.csv$Alignment_ID.
##
## LC tiers: Level 1 (target, confirmed) / Level 2 (2a, library
##   annotation) / Level 3 (NAP, tentative spectral-match candidate,
##   mz_match=TRUE with a Top1_SMILES, excluding anything already in
##   Level 1/2) / Not annotated.
##
## GC: native "ID level" column already covers 1-5, no NAP needed.
############################################################

library(data.table)
library(dplyr)
library(readxl)
library(ggplot2)
library(scales)

out_dir <- "path/to/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setwd(out_dir)

bar_chart <- function(counts_df, tier_col, n_col, total, colors, title){
  ggplot(counts_df, aes(x = .data[[tier_col]], y = .data[[n_col]], fill = .data[[tier_col]])) +
    geom_col(width = 0.65) +
    geom_text(aes(label = paste0(comma(.data[[n_col]]), "\n(", round(100*.data[[n_col]]/total,1), "%)")),
              vjust = -0.15, size = 3.4) +
    scale_fill_manual(values = colors) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15))) +
    labs(title = title, x = NULL, y = "Number of features") +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none", plot.title = element_text(face = "bold"))
}

# ================================================================
# LC-HRMS
# ================================================================
cat("=== LC-HRMS ===\n")

nap <- fread("NAP_xref_all_features.csv")
n_lc_total <- nrow(nap)
cat("Total LC features (filtered universe, from NAP_xref_all_features.csv):", n_lc_total, "\n")

# ---- Confirmed targets/annotations (Level 1 / 2) ----
targets_file <- "LC_Level1_2_polished_FINAL.xlsx"
targets <- read_excel(targets_file, sheet = 1)
targets$Alignment_ID <- sub("_SU$", "", targets$`Feature ID`)
targets$Tier <- case_when(
  targets$`Confidence level` == 1    ~ "Level 1",
  targets$`Confidence level` == "2a" ~ "Level 2",
  TRUE ~ NA_character_
)
if (any(is.na(targets$Tier))) {
  warning("Unrecognized Confidence level value(s) found: ",
          paste(unique(targets$`Confidence level`[is.na(targets$Tier)]), collapse = ", "))
}
cat("Confirmed targets/annotations:", nrow(targets),
    "| Level 1:", sum(targets$Tier == "Level 1"),
    "| Level 2:", sum(targets$Tier == "Level 2"), "\n")

confirmed_ids <- targets$Alignment_ID

# ---- NAP tentative (Level 3), excluding anything already confirmed ----
nap_tentative <- nap %>%
  filter(mz_match == TRUE, !is.na(Top1_SMILES), Top1_SMILES != "",
         !(Alignment_ID %in% confirmed_ids))
cat("NAP tentative (Level 3, excluding confirmed overlap):", nrow(nap_tentative), "\n")
nap_ids <- nap_tentative$Alignment_ID

# ---- Assign every LC feature (of the 155,448) to one tier ----
lc_ids <- data.table(Alignment_ID = nap$Alignment_ID)
lc_ids[, Tier := "Not annotated"]
lc_ids[Alignment_ID %in% confirmed_ids[targets$Tier == "Level 1"], Tier := "Level 1"]
lc_ids[Alignment_ID %in% confirmed_ids[targets$Tier == "Level 2"], Tier := "Level 2"]
lc_ids[Alignment_ID %in% nap_ids, Tier := "Level 3 (NAP)"]

lc_counts <- lc_ids[, .N, by = Tier]
setnames(lc_counts, "N", "n")

# Use the full confirmed-target counts (115 / 204) rather than only those
# matched within the NAP-derived universe, adjusting Not annotated down
# to keep the total at n_lc_total.
n_level1 <- sum(targets$Tier == "Level 1")
n_level2 <- sum(targets$Tier == "Level 2")
n_level3 <- lc_counts$n[lc_counts$Tier == "Level 3 (NAP)"]
n_not_annotated <- n_lc_total - n_level1 - n_level2 - n_level3

lc_counts <- data.table(
  Tier = c("Level 1", "Level 2", "Level 3 (NAP)", "Not annotated"),
  n = c(n_level1, n_level2, n_level3, n_not_annotated)
)
lc_counts$Tier <- factor(lc_counts$Tier,
  levels = c("Level 1", "Level 2", "Level 3 (NAP)", "Not annotated"))
lc_counts <- lc_counts[order(lc_counts$Tier), ]

cat("\nLC-HRMS tier counts:\n"); print(lc_counts)
stopifnot(sum(lc_counts$n) == n_lc_total)

write.csv(lc_counts, "Annotation_ConfidenceLevel_LC.csv", row.names = FALSE)

tier_colors <- c(
  "Level 1"       = "#5A2D75",
  "Level 2"       = "#9476A5",
  "Level 3 (NAP)" = "#CEC0D6",
  "Not annotated" = "#D9D9D6"
)

p_lc <- bar_chart(as.data.frame(lc_counts), "Tier", "n", n_lc_total, tier_colors,
                   paste0("LC-HRMS annotation confidence level distribution (n = ", comma(n_lc_total), ")"))
print(p_lc)
ggsave("Annotation_ConfidenceLevel_LC_bar.png", p_lc, width = 7.5, height = 6, dpi = 300, bg = "white")
cat("Saved: Annotation_ConfidenceLevel_LC_bar.png\n\n")

# ================================================================
# GC-HRMS
# ================================================================
cat("=== GC-HRMS ===\n")

gc_file <- "GC_Level1_5_polished_FULL_final.xlsx"

gc_raw <- read_excel(gc_file)
cat("Column names include 'ID level'?", "ID level" %in% names(gc_raw), "\n")

gc_raw$`ID level` <- as.numeric(gc_raw$`ID level`)

n_gc_total <- nrow(gc_raw)
cat("Total GC features:", n_gc_total, "\n")

gc_counts <- gc_raw %>%
  mutate(Tier = paste0("Level ", `ID level`)) %>%
  count(Tier, name = "n")

gc_counts$Tier <- factor(gc_counts$Tier,
  levels = paste0("Level ", sort(unique(gc_raw$`ID level`))))
gc_counts <- gc_counts[order(gc_counts$Tier), ]

cat("\nGC-HRMS tier counts:\n"); print(gc_counts)
stopifnot(sum(gc_counts$n) == n_gc_total)

write.csv(gc_counts, "Annotation_ConfidenceLevel_GC.csv", row.names = FALSE)

gc_tier_colors <- c(
  "Level 1" = "#A85A1F", "Level 2" = "#BE8357", "Level 3" = "#D4AC8F",
  "Level 4" = "#E5CEBC", "Level 5" = "#F4EAE2"
)

p_gc <- bar_chart(as.data.frame(gc_counts), "Tier", "n", n_gc_total, gc_tier_colors,
                   paste0("GC-HRMS annotation confidence level distribution (n = ", comma(n_gc_total), ")"))
print(p_gc)
ggsave("Annotation_ConfidenceLevel_GC_bar.png", p_gc, width = 7.5, height = 6, dpi = 300, bg = "white")
cat("Saved: Annotation_ConfidenceLevel_GC_bar.png\n")

cat("\n=== DONE ===\n")

