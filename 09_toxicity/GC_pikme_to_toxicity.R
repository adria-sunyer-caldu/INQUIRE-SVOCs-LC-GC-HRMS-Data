# ============================================================================
# INQUIRE — GC-HRMS: DTXSID retrieval, pikme hazard-data filtering, and
# toxicity/EDC characterization, all in one pipeline.
#
# Restricted to ID level 1-2 (confirmed/probable structure) — level 3
# (tentative candidate, real spectral match but RI outside confidence
# interval) tested as a sensitivity check and excluded from the primary
# analysis; see summary doc for the comparison.
# Paired (per-household, main samplers only) indoor/outdoor ratio, all confirmed
# compounds used for every statistic/plot, no concentration prioritization.
# CATMoS GHS grouping dropped (n=2, uninformative).
#
# Steps in this script:
#   1) Retrieve DTXSIDs for the GC-annotated compounds via EPA's ctxR
#      package, using InChIKey (falling back to CAS for anything that
#      doesn't resolve).
#   2) Filter the pikme hazard database (pikme_all.parquet) down to just
#      those DTXSIDs.
#   3) Merge with GC concentration data, OPERA predictions, and pikme
#      hazard scores; compute paired indoor/outdoor ratios; run the
#      statistics and produce the final figures.
# ============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(arrow)   # for reading pikme_all.parquet

if (!requireNamespace("ctxR", quietly = TRUE)) install.packages("ctxR", repos = "https://cloud.r-project.org")
library(ctxR)

# ---- Set your working/data directory here ----
data_dir <- "path/to/data"
setwd(data_dir)

# ============================================================================
# STEP 1: Retrieve DTXSIDs via EPA CTX (ctxR), using InChIKey as the lookup
# key, with CAS as a fallback for anything that doesn't resolve. Using ctxR
# instead of hand-rolled HTTP calls because it has the correct, EPA-
# maintained endpoint and header logic built in.
# ============================================================================

# ---- EPA CTX API key ----
# Get your own free key at https://ctx.ct2b.epa.gov/#/apikey (or the
# current CompTox API registration page). Set it as an environment
# variable named CTX_API_KEY before running this script — never hardcode
# the key itself here, since this script is version-controlled.
api_key <- Sys.getenv("CTX_API_KEY")
if (api_key == "") {
  stop("CTX_API_KEY environment variable is not set. Set it with Sys.setenv(CTX_API_KEY = \"your-key-here\") or in your .Renviron file before running this script.")
}
register_ctx_api_key(key = api_key)
cat("API key registered:", has_ctx_key(), "\n\n")

compounds <- read.csv("GC_159_compounds_for_DTXSID_lookup.csv", stringsAsFactors = FALSE)
cat("Compounds to look up:", nrow(compounds), "\n")

# ---- Single-compound connectivity test before the full batch ----
cat("Testing ctxR connection with one InChIKey:", compounds$InChiKey[1], "\n")
test_result <- tryCatch(
  chemical_equal(word = compounds$InChiKey[1], verbose = TRUE),
  error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL }
)
print(test_result)
if (is.null(test_result) || (is.data.frame(test_result) && nrow(test_result) == 0)) {
  stop("Connectivity test failed — check your API key and network connection before continuing.")
}

# ---- Full batch lookup by InChIKey ----
inchikeys <- compounds$InChiKey[compounds$InChiKey != "" & !is.na(compounds$InChiKey)]
batch_result <- chemical_equal_batch(word_list = inchikeys, verbose = TRUE)
cat("\nMatched", nrow(batch_result), "of", length(inchikeys), "InChIKeys queried\n")

# batch_result is a flat table, one row per matched compound.
# searchValue = the InChIKey that was queried; dtxsid = the result.
dtxsid_map <- as.data.frame(batch_result)[, c("searchValue", "dtxsid")]
names(dtxsid_map) <- c("InChiKey", "dtxsid")
dtxsid_map <- dtxsid_map[!duplicated(dtxsid_map$InChiKey), ]  # in case of duplicate hits

compounds <- compounds %>% left_join(dtxsid_map, by = "InChiKey")
cat("\nDTXSIDs found via InChIKey:", sum(!is.na(compounds$dtxsid)), "/", nrow(compounds), "\n")

# ---- Retry pass: use CAS for anything that didn't resolve via InChIKey ----
failed <- compounds %>% filter(is.na(dtxsid))
if (nrow(failed) > 0) {
  cat(nrow(failed), "compounds did not resolve via InChIKey",
      "(these likely need CAS-based lookup or don't exist in DSSTox)\n")

  retry_cas <- failed %>% filter(!is.na(CAS) & CAS != "")
  cat("\nRetrying", nrow(retry_cas), "of", nrow(failed), "failures via CAS...\n")

  if (nrow(retry_cas) > 0) {
    cas_batch_result <- chemical_equal_batch(word_list = retry_cas$CAS, verbose = TRUE)
    cat("CAS retry matched:", nrow(cas_batch_result), "of", nrow(retry_cas), "\n")

    if (nrow(cas_batch_result) > 0) {
      cas_map <- as.data.frame(cas_batch_result)[, c("searchValue", "dtxsid")]
      names(cas_map) <- c("CAS", "dtxsid_from_cas")
      cas_map <- cas_map[!duplicated(cas_map$CAS), ]

      compounds <- compounds %>%
        left_join(cas_map, by = "CAS") %>%
        mutate(dtxsid = ifelse(is.na(dtxsid), dtxsid_from_cas, dtxsid)) %>%
        select(-dtxsid_from_cas)

      cat("\nDTXSIDs found after CAS retry:", sum(!is.na(compounds$dtxsid)), "/", nrow(compounds), "\n")
    }
  }

  still_failed <- compounds %>% filter(is.na(dtxsid))
  if (nrow(still_failed) > 0) {
    write.csv(still_failed, "GC_dtxsid_lookup_failures.csv", row.names = FALSE)
    cat(nrow(still_failed), "compounds still unresolved after both InChIKey and CAS lookup",
        "— see GC_dtxsid_lookup_failures.csv\n")
  }
}

write.csv(compounds, "Toxicity data GC.csv", row.names = FALSE)
cat("Saved: Toxicity data GC.csv (includes dtxsid column)\n\n")

# ============================================================================
# STEP 2: Filter the pikme hazard database to just the DTXSIDs found above
# ============================================================================
dtxsids <- unique(compounds$dtxsid[!is.na(compounds$dtxsid) & compounds$dtxsid != ""])
cat("Unique valid DTXSIDs to look up in pikme:", length(dtxsids), "\n")

pikme_data <- read_parquet("pikme_all.parquet")
pikme_filtered <- pikme_data %>% filter(dtxsid %in% dtxsids)

write.csv(pikme_filtered, "pikme_filtered_GC.csv", row.names = FALSE)
cat("Saved: pikme_filtered_GC.csv\n")
cat("DTXSIDs requested:", length(dtxsids), "| DTXSIDs matched in pikme:", nrow(pikme_filtered), "\n\n")

# ============================================================================
# STEP 3: Toxicity / EDC characterization
# ============================================================================

# ---------------------------
# 3a) Load GC concentration data + metadata (ID level, SMILES, Feature_ID)
# ---------------------------
gc_raw <- read_excel("GC_Level1_2_polished_FINAL.xlsx")

meta_cols <- c("RT","RT2","Feature_ID","Database match","Formula","InChiKey","Smiles","CAS",
               "IUPAC name","Method","m/z","Ion","MS/MS","MW","Exact Mass","ID level",
               "Componant info","RT/RI LC-MS","RT/RI UoA","RI","LoD","LoQ","Uncertainty",
               "Quantification Method","Notes on Batch Effects","Concentration Type")
meta_cols <- intersect(meta_cols, names(gc_raw))

gc_confirmed <- gc_raw %>% filter(`ID level` %in% c(1, 2))
cat("Restricted to ID level 1-2:", nrow(gc_confirmed), "of", nrow(gc_raw), "total rows\n")
cat("(", sum(gc_raw$`ID level` == 3, na.rm=TRUE), "level-3 tentative candidates excluded — see summary doc for sensitivity check)\n")

sample_cols_all <- setdiff(names(gc_raw), meta_cols)
sample_cols <- sample_cols_all[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", sample_cols_all)]
cat("Indoor/outdoor sample columns kept:", length(sample_cols), "of", length(sample_cols_all), "\n")

gc_confirmed[sample_cols] <- lapply(gc_confirmed[sample_cols], as.numeric)

# ---------------------------
# 3b) PAIRED indoor/outdoor ratio (per household, main samplers only)
# ---------------------------
indoor_outdoor <- ifelse(grepl("_Indoor_", sample_cols), "Indoor", "Outdoor")
house_id <- sub("_(Indoor|Outdoor)_", "_", sample_cols)   # {Country}_{code}

pair_table <- data.frame(sample = sample_cols, house = house_id, io = indoor_outdoor, stringsAsFactors = FALSE)
houses_with_both <- pair_table %>% group_by(house) %>% summarise(n_io = n_distinct(io)) %>% filter(n_io == 2) %>% pull(house)
cat("Houses with both indoor and outdoor samples:", length(houses_with_both), "of", n_distinct(house_id), "\n")

paired_ratio_per_feature <- function(row_values) {
  ratios <- sapply(houses_with_both, function(h) {
    in_col  <- pair_table$sample[pair_table$house == h & pair_table$io == "Indoor"]
    out_col <- pair_table$sample[pair_table$house == h & pair_table$io == "Outdoor"]
    in_val  <- sum(row_values[in_col], na.rm = TRUE)
    out_val <- sum(row_values[out_col], na.rm = TRUE)
    if ((in_val + out_val) == 0) return(NA_real_)
    in_val / (in_val + out_val)
  })
  mean(ratios, na.rm = TRUE)
}

gc_confirmed <- gc_confirmed %>%
  mutate(
    Mean_Conc = rowMeans(across(all_of(sample_cols)), na.rm = TRUE),
    Detection_Freq = rowSums(across(all_of(sample_cols), ~ . > 0), na.rm = TRUE) / length(sample_cols)
  )
sample_matrix <- gc_confirmed[, sample_cols]
gc_confirmed$Paired_Indoor_Ratio <- apply(sample_matrix, 1, paired_ratio_per_feature)
gc_confirmed <- gc_confirmed %>% filter(!is.na(Paired_Indoor_Ratio))

# ---------------------------
# 3c) Merge crosswalk (Feature_ID <-> dtxsid) + OPERA + pikme t_human_score
# ---------------------------
crosswalk <- compounds %>% select(Feature_ID, dtxsid)   # from Step 1, already in memory
opera <- read.csv("OPERA_GC_results.csv", stringsAsFactors = FALSE) %>%
  rename(Feature_ID = MoleculeID) %>%
  select(Feature_ID, LogKOA_pred, AD_KOA, CERAPP_Bind_pred, AD_CERAPP_Bind,
         CoMPARA_Bind_pred, AD_CoMPARA_Bind, CATMoS_LD50_pred, CATMoS_GHS_pred, AD_CATMoS)
pikme <- pikme_filtered %>%   # from Step 2, already in memory
  select(dtxsid, t_human_score) %>% mutate(t_human_score = as.numeric(t_human_score))

results <- gc_confirmed %>%
  left_join(crosswalk, by = "Feature_ID") %>%
  left_join(opera, by = "Feature_ID") %>%
  left_join(pikme, by = "dtxsid") %>%
  mutate(
    EDC_flag = ifelse(CERAPP_Bind_pred == 1 | CoMPARA_Bind_pred == 1, 1, 0),
    AD_EDC_lab = ifelse(AD_CERAPP_Bind == 1 | AD_CoMPARA_Bind == 1, "In domain", "Out of domain")
  )

cat("\nFinal GC working set:", nrow(results), "confirmed (ID level 1-2) compounds\n")
write.csv(results, "GC_toxicity_results_full.csv", row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE)) writexl::write_xlsx(results, "GC_Toxicity_FULL_RESULTS.xlsx")

# ============================================================================
# STEP 4: BASELINE-RELATIVE STATISTICS + WILCOXON CONFIRMATION
# ============================================================================
n_all <- sum(!is.na(results$Paired_Indoor_Ratio))
n_indoor_all <- sum(results$Paired_Indoor_Ratio > 0.5, na.rm = TRUE)
baseline_pct <- 100 * n_indoor_all / n_all

bt_all <- binom.test(n_indoor_all, n_all, p = 0.5)
wt_all <- wilcox.test(results$Paired_Indoor_Ratio, mu = 0.5, alternative = "greater")

cat("\n==== GC whole-dataset baseline ====\n")
cat("Indoor-dominant:", n_indoor_all, "/", n_all, "=", round(baseline_pct, 1), "%\n")
cat("Binomial sign test vs 50%: p =", format.pval(bt_all$p.value, digits = 3), "\n")
cat("Wilcoxon signed-rank vs 0.5: V =", wt_all$statistic, ", p =", format.pval(wt_all$p.value, digits = 3), "\n\n")

test_group_vs_baseline <- function(label, subset_df) {
  n_g <- sum(!is.na(subset_df$Paired_Indoor_Ratio))
  if (n_g == 0) { cat(sprintf("%-18s n=0, skipped\n", label)); return(NULL) }
  n_ind <- sum(subset_df$Paired_Indoor_Ratio > 0.5, na.rm = TRUE)
  pct <- 100 * n_ind / n_g
  bt <- binom.test(n_ind, n_g, p = 0.5)
  wt <- tryCatch(wilcox.test(subset_df$Paired_Indoor_Ratio, mu = 0.5, alternative = "greater"), error = function(e) NULL)
  wp <- if (!is.null(wt)) format.pval(wt$p.value, digits = 3) else "n/a"
  cat(sprintf("%-18s n=%3d  indoor=%3d/%3d (%5.1f%%)  binomial p=%-9s wilcoxon p=%s\n",
      label, n_g, n_ind, n_g, pct, format.pval(bt$p.value, digits = 3), wp))
  data.frame(Group = label, n = n_g, n_indoor = n_ind, pct_indoor = pct,
             binomial_p = format.pval(bt$p.value, digits = 3), wilcoxon_p = wp)
}

cat("==== GC: toxic/EDC subsets tested against 50%, on their own (not vs baseline population) ====\n")
summ <- bind_rows(
  test_group_vs_baseline("t_human_score>=4", results %>% filter(t_human_score >= 4)),
  test_group_vs_baseline("EDC-flagged", results %>% filter(EDC_flag == 1))
)
write.csv(summ, "GC_toxic_subset_indoor_tests.csv", row.names = FALSE)

abs_counts <- bind_rows(
  data.frame(Group = "t_human_score>=4", Indoor = sum(results$t_human_score>=4 & results$Paired_Indoor_Ratio>0.5, na.rm=TRUE),
             Outdoor = sum(results$t_human_score>=4 & results$Paired_Indoor_Ratio<=0.5, na.rm=TRUE)),
  data.frame(Group = "EDC-flagged", Indoor = sum(results$EDC_flag==1 & results$Paired_Indoor_Ratio>0.5, na.rm=TRUE),
             Outdoor = sum(results$EDC_flag==1 & results$Paired_Indoor_Ratio<=0.5, na.rm=TRUE))
)
write.csv(abs_counts, "GC_absolute_indoor_outdoor_counts.csv", row.names=FALSE)

confirmatory_stats <- data.frame(
  Test = c("Binomial sign test (indoor-dominant vs 50%) — ALL confirmed", "Wilcoxon signed-rank (ratio vs 0.5) — ALL confirmed"),
  Statistic = c(paste0(n_indoor_all, "/", n_all), paste0("V=", wt_all$statistic)),
  p_value = c(format.pval(bt_all$p.value, digits = 3), format.pval(wt_all$p.value, digits = 3)),
  Interpretation = c("Highly significant indoor skew, whole dataset", "Confirms indoor skew (magnitude-sensitive), whole dataset")
)
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(Whole_dataset_stats = confirmatory_stats, Toxic_subset_tests = summ, Absolute_counts = abs_counts),
    "GC_Toxicity_STATISTICS_SUMMARY.xlsx"
  )
}

# ---------------------------
# 4b) GC_Results_Summary.xlsx — same structure as the LC script's
# Targets_Toxicity_STATISTICS_SUMMARY.xlsx, for easy side-by-side comparison
# ---------------------------
whole_dataset_stats_LCstyle <- data.frame(
  Metric = c("n compounds", "n indoor-dominant", "pct indoor-dominant",
             "Binomial p vs 50%", "Wilcoxon p vs 0.5"),
  Value = c(n_all, n_indoor_all, round(baseline_pct, 1),
            format.pval(bt_all$p.value, digits = 3), format.pval(wt_all$p.value, digits = 3))
)

toxic_edc_subset_stats_LCstyle <- bind_rows(
  results %>% filter(t_human_score >= 4) %>%
    summarise(Subset = "t_human_score>=4", n = n(),
              n_indoor = sum(Paired_Indoor_Ratio > 0.5, na.rm = TRUE),
              pct_indoor = round(100 * n_indoor / n, 1),
              Binomial_p_vs_50pct = format.pval(binom.test(n_indoor, n, p = 0.5)$p.value, digits = 3),
              Wilcoxon_p_vs_0.5 = format.pval(tryCatch(wilcox.test(Paired_Indoor_Ratio, mu = 0.5, alternative = "greater")$p.value, error = function(e) NA), digits = 3)),
  results %>% filter(EDC_flag == 1) %>%
    summarise(Subset = "EDC-flagged", n = n(),
              n_indoor = sum(Paired_Indoor_Ratio > 0.5, na.rm = TRUE),
              pct_indoor = round(100 * n_indoor / n, 1),
              Binomial_p_vs_50pct = format.pval(binom.test(n_indoor, n, p = 0.5)$p.value, digits = 3),
              Wilcoxon_p_vs_0.5 = format.pval(tryCatch(wilcox.test(Paired_Indoor_Ratio, mu = 0.5, alternative = "greater")$p.value, error = function(e) NA), digits = 3))
)

full_compound_data_LCstyle <- results %>% select(-all_of(sample_cols))

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      "1_Whole_dataset_stats" = whole_dataset_stats_LCstyle,
      "2_Toxic_EDC_subset_stats" = toxic_edc_subset_stats_LCstyle,
      "3_Absolute_counts" = abs_counts,
      "4_Full_compound_data" = full_compound_data_LCstyle
    ),
    "GC_Results_Summary.xlsx"
  )
  cat("\nSaved: GC_Results_Summary.xlsx (LC-matching 4-sheet format)\n")
}

# ============================================================================
# STEP 5: FIGURES — all titled "GC-HRMS" explicitly, matching LC-HRMS
# script aesthetics exactly
# ============================================================================

## Fig 0: whole-dataset baseline result, its own figure
whole_df <- data.frame(Group = "All confirmed GC compounds", n = n_all, pct_indoor = baseline_pct,
                        binomial_p = format.pval(bt_all$p.value, digits = 3),
                        wilcoxon_p = format.pval(wt_all$p.value, digits = 3))
whole_df$label <- paste0(round(whole_df$pct_indoor, 1), "% indoor-dominant (n=", whole_df$n, ")\n",
                          "Binomial p=", whole_df$binomial_p, "   Wilcoxon p=", whole_df$wilcoxon_p)
p0 <- ggplot(whole_df, aes(x = Group, y = pct_indoor)) +
  geom_col(fill = "#0F6E56", alpha = 0.85, width = 0.5) +
  geom_text(aes(label = label), hjust = -0.03, size = 3.6, lineheight = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#D62828") +
  coord_flip(clip = "off") + ylim(0, 145) +
  labs(title = "GC-HRMS: whole confirmed-compound dataset is indoor-dominant",
       subtitle = "All confirmed (ID level 1-2) sVOCs, tested against 50%",
       x = NULL, y = "% indoor-dominant") +
  theme_minimal(base_size = 13) + theme(plot.margin = margin(5, 130, 5, 5), axis.text.y = element_blank())
print(p0)
ggsave("Fig_GC_whole_dataset_indoor_pct.png", p0, width = 10, height = 3.5, dpi = 300)

## Fig 1: toxic/EDC subsets, % indoor + significance, tested against 50% directly
summ$Group <- factor(summ$Group, levels = rev(summ$Group))
summ$sig_label <- paste0(round(summ$pct_indoor, 1), "%\n(n=", summ$n, ")\nbinom p=", summ$binomial_p, "\nwilcox p=", summ$wilcoxon_p)
p1 <- ggplot(summ, aes(x = Group, y = pct_indoor)) +
  geom_col(fill = "#A85A1F", alpha = 0.85) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#D62828") +
  coord_flip(clip = "off") + ylim(0, 100) +
  labs(title = "GC-HRMS: are toxic/EDC compounds themselves indoor-skewed?",
       subtitle = "Each subset tested on its own against 50% (dashed line) — not compared to the rest of the dataset; ID level 1-2 only",
       x = NULL, y = "% indoor-dominant within the subset (Paired_Indoor_Ratio > 0.5)") +
  theme_minimal(base_size = 12) + theme(plot.margin = margin(5, 90, 5, 5))
print(p1)
ggsave("Fig_GC_toxic_subset_indoor_pct.png", p1, width = 12, height = 6, dpi = 600)

## Fig 2: t_human_score volcano
pd <- results %>%
  filter(!is.na(Paired_Indoor_Ratio)) %>%
  mutate(Log2FC = log2((Paired_Indoor_Ratio + 1e-3) / (1 - Paired_Indoor_Ratio + 1e-3)),
         DetFreqPct = Detection_Freq * 100,
         Hazard_cat = case_when(
           !is.na(t_human_score) & t_human_score == 5 ~ "Score 5",
           !is.na(t_human_score) & t_human_score == 4 ~ "Score 4",
           TRUE ~ "Unscored/lower"
         ))
pd_haz45 <- pd %>% filter(Hazard_cat %in% c("Score 4", "Score 5"))
pd_haz45 <- pd_haz45 %>%
  mutate(priority_score = DetFreqPct * pmax(Log2FC, 0)) %>%
  arrange(desc(Hazard_cat == "Score 5"), desc(priority_score))
cat("\n=== GC Regulatory priority ranking ===\n")
print(pd_haz45 %>% select(any_of(c("Name", "Database match")), Hazard_cat, DetFreqPct, Log2FC, priority_score), n = Inf)

## Threshold is a hard filter on the data itself (not a soft axis zoom), so
## nothing below 2% DF or outside the x-range can ever be plotted or
## labeled, regardless of ggplot's automatic axis expansion.
pd_haz45_plot <- pd_haz45 %>%
  filter(DetFreqPct >= 2, Log2FC >= -2.5, Log2FC <= 5)

p2 <- ggplot(pd_haz45_plot, aes(Log2FC, DetFreqPct, color = Hazard_cat)) +
  geom_point(size = 3, alpha = 0.75) +
  geom_text_repel(aes(label = `Database match`), color = "black", size = 4,
                   box.padding = 0.5, point.padding = 0.3, force = 1, max.iter = 20000,
                   max.overlaps = Inf, segment.size = 0.3, segment.color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Score 5" = "#d62828", "Score 4" = "#fcbf49"), name = "t_human_score") +
  coord_cartesian(xlim = c(-2.5, 5)) +
  theme_minimal(base_size = 11) +
  labs(title = "GC-HRMS: paired indoor/outdoor ratio vs. human hazard score",
       subtitle = paste0("Score 4-5 compounds only, n=", nrow(pd_haz45), " of ", nrow(results), " confirmed; all visible compounds labeled"),
       x = "Paired indoor/outdoor ratio (log2)", y = "Detection frequency (%)")
print(p2)
ggsave("Volcano_GC_hazard_score.png", p2, width = 8, height = 9, dpi = 600)

cat("\nGC Score 5 labeled (", sum(pd_haz45_plot$Hazard_cat == "Score 5"), "):\n")
print(pd_haz45_plot$`Database match`[pd_haz45_plot$Hazard_cat == "Score 5"])
cat("\nGC Score 4 labeled (", sum(pd_haz45_plot$Hazard_cat == "Score 4"), "):\n")
print(pd_haz45_plot$`Database match`[pd_haz45_plot$Hazard_cat == "Score 4"])

## Fig 3: EDC volcano, colored by indoor/outdoor direction among EDC-flagged
pd2 <- results %>%
  filter(!is.na(EDC_flag) & !is.na(Paired_Indoor_Ratio)) %>%
  mutate(Log2FC = log2((Paired_Indoor_Ratio + 1e-3) / (1 - Paired_Indoor_Ratio + 1e-3)),
         DetFreqPct = Detection_Freq * 100,
         EDC_direction = case_when(
           EDC_flag == 1 & Paired_Indoor_Ratio > 0.5  ~ "EDC, indoor-skewed",
           EDC_flag == 1 & Paired_Indoor_Ratio <= 0.5 ~ "EDC, outdoor-skewed",
           TRUE ~ "Non-EDC"
         ))
n_edc_indoor <- sum(pd2$EDC_direction == "EDC, indoor-skewed")
n_edc_outdoor <- sum(pd2$EDC_direction == "EDC, outdoor-skewed")
p3 <- ggplot(pd2, aes(Log2FC, DetFreqPct, color = EDC_direction, alpha = AD_EDC_lab)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("EDC, indoor-skewed" = "#007e5d", "EDC, outdoor-skewed" = "#e8c547", "Non-EDC" = "grey80"),
                      name = "EDC status & direction") +
  scale_alpha_manual(values = c("In domain" = 0.85, "Out of domain" = 0), guide = "none") +
  theme_minimal(base_size = 11) +
  labs(title = "GC-HRMS: paired indoor/outdoor ratio vs. DF; colored by EDC activity",
       subtitle = paste0(n_edc_indoor, " EDC-flagged indoor-skewed, ", n_edc_outdoor, " outdoor-skewed"),
       x = "Log2 (Paired Indoor Ratio / (1 - Paired Indoor Ratio))", y = "Detection frequency (%)")
print(p3)
ggsave("Volcano_GC_EDC_ALL.png", p3, width = 9.5, height = 6.5, dpi = 600)

## Fig 4: LogKOA vs concentration, DF as color/size (matching LC combined-legend style)
pd_koa <- results %>% filter(!is.na(LogKOA_pred)) %>% mutate(AD_lab = ifelse(AD_KOA == 1, "In domain", "Out of domain"))
p4 <- ggplot(pd_koa, aes(x = Mean_Conc, y = LogKOA_pred)) +
  geom_point(aes(color = Detection_Freq, size = Detection_Freq, alpha = AD_lab)) +
  scale_x_log10() +
  scale_color_gradient(low = "#dad7cd", high = "#21b0fe",
                       name = "Detection Frequency",
                       limits = c(0, 1)) +
  scale_size(range = c(1.5, 8),
             name = "Detection Frequency",
             limits = c(0, 1)) +
  guides(color = guide_legend(), size = guide_legend()) +
  scale_alpha_manual(values = c("In domain" = 0.9, "Out of domain" = 0), guide = "none") +
  theme_minimal(base_size = 12) +
  labs(title = "GC-HRMS: LogKOA vs. mean area",
       subtitle = paste0("n=", nrow(pd_koa), ", level 1-2 only"),
       x = "Mean Area", y = "LogKOA (Air Partitioning)")
print(p4)
ggsave("Scatter_GC_logKoa_conc_DF.png", p4, width = 9.5, height = 7, dpi = 600)

## Fig 5: absolute indoor vs outdoor counts, hazard/EDC groups
abs_long <- abs_counts %>% pivot_longer(cols = c(Indoor, Outdoor), names_to = "Direction", values_to = "n")
p5 <- ggplot(abs_long, aes(x = Group, y = n, fill = Direction)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c(Indoor = "#B03A2E", Outdoor = "#1F618D")) +
  labs(title = "GC-HRMS: absolute count of hazardous/EDC-flagged compounds, indoor vs outdoor",
       subtitle = "More hazardous compounds found indoor than outdoor in raw terms",
       x = NULL, y = "Number of compounds") +
  theme_minimal(base_size = 12)
print(p5)
ggsave("Fig_GC_absolute_counts_hazard.png", p5, width = 8.5, height = 6, dpi = 600)

cat("\n==== DONE — GC-HRMS pipeline complete (DTXSID retrieval -> pikme -> toxicity analysis) ====\n")
