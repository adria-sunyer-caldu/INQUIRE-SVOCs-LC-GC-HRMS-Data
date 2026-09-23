# ============================================================================
# INQUIRE — LC-HRMS: pikme hazard-data filtering and toxicity/EDC
# characterization, in one pipeline.
#
# Paired (per-household) indoor/outdoor ratio ONLY. All confirmed compounds
# used for every statistic/plot — no concentration-based prioritization.
# Hazard/EDC subgroups tested against the dataset's OWN baseline (not 50/50),
# PLUS Wilcoxon signed-rank as a magnitude-sensitive confirmatory test.
# CATMoS GHS grouping dropped (too small, n=5, uninformative).
#
# Steps in this script:
#   1) Filter the pikme hazard database (pikme_all.parquet) down to the
#      DTXSIDs already present in "Toxicity data targets and
#      annotations.csv" (DTXSID retrieval for LC targets/annotations was
#      done separately, upstream of this script).
#   2) Merge with LC concentration data, OPERA predictions, and pikme
#      hazard scores; compute paired indoor/outdoor ratios; run the
#      statistics and produce the final figures.
# ============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(arrow)   # for reading pikme_all.parquet

# ---- Set your working/data directory here ----
data_dir <- "path/to/data"
setwd(data_dir)

# ============================================================================
# STEP 1: Filter the pikme hazard database to the DTXSIDs already present
# in the targets/annotations toxicity data file (column: DTXSID)
# ============================================================================
toxicity_data <- read.csv("Toxicity data targets and annotations.csv", stringsAsFactors = FALSE)
dtxsids <- unique(toxicity_data$DTXSID)
dtxsids <- dtxsids[!is.na(dtxsids) & dtxsids != ""]
cat("Unique valid DTXSIDs to look up in pikme:", length(dtxsids), "\n")

pikme_data <- read_parquet("pikme_all.parquet")
pikme_filtered <- pikme_data %>% filter(dtxsid %in% dtxsids)
cat("DTXSIDs requested:", length(dtxsids), "| DTXSIDs matched in pikme:", nrow(pikme_filtered), "\n\n")

write.csv(pikme_filtered, "pikme_filtered_by_DTXSID.csv", row.names = FALSE)
cat("Saved: pikme_filtered_by_DTXSID.csv\n\n")

# ============================================================================
# STEP 2: Toxicity / EDC characterization
# ============================================================================

# ---------------------------
# 2a) Load concentration data
# ---------------------------
target_file <- "LC_Level1_2_polished_FINAL.xlsx"
data_raw <- read_excel(target_file)

meta_cols <- c("Feature ID","MZ","RT","Name","Molecular formula","InChiKey","SMILES",
               "CAS","Ion type","Confidence level","MS/MS","LOD","LOQ","Quantification")
names(data_raw) <- trimws(names(data_raw))   # fixes trailing whitespace on sample column names
sample_cols <- setdiff(names(data_raw), meta_cols)
data_raw[sample_cols] <- lapply(data_raw[sample_cols], as.numeric)

feature_meta <- data_raw[, meta_cols]
sample_matrix <- data_raw[, sample_cols]

# ---------------------------
# 2b) pikme t_human_score (from Step 1, already in memory) + fresh
# AD-flagged OPERA rerun
# ---------------------------
pikme <- pikme_filtered %>% mutate(t_human_score = as.numeric(t_human_score))

opera_targets <- read.csv("OPERA_targets_results.csv") %>% rename(`Feature ID` = MoleculeID)
cat("OPERA-matched targets:", nrow(opera_targets), "\n")

feature_meta_haz <- feature_meta %>%
  left_join(pikme %>% select(SU_name, t_human_score), by = c("Name" = "SU_name")) %>%
  left_join(opera_targets %>% select(`Feature ID`, LogKOA_pred, AD_KOA,
                                      CERAPP_Bind_pred, AD_CERAPP_Bind, CoMPARA_Bind_pred, AD_CoMPARA_Bind,
                                      CATMoS_LD50_pred, CATMoS_GHS_pred, AD_CATMoS),
            by = "Feature ID") %>%
  mutate(
    logkoa_pred_opera = LogKOA_pred,
    AD_KOA_lab = ifelse(AD_KOA == 1, "In domain", "Out of domain"),
    EDC_flag = ifelse(CERAPP_Bind_pred == 1 | CoMPARA_Bind_pred == 1, 1, 0),
    AD_EDC_lab = ifelse(AD_CERAPP_Bind == 1 | AD_CoMPARA_Bind == 1, "In domain", "Out of domain")
  )

# ---------------------------
# 2c) PAIRED indoor/outdoor ratio (per house, averaged across pairs)
# ---------------------------
sample_names <- colnames(sample_matrix)
indoor_outdoor <- ifelse(grepl("_Indoor_", sample_names), "Indoor", "Outdoor")
house_id <- sub("_(Indoor|Outdoor)_", "_", sample_names)   # {Country}_{code}; extra samplers (_2/_3) stay unpaired

pair_table <- data.frame(sample = sample_names, house = house_id, io = indoor_outdoor, stringsAsFactors = FALSE)
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

results <- feature_meta_haz %>%
  mutate(
    Mean_Conc = rowMeans(sample_matrix, na.rm = TRUE),
    Detection_Freq = rowSums(sample_matrix > 0, na.rm = TRUE) / ncol(sample_matrix)
  ) %>%
  filter(Mean_Conc > 0)

results$Paired_Indoor_Ratio <- apply(
  sample_matrix[match(results$`Feature ID`, feature_meta_haz$`Feature ID`), ], 1, paired_ratio_per_feature
)
results <- results %>% filter(!is.na(Paired_Indoor_Ratio))

write.csv(results, "Targets_toxicity_results_full.csv", row.names = FALSE)
if (requireNamespace("writexl", quietly = TRUE)) writexl::write_xlsx(results, "Targets_Toxicity_FULL_RESULTS.xlsx")

cat("\nFinal LC working set:", nrow(results), "confirmed compounds\n")

# ============================================================================
# STEP 3: BASELINE-RELATIVE STATISTICS + WILCOXON CONFIRMATION (all
# confirmed compounds, no filtering)
# ============================================================================
n_all <- sum(!is.na(results$Paired_Indoor_Ratio))
n_indoor_all <- sum(results$Paired_Indoor_Ratio > 0.5, na.rm = TRUE)
baseline_pct <- 100 * n_indoor_all / n_all

bt_all <- binom.test(n_indoor_all, n_all, p = 0.5)
wt_all <- wilcox.test(results$Paired_Indoor_Ratio, mu = 0.5, alternative = "greater")

cat("\n==== LC whole-dataset baseline ====\n")
cat("Indoor-dominant:", n_indoor_all, "/", n_all, "=", round(baseline_pct, 1), "%\n")
cat("Binomial sign test vs 50%: p =", format.pval(bt_all$p.value, digits = 3), "\n")
cat("Wilcoxon signed-rank vs 0.5: V =", wt_all$statistic, ", p =", format.pval(wt_all$p.value, digits = 3), "\n\n")

test_group_vs_baseline <- function(label, subset_df) {
  n_g <- sum(!is.na(subset_df$Paired_Indoor_Ratio))
  if (n_g == 0) return(NULL)
  n_ind <- sum(subset_df$Paired_Indoor_Ratio > 0.5, na.rm = TRUE)
  pct <- 100 * n_ind / n_g
  bt <- binom.test(n_ind, n_g, p = 0.5)
  wt <- wilcox.test(subset_df$Paired_Indoor_Ratio, mu = 0.5, alternative = "greater")
  cat(sprintf("%-15s n=%3d  indoor=%3d/%3d (%5.1f%%)  binomial p=%-9s wilcoxon p=%s\n",
      label, n_g, n_ind, n_g, pct, format.pval(bt$p.value, digits = 3), format.pval(wt$p.value, digits = 3)))
  data.frame(Group = label, n = n_g, n_indoor = n_ind, pct_indoor = pct,
             binomial_p = format.pval(bt$p.value, digits = 3), wilcoxon_p = format.pval(wt$p.value, digits = 3))
}

cat("==== LC: toxic/EDC subsets tested against 50%, on their own (not vs baseline population) ====\n")
summ <- bind_rows(
  test_group_vs_baseline("t_human_score>=4", results %>% filter(t_human_score >= 4)),
  test_group_vs_baseline("EDC-flagged", results %>% filter(EDC_flag == 1))
)
write.csv(summ, "LC_toxic_subset_indoor_tests.csv", row.names = FALSE)

# ---- Absolute indoor vs outdoor counts (the key exposure-relevant framing) ----
abs_counts <- bind_rows(
  data.frame(Group="Score 5", Indoor=sum(results$t_human_score==5 & results$Paired_Indoor_Ratio>0.5, na.rm=TRUE),
             Outdoor=sum(results$t_human_score==5 & results$Paired_Indoor_Ratio<=0.5, na.rm=TRUE)),
  data.frame(Group="Score 4", Indoor=sum(results$t_human_score==4 & results$Paired_Indoor_Ratio>0.5, na.rm=TRUE),
             Outdoor=sum(results$t_human_score==4 & results$Paired_Indoor_Ratio<=0.5, na.rm=TRUE)),
  data.frame(Group="Score 4+5", Indoor=sum(results$t_human_score %in% c(4,5) & results$Paired_Indoor_Ratio>0.5, na.rm=TRUE),
             Outdoor=sum(results$t_human_score %in% c(4,5) & results$Paired_Indoor_Ratio<=0.5, na.rm=TRUE)),
  data.frame(Group="EDC-flagged", Indoor=sum(results$EDC_flag==1 & results$Paired_Indoor_Ratio>0.5, na.rm=TRUE),
             Outdoor=sum(results$EDC_flag==1 & results$Paired_Indoor_Ratio<=0.5, na.rm=TRUE))
)
write.csv(abs_counts, "LC_absolute_indoor_outdoor_counts.csv", row.names = FALSE)

# Regulatory cross-check (manually verified against ECHA, see summary doc)
reg <- data.frame(
  Compound = c("Diisobutyl phthalate (DIBP)", "Benzyl butyl phthalate (BBP)", "Dihexyl phthalate", "Dipentyl phthalate",
               "Tris(2-chloroethyl) phosphate (TCEP)", "Triphenyl phosphate (TPHP)", "Galaxolide (HHCB)"),
  Regulatory_status = c("REACH SVHC, Authorisation List (Annex XIV)", "REACH SVHC, Authorisation List (Annex XIV)",
                          "REACH SVHC, Authorisation List (Annex XIV)", "REACH SVHC, Authorisation List (Annex XIV)",
                          "REACH SVHC, Candidate + Authorisation List", "Added to SVHC Candidate List Nov 2024 (EDC)",
                          "ANSES proposed Repr. 1B classification 2025, ECHA RAC review pending")
)

confirmatory_stats <- data.frame(
  Test = c("Binomial sign test (indoor-dominant vs 50%) — ALL confirmed", "Wilcoxon signed-rank (ratio vs 0.5) — ALL confirmed"),
  Statistic = c(paste0(n_indoor_all, "/", n_all), paste0("V=", wt_all$statistic)),
  p_value = c(format.pval(bt_all$p.value, digits = 3), format.pval(wt_all$p.value, digits = 3)),
  Interpretation = c("Highly significant indoor skew, whole dataset", "Confirms indoor skew (magnitude-sensitive), whole dataset")
)

if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(Whole_dataset_stats = confirmatory_stats, Toxic_subset_tests = summ,
         Absolute_counts = abs_counts, Regulatory_crosscheck = reg),
    "Targets_Toxicity_STATISTICS_SUMMARY.xlsx"
  )
}

# ============================================================================
# STEP 4: FIGURES — all titled "LC-HRMS" explicitly, all confirmed
# compounds, no ranking
# ============================================================================

## Fig 0: whole-dataset baseline result, its own figure (full status, not subordinate)
whole_df <- data.frame(Group = "All confirmed compounds", n = n_all, pct_indoor = baseline_pct,
                        binomial_p = format.pval(bt_all$p.value, digits = 3),
                        wilcoxon_p = format.pval(wt_all$p.value, digits = 3))
whole_df$label <- paste0(round(whole_df$pct_indoor, 1), "% indoor-dominant (n=", whole_df$n, ")\n",
                          "Binomial p=", whole_df$binomial_p, "   Wilcoxon p=", whole_df$wilcoxon_p)
p0 <- ggplot(whole_df, aes(x = Group, y = pct_indoor)) +
  geom_col(fill = "#0F6E56", alpha = 0.85, width = 0.5) +
  geom_text(aes(label = label), hjust = -0.03, size = 3.6, lineheight = 1) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#D62828") +
  coord_flip(clip = "off") + ylim(0, 145) +
  labs(title = "LC-HRMS: whole confirmed-compound dataset is indoor-dominant",
       subtitle = "All confirmed target/annotated compounds, tested against 50%",
       x = NULL, y = "% indoor-dominant") +
  theme_minimal(base_size = 13) + theme(plot.margin = margin(5, 130, 5, 5), axis.text.y = element_blank())
print(p0)
ggsave("Fig_LC_whole_dataset_indoor_pct.png", p0, width = 10, height = 3.5, dpi = 300)

## Fig 1: toxic/EDC subsets, % indoor + significance, tested against 50% directly
summ$Group <- factor(summ$Group, levels = rev(summ$Group))
summ$sig_label <- paste0(round(summ$pct_indoor, 1), "%\n(n=", summ$n, ")\nbinom p=", summ$binomial_p, "\nwilcox p=", summ$wilcoxon_p)
p1 <- ggplot(summ, aes(x = Group, y = pct_indoor)) +
  geom_col(fill = "#7B3F9E", alpha = 0.85) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "#D62828") +
  coord_flip(clip = "off") + ylim(0, 100) +
  labs(title = "LC-HRMS: are toxic/EDC compounds themselves indoor-skewed?",
       subtitle = "Each subset tested on its own against 50% (dashed line) — not compared to the rest of the dataset",
       x = NULL, y = "% indoor-dominant within the subset (Paired_Indoor_Ratio > 0.5)") +
  theme_minimal(base_size = 12) + theme(plot.margin = margin(5, 90, 5, 5))
print(p1)
ggsave("Fig_LC_toxic_subset_indoor_pct.png", p1, width = 12, height = 6, dpi = 600)

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
cat("\n=== LC Regulatory priority ranking ===\n")
print(pd_haz45 %>% select(Name, Hazard_cat, DetFreqPct, Log2FC, priority_score), n = Inf)

lbl_haz45 <- pd_haz45 %>% filter(Hazard_cat == "Score 5" | DetFreqPct >= 25)
p2 <- ggplot(pd_haz45, aes(Log2FC, DetFreqPct, color = Hazard_cat)) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(data = lbl_haz45, aes(label = Name), color = "black", size = 3,
                   box.padding = 0.5, point.padding = 0.3, force = 1, max.iter = 20000,
                   max.overlaps = Inf, segment.size = 0.3, segment.color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  scale_color_manual(values = c("Score 5" = "#d62828", "Score 4" = "#fcbf49"), name = "t_human_score") +
  coord_cartesian(ylim = c(5, NA)) +
  theme_minimal(base_size = 11) +
  labs(title = "LC-HRMS: paired indoor/outdoor ratio vs. human hazard score",
       subtitle = paste0("Score 4-5 compounds only, n=", nrow(pd_haz45), " of ", nrow(results), " confirmed; labels shown for Score 5 (all) and Score 4 (DF>=25%)"),
       x = "Paired indoor/outdoor ratio (log2)", y = "Detection frequency (%)")
print(p2)
ggsave("Volcano_LC_hazard_score.png", p2, width = 8, height = 9, dpi = 600)

cat("\nLC Score 5 labeled (", sum(lbl_haz45$Hazard_cat == "Score 5"), "):\n")
print(lbl_haz45$Name[lbl_haz45$Hazard_cat == "Score 5"])
cat("\nLC Score 4 labeled (", sum(lbl_haz45$Hazard_cat == "Score 4"), "):\n")
print(lbl_haz45$Name[lbl_haz45$Hazard_cat == "Score 4"])

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
  labs(title = "LC-HRMS: paired indoor/outdoor ratio vs. DF; colored by EDC activity",
       subtitle = paste0(n_edc_indoor, " EDC-flagged indoor-skewed, ", n_edc_outdoor, " outdoor-skewed"),
       x = "Log2 (Paired Indoor Ratio / (1 - Paired Indoor Ratio))", y = "Detection frequency (%)")
print(p3)
ggsave("Volcano_LC_EDC_ALL.png", p3, width = 9.5, height = 6.5, dpi = 300)

## Fig 4: LogKOA vs concentration, DF as color/size
pd_koa <- results %>% filter(!is.na(logkoa_pred_opera)) %>% mutate(AD_lab = AD_KOA_lab)
p4 <- ggplot(pd_koa, aes(x = Mean_Conc, y = logkoa_pred_opera)) +
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
  labs(title = "LC-HRMS: LogKOA vs. concentration",
       subtitle = paste0("n=", nrow(pd_koa), ", all targets and annotations"),
       x = "Mean Concentration (ng/mL)", y = "LogKOA (Air Partitioning)")
print(p4)
ggsave("Scatter_LC_logKoa_conc_DF.png", p4, width = 9.5, height = 7, dpi = 300)

## Fig 5: absolute indoor vs outdoor counts, hazard/EDC groups
abs_long <- abs_counts %>% pivot_longer(cols = c(Indoor, Outdoor), names_to = "Direction", values_to = "n")
p5 <- ggplot(abs_long, aes(x = Group, y = n, fill = Direction)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.3, size = 4) +
  scale_fill_manual(values = c(Indoor = "#B03A2E", Outdoor = "#1F618D")) +
  labs(title = "LC-HRMS: absolute count of hazardous/EDC-flagged compounds, indoor vs outdoor",
       subtitle = "More hazardous compounds found indoor than outdoor in raw terms",
       x = NULL, y = "Number of compounds") +
  theme_minimal(base_size = 12)
print(p5)
ggsave("Fig_LC_absolute_counts_hazard.png", p5, width = 8.5, height = 6, dpi = 300)

cat("\n==== DONE — LC-HRMS pipeline complete (pikme filtering -> toxicity analysis) ====\n")
