# ================================
# Chemical profiling analysis — NAP-annotated features
# Part A: original plot set (Indoor vs Outdoor overlay)
# Part B: same plots, faceted Indoor/Outdoor where it adds clarity
# Filter: excludes formulas containing As, Ge, Te.
# ElementClass resolves individual halogens (Cl/F/Br); >1 type -> MultiHal.
# Van Krevelen capped to conventional bounds (O/C<=1.2, H/C<=3).
# Statistics: Wilcoxon for DBE, O/C, H/C, AI_mod, saved with
# Benjamini-Hochberg adjusted p-values in Stats_summary.csv.
# ================================

library(dplyr)
library(stringr)
library(ggplot2)
library(readr)
library(tidyr)
library(rcdk)

setwd("path/to/data")

# ------------------------
# Load data
# ------------------------
data <- read_csv("NAP_xref_all_features.csv")
cat("Total features in NAP xref:", nrow(data), "\n")

data <- data %>% filter(mz_match == TRUE, !is.na(Top1_SMILES), Top1_SMILES != "")
cat("Annotated (mz_match=TRUE, has SMILES):", nrow(data), "\n")

# ------------------------
# Indoor vs Outdoor — paired ratio Category only
# ------------------------
data <- data %>%
  mutate(
    SampleType = case_when(
      Category %in% c("Strong_Indoor", "Indoor")  ~ "Indoor",
      Category %in% c("Strong_Outdoor", "Outdoor") ~ "Outdoor",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SampleType))
cat("Indoor/Outdoor assigned:", nrow(data), "\n")
print(table(data$SampleType))

# ------------------------
# SMILES -> molecular formula (rcdk)
# ------------------------
cat("\nComputing formulas from SMILES...\n")
unique_smiles <- unique(data$Top1_SMILES)
formula_lookup <- character(length(unique_smiles))
names(formula_lookup) <- unique_smiles

for (i in seq_along(unique_smiles)) {
  smi <- unique_smiles[i]
  formula_lookup[i] <- tryCatch({
    mol <- parse.smiles(smi)[[1]]
    if (is.null(mol)) return(NA_character_)
    convert.implicit.to.explicit(mol)
    get.mol2formula(mol, charge = 0)@string
  }, error = function(e) NA_character_)
  if (i %% 2000 == 0) cat("  ", i, "/", length(unique_smiles), "\n")
}

data$ConsensusFormula <- formula_lookup[data$Top1_SMILES]
n_failed <- sum(is.na(data$ConsensusFormula))
cat("Formula parse failures:", n_failed, "/", nrow(data), "\n")
data <- data %>% filter(!is.na(ConsensusFormula))

# ------------------------
# Remove implausible elements (As, Ge, Te only)
# ------------------------
implausible_pattern <- "(As|Ge|Te)([0-9]|[A-Z]|$)"
n_before <- nrow(data)
data <- data %>% filter(!str_detect(ConsensusFormula, implausible_pattern))
cat("Removed", n_before - nrow(data), "features with As/Ge/Te (kept", nrow(data), ")\n")

# ------------------------
# Parse elemental composition
# ------------------------
parse_element <- function(formula, element){
  val <- str_extract(formula, paste0(element,"[0-9]*"))
  val <- str_replace(val, element, "")
  val[val==""] <- 1
  val <- as.numeric(val)
  val[is.na(val)] <- 0
  return(val)
}

data <- data %>%
  mutate(
    C = parse_element(ConsensusFormula,"C"),
    H = parse_element(ConsensusFormula,"H"),
    O = parse_element(ConsensusFormula,"O"),
    N = parse_element(ConsensusFormula,"N"),
    S = parse_element(ConsensusFormula,"S"),
    Cl = parse_element(ConsensusFormula,"Cl"),
    Br = parse_element(ConsensusFormula,"Br"),
    F = parse_element(ConsensusFormula,"F")
  )

# ------------------------
# FormulaGroup (exact, kept for the record) + ElementClass (coarse, for plotting)
# Halogen resolved individually: Cl / F / Br / MultiHal (>1 type present)
# ------------------------
data <- data %>%
  rowwise() %>%
  mutate(
    FormulaGroup = paste0(
      if(N>0) ifelse(N==1,"N",paste0("N",N)) else "",
      if(O>0) ifelse(O==1,"O",paste0("O",O)) else "",
      if(S>0) ifelse(S==1,"S",paste0("S",S)) else "",
      if(Cl>0) ifelse(Cl==1,"Cl",paste0("Cl",Cl)) else "",
      if(Br>0) ifelse(Br==1,"Br",paste0("Br",Br)) else "",
      if(F>0) ifelse(F==1,"F",paste0("F",F)) else ""
    ),
    HalogenType = {
      n_types <- sum(Cl>0, Br>0, F>0)
      if (n_types == 0) "None"
      else if (n_types > 1) "MultiHal"
      else if (Cl>0) "Cl"
      else if (Br>0) "Br"
      else "F"
    },
    ElementClass = paste0(
      "CH",
      if(N>0) "N" else "",
      if(O>0) "O" else "",
      if(S>0) "S" else "",
      if(HalogenType != "None") HalogenType else ""
    )
  ) %>%
  ungroup() %>%
  filter(FormulaGroup != "")   # drops pure-CH (no heteroatoms) as before

group_counts <- data %>% count(FormulaGroup)
valid_groups <- group_counts %>% filter(n >= 10) %>% pull(FormulaGroup)
data <- data %>% filter(FormulaGroup %in% valid_groups)

# ------------------------
# DBE, O/C, H/C, KMD, Halogenated (binary, legacy), AI_mod
# ------------------------
data <- data %>%
  mutate(
    DBE = C - H/2 + N/2 + 1,
    O_C = ifelse(C>0, O/C, NA),
    H_C = ifelse(C>0, H/C, NA),
    KendrickMass = MZ * (14.00000 / 14.01565),
    NominalKM = round(KendrickMass),
    KMD = NominalKM - KendrickMass,
    Halogenated = ifelse(HalogenType != "None", "Halogenated", "Non-halogenated"),
    denom = C - 0.5*O - S - N,
    AI_mod = ifelse(denom <= 0, NA, (1 + C - 0.5*O - S - 0.5*(H + N))/denom)
  )

colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")

# ================================================================
# PART A — ORIGINAL PLOT SET
# ================================================================

# ---- Plot 1: ElementClass Distribution ----
class_counts <- data %>%
  group_by(ElementClass, SampleType) %>%
  summarise(Count = n(), .groups="drop") %>%
  tidyr::complete(ElementClass, SampleType, fill=list(Count=0))

class_order <- class_counts %>%
  pivot_wider(names_from=SampleType, values_from=Count, values_fill=0) %>%
  mutate(FoldChange = Outdoor - Indoor) %>%
  arrange(desc(FoldChange)) %>%
  pull(ElementClass)
class_counts$ElementClass <- factor(class_counts$ElementClass, levels=class_order)

p1 <- ggplot(class_counts, aes(x=ElementClass, y=Count, fill=SampleType)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
  scale_fill_manual(values=colors) +
  scale_y_continuous(trans = scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10, 100, 1000, 5000),
                     expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(text=element_text(size=14), axis.text.x=element_text(angle=45, hjust=1), panel.grid=element_blank()) +
  labs(title="Element Class Distribution (NAP-annotated, log scale)", x="Element Class", y="Count (log scale)")
print(p1)
ggsave("ElementClass_distribution.png", p1, width=10, height=6, dpi=600)

# ---- Plot 3: Van Krevelen (overlay, outliers capped) ----
vk_bound_O <- 1.2; vk_bound_H <- 3
n_vk_before <- nrow(data)
data_vk <- data %>% filter(O_C <= vk_bound_O, H_C <= vk_bound_H)
cat("\nVan Krevelen: removed", n_vk_before - nrow(data_vk), "points outside O/C<=1.2, H/C<=3 (",
    round(100*(n_vk_before - nrow(data_vk))/n_vk_before,2), "% )\n")

p3 <- ggplot(data_vk, aes(O_C, H_C, color=SampleType)) +
  geom_point(size=1.2, alpha=0.5) +
  scale_color_manual(values=colors) +
  theme_bw() +
  labs(title="Van Krevelen Diagram (outliers removed: O/C<=1.2, H/C<=3)", x="O/C", y="H/C")
print(p3)
ggsave("Van_Krevelen.png", p3, width=10, height=6, dpi=600)

# ---- Plot 4: DBE Distribution ----
data <- data %>% filter(DBE != -1) %>% mutate(DBE_factor = factor(round(DBE)))
dbe_counts <- data %>% group_by(DBE_factor, SampleType) %>% summarise(Count=n(), .groups="drop")
p4 <- ggplot(dbe_counts, aes(x=DBE_factor, y=Count, fill=SampleType)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
  scale_fill_manual(values=colors) +
  theme_bw() + labs(title="DBE Distribution", x="DBE", y="Count") +
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
print(p4)
ggsave("DBE_distribution.png", p4, width=8, height=6, dpi=600)

# ---- Plot 5: AI_mod Distribution ----
data_ai <- data %>% filter(AI_mod >= -1 & AI_mod <= 1) %>% mutate(AI_mod_factor = factor(round(AI_mod,1)))
ai_counts <- data_ai %>% group_by(AI_mod_factor, SampleType) %>% summarise(Count=n(), .groups="drop")
p5 <- ggplot(ai_counts, aes(x=AI_mod_factor, y=Count, fill=SampleType)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
  scale_fill_manual(values=colors) +
  theme_bw() + labs(title="AI_mod Distribution", x="AI_mod", y="Count") +
  theme(panel.grid=element_blank(), axis.text.x=element_text(angle=45, hjust=1))
print(p5)
ggsave("AI_mod_distribution.png", p5, width=8, height=6, dpi=300)

# ---- Plot 6: KMD (overlay) ----
p6 <- ggplot(data, aes(NominalKM, KMD, color=SampleType)) +
  geom_point(alpha=0.5, size=0.8) +
  scale_color_manual(values=colors) +
  theme_bw() +
  labs(title="Kendrick Mass Defect Plot", x="Nominal Kendrick Mass", y="KMD")
print(p6)
ggsave("KMD_plot.png", p6, width=7, height=6, dpi=600)

# ---- Plot 7: Halogen Enrichment (split by type: F/Cl/Br/MultiHal/None) ----
halo_counts <- data %>%
  count(HalogenType, SampleType) %>%
  mutate(HalogenType = factor(HalogenType, levels=c("None","F","Cl","Br","MultiHal")))

p7 <- ggplot(halo_counts, aes(HalogenType, n, fill=SampleType)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=colors) +
  theme_bw() +
  labs(title="Halogen Type Enrichment", x="Halogen Type", y="Count")
print(p7)
ggsave("Halogen_enrichment.png", p7, width=8, height=6, dpi=300)

# ================================================================
# STATISTICS — DBE, O/C, H/C, AI_mod (Wilcoxon), with BH-adjusted p
# ================================================================
wilcox_dbe    <- wilcox.test(DBE ~ SampleType, data=data)
wilcox_oc     <- wilcox.test(O_C ~ SampleType, data=data)
wilcox_hc     <- wilcox.test(H_C ~ SampleType, data=data)
wilcox_aimod  <- wilcox.test(AI_mod ~ SampleType, data=data)

cat("\nWilcoxon tests:\n")
cat("DBE:\n");    print(wilcox_dbe)
cat("O/C:\n");    print(wilcox_oc)
cat("H/C:\n");    print(wilcox_hc)
cat("AI_mod:\n"); print(wilcox_aimod)

raw_p <- c(DBE=wilcox_dbe$p.value, O_C=wilcox_oc$p.value,
           H_C=wilcox_hc$p.value, AI_mod=wilcox_aimod$p.value)
adj_p <- p.adjust(raw_p, method="BH")
cat("\nBenjamini-Hochberg adjusted p-values:\n"); print(adj_p)

desc_stats <- data %>%
  group_by(SampleType) %>%
  summarise(
    n = n(),
    DBE_median = median(DBE, na.rm=TRUE), DBE_IQR = IQR(DBE, na.rm=TRUE),
    O_C_median = median(O_C, na.rm=TRUE), O_C_IQR = IQR(O_C, na.rm=TRUE),
    H_C_median = median(H_C, na.rm=TRUE), H_C_IQR = IQR(H_C, na.rm=TRUE),
    AI_mod_median = median(AI_mod, na.rm=TRUE), AI_mod_IQR = IQR(AI_mod, na.rm=TRUE),
    .groups="drop"
  )
print(desc_stats)

stats_summary <- desc_stats %>%
  mutate(
    DBE_wilcox_W = wilcox_dbe$statistic,    DBE_wilcox_p = raw_p["DBE"],    DBE_wilcox_p_BH = adj_p["DBE"],
    O_C_wilcox_W = wilcox_oc$statistic,     O_C_wilcox_p = raw_p["O_C"],    O_C_wilcox_p_BH = adj_p["O_C"],
    H_C_wilcox_W = wilcox_hc$statistic,     H_C_wilcox_p = raw_p["H_C"],    H_C_wilcox_p_BH = adj_p["H_C"],
    AI_mod_wilcox_W = wilcox_aimod$statistic, AI_mod_wilcox_p = raw_p["AI_mod"], AI_mod_wilcox_p_BH = adj_p["AI_mod"]
  )

if("Compound_Class" %in% colnames(data)){
  cc_table <- table(data$Compound_Class, data$SampleType)
  chisq_res <- chisq.test(cc_table)
  cat("\nChi-square test (Compound_Class x SampleType):\n")
  print(chisq_res)
  write_csv(as.data.frame.matrix(cc_table) %>% tibble::rownames_to_column("Compound_Class"),
            "Stats_CompoundClass_contingency.csv")
  cat("Chi-sq X2=", round(chisq_res$statistic,2), " df=", chisq_res$parameter,
      " p=", format(chisq_res$p.value, scientific=TRUE), "\n")
} else {
  cat("Compound_Class column not found, skipping chi-square.\n")
}

write_csv(stats_summary, "Stats_summary.csv")
cat("\nSaved: Stats_summary.csv (medians/IQRs + Wilcoxon raw & BH-adjusted p per variable)\n")

write_csv(data, "Chemical_characterization_results_filtered.csv")

# ================================================================
# PART B — FACETED (Van Krevelen, KMD)
# ================================================================
p3_facet <- ggplot(data_vk, aes(O_C, H_C, color=SampleType)) +
  geom_point(size=1.2, alpha=0.5) +
  scale_color_manual(values=colors) +
  facet_wrap(~ SampleType) +
  theme_bw() +
  theme(legend.position="none", strip.background=element_rect(fill="grey95"), strip.text=element_text(face="bold")) +
  labs(title="Van Krevelen Diagram — Indoor vs Outdoor (faceted, outliers removed)", x="O/C", y="H/C")
print(p3_facet)
ggsave("Van_Krevelen_faceted.png", p3_facet, width=12, height=6, dpi=600)

p6_facet <- ggplot(data, aes(NominalKM, KMD, color=SampleType)) +
  geom_point(alpha=0.5, size=0.8) +
  scale_color_manual(values=colors) +
  facet_wrap(~ SampleType) +
  geom_hline(yintercept = 0, color="grey60", linewidth=0.4, linetype="dashed") +
  theme_bw() +
  theme(legend.position="none", strip.background=element_rect(fill="grey95"), strip.text=element_text(face="bold")) +
  labs(title="Kendrick Mass Defect — Indoor vs Outdoor (faceted)", x="Nominal Kendrick Mass", y="KMD")
print(p6_facet)
ggsave("KMD_plot_faceted.png", p6_facet, width=10, height=6, dpi=600)

cat("\nDone.\n")
