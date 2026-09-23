# ============================================================================
# INQUIRE — LC-HRMS: chemical correlation analysis
#
# Full pipeline for compound-compound relationships among confirmed
# target/annotated LC-HRMS compounds:
#   1) Pairwise Spearman correlations (concentration magnitude) + top-40
#      scatter plots
#   2) Correlation network built from (1) — nodes = compounds, edges =
#      FDR-significant pairs at |rho| >= threshold
#   3) Pairwise Fisher's exact co-occurrence tests (presence/absence) —
#      complements (1): valid at very low n, where Spearman breaks down
#   4) Rare-compound co-occurrence network built from (3)
#   5) Sample-level detection heatmap for the rare compounds in (4) —
#      shows WHICH samples drive each co-occurrence link
#
# Each section reads the previous section's own CSV output rather than
# recomputing, so results always match whatever was last produced — run
# top to bottom, in order.
# ============================================================================

library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(igraph)
library(ggraph)
library(pheatmap)

# ---- Set your working/data directory here ----
setwd("path/to/data")

lc_file <- "LC_Level1_2_polished_FINAL.xlsx"

meta_cols <- c("Feature ID","MZ","RT","Name","Molecular formula","InChiKey","SMILES",
               "CAS","Ion type","Confidence level","MS/MS","LOD","LOQ","Quantification")

# ---- Shared: load LC target/annotated concentrations once ----
lc_raw <- read_excel(lc_file)
names(lc_raw) <- trimws(names(lc_raw))
sample_cols <- setdiff(names(lc_raw), meta_cols)
lc_raw[sample_cols] <- lapply(lc_raw[sample_cols], as.numeric)

conc_mat <- as.matrix(lc_raw[, sample_cols])
rownames(conc_mat) <- lc_raw$Name
n_compounds <- nrow(conc_mat)
n_samples   <- ncol(conc_mat)
cat("LC: ", n_compounds, "compounds x", n_samples, "samples\n")

# ============================================================================
# SECTION 1: Pairwise Spearman correlations
# Tests all pairwise combinations of target/annotated compounds, ranks by
# |rho| among FDR-significant pairs, and plots the top 40 as scatter plots.
#
# Method notes (read before trusting the output):
#  - Spearman rank correlation: robust to skewed concentration data and
#    non-detects, no log-transform/normality assumption required.
#  - Samples pooled across Indoor + Outdoor (all paired samples together) —
#    this tests general co-occurrence/shared-source patterns, NOT
#    deployment-specific correlation.
#  - MIN_DF: both compounds in a pair must be detected (>0) in at least this
#    fraction of samples, or the pair is skipped. At MIN_DF = 0.10 (~41 of
#    410 samples), pairs backed by very few detections can still produce an
#    unstable rho — check the `n` column before trusting any single pair,
#    especially further down the top-40 list.
#  - BH-FDR correction applied across ALL pairwise tests. Only pairs with
#    padj < ALPHA are eligible for the "top 40" ranking.
# ============================================================================

MIN_DF   <- 0.10   # minimum detection frequency required for BOTH compounds in a pair
ALPHA    <- 0.05   # BH-FDR threshold for a pair to be eligible for top-N ranking
TOP_N    <- 40

det_freq <- rowSums(conc_mat > 0, na.rm = TRUE) / n_samples
keep <- names(det_freq)[det_freq >= MIN_DF]
cat("Compounds passing DF >=", MIN_DF, ":", length(keep), "of", n_compounds, "\n")

conc_mat_f <- conc_mat[keep, , drop = FALSE]

pairs <- combn(rownames(conc_mat_f), 2, simplify = FALSE)
cat("Testing", length(pairs), "compound pairs...\n")

results <- vector("list", length(pairs))
for (i in seq_along(pairs)) {
  a <- pairs[[i]][1]; b <- pairs[[i]][2]
  xa <- conc_mat_f[a, ]; xb <- conc_mat_f[b, ]
  # Only correlate samples where BOTH compounds were actually detected
  # (>0), not just non-missing. With non-detects encoded as 0 rather than
  # NA, using !is.na() alone includes every sample regardless of detection
  # — for sparse compounds this lets shared absence (both = 0) drive most
  # of the correlation, inflating rho without reflecting real co-occurrence.
  ok <- !is.na(xa) & !is.na(xb) & xa > 0 & xb > 0
  n_ok <- sum(ok)
  if (n_ok < 10) {
    results[[i]] <- data.frame(Compound_A = a, Compound_B = b, rho = NA, p = NA, n = n_ok)
    next
  }
  ct <- suppressWarnings(cor.test(xa[ok], xb[ok], method = "spearman"))
  results[[i]] <- data.frame(Compound_A = a, Compound_B = b,
                              rho = unname(ct$estimate), p = ct$p.value, n = n_ok)
}
cor_results <- bind_rows(results) %>% filter(!is.na(rho))
cor_results$padj <- p.adjust(cor_results$p, method = "BH")

write.csv(cor_results, "LC_all_pairwise_correlations.csv", row.names = FALSE)
cat("Saved: LC_all_pairwise_correlations.csv (", nrow(cor_results), "pairs tested)\n")

top_pairs <- cor_results %>%
  filter(padj < ALPHA) %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = TOP_N)

cat("\nTop", TOP_N, "strongest significant correlations (LC):\n")
print(top_pairs)
write.csv(top_pairs, "LC_top40_correlations.csv", row.names = FALSE)

if (nrow(top_pairs) == 0) {
  warning("No pairs reached padj < ", ALPHA, " with n >= 10 and DF >= ", MIN_DF,
          ". Loosen MIN_DF/ALPHA and re-run.")
}

# ---- Build one scatter plot per pair, print + save each individually.
# Loops over however many pairs actually passed the filter (handles < 40
# pairs automatically, unlike a hardcoded p1..p40 block). ----
plot_list <- list()
for (i in seq_len(nrow(top_pairs))) {
  a <- top_pairs$Compound_A[i]; b <- top_pairs$Compound_B[i]
  # Same detection filter as the correlation test itself — only points
  # where both compounds were actually detected, not all samples.
  df <- data.frame(x = conc_mat_f[a, ], y = conc_mat_f[b, ])
  df <- df[!is.na(df$x) & !is.na(df$y) & df$x > 0 & df$y > 0, ]

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.5, color = "#7B3F9E") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    labs(title = paste0(a, " vs. ", b),
         subtitle = sprintf("Spearman rho = %.2f | padj = %.2e | n = %d",
                             top_pairs$rho[i], top_pairs$padj[i], top_pairs$n[i]),
         x = paste0(a, " (log10 conc.)"), y = paste0(b, " (log10 conc.)")) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 8))

  plot_list[[i]] <- p
}
cat("\nBuilt", length(plot_list), "scatter plots.\n")

for (i in seq_along(plot_list)) {
  a <- top_pairs$Compound_A[i]; b <- top_pairs$Compound_B[i]
  fname <- sprintf("LC_corr_%02d_%s_vs_%s.png", i,
                    gsub("[^A-Za-z0-9]", "", a), gsub("[^A-Za-z0-9]", "", b))
  print(plot_list[[i]])
  ggsave(fname, plot_list[[i]], width = 5, height = 4, dpi = 600)
}
cat("Saved", length(plot_list), "individual scatter plots (LC_corr_NN_*.png)\n\n")

# ============================================================================
# SECTION 2: Compound co-occurrence network (from Section 1's correlations)
# Nodes = compounds, edges = FDR-significant pairs at |rho| >= RHO_THRESH.
#
# Method notes:
#  - RHO_THRESH = 0.5: a pre-specified effect-size cutoff (moderate-to-
#    strong correlation), not a round number chosen to fit a figure —
#    this is the defensible way to decide "how many pairs to show" instead
#    of an arbitrary top-N.
#  - Only pairs with padj < ALPHA are eligible as edges.
#  - Nodes are colored by connected component (cluster).
#  - The largest connected component can end up dense/illegible at this
#    threshold — GIANT_CLUSTER_THRESH applies a stricter threshold ONLY
#    within the largest cluster, to help it resolve into sub-structure,
#    without changing the rest of the network.
# ============================================================================

RHO_THRESH           <- 0.5    # minimum |rho| for an edge
GIANT_CLUSTER_THRESH <- 0.65   # NULL = off

cor_results_net <- read_csv("LC_all_pairwise_correlations.csv", show_col_types = FALSE)

edges <- cor_results_net %>%
  filter(padj < ALPHA, abs(rho) >= RHO_THRESH)

cat("Edges at |rho| >=", RHO_THRESH, ", padj <", ALPHA, ":", nrow(edges), "\n")

g <- graph_from_data_frame(
  edges %>% select(Compound_A, Compound_B, rho, padj, n),
  directed = FALSE
)

cat("Nodes (unique compounds):", vcount(g), "\n")

comp <- components(g)
cat("Connected clusters:", comp$no, "\n")
cluster_sizes <- sort(table(comp$membership), decreasing = TRUE)
print(cluster_sizes)

V(g)$cluster <- as.factor(comp$membership)

if (!is.null(GIANT_CLUSTER_THRESH)) {
  giant_id <- as.integer(names(cluster_sizes)[1])
  giant_nodes <- names(comp$membership)[comp$membership == giant_id]

  cat("\nRe-thresholding giant cluster (", length(giant_nodes),
      "nodes) at |rho| >=", GIANT_CLUSTER_THRESH, "\n")

  edges_final <- edges %>%
    filter(!(Compound_A %in% giant_nodes & Compound_B %in% giant_nodes) |
             abs(rho) >= GIANT_CLUSTER_THRESH)

  g <- graph_from_data_frame(
    edges_final %>% select(Compound_A, Compound_B, rho, padj, n),
    directed = FALSE
  )
  comp <- components(g)
  V(g)$cluster <- as.factor(comp$membership)
  cat("After re-thresholding: ", vcount(g), "nodes,", ecount(g), "edges,",
      comp$no, "clusters\n")
}

E(g)$sign <- ifelse(E(g)$rho > 0, "Positive", "Negative")

set.seed(42)
p_network <- ggraph(g, layout = "fr") +
  geom_edge_link(aes(color = sign, width = abs(rho)), alpha = 0.4) +
  scale_edge_color_manual(values = c("Positive" = "#B03A2E", "Negative" = "#1F618D"),
                          name = "Correlation") +
  scale_edge_width(range = c(0.3, 2), name = "|rho|") +
  geom_node_point(aes(color = cluster), size = 8, show.legend = FALSE) +
  geom_node_text(aes(label = name), size = 3.2, repel = TRUE, max.overlaps = 30) +
  labs(title = "LC-HRMS compound co-occurrence network",
       subtitle = sprintf("|Spearman rho| >= %.2f, FDR-significant, jointly-detected samples only | %d edges, %d compounds, %d clusters",
                          RHO_THRESH, ecount(g), vcount(g), comp$no)) +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))

print(p_network)
ggsave("LC_correlation_network.png", p_network, width = 16, height = 12, dpi = 600)
cat("\nSaved: LC_correlation_network.png\n")

cluster_table <- data.frame(Compound = names(comp$membership),
                             Cluster = comp$membership) %>%
  arrange(Cluster, Compound)
write.csv(cluster_table, "LC_correlation_network_clusters.csv", row.names = FALSE)
cat("Saved: LC_correlation_network_clusters.csv\n\n")

# ============================================================================
# SECTION 3: Pairwise Fisher's exact co-occurrence tests
# Tests presence/absence co-occurrence (NOT concentration correlation —
# that's Section 1). Fisher's exact stays valid at very low n, which is
# exactly where Spearman breaks down (e.g. two rare pesticides detected
# together in only a handful of samples can still give a solid, testable
# co-occurrence signal here).
#
# Method notes (read before trusting the output):
#  - One-sided Fisher's exact test (alternative = "greater") on each
#    compound pair's 2x2 detection contingency table (both / A-only /
#    B-only / neither), across ALL samples — this is a presence/absence
#    question, not a magnitude one, so unlike Section 1 there's no
#    "joint-detection-only" filtering here; absence is exactly what's
#    being tested, not excluded.
#  - MIN_DETECTIONS: a compound must be detected (>0) in at least this many
#    samples (absolute count, not %) to be included at all.
#  - Odds ratio can be Infinite (perfect separation). Reported as Inf, not
#    treated as an error.
#  - BH-FDR correction applied across ALL pairwise tests.
#  - MAX_DF_FOR_RANKING: for the ranked "top 40" table specifically (not
#    the full CSV), both compounds in a pair must ALSO have a DF at or
#    below this ceiling. Raw p-value ranking is otherwise dominated by
#    common compound pairs — with large n, even a modest odds ratio
#    produces a tiny p-value, burying rare-but-tightly-linked pairs far
#    down the list. 15% sits just above the lowest-frequency compound in
#    a known true-positive case (a pesticide transformation-product pair,
#    ~2-10% DF each). Does NOT affect the full CSV.
# ============================================================================

MIN_DETECTIONS     <- 3     # minimum absolute number of detections (>0) required per compound
MAX_DF_FOR_RANKING <- 0.15  # DF ceiling (both compounds) for the ranked top-40 table only

det_mat <- (conc_mat > 0) & !is.na(conc_mat)
n_detections <- rowSums(det_mat)

keep_fisher <- names(n_detections)[n_detections >= MIN_DETECTIONS]
cat("Compounds passing >=", MIN_DETECTIONS, "detections:", length(keep_fisher), "of", n_compounds, "\n")

det_mat_f <- det_mat[keep_fisher, , drop = FALSE]

pairs_fisher <- combn(rownames(det_mat_f), 2, simplify = FALSE)
cat("Testing", length(pairs_fisher), "compound pairs...\n")

results_fisher <- vector("list", length(pairs_fisher))
for (i in seq_along(pairs_fisher)) {
  a <- pairs_fisher[[i]][1]; b <- pairs_fisher[[i]][2]
  xa <- det_mat_f[a, ]; xb <- det_mat_f[b, ]

  both    <- sum(xa & xb)
  a_only  <- sum(xa & !xb)
  b_only  <- sum(!xa & xb)
  neither <- sum(!xa & !xb)

  tbl <- matrix(c(both, a_only, b_only, neither), nrow = 2)
  ft <- fisher.test(tbl, alternative = "greater")

  results_fisher[[i]] <- data.frame(
    Compound_A = a, Compound_B = b,
    n_both = both, n_A_only = a_only, n_B_only = b_only, n_neither = neither,
    odds_ratio = unname(ft$estimate), p = ft$p.value
  )
}
fisher_results <- bind_rows(results_fisher)
fisher_results$padj <- p.adjust(fisher_results$p, method = "BH")

fisher_results <- fisher_results %>%
  mutate(
    total   = n_both + n_A_only + n_B_only + n_neither,
    DF_A    = (n_both + n_A_only) / total,
    DF_B    = (n_both + n_B_only) / total
  ) %>%
  arrange(p)

write.csv(fisher_results, "LC_all_pairwise_fisher_cooccurrence.csv", row.names = FALSE)
cat("Saved: LC_all_pairwise_fisher_cooccurrence.csv (", nrow(fisher_results), "pairs tested)\n")

n_sig <- sum(fisher_results$padj < ALPHA, na.rm = TRUE)
n_inf <- sum(is.infinite(fisher_results$odds_ratio))
cat("\nPairs significant at padj <", ALPHA, ":", n_sig, "\n")
cat("Pairs with infinite OR (perfect co-occurrence separation):", n_inf, "\n")

top_pairs_fisher <- fisher_results %>%
  filter(padj < ALPHA, DF_A <= MAX_DF_FOR_RANKING, DF_B <= MAX_DF_FOR_RANKING) %>%
  arrange(p) %>%
  slice_head(n = TOP_N)

cat("\nTop", TOP_N, "rare-pair co-occurrences (DF <=", MAX_DF_FOR_RANKING, "each, padj <", ALPHA, "):\n")
print(top_pairs_fisher %>% select(Compound_A, Compound_B, n_both, DF_A, DF_B, odds_ratio, p, padj))

write.csv(top_pairs_fisher, "LC_top40_fisher_cooccurrence.csv", row.names = FALSE)
cat("Saved: LC_top40_fisher_cooccurrence.csv\n\n")

# ============================================================================
# SECTION 4: Rare-compound co-occurrence network (from Section 3's Fisher
# results). Companion to Section 2's Spearman network — this one
# visualizes presence/absence co-occurrence among RARE compound pairs.
#
# Method notes:
#  - Edges = pairs with padj < ALPHA AND both compounds' DF <= MAX_DF, same
#    rarity ceiling used for the top-40 table — NOT the full unfiltered
#    Fisher CSV. Without this filter, the network would be dominated by
#    near-ubiquitous compounds that trivially "co-occur" because both are
#    common — that's not the rare-pair structure this network is meant to
#    show.
#  - Fisher's exact here is one-sided (positive co-occurrence only), so
#    unlike the Spearman network there is no negative-correlation edge
#    color — all edges represent the same kind of relationship.
#  - Edge width = log10(odds ratio). Infinite odds ratios (perfect
#    separation) are capped at CAPPED_OR (1.5x the max finite OR in the
#    filtered set) purely for plotting; the underlying CSV/cluster table
#    still records them as Inf. A capped edge is visually flagged (dashed).
#  - Nodes colored by connected component (cluster), for visual/
#    interpretive consistency with Section 2's network.
# ============================================================================

MAX_DF <- 0.15   # same rarity ceiling as the top-40 table
GIANT_CLUSTER_THRESH_FISHER <- NULL  # e.g. a higher OR floor, to re-split a dense cluster if needed

fisher_results_net <- read_csv("LC_all_pairwise_fisher_cooccurrence.csv", show_col_types = FALSE)

edges_fisher <- fisher_results_net %>%
  filter(padj < ALPHA, DF_A <= MAX_DF, DF_B <= MAX_DF)

cat("Edges at padj <", ALPHA, ", DF <=", MAX_DF, "(both compounds):", nrow(edges_fisher), "\n")

if (nrow(edges_fisher) == 0) {
  stop("No pairs pass the current ALPHA/MAX_DF filter — loosen MAX_DF and re-run.")
}

finite_or <- edges_fisher$odds_ratio[is.finite(edges_fisher$odds_ratio)]
CAPPED_OR <- 1.5 * max(finite_or, na.rm = TRUE)
cat("Finite OR range:", min(finite_or), "-", max(finite_or), "| Capping Inf edges at", round(CAPPED_OR), "\n")

edges_fisher <- edges_fisher %>%
  mutate(
    OR_was_infinite = is.infinite(odds_ratio),
    OR_for_plot = ifelse(is.infinite(odds_ratio), CAPPED_OR, odds_ratio)
  )

g_fisher <- graph_from_data_frame(
  edges_fisher %>% select(Compound_A, Compound_B, odds_ratio, OR_for_plot, OR_was_infinite, padj, n_both),
  directed = FALSE
)

cat("Nodes (unique compounds):", vcount(g_fisher), "\n")

comp_fisher <- components(g_fisher)
cat("Connected clusters:", comp_fisher$no, "\n")
cluster_sizes_fisher <- sort(table(comp_fisher$membership), decreasing = TRUE)
print(cluster_sizes_fisher)

V(g_fisher)$cluster <- as.factor(comp_fisher$membership)

if (!is.null(GIANT_CLUSTER_THRESH_FISHER)) {
  giant_id <- as.integer(names(cluster_sizes_fisher)[1])
  giant_nodes <- names(comp_fisher$membership)[comp_fisher$membership == giant_id]

  cat("\nRe-thresholding giant cluster (", length(giant_nodes),
      "nodes) at OR >=", GIANT_CLUSTER_THRESH_FISHER, "\n")

  edges_final_fisher <- edges_fisher %>%
    filter(!(Compound_A %in% giant_nodes & Compound_B %in% giant_nodes) |
             OR_for_plot >= GIANT_CLUSTER_THRESH_FISHER)

  g_fisher <- graph_from_data_frame(
    edges_final_fisher %>% select(Compound_A, Compound_B, odds_ratio, OR_for_plot, OR_was_infinite, padj, n_both),
    directed = FALSE
  )
  comp_fisher <- components(g_fisher)
  V(g_fisher)$cluster <- as.factor(comp_fisher$membership)
  cat("After re-thresholding: ", vcount(g_fisher), "nodes,", ecount(g_fisher), "edges,",
      comp_fisher$no, "clusters\n")
}

E(g_fisher)$edge_style <- ifelse(E(g_fisher)$OR_was_infinite, "Infinite OR (capped)", "Finite OR")

set.seed(42)
p_network_fisher <- ggraph(g_fisher, layout = "fr") +
  geom_edge_link(aes(width = log10(OR_for_plot), linetype = edge_style),
                  color = "#B03A2E", alpha = 0.5) +
  scale_edge_width(range = c(0.3, 2.5), name = "log10(OR)") +
  scale_edge_linetype_manual(values = c("Finite OR" = "solid", "Infinite OR (capped)" = "22"),
                              name = NULL) +
  geom_node_point(aes(color = cluster), size = 4, show.legend = FALSE) +
  geom_node_text(aes(label = name), size = 3.2, repel = TRUE, max.overlaps = 30) +
  labs(title = "LC-HRMS rare-compound co-occurrence network (Fisher's exact)",
       subtitle = sprintf("padj < %.2f, both compounds DF <= %.0f%% | %d edges, %d compounds, %d clusters",
                           ALPHA, MAX_DF * 100, ecount(g_fisher), vcount(g_fisher), comp_fisher$no)) +
  theme_void(base_size = 11) +
  theme(plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))

print(p_network_fisher)
ggsave("LC_fisher_network.png", p_network_fisher, width = 16, height = 12, dpi = 300)
cat("\nSaved: LC_fisher_network.png\n")

cluster_table_fisher <- data.frame(Compound = names(comp_fisher$membership),
                             Cluster = comp_fisher$membership) %>%
  arrange(Cluster, Compound)
write.csv(cluster_table_fisher, "LC_fisher_network_clusters.csv", row.names = FALSE)
cat("Saved: LC_fisher_network_clusters.csv\n\n")

# ============================================================================
# SECTION 5: Sample-level detection heatmap for rare, co-occurring compounds
# Complements Section 4: the network shows WHICH pairs are statistically
# linked (OR/padj), but not WHICH SAMPLES drive that link. This shows that
# directly — rows = compounds, columns = samples, colored by detected/
# not-detected, rows clustered by JACCARD distance (the correct distance
# metric for binary presence/absence data — it measures shared-sample
# overlap directly, not just statistical association). If a cluster is
# real and not a statistical artifact, its compounds should show red
# blocks concentrated in the SAME handful of sample columns.
#
# Method notes:
#  - Restricted to the same rare-compound set as Section 4 (compounds with
#    DF <= MAX_DF appearing in at least one significant pair) — plotting
#    every compound x every sample would be unreadable and pointless
#    (most compounds are common, not what this figure is about).
#  - Row order = hierarchical clustering on Jaccard distance
#    (dist(..., method = "binary") in R IS Jaccard distance on 0/1 data).
#    Column order = hierarchical clustering the same way, so samples that
#    share similar rare-compound detection profiles group together too.
#  - Column annotation by country, if country can be parsed from the
#    sample name — lets you see directly whether a cluster is
#    geographically concentrated rather than spread evenly.
#  - This is presence/absence only (matches the Fisher test's data type) —
#    NOT concentration.
# ============================================================================

SHOW_COLUMN_DENDROGRAM <- TRUE  # set FALSE to hide the (dense) column tree
                                  # and show only the country color strip —
                                  # columns stay Jaccard-ordered either way

det_mat_full <- (conc_mat > 0) & !is.na(conc_mat)   # full presence/absence matrix, all compounds

rare_sig <- fisher_results_net %>%
  filter(padj < ALPHA, DF_A <= MAX_DF, DF_B <= MAX_DF)

rare_compounds <- unique(c(rare_sig$Compound_A, rare_sig$Compound_B))
cat("Rare compounds in at least one significant pair:", length(rare_compounds), "\n")

det_mat_heatmap <- det_mat_full[rownames(det_mat_full) %in% rare_compounds, , drop = FALSE]
storage.mode(det_mat_heatmap) <- "numeric"  # pheatmap wants numeric, not logical

# Drop samples where NONE of these rare compounds were detected — otherwise
# most columns are entirely empty and just add visual noise
any_detected <- colSums(det_mat_heatmap) > 0
det_mat_active <- det_mat_heatmap[, any_detected, drop = FALSE]
cat("Samples with at least one rare-compound detection:", ncol(det_mat_active), "of", ncol(det_mat_heatmap), "\n")

# Country annotation (parsed from sample name, e.g. "NL_Indoor_AA" -> "NL")
sample_country <- str_match(colnames(det_mat_active), "^([A-Z]{2})_")[, 2]
col_annotation <- data.frame(Country = sample_country)
rownames(col_annotation) <- colnames(det_mat_active)

country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C"
)
ann_colors <- list(Country = country_colors[names(country_colors) %in% unique(sample_country)])

# Highlight known named clusters as a row annotation (from the discussion
# text / top-40 table) — everything else labeled "Other"
triazine <- c("Prometon","Terbutylazine-2-hydroxy","Terbutylazine-desethyl-2-hydroxy","Terbutryn")
peg      <- c("Hexaethylene glycol","Decaethylene glycol","Undecaethylene glycol")
b6       <- c("Pyridoxine","6-(hydroxymethyl)pyridin-3-ol","3-Pyridinemethanol")
phthalate_sub <- c("Bis(4-methyl-pentyl) phthalate","Dihexyl phthalate","Monoisobutyl phthalate",
                    "Triethylene glycol bis(2-ethylhexanoate)","Phenylacetic acid")
biogenic <- c("Histamine","Glycerophosphocholine","L-Tyrosine","Methyl 4-hydroxycinnamate","N(6)-OH-Me-Adenosine")

group_lookup <- rep("Other", nrow(det_mat_active))
names(group_lookup) <- rownames(det_mat_active)
group_lookup[names(group_lookup) %in% triazine]      <- "Triazine herbicides"
group_lookup[names(group_lookup) %in% peg]            <- "PEG oligomers"
group_lookup[names(group_lookup) %in% b6]              <- "Vitamin B6-related"
group_lookup[names(group_lookup) %in% phthalate_sub]   <- "Phthalate sub-cluster"
group_lookup[names(group_lookup) %in% biogenic]        <- "Biogenic amine/metabolite"

row_annotation <- data.frame(Group = group_lookup)
rownames(row_annotation) <- rownames(det_mat_active)

group_colors <- c(
  "Triazine herbicides"       = "#B03A2E",
  "PEG oligomers"             = "#1F618D",
  "Vitamin B6-related"        = "#F5B041",
  "Phthalate sub-cluster"     = "#7B3F9E",
  "Biogenic amine/metabolite" = "#28A745",
  "Other"                     = "grey85"
)
ann_colors$Group <- group_colors

# Truncate long row labels for display only — the full compound name is
# preserved in the exported CSVs below, only the heatmap's own row labels
# are shortened.
display_names <- rownames(det_mat_active)
is_long <- nchar(display_names) > 40
display_names[is_long] <- paste0(substr(display_names[is_long], 1, 37), "...")
rownames(det_mat_active)  <- display_names
rownames(row_annotation)  <- display_names

p_heatmap <- pheatmap(
  det_mat_active,
  color = c("white", "#B03A2E"),
  legend_breaks = c(0, 1),
  legend_labels = c("Not detected", "Detected"),
  clustering_distance_rows = "binary",   # Jaccard distance on 0/1 data
  clustering_distance_cols = "binary",
  clustering_method = "ward.D2",
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  show_colnames = FALSE,      # too many samples to label individually
  treeheight_col = if (SHOW_COLUMN_DENDROGRAM) 50 else 0,
  fontsize_row = 6,
  main = sprintf("Sample-level detection of rare co-occurring LC-HRMS compounds (DF <= %.0f%%, padj < %.2f in >=1 pair)\nRows/columns clustered by Jaccard distance",
                  MAX_DF * 100, ALPHA),
  filename = "LC_rare_compound_sample_heatmap.png",
  width = 18, height = 11, dpi = 300
)

cat("\nSaved: LC_rare_compound_sample_heatmap.png\n")

write.csv(det_mat_active, "LC_rare_compound_detection_matrix.csv")
write.csv(data.frame(Compound = rownames(row_annotation), Group = row_annotation$Group),
          "LC_rare_compound_groups.csv", row.names = FALSE)
write.csv(data.frame(SampleID = rownames(col_annotation), Country = col_annotation$Country),
          "LC_rare_compound_sample_countries.csv", row.names = FALSE)
cat("Saved: LC_rare_compound_detection_matrix.csv, LC_rare_compound_groups.csv, LC_rare_compound_sample_countries.csv\n\n")

cat("==== DONE — LC chemical correlation analysis complete ====\n")
