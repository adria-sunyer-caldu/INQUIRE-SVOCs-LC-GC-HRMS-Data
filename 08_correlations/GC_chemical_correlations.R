# ============================================================================
# INQUIRE — GC-HRMS: chemical correlation analysis
# Direct GC counterpart of LC_chemical_correlations.R — same methods, same
# thresholds, restricted to confirmed compounds (ID level 1-2), matching
# the confidence-tier convention used elsewhere in the GC pipeline.
#
# Full pipeline for compound-compound relationships among confirmed GC-HRMS
# compounds:
#   1) Pairwise Spearman correlations (concentration magnitude) + top-40
#      scatter plots
#   2) Correlation network built from (1) — Louvain community detection
#      (see Section 2 notes for why, unlike the LC network, which uses
#      connected components)
#   3) Pairwise Fisher's exact co-occurrence tests (presence/absence) —
#      complements (1): valid at very low n, where Spearman breaks down
#   4) Sample-level detection heatmap for the rare compounds from (3) —
#      shows WHICH samples drive each co-occurrence link
#
# GC has no separate rare-compound Fisher NETWORK step (unlike LC) — the
# sample heatmap in Section 4 reads directly from Section 3's output.
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
library(patchwork)
library(pheatmap)

# ---- Set your working/data directory here ----
setwd("path/to/data")

gc_file <- "GC_Level1_2_polished_FINAL.xlsx"

# ---- Shared: load GC data, restrict to confirmed compounds (ID level 1-2),
# build a de-duplicated compound name column, once ----
# NOTE: read with no sheet= argument — the polished file has a single sheet.
gc_raw <- read_excel(gc_file)

gc_raw$`ID level` <- suppressWarnings(as.numeric(gc_raw$`ID level`))
gc_confirmed <- gc_raw %>% filter(`ID level` %in% c(1, 2))
cat("GC: ", nrow(gc_confirmed), "confirmed compounds (ID level 1-2)\n")

sample_cols <- names(gc_confirmed)[str_detect(names(gc_confirmed), "^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$")]
gc_confirmed[sample_cols] <- lapply(gc_confirmed[sample_cols], as.numeric)

# Use IUPAC name if present, else Database match, else Feature_ID, as the row label
compound_name <- ifelse(!is.na(gc_confirmed$`IUPAC name`) & gc_confirmed$`IUPAC name` != "",
                         gc_confirmed$`IUPAC name`,
                         ifelse(!is.na(gc_confirmed$`Database match`) & gc_confirmed$`Database match` != "",
                                gc_confirmed$`Database match`,
                                gc_confirmed$Feature_ID))

# De-duplicate names if needed (append Feature_ID suffix on collision) —
# same precaution as the LC file's duplicate-InChIKey issue; check
# manually if this fires, since it may indicate the same compound
# annotated twice rather than a genuine naming collision.
dup <- duplicated(compound_name) | duplicated(compound_name, fromLast = TRUE)
if (any(dup)) {
  cat("WARNING:", sum(dup), "duplicate compound names found — check these before trusting results:\n")
  print(unique(compound_name[dup]))
}
compound_name[dup] <- paste0(compound_name[dup], " [", gc_confirmed$Feature_ID[dup], "]")

conc_mat <- as.matrix(gc_confirmed[, sample_cols])
rownames(conc_mat) <- compound_name
n_compounds <- nrow(conc_mat)
n_samples   <- ncol(conc_mat)
cat("GC: ", n_compounds, "compounds x", n_samples, "samples\n")

ALPHA <- 0.05   # BH-FDR threshold, used throughout every section below

# ============================================================================
# SECTION 1: Pairwise Spearman correlations
# Direct GC counterpart of Section 1 in LC_chemical_correlations.R — same
# method, same fixes (joint-detection-only filtering), same top-N structure.
#
# Method notes (identical logic to the LC script, read there for detail):
#  - Spearman rank correlation, samples pooled Indoor + Outdoor.
#  - ok <- !is.na(xa) & !is.na(xb) & xa > 0 & xb > 0 — only jointly-detected
#    samples are correlated (NOT just non-missing), same fix applied to the
#    LC script after non-detects encoded as 0 (not NA) were found to
#    inflate rho via shared absence.
#  - MIN_DF: fraction of GC samples both compounds must be detected in.
#    GC's confirmed set is far smaller (~106 compounds) than LC's (319) —
#    expect far fewer pairs tested and possibly fewer/no FDR-significant
#    results. Check console output before assuming this dataset behaves
#    like LC's.
# ============================================================================

MIN_DF <- 0.10   # minimum detection frequency required for BOTH compounds in a pair
TOP_N  <- 40

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

write.csv(cor_results, "GC_all_pairwise_correlations.csv", row.names = FALSE)
cat("Saved: GC_all_pairwise_correlations.csv (", nrow(cor_results), "pairs tested)\n")

top_pairs <- cor_results %>%
  filter(padj < ALPHA) %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = TOP_N)

cat("\nTop", TOP_N, "strongest significant correlations (GC):\n")
print(top_pairs)
write.csv(top_pairs, "GC_top40_correlations.csv", row.names = FALSE)

if (nrow(top_pairs) == 0) {
  warning("No pairs reached padj < ", ALPHA, " with n >= 10 and DF >= ", MIN_DF,
          ". Loosen MIN_DF/ALPHA and re-run. GC's much smaller confirmed set ",
          "(vs LC's) may simply not have enough power for this — check",
          " before assuming a script problem.")
}

# ---- Build one scatter plot per pair, print + save each individually.
# Loops over however many pairs actually passed the filter (handles < 40
# pairs automatically, unlike a hardcoded p1..p40 block — very likely for
# GC given its much smaller confirmed compound set). ----
plot_list <- list()
for (i in seq_len(nrow(top_pairs))) {
  a <- top_pairs$Compound_A[i]; b <- top_pairs$Compound_B[i]
  df <- data.frame(x = conc_mat_f[a, ], y = conc_mat_f[b, ])
  df <- df[!is.na(df$x) & !is.na(df$y) & df$x > 0 & df$y > 0, ]

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(alpha = 0.5, color = "#B5651D") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.5, linetype = "dashed") +
    scale_x_log10() + scale_y_log10() +
    labs(title = paste0(a, " vs. ", b),
         subtitle = sprintf("Spearman rho = %.2f | padj = %.2e | n = %d",
                             top_pairs$rho[i], top_pairs$padj[i], top_pairs$n[i]),
         x = paste0(a, " (log10 area)"), y = paste0(b, " (log10 area)")) +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(size = 9, face = "bold"),
          plot.subtitle = element_text(size = 8))

  plot_list[[i]] <- p
}
cat("\nBuilt", length(plot_list), "scatter plots.\n")

for (i in seq_along(plot_list)) {
  a <- top_pairs$Compound_A[i]; b <- top_pairs$Compound_B[i]
  fname <- sprintf("GC_corr_%02d_%s_vs_%s.png", i,
                    gsub("[^A-Za-z0-9]", "", a), gsub("[^A-Za-z0-9]", "", b))
  print(plot_list[[i]])
  ggsave(fname, plot_list[[i]], width = 5, height = 4, dpi = 600)
}
cat("Saved", length(plot_list), "individual scatter plots (GC_corr_NN_*.png)\n\n")

# ============================================================================
# SECTION 2: Compound co-occurrence network (Louvain community detection)
# Built from Section 1's correlations.
#
# Method notes:
#  - Unlike the LC network, raw |rho| re-thresholding (the
#    GIANT_CLUSTER_THRESH connected-components approach used for LC) does
#    NOT resolve the GC network into distinct clusters — checked up to
#    |rho| >= 0.85 and one ~35-node component still persists. So this
#    section uses Louvain community detection (modularity-based, weighted
#    by |rho|) instead of connected components to assign cluster color —
#    finds substructure WITHOUT deleting real edges.
#  - RHO_THRESH = 0.5 kept identical to the LC/GC pairwise correlation
#    sections for consistency.
# ============================================================================

RHO_THRESH <- 0.5    # minimum |rho| for an edge (used for community detection)
FACET_BY_CLUSTER <- TRUE  # TRUE = one panel per Louvain community (recommended
                           # given only 3 communities — each gets its own
                           # layout, so within-cluster structure is legible
                           # without needing an extra edge-strength cutoff).
                           # FALSE = old single-panel behavior using
                           # PLOT_RHO_THRESH below.
PLOT_RHO_THRESH <- 0.85   # only used if FACET_BY_CLUSTER = FALSE

cor_results_net <- read_csv("GC_all_pairwise_correlations.csv", show_col_types = FALSE)

edges <- cor_results_net %>%
  filter(padj < ALPHA, abs(rho) >= RHO_THRESH)

cat("Edges at |rho| >=", RHO_THRESH, ", padj <", ALPHA, ":", nrow(edges), "\n")

g <- graph_from_data_frame(
  edges %>% select(Compound_A, Compound_B, rho, padj, n),
  directed = FALSE
)

cat("Nodes (unique compounds):", vcount(g), "\n")
cat("Raw connected components:", components(g)$no, "\n")  # expect 1 — this is why Louvain is used

E(g)$weight <- abs(E(g)$rho)
set.seed(42)
lou <- cluster_louvain(g, weights = E(g)$weight)

V(g)$cluster <- as.factor(membership(lou))
cat("Louvain communities:", length(lou), "\n")
print(sort(table(membership(lou)), decreasing = TRUE))
cat("Modularity:", modularity(lou), "\n")

if (FACET_BY_CLUSTER) {
  mem <- membership(lou)                    # named integer vector: compound -> cluster id
  cluster_ids <- sort(unique(as.integer(mem)))
  panels <- list()

  # Distinct fixed color per cluster panel (not mapped via aes, since each
  # panel only contains one cluster — a manual palette keeps panels visually
  # distinguishable from each other at a glance)
  cluster_palette <- c("#7B3F9E", "#D9782D", "#2E8B57", "#25998F", "#F36E98", "#9862A2")
  cluster_colors <- setNames(cluster_palette[seq_along(cluster_ids)], as.character(cluster_ids))

  for (cid in cluster_ids) {
    # Using membership() (plain integer vector keyed by compound name)
    # instead of comparing V(g)$cluster factors directly — factor == factor
    # comparisons can silently misbehave if level sets/order ever diverge
    # between cluster_ids and V(g)$cluster, giving wrong or empty subsets
    # with no error.
    nodes_in_cluster <- names(mem)[mem == cid]
    sub_g <- induced_subgraph(g, vids = which(V(g)$name %in% nodes_in_cluster))
    E(sub_g)$sign <- ifelse(E(sub_g)$rho > 0, "Positive", "Negative")

    cat(sprintf("Cluster %d: %d nodes, %d edges\n", cid, vcount(sub_g), ecount(sub_g)))

    set.seed(42)
    panels[[as.character(cid)]] <- ggraph(sub_g, layout = "fr") +
      geom_edge_link(aes(color = sign, width = abs(rho)), alpha = 0.5) +
      scale_edge_color_manual(values = c("Positive" = "#B03A2E", "Negative" = "#1F618D"),
                              name = "Correlation") +
      scale_edge_width(range = c(0.3, 2), name = "|rho|") +
      geom_node_point(color = cluster_colors[[as.character(cid)]], size = 6) +
      geom_node_text(aes(label = name), size = 5, repel = TRUE, max.overlaps = 30) +
      labs(title = sprintf("Cluster %d (n=%d compounds)", cid, vcount(sub_g))) +
      theme_void(base_size = 10) +
      theme(plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            plot.title = element_text(face = "bold", hjust = 0.5, size = 11,
                                       color = cluster_colors[[as.character(cid)]]),
            legend.position = "none")
  }

  p_network <- wrap_plots(panels, ncol = 2) +
    plot_annotation(
      title = "GC-HRMS compound co-occurrence network, faceted by Louvain community",
      subtitle = sprintf("|Spearman rho| >= %.2f, FDR-significant, jointly-detected samples only | %d compounds, %d communities (cross-cluster edges not shown)",
                         RHO_THRESH, vcount(g), length(lou)),
      theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))
    )

  print(p_network)
  ggsave("GC_correlation_network_louvain_faceted.png", p_network, width = 18, height = 14, dpi = 600)
  cat("\nSaved: GC_correlation_network_louvain_faceted.png\n")

} else {

  plot_edges <- edges %>% filter(abs(rho) >= PLOT_RHO_THRESH)
  cat("\nEdges retained for plotting at |rho| >=", PLOT_RHO_THRESH, ":", nrow(plot_edges), "\n")

  g_plot <- graph_from_data_frame(
    plot_edges %>% select(Compound_A, Compound_B, rho, padj, n),
    directed = FALSE,
    vertices = data.frame(name = V(g)$name, cluster = V(g)$cluster)
  )

  # Drop isolated nodes (0 edges at PLOT_RHO_THRESH) — otherwise these scatter
  # across the plot as unconnected singletons and add clutter with no
  # information (their correlations exist only below PLOT_RHO_THRESH, and
  # are still captured in the full cluster_table export below).
  n_isolated <- sum(degree(g_plot) == 0)
  cat("Isolated nodes after thinning (0 edges at plot threshold):", n_isolated, "\n")
  g_plot <- delete_vertices(g_plot, which(degree(g_plot) == 0))
  cat("Nodes remaining after dropping isolates:", vcount(g_plot), "\n")

  g <- g_plot

  E(g)$sign <- ifelse(E(g)$rho > 0, "Positive", "Negative")

  set.seed(42)
  p_network <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(color = sign, width = abs(rho)), alpha = 0.4) +
    scale_edge_color_manual(values = c("Positive" = "#B03A2E", "Negative" = "#1F618D"),
                            name = "Correlation") +
    scale_edge_width(range = c(0.3, 2), name = "|rho|") +
    geom_node_point(aes(color = cluster), size = 8, show.legend = FALSE) +
    geom_node_text(aes(label = name), size = 5, repel = TRUE, max.overlaps = 30) +
    labs(title = "GC-HRMS compound co-occurrence network",
         subtitle = sprintf("Edges shown at |Spearman rho| >= %.2f (communities computed at >= %.2f), FDR-significant, isolated nodes hidden | %d edges, %d compounds, %d Louvain communities",
                            PLOT_RHO_THRESH, RHO_THRESH, ecount(g), vcount(g), length(lou))) +
    theme_void(base_size = 11) +
    theme(plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA)) +
    theme(plot.title = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"))

  print(p_network)
  ggsave("GC_correlation_network_louvain.png", p_network, width = 16, height = 12, dpi = 600)
  cat("\nSaved: GC_correlation_network_louvain.png\n")
}

# Cluster membership table (full RHO_THRESH-graph membership, not the
# thinned plot graph)
cluster_table <- data.frame(Compound = names(membership(lou)),
                             Cluster = as.integer(membership(lou))) %>%
  arrange(Cluster, Compound)
write.csv(cluster_table, "GC_correlation_network_louvain_clusters.csv", row.names = FALSE)
cat("Saved: GC_correlation_network_louvain_clusters.csv\n\n")

# ============================================================================
# SECTION 3: Pairwise Fisher's exact co-occurrence tests
# Direct GC counterpart of Section 3 in LC_chemical_correlations.R — same
# method, same MAX_DF_FOR_RANKING logic.
#
# Method notes (identical logic to the LC script, read there for detail):
#  - One-sided Fisher's exact test (alternative = "greater") on each pair's
#    2x2 detection contingency table, across ALL samples.
#  - MIN_DETECTIONS = 3 (absolute count, not %).
#  - MAX_DF_FOR_RANKING = 0.15 — same 15% rarity ceiling used for LC. GC has
#    far fewer samples and a much smaller confirmed compound set than LC —
#    re-check whether 15% is still a sensible ceiling here rather than
#    assuming it transfers directly.
#  - BH-FDR correction applied across all pairwise tests.
#  - Odds ratio can be Infinite (perfect separation) — reported as Inf,
#    not an error.
# ============================================================================

MIN_DETECTIONS     <- 3     # minimum absolute number of detections (>0) required per compound
MAX_DF_FOR_RANKING  <- 0.15  # DF ceiling (both compounds) for the ranked top-40 table only

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
    total = n_both + n_A_only + n_B_only + n_neither,
    DF_A  = (n_both + n_A_only) / total,
    DF_B  = (n_both + n_B_only) / total
  ) %>%
  arrange(p)

write.csv(fisher_results, "GC_all_pairwise_fisher_cooccurrence.csv", row.names = FALSE)
cat("Saved: GC_all_pairwise_fisher_cooccurrence.csv (", nrow(fisher_results), "pairs tested)\n")

n_sig <- sum(fisher_results$padj < ALPHA, na.rm = TRUE)
n_inf <- sum(is.infinite(fisher_results$odds_ratio))
cat("\nPairs significant at padj <", ALPHA, ":", n_sig, "\n")
cat("Pairs with infinite OR (perfect co-occurrence separation):", n_inf, "\n")

top_pairs_fisher <- fisher_results %>%
  filter(padj < ALPHA, DF_A <= MAX_DF_FOR_RANKING, DF_B <= MAX_DF_FOR_RANKING) %>%
  arrange(p) %>%
  slice_head(n = TOP_N)

cat("\nPairs passing DF <=", MAX_DF_FOR_RANKING, "ceiling (both compounds), padj <", ALPHA, ":",
    sum(fisher_results$padj < ALPHA & fisher_results$DF_A <= MAX_DF_FOR_RANKING &
        fisher_results$DF_B <= MAX_DF_FOR_RANKING, na.rm = TRUE), "\n")

cat("\nTop", TOP_N, "rare-pair co-occurrences (DF <=", MAX_DF_FOR_RANKING, "each, padj <", ALPHA, "):\n")
print(top_pairs_fisher %>% select(Compound_A, Compound_B, n_both, DF_A, DF_B, odds_ratio, p, padj))

write.csv(top_pairs_fisher, "GC_top40_fisher_cooccurrence.csv", row.names = FALSE)
cat("Saved: GC_top40_fisher_cooccurrence.csv\n\n")

# ============================================================================
# SECTION 4: Sample-level detection heatmap for rare, co-occurring compounds
# Direct GC counterpart of Section 5 in LC_chemical_correlations.R. GC has
# no separate rare-compound Fisher NETWORK step (unlike LC) — this reads
# directly from Section 3's Fisher results.
#
# Method notes:
#  - Restricted to compounds appearing in >=1 significant pair from
#    Section 3 (DF <= MAX_DF).
#  - Row/column order = hierarchical clustering on Jaccard distance
#    (dist(..., method="binary")).
#  - Country parsed from GC sample names (pseudonymized pattern
#    "{Country}_{Indoor|Outdoor}_{code}"; SL -> SI remap kept as a safeguard).
#  - GROUP_LOOKUP below reflects named clusters identified by inspecting
#    GC_top40_fisher_cooccurrence.csv (padj-ranked pairs) — chemically
#    plausible but not yet checked for indoor/outdoor deployment
#    consistency (the step that validated/corrected LC's named groups,
#    e.g. Urocanic acid's reversed pattern). Update this list if new named
#    clusters emerge from a fresh run.
# ============================================================================

MAX_DF <- 0.15   # same rarity ceiling as Section 3
SHOW_COLUMN_DENDROGRAM <- TRUE

det_mat_full <- (conc_mat > 0) & !is.na(conc_mat)

fisher_results_heatmap <- read.csv("GC_all_pairwise_fisher_cooccurrence.csv")
rare_sig <- fisher_results_heatmap %>%
  filter(padj < ALPHA, DF_A <= MAX_DF, DF_B <= MAX_DF)

rare_compounds <- unique(c(rare_sig$Compound_A, rare_sig$Compound_B))
cat("Rare compounds in at least one significant pair:", length(rare_compounds), "\n")

if (length(rare_compounds) == 0) {
  stop("No compounds passed the significance/DF filter — check Section 3's ",
       "output before running this section; GC's smaller sample size may need a looser MAX_DF.")
}

det_mat_heatmap <- det_mat_full[rownames(det_mat_full) %in% rare_compounds, , drop = FALSE]
storage.mode(det_mat_heatmap) <- "numeric"

any_detected <- colSums(det_mat_heatmap) > 0
det_mat_active <- det_mat_heatmap[, any_detected, drop = FALSE]
cat("Samples with at least one rare-compound detection:", ncol(det_mat_active), "of", ncol(det_mat_heatmap), "\n")

sample_country_raw <- str_match(colnames(det_mat_active), "^([A-Z]{2})_")[, 2]
sample_country <- ifelse(sample_country_raw == "SL", "SI", sample_country_raw)
col_annotation <- data.frame(Country = sample_country)
rownames(col_annotation) <- colnames(det_mat_active)

country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C"
)
ann_colors <- list(Country = country_colors[names(country_colors) %in% unique(sample_country)])

GROUP_LOOKUP <- list(
  # "Group Name" = c("Compound A", "Compound B", ...)
  "Glycol ether solvents" = c(
    "2-(2-Methoxypropoxy)propan-1-ol",
    "2-Ethylhexan-1-ol",
    "2-Phenoxyethan-1-ol",
    "1-[(1-Methoxypropan-2-yl)oxy]propan-2-ol",
    "1-[1-(1-butoxypropan-2-yloxy)propan-2-yloxy]propan-2-ol"
  ),
  "Salicylate esters" = c(
    "Benzyl 2-hydroxybenzoate",
    "Cyclohexyl 2-hydroxybenzoate",
    "decan-2-yl benzoate",
    "Pentyl 2-hydroxybenzoate"
  )
)

group_lookup <- rep("Other", nrow(det_mat_active))
names(group_lookup) <- rownames(det_mat_active)
for (grp_name in names(GROUP_LOOKUP)) {
  group_lookup[names(group_lookup) %in% GROUP_LOOKUP[[grp_name]]] <- grp_name
}

row_annotation <- data.frame(Group = group_lookup)
rownames(row_annotation) <- rownames(det_mat_active)

group_colors <- c("Other" = "grey85")
if (length(GROUP_LOOKUP) > 0) {
  palette_extra <- c("#B03A2E", "#1F618D", "#F5B041", "#7B3F9E", "#28A745", "#0AA0BF")
  for (i in seq_along(GROUP_LOOKUP)) {
    group_colors[names(GROUP_LOOKUP)[i]] <- palette_extra[((i - 1) %% length(palette_extra)) + 1]
  }
}
ann_colors$Group <- group_colors

display_names <- rownames(det_mat_active)
is_long <- nchar(display_names) > 40
display_names[is_long] <- paste0(substr(display_names[is_long], 1, 37), "...")
rownames(det_mat_active)  <- display_names
rownames(row_annotation)  <- display_names

p_heatmap <- pheatmap(
  det_mat_active,
  color = c("white", "#B03A2E"),
  border_color = NA,   # GC's matrix has far fewer rows/cols than LC's, so
                       # cells render much larger on the same canvas size —
                       # pheatmap's default grey60 cell border becomes
                       # visually dominant and reads as a grey background.
  legend_breaks = c(0, 1),
  legend_labels = c("Not detected", "Detected"),
  clustering_distance_rows = "binary",
  clustering_distance_cols = "binary",
  clustering_method = "ward.D2",
  annotation_row = row_annotation,
  annotation_col = col_annotation,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  treeheight_col = if (SHOW_COLUMN_DENDROGRAM) 50 else 0,
  fontsize_row = 6,
  main = sprintf("Sample-level detection of rare co-occurring GC-HRMS compounds (DF <= %.0f%%, padj < %.2f in >=1 pair)\nRows/columns clustered by Jaccard distance",
                  MAX_DF * 100, ALPHA),
  filename = "GC_rare_compound_sample_heatmap.png",
  width = 18, height = 11, dpi = 300
)

cat("\nSaved: GC_rare_compound_sample_heatmap.png\n")

write.csv(det_mat_active, "GC_rare_compound_detection_matrix.csv")
write.csv(data.frame(Compound = rownames(row_annotation), Group = row_annotation$Group),
          "GC_rare_compound_groups.csv", row.names = FALSE)
write.csv(data.frame(SampleID = rownames(col_annotation), Country = col_annotation$Country),
          "GC_rare_compound_sample_countries.csv", row.names = FALSE)
cat("Saved: GC_rare_compound_detection_matrix.csv, GC_rare_compound_groups.csv, GC_rare_compound_sample_countries.csv\n\n")

cat("==== DONE — GC chemical correlation analysis complete ====\n")
