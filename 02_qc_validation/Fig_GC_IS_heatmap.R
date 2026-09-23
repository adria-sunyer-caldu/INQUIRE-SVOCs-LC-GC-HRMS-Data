# =============================================================================
# Fig_GC_IS_heatmap.R
#
# GC-HRMS counterpart of Fig_IS_heatmap_FINAL.R (LC). Same layout, colors,
# title, axis labels and RSD rule as the LC script.
#
# Heatmap: RSD% per internal standard (row) x batch (column), using the
# 4CPS_PDMS QC replicate injections only (QC_NORM equivalent for GC).
# 8 internal standards: 13C12 PBDE-183 was removed from the input file,
# being gap-filled in 50% of the injections (below the >10% threshold used
# for LC-HRMS) and never the closest-eluting IS for any feature.
#
# BATCH ASSIGNMENT (replaces the earlier proportional allocation):
# In IS_Area, columns follow the acquisition sequence: each country is one
# contiguous block and the 23 QC injections (4CPS_PDMS_01-23, numbered in
# injection order) sit inside those blocks. Each QC is assigned to the batch
# of the household sample injected immediately BEFORE it. Batch IDs are
# numbered 1-8 in run order (country is not shown in the figure).
# Two QCs sit exactly at a batch boundary (after the last sample of one
# country, before the first sample of the next): 4CPS_PDMS_09 and
# 4CPS_PDMS_15. Under the "preceding sample" rule they go to the earlier
# batch. Override in `boundary_override` below if the true sequence differs.
#
# GAP-FILLED VALUES: in IS_Area each IS has one repeated constant equal to
# its minimum value. These are MZmine gap-fill values (confirmed by NILU),
# not measured areas, so they are EXCLUDED from the RSD rather than counted
# as signal. A cell with fewer than two measured QC injections gets
# RSD = NA -> grey tile.
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)

base_dir <- "path/to/data"
out_dir  <- base_dir
xlsx_path <- file.path(base_dir, "20260823_INQUIRE_NILU_IS-RT-FB_data_8IS.xlsx")

boundary_override <- c()   # e.g. c("4CPS_PDMS_09" = "4") to move a boundary QC

# ---- Read IS areas (wide: IS x injection) ----
is_raw <- read_excel(xlsx_path, sheet = "IS_Area", .name_repair = "minimal")
all_cols <- names(is_raw)[-(1:4)]   # Name, Base_Peak, RT, RI come first

is_hh <- grepl("^[A-Z]{2}_HH_[0-9]+_(IS[0-9]|OS1)$", all_cols)
is_qc <- grepl("^4CPS_PDMS_[0-9]+$", all_cols)
dropped <- all_cols[!is_hh & !is_qc]
if (length(dropped) > 0) {
  cat(sprintf("NOTE: %d column(s) with no household/QC name dropped: %s\n",
              length(dropped), paste0("'", dropped, "'", collapse = ", ")))
}

# ---- Batch per column from acquisition sequence ----
seq_df <- data.frame(Column = all_cols, Position = seq_along(all_cols),
                     Country = ifelse(is_hh, substr(all_cols, 1, 2), NA),
                     stringsAsFactors = FALSE) %>%
  filter(is_hh | is_qc)
# run order of countries = order of first appearance
country_run_order <- unique(na.omit(seq_df$Country))
# QC inherits the country of the preceding household sample
seq_df$Batch_country <- seq_df$Country
for (i in seq_len(nrow(seq_df))) {
  if (is.na(seq_df$Batch_country[i])) {
    prev <- seq_df$Batch_country[seq_len(i - 1)]
    prev <- prev[!is.na(prev)]
    seq_df$Batch_country[i] <- if (length(prev)) tail(prev, 1) else seq_df$Country[!is.na(seq_df$Country)][1]
  }
}
seq_df$Batch_ID <- as.character(match(seq_df$Batch_country, country_run_order))
for (qc in names(boundary_override)) seq_df$Batch_ID[seq_df$Column == qc] <- boundary_override[[qc]]

cat("Batch run order (Batch_ID = country):\n")
print(setNames(seq_along(country_run_order), country_run_order))
cat("\nQC injections per batch:\n")
print(table(seq_df$Batch_ID[grepl("^4CPS_PDMS", seq_df$Column)]))

# ---- Long format, placeholder -> 0 ----
long <- lapply(seq_len(nrow(is_raw)), function(j) {
  v <- as.numeric(unlist(is_raw[j, seq_df$Column]))
  mn <- min(v, na.rm = TRUE)
  gap <- v == mn & sum(v == mn, na.rm = TRUE) >= 2   # repeated minimum = MZmine gap-fill
  v[gap] <- NA
  data.frame(Name = is_raw$Name[j], Ionization = "EI",
             Column = seq_df$Column, Batch_ID = seq_df$Batch_ID, Value = v,
             stringsAsFactors = FALSE)
}) %>% bind_rows()

qc_long <- long %>% filter(grepl("^4CPS_PDMS", Column))

rsd_by_batch <- qc_long %>%
  group_by(Name, Ionization, Batch_ID) %>%
  summarise(Mean = mean(Value, na.rm = TRUE), SD = sd(Value, na.rm = TRUE),
            n_qc = n(), n_measured = sum(!is.na(Value)), .groups = "drop") %>%
  mutate(RSD_pct = ifelse(n_measured >= 2 & !is.na(Mean) & Mean > 0, SD / Mean * 100, NA))

rsd_by_batch <- rsd_by_batch %>%
  mutate(Name_label = paste0(Name, " (", Ionization, ")"))

name_order <- rsd_by_batch %>%
  group_by(Name_label, Ionization) %>% summarise(mean_rsd = mean(RSD_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(Ionization, desc(mean_rsd)) %>% pull(Name_label)
rsd_by_batch$Name_label <- factor(rsd_by_batch$Name_label, levels = rev(name_order))
rsd_by_batch$Batch_ID <- factor(rsd_by_batch$Batch_ID,
                                 levels = as.character(sort(as.numeric(unique(rsd_by_batch$Batch_ID)))))

p_heatmap <- ggplot(rsd_by_batch, aes(x = Batch_ID, y = Name_label, fill = RSD_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(RSD_pct), paste0("n=", n_measured), as.character(round(RSD_pct, 0)))),
            size = 4, color = "black") +
  facet_wrap(~Ionization, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "firebrick", na.value = "grey90", name = "RSD %") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y = element_text(size = 10), panel.grid = element_blank()) +
  labs(x = "Batch ID", y = NULL,
       title = "RSD% per internal standard (QC_NORM replicates only), by batch")
print(p_heatmap)
ggsave(file.path(out_dir, "Fig_GC_IS_heatmap_QCnorm.png"), p_heatmap, width = 9, height = 5, dpi = 600)
cat("Saved: Fig_GC_IS_heatmap_QCnorm.png\n")

# Source Data (country kept here for traceability only; not shown in figure)
batch_key <- data.frame(Batch_ID = as.character(seq_along(country_run_order)),
                        Batch_country = country_run_order)
write.csv(rsd_by_batch %>% mutate(Batch_ID = as.character(Batch_ID)) %>% left_join(batch_key, by = "Batch_ID"),
          file.path(out_dir, "GC_IS_QCnorm_RSD_by_batch.csv"), row.names = FALSE)
write.csv(seq_df, file.path(out_dir, "GC_injection_sequence_batch_assignment.csv"), row.names = FALSE)
cat("Saved underlying data: GC_IS_QCnorm_RSD_by_batch.csv, GC_injection_sequence_batch_assignment.csv\n")
