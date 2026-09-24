# =============================================================================
# Fig_GC_IS_heatmap.R
#
# GC-HRMS internal standard reproducibility. Heatmap: RSD% per internal
# standard (row) x country batch (column, x axis shows country codes),
# from the 4CPS_PDMS QC replicate injections.
#
# BATCH ASSIGNMENT: In IS_Area, columns follow the acquisition sequence,
# with each country forming one contiguous block. Each QC injection is
# assigned the country of the household sample injected immediately before
# it. Override a specific QC via `boundary_override` below if needed.
#
# GAP-FILLED VALUES: the per-IS repeated minimum constant in IS_Area is an
# MZmine gap-fill value, not a measured area, and is excluded from the RSD.
# A cell with fewer than two measured QC injections shows n=<count> instead
# of an RSD and is shaded grey.
#
# 8 internal standards: 13C12 PBDE-183 was removed from the input file
# (gap-filled in 50% of injections, and never the closest-eluting IS for
# any feature, since it elutes after every feature in the dataset).
#
# Source data output matches the "Suppl. Figure 41" sheet format:
# Name, Ionization, Batch_ID, Mean_area, SD_area, n_QC_injections,
# n_measured, RSD_pct.
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)

base_dir  <- "path/to/data"
out_dir   <- base_dir
xlsx_path <- file.path(base_dir, "20260823_INQUIRE_NILU_IS-RT-FB_data_8IS.xlsx")

boundary_override <- c()   # e.g. c("4CPS_PDMS_09" = "IT") to move a boundary QC

# ---- Read IS areas ----
is_raw   <- read_excel(xlsx_path, sheet = "IS_Area", .name_repair = "minimal")
all_cols <- names(is_raw)[-(1:4)]

is_hh <- grepl("^[A-Z]{2}_HH_[0-9]+_(IS[0-9]|OS1)$", all_cols)
is_qc <- grepl("^4CPS_PDMS_[0-9]+$", all_cols)
dropped <- all_cols[!is_hh & !is_qc]
if (length(dropped) > 0)
  cat(sprintf("NOTE: %d column(s) dropped: %s\n", length(dropped),
              paste0("'", dropped, "'", collapse = ", ")))

# ---- Batch assignment (preceding sample rule -> country code) ----
seq_df <- data.frame(Column   = all_cols,
                     Position = seq_along(all_cols),
                     Country  = ifelse(is_hh, substr(all_cols, 1, 2), NA),
                     stringsAsFactors = FALSE) %>%
  filter(is_hh | is_qc)

country_run_order <- unique(na.omit(seq_df$Country))

seq_df$Batch_country <- seq_df$Country
for (i in seq_len(nrow(seq_df))) {
  if (is.na(seq_df$Batch_country[i])) {
    prev <- seq_df$Batch_country[seq_len(i - 1)]
    prev <- prev[!is.na(prev)]
    seq_df$Batch_country[i] <- if (length(prev)) tail(prev, 1) else
      seq_df$Country[!is.na(seq_df$Country)][1]
  }
}
seq_df$Batch_ID <- as.character(match(seq_df$Batch_country, country_run_order))
for (qc in names(boundary_override))
  seq_df$Batch_ID[seq_df$Column == qc] <- boundary_override[[qc]]

cat("Country run order:\n")
print(setNames(seq_along(country_run_order), country_run_order))
cat("\nQC injections per country:\n")
print(table(seq_df$Batch_country[grepl("^4CPS_PDMS", seq_df$Column)]))

# ---- Long format; gap-filled -> NA ----
long <- lapply(seq_len(nrow(is_raw)), function(j) {
  v  <- as.numeric(unlist(is_raw[j, seq_df$Column]))
  mn <- min(v, na.rm = TRUE)
  v[v == mn & sum(v == mn, na.rm = TRUE) >= 2] <- NA
  data.frame(Name         = is_raw$Name[j],
             Ionization   = "EI",
             Column       = seq_df$Column,
             Batch_ID     = seq_df$Batch_ID,
             Batch_country = seq_df$Batch_country,
             Value        = v,
             stringsAsFactors = FALSE)
}) %>% bind_rows()

qc_long <- long %>% filter(grepl("^4CPS_PDMS", Column))

# ---- RSD per IS per country ----
rsd_by_batch <- qc_long %>%
  group_by(Name, Ionization, Batch_country) %>%
  summarise(Mean_area       = mean(Value, na.rm = TRUE),
            SD_area         = sd(Value,   na.rm = TRUE),
            n_QC_injections = n(),
            n_measured      = sum(!is.na(Value)),
            .groups = "drop") %>%
  mutate(RSD_pct = ifelse(n_measured >= 2 & !is.na(Mean_area) & Mean_area > 0,
                          SD_area / Mean_area * 100, NA),
         Name_label = paste0(Name, " (EI)"))

# Factor ordering
name_order <- rsd_by_batch %>%
  group_by(Name_label) %>%
  summarise(mean_rsd = mean(RSD_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_rsd)) %>% pull(Name_label)
rsd_by_batch$Name_label   <- factor(rsd_by_batch$Name_label, levels = rev(name_order))
rsd_by_batch$Batch_country <- factor(rsd_by_batch$Batch_country, levels = country_run_order)

# ---- Figure ----
p_heatmap <- ggplot(rsd_by_batch,
                    aes(x = Batch_country, y = Name_label, fill = RSD_pct)) +
  geom_tile(color = "white") +
  geom_text(aes(label = ifelse(is.na(RSD_pct),
                               paste0("n=", n_measured),
                               as.character(round(RSD_pct, 0)))),
            size = 4, color = "black") +
  facet_wrap(~Ionization, scales = "free_y") +
  scale_fill_gradient(low = "white", high = "firebrick",
                      na.value = "grey90", name = "RSD %") +
  theme_minimal(base_size = 10) +
  theme(axis.text.y  = element_text(size = 10),
        panel.grid   = element_blank()) +
  labs(x     = "Country",
       y     = NULL,
       title = "RSD% per internal standard (QC replicates only), by country batch")

print(p_heatmap)
ggsave(file.path(out_dir, "Fig_GC_IS_heatmap_QCnorm.png"),
       p_heatmap, width = 9, height = 5, dpi = 600)
cat("Saved: Fig_GC_IS_heatmap_QCnorm.png\n")

# ---- Source data — matches Suppl. Figure 41 sheet format exactly ----
source_data_41 <- rsd_by_batch %>%
  rename(Batch_ID = Batch_country) %>%
  select(Name, Ionization, Batch_ID, Mean_area, SD_area,
         n_QC_injections, n_measured, RSD_pct)

write.csv(source_data_41,
          file.path(out_dir, "GC_IS_QCnorm_RSD_by_country.csv"),
          row.names = FALSE)
cat("Saved: GC_IS_QCnorm_RSD_by_country.csv\n")
cat("  -> Replace sheet 'Suppl. Figure 41' in Source_Data_INQUIRE_FINAL.xlsx with this file.\n")
