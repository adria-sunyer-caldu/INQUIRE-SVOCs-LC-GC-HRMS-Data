library(ggplot2)
library(dplyr)

# ==========================================================
# Lollipop/forest-plot comparing indoor/outdoor total chemical
# burden fold change by country, LC and GC side by side on the
# same log-x axis, significance marked.
#
# Reads directly from the two by-country statistics tables
# produced by the total chemical burden scripts (LC:
# TotalBurden_Statistics_ByCountry.csv, GC:
# GC_TotalBurden_Statistics_ByCountry.csv) rather than hardcoding
# values, so this plot always reflects the current data.
# ==========================================================

# ---- Set your data directory here ----
data_dir <- "path/to/data"

lc_stats <- read.csv(file.path(data_dir, "TotalBurden_Statistics_ByCountry.csv"), stringsAsFactors = FALSE)
gc_stats <- read.csv(file.path(data_dir, "GC_TotalBurden_Statistics_ByCountry.csv"), stringsAsFactors = FALSE)

fc_data <- bind_rows(
  lc_stats %>% transmute(Country, Platform = "LC-HRMS", FC = Median_FC, Sig),
  gc_stats %>% transmute(Country, Platform = "GC-HRMS", FC = Median_FC, Sig)
)

cat("Fold-change data loaded from source CSVs:\n")
print(fc_data)

# Order countries by LC fold-change (descending)
country_order <- lc_stats$Country[order(-lc_stats$Median_FC)]
fc_data$Country <- factor(fc_data$Country, levels = rev(country_order))

# Manual y-offset instead of position_dodge (dodge breaks geom_segment when
# endpoints are far apart on x - it dodges start/end independently, producing
# diagonal lines instead of horizontal lollipops)
fc_data$y_num <- as.numeric(fc_data$Country)
fc_data$y_offset <- ifelse(fc_data$Platform == "LC-HRMS", fc_data$y_num + 0.2, fc_data$y_num - 0.2)

p_lollipop <- ggplot(fc_data, aes(x = FC, y = y_offset, color = Platform)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
  geom_segment(aes(x = 1, xend = FC, y = y_offset, yend = y_offset),
               linewidth = 1.3, alpha = 0.7) +
  geom_point(size = 5.5) +
  geom_text(aes(label = Sig, x = FC * ifelse(FC >= 1, 1.08, 1/1.08),
                hjust = ifelse(FC >= 1, 0, 1)),
            size = 5.5, color = "grey20", show.legend = FALSE, fontface = "bold") +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 10),
                labels = c("0.5x", "1x", "2x", "5x", "10x"),
                # Limits set just below the smallest FC in the data and just
                # above the largest, rather than a fixed symmetric range, so
                # the plot always fits whatever the current data actually is.
                limits = c(min(fc_data$FC) * 0.9, max(fc_data$FC) * 1.1)) +
  scale_y_continuous(breaks = seq_along(levels(fc_data$Country)),
                     labels = levels(fc_data$Country)) +
  scale_color_manual(values = c("LC-HRMS" = "#5A2D75", "GC-HRMS" = "#A85A1F")) +
  labs(
    title = "Indoor/Outdoor Fold Change by Country",
    subtitle = "Median fold change (log scale). Significance: * p<0.05, ** p<0.01, *** p<0.001 (BH-corrected paired Wilcoxon)",
    x = "Indoor / Outdoor fold change",
    y = NULL,
    color = NULL
  ) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(size = 18, face = "bold"),
    axis.text.x        = element_text(size = 15),
    axis.title.x       = element_text(size = 17),
    plot.title          = element_text(size = 21, face = "bold"),
    plot.subtitle       = element_text(size = 13, color = "grey40"),
    legend.position      = "top",
    legend.text          = element_text(size = 15),
    legend.title          = element_blank()
  )
print(p_lollipop)

setwd(data_dir)
ggsave("Fig_FoldChange_Lollipop.png", p_lollipop, dpi = 600, width = 8, height = 9)
cat("Saved: Fig_FoldChange_Lollipop.png\n")
