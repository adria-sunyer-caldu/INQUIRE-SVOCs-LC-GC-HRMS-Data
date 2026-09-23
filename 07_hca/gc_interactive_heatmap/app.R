############################################################
## Interactive Hierarchical Cluster Heatmap – INQUIRE GC-HRMS SVOCs
## Restricted to confirmed compounds (ID level 1-2), main indoor/outdoor samplers only.
## Sample codes are pseudonymized: {Country}_{Indoor|Outdoor}_{code}, e.g. PT_Indoor_AK.
## Mirrors the final working LC-HRMS interactive app (app.R, LC repo).
############################################################

# ---- Load required packages ----
# NOTE: for ShinyProxy deployment, packages must be pre-installed in the Docker
# image (see Dockerfile) — do not call install.packages()/BiocManager::install()
# inside app.R, the container has no write access to the library at runtime.
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(shiny)
library(InteractiveComplexHeatmap)

# ---- Load dataset ----
# File must sit alongside app.R in the same directory that gets copied into the
# Docker image (see Dockerfile) — relative path, no local machine path.
gc_file <- "GC_Level1_2_polished_FINAL.csv"

# ---- Load and filter to confirmed compounds (ID level 1-2) ----
gc_raw <- read.csv(gc_file, fileEncoding = "UTF-8", header = TRUE, check.names = FALSE)
gc_raw <- gc_raw[, !(is.na(names(gc_raw)) | names(gc_raw) == "")]  # drop stray blank/unnamed columns

meta_cols <- c("RT","RT2","Feature_ID","Database match","Formula","InChiKey","Smiles","CAS",
               "IUPAC name","Method","m/z","Ion","MS/MS","MW","Exact Mass","ID level",
               "Componant info","RT/RI LC-MS (Other Index), (optional)","RT/RI UoA","RI",
               "LoD","LoQ","Uncertainty","Quantification Method (optional)",
               "Notes on Batch Effects (optional)","Concentration Type")
meta_cols <- intersect(meta_cols, names(gc_raw))

gc_confirmed <- gc_raw %>% filter(`ID level` %in% c(1, 2))
cat("Confirmed (ID level 1-2) features:", nrow(gc_confirmed), "\n")

# Use Database match (chemical name) as the row label, fall back to Feature_ID if missing
gc_confirmed$RowLabel <- ifelse(!is.na(gc_confirmed$`Database match`) & gc_confirmed$`Database match` != "",
                                 gc_confirmed$`Database match`, gc_confirmed$Feature_ID)
if (any(duplicated(gc_confirmed$RowLabel))) gc_confirmed$RowLabel <- make.unique(gc_confirmed$RowLabel)
rownames(gc_confirmed) <- gc_confirmed$RowLabel

# ---- Sample columns: main indoor/outdoor samplers only (extra CZ _2/_3 samplers excluded) ----
sample_cols_all <- setdiff(names(gc_raw), meta_cols)
sample_cols <- sample_cols_all[grepl("^[A-Z]{2}_(Indoor|Outdoor)_[A-Z]{2}$", sample_cols_all)]
cat("Indoor/outdoor sample columns:", length(sample_cols), "\n")

gc_confirmed[sample_cols] <- lapply(gc_confirmed[sample_cols], as.numeric)
mat <- as.matrix(gc_confirmed[, sample_cols])
rownames(mat) <- gc_confirmed$RowLabel

# ---- Detection frequency (row metadata) ----
detection_frequency <- rowSums(mat > 0, na.rm = TRUE) / ncol(mat) * 100

# ---- Z-score normalization per feature ----
cal_z_score <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
mat_norm <- t(apply(mat, 1, cal_z_score))
mat_norm[is.nan(mat_norm)] <- 0   # guards against zero-variance rows

# ---- Column ordering: Country, then Indoor/Outdoor within country ----
# Fixed country order — MUST match HCA_Heatmap_LC.R exactly for cross-platform comparability.
country_order <- c("IT", "NL", "UK", "SE", "EE", "PT", "CZ", "SI")

country_from_col <- sub("_.*", "", sample_cols)
io_from_col <- sub("^[A-Z]{2}_(Indoor|Outdoor)_.*$", "\\1", sample_cols)
col_df <- data.frame(sample = sample_cols, country = country_from_col, io = io_from_col, stringsAsFactors = FALSE)
col_df$country <- factor(col_df$country, levels = country_order)
col_order <- col_df %>% arrange(country, match(io, c("Indoor", "Outdoor")), sample) %>% pull(sample)   # Indoor before Outdoor within each country block
mat_norm <- mat_norm[, col_order]

Country <- col_df$country[match(col_order, col_df$sample)]
IndoorOutdoor <- col_df$io[match(col_order, col_df$sample)]

# ---- Custom colors (locked project palette; same as LC) ----
country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C"
)
country_colors <- country_colors[names(country_colors) %in% unique(Country)]  # keep only countries present
inout_colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")

col_ha <- HeatmapAnnotation(
  Country = Country,
  Sampling = IndoorOutdoor,
  col = list(Country = country_colors, Sampling = inout_colors),
  annotation_legend_param = list(Sampling = list(title = "Sampling"))
)

# ---- Row annotation: ID level (confidence), Compound_Class, Detection Frequency ----
class_file <- "GC_Compound_Class.csv"
if (file.exists(class_file)) {
  gc_class <- read.csv(class_file, stringsAsFactors = FALSE)
  gc_confirmed <- gc_confirmed %>% left_join(gc_class %>% select(Feature_ID, Compound_Class), by = "Feature_ID")
  gc_confirmed$Compound_Class[is.na(gc_confirmed$Compound_Class)] <- "Other"
  cat("Compound_Class merged from", class_file, "\n")
} else {
  gc_confirmed$Compound_Class <- "Not classified"
  cat("WARNING:", class_file, "not found — run GC_SMARTS_classifier.py first. Class annotation will be blank.\n")
}

class_levels <- sort(unique(gc_confirmed$Compound_Class))
class_palette <- c("#E69F00","#56B4E9","#001219","#c9184a","#0072B2","#D55E00","#CC79A7","#999999",
                    "#9C179E","#3B528B","#FDE725","#1F968B","#DE8C00","#482878","#F94144","#43AA8B",
                    "#F3722C","#277DA1","#90BE6D","#577590")
class_colors <- setNames(class_palette[seq_along(class_levels) %% length(class_palette) + 1], class_levels)

row_ha <- rowAnnotation(
  Confidence = paste0("Level ", gc_confirmed$`ID level`[match(rownames(mat_norm), gc_confirmed$RowLabel)]),
  Compound_Class = gc_confirmed$Compound_Class[match(rownames(mat_norm), gc_confirmed$RowLabel)],
  DetectionFreq = detection_frequency[match(rownames(mat_norm), names(detection_frequency))],
  col = list(
    Confidence = c("Level 1" = "black", "Level 2" = "deepskyblue"),
    Compound_Class = class_colors,
    DetectionFreq = colorRamp2(c(0, 100), c("white", "darkgreen"))
  )
)

# ---- Heatmap colors ----
mycols <- colorRamp2(c(-10, -5, 0, 5, 10), c("darkblue", "blue", "white", "red", "darkred"))

# ---- Create the main heatmap ----
ht <- Heatmap(
  mat_norm,
  name = "Z-score",
  cluster_columns = FALSE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "average",
  top_annotation = col_ha,
  right_annotation = row_ha,
  column_title = NULL,
  row_title = "Chemical Features (GC-HRMS, confirmed, ID level 1-2)",
  show_row_names = FALSE,
  show_column_names = FALSE,
  col = mycols
)

# ---- Click action ----
click_action <- function(df, output) {
  output$info <- renderPrint({
    str(df)
  })
}

# ---- Shiny UI ----
ui <- fluidPage(
  titlePanel("GC-HRMS Interactive Heatmap of Confidently Annotated Chemicals"),
  p("The INQUIRE study is a large-scale exposomics monitoring campaign of semi-volatile organic compounds (SVOCs) in indoor and outdoor air, sampled across 8 European countries (CZ, EE, IT, NL, PT, SE, SI, UK). This app shows the confirmed and confidently annotated GC-HRMS chemical features detected across the study. Rows are chemical features (z-scored intensity), annotated by confidence level, compound class, and detection frequency; columns are samples, annotated by country and indoor/outdoor sampling. Click and drag on the original heatmap to zoom into a region in the sub-heatmap panel on the right. You may need to widen your browser window first."),
  InteractiveComplexHeatmapOutput(
    width1 = 800, height1 = 800,
    width2 = 800, height2 = 800
  ),
  verbatimTextOutput("info")
)

# ---- Shiny server ----
server <- function(input, output, session) {
  makeInteractiveComplexHeatmap(
    input, output, session,
    ht,
    click_action = click_action
  )
}

# ---- Run the Shiny app ----
shinyApp(ui, server)
