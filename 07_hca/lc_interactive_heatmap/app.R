############################################################
## Interactive Hierarchical Cluster Heatmap – INQUIRE SVOCs
## Author: Adapted from Stefano Papazian's workflow
## Date: 2025-10-09
############################################################

# ---- Load required packages ----
# NOTE: for ShinyProxy deployment, packages must be pre-installed in the Docker
# image (see Dockerfile) — do not call install.packages()/BiocManager::install()
# inside app.R, the container has no write access to the library at runtime.
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(colorspace)
library(shiny)
library(InteractiveComplexHeatmap)

# ---- Load dataset ----
# File must sit alongside app.R in the same directory that gets copied into the
# Docker image (see Dockerfile) — relative path, no local machine path.
file_path <- "INQUIRE_HCA_Heatmap_Areas_FINAL.csv"
SVOCs <- read.csv(file_path, fileEncoding = "UTF-8", header = TRUE, check.names = FALSE)

if(any(duplicated(SVOCs$Annotation))){
  print(SVOCs$Annotation[duplicated(SVOCs$Annotation)])
  stop("Duplicate Annotation names found — please fix them first.")
}

row.names(SVOCs) <- SVOCs$Annotation
df_metadata <- SVOCs %>% dplyr::select(Ionization, Confidence_level, class, subclass, detection_frequency)
df_SVOCs <- SVOCs %>% dplyr::select(-(1:7))
SVOCs_mat <- as.matrix(df_SVOCs)

# ---- Z-score normalization per feature ----
cal_z_score <- function(x) (x - mean(x)) / sd(x)
SVOCs_norm <- t(apply(SVOCs_mat, 1, cal_z_score))

# ---- Sample ordering ----
# Sample codes are pseudonymized: {Country}_{Indoor|Outdoor}_{code}, e.g. PT_Indoor_AK.
# Indoor and outdoor samples from the same household share the same code.
# Extra Czech indoor samplers carry a _2 / _3 suffix (e.g. CZ_Indoor_AL_2).
# Order: fixed country order, Indoor block then Outdoor block, codes alphabetical.
country_order <- c("IT", "NL", "UK", "SE", "EE", "PT", "CZ", "SI")
samp <- colnames(SVOCs_norm)
samp_country <- sub("_.*", "", samp)
samp_io <- sub("^[A-Z]{2}_(Indoor|Outdoor)_.*$", "\\1", samp)
stopifnot(all(samp_country %in% country_order), all(samp_io %in% c("Indoor", "Outdoor")))
sample_order <- samp[order(match(samp_country, country_order),
                           match(samp_io, c("Indoor", "Outdoor")),
                           samp)]
SVOCs_norm <- SVOCs_norm[, sample_order]

Country <- sub("_.*", "", colnames(SVOCs_norm))
IndoorOutdoor <- sub("^[A-Z]{2}_(Indoor|Outdoor)_.*$", "\\1", colnames(SVOCs_norm))

# ---- Custom colors ----
country_colors <- c(
  "IT" = "#25998F", "NL" = "#F36E98", "UK" = "#0AA0BF",
  "SE" = "#78B177", "EE" = "#F05006", "PT" = "#F6114A",
  "CZ" = "#9862A2", "SI" = "#FCA00C", "SL" = "#FCA00C"
)
inout_colors <- c("Indoor" = "#B03A2E", "Outdoor" = "#1F618D")
class_colors <- c("#000000", "#d00000", "#fffcf9")
subclass_colors <- c(
  "#E69F00", "#56B4E9", "#001219", "#c9184a", "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#9C179E", "#3B528B", "#FDE725", "#1F968B", "#DE8C00", "#482878", "#F94144", "#43AA8B",
  "#F3722C", "#277DA1", "#90BE6D", "#577590", "#dcdfe5"
)

# ---- Column annotation (Indoor/Outdoor + Country colors) ----
col_ha <- HeatmapAnnotation(
  Country = Country,
  Sampling = IndoorOutdoor,
  col = list(
    Country = country_colors,
    Sampling = inout_colors
  ),
  annotation_legend_param = list(
    Sampling = list(title = "Sampling")
  )
)

# ---- Row annotation ----
class_levels <- sort(unique(df_metadata$class))
subclass_levels <- sort(unique(df_metadata$subclass))

row_ha <- rowAnnotation(
  Ionization = df_metadata$Ionization,
  Confidence = df_metadata$Confidence_level,
  Class = df_metadata$class,
  Subclass = df_metadata$subclass,
  DetectionFreq = df_metadata$detection_frequency,
  col = list(
    Ionization = c("pos" = "blue", "neg" = "red"),
    Confidence = c("Level 1" = "black", "Level 2" = "deepskyblue"),
    DetectionFreq = colorRamp2(c(0, 100), c("white", "darkgreen")),
    Class = structure(class_colors[seq_along(class_levels)], names = class_levels),
    Subclass = structure(subclass_colors[seq_along(subclass_levels)], names = subclass_levels)
  )
)

# ---- Heatmap colors ----
mycols <- colorRamp2(
  breaks = c(-10, -5, 0, 5, 10),
  colors = c("darkblue", "blue", "white", "red", "darkred")
)

# ---- Create the main heatmap ----
ht <- Heatmap(
  SVOCs_norm,
  name = "Z-score",
  cluster_columns = FALSE,
  clustering_distance_rows = "pearson",
  clustering_method_rows = "average",
  top_annotation = col_ha,
  right_annotation = row_ha,
  column_title = NULL,
  row_title = "Chemical Features",
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
  titlePanel("LC-HRMS Interactive Heatmap of Target and Confidently Annotated Chemicals"),
  p("Interactive Supplementary Data for the INQUIRE exposomics manuscript. Hierarchical clustering heatmap of confirmed and confidently annotated LC-HRMS chemical features (indoor and outdoor air) from the INQUIRE exposomics study, across 8 European countries. Rows are chemical features (z-scored intensity), annotated by ionization mode, confidence level, chemical class, subclass, and detection frequency; columns are samples, annotated by country and indoor/outdoor sampling. Click and drag on the original heatmap to zoom into a region in the sub-heatmap panel on the right. You may need to widen your browser window first."),
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
