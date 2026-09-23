# INQUIRE LC-HRMS Heatmap ShinyApp

The INQUIRE study is a large-scale exposomics monitoring campaign of semi-volatile organic compounds (SVOCs) in indoor and outdoor air, sampled across 8 European countries (CZ, EE, IT, NL, PT, SE, SI, UK). This repository contains the code used to generate an interactive, browsable heatmap of the confirmed and confidently annotated **LC-HRMS** chemical features detected across the study, deployed live on [SciLifeLab Serve](https://serve.scilifelab.se/) as supplementary material for the manuscript.

**This repository covers LC-HRMS data only.** A separate companion repository/app exists for GC-HRMS data.

## Live app (LC-HRMS)

https://lc-heatmap.serve.scilifelab.se/app/lc-heatmap

## Contents

- `app.R` — Shiny app source code for the **LC-HRMS** heatmap (data loading, z-score normalization, ComplexHeatmap construction, InteractiveComplexHeatmap server/UI)
- `Dockerfile` — container build instructions used for deployment on SciLifeLab Serve
- `INQUIRE_HCA_Heatmap_Areas_FINAL.csv` — underlying **LC-HRMS** feature data (confirmed and annotated compounds, Level 1–2, indoor/outdoor paired samples, pseudonymized sample codes). Values are chromatographic peak areas (not concentrations).

## What it shows

**LC-HRMS** chemical features only. Rows are chemical features (z-scored peak area), annotated by ionization mode, confidence level, chemical class, subclass, and detection frequency. Columns are samples, annotated by country and indoor/outdoor sampling. Click and drag on the original heatmap to zoom into a region in the sub-heatmap panel on the right.

## Running locally

```r
# Requires: shiny, ComplexHeatmap, InteractiveComplexHeatmap, circlize, dplyr, colorspace
shiny::runApp("app.R")
```

Or via Docker:

```bash
docker build --platform linux/amd64 -t inquire-lc-heatmap .
docker run --rm -p 3838:3838 inquire-lc-heatmap
```

Then open `http://localhost:3838/`.

## Citation

Supplementary material for the INQUIRE exposomics manuscript (indoor/outdoor air chemical exposure across 8 European countries). This repository/app corresponds to the **LC-HRMS** dataset.
