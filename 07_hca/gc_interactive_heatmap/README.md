# INQUIRE GC-HRMS Heatmap ShinyApp

The INQUIRE study is a large-scale exposomics monitoring campaign of semi-volatile organic compounds (SVOCs) in indoor and outdoor air, sampled across 8 European countries (CZ, EE, IT, NL, PT, SE, SI, UK). This repository contains the code used to generate an interactive, browsable heatmap of the confirmed and confidently annotated **GC-HRMS** chemical features detected across the study, deployed live on [SciLifeLab Serve](https://serve.scilifelab.se/) as supplementary material for the manuscript.

**This repository covers GC-HRMS data only.** A separate companion repository/app exists for LC-HRMS data.

## Live app (GC-HRMS)

https://gc-heatmap.serve.scilifelab.se/app/gc-heatmap

## Contents

- `app.R` — Shiny app source code for the **GC-HRMS** heatmap (data loading, filtering to confirmed compounds ID level 1–2, z-score normalization, ComplexHeatmap construction, InteractiveComplexHeatmap server/UI)
- `Dockerfile` — container build instructions used for deployment on SciLifeLab Serve
- `GC_Level1_2_polished_FINAL.csv` — underlying **GC-HRMS** feature intensity data, restricted to confirmed compounds (ID level 1–2), indoor/outdoor paired samples (main samplers only), with pseudonymized sample codes
- `GC_Compound_Class.csv` — compound class annotations (from SMARTS-based classification), optional but expected alongside app.R

## What it shows

**GC-HRMS** chemical features only, restricted to confirmed compounds (ID level 1–2). Rows are chemical features (z-scored intensity), annotated by confidence level, compound class, and detection frequency. Columns are samples, annotated by country and indoor/outdoor sampling. Click and drag on the original heatmap to zoom into a region in the sub-heatmap panel on the right.

## Running locally

```r
# Requires: shiny, ComplexHeatmap, InteractiveComplexHeatmap, circlize, dplyr, colorspace, readxl
shiny::runApp("app.R")
```

Or via Docker:

```bash
docker build --platform linux/amd64 -t inquire-gc-heatmap .
docker run --rm -p 3838:3838 inquire-gc-heatmap
```

Then open `http://localhost:3838/`.

## Citation

Supplementary material for the INQUIRE exposomics manuscript (indoor/outdoor air chemical exposure across 8 European countries). This repository/app corresponds to the **GC-HRMS** dataset.
