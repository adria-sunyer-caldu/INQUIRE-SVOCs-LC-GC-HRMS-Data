# INQUIRE-SVOCs-LC-GC-HRMS-Data

Target and nontarget LC-HRMS/GC-HRMS analysis of semi-volatile organic
compounds (SVOCs) in indoor/outdoor air across Europe, using PDMS foam
passive samplers — INQUIRE project.

This repository contains the full analysis pipeline behind the associated
manuscript (citation to be added upon publication), organized by pipeline
stage so each folder can be run largely independently once its inputs are
available.

## Data availability

Two files are deposited on Zenodo (**https://doi.org/10.5281/zenodo.22113617**), one per
platform, each holding the full nontarget feature list:

| Zenodo file | Content |
|---|---|
| `LC_All_Features_FULL_final_zenodo.xlsx` | LC-HRMS, all aligned nontarget features, sample columns only. Sheet `Full_dataset_INQUIRE_nontarget_`, MS-DIAL layout: 3 metadata rows (Class, Sample Type, Injection Order), header in row 4, features from row 5, sample columns from column 33 |
| `GC_Level1_5_polished_FULL_zenodo.xlsx` | GC-HRMS, all nontarget features (ID levels 1-5), one row per feature, header in row 1, 13 metadata columns followed by the 419 sample columns |

Everything else a script reads is either written by an earlier script in
this repository, or a subset/reshaping of one of those two datasets, or a
file produced by an external tool. Nothing else is deposited. The table
below says which is which; scripts whose inputs are marked "not deposited"
cannot be rerun from the Zenodo deposit alone, and those inputs are
available from the corresponding author on reasonable request.

| Script input | Source |
|---|---|
| `LC_All_Features_FULL_final_zenodo.xlsx` | Zenodo (LC) |
| `GC_Level1_5_polished_FULL_zenodo.xlsx` | Zenodo (GC) |
| `LC_Level1_2_polished_FINAL.xlsx`, `GC_Level1_2_polished_FINAL.xlsx` | Written by `01_datasets_preparation/Polish_LC_GC_Level1_2_datasets.R` (as `*_polished.xlsx`; renamed after manual review of the duplicate-InChIKey flags) |
| `GC_Level1_5_polished_FULL_final.xlsx` | Written by the same script (`GC_Level1_5_polished_FULL.xlsx`): annotated GC features only, ID levels 1-5 |
| `Full dataset INQUIRE nontarget without removing adducts.csv` | LC MS-DIAL alignment before adduct/in-source-fragment removal; not deposited |
| `AdductFinder_removed_features_*.xlsx`, `InSourceFinder_removed_features_*.xlsx` | MS-DIAL AdductFinder / InSourceFinder outputs, POS and NEG; not deposited |
| `Targets and annotations concentrations.xlsx` | LC target and annotation table with calibrated concentrations, restricted to confirmed compounds; not deposited |
| `RT_targets_annotations_POS_and_NEG_FINAL.txt`, `Blanks_targets_annotations_POS_and_NEG_FINAL.txt` | Retention-time and field-blank exports for the same LC features, merged and renamed by `01_datasets_preparation/Rename_sample_columns_RT_Blanks.R`; not deposited |
| `PCA-like dataset_FINAL.csv`, `Targets and annotations_FINAL.csv` | The LC dataset transposed to one row per sample, with `SampleID`, `Country` and `Deployment` columns: all features, and confirmed targets/annotations only; not deposited |
| `NAP_xref_all_features.csv` | GNPS Network Annotation Propagation results for the LC features; not deposited |
| `OPERA_targets_results.csv`, `OPERA_GC_results.csv` | OPERA property predictions for the confirmed LC and GC compounds; not deposited |
| `Toxicity data targets and annotations.csv`, `GC_159_compounds_for_DTXSID_lookup.csv` | Compound lists used for the hazard lookups; not deposited |
| `20260823_INQUIRE_NILU_IS-RT-FB_data_8IS.xlsx` | GC-HRMS QA/QC export from NILU (internal standard areas, retention times, field blanks), with 13C12 PBDE-183 removed; not deposited |
| `pikme_all.parquet` | The pikme hazard database, a third-party dataset obtained separately |

Scripts expect their inputs in one local working directory. Each script has
a `data_dir <- "path/to/data"` (or similarly named) variable near the top
that you should point at wherever you have placed the data.

## Sample codes

Sample names in the deposited datasets are pseudonymized as
`{Country}_{Indoor|Outdoor}_{code}` (e.g. `PT_Indoor_AK`). The indoor and
outdoor samples of the same household share the same two-letter code. Extra
indoor samplers in five Czech households carry a `_2` or `_3` suffix
(e.g. `CZ_Indoor_AL_2`) and have no outdoor partner. All analysis scripts
from `02_qc_validation` onward read these codes. The scripts in
`01_datasets_preparation` run on the raw instrument exports, which use the
original internal sample names and are available upon reasonable request;
sample names were pseudonymized after that step.

## Repository structure and run order

Folders are numbered in roughly the order they depend on each other's
outputs. Within `09_toxicity`, the pikme-based scripts also require a
separate hazard database (pikme) subset — see comments inside those
scripts for details.

| Folder | Contents | Depends on |
|---|---|---|
| `01_datasets_preparation` | Rename/merge raw RT & field-blank files (POS+NEG); remove adduct/in-source-fragment features; produce the final polished Level 1–2 (confirmed) and Level 1–5 (all annotation levels) LC/GC datasets | Raw MS-DIAL exports (Zenodo) |
| `02_qc_validation` | RT stability, field blank levels, internal-standard normalization QC for both platforms (RSD heatmap + before/after normalization; LC scripts and their GC counterparts, `Fig_GC_IS_heatmap.R` and `Fig_GC_Before_After_Normalization.R`) | `01_datasets_preparation` outputs; GC QA/QC export |
| `03_total_burden` | Total chemical burden per sample/household, LC and GC, by country; fold-change and raincloud summary figures | `01_datasets_preparation` outputs |
| `04_target_analysis` | Compound-class indoor/outdoor enrichment, LC-vs-GC identity overlap, merged class-panel detection-frequency figures, named-compound violin plots | `01_datasets_preparation` outputs |
| `05_nontarget_analysis` | Full nontarget feature prevalence/prioritization (paired indoor/outdoor ratio), Kendrick mass defect, NAP-based molecular characterization (Van Krevelen, DBE, halogenation) | `01_datasets_preparation` outputs |
| `06_pca` | PCA of LC/GC data, both for confirmed targets/annotations only and for the full nontarget feature set | `01_datasets_preparation` outputs |
| `07_hca` | Interactive hierarchical clustering heatmaps (Shiny apps), LC and GC. Snapshot of the deployed apps; live versions hosted on SciLifeLab Serve | `01_datasets_preparation` outputs |
| `08_correlations` | Pairwise Spearman correlation and Fisher's exact co-occurrence networks between compounds, plus sample-level detection heatmaps for rare co-occurring compounds | `01_datasets_preparation` outputs |
| `09_toxicity` | Hazard/EDC characterization using the pikme database (DTXSID retrieval, hazard-score merging, indoor/outdoor hazard statistics and figures) | `01_datasets_preparation` outputs + pikme database subset (Zenodo) |
| `10_annotation_confidence` | Annotation confidence-level (Schymanski/Koelmel scale) distribution, LC and GC | `01_datasets_preparation` outputs + NAP cross-reference file |
| `utils` | Standalone Python helpers for retrieving InChIKeys/CAS/SMILES/molecular formula from PubChem, used during compound identification | none (run independently as needed) |

## Interactive visualizations

Live, interactive versions of the hierarchical clustering heatmaps:
- LC-HRMS: https://lc-heatmap.serve.scilifelab.se/app/lc-heatmap
- GC-HRMS: https://gc-heatmap.serve.scilifelab.se/app/gc-heatmap

## Requirements

**R** (≥ 4.x recommended). Key packages used across scripts:
`data.table`, `dplyr`, `tidyr`, `stringr`, `ggplot2`, `readxl`, `writexl`,
`scales`, `patchwork`, `igraph`, `ggraph`, `pheatmap`, `ComplexHeatmap`,
`InteractiveComplexHeatmap`, `circlize`, `arrow`, `ctxR`, `shiny`.
Install via `install.packages()`; a few scripts install missing packages
automatically on first run.

**Python** (≥ 3.9) for the `utils` scripts: `pandas`, `requests`, `openpyxl`.

## Notes on reproducibility

- Every script's working/data directory is a placeholder
  (`"path/to/data"`) — edit this at the top of each script before running.
- Scripts are designed to be run top-to-bottom within their file; several
  read a previous script's CSV output rather than recomputing it, so
  running out of order within a folder may fail with a missing-file error.
- One script (`retrieve_DTXSID_GC` step, folder `09_toxicity`) requires a
  free EPA CTX API key, supplied via the environment variable
  `CTX_API_KEY` (never hardcoded in the script itself).

## License

All content in this repository — code and data alike — is released under
the [Creative Commons Attribution 4.0 International License (CC-BY 4.0)](https://creativecommons.org/licenses/by/4.0/).
You are free to share and adapt this material for any purpose, including
commercially, as long as you give appropriate credit.

**Attribution:** Adrià Sunyer Caldú, Stockholm University — INQUIRE project.

Data deposited separately on Zenodo is likewise released under CC-BY 4.0.

## Citation

If you use this code or data, please give credit to:
Adrià Sunyer Caldú, Stockholm University — [paper citation to be added upon publication]
Code archive DOI (all versions): https://doi.org/10.5281/zenodo.22924214

## Contact

Adrià Sunyer Caldú — Stockholm University
