# === 1. Load packages ===
library(readxl)
library(data.table)

# === 2. Define file paths ===
dir_path <- "path/to/data"

file_adduct_POS   <- file.path(dir_path, "AdductFinder_removed_features_20251024_0940_POS.xlsx")
file_insource_POS <- file.path(dir_path, "InSourceFinder_removed_features_20251024_1058_POS.xlsx")
file_adduct_NEG   <- file.path(dir_path, "AdductFinder_removed_features_20251024_1038_NEG.xlsx")
file_insource_NEG <- file.path(dir_path, "InSourceFinder_removed_features_20251024_1058_NEG.xlsx")
big_csv_file      <- file.path(dir_path, "Full dataset INQUIRE nontarget without removing adducts.csv")
output_file       <- file.path(dir_path, "Full_dataset_INQUIRE_nontarget_cleaned.csv")

# === 3. Read features to remove ===
adduct_POS   <- read_excel(file_adduct_POS, sheet = 1)
insource_POS <- read_excel(file_insource_POS, sheet = 1)
adduct_NEG   <- read_excel(file_adduct_NEG, sheet = 1)
insource_NEG <- read_excel(file_insource_NEG, sheet = 1)

# === 4. Extract IDs ===
feature_ids_POS  <- adduct_POS$Feature_ID
insource_ids_POS <- insource_POS[["In Source ID"]]
feature_ids_NEG  <- adduct_NEG$Feature_ID
insource_ids_NEG <- insource_NEG[["In Source ID"]]

all_remove_ids <- unique(c(feature_ids_POS, insource_ids_POS, feature_ids_NEG, insource_ids_NEG))
all_remove_ids <- all_remove_ids[!is.na(all_remove_ids)]
cat("Total features to remove:", length(all_remove_ids), "\n")

# === 5. Read metadata and header lines ===
meta_lines  <- readLines(big_csv_file, n = 3)
header_line <- readLines(big_csv_file, n = 4)[4]
cat("Header line read successfully.\n")

# === 6. Read main data skipping first 4 rows ===
big_data <- fread(big_csv_file, skip = 4, header = FALSE)
colnames(big_data) <- strsplit(header_line, ",")[[1]]
cat("Total rows before cleaning:", nrow(big_data), "\n")

# === 7. Remove rows matching Alignment IDs ===
cleaned_data <- big_data[!(big_data$`Alignment ID` %in% all_remove_ids), ]
cat("Total rows after cleaning:", nrow(cleaned_data), "\n")

# === 8. Write cleaned CSV ===
# First, write metadata and header
writeLines(c(meta_lines, header_line), con = output_file)

# Then append cleaned data (without column names)
fwrite(cleaned_data, output_file, append = TRUE, col.names = FALSE, sep = ",")

cat("✅ Cleaned CSV saved at:", output_file, "\n")
