# main
rm(list = ls())

# reading and explo ####
source("./reqs.R")

# Load data from second sheet (clean log2 TMT data)
excel_path <- ("../CRC_ProteoResistance_data/CPMSF_GPPF-PM-43_CPPF-PM-43_results_20240119.xlsx")

sheet_names <- readxl::excel_sheets(excel_path)

# read all the sheets
all_sheets <- sheet_names %>% 
  set_names() %>% 
  map(~ read_excel(excel_path, sheet = .x) %>%  janitor::clean_names())

print(names(all_sheets))

# Summarize each sheet: name, dimensions, column names (first 10), and head
for (sheet_name in names(all_sheets)) {
  cat("\n=== Sheet:", sheet_name, "===\n")
  df <- all_sheets[[sheet_name]]
  cat("Dimensions:", dim(df)[1], "rows x", dim(df)[2], "columns\n")
  cat("Column names (first 10):", paste0(colnames(df)[1:min(10, ncol(df))], collapse = ", "), "\n")
  cat("Head of data:\n")
  print(head(df, 3))
  cat("\n-----------------------------\n")
}

# laod experess ####
# Load expression data
expr_data <- all_sheets[["Quan_data"]]

# Extract protein annotations
row_annot <- dplyr::select(expr_data, protein_id, gene_name, description)

# Extract expression matrix
expr_matrix <- expr_data %>%
  dplyr::select(-gene_name, -description) %>%
  column_to_rownames("protein_id") %>%
  as.matrix() %>%
  apply(2, as.numeric)  # convert character columns to numeric

# Load sample metadata
labelling_scheme <- all_sheets[["Labelling_Scheme"]]

# Fix sample names to match expression matrix column names
sample_names <- colnames(expr_matrix)

# Use cleaned-up column names as rownames for sample metadata
sample_meta <- tibble(
  sample = sample_names,
  condition = case_when(
    str_detect(sample, "HCT116_parental") ~ "HCT116_parental",
    str_detect(sample, "SW620_parental") ~ "SW620_parental",
    str_detect(sample, "HCT116_TH9616") ~ "HCT116_TH9616_r",
    str_detect(sample, "SW620_TH9619") ~ "SW620_TH9619_r",
    str_detect(sample, "HCT116_MTX") ~ "HCT116_MTX_r",
    str_detect(sample, "SW620_MTX") ~ "SW620_MTX_r",
    TRUE ~ "unknown"
  )
) %>%
  column_to_rownames("sample")

# Create SummarizedExperiment object
se <- SummarizedExperiment(
  assays = list(log2_Intensity = expr_matrix),
  rowData = row_annot,
  colData = DataFrame(sample_meta)
)

# Optional: Wrap into QFeatures object
qf <- QFeatures(list(proteins = se))





