# ------------------------- Setup -------------------------
rm(list = ls())
source("./reqs.R")

# ------------------------- Load Data -------------------------
excel_path <- "../CRC_ProteoResistance_data/CPMSF_GPPF-PM-43_CPPF-PM-43_results_20240119.xlsx"
sheet_names <- readxl::excel_sheets(excel_path)
all_sheets <- set_names(sheet_names) %>% 
  map(~ read_excel(excel_path, sheet = .x) %>% janitor::clean_names())

expr_data <- all_sheets[["Quan_data"]]

# Remove duplicated or NA protein IDs BEFORE creating matrix
expr_data <- expr_data %>%
  filter(!is.na(protein_id)) %>%
  distinct(protein_id, .keep_all = TRUE)

# Extract matrix
expr_matrix <- expr_data %>%
  select(-gene_name, -description) %>%
  column_to_rownames("protein_id") %>%
  as.matrix() %>%
  apply(2, as.numeric)

# Row annotation
row_annot <- expr_data %>% select(protein_id, gene_name, description)

# Sample metadata
sample_names <- colnames(expr_matrix)
sample_meta <- tibble(
  sample = sample_names,
  condition = case_when(
    str_detect(sample, "HCT116_parental") ~ "hct116_parental",
    str_detect(sample, "TH9616") ~ "hct116_th9616_r",
    str_detect(sample, "MTX") ~ "hct116_mtx_r",
    TRUE ~ "unknown"
  )
) %>% column_to_rownames("sample")

# ------------------------- SummarizedExperiment + QFeatures -------------------------
se <- SummarizedExperiment(
  assays = list(log2_Intensity = expr_matrix),
  rowData = row_annot,
  colData = DataFrame(sample_meta)
)

# Subset for HCT116 only
se_hct116 <- se[, grepl("^hct116", colData(se)$condition)]

# Create QFeatures object
qf <- QFeatures(list(proteins = se_hct116))

# ------------------------- Normalize (Center Median) -------------------------
qf <- normalize(qf, i = "proteins", name = "norm", method = "center.median")
