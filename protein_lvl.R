rm(list = ls())

# Load packages
library(QFeatures)
library(SummarizedExperiment)
library(readxl)
library(tidyverse)
library(msqrob2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(ggfortify)
library(pheatmap)

# Load data
excel_path <- "../CRC_ProteoResistance_data/CPMSF_GPPF-PM-43_CPPF-PM-43_results_20240119.xlsx"
quan_df <- read_excel(excel_path, sheet = "Quan_data")

# Standardize first 3 columns
colnames(quan_df)[1:3] <- c("protein_id", "gene_name", "description")

# Fix TH9616 → TH9619 in column names
fixed_colnames <- colnames(quan_df)
fixed_colnames <- str_replace_all(fixed_colnames, "TH9616", "TH9619")
colnames(quan_df) <- fixed_colnames

# Convert intensity columns to numeric and fix 0s or NAs
for (i in 4:ncol(quan_df)) {
  quan_df[[i]] <- suppressWarnings(as.numeric(quan_df[[i]]))
  quan_df[[i]][is.na(quan_df[[i]]) | quan_df[[i]] <= 0] <- 1e-3
}

# Build assay matrix and row metadata
assay_data <- as.matrix(quan_df[, 4:ncol(quan_df)])
rownames(assay_data) <- make.unique(quan_df$protein_id)
row_metadata <- quan_df[, 1:3]
rownames(row_metadata) <- rownames(assay_data)

# Construct colData metadata (cell line, treatment)
sample_names <- colnames(assay_data)
cell_line <- case_when(
  str_detect(sample_names, "HCT116") ~ "HCT116",
  str_detect(sample_names, "SW620") ~ "SW620",
  TRUE ~ NA_character_
)
treatment <- case_when(
  str_detect(sample_names, "parental") ~ "Parental",
  str_detect(sample_names, "TH9619") ~ "TH9619_R",
  str_detect(sample_names, "MTX") ~ "MTX_R",
  TRUE ~ "Unknown"
)
group <- paste(cell_line, treatment, sep = ".")

col_metadata <- DataFrame(sample = sample_names,
                          cell_line = cell_line,
                          treatment = treatment,
                          group = group,
                          row.names = sample_names)

# Construct SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(intensity = assay_data),
  rowData = row_metadata,
  colData = col_metadata
)

# Wrap in QFeatures
qf <- QFeatures(list(proteins = se))

# Transform and normalize
qf <- logTransform(qf, i = "proteins", name = "log_proteins")
qf <- normalize(qf, i = "log_proteins", name = "norm_proteins", method = "center.median")

# Filter out rows with too many missing values
qf <- filterNA(qf, i = "norm_proteins", pNA = 0.5)

# Preview summary
print(qf)
print(table(colData(qf[["norm_proteins"]])$group))

# save env ...
save.image(file = "../CRC_ProteoResistance_data/CRC_ProteoResistance_workspace.RData")

# Later, reload it using:
# load("CRC_ProteoResistance_workspace.RData")

