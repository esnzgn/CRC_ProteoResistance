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
library(limma)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(stringr)
library(UpSetR)
library(DT)


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
load("../CRC_ProteoResistance_data/CRC_ProteoResistance_workspace.RData")

# Later, reload it using:
# load("CRC_ProteoResistance_workspace.RData")

prot <- getWithColData(qf, "norm_proteins")

cell_lines <- c("HCT116", "SW620")

de_results <- list()
deg_filtered <- list()

for (cl in cell_lines) {
  message("Processing: ", cl)
  
  idx <- colData(prot)$cell_line == cl
  se_sub <- prot[, idx]
  se_sub$group <- droplevels(factor(se_sub$group))
  
  treatments <- unique(se_sub$treatment)
  treatments <- treatments[treatments != "Parental"]
  
  for (tr in treatments) {
    comp_name <- paste0(cl, "_", tr, "_vs_Parental")
    message("  Contrast: ", comp_name)
    
    # Subset relevant samples
    selected <- se_sub[, se_sub$treatment %in% c("Parental", tr)]
    selected$treatment <- droplevels(factor(selected$treatment))
    
    # Design matrix
    design <- model.matrix(~ 0 + selected$treatment)
    colnames(design) <- levels(selected$treatment)  # e.g., "Parental", "TH9619_R"
    
    # Linear modeling
    fit <- lmFit(assay(selected), design)
    contrast_str <- paste0(tr, " - Parental")
    contrast_matrix <- makeContrasts(contrasts = contrast_str, levels = design)
    
    fit2 <- contrasts.fit(fit, contrast_matrix)
    fit2 <- eBayes(fit2)
    
    # Full result table
    top_res <- topTable(fit2, coef = 1, number = Inf) %>%
      rownames_to_column("protein_id") %>%
      mutate(comparison = comp_name)
    
    # Apply STRICT filtering
    degs <- top_res %>%
      filter(adj.P.Val < 0.0005, abs(logFC) > 3, AveExpr > 3)
    
    message("    DEGs found: ", nrow(degs))
    
    # Store results
    de_results[[comp_name]] <- top_res
    deg_filtered[[comp_name]] <- degs
  }
}

deg_counts <- sapply(deg_filtered, nrow)
deg_summary <- data.frame(
  Comparison = names(deg_counts),
  DEG_Count = deg_counts
)
print(deg_summary)

pdf("./output/plots/volcano_heat.pdf")

for (comp in names(deg_filtered)) {
  df <- de_results[[comp]]  # use full results for volcano plot
  df$Significant <- with(df, ifelse(adj.P.Val < 0.001 & abs(logFC) > 2 & AveExpr > 2, "Yes", "No"))
  
  p <- ggplot(df, aes(x = logFC, y = -log10(adj.P.Val), color = Significant)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "blue") +
    geom_hline(yintercept = -log10(0.001), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("grey60", "firebrick")) +
    labs(title = paste("Volcano Plot:", comp), x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
    theme_minimal()
  
  print(p)
}

for (comp in names(deg_filtered)) {
  top30 <- deg_filtered[[comp]] %>% 
    arrange(adj.P.Val) %>%
    slice(1:30) %>%
    pull(protein_id)
  
  mat <- assay(qf[["norm_proteins"]])[top30, , drop = FALSE]
  ann <- as.data.frame(colData(qf[["norm_proteins"]])[, c("cell_line", "treatment")])
  
  pheatmap(mat, 
           scale = "row", 
           annotation_col = ann,
           main = paste("Top 30 DEGs -", comp))
}

dev.off()

# go_results <- list()
# for (comp in names(deg_filtered)) {
#   message("Running GO BP enrichment for: ", comp)
#   
#   # Use top 100 DEGs based on adjusted P-value
#   top_genes <- deg_filtered[[comp]] %>%
#     arrange(adj.P.Val) %>%
#     head(100)
#   
#   # Extract UniProt IDs (middle part between pipes)
#   uniprot_ids <- str_extract(top_genes$protein_id, "(?<=\\|)[^|]+(?=\\|)")
#   
#   # Perform GO enrichment
#   ego <- enrichGO(
#     gene = uniprot_ids,
#     OrgDb = org.Hs.eg.db,
#     keyType = "UNIPROT",
#     ont = "BP",
#     pAdjustMethod = "BH",
#     readable = TRUE
#   )
#   
#   # Save and plot
#   go_results[[comp]] <- ego
#   
#   if (!is.null(ego) && nrow(ego) > 0) {
#     print(dotplot(ego, showCategory = 15) + ggtitle(comp))
#   } else {
#     message("No enriched GO terms found for: ", comp)
#   }
# }
# 
# 
# 
# kegg_results <- list()
# 
# for (comp in names(deg_filtered)) {
#   message("Processing KEGG for: ", comp)
#   
#   # Extract gene names
#   genes <- deg_filtered[[comp]]$protein_id
#   gene_symbols <- deg_filtered[[comp]]$gene_name  # Add this if available
#   
#   # Convert to Entrez IDs
#   entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", 
#                      toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#   
#   # Run enrichment
#   kegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
#                      organism = 'hsa', pvalueCutoff = 0.05)
#   
#   kegg_results[[comp]] <- kegg
# }



for (comp in names(deg_filtered)) {
  message("Rendering top 20 DEGs for: ", comp)
  
  dat <- deg_filtered[[comp]] %>%
    arrange(adj.P.Val) %>%
    head(20)
  
  print(DT::datatable(
    dat,
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: left;',
      paste("Top 20 DEGs for", comp)
    ),
    options = list(pageLength = 10, scrollX = TRUE),
    rownames = FALSE
  ))
}


# Create a binary matrix of DEGs
all_prots <- unique(unlist(lapply(deg_filtered, function(df) df$protein_id)))
binary_mat <- sapply(deg_filtered, function(df) all_prots %in% df$protein_id)
rownames(binary_mat) <- all_prots

# Plot overlap
upset(as.data.frame(binary_mat), nsets = length(deg_filtered), 
      order.by = "freq")


saveRDS(deg_filtered, "./output/final_deg.rds")

