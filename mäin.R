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

# QC ####
# Count missing per sample
# Extract the expression matrix correctly
expr_matrix <- assay(qf[["proteins"]], "log2_Intensity")

# Now compute missing values per sample (column-wise)
missing_counts <- colSums(is.na(expr_matrix))

# View missing data
print(missing_counts)

# Optional: visualize
# Prepare the data frame for plotting
missing_df <- data.frame(
  Sample = names(missing_counts),
  MissingValues = missing_counts
)

# Plot using ggplot
ggplot(missing_df, aes(x = reorder(Sample, MissingValues), y = MissingValues)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "Missing values per sample",
       x = "Sample",
       y = "Count of missing values") +
  theme(axis.text.y = element_text(size = 8))

protein_missing_counts <- rowSums(is.na(expr_matrix))
summary(protein_missing_counts)
table(protein_missing_counts)

ggplot(data.frame(MissingCount = protein_missing_counts), aes(x = MissingCount)) +
  geom_histogram(binwidth = 1, fill = "steelblue") +
  theme_minimal() +
  labs(title = "Missing value distribution across proteins",
       x = "Number of missing values",
       y = "Number of proteins")

threshold <- ceiling(ncol(expr_matrix) * 2 / 3)
expr_matrix_filtered <- expr_matrix[rowSums(!is.na(expr_matrix)) >= threshold, ]
dim(expr_matrix_filtered)  # Check how many proteins remain

expr_matrix_imputed <- impute.knn(expr_matrix_filtered)$data

prcomp_res <- prcomp(t(expr_matrix_imputed), scale. = TRUE)

autoplot(prcomp_res, data = as.data.frame(colData(qf[["proteins"]])), colour = "condition") +
  theme_minimal() +
  labs(title = "PCA of Samples", x = "PC1", y = "PC2")

# qf checks
# Check number of assays (quantification levels)
length(qf)

# Check names of each assay
names(qf)

# Check the class of each assay
sapply(qf, class)

# Or see structure
str(qf, max.level = 1)

# checks the excell file again
names(all_sheets)


# EDA ####
expr_long <- melt(expr_matrix_imputed)
colnames(expr_long) <- c("Protein", "Sample", "Log2Intensity")

ggplot(expr_long, aes(x = Sample, y = Log2Intensity)) +
  geom_boxplot() +
  theme_minimal(base_size = 10) +
  coord_flip() +
  labs(title = "Log2 Intensities per Sample")

sample_cor <- cor(expr_matrix_imputed, method = "pearson")
pheatmap(sample_cor, main = "Sample Correlation Heatmap")

# dep ####
# Make sure protein IDs exist in the annotation
head(row_annot$protein_id)

# Set rownames of both to protein_id
rownames(expr_matrix) <- row_annot$protein_id
rownames(row_annot) <- row_annot$protein_id

library(SummarizedExperiment)

se <- SummarizedExperiment(
  assays = list(log2_Intensity = expr_matrix),
  rowData = row_annot,
  colData = DataFrame(sample_meta)
)

validObject(se)

unique(colData(se)$condition)
table(colData(se)$condition)


colData = DataFrame(sample_meta)
rownames(sample_meta)
colnames(se)

sample_meta$condition <- c(
  "hct116_parental", "sw620_parental", "hct116_parental",
  "sw620_parental", "hct116_parental", "sw620_parental",
  "hct116_th9616_r", "sw620_th9619_r", "hct116_th9616_r",
  "sw620_th9619_r", "hct116_th9616_r", "sw620_th9619_r",
  "hct116_mtx_r", "sw620_mtx_r", "hct116_mtx_r",
  "sw620_mtx_r", "hct116_mtx_r", "sw620_mtx_r"
)

colData(se)$condition <- sample_meta$condition

table(colData(se)$condition)


# Subset to HCT116 samples only
se_hct116 <- se[, grepl("^hct116", colData(se)$condition)]

# Check sample conditions
table(colData(se_hct116)$condition)

se_hct116 <- normalizeD(se_hct116, i = "log2_Intensity", name = "norm")

# Manually normalize by median-centering (column-wise)
assay(se_hct116, "norm") <- sweep(
  assay(se_hct116, "log2_Intensity"),
  2,
  apply(assay(se_hct116, "log2_Intensity"), 2, median, na.rm = TRUE),
  FUN = "-"
)

# Confirm assay was added
assays(se_hct116)

se_hct116 <- aggregateFeatures(
  se_hct116,
  i = "norm",
  fcol = "protein_id", # or another suitable feature column
  name = "protein",
  na.rm = TRUE
)

se_hct116 <- msqrob2::msqrob(object = se_hct116, formula = ~condition, overwrite = TRUE)


msqrob2::getCoef(rowData(se_hct116)[["protein"]]$msqrobModels[[1]])


objs <- mget(ls("package:ggplot2"), envir = as.environment("package:ggplot2"))
funcs <- objs[sapply(objs, is.function)]
names(funcs)


# Check structure of rowData
head(rowData(se_hct116))

# Access one model correctly
model1 <- rowData(se_hct116)$msqrobModels[[1]]

# Check if model is not NULL
if (!is.null(model1)) {
  coef1 <- msqrob2::getCoef(model1)
  print(coef1)
}


coefs <- lapply(rowData(se_hct116)$msqrobModels, function(model) {
  if (!is.null(model)) msqrob2::getCoef(model) else NA
})
names(coefs) <- rownames(se_hct116)

coef_df <- do.call(rbind, lapply(seq_along(coefs), function(i) {
  coef_row <- coefs[[i]]
  if (!is.null(coef_row)) {
    data.frame(protein = rownames(se_hct116)[i], t(coef_row))
  } else {
    data.frame(protein = rownames(se_hct116)[i], NA)
  }
}))


coef_df <- lapply(seq_along(coefs), function(i) {
  coef_row <- coefs[[i]]
  if (!is.null(coef_row)) {
    # Convert to data frame and add protein name
    as.data.frame(t(coef_row)) %>%
      mutate(protein = rownames(se_hct116)[i])
  } else {
    # Return NA row with protein name
    data.frame(protein = rownames(se_hct116)[i], NA)
  }
}) %>%
  bind_rows()

table(sapply(coefs, is.null))

which(sapply(coefs, is.null))[1:10]

# Get all unique coefficient names
all_coef_names <- unique(unlist(lapply(coefs, function(x) if (!is.null(x)) names(x))))

# Build data frame safely
coef_df <- lapply(seq_along(coefs), function(i) {
  coef_row <- coefs[[i]]
  if (!is.null(coef_row) && length(coef_row) == length(all_coef_names)) {
    out <- as.data.frame(t(coef_row))
  } else if (!is.null(coef_row)) {
    # If coef_row is not null but the length is wrong, fill with NAs and give a warning
    warning(sprintf("Protein %s has coefficient length %d ≠ expected %d. Filling with NAs.",
                    rownames(se_hct116)[i], length(coef_row), length(all_coef_names)))
    out <- as.data.frame(setNames(as.list(rep(NA, length(all_coef_names))), all_coef_names))
  } else {
    out <- as.data.frame(setNames(as.list(rep(NA, length(all_coef_names))), all_coef_names))
  }
  out$protein <- rownames(se_hct116)[i]
  out
}) %>% bind_rows()

if (!is.null(coef_row) && length(coef_row) == length(all_coef_names)) {
  coef_row <- setNames(coef_row, all_coef_names)
  out <- as.data.frame(t(coef_row))
}

coef_df <- lapply(seq_along(coefs), function(i) {
  coef_row <- coefs[[i]]
  
  if (!is.null(coef_row) && length(coef_row) == length(all_coef_names)) {
    coef_row <- setNames(coef_row, all_coef_names)
    out <- as.data.frame(t(coef_row))
  } else {
    if (!is.null(coef_row)) {
      warning(sprintf("Protein %s has coefficient length %d ≠ expected %d. Filling with NAs.",
                      rownames(se_hct116)[i], length(coef_row), length(all_coef_names)))
    }
    out <- as.data.frame(setNames(as.list(rep(NA, length(all_coef_names))), all_coef_names))
  }
  
  out$protein <- rownames(se_hct116)[i]
  out
}) %>% bind_rows()

# Optional cleanup of accidental columns like "V1"
coef_df <- coef_df[, !(names(coef_df) %in% c("V1", "X.Intercept."))]

sum(complete.cases(coef_df[, all_coef_names]))

coef_df %>%
  arrange(desc(abs(conditionhct116_th9616_r))) %>%
  head(10)

coef_df_long <- coef_df %>%
  tidyr::pivot_longer(cols = all_coef_names, names_to = "term", values_to = "estimate")

ggplot(coef_df_long, aes(x = estimate)) +
  geom_histogram(bins = 50) +
  facet_wrap(~ term, scales = "free_x") +
  theme_minimal() +
  labs(title = "Distribution of model coefficients", x = "Estimate", y = "Count")

L <- makeContrast(
  "conditionhct116_th9616_r - conditionhct116_parental = 0",
  parameterNames = c("conditionhct116_parental", "conditionhct116_th9616_r")
)
se_hct116 <- hypothesisTest(se_hct116, contrast = L)

# Extract name of contrast you used
colnames(L)
# [1] "conditionhct116_th9616_r - conditionhct116_parental"

# Use it to get the DE results:
res <- rowData(se_hct116)[["protein"]][[colnames(L)]]

# Check rowData structure
head(rowData(se_hct116))

# Extract result for TH9616-R vs Parental
res <- rowData(se_hct116)[["conditionhct116_th9616_r - conditionhct116_parental"]]
str(res)  # should now show a data.frame with logFC, pval, adjPval

ggplot(res, aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(
    title = "TH9616-R vs Parental",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  )

L_mtx <- makeContrast(
  "conditionhct116_mtx_r - conditionhct116_parental = 0",
  parameterNames = c("conditionhct116_parental", "conditionhct116_mtx_r")
)
se_hct116 <- hypothesisTest(se_hct116, contrast = L_mtx, overwrite = TRUE)

se_hct116 <- hypothesisTest(
  object = se_hct116,
  contrast = L_mtx,
  overwrite = TRUE,
  resultsColumnNamePrefix = "MTX_vs_Parental"
)

res_mtx <- rowData(se_hct116)[["MTX_vs_Parentalconditionhct116_mtx_r"]]

summary(res_mtx$logFC)
summary(res_mtx$pval)

ggplot(res_mtx, aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(alpha = 0.7) +
  theme_minimal() +
  scale_color_manual(values = c("black", "red")) +
  labs(
    title = "MTX-R vs Parental",
    x = "log2 Fold Change",
    y = "-log10(p-value)"
  )


# Filter and preserve rownames as a column
top_hits <- res_mtx %>%
  dplyr::filter(adjPval < 0.05) %>%
  dplyr::arrange(desc(abs(logFC))) %>%
  dplyr::mutate(protein = rownames(.) )

head(top_hits, 10)


# Extract UniProt Accession
top_hits <- top_hits %>%
  mutate(
    accession = sub(".*\\|(.*?)\\|.*", "\\1", protein),
    gene_name = sub(".*\\|.*\\|(.*?)_HUMAN", "\\1", protein)
  )

fetch_uniprot_info <- function(accession) {
  url <- paste0("https://rest.uniprot.org/uniprotkb/", accession, ".json")
  res <- GET(url)
  
  if (status_code(res) == 200) {
    info <- fromJSON(content(res, as = "text", encoding = "UTF-8"))
    
    # Safely extract full name (check if recommendedName exists)
    full_name <- tryCatch({
      if (!is.null(info$proteinDescription$recommendedName$fullName$value)) {
        info$proteinDescription$recommendedName$fullName$value
      } else if (!is.null(info$proteinDescription$submissionNames[[1]]$fullName$value)) {
        info$proteinDescription$submissionNames[[1]]$fullName$value
      } else {
        NA
      }
    }, error = function(e) NA)
    
    # Safely extract gene name
    gene <- tryCatch({
      if (!is.null(info$genes[[1]]$geneName$value)) {
        info$genes[[1]]$geneName$value
      } else {
        NA
      }
    }, error = function(e) NA)
    
    data.frame(
      accession = accession,
      full_name = full_name,
      gene = gene,
      stringsAsFactors = FALSE
    )
  } else {
    message("Failed to fetch for: ", accession)
    data.frame(accession = accession, full_name = NA, gene = NA)
  }
}

# Annotate a subset (e.g. top 10) to avoid throttling
annotations <- lapply(head(top_hits$accession, 200), fetch_uniprot_info) %>% bind_rows()
print(annotations)

# Merge with top hits
google  <- left_join(top_hits, annotations, by = "accession")

head(google)

write.csv(google, "MTX_vs_Parental_annotated_hits.csv", row.names = FALSE)

ggplot(google, aes(x = logFC, y = -log10(pval), color = adjPval < 0.05)) +
  geom_point(alpha = 0.7) +
  geom_text(data = subset(google, adjPval < 0.05 & abs(logFC) > 1.5),
            aes(label = gene_name), vjust = 1.5, hjust = 0.5, size = 3) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(title = "MTX-R vs Parental (Annotated Top Hits)",
       x = "log2 Fold Change", y = "-log10(p-value)")

sig_genes <- na.omit(unique(google$gene_name[google$adjPval < 0.05]))
write.table(sig_genes, "sig_genes_mtx_vs_parental.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#enrichment ####
# Your significant gene symbols (top hits)
gene_symbols <- top_hits$gene_name

# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# GO enrichment (Biological Process)
ego_bp <- enrichGO(gene = entrez_ids$ENTREZID,
                   OrgDb = org.Hs.eg.db,
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE)

# KEGG pathway enrichment
ekegg <- enrichKEGG(gene = entrez_ids$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)

# Plot top GO BP terms
dotplot(ego_bp, showCategory = 10, title = "GO Biological Process Enrichment") + theme_minimal()

# Plot top KEGG pathways
dotplot(ekegg, showCategory = 10, title = "KEGG Pathway Enrichment") + theme_minimal()

# View results
head(as.data.frame(ego_bp))
head(as.data.frame(ekegg))

# Assuming you have this mapping already from top_hits:
# top_hits$gene_name and top_hits$logFC

# Convert gene names to Entrez IDs
gene_symbols <- top_hits$gene_name
logfc_vector <- top_hits$logFC
names(logfc_vector) <- gene_symbols

# Map SYMBOL to ENTREZID
entrez_map <- bitr(names(logfc_vector), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_logfc <- merge(entrez_map, data.frame(SYMBOL = names(logfc_vector), logFC = logfc_vector), by = "SYMBOL")

# Make a named vector for pathview
gene_data <- entrez_logfc$logFC
names(gene_data) <- entrez_logfc$ENTREZID

# Example: pathway in cancer (hsa05200)
pathview(
  gene.data  = gene_data,
  pathway.id = "hsa05200",  # Use any valid KEGG pathway ID (see below for tips)
  species    = "hsa",
  gene.idtype = "entrez",
  limit      = list(gene = max(abs(gene_data))),  # auto scale
  out.suffix = "MTX_resistance"
)

head(ekegg@result[, c("ID", "Description")], 10)


# endhance pathview outp
# Your gene data
gene_list <- setNames(top_hits$logFC, top_hits$accession)

# Convert to Entrez IDs
entrez <- bitr(names(gene_list), fromType="UNIPROT", toType="ENTREZID", OrgDb=org.Hs.eg.db)
gene_entrez <- setNames(gene_list[entrez$UNIPROT], entrez$ENTREZID)

# Run pathview only with mapped genes (suppress others)
pathview(gene.data = gene_entrez,
         pathway.id = "hsa05200",
         species = "hsa",
         limit = list(gene=max(abs(gene_entrez))),
         gene.idtype = "entrez",
         kegg.native = TRUE,
         out.suffix = "MTX_resistance",
         low = list(gene="green"),
         mid = list(gene="gray"),
         high = list(gene="red"),
         node.sum = "max")

# Download KGML file
kegg_url <- "https://rest.kegg.jp/get/hsa05200/kgml"
download.file(kegg_url, "hsa05200.xml")

# Parse and convert to graph
kegg_graph <- parseKGML("hsa05200.xml")
gR <- KEGGpathway2Graph(kegg_graph, genesOnly = TRUE)

# Subset graph to only your genes
g_igraph <- igraph.from.graphNEL(gR)
nodes_of_interest <- V(g_igraph)$name[V(g_igraph)$name %in% entrez_ids]
sub_g <- induced_subgraph(g_igraph, vids = nodes_of_interest)

# Plot reduced network
plot(sub_g,
     vertex.label = V(sub_g)$name,
     main = "Top Genes in Pathway hsa05200",
     vertex.color = "lightblue", edge.arrow.size = 0.5)

# Load and parse pathway
kgml_file <- "hsa05200.xml"
pathway <- parseKGML(kgml_file)
graph <- KEGGpathway2Graph(pathway, genesOnly = TRUE)

# Check node names
graph::nodes(graph)[1:10]

# Strip the "hsa:" prefix from the beginning of the string
gene_string <- sub("^hsa:", "", genes_kegg[1])

# Now safely evaluate the string into an actual character vector
real_genes <- eval(parse(text = gene_string))

# Prefix the gene symbols with "hsa:" to match node names in the graph
genes_kegg_cleaned <- paste0("hsa:", real_genes)

# Convert gene symbols to ENTREZ IDs
entrez_ids <- mapIds(
  org.Hs.eg.db,
  keys = real_genes,
  column = "ENTREZID",
  keytype = "SYMBOL",
  multiVals = "first"
)

# Drop NAs (genes that couldn't be mapped)
entrez_ids <- na.omit(entrez_ids)

# Add "hsa:" prefix to match graph format
genes_kegg_matched <- paste0("hsa:", entrez_ids)

# Now find matching nodes
nodes_of_interest <- V(g_igraph)$name[V(g_igraph)$name %in% genes_kegg_matched]

# Check how many matched
length(nodes_of_interest)

real_genes[is.na(entrez_ids)]

subgraph_interest <- induced_subgraph(g_igraph, vids = nodes_of_interest)

g_tbl <- as_tbl_graph(subgraph_interest)

ggraph(g_tbl, layout = "kk") +
  geom_edge_link(alpha = 0.3) +
  geom_node_point(color = "steelblue", size = 3) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.5) +
  theme_void()

V(g_igraph)$matched <- V(g_igraph)$name %in% nodes_of_interest

plot(
  g_igraph,
  vertex.color = ifelse(V(g_igraph)$matched, "red", "gray"),
  vertex.size = 3,
  vertex.label = NA,
  edge.arrow.size = 0.3
)

# 1. Extract vertex names (Entrez-style KEGG IDs)
entrez_kegg_ids <- V(g_igraph)$name

# 2. Remove "hsa:" prefix
entrez_ids_clean <- gsub("hsa:", "", entrez_kegg_ids)

# 3. Map back to gene symbols using org.Hs.eg.db
symbol_map <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids_clean,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# 4. Replace NA with original IDs (fallback)
symbol_labels <- ifelse(is.na(symbol_map), entrez_kegg_ids, symbol_map)

# 5. Assign as vertex labels
V(g_igraph)$label <- symbol_labels

plot(g_igraph,
     vertex.size = 5,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans")

# Create a label vector: only assign gene symbols to nodes of interest
V(g_igraph)$label <- ifelse(
  V(g_igraph)$name %in% nodes_of_interest,
  symbol_map[gsub("hsa:", "", V(g_igraph)$name)],
  NA  # Leave other labels blank
)

# Optional: tweak the appearance for better clarity
plot(g_igraph,
     vertex.size = 5,
     vertex.label.cex = 0.7,
     vertex.label.color = "black",
     vertex.label.family = "sans",
     vertex.label.dist = 0.5,
     vertex.color = ifelse(V(g_igraph)$name %in% nodes_of_interest, "orange", "gray"),
     edge.color = "gray")

V(g_igraph)$size <- ifelse(V(g_igraph)$name %in% nodes_of_interest, 8, 3)
V(g_igraph)$frame.color <- ifelse(V(g_igraph)$name %in% nodes_of_interest, "black", NA)

layout <- layout_with_fr(g_igraph)
plot(g_igraph, layout = layout, ...)



# Confirm assay was added values of those not expressingshalonf
assays(se_hct116)

se_hct116 <- aggregateFeatures(
  se_hct116,
  i = "norm",
  fcol = "protein_id", # or another suitable feature column
  name = "protein",
  na.rm = TRUE
)
