# main
rm(list = ls())

source("./reqs.R")

# Step 1: Load data
raw_data <- read_excel("../CRC_ProteoResistance_data/CPMSF_GPPF-PM-43_CPPF-PM-43_results_20240119.xlsx", sheet = 2)

# Basic filtering
raw_data <- raw_data %>%
  filter(`Gene Name` != "NA") %>%
  filter(as.numeric(`Amount PSMs`) >= 2) %>%
  distinct(`Gene Name`, .keep_all = TRUE)

expr_cols <- grep("tmt18plex_\\d+", names(raw_data), value = TRUE)
expr <- raw_data[, expr_cols]
row.names(expr) <- raw_data$`Gene Name`


# Step 2: Create QFeatures object
expr_mat <- raw_data %>% select(all_of(expr_cols)) %>% as.matrix()
rownames(expr_mat) <- raw_data$`Gene Name`

se <- SummarizedExperiment(assays = list(counts = expr_mat))
colData(se) <- DataFrame(metadata)
qf <- QFeatures(list(proteins = se))


# Step 3: Prepare data for MSqRob2TMT
# This will extract peptide-level QFeatures and generate a structure ready for msqrob2tmt
qf_tmt <- prepareQFeatures(qf, i = 1, assayName = "peptides")


