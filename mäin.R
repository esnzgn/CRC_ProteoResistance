# main
rm(list = ls())

source("./reqs.R")

# Step 1: Load data
raw_data <- read_excel("../CRC_ProteoResistance_data/CPMSF_GPPF-PM-43_CPPF-PM-43_results_20240119.xlsx")

# length(unique(raw_data$`Gene Name`))

# Extract metadata
expr_cols <- grep("tmt18plex", colnames(raw_data), value = TRUE)
metadata <- tibble(
  name = expr_cols,
  condition = rep(c("HCT116_parental", "SW620_parental", "HCT116_TH9619", "SW620_TH9619", "HCT116_MTX", "SW620_MTX"), each = 3)
)

# Step 2: Create QFeatures object
expr_mat <- raw_data %>% select(all_of(expr_cols)) %>% as.matrix()
rownames(expr_mat) <- raw_data$`Gene Name`

se <- SummarizedExperiment(assays = list(counts = expr_mat))
colData(se) <- DataFrame(metadata)
qf <- QFeatures(list(proteins = se))

# Step 3: Prepare data for MSqRob2TMT
# This will extract peptide-level QFeatures and generate a structure ready for msqrob2tmt
qf_tmt <- prepareQFeatures(qf, i = 1, assayName = "peptides")


