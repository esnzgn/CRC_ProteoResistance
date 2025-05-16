# req
# install.packages(c("EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("ggfortify", "msqrob2", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
# Load libraries
# Install only if not already installed
# if (!requireNamespace("pathview", quietly = TRUE)) {
#   BiocManager::install("pathview")
# }
# Required Bioconductor packages
# BiocManager::install(c("KEGGREST"))
# install.packages(c("igraph", "ggraph", "tidygraph"))
library("KEGGREST")
library("ggraph")
library("tidygraph")
library("igraph")
library(pathview)
library(QFeatures)
library(SummarizedExperiment)
library(readxl)
library(tidyverse)
library(msqrob2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(janitor)
library(dplyr)
library(ggplot2)
library(impute)
library(ggfortify)
library(reshape2)
library(pheatmap)
library(SummarizedExperiment)
library(tidyr)
library(httr)
library(jsonlite)
library(enrichplot)


