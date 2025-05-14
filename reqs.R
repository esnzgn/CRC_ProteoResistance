# req
# install.packages(c("EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("ggfortify", "msqrob2", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
# Load libraries
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


