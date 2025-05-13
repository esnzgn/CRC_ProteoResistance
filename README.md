# CRC_ProteoResistance

**CRC_ProteoResistance** is an R-based analysis project focused on proteomic changes in colorectal cancer cell lines under drug resistance conditions to the MTHFD1/2 inhibitor **TH9619** and **methotrexate (MTX)**. The study involves TMT-based quantification of six cell lines (parental and resistant lines of HCT116 and SW620).

## 📦 Dataset Summary

- Tandem Mass Tag (TMT) 18-plex quantitative proteomics
- Cell lines:
  - HCT116 (Parental, TH9619-R, MTX-R)
  - SW620 (Parental, TH9619-R, MTX-R)

## 🧪 Main Objectives

- Differential expression analysis between resistant and parental lines
- Visualization using volcano plots, heatmaps, and PCA
- Functional enrichment (GO, KEGG) of resistance-associated proteins
- Identification of unique protein markers of resistance to TH9619

## 🔬 Analysis Pipeline

- [ ] Data loading and filtering
- [ ] Differential expression analysis using limma
- [ ] Visualization:
  - [ ] Volcano plots
  - [ ] PCA plots
  - [ ] Heatmaps of top proteins
- [ ] Functional enrichment (GO/KEGG)
- [ ] Reporting significant resistance-related hits

## 📁 Folder Overview

- `data/`: Raw input Excel files
- `scripts/`: R scripts for analysis and plotting
- `results/`: Output tables and plots

## 📚 Dependencies

```r
# Install required packages
install.packages(c("tidyverse", "readxl", "pheatmap", "ggfortify"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("limma", "EnhancedVolcano", "clusterProfiler", "org.Hs.eg.db"))
