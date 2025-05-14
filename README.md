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

- [x] Data loading and filtering
- [ ] Contaminant Filtering (e.g., GeneID == "NA", low PSMs, zero values)
- [ ] Data Wrangling into QFeatures
- [ ] Normalization & Quality Control
- [ ] Sample Annotation & Experimental Design Setup
- [ ] Peptide-to-Protein Summarization
- [ ] Statistical Analysis (MSqRob2 or MSqRob2TMT)
- [ ] DEG Detection Across Contrasts
- [ ] Visualization (Volcano, Heatmap, PCA)
- [ ] Functional Enrichment (GO/KEGG using clusterProfiler)

## 📁 Folder Overview

- `data/`: Raw input Excel file(s) from ../ same name but _data sufixed
- `scripts/`: R scripts for analysis and plotting on the root
- `results/`: Output tables and plots

## 📚 Dependencies

in reqs.R file
