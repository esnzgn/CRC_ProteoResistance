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

## 🔬 Analysis Pipeline and rep

- [x] Data loading
-   [x] Data cleaning
-   [x] extracted and converted the expression matrix
-   [x] Built the sample metadata (colData)
-   [x] Built the feature metadata (rowData)
-   [x] Created a SummarizedExperiment object
-   [x] Wrapped everything in a QFeatures object for downstream quantitative proteomics analysis 🎉
-   [x] Preprocessing and Quality Control (QC) Summary
-   [x] Missing Value Assessment:
  -   [x] All 18 samples had exactly 8 missing values, which is low and uniform.
-   [x] Over 10,771 proteins had no missing values.
-   [x] Only 8 proteins had 18 missing values (likely dropped after filtering).
-   [x] Filtering and Imputation:
  -   [x] You applied a 2/3 presence threshold: proteins retained if observed in ≥12 samples.
-   [x] Resulted in a matrix of 10,771 proteins × 18 samples.
-   [x] Used impute.knn() from the impute package to handle remaining NAs.
-   [x] PCA & EDA:
  -   [x] PCA plot shows clear sample separation, though the condition label needs refinement (currently "unknown").
-   [x] Boxplots of log2 intensities per sample show consistent distributions.
-   [x] Sample correlation heatmap reveals clear clustering, indicating batch or condition-related grouping. 
-   [x] 10779 rows: proteins (with names like "sp|A0A0B4J2F0|PIOS1_HUMAN")
-   [x] 18 columns: samples (with names like "hct116_parental_..._tmt18plex_126")
-   [x] rowData: includes protein_id, gene_name, description
-   [x] colData: includes condition for downstream modeling
-   [x] alidObject(se) returns TRUE → no internal issues
- [x] Peptide-to-Protein Summarization
- [x] Statistical Analysis (MSqRob2 or MSqRob2TMT)
- [x] DEG Detection Across Contrasts
- [x] Visualization (Volcano, Heatmap, PCA)
- [x] Functional Enrichment (GO/KEGG using clusterProfiler)

## 📁 Folder Overview

- `data/`: Raw input Excel file(s) from ../ same name but _data sufixed
- `scripts/`: R scripts for analysis and plotting on the root
- `results/`: Output tables and plots

## 📚 Dependencies

in reqs.R file
