# Breast-Cancer_Seurat-Analysis
Single-cell RNA-seq analysis of breast cancer data using Seurat.
This project performs single-cell RNA-seq analysis of **normal breast tissue samples** from the GSE161529 dataset using **Seurat**.  
The primary objective is to identify cell populations, perform clustering, and determine marker genes.

## Project Objective
- Perform a comprehensive single-cell RNA sequencing analysis of normal human breast tissue using the GSE161529 dataset. 
- To identify transcriptionally distinct cell populations.
- Annotate major breast cell types based on known marker genes.
- Explore cellular heterogeneity within normal breast tissue.
- Generate a high-quality annotated reference dataset representing the cellular composition of normal breast tissue.


# Dataset
- **GEO Accession:** GSE161529  
- **Tissue:** Normal breast tissue  
- **Organism:** Homo sapiens  
- **Platform**:	GPL18573	Illumina NextSeq 500 (Homo sapiens)
- **Data Type** : Single-cell RNA seq

# Analysis Workflow
  - Data Download and Pre-Processing
  - Quality Control and filtering
  - Normalization and Scaling
  - Dimensionality Reduction (PCA, UMAP)
  - Clustering and Marker Identification
  - Visualization and annotation of Result
  
# Tools and Environment
- Operating Environment: Visual Studio Code
- Programming Language: R
- R Version: 4.5.2
- Seurat Version: 5.2.0
- Additional Packages: ggplot2, dplyr, patchwork, tidyverse, Matrix.    
  
