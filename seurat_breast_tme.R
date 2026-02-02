# ============================================================================
# Single-Cell RNA-seq Analysis Pipeline
# Dataset: GSE161529 (Normal samples)
# Description: Quality control, normalization, clustering, and marker 
#              identification for scRNA-seq data using Seurat
# ============================================================================

# ============================================================================
# SETUP
# ============================================================================
library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(Matrix)
library(patchwork)
library(readr)
library(tidyverse)


# Configuration - EDIT THESE PATHS FOR YOUR SYSTEM
CONFIG <- list(
  paths = list(
    base_dir = "YOUR_PATH_HERE/Github project",
    data_dir = "YOUR_PATH_HERE/Github project/GSE161529_RAW_NORMAL",
    plot_dir = "YOUR_PATH_HERE/Github project/PLOTS",
    output_dir = "YOUR_PATH_HERE/Github project"
  ),
  qc_thresholds = list(
    min_features = 200,    # Minimum genes per cell
    max_features = 7000,   # Maximum genes per cell (doublet filter)
    max_counts = 50000,    # Maximum UMI counts (doublet filter)
    max_mito = 15          # Maximum mitochondrial percentage
  ),
  clustering = list(
    pca_dims = 15,         # Number of PCA dimensions to use
    resolution = 0.5       # Clustering resolution
  ),
  variable_features = 2000 # Number of variable features to identify
)

# Set working directory
setwd(CONFIG$paths$data_dir)

# Create output directories if they don't exist
if (!dir.exists(CONFIG$paths$plot_dir)) {
  dir.create(CONFIG$paths$plot_dir, recursive = TRUE)
}

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

# Function to save plots consistently
save_plot <- function(plot, filename, width = 20, height = 10, dpi = 300) {
  ggsave(
    filename = filename,
    plot = plot,
    path = CONFIG$paths$plot_dir,
    width = width,
    height = height,
    dpi = dpi
  )
  message("Saved plot: ", filename)
}

# Function to load a single sample
load_sample <- function(matrix_file, barcode_file, gene_names, sample_name) {
  message("  Loading matrix and barcodes...")
  
  # Read matrix and barcodes
  mat <- readMM(gzfile(matrix_file))
  barcodes <- read.delim(
    gzfile(barcode_file), 
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  # Assign row and column names
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes$V1
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = CONFIG$qc_thresholds$min_features
  )
  
  # Add sample metadata
  seurat_obj$sample <- sample_name
  
  message("  Created Seurat object: ", ncol(seurat_obj), " cells, ", 
          nrow(seurat_obj), " genes")
  
  return(seurat_obj)
}

# ============================================================================
# 1. DATA LOADING
# ============================================================================

message("\n=== DATA LOADING ===\n")

# Load feature file
features_file <- file.path(CONFIG$paths$data_dir, "features.tsv")
if (!file.exists(features_file)) {
  stop("Features file not found: ", features_file)
}

features <- read.delim(
  gzfile(features_file),
  header = FALSE,
  stringsAsFactors = FALSE
)

# Make gene names unique
gene_names <- make.unique(features$V2)
message("Loaded ", length(gene_names), " gene names")

# Find matrix and barcode files
matrix_files <- list.files(
  pattern = "matrix\\.mtx\\.gz$", 
  full.names = TRUE
)
barcode_files <- list.files(
  pattern = "barcodes\\.tsv\\.gz$", 
  full.names = TRUE
)

# Sort files to ensure matching
matrix_files <- sort(matrix_files)
barcode_files <- sort(barcode_files)

# Validate file matching
if (length(matrix_files) == 0) {
  stop("No matrix files found matching pattern 'matrix.mtx.gz'")
}

if (length(matrix_files) != length(barcode_files)) {
  stop("Mismatch: ", length(matrix_files), " matrix files but ", 
       length(barcode_files), " barcode files")
}

message("Found ", length(matrix_files), " samples to process")
message("Sample files:")
for (i in seq_along(matrix_files)) {
  message("  ", i, ". ", basename(matrix_files[i]))
}

# Load all samples into Seurat objects
seurat_list <- vector("list", length = length(matrix_files))

for (i in seq_along(matrix_files)) {
  sample_name <- sub("-matrix.*", "", basename(matrix_files[i]))
  message("\nProcessing sample ", i, " of ", length(matrix_files), 
          ": ", sample_name)
  
  seurat_list[[i]] <- load_sample(
    matrix_file = matrix_files[i],
    barcode_file = barcode_files[i],
    gene_names = gene_names,
    sample_name = sample_name
  )
}

# Name the list elements
names(seurat_list) <- sapply(seurat_list, function(x) unique(x$sample))

message("\nSuccessfully loaded ", length(seurat_list), " samples")

# ============================================================================
# 2. MERGING SEURAT OBJECTS
# ============================================================================

message("\n=== MERGING SAMPLES ===\n")

seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE161529"
)

# Validate merge
message("Merged object contains:")
message("  Total cells: ", ncol(seurat_merged))
message("  Total genes: ", nrow(seurat_merged))
message("  Samples: ", length(unique(seurat_merged$sample)))
message("\nCells per sample:")
print(table(seurat_merged$sample))

# ============================================================================
# 3. QUALITY CONTROL METRICS
# ============================================================================

message("\n=== CALCULATING QC METRICS ===\n")

# Calculate mitochondrial and ribosomal percentages
seurat_merged[["percent.mt"]] <- PercentageFeatureSet(
  seurat_merged, 
  pattern = "^MT-"
)
seurat_merged[["percent.ribo"]] <- PercentageFeatureSet(
  seurat_merged,
  pattern = "^RP[SL]"
)

message("QC metrics calculated")

# Visualize QC metrics
message("Generating QC plots...")

plots <- VlnPlot(
  seurat_merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.6,
  combine = FALSE
)

# Adjust x-axis labels
plots <- lapply(plots, function(p) {
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
})

# Save individual QC plots
save_plot(plots[[1]], "qc_nFeature_RNA.png", width = 15, height = 10)
save_plot(plots[[2]], "qc_nCount_RNA.png", width = 15, height = 10)
save_plot(plots[[3]], "qc_percent_mt.png", width = 15, height = 10)

# Feature scatter plots
plot_count_vs_mt <- FeatureScatter(
  seurat_merged, 
  feature1 = "nCount_RNA", 
  feature2 = "percent.mt", 
  raster = FALSE
)

plot_count_vs_feature <- FeatureScatter(
  seurat_merged, 
  feature1 = "nCount_RNA", 
  feature2 = "nFeature_RNA", 
  raster = FALSE
)

plot_combined <- plot_count_vs_mt + plot_count_vs_feature

save_plot(plot_count_vs_mt, "scatter_count_vs_mt.png", width = 15, height = 10)
save_plot(plot_count_vs_feature, "scatter_count_vs_feature.png", width = 15, height = 10)
save_plot(plot_combined, "scatter_combined.png", width = 20, height = 10)

# ============================================================================
# 4. FILTERING
# ============================================================================

message("\n=== FILTERING CELLS ===\n")

n_before <- ncol(seurat_merged)

seurat_merged_qc <- subset(
  seurat_merged,
  subset =
    nFeature_RNA > CONFIG$qc_thresholds$min_features &
    nFeature_RNA < CONFIG$qc_thresholds$max_features &
    nCount_RNA < CONFIG$qc_thresholds$max_counts &
    percent.mt < CONFIG$qc_thresholds$max_mito
)

n_after <- ncol(seurat_merged_qc)
n_filtered <- n_before - n_after
pct_filtered <- round((n_filtered / n_before) * 100, 1)

message("Filtering complete:")
message("  Before: ", n_before, " cells")
message("  After: ", n_after, " cells")
message("  Filtered: ", n_filtered, " cells (", pct_filtered, "%)")

# ============================================================================
# 5. NORMALIZATION
# ============================================================================

message("\n=== NORMALIZATION ===\n")

seurat_merged_qc <- NormalizeData(seurat_merged_qc)
message("Data normalized")

# ============================================================================
# 6. IDENTIFY HIGHLY VARIABLE FEATURES
# ============================================================================

message("\n=== IDENTIFYING VARIABLE FEATURES ===\n")

seurat_merged_qc <- FindVariableFeatures(
  seurat_merged_qc, 
  selection.method = "vst", 
  nfeatures = CONFIG$variable_features
)

# Get top 10 most variable genes
top10 <- head(VariableFeatures(seurat_merged_qc), 10)
message("Top 10 variable genes: ", paste(top10, collapse = ", "))

# Plot variable features with labels
variable_features_plot <- VariableFeaturePlot(seurat_merged_qc)
labeled_plot <- LabelPoints(
  plot = variable_features_plot, 
  points = top10, 
  repel = TRUE
)

save_plot(labeled_plot, "variable_features.png", width = 20, height = 10)

# ============================================================================
# 7. SCALING DATA
# ============================================================================

message("\n=== SCALING DATA ===\n")

seurat_merged_qc <- ScaleData(
  seurat_merged_qc,
  features = VariableFeatures(seurat_merged_qc)
)
message("Data scaled for ", length(VariableFeatures(seurat_merged_qc)), 
        " variable features")

# ============================================================================
# 8. DIMENSIONALITY REDUCTION (PCA)
# ============================================================================

message("\n=== RUNNING PCA ===\n")

seurat_merged_qc <- RunPCA(
  seurat_merged_qc, 
  features = VariableFeatures(object = seurat_merged_qc)
)
message("PCA complete")

# Visualize PCA
viz_dim_plot <- VizDimLoadings(
  seurat_merged_qc, 
  dims = 1:2, 
  reduction = "pca"
)
save_plot(viz_dim_plot, "pca_loadings.png", width = 20, height = 10)

dim_plot <- DimPlot(seurat_merged_qc, reduction = "pca", raster = FALSE) + 
  NoLegend()
save_plot(dim_plot, "pca_dimplot.png", width = 20, height = 10)

# Heatmap of PC genes
message("Generating PCA heatmap (this may take a while)...")
heatmap <- DimHeatmap(
  seurat_merged_qc, 
  dims = 1:15, 
  cells = 500, 
  balanced = TRUE
)

# Save heatmap to PDF
pdf(
  file = file.path(CONFIG$paths$plot_dir, "pca_heatmap.pdf"),
  width = 12,
  height = 8
)
print(heatmap)
dev.off()
message("Saved plot: pca_heatmap.pdf")

# Elbow plot for determining dimensions
elbow_plot <- ElbowPlot(seurat_merged_qc)
save_plot(elbow_plot, "pca_elbow.png", width = 20, height = 10)

# ============================================================================
# 9. CLUSTERING
# ============================================================================

message("\n=== CLUSTERING CELLS ===\n")

seurat_merged_qc <- FindNeighbors(
  seurat_merged_qc, 
  dims = 1:CONFIG$clustering$pca_dims
)

seurat_merged_qc <- FindClusters(
  seurat_merged_qc, 
  resolution = CONFIG$clustering$resolution
)

n_clusters <- length(unique(Idents(seurat_merged_qc)))
message("Identified ", n_clusters, " clusters")
message("Cluster distribution:")
print(table(Idents(seurat_merged_qc)))

# ============================================================================
# 10. UMAP VISUALIZATION
# ============================================================================

message("\n=== RUNNING UMAP ===\n")

seurat_merged_qc <- RunUMAP(
  seurat_merged_qc, 
  dims = 1:CONFIG$clustering$pca_dims
)

umap_plot <- DimPlot(
  seurat_merged_qc,
  reduction = "umap",
  raster = FALSE,
  label = TRUE,
  repel = TRUE
)
save_plot(umap_plot, "umap_clusters.png", width = 20, height = 10)



# ============================================================================
# 11. SAVE PROCESSED OBJECT
# ============================================================================

message("\n=== SAVING PROCESSED DATA ===\n")

output_file <- file.path(CONFIG$paths$output_dir, "seurat_merged_qc.rds")
saveRDS(seurat_merged_qc, file = output_file)
message("Saved Seurat object: ", output_file)

# ============================================================================
# 12. FIND MARKER GENES
# ============================================================================

message("\n=== FINDING MARKER GENES ===\n")

# Join layers if needed (for Seurat v5)
seurat_merged_qc <- JoinLayers(seurat_merged_qc)

# Find markers for all clusters
message("Finding differentially expressed genes (this may take a while)...")
seurat_markers <- FindAllMarkers(
  seurat_merged_qc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

message("Found ", nrow(seurat_markers), " marker genes across all clusters")

# Save all markers
markers_file <- file.path(CONFIG$paths$output_dir, "seurat_markers.rds")
saveRDS(seurat_markers, file = markers_file)
message("Saved marker genes: ", markers_file)

# Marker identification will be added in the next commit
# TO DO: Identify top markers genes per cluster





