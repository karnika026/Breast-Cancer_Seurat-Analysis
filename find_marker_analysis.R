# ============================================================================
# Cluster Annotation and Visualization
# Dataset: GSE161529 (Normal Breast Tissue)
# Description: Annotate cell type clusters and create labeled UMAP visualization
# ============================================================================

# ============================================================================
# SETUP
# ============================================================================

# Load required libraries
library(Seurat)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Configuration - UPDATE THESE PATHS FOR YOUR SYSTEM
CONFIG <- list(
  paths = list(
    seurat_object = "path/to/seurat_merged_qc.rds",
    markers_file = "path/to/seurat_markers.rds",
    plot_dir = "path/to/PLOTS"
  ),
  marker_filters = list(
    p_val_threshold = 0.05,   # Adjusted p-value threshold
    logfc_threshold = 0.5,    # Log fold-change threshold
    top_n = 10                # Number of top markers per cluster
  )
)

# Create output directory if it doesn't exist
if (!dir.exists(CONFIG$paths$plot_dir)) {
  dir.create(CONFIG$paths$plot_dir, recursive = TRUE)
  message("Created output directory: ", CONFIG$paths$plot_dir)
}

# ============================================================================
# 1. LOAD DATA
# ============================================================================

message("\n=== LOADING DATA ===\n")

# Validate file paths
if (!file.exists(CONFIG$paths$seurat_object)) {
  stop("Seurat object not found: ", CONFIG$paths$seurat_object)
}

if (!file.exists(CONFIG$paths$markers_file)) {
  stop("Markers file not found: ", CONFIG$paths$markers_file)
}

# Load Seurat object and markers
seurat_merged_qc <- readRDS(CONFIG$paths$seurat_object)
seurat_markers <- readRDS(CONFIG$paths$markers_file)

message("Loaded Seurat object:")
message("  Cells: ", ncol(seurat_merged_qc))
message("  Genes: ", nrow(seurat_merged_qc))
message("  Clusters: ", length(unique(Idents(seurat_merged_qc))))

message("\nLoaded marker genes:")
message("  Total markers: ", nrow(seurat_markers))

# ============================================================================
# 2. IDENTIFY TOP MARKERS PER CLUSTER
# ============================================================================

message("\n=== IDENTIFYING TOP MARKERS ===\n")

# Filter markers by statistical significance and fold-change
top_markers <- seurat_markers %>%
  filter(
    p_val_adj < CONFIG$marker_filters$p_val_threshold,
    avg_log2FC > CONFIG$marker_filters$logfc_threshold
  ) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = CONFIG$marker_filters$top_n) %>%
  ungroup()

message("Identified top markers:")
message("  Clusters with markers: ", length(unique(top_markers$cluster)))
message("  Total top markers: ", nrow(top_markers))

# Display top markers summary
message("\nTop marker genes per cluster (first 3 per cluster):")
top3_summary <- top_markers %>%
  group_by(cluster) %>%
  slice_head(n = 3) %>%
  select(cluster, gene, avg_log2FC, p_val_adj)
print(top3_summary, n = Inf)

# Save top markers to CSV for reference
top_markers_file <- file.path(
  dirname(CONFIG$paths$plot_dir),
  "top_markers_filtered.csv"
)
write.csv(top_markers, top_markers_file, row.names = FALSE)
message("\nSaved top markers to: ", top_markers_file)

# ============================================================================
# 3. CLUSTER ANNOTATION
# ============================================================================

message("\n=== ANNOTATING CLUSTERS ===\n")

# Cell type annotations based on marker gene analysis
# These annotations are based on known marker genes for breast tissue cell types
cluster_annotation <- c(
  "0" = "Basal / myoepithelial cells",
  "1" = "Secretory luminal epithelial cells",
  "2" = "Luminal epithelial (secretory progenitor)",
  "3" = "Secretory luminal epithelial",
  "4" = "Luminal epithelial",
  "5" = "Endothelial cells",
  "6" = "Fibroblasts / stromal cells",
  "7" = "Secretory luminal epithelial",
  "8" = "Myoepithelial / basal epithelial",
  "9" = "Pericytes / vascular support",
  "10" = "Fibroblasts / ECM-remodeling",
  "11" = "Luminal epithelial / hormone-responsive",
  "12" = "Macrophages / monocytes",
  "13" = "T cells (lymphocytes)",
  "14" = "Macrophages / monocytes (activated/resident)",
  "15" = "Luminal progenitor epithelial cells",
  "16" = "Stromal / immune-interacting fibroblasts",
  "17" = "Activated fibroblasts / immune-interacting stroma",
  "18" = "Activated fibroblasts / ECM-remodeling stromal cells"
)

# Apply annotations to Seurat object
seurat_merged_qc <- RenameIdents(seurat_merged_qc, cluster_annotation)

# Store annotations in metadata for future reference
seurat_merged_qc$cell_type <- Idents(seurat_merged_qc)

message("Annotated ", length(cluster_annotation), " clusters")
message("\nCluster distribution:")
print(table(seurat_merged_qc$cell_type))

# ============================================================================
# 4. VISUALIZATION
# ============================================================================

message("\n=== CREATING VISUALIZATIONS ===\n")

# Get number of unique cell types
n_clusters <- length(unique(Idents(seurat_merged_qc)))

# Generate color palette
# Use RColorBrewer Paired palette and extend if needed
cluster_colors <- brewer.pal(n = max(3, min(12, n_clusters)), name = "Paired")

# If more colors needed, interpolate
if (n_clusters > length(cluster_colors)) {
  cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
}

message("Using ", length(cluster_colors), " colors for ", n_clusters, " cell types")

# Create annotated UMAP plot
umap_annotated <- DimPlot(
  seurat_merged_qc,
  reduction = "umap",
  label = TRUE,
  label.size = 3.5,
  repel = TRUE,
  cols = cluster_colors,
  raster = FALSE
) +
  ggtitle("Annotated Normal Breast scRNA-seq Clusters") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

# Save annotated UMAP
output_file <- file.path(CONFIG$paths$plot_dir, "annotated_breast_umap.png")
ggsave(
  filename = output_file,
  plot = umap_annotated,
  width = 14,
  height = 10,
  dpi = 300
)
message("Saved annotated UMAP: ", output_file)

# Create UMAP without labels (cleaner for publications)
umap_no_labels <- DimPlot(
  seurat_merged_qc,
  reduction = "umap",
  label = FALSE,
  cols = cluster_colors,
  raster = FALSE
) +
  ggtitle("Annotated Normal Breast scRNA-seq Clusters") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  )

output_file_no_labels <- file.path(
  CONFIG$paths$plot_dir, 
  "annotated_breast_umap_no_labels.png"
)
ggsave(
  filename = output_file_no_labels,
  plot = umap_no_labels,
  width = 14,
  height = 10,
  dpi = 300
)
message("Saved UMAP (no labels): ", output_file_no_labels)

# Create split UMAP by sample
umap_by_sample <- DimPlot(
  seurat_merged_qc,
  reduction = "umap",
  split.by = "sample",
  label = FALSE,
  cols = cluster_colors,
  raster = FALSE,
  ncol = 3
) +
  ggtitle("Cell Types Across Samples") +
  theme_minimal(base_size = 12)

output_file_sample <- file.path(
  CONFIG$paths$plot_dir, 
  "annotated_umap_by_sample.png"
)
ggsave(
  filename = output_file_sample,
  plot = umap_by_sample,
  width = 18,
  height = 12,
  dpi = 300
)
message("Saved UMAP by sample: ", output_file_sample)

# ============================================================================
# 5. SAVE ANNOTATED OBJECT
# ============================================================================

message("\n=== SAVING ANNOTATED DATA ===\n")

# Save annotated Seurat object
annotated_object_file <- file.path(
  dirname(CONFIG$paths$seurat_object),
  "seurat_merged_qc_annotated.rds"
)
saveRDS(seurat_merged_qc, file = annotated_object_file)
message("Saved annotated Seurat object: ", annotated_object_file)

# Create summary data frame
annotation_summary <- data.frame(
  cluster = names(cluster_annotation),
  cell_type = cluster_annotation,
  n_cells = as.vector(table(seurat_merged_qc$cell_type)[cluster_annotation])
)

# Add marker gene information
marker_summary <- top_markers %>%
  group_by(cluster) %>%
  summarize(
    top_markers = paste(head(gene, 5), collapse = ", "),
    .groups = "drop"
  )

annotation_summary <- annotation_summary %>%
  left_join(marker_summary, by = "cluster")

# Save annotation summary
summary_file <- file.path(
  dirname(CONFIG$paths$plot_dir),
  "cluster_annotations_summary.csv"
)
write.csv(annotation_summary, summary_file, row.names = FALSE)
message("Saved annotation summary: ", summary_file)

# ============================================================================
# ANALYSIS COMPLETE
# ============================================================================

message("\n=== ANNOTATION COMPLETE ===\n")
message("Summary:")
message("  Total cells: ", ncol(seurat_merged_qc))
message("  Cell types identified: ", n_clusters)
message("  Output files created:")
message("    - Annotated UMAP plots")
message("    - Annotated Seurat object")
message("    - Top markers CSV")
message("    - Annotation summary CSV")
message("\nAll files saved to: ", dirname(CONFIG$paths$plot_dir))
