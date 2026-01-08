library(Seurat)
library(SeuratObject)
library(patchwork)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

seurat_markers <- readRDS("C:/Users/kas2071/Desktop/Github project/seurat_markers.rds")
head(seurat_markers)
tail(seurat_markers)
# Identifying top markers per cluster
top_markers <- seurat_markers %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)
head(top_markers)
tail(top_markers)
top_markers

 
 #Cluster Annotation
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
names(cluster_annotation) <- 0:18
seurat_merged_qc <- RenameIdents(seurat_merged_qc, cluster_annotation)
n_clusters <- length(unique(Idents(seurat_merged_qc)))
cluster_colors <- brewer.pal(n = max(3, n_clusters), name = "Paired") # "Paired" gives visually distinct colors
if(n_clusters > length(cluster_colors)) {
  cluster_colors <- colorRampPalette(cluster_colors)(n_clusters)
}
umap_plot <- DimPlot(
  seurat_merged_qc,
  reduction = "umap",
  label = TRUE,
  label.size = 4,   # adjust text size
  repel = TRUE,     # avoids overlapping labels
  cols = cluster_colors,
  raster = FALSE
) +
  ggtitle("Annotated Normal Breast scRNA-seq Clusters") +
  theme_minimal(base_size = 14)



ggsave(
  filename = "Annotated_Breast_UMAP.png",
  plot = umap_plot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 12,
  height = 10,
  dpi = 300
)
