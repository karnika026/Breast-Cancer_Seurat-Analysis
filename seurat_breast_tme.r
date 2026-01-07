library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(Matrix)
library(patchwork)
library(readr)
library(tidyverse)
setwd("C:/Users/kas2071/Desktop/Github project/GSE161529_RAW_NORMAL")
features_file <- "C:/Users/kas2071/Desktop/Github project/GSE161529_RAW_NORMAL/features.tsv"
matrix_files  <- list.files(pattern = "matrix\\.mtx\\.gz$", full.names = TRUE)
barcode_files <- list.files(pattern = "barcodes\\.tsv\\.gz$", full.names = TRUE)


matrix_files <- sort(matrix_files)
barcode_files <- sort(barcode_files)
length(matrix_files)
length(barcode_files)
basename(matrix_files)[1]
basename(barcode_files)[1]
seurat_list <- list()
features <- read.delim(gzfile(features_file),
                       header = FALSE,
                       stringsAsFactors = FALSE)

gene_names <- make.unique(features$V2)

for (i in seq_along(matrix_files)) {

  mat <- readMM(gzfile(matrix_files[i]))
  barcodes <- read.delim(gzfile(barcode_files[i]),  # NOTE: features.tsv contains duplicated gene symbols.# Seurat requires unique rownames, so we use make.unique().
                         header = FALSE,
                         stringsAsFactors = FALSE)

  rownames(mat) <- gene_names
  colnames(mat) <- barcodes$V1

 sample_name <- sub("-matrix.*", "", basename(matrix_files[i]))

  seurat_obj <- CreateSeuratObject(
    counts = mat,
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )

  seurat_obj$sample <- sample_name
  seurat_list[[sample_name]] <- seurat_obj
}

length(seurat_list)
seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE161529"
)
class(seurat_merged)
seurat_merged
seurat_merged[["percent.mt"]] <- PercentageFeatureSet(seurat_merged, pattern = "^MT-")
seurat_merged[["percent.ribo"]] <- PercentageFeatureSet(seurat_merged,pattern = "^RP[SL]")
plots <- VlnPlot(
  seurat_merged,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  pt.size = 0.6,
  combine = FALSE
)
plots <- lapply(plots, function(p) {
  p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}) 
patchwork::wrap_plots(plots, ncol = 3)
print(plots[[1]])
print(plots[[2]])
print(plots[[3]])
ggsave(
  filename = "percent.mt.png",
  plot = plots[[3]],
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
plot1 <- FeatureScatter(seurat_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", raster=FALSE)
plot2 <- FeatureScatter(seurat_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster=FALSE)
plot = plot1 + plot2
print(plot1)
ggsave(filename = "plot1.png",
  plot = plot1,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(filename = "plot2.png",
  plot = plot2,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(filename = "plot.png",
  plot = plot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
seurat_merged_qc <- subset(
  seurat_merged,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 7000 &
    nCount_RNA < 50000 &
    percent.mt < 15
)
dim(seurat_merged)
dim(seurat_merged_qc)

seurat_merged_qc <- NormalizeData(seurat_merged_qc)
seurat_merged_qc <- FindVariableFeatures(seurat_merged_qc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_merged_qc), 10)

# plot variable features with and without labels
Variablefeaturesplot <- VariableFeaturePlot(seurat_merged_qc)
labelpointplot <- LabelPoints(plot = Variablefeaturesplot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)
Variableplot = Variablefeaturesplot + labelpointplot
print(Variableplot)
ggsave(filename = "Variableplot.png",
  plot = Variableplot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)

seurat_merged_qc <- ScaleData(seurat_merged_qc,features = VariableFeatures(seurat_merged_qc))
str(seurat_merged_qc)
class(seurat_merged_qc)
seurat_merged_qc
head(seurat_merged_qc@meta.data)
tail(seurat_merged_qc@meta.data)

#Run PCA
seurat_merged_qc <- RunPCA(seurat_merged_qc, features = VariableFeatures(object = seurat_merged_qc))
Vizdimplot <-VizDimLoadings(seurat_merged_qc, dims = 1:2, reduction = "pca")
ggsave(filename = "Vizdimplot.png",
  plot = Vizdimplot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
dimplot <- DimPlot(seurat_merged_qc, reduction = "pca", raster = FALSE) + NoLegend()
ggsave(filename = "dimplot.png",
  plot = dimplot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
heatmap <- DimHeatmap(seurat_merged_qc, dims = 1:15, cells = 500, balanced = TRUE)
pdf(
  file = "C:/Users/kas2071/Desktop/Github project/PLOTS/heatmap.pdf",
  width = 12,
  height = 8
)
dev.off()
#Dimensional heatmaps were saved using base R graphics devices due to ComplexHeatmap-based rendering in Seurat.

elbowplot <- ElbowPlot(seurat_merged_qc)
ggsave(filename = "elbowplot.png",
  plot = elbowplot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
seurat_merged_qc <- FindNeighbors(seurat_merged_qc, dims = 1:15)
seurat_merged_qc <- FindClusters(seurat_merged_qc, resolution = 0.5)
head(Idents(seurat_merged_qc), 5)
seurat_merged_qc <- RunUMAP(seurat_merged_qc, dims = 1:15)
umap_plot <- DimPlot(
  seurat_merged_qc,
  reduction = "umap",
  raster = FALSE,
  label = TRUE,
  repel = TRUE
)
ggsave(filename = "umap_plot.png",
  plot = umap_plot,
  path = "C:/Users/kas2071/Desktop/Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)

saveRDS(
  seurat_merged_qc,
  file = "C:/Users/kas2071/Desktop/Github project/seurat_merged_qc.rds"
)
seurat_merged_qc <- readRDS(
  "C:/Users/kas2071/Desktop/Github project/seurat_merged_qc.rds"
)
class(seurat_merged_qc)
Idents(seurat_merged_qc)
seurat_merged_qc <- JoinLayers(seurat_merged_qc)
seurat_markers <- FindAllMarkers(
  seurat_merged_qc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

DefaultAssay(seurat_merged_qc)
Assays(seurat_merged_qc)
class(seurat_markers)
colnames(seurat_markers)
head(seurat_markers)
tail(seurat_markers)
saveRDS(
  seurat_markers,
  file = "C:/Users/kas2071/Desktop/Github project/seurat_markers.rds"
)
