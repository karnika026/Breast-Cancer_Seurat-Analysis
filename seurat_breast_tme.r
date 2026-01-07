library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(Matrix)
library(patchwork)
library(readr)
library(tidyverse)
setwd("....../Github project/GSE161529_RAW_NORMAL")

#Loading data
features_file <- "......./Github project/GSE161529_RAW_NORMAL/features.tsv"
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
  barcodes <- read.delim(gzfile(barcode_files[i]), header = FALSE,
                         stringsAsFactors = FALSE)
  rownames(mat) <- gene_names
  colnames(mat) <- barcodes$V1
sample_name <- sub("-matrix.*", "", basename(matrix_files[i]))

#Creating Seurat Object
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


#Merging Seurat Object
seurat_merged <- merge(
  seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  project = "GSE161529"
)
class(seurat_merged)


#Quality control 
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
  path = "...../Github project/PLOTS",
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
  path = "......./Github project/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(filename = "plot2.png",
  plot = plot2,
  path = "....../Github project/PLOTS",
  width = 15,
  height = 10,
  dpi = 300
)
ggsave(filename = "plot.png",
  plot = plot,
  path = "....../Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)


# Quality Control and Filtering

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
  path = "....../Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)

#Scaling Data
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
  path = "......./Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
dimplot <- DimPlot(seurat_merged_qc, reduction = "pca", raster = FALSE) + NoLegend()
ggsave(filename = "dimplot.png",
  plot = dimplot,
  path = "......./Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)
heatmap <- DimHeatmap(seurat_merged_qc, dims = 1:15, cells = 500, balanced = TRUE)
pdf(file = "......../Github project/PLOTS/heatmap.pdf",width = 12,height = 8)
dev.off()

#Dimensional heatmaps were saved using base R graphics devices due to ComplexHeatmap-based rendering in Seurat.


#Elbow Plot
elbowplot <- ElbowPlot(seurat_merged_qc)
ggsave(filename = "elbowplot.png",
  plot = elbowplot,
  path = "......./Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)

#Clustering the cells
seurat_merged_qc <- FindNeighbors(seurat_merged_qc, dims = 1:15)
seurat_merged_qc <- FindClusters(seurat_merged_qc, resolution = 0.5)
head(Idents(seurat_merged_qc), 5)

#Running UMAP
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
  path = "......../Github project/PLOTS",
  width = 20,
  height = 10,
  dpi = 300
)

saveRDS(seurat_merged_qc,file = "........./Github project/seurat_merged_qc.rds")

#FindMarkers function 
seurat_merged_qc <- JoinLayers(seurat_merged_qc)
seurat_markers <- FindAllMarkers(
  seurat_merged_qc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

class(seurat_markers)
colnames(seurat_markers)
head(seurat_markers)
tail(seurat_markers)
saveRDS(seurat_markers,file = "......../Github project/seurat_markers.rds")

# Marker identification will be added in the next commit
# TO DO: Identify top markers genes per cluster


