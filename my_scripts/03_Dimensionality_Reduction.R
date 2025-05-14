## Dimensionality Reduction

library(Seurat)

seu <- readRDS("seu_day1-3_Caroline.rds")

head(seu@meta.data)

seu$orig.ident <- seu@meta.data$sample
table(seu$orig.ident)

## PCA

seu <- Seurat::RunPCA(seu)
Seurat::DimPlot(seu, reduction = "pca", group.by = "orig.ident")
Seurat::FeaturePlot(seu, reduction = "pca", features = "percent.mito")


Seurat::FeaturePlot(seu, reduction = "pca", features = "C3")

Seurat::DimHeatmap(seu, dims = 1:12, cells = 500, balanced = TRUE)

Seurat::ElbowPlot(seu, ndims = 40)

## Using too many PCs doesn't matter to much, using too few will give you less information to work with

seu <- Seurat::RunUMAP(seu, dims = 1:25)

Seurat::DimPlot(seu, reduction = "umap", group.by = "orig.ident")

## Change the number of neighbors 
## Defualt is 30

seu <- Seurat::RunUMAP(seu, dims = 1:25, n.neighbors = 5)
Seurat::DimPlot(seu, reduction = "umap",  group.by = "orig.ident")


## Change the numner of PCs


seu <- Seurat::RunUMAP(seu, dims = 1:5)
Seurat::DimPlot(seu, reduction = "umap",  group.by = "orig.ident")
## Way less clustered 

seu <- Seurat::RunUMAP(seu, dims = 1:50)
Seurat::DimPlot(seu, reduction = "umap",  group.by = "orig.ident")
## A little more spread out



## Return to before

seu <- Seurat::RunUMAP(seu, dims = 1:25)
Seurat::DimPlot(seu, reduction = "umap",  group.by = "orig.ident")


saveRDS(seu, "seu_dr.rds")
