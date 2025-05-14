# Clustering 
library(Seurat)
seu <- readRDS("seu_day2-2.rds")

seu <- Seurat::FindNeighbors(seu, dims = 1:25, reduction = "integrated.cca")


seu <- Seurat::FindClusters(seu, resolution = seq(0.1, 0.8, by=0.1))
head(seu@meta.data)


library(clustree)

clustree::clustree(seu@meta.data[,grep("RNA_snn_res", colnames(seu@meta.data))],
                   prefix = "RNA_snn_res.")


Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.1")
Seurat::DimPlot(seu, group.by = "RNA_snn_res.0.3")

