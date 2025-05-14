# Integration 

library(Seurat)

seu <- readRDS("seu_dr.rds")

Seurat::DimPlot(seu, reduction = "umap", 
                group.by = "orig.ident")



### Split by samples and integrate, then rejoin

seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)

seu <- Seurat::IntegrateLayers(object = seu, method = CCAIntegration,
                               orig.reduction = "pca",
                               new.reduction = "integrated.cca",
                               verbose = FALSE)

# re-join layers after integration
seu[["RNA"]] <- JoinLayers(seu[["RNA"]])


### Run UMAP again with appropriate "reduciton"


seu <- RunUMAP(seu, dims = 1:25, reduction = "integrated.cca")

Seurat::DimPlot(seu, reduction = "umap", group.by = "orig.ident")


# Check current identities
table(Idents(seu))

# If needed, set identities to your samples
Idents(seu) <- "orig.ident"

# Then you can plot without group.by
Seurat::DimPlot(seu, reduction = "umap")

saveRDS(seu, "seu_day2-2.rds")
