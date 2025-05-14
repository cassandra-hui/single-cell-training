## Trajectory

##############
# Keep in mind to ONLY do a trajectory analysis on data with a time aspect.
# You can force your data to make a trajectory but if you are not looking at 
# developmental stages, time points, or something else similar then it is not
# "biologially real".
################




library(SingleCellExperiment)
library(scater)
library(slingshot)
library(ggplot2)
library(ggbeeswarm)
library(Seurat)


deng_SCE <- readRDS("deng-reads.rds")
deng_SCE

# Change level from defualy alhabetical to time order
deng_SCE$cell_type2 <- factor(deng_SCE$cell_type2,
                              levels = c("zy",
                                         "early2cell",
                                         "mid2cell",
                                         "late2cell",
                                         "4cell",
                                         "8cell",
                                         "16cell",
                                         "earlyblast",
                                         "midblast",
                                         "lateblast"))


# Run a PCA
deng_SCE <- scater::runPCA(deng_SCE, ncomponents = 50)

# Asses PCA and store results
pca <- SingleCellExperiment::reducedDim(deng_SCE, "PCA")
head(pca)

# Add PCA data to SCE object
deng_SCE$PC1 <- pca[, 1]
deng_SCE$PC2 <- pca[, 2]

ggplot(as.data.frame(colData(deng_SCE)), aes(x = PC1, y = PC2, color = cell_type2)) +
  geom_point(size=2, shape=20) +
  theme_classic() +
  xlab("PC1") + ylab("PC2") + ggtitle("PC biplot")


deng_SCE$pseudotime_PC1 <- rank(deng_SCE$PC1)  # rank cells by their PC1 score

ggplot(as.data.frame(colData(deng_SCE)), aes(x = pseudotime_PC1, y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("PC1") + ylab("Timepoint") +
  ggtitle("Cells ordered by first principal component")

# Using slingshot for trajectory
sce <- slingshot::slingshot(deng_SCE, reducedDim = 'PCA')




## Custom fucnstion to plot the PCA in slingshot object
PCAplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
  # set palette for factorial variables
  palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  # set palette for numeric variables
  paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  # extract pca from SingleCellExperiment object
  pca <- SingleCellExperiment::reducedDims(sce)$PCA
  
  if(is.null(variable)){
    col <- "black"
  }
  if(is.character(variable)){
    variable <- as.factor(variable)
  }
  if(is.factor(variable)){
    colpal <- palf(length(levels(variable)))
    colors <- colpal[variable]
  }
  if(is.numeric(variable)){
    colpal <- paln(50)
    colors <- colpal[cut(variable,breaks=50)]
  }
  
  # draw the plot
  plot(pca, bg = colors, pch = 21)
  # draw lines
  if(draw_lines){
    lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
  }
  # add legend
  if(legend & is.factor(variable)){
    legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)
    
  }
}



UMAPplot_slingshot <- function(sce, draw_lines = TRUE, variable = NULL, legend = FALSE, ...){
  # set palette for factorial variables
  palf <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))
  # set palette for numeric variables
  paln <- colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))
  # extract pca from SingleCellExperiment object
  pca <- SingleCellExperiment::reducedDims(sce)$UMAP
  
  if(is.null(variable)){
    col <- "black"
  }
  if(is.character(variable)){
    variable <- as.factor(variable)
  }
  if(is.factor(variable)){
    colpal <- palf(length(levels(variable)))
    colors <- colpal[variable]
  }
  if(is.numeric(variable)){
    colpal <- paln(50)
    colors <- colpal[cut(variable,breaks=50)]
  }
  
  # draw the plot
  plot(pca, bg = colors, pch = 21)
  # draw lines
  if(draw_lines){
    lines(slingshot::SlingshotDataSet(sce), lwd = 2, ... )
  }
  # add legend
  if(legend & is.factor(variable)){
    legend("bottomright", pt.bg = colpal,legend = levels(variable),pch=21)
    
  }
}



PCAplot_slingshot(sce, variable = sce$slingPseudotime_1, draw_lines = TRUE)

ggplot(as.data.frame(colData(deng_SCE)), aes(x = sce$slingPseudotime_1,
                                             y = cell_type2,
                                             colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")



############################################
## Using Seurat
############################################

# Create object
gcdata <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(deng_SCE),
                                     project = "slingshot")
# Noramilze
gcdata <- Seurat::NormalizeData(object = gcdata,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
# Variable Features
gcdata <- Seurat::FindVariableFeatures(object = gcdata,
                                       mean.function = ExpMean,
                                       dispersion.function = LogVMR)
# Scale
gcdata <- Seurat::ScaleData(object = gcdata,
                            do.center = T,
                            do.scale = F)
# PCA
gcdata <- Seurat::RunPCA(object = gcdata,
                         pc.genes = gcdata@var.genes)
# Find Neighbors
gcdata <- Seurat::FindNeighbors(gcdata,
                                reduction = "pca",
                                dims = 1:5)

# clustering with resolution of 0.6
gcdata <- Seurat::FindClusters(object = gcdata,
                               resolution = 0.6)

# Add clusters to SCE object
deng_SCE$Seurat_clusters <- as.character(Idents(gcdata))  # go from factor to character

# Run Slingshot tragectory
sce <- slingshot::slingshot(deng_SCE,
                            clusterLabels = 'Seurat_clusters',
                            reducedDim = 'PCA',
                            start.clus = "2")


SlingshotDataSet(sce)
# class: SlingshotDataSet 
# 
# Samples Dimensions
# 268         50
# 
# lineages: 2 
# Lineage1: 2  4  0  5  3  
# Lineage2: 2  4  1  
# 
# curves: 2 
# Curve1: Length: 425.88	Samples: 234.62
# Curve2: Length: 340.99	Samples: 132.38


PCAplot_slingshot(sce, variable = sce$slingPseudotime_2)


ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_1 = sce$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = cell_type2,
           colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

ggplot(data.frame(cell_type2 = deng_SCE$cell_type2,
                  slingPseudotime_2 = sce$slingPseudotime_2),
       aes(x = slingPseudotime_2, y = cell_type2,
           colour = cell_type2)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")


PCAplot_slingshot(sce,
                  variable = deng_SCE$Seurat_clusters,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)

PCAplot_slingshot(sce,
                  variable = deng_SCE$cell_type2,
                  type = 'lineages',
                  col = 'black',
                  legend = TRUE)



sce <- slingshot::slingshot(deng_SCE,
                            clusterLabels = 'Seurat_clusters',
                            reducedDim = 'PCA',
                            end.clus = c("0", "3", "5")) ## check which would be the best according to bio




###################################
# Testing on Caroline's data
###################################

seu <- readRDS("my_objects/seu_day2-4.rds")
seu.sce <- as.SingleCellExperiment(seu) # Switch to SCE


sling <- slingshot::slingshot(seu.sce, reducedDim = 'PCA')

SlingshotDataSet(sling)


PCAplot_slingshot(sling, variable = sling$slingPseudotime_1)


seu.sce$sample

ggplot(data.frame(sample = seu.sce$sample,
                  slingPseudotime_1 = sling$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = sample,
           colour = sample)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

## Try UMAP
reducedDims(seu.sce)



sling2 <- slingshot::slingshot(seu.sce, reducedDim = 'UMAP')

SlingshotDataSet(sling2)


UMAPplot_slingshot(sling2, variable = sling2$slingPseudotime_1)


seu.sce$sample

ggplot(data.frame(sample = seu.sce$sample,
                  slingPseudotime_1 = sling2$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = sample,
           colour = sample)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")


###################################
# Testing on Sal's data
###################################

seu <- readRDS("../../Projects/Sal_Baker/scRNAseq_2024/obj.rds")
seu.sce <- as.SingleCellExperiment(seu) # Switch to SCE


sling <- slingshot::slingshot(seu.sce, reducedDim = 'PCA')

SlingshotDataSet(sling)


PCAplot_slingshot(sling, variable = sling$slingPseudotime_1)


seu.sce$sample

ggplot(data.frame(sample = seu.sce$sample,
                  slingPseudotime_1 = sling$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = sample,
           colour = sample)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

## Try UMAP
reducedDims(seu.sce)



sling2 <- slingshot::slingshot(seu.sce, reducedDim = 'UMAP')

SlingshotDataSet(sling2)


UMAPplot_slingshot(sling2, variable = sling2$slingPseudotime_1)

seu.sce$sample

ggplot(data.frame(sample = seu.sce$sample,
                  slingPseudotime_1 = sling2$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = sample,
           colour = sample)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")


###################################
# Testing on Misty's data
###################################

seu <- readRDS("../../Projects/Misty_Riddle/scRNAseq_2024/Seurat_Obj_UMAP_02.03.2025.RDS")
seu.sce <- as.SingleCellExperiment(seu) # Switch to SCE


sling <- slingshot::slingshot(seu.sce, reducedDim = 'PCA')

SlingshotDataSet(sling)


PCAplot_slingshot(sling, variable = sling$slingPseudotime_1)


seu.sce$sample

ggplot(data.frame(sample = seu.sce$sample,
                  slingPseudotime_1 = sling$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = sample,
           colour = sample)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")

## Try UMAP
reducedDims(seu.sce)



sling2 <- slingshot::slingshot(seu.sce, reducedDim = 'UMAP')

SlingshotDataSet(sling2)

UMAPplot_slingshot(sling2, variable = sling2$slingPseudotime_1)

# Open a PNG device before plotting
png("Pseudotime_Trajectory.png", width = 1200, height = 1000, res = 300)
UMAPplot_slingshot(sling2, variable = sling2$slingPseudotime_1)
# Close the PNG device to save the file
dev.off()

# Try PDF which handles margins better
pdf("Pseudotime_Trajectory.pdf", width = 8, height = 6)
UMAPplot_slingshot(sling2, variable = sling2$slingPseudotime_1)
dev.off()



seu.sce$sample

ggplot(data.frame(sample = seu.sce$sample,
                  slingPseudotime_1 = sling2$slingPseudotime_1),
       aes(x = slingPseudotime_1, y = sample,
           colour = sample)) +
  ggbeeswarm::geom_quasirandom(groupOnX = FALSE) +
  theme_classic() +
  xlab("Slingshot pseudotime") + ylab("Timepoint") +
  ggtitle("Cells ordered by Slingshot pseudotime")


